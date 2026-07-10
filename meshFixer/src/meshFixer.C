#include "meshFixer.H"
#include "fvm.H"
#include "fvc.H"
#include "syncTools.H"
#include "polyMesh.H"
#include "volPointInterpolation.H"
#include "mathematicalConstants.H"
#include "unitConversion.H"
#include <iomanip>

namespace Foam
{

meshFixer::meshFixer(fvMesh& mesh)
:
    mesh_(mesh),
    runTime_(const_cast<Time&>(mesh.time()))
{
    IOdictionary dict
    (
        IOobject
        (
            "meshFixerDict",
            runTime_.system(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const dictionary& smootherDict = dict.subDict("smoother");
    totalIterations_ = smootherDict.lookupOrDefault<label>("totalIterations", 10);
    writeEvery_ = smootherDict.lookupOrDefault<label>("writeEvery", 10);
    if(writeEvery_ == 0){
        writeEvery_ = 99999;
    }
    writeIntermediate_ = smootherDict.lookupOrDefault<bool>("writeIntermediate", true);  
    writeDisplacement_ = smootherDict.lookupOrDefault<bool>("writeDisplacement", false); 
    localRelaxation_   = smootherDict.lookupOrDefault<scalar>("localRelaxation", 1.0);

    qualityControlDict_ = dict.subDict("meshQualityControls");
    const dictionary& localRelaxDict = smootherDict.subDict("localAdaptation");
    movementRelaxIters_ = localRelaxDict.lookupOrDefault<label>("relaxationIterations", 1000);
    movementRelaxation_ = localRelaxDict.lookupOrDefault<scalar>("relaxationFactor", 0.5);
    relativeQualityControlDict_ = localRelaxDict.subDict("meshQualityControls");
    const dictionary& wiggleDict = smootherDict.subDict("wiggleSmooth");
    wiggleQualityControlDict_ = wiggleDict.subDict("meshQualityControls");

    const dictionary& moDict = smootherDict.subDict("multiObjective");
    wVol_        = moDict.lookupOrDefault<scalar>("weightVolRegularizer", 0.1);
    wLap_        = moDict.lookupOrDefault<scalar>("weightLaplacian", 0.05);
    wExp_        = moDict.lookupOrDefault<scalar>("weightExpansion", 0.5);
    maxExpRatio_ = moDict.lookupOrDefault<scalar>("maxExpansionRatio", 1.5);
    maxARShield_ = moDict.lookupOrDefault<scalar>("maxAspectRatioShield", 5.0);

    const dictionary& qualDict = smootherDict.subDict("quality");
    qGood_ = qualDict.lookupOrDefault<scalar>("qGood", 60.0);
    qBad_  = qualDict.lookupOrDefault<scalar>("qBad", 75.0);
    sGood_ = qualDict.lookupOrDefault<scalar>("sGood", 2.0); 
    sBad_  = qualDict.lookupOrDefault<scalar>("sBad", 4.0);  
    orthoSpring_ = qualDict.lookupOrDefault<scalar>("orthoMaxSpring", 60.0); 
    skewSpring_ = qualDict.lookupOrDefault<scalar>("skewMaxSpring", 3.0);  

    const dictionary& diffDict = smootherDict.subDict("diffusion");
    diffConstant_ = diffDict.lookupOrDefault<scalar>("baseConstant", 1e-3);

    nSmoothLayers_  = smootherDict.lookupOrDefault<label>("nSmoothLayers", 2);
    laplacianIters_ = smootherDict.lookupOrDefault<label>("laplacianIters", 5);
    localLaplacianRelaxation_ = smootherDict.lookupOrDefault<scalar>("laplacianRelaxation", 0.5); 

    wiggleFactor = wiggleDict.lookupOrDefault<scalar>("relStepLength", 0.1);
    wiggleInitialized = false;
    gradientDelta = wiggleDict.lookupOrDefault<scalar>("angleStep", 15)/180*Foam::constant::mathematical::pi;
    randomLambLimit = wiggleDict.lookupOrDefault<scalar>("randomLambLimit", 0.1);
    //rndGen = uniformGeneratorOp();
    wiggleSurroundingLayers_ = wiggleDict.lookupOrDefault<label>("surroundingLayers", 0);
    complexReward_ = wiggleDict.lookupOrDefault<bool>("complexReward", false);
    selfFactor_ = wiggleDict.lookupOrDefault<scalar>("selfFactor", 2);

    const dictionary& workflowDict = dict.subDict("workflow");
    discreteStep_ = workflowDict.lookupOrDefault<bool>("discreteAdjustment", true);  
    laplacianStep_ = workflowDict.lookupOrDefault<bool>("laplacianSmooth", true);  
    wiggleStep_ = workflowDict.lookupOrDefault<bool>("wiggleSmooth", true);  
    verboseWiggle_ = false;

}

void meshFixer::identifyBoundaries()
{
    isBoundaryPoint_.setSize(mesh_.nPoints(), false);
    const polyBoundaryMesh& bMesh = mesh_.boundaryMesh();

    forAll(bMesh, patchi)
    {
        const polyPatch& p = bMesh[patchi];
        if (!p.coupled() && p.size() > 0)
        {
            forAll(p, localFacei)
            {
                const face& f = p[localFacei];
                forAll(f, i)
                {
                    isBoundaryPoint_[f[i]] = true;
                }
            }
        }
    }
    syncTools::syncPointList(mesh_, isBoundaryPoint_, orEqOp<bool>(), false);
}

// ------------------------------------------------------------------------- //
// Mesh Metric Functions (Static Arrays)
// ------------------------------------------------------------------------- //

scalarField meshFixer::calculateFaceNonOrtho() const
{
    scalarField nonOrtho(mesh_.nFaces(), 0.0);
    const surfaceVectorField& Sf = mesh_.Sf();
    const volVectorField& C = mesh_.C();

    for (label facei = 0; facei < mesh_.nInternalFaces(); ++facei)
    {
        label own = mesh_.faceOwner()[facei];
        label nei = mesh_.faceNeighbour()[facei];
        vector d = C[nei] - C[own];
        vector n = Sf[facei] / (mag(Sf[facei]) + VSMALL);
        
        scalar cosAngle = (d & n) / (mag(d) + VSMALL);
        nonOrtho[facei] = acos(min(1.0, max(-1.0, cosAngle))) * 180.0 / Foam::constant::mathematical::pi;
    }
    return nonOrtho;
}

// ------------------------------------------------------------------------- //
// Helper: Calculate Skewness (Corrected for Boundary Faces)
// ------------------------------------------------------------------------- //
scalar meshFixer::calculateFaceSkewness
(
    const pointField& p,
    const label facei
) const
{
    const face& f = mesh_.faces()[facei];
    vector s = f.areaNormal(p);
    scalar magS = mag(s);
    
    if (magS < VSMALL) return GREAT;

    const label own = mesh_.faceOwner()[facei];
    scalar ownVol; vector ownCc;
    calculateCellVolumeAndCentre(p, own, ownVol, ownCc);
    
    vector fc = f.centre(p);
    vector dOwn = fc - ownCc;

    if (mesh_.isInternalFace(facei))
    {
        // --- Internal Face Skewness ---
        const label nei = mesh_.faceNeighbour()[facei];
        scalar neiVol; vector neiCc;
        calculateCellVolumeAndCentre(p, nei, neiVol, neiCc);
        
        vector d = neiCc - ownCc;
        
        // Intersect cell-to-cell line with face plane
        scalar t = (s & dOwn) / ((s & d) + VSMALL);
        vector pInt = ownCc + t * d;

        return mag(fc - pInt) / (mag(d) + VSMALL);
    }
    else
    {
        // --- Boundary Face Skewness ---
        // Project owner cell center onto the face plane along the face normal
        vector n = s / magS;
        vector pInt = ownCc + n * (dOwn & n); // Projection point
        
        return mag(fc - pInt) / (mag(dOwn) + VSMALL);
    }
}

scalar meshFixer::calculateCellAspectRatio(label celli) const
{
    const labelList& cFaces = mesh_.cells()[celli];
    const faceList& faces = mesh_.faces();
    const pointField& pts = mesh_.points();

    scalar maxLenSqr = -VGREAT;
    scalar minLenSqr = VGREAT;

    for (label facei : cFaces)
    {
        const face& f = faces[facei];
        for (label i = 0; i < f.size(); ++i)
        {
            label pt1 = f[i];
            label pt2 = f[(i + 1) % f.size()];
            scalar lenSqr = magSqr(pts[pt1] - pts[pt2]);
            
            maxLenSqr = max(maxLenSqr, lenSqr);
            minLenSqr = min(minLenSqr, lenSqr);
        }
    }
    
    return sqrt(maxLenSqr / (minLenSqr + VSMALL));
}

scalar meshFixer::getMinEdgeLength(label ptI) const
{
    const labelList& pEdges = mesh_.pointEdges()[ptI];
    scalar minLen = VGREAT;
    for (label edgeI : pEdges)
    {
        const edge& e = mesh_.edges()[edgeI];
        minLen = min(minLen, mag(mesh_.points()[e.start()] - mesh_.points()[e.end()]));
    }
    return minLen;
}

// ------------------------------------------------------------------------- //
// Dynamic Mesh Metric Functions (For the Bisection Engine)
// ------------------------------------------------------------------------- //

void meshFixer::calculateCellVolumeAndCentre
(
    const pointField& p,
    const label celli,
    scalar& volume,
    vector& centre 
) const
{
    const cell& cFaces = mesh_.cells()[celli];
    
    vector cEst = Zero;
    for (const label facei : cFaces)
    {
        cEst += mesh_.faces()[facei].centre(p);
    }
    cEst /= cFaces.size();

    centre = Zero;
    volume = 0.0;

    for (const label facei : cFaces)
    {
        const face& f = mesh_.faces()[facei];
        vector fc = f.centre(p);
        vector fa = f.areaNormal(p);
        
        scalar pyr3Vol = fa & (fc - cEst);
        if (mesh_.faceOwner()[facei] != celli)
        {
            pyr3Vol = -pyr3Vol;
        }

        volume += pyr3Vol;
        centre += pyr3Vol * (0.75 * fc + 0.25 * cEst);
    }

    if (mag(volume) > VSMALL)
    {
        centre /= volume;
    }
    else
    {
        centre = cEst;
    }
    volume /= 3.0;
}

scalar meshFixer::calculateFaceNonOrtho
(
    const pointField& p,
    const label facei
) const
{
    const face& f = mesh_.faces()[facei];
    vector s = f.areaNormal(p);
    scalar magS = mag(s);
    
    if (magS < VSMALL) return 0.0; // Perfectly orthogonal if no area

    const label own = mesh_.faceOwner()[facei];
    scalar ownVol; vector ownCc;
    calculateCellVolumeAndCentre(p, own, ownVol, ownCc);

    vector d(Zero);
    if (mesh_.isInternalFace(facei))
    {
        const label nei = mesh_.faceNeighbour()[facei];
        scalar neiVol; vector neiCc;
        calculateCellVolumeAndCentre(p, nei, neiVol, neiCc);
        d = neiCc - ownCc;
    }
    else
    {
        // EXACT OPENFOAM MATCH:
        // Original motionSmoother ignores non-ortho on uncoupled boundaries.
        // If you want the exact same results, return 0.0 for uncoupled boundaries.
        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
        const label patchi = pbm.whichPatch(facei);
        if (!pbm[patchi].coupled())
        {
            return 0.0; 
        }
        
        // If it is coupled, you would theoretically need the neighbour cell center.
        // Since we are doing local morphing without MPI syncing centers:
        d = f.centre(p) - ownCc;
    }

    scalar magD = mag(d);
    if (magD < VSMALL) return 0.0;

    // The dot product logic exactly as implemented in polyMeshGeometry.C
    scalar dDotS = (d & s) / (magD * magS + VSMALL);
    dDotS = max(-1.0, min(1.0, dDotS)); 
    
    return radToDeg(acos(dDotS));
}

// ------------------------------------------------------------------------- //
// Morphing & Solver Logic
// ------------------------------------------------------------------------- //

void meshFixer::calcLocalCorrections(pointVectorField& localDisp)
{
    vectorField& ptDisp = localDisp.primitiveFieldRef();
    scalarField sumWeights(mesh_.nPoints(), 0.0);
    scalarField minEdgeLength(mesh_.nPoints(), 1e20);
    
    isBadPoint_.setSize(mesh_.nPoints(), false);
    isBadPoint_ = false;

    scalarField nonOrthoAngle = calculateFaceNonOrtho();
    scalarField skewness = calculateFaceSkewness();



    for (label facei = 0; facei < mesh_.nInternalFaces(); ++facei)
    {
        if (true) 
        {

            bool badFace = (nonOrthoAngle[facei] > qBad_ || skewness[facei] > sBad_);

            label own = mesh_.faceOwner()[facei];
            label nei = mesh_.faceNeighbour()[facei];
            
            point idealCenter = 0.5*(mesh_.C()[nei]+mesh_.C()[own]);
            vector idealNormal = mesh_.C()[nei] - mesh_.C()[own];
            idealNormal = idealNormal / mag(idealNormal);

            const face& f = mesh_.faces()[facei];

            //calculate face center
            point faceCenter(0,0,0);
            label pointCount = 0;
            scalar minEdgeLengthFace = 1e20;
            
            forAll(f, i)
            {
                label ptI = f[i];
                faceCenter += mesh_.points()[ptI];
                pointCount++;
                //calculate min edge length

                Foam::edge edge_ = f.edge(i);
                scalar edgeLength = mag(
                    mesh_.points()[edge_.a()]
                    - mesh_.points()[edge_.b()]
                );
                minEdgeLengthFace = min(minEdgeLengthFace,edgeLength);
            }

            faceCenter = faceCenter/scalar(pointCount);

            vector uniformMovement = idealCenter - faceCenter;

            forAll(f, i)
            {
                label ptI = f[i];
                if (!isBoundaryPoint_[ptI])
                {
                    vector pointMove = uniformMovement;
                    
                    vector pointToCenter = faceCenter - mesh_.points()[ptI];
                    vector pointToPlane = idealNormal * (idealNormal & pointToCenter);

                    pointMove += pointToPlane;
                    
                    ptDisp[ptI] += pointMove;
                    sumWeights[ptI] += 1.0;
                    minEdgeLength[ptI] = min(minEdgeLength[ptI],minEdgeLengthFace);
                    isBadPoint_[ptI] = isBadPoint_[ptI] || badFace;
                }
            }
        }
    }

    syncTools::syncPointList(mesh_, ptDisp, plusEqOp<vector>(), vector::zero);
    syncTools::syncPointList(mesh_, sumWeights, plusEqOp<scalar>(), 0.0);
    syncTools::syncPointList(mesh_, isBadPoint_, orEqOp<bool>(), false);
    syncTools::syncPointList(mesh_, minEdgeLength, minEqOp<scalar>(), 0.0);

    scalar maxDisplacement = 0.0;

    forAll(ptDisp, i)
    {
        if (isBadPoint_[i])
        {
            ptDisp[i] = localRelaxation_ * (ptDisp[i] / sumWeights[i]);
            scalar edgeLimit = 0.3;
            if(mag(ptDisp[i]) > edgeLimit*minEdgeLength[i])
            {
                ptDisp[i] = ptDisp[i] * (edgeLimit*minEdgeLength[i]/mag(ptDisp[i]));
            }
            maxDisplacement = max(maxDisplacement,mag(ptDisp[i]));
        }
    }
    reduce(maxDisplacement,maxOp<scalar>());
    Info << "Max displacement from localRelaxer " << maxDisplacement << endl;
}

void meshFixer::solveGlobalDiffusion(const pointVectorField& localDisp)
{
    volVectorField cellSource
    (
        IOobject("cellSource", mesh_.time().timeName(), mesh_, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh_,
        dimensionedVector(dimLength, vector::zero)
    );

    const labelListList& cellPoints = mesh_.cellPoints();
    forAll(cellPoints, celli)
    {
        vector avgDisp = vector::zero;
        const labelList& cPts = cellPoints[celli];
        for (label ptI : cPts)
        {
            avgDisp += localDisp[ptI];
        }
        cellSource[celli] = avgDisp / cPts.size();
    }

    volScalarField gamma
    (
        IOobject("gamma", mesh_.time().timeName(), mesh_, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh_,
        dimensionedScalar("g", dimArea, 0.0)
    );

    volScalarField cellQualityDamping
    (
        IOobject("cellQualityDamping", mesh_.time().timeName(), mesh_, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh_,
        dimensionedScalar("qDamp", dimless, 1.0)
    );

    scalarField nonOrthoAngle = calculateFaceNonOrtho();
    scalarField skewness = calculateFaceSkewness();

    forAll(mesh_.cells(), celli)
    {
        gamma[celli] = mesh_.V()[celli] * diffConstant_;

        scalar maxCellAngle = 0.0;
        scalar maxCellSkew = 0.0; 
        const labelList& cFaces = mesh_.cells()[celli];
        for (label facei : cFaces)
        {
            if (facei < mesh_.nInternalFaces())
            {
                maxCellAngle = max(maxCellAngle, nonOrthoAngle[facei]);
                maxCellSkew = max(maxCellSkew, skewness[facei]); 
            }
        }

        scalar qTildeOrtho = max(0.0, min(1.0, (qBad_ - maxCellAngle) / (qBad_ - qGood_ + VSMALL)));
        scalar qTildeSkew  = max(0.0, min(1.0, (sBad_ - maxCellSkew)  / (sBad_ - sGood_ + VSMALL)));
        
        scalar qTilde = min(qTildeOrtho, qTildeSkew);
        cellQualityDamping[celli] = 1.0 + (qTilde * 1e3); 
    }

    volVectorField cellDisp
    (
        IOobject("cellDisp", mesh_.time().timeName(), mesh_, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh_,
        dimensionedVector(dimLength, vector::zero),
        "fixedValue" 
    );

    forAll(cellDisp.boundaryField(), patchi)
    {
        cellDisp.boundaryFieldRef()[patchi] == vector::zero;
    }

    fvVectorMatrix dispEqn
    (
        fvm::Sp(cellQualityDamping, cellDisp)
      - fvm::laplacian(gamma, cellDisp)
     == cellSource
    );

    dispEqn.solve();

    pointVectorField finalPointDisp
    (
        IOobject
        (
            "pointDisplacement",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        volPointInterpolation::New(mesh_).interpolate(cellDisp)
    );

    forAll(finalPointDisp, ptI)
    {
        if (isBoundaryPoint_[ptI])
        {
            finalPointDisp[ptI] = vector::zero;
        }
        else if (isBadPoint_[ptI])
        {
            finalPointDisp[ptI] = localDisp[ptI];
        }
    }

    if (writeDisplacement_)
    {
        finalPointDisp.write();
    }

    // You can choose which morph algorithm to use here. We will use the new relative one.
    //applyRelativeGuardedMorph(finalPointDisp);
    applyQualityGuardedMorph(applyRelativeGuardedMorph(finalPointDisp,false));
}

void meshFixer::applyQualityGuardedMorph(const pointVectorField& finalPointDisp)
{
    Info<< "   Applying legacy quality-guarded morphing..." << endl;

    labelList emptyLabels;
    indirectPrimitivePatch dummyPatch
    (
        IndirectList<face>(mesh_.faces(), emptyLabels), 
        mesh_.points()
    );

    motionSmoother meshMover
    (
        mesh_,
        dummyPatch,
        emptyLabels,
        finalPointDisp, 
        qualityControlDict_
    );

    labelList checkFaces(identity(mesh_.nFaces()));
    List<labelPair> baffles;
    label nInitErrors = 0; 
    bool meshOk = false;
    label nScaleIterations = 10;

    for (label iter = 0; iter < 2 * nScaleIterations; iter++)
    {
        if (iter == nScaleIterations)
        {
            meshMover.setErrorReduction(0.0);
        }

        meshOk = meshMover.scaleMesh
        (
            checkFaces, 
            baffles, 
            true, 
            nInitErrors
        );

        if (meshOk)
        {
            Info<< "   -> Displacement applied fully without quality violations." << endl;
            break;
        }
    }

    if (!meshOk)
    {
        Info<< "   -> Notice: Displacement scaled back to prevent bad cells." << endl;
    }

    meshMover.correct();
}

pointVectorField meshFixer::applyRelativeGuardedMorph(const pointVectorField& finalPointDisp, bool move)
{
    const pointField& p0 = mesh_.points();
    const vectorField& disp = finalPointDisp.primitiveField();
    
    label pendulumState = 0;
    scalar limit = 0.01;
    scalar lamFactor = 1;
    scalar lamSign = 1;

    // --- 1. Read BOTH Quality Metrics Limits ---
    
    // Default Limits (from localAdaptation)
    const scalar minPyrVolDef  = relativeQualityControlDict_.lookupOrDefault<scalar>("minVol", 1e-13);
    const scalar maxConcaveDef = relativeQualityControlDict_.lookupOrDefault<scalar>("maxConcave", 80.0);
    const scalar maxIntSkewDef = relativeQualityControlDict_.lookupOrDefault<scalar>("maxInternalSkewness", 4.0);
    const scalar maxBounSkewDef= relativeQualityControlDict_.lookupOrDefault<scalar>("maxBoundarySkewness", 20.0);
    const scalar minWeightDef  = relativeQualityControlDict_.lookupOrDefault<scalar>("minFaceWeight", 0.0001);
    const scalar minVolRatioDef= relativeQualityControlDict_.lookupOrDefault<scalar>("minVolRatio", 0.01);
    const scalar minTwistDef   = relativeQualityControlDict_.lookupOrDefault<scalar>("minTwist", 0.05);
    const scalar qBadDef       = relativeQualityControlDict_.lookupOrDefault<scalar>("maxNonOrtho", 70.0);

    // Stricter Limits (from wiggleSmooth)
    const scalar minPyrVolWig  = wiggleQualityControlDict_.lookupOrDefault<scalar>("minVol", 1e-13);
    const scalar maxConcaveWig = wiggleQualityControlDict_.lookupOrDefault<scalar>("maxConcave", 80.0);
    const scalar maxIntSkewWig = wiggleQualityControlDict_.lookupOrDefault<scalar>("maxInternalSkewness", 4.0);
    const scalar maxBounSkewWig= wiggleQualityControlDict_.lookupOrDefault<scalar>("maxBoundarySkewness", 20.0);
    const scalar minWeightWig  = wiggleQualityControlDict_.lookupOrDefault<scalar>("minFaceWeight", 0.0001);
    const scalar minVolRatioWig= wiggleQualityControlDict_.lookupOrDefault<scalar>("minVolRatio", 0.01);
    const scalar minTwistWig   = wiggleQualityControlDict_.lookupOrDefault<scalar>("minTwist", 0.05);
    const scalar qBadWig       = wiggleQualityControlDict_.lookupOrDefault<scalar>("maxNonOrtho", 70.0);

    const scalar minPyrVolTight   = max(minPyrVolDef, minPyrVolWig);
    const scalar maxConcaveTight  = min(maxConcaveDef, maxConcaveWig);
    const scalar maxIntSkewTight  = min(maxIntSkewDef, maxIntSkewWig);
    const scalar maxBounSkewTight = min(maxBounSkewDef, maxBounSkewWig);
    const scalar minWeightTight   = max(minWeightDef, minWeightWig);
    const scalar minVolRatioTight = max(minVolRatioDef, minVolRatioWig);
    const scalar minTwistTight    = max(minTwistDef, minTwistWig);
    const scalar qBadTight        = min(qBadDef, qBadWig);

    // 2. Identify Active Topology
    labelHashSet activePointsSet;
    forAll(disp, pointi)
    {
        if (mag(disp[pointi]) > VSMALL)
        {
            activePointsSet.insert(pointi);
        }
    }

    // Expand to connected faces and cells
    labelHashSet activeFaces;
    labelHashSet activeCells;
    const labelListList& pointFaces = mesh_.pointFaces();
    
    for (const label pointi : activePointsSet)
    {
        for (const label facei : pointFaces[pointi])
        {
            activeFaces.insert(facei);
            activeCells.insert(mesh_.faceOwner()[facei]);
            if (mesh_.isInternalFace(facei))
            {
                activeCells.insert(mesh_.faceNeighbour()[facei]);
            }
        }
    }

    label initialFaceCount = activeFaces.size();

    //All faces that belong to cells that are Part of faces that move must also be considered (cell center moves)
    for (const label celli : activeCells)
    {
        const cell& cFaces = mesh_.cells()[celli];
        for (const label facei : cFaces)
        {
            if (mesh_.isInternalFace(facei))
            {
                activeFaces.insert(facei);
            }
        }
    }

    // --- Calculate newly added faces and reduce across all cores ---
    label localAddedFaces = activeFaces.size() - initialFaceCount;
    label globalAddedFaces = returnReduce(localAddedFaces, sumOp<label>());

    Info << "   -> Faces added due to active cell expansion: " << globalAddedFaces << endl;

    // 3. Pre-calculate Old State & Per-Face Limits
    Map<scalar> orthoOld(activeFaces.size());
    Map<scalar> skewOld(activeFaces.size());
    Map<scalar> pyrOwnOld(activeFaces.size());
    Map<scalar> pyrNeiOld(activeFaces.size());
    Map<scalar> concOld(activeFaces.size());
    Map<scalar> weightOld(activeFaces.size());
    Map<scalar> ratioOld(activeFaces.size());
    Map<scalar> twistOld(activeFaces.size());

    // Step 3a: Cache initial states first so we can use them for limit sub-checks
    for (const label facei : activeFaces)
    {
        orthoOld.insert(facei, calculateFaceNonOrtho(p0, facei));
        skewOld.insert(facei, calculateFaceSkewness(p0, facei));
        pyrOwnOld.insert(facei, calculateFacePyramidVolume(p0, facei, true));
        if (mesh_.isInternalFace(facei)) {
            pyrNeiOld.insert(facei, calculateFacePyramidVolume(p0, facei, false));
        } else {
            pyrNeiOld.insert(facei, GREAT); // Boundaries don't have a neighbour pyramid
        }
        concOld.insert(facei, calculateFaceConcavity(p0, facei));
        weightOld.insert(facei, calculateInterpolationWeight(p0, facei));
        ratioOld.insert(facei, calculateVolumeRatio(p0, facei));
        twistOld.insert(facei, calculateFaceTwist(p0, facei));
    }

    // Per-face Dynamic Limit Maps
    Map<scalar> limPyrVol(activeFaces.size());
    Map<scalar> limConcave(activeFaces.size());
    Map<scalar> limIntSkew(activeFaces.size());
    Map<scalar> limBounSkew(activeFaces.size());
    Map<scalar> limWeight(activeFaces.size());
    Map<scalar> limRatio(activeFaces.size());
    Map<scalar> limTwist(activeFaces.size());
    Map<scalar> limOrtho(activeFaces.size());
    
    // Step 3b: Assign limits dynamically based on culpability and layer buffering
    for (const label facei : activeFaces)
    {
        // Check if any point on this face is currently being wiggled
        bool useWiggle = false;
        const face& f = mesh_.faces()[facei];
        forAll(f, i)
        {
            if (wigglePoints.found(f[i]))
            {
                useWiggle = true;
                break;
            }
        }

        // Check if the face is completely clean against the STRICT limits
        bool isCleanFace = !(
            (orthoOld[facei] > qBadWig) ||
            (mesh_.isInternalFace(facei) && skewOld[facei] > maxIntSkewWig) ||
            (!mesh_.isInternalFace(facei) && skewOld[facei] > maxBounSkewWig) ||
            (pyrOwnOld[facei] <= minPyrVolWig) ||
            (pyrNeiOld[facei] <= minPyrVolWig) ||
            (concOld[facei] > maxConcaveWig) ||
            (weightOld[facei] < minWeightWig) ||
            (ratioOld[facei] < minVolRatioWig) ||
            (twistOld[facei] < minTwistWig)
        );

        bool isDirtyFace = (
            (orthoOld[facei] > qBadDef) ||
            (mesh_.isInternalFace(facei) && skewOld[facei] > maxIntSkewDef) ||
            (!mesh_.isInternalFace(facei) && skewOld[facei] > maxBounSkewDef) ||
            (pyrOwnOld[facei] <= minPyrVolDef) ||
            (pyrNeiOld[facei] <= minPyrVolDef) ||
            (concOld[facei] > maxConcaveDef) ||
            (weightOld[facei] < minWeightDef) ||
            (ratioOld[facei] < minVolRatioDef) ||
            (twistOld[facei] < minTwistDef)
        );

        if (isCleanFace)
        {
            // It's a healthy buffer face (added via layers). 
            // Enforce strict limits on EVERYTHING so it doesn't break during the morph.
            limPyrVol.insert(facei, max(minPyrVolDef, minPyrVolWig));
            limConcave.insert(facei, min(maxConcaveDef, maxConcaveWig));
            limIntSkew.insert(facei, min(maxIntSkewDef, maxIntSkewWig));
            limBounSkew.insert(facei, min(maxBounSkewDef, maxBounSkewWig));
            limWeight.insert(facei, max(minWeightDef, minWeightWig));
            limRatio.insert(facei, max(minVolRatioDef, minVolRatioWig));
            limTwist.insert(facei, max(minTwistDef, minTwistWig));
            limOrtho.insert(facei, min(qBadDef, qBadWig));
        }
        else if (
            !isCleanFace 
            && 
            !isDirtyFace
        ) //must be a wiggle face
        {
            limPyrVol.insert(facei,   min(minPyrVolTight, min(pyrOwnOld[facei],pyrNeiOld[facei])));
            limConcave.insert(facei,  max(maxConcaveTight, concOld[facei]));
            limIntSkew.insert(facei,  max(maxIntSkewTight, skewOld[facei]));
            limBounSkew.insert(facei, max(maxBounSkewTight, skewOld[facei]));
            limWeight.insert(facei,   min(minWeightTight, weightOld[facei]));
            limRatio.insert(facei,    min(minVolRatioTight, ratioOld[facei]));
            limTwist.insert(facei,    min(minTwistTight, twistOld[facei]));
            limOrtho.insert(facei,    max(qBadTight, orthoOld[facei]));
        }
        else
        {
            // Face has violations but does not touch a wiggle point (e.g. from general laplacian smoothing)
            limPyrVol.insert(facei,   min(minPyrVolDef, min(pyrOwnOld[facei],pyrNeiOld[facei])));
            limConcave.insert(facei,  max(maxConcaveDef, concOld[facei]));
            limIntSkew.insert(facei,  max(maxIntSkewDef, skewOld[facei]));
            limBounSkew.insert(facei, max(maxBounSkewDef, skewOld[facei]));
            limWeight.insert(facei,   min(minWeightDef, weightOld[facei]));
            limRatio.insert(facei,    min(minVolRatioDef, ratioOld[facei]));
            limTwist.insert(facei,    min(minTwistDef, twistOld[facei]));
            limOrtho.insert(facei,    max(qBadDef, orthoOld[facei]));
        }
    }

    // 4. Local Lambda Setup
    scalarField lambda(mesh_.nPoints(), 1.0);
    const label maxIters = movementRelaxIters_; 
    const scalar lambdaReduction = movementRelaxation_;
    
    pointField pTest = p0;
    bool completelyValid = false;

    Info << "   -> Starting Local Guarded Morph on " 
         << returnReduce(activeFaces.size(), sumOp<label>()) << " active faces." << endl;

    // 5. Local Rollback Loop
    for (label iter = 0; iter < maxIters; ++iter)
    {
        syncTools::syncPointList(mesh_, lambda, minMagSqrEqOp<scalar>(), 1.0);

        forAll(pTest, ptI)
        {
            pTest[ptI] = p0[ptI] + (lambda[ptI] * disp[ptI]);
        }

        bool iterValid = true;
        
        // --- TELEMETRY COUNTERS ---
        label nFailVol = 0, nFailArea = 0, nFailQual = 0;
        label nPassAbs = 0, nPassImp = 0;

        label nFailOrtho = 0, nFailSkew = 0, nFailPyr = 0;
        label nFailConcave = 0, nFailWeight = 0, nFailRatio = 0, nFailTwist = 0;

        // --- TIER 1 EVALUATION (Cells) ---
        labelHashSet failedCells;
        labelHashSet failedAssociationCells;
        for (const label celli : activeCells)
        {
            scalar volTest; vector ccTest;
            calculateCellVolumeAndCentre(pTest, celli, volTest, ccTest);
            
            if (volTest <= VSMALL)
            {
                iterValid = false;
                nFailVol++;
                failedCells.insert(celli);
            }
        }

        // --- TIER 1, 2, and 3 EVALUATION (Faces) ---
        for (const label facei : activeFaces)
        {
            bool faceFails = false;

            label ownerCellI = mesh_.faceOwner()[facei];
            label neighbourCellI = mesh_.faceNeighbour()[facei];

            if (failedCells.found(ownerCellI) || failedCells.found(neighbourCellI))
            {
                faceFails = true;
            }

            // Tier 1 (Face Area)
            scalar areaTest = mag(mesh_.faces()[facei].areaNormal(pTest));
            if (areaTest <= VSMALL)
            {
                faceFails = true;
                nFailArea++;
            }
            else if (!faceFails)
            {
                // Calculate new dynamic metrics
                scalar orthoNew = calculateFaceNonOrtho(pTest, facei);
                scalar skewNew  = calculateFaceSkewness(pTest, facei);
                scalar pyrVolOwnNew = calculateFacePyramidVolume(pTest, facei, true);
                scalar concNew  = calculateFaceConcavity(pTest, facei);
                scalar weightNew = calculateInterpolationWeight(pTest, facei);
                scalar ratioNew = calculateVolumeRatio(pTest, facei);
                scalar twistNew = calculateFaceTwist(pTest, facei);
                
                scalar pyrVolNeiNew = GREAT;
                if (mesh_.isInternalFace(facei)) {
                    pyrVolNeiNew = calculateFacePyramidVolume(pTest, facei, false);
                }

                // Check against Per-Face dynamic limits
                scalar targetSkewLimit = mesh_.isInternalFace(facei) ? limIntSkew[facei] : limBounSkew[facei];

                bool orthoPassAbs  = (orthoNew <= limOrtho[facei]);
                bool skewPassAbs   = (skewNew <= targetSkewLimit);
                bool pyrOwnPassAbs = (pyrVolOwnNew > limPyrVol[facei]);
                bool pyrNeiPassAbs = (pyrVolNeiNew > limPyrVol[facei]);
                bool concPassAbs   = (concNew <= limConcave[facei]);
                bool weightPassAbs = (weightNew >= limWeight[facei]);
                bool ratioPassAbs  = (ratioNew >= limRatio[facei]);
                bool twistPassAbs  = (twistNew >= limTwist[facei]);

                bool allAbsPass = orthoPassAbs && skewPassAbs && pyrOwnPassAbs && 
                                  pyrNeiPassAbs && concPassAbs && weightPassAbs && 
                                  ratioPassAbs && twistPassAbs;

                bool faceFailsQual = false;

                if (!allAbsPass)
                {
                    // Strict Improvement Check: Must improve on ALL failed limits
                    bool allFailedImproved = true;

                    if (!orthoPassAbs  && (orthoNew >= orthoOld[facei]))       allFailedImproved = false;
                    if (!skewPassAbs   && (skewNew >= skewOld[facei]))         allFailedImproved = false;
                    if (!pyrOwnPassAbs && (pyrVolOwnNew <= pyrOwnOld[facei]))  allFailedImproved = false;
                    if (!pyrNeiPassAbs && (pyrVolNeiNew <= pyrNeiOld[facei]))  allFailedImproved = false;
                    if (!concPassAbs   && (concNew >= concOld[facei]))         allFailedImproved = false;
                    if (!weightPassAbs && (weightNew <= weightOld[facei]))     allFailedImproved = false;
                    if (!ratioPassAbs  && (ratioNew <= ratioOld[facei]))       allFailedImproved = false;
                    if (!twistPassAbs  && (twistNew <= twistOld[facei]))       allFailedImproved = false;

                    if (!allFailedImproved)
                    {
                        faceFailsQual = true;
                    }
                }

                // Tally specific failures (for the printout)
                if (!orthoPassAbs) nFailOrtho++;
                if (!skewPassAbs)  nFailSkew++;
                if (!pyrOwnPassAbs || !pyrNeiPassAbs) nFailPyr++;
                if (!concPassAbs)  nFailConcave++;
                if (!weightPassAbs)nFailWeight++;
                if (!ratioPassAbs) nFailRatio++;
                if (!twistPassAbs) nFailTwist++;

                if (faceFailsQual)
                {
                    faceFails = true;
                    nFailQual++;
                }
                else
                {
                    if (allAbsPass) nPassAbs++;
                    else nPassImp++;
                }
            }

            if (faceFails)
            {
                failedAssociationCells.insert(ownerCellI);
                failedAssociationCells.insert(neighbourCellI);
                
                /*
                iterValid = false;
                const face& f = mesh_.faces()[facei];
                forAll(f, i)
                {
                    label ptI = f[i];
                    if (activePointsSet.found(ptI))
                    {
                        lambda[ptI] *= lamSign * lambdaReduction;
                    }
                }*/
            }
        }

        labelHashSet splitPoints;

        for (const label facei : activeFaces) 
        {
            //now for all faces that somehow are part of an associated cell we roll back
            bool associated = false;

            label ownerCellI = mesh_.faceOwner()[facei];
            label neighbourCellI = mesh_.faceNeighbour()[facei];

            if(failedAssociationCells.found(ownerCellI)){
                associated = true;
            }

            if(failedAssociationCells.found(neighbourCellI)){
                associated = true;
            }

            if(
                associated
            ){
                const face& f = mesh_.faces()[facei];
                forAll(f, i)
                {
                    label ptI = f[i];
                    if (activePointsSet.found(ptI))
                    {
                        splitPoints.insert(ptI);
                    }
                }
            }
        }

        for (const label pointI : splitPoints)
        {
            lambda[pointI] *= lamSign * lambdaReduction;
        } 

        lamFactor *= lambdaReduction;
        if(lamFactor < limit)
        {
            if(pendulumState == 0){
                lamSign = -1;
                pendulumState = 1;
                Info << "      -> Trying displacement flip" << endl;
            }else if(pendulumState == 1){
                pendulumState = 2;
                Info << "      -> Reverting displacement flip" << endl;
            }else if(pendulumState == 3){
                lamSign = 1;
                pendulumState = 4;
                Info << "      -> Back to normal operation" << endl;
            }
        }
        
        // Parallel Reductions
        reduce(nFailVol, sumOp<label>());
        reduce(nFailArea, sumOp<label>());
        reduce(nFailQual, sumOp<label>());
        reduce(nPassAbs, sumOp<label>());
        reduce(nPassImp, sumOp<label>());
        
        reduce(nFailOrtho, sumOp<label>());
        reduce(nFailSkew, sumOp<label>());
        reduce(nFailPyr, sumOp<label>());
        reduce(nFailConcave, sumOp<label>());
        reduce(nFailWeight, sumOp<label>());
        reduce(nFailRatio, sumOp<label>());
        reduce(nFailTwist, sumOp<label>());

        // Verbose Printout on Iteration 0 (Unscaled Displacement)
        if (iter == 0 && (nFailQual > 0 || nFailVol > 0 || nFailArea > 0))
        {
            Info << "      Initial displacement causes the following errors (Dynamic Default/Strict Limits):" << nl
                 << "         non-orthogonality (> limit)                            : " << nFailOrtho << nl
                 << "         faces with face pyramid volume (< limit)               : " << nFailPyr << nl
                 << "         faces with concavity (> limit)                         : " << nFailConcave << nl
                 << "         faces with skewness (> limit)                          : " << nFailSkew << nl
                 << "         faces with interpolation weights (< limit)             : " << nFailWeight << nl
                 << "         faces with volume ratio of neighbour cells (< limit)   : " << nFailRatio << nl
                 << "         faces with face twist (< limit)                        : " << nFailTwist << nl
                 << "         [Critical] Cell Volume <= 0                            : " << nFailVol << nl
                 << "         [Critical] Face Area <= 0                              : " << nFailArea << endl;
        }

        // Compact iteration tracker
        Info << "      Iter " << iter << " | "
             << "Pass(Abs/Imp): " << nPassAbs << "/" << nPassImp << " || "
             << "Fail[C:" << nFailVol << " A:" << nFailArea << " | "
             << "O:" << nFailOrtho << " S:" << nFailSkew 
             << " P:" << nFailPyr << " C:" << nFailConcave 
             << " W:" << nFailWeight << " R:" << nFailRatio 
             << " T:" << nFailTwist << "]" << endl;

        // Break if clean
        if ((nFailVol + nFailArea + nFailQual) == 0)
        {
            completelyValid = true;
            break;
        }
    }
    
    Info << "Relaxation done." << endl;

    // 6. Final Sync and Commit
    syncTools::syncPointList(mesh_, lambda, minMagSqrEqOp<scalar>(), 1.0);

    pointVectorField adaptedDisplacement(finalPointDisp);
    forAll(pTest, ptI)
    {
        pTest[ptI] = p0[ptI] + (lambda[ptI] * disp[ptI]);
        adaptedDisplacement[ptI] = lambda[ptI] * disp[ptI];
    }

    syncTools::syncPointList(mesh_, pTest, minEqOp<point>(), point(GREAT, GREAT, GREAT));

    // NaN / Inf Check
    forAll(pTest, ptI)
    {
        if (!finite(mag(pTest[ptI])))
        {
            FatalErrorInFunction
                << "NaN or Inf detected in pTest at point " << ptI 
                << " coord: " << pTest[ptI] 
                << abort(FatalError);
        }
    }

    Info << "Attempting to move mesh with " << returnReduce(pTest.size(), sumOp<label>()) << " points." << endl;

    // Logging the Lambda Range
    scalar minActiveLambda = 1.0;
    scalar maxActiveLambda = -1.0;
    for (const label ptI : activePointsSet)
    {
        minActiveLambda = min(minActiveLambda, lambda[ptI]);
        maxActiveLambda = max(maxActiveLambda, lambda[ptI]);
    }
    
    reduce(minActiveLambda, minOp<scalar>());
    reduce(maxActiveLambda, maxOp<scalar>());

    if (completelyValid)
    {
        Info << "   -> Guarded displacement achieved. Active lambda range: [" 
             << minActiveLambda << " to " << maxActiveLambda << "]" << endl;
    }
    else
    {
        Info << "   -> Notice: Max rollback iterations hit. Mesh safety preserved via severe local damping. Active lambda range: [" 
             << minActiveLambda << " to " << maxActiveLambda << "]" << endl;
    }

    // Move the mesh
    if(move){
        mesh_.movePoints(pTest);
    }
    
    if(wiggleInitialized){
        Info << "Transfering wiggle lambda with strict limit check..." << endl;
        
        DynamicList<scalar> newLambda(wigglePoints.size());
        labelList badPointIDs = wigglePoints.toc();
        
        forAll(badPointIDs, ind)
        {
            label pointI = badPointIDs[ind];
            
            // Check if this specific point was originally "clean" 
            // by checking if it touches any face that violates the strict limits
            bool isCleanPoint = true;
            const labelList& pFaces = mesh_.pointFaces()[pointI];
            
            for (const label facei : pFaces)
            {
                // Re-evaluate the metric violation for this face using the strict wiggle limits
                // Note: We use the final displaced mesh (pTest)
                scalar ortho = calculateFaceNonOrtho(pTest, facei);
                scalar skew  = calculateFaceSkewness(pTest, facei);
                
                // If it violates ANY strict wiggle limit, it's NOT a clean point
                if (ortho > qBadWig || skew > (mesh_.isInternalFace(facei) ? maxIntSkewWig : maxBounSkewWig))
                {
                    isCleanPoint = false;
                    break;
                }
            }

            // If it is clean, enforce lambda 0, otherwise store the applied lambda
            scalar localLamb = isCleanPoint ? 0.0 : lambda[pointI];
            newLambda.append(localLamb);
        }
        wiggleLambda = newLambda;
    }

    return adaptedDisplacement;
}

void meshFixer::solveLocalLaplacian(const pointMesh& pMesh)
{
    scalarField nonOrtho = calculateFaceNonOrtho();
    scalarField skewness = calculateFaceSkewness();
    boolList isBadFace(mesh_.nFaces(), false);
    boolList isBadCell(mesh_.nCells(), false);

    for (label f = 0; f < mesh_.nInternalFaces(); ++f) {
        if (nonOrtho[f] > orthoSpring_ || skewness[f] > skewSpring_) {
            isBadFace[f] = true;
            isBadCell[mesh_.faceOwner()[f]] = true;
            isBadCell[mesh_.faceNeighbour()[f]] = true;
        }
    }

    boolList bufferCells = isBadCell;
    for (label layer = 0; layer < nSmoothLayers_; ++layer) {
        boolList nextBuffer = bufferCells;
        forAll(mesh_.cells(), c) {
            if (bufferCells[c]) {
                for (label f : mesh_.cells()[c]) {
                    if (mesh_.isInternalFace(f)) {
                        nextBuffer[mesh_.faceOwner()[f]] = true;
                        nextBuffer[mesh_.faceNeighbour()[f]] = true;
                    }
                }
            }
        }
        bufferCells = nextBuffer;
    }

    boolList activePoints(mesh_.nPoints(), false);
    forAll(activePoints, ptI) {
        const labelList& pCells = mesh_.pointCells()[ptI];
        bool allInBuffer = true;
        for (label c : pCells) { if (!bufferCells[c]) allInBuffer = false; }
        if (allInBuffer && !isBoundaryPoint_[ptI]) activePoints[ptI] = true;
    }

    pointVectorField laplaceDisp
    (
        IOobject("laplaceDisp", mesh_.time().timeName(), mesh_, IOobject::NO_READ, IOobject::NO_WRITE),
        pMesh, dimensionedVector(dimLength, vector::zero)
    );
    
    for (label iter = 0; iter < laplacianIters_; ++iter)
    {
        vectorField iterNetForce(mesh_.nPoints(), vector::zero);
        scalarField iterPointsAffecting(mesh_.nPoints(), 0.0);
        
        const pointField& pts = mesh_.points();

        for (label ptI = 0; ptI < mesh_.nPoints(); ++ptI) {
            if (!activePoints[ptI]) continue;

            vector localNetForce = vector::zero;
            scalar localPointsAffecting = 0.0;
            bool isOnBadFace = isBadPoint_[ptI]; 

            point currentPos = pts[ptI] + laplaceDisp[ptI];

            const labelList& pEdges = mesh_.pointEdges()[ptI];
            for (label eI : pEdges) {
                const edge& e = mesh_.edges()[eI];
                label nPt = (e.start() == ptI) ? e.end() : e.start();

                point neighborPos = pts[nPt] + laplaceDisp[nPt];
                localNetForce += (neighborPos - currentPos);
                localPointsAffecting += 1.0;
            }

            if (isOnBadFace && false) {
                const labelList& pFaces = mesh_.pointFaces()[ptI];
                scalar strengthMultiplier = 1.0;

                for (label f : pFaces) {
                    if (isBadFace[f] && mesh_.isInternalFace(f)) {
                        point fCenter = 0.5 * (mesh_.C()[mesh_.faceOwner()[f]] + mesh_.C()[mesh_.faceNeighbour()[f]]);
                        localNetForce += strengthMultiplier * (fCenter - currentPos);
                        localPointsAffecting += 1.0;
                    }
                }
            }

            iterNetForce[ptI] = localNetForce;
            iterPointsAffecting[ptI] = localPointsAffecting;
        }

        syncTools::syncPointList(mesh_, iterNetForce, plusEqOp<vector>(), vector::zero);
        syncTools::syncPointList(mesh_, iterPointsAffecting, plusEqOp<scalar>(), 0.0);
        
        vectorField currentIterDisp(mesh_.nPoints(), vector::zero);

        for (label ptI = 0; ptI < mesh_.nPoints(); ++ptI) {
            if (!activePoints[ptI]) continue;

            vector netForce = iterNetForce[ptI];
            scalar pointsAffecting = iterPointsAffecting[ptI];
            
            if (pointsAffecting > VSMALL && mag(netForce) > VSMALL) {
                vector dir = netForce / mag(netForce);
                
                scalar moveMax = getMinEdgeLength(ptI) * 2;
                scalar moveForce = mag(netForce) / pointsAffecting;
                scalar moveRes = localLaplacianRelaxation_ * min(moveForce, moveMax);
                
                currentIterDisp[ptI] = dir * moveRes;
            }
        }

        laplaceDisp.primitiveFieldRef() += currentIterDisp;
    }

    Info<< "   Applying iterative diffusion (" << laplacianIters_ << " iterations)..." << endl;
    
    // Fixed the truncated call here:
    //applyRelativeGuardedMorph(laplaceDisp);
    applyQualityGuardedMorph(applyRelativeGuardedMorph(laplaceDisp,false));
}

void meshFixer::fix()
{
    Info<< "\nStarting meshFixer over " << totalIterations_ << " iterations." << endl;
    identifyBoundaries();

    const pointMesh& pMesh = pointMesh::New(mesh_);

    haveBadFaces_ = true;

    for (label iter = 0; iter < totalIterations_; ++iter)
    {
        if (writeIntermediate_)
        {
            runTime_++;
        }

        Info<< "Iteration " << iter + 1 << " / " << totalIterations_ 
            << " (Time: " << runTime_.timeName() << ")" << endl;

        pointVectorField localDisplacement
        (
            IOobject("localDisp", mesh_.time().timeName(), mesh_, IOobject::NO_READ, IOobject::NO_WRITE),
            pMesh,
            dimensionedVector(dimLength, vector::zero)
        );

        if(discreteStep_)
        {
            calcLocalCorrections(localDisplacement);
            solveGlobalDiffusion(localDisplacement); 
        }
        if(laplacianStep_)
        {
            solveLocalLaplacian(pMesh);
        }
        if(wiggleStep_)
        {
            wiggleSmooth(pMesh);
        }





        if (writeIntermediate_)
        {
            mesh_.write();
        }

        if(
            iter % writeEvery_ == 0
            &&
            iter != 0
        )
        {
            runTime_++;
            mesh_.write();
        }

        if(haveBadFaces_ == false){
            Info << "No bad faces remaining. Exiting smoothing loop early at iteration " << iter + 1 << "." << endl;
            break;
        }
    }
    
    Info<< "Fix sequence complete. Writing final mesh." << endl;
    if (!writeIntermediate_) 
    {
        runTime_++;
        mesh_.write(); 
    }
}

void meshFixer::wiggleSmooth(const pointMesh& pMesh)
{
    // --- 1. Helper Code: Read Limits from relativeQualityControlDict_ ---
    const scalar maxNonOrtho = wiggleQualityControlDict_.lookupOrDefault<scalar>("maxNonOrtho", 70);
    const scalar minPyrVol  = wiggleQualityControlDict_.lookupOrDefault<scalar>("minVol", 1e-13);
    const scalar maxConcave = wiggleQualityControlDict_.lookupOrDefault<scalar>("maxConcave", 80.0);
    const scalar maxIntSkew = wiggleQualityControlDict_.lookupOrDefault<scalar>("maxInternalSkewness", 4.0);
    // Boundary skew not strictly needed here since loop only covers internal faces, but good for completeness
    const scalar maxBounSkew = wiggleQualityControlDict_.lookupOrDefault<scalar>("maxBoundarySkewness", 20.0);
    const scalar minWeight   = wiggleQualityControlDict_.lookupOrDefault<scalar>("minFaceWeight", 0.0001);
    const scalar minVolRatio = wiggleQualityControlDict_.lookupOrDefault<scalar>("minVolRatio", 0.01);
    const scalar minTwist    = wiggleQualityControlDict_.lookupOrDefault<scalar>("minTwist", 0.05);

    // --- 1. Read Limits from relativeQualityControlDict_ ---
    const scalar maxNonOrthoLimit   = relativeQualityControlDict_.lookupOrDefault<scalar>("maxNonOrtho", 70.0);
    const scalar minPyrVolLimit     = relativeQualityControlDict_.lookupOrDefault<scalar>("minVol", 1e-13);
    const scalar maxConcaveLimit    = relativeQualityControlDict_.lookupOrDefault<scalar>("maxConcave", 80.0);
    const scalar maxIntSkewLimit    = relativeQualityControlDict_.lookupOrDefault<scalar>("maxInternalSkewness", 4.0);
    const scalar maxBounSkewLimit   = relativeQualityControlDict_.lookupOrDefault<scalar>("maxBoundarySkewness", 20.0);
    const scalar minWeightLimit     = relativeQualityControlDict_.lookupOrDefault<scalar>("minFaceWeight", 0.0001);
    const scalar minVolRatioLimit   = relativeQualityControlDict_.lookupOrDefault<scalar>("minVolRatio", 0.01);
    const scalar minTwistLimit      = relativeQualityControlDict_.lookupOrDefault<scalar>("minTwist", 0.05);


    // Current points of the mesh
    const pointField& p = mesh_.points();

    // --- 2. Initialize Boundary Mask ---
    if(!wiggleInitialized)
    {
        wiggleInitialized = true;
        for (label facei = mesh_.nInternalFaces(); facei < mesh_.nFaces(); ++facei)
        {
            const face& f = mesh_.faces()[facei];
            forAll(f, i)
            {
                label ptI = f[i];
                boundaryPoints.insert(ptI);
            }
        }  
    }
    
    // --- 3. Identify Bad Points based on Limit Quality Metrics ---
    labelHashSet currentBadPoints;
    labelHashSet currentBranchPoints;

    for (label facei = 0; facei < mesh_.nInternalFaces(); ++facei)
    {
        // Calculate all metrics dynamically
        scalar ortho    = calculateFaceNonOrtho(p, facei);
        scalar skew     = calculateFaceSkewness(p, facei);
        scalar pyrOwn   = calculateFacePyramidVolume(p, facei, true);
        scalar pyrNei   = calculateFacePyramidVolume(p, facei, false);
        scalar conc     = calculateFaceConcavity(p, facei);
        scalar weight   = calculateInterpolationWeight(p, facei);
        scalar ratio    = calculateVolumeRatio(p, facei);
        scalar twist    = calculateFaceTwist(p, facei);

        // A face is bad if it violates ANY of the limits
        // A face is bad if it violates ANY of the defined limits
        bool badFace = (ortho > maxNonOrthoLimit) ||
                       (skew > maxIntSkewLimit) ||
                       (pyrOwn <= minPyrVolLimit) ||
                       (pyrNei <= minPyrVolLimit) ||
                       (conc > maxConcaveLimit) ||
                       (weight < minWeightLimit) ||
                       (ratio < minVolRatioLimit) ||
                       (twist < minTwistLimit);

        if(badFace)
        {
            const face& f = mesh_.faces()[facei];
            forAll(f, i)
            {
                label ptI = f[i];
                if(!boundaryPoints.found(ptI))
                {
                    currentBadPoints.insert(ptI);
                }
            }
        }
    }

    // --- 3. Identify Bad Points based on Wiggle Quality Metrics ---
    for (label facei = 0; facei < mesh_.nInternalFaces(); ++facei)
    {
        // Calculate all metrics dynamically
        scalar ortho    = calculateFaceNonOrtho(p, facei);
        scalar skew     = calculateFaceSkewness(p, facei);
        scalar pyrOwn   = calculateFacePyramidVolume(p, facei, true);
        scalar pyrNei   = calculateFacePyramidVolume(p, facei, false);
        scalar conc     = calculateFaceConcavity(p, facei);
        scalar weight   = calculateInterpolationWeight(p, facei);
        scalar ratio    = calculateVolumeRatio(p, facei);
        scalar twist    = calculateFaceTwist(p, facei);

        // A face is bad if it violates ANY of the limits
        bool badFace = (ortho > maxNonOrtho) ||
                       (skew > maxIntSkew) ||
                       (pyrOwn <= minPyrVol) ||
                       (pyrNei <= minPyrVol) ||
                       (conc > maxConcave) ||
                       (weight < minWeight) ||
                       (ratio < minVolRatio) ||
                       (twist < minTwist);

        if(badFace)
        {
            const face& f = mesh_.faces()[facei];
            forAll(f, i)
            {
                label ptI = f[i];
                if(!boundaryPoints.found(ptI))
                {
                    if(!currentBadPoints.found(ptI))
                    {
                        currentBadPoints.insert(ptI);
                        currentBranchPoints.erase(ptI);
                        //we were bad enough ourselfs to wiggle so were no longer a branch
                    }
                    
                }
            }
        }
    }

    // --- 3b. Expand Bad Points Set by surroundingLayers ---
    for (label layer = 0; layer < wiggleSurroundingLayers_; ++layer)
    {
        // Extract current points to an array so we don't infinitely loop as we insert new ones
        labelList currentPointsArray = currentBadPoints.toc();
        
        forAll(currentPointsArray, i)
        {
            label ptI = currentPointsArray[i];
            
            // Get all edges connected to this point
            const labelList& pEdges = mesh_.pointEdges()[ptI];
            
            forAll(pEdges, edgeIdx)
            {
                const edge& e = mesh_.edges()[pEdges[edgeIdx]];
                
                // Get the point on the other side of the edge
                label adjacentPt = (e.start() == ptI) ? e.end() : e.start();

                // Add to the wiggle set if it's not a boundary point
                if (!boundaryPoints.found(adjacentPt))
                {
                    if(!currentBadPoints.found(adjacentPt))
                    {
                        currentBadPoints.insert(adjacentPt);
                        currentBranchPoints.insert(adjacentPt);
                    }
                    
                }
            }
        }
    }

    labelList badPointIDs = currentBadPoints.toc();
    
    // Perform MPI reductions to sum across all processor cores
    label globalBadPoints = returnReduce(currentBadPoints.size(), sumOp<label>());
    label globalBranchPoints = returnReduce(currentBranchPoints.size(), sumOp<label>());
    
    // Calculate the original "core" points that were bad on their own merit
    label globalCorePoints = globalBadPoints - globalBranchPoints;

    Info << "Wiggeling Points (including " << wiggleSurroundingLayers_ << " surrounding layers): " 
         << globalBadPoints << nl
         << "    -> Core (Failing) Points : " << globalCorePoints << nl
         << "    -> Branch (Buffer) Points: " << globalBranchPoints << endl;

    if(globalBadPoints == 0)
    {
        haveBadFaces_ = false;
    }else{
        haveBadFaces_ = true;
    }

    scalarList minEdge(
        currentBadPoints.size()
    );

    {
        label i = 0;
        forAll(badPointIDs,ind)
        {
            label pointI = badPointIDs[ind];
            minEdge[i] = getMinEdgeLength(pointI);
            ++i;
        };
    }

    //size all the needed containers
    wiggleDisplacement.reserve(
        currentBadPoints.size()
    );
    DynamicList<scalar> thetaList_(currentBadPoints.size());
    DynamicList<scalar> phiList_(currentBadPoints.size());
    Map<label> pointToId;
    DynamicList<vector> nextDisplacement;
    nextDisplacement.reserve(currentBadPoints.size());
            

    label i = 0;
    forAll(badPointIDs,ind)
    {
        label pointI = badPointIDs[ind];
        pointToId.insert(pointI,i);
     //first we check what information we have for this point
     //1. know nothing -> random step
     //2. know only 1 last step -> random step
     //3. know more but lambda to bac -> random step
     //4. know 2 and good -> gradient step
        scalar theta_ = 0;
        scalar phi_ = 0;

        bool gradientStep = false;
        label containedIn = 0;

        if(wigglePoints.found(pointI))
        {
            if(verboseWiggle_){
                Info << "Found in 0"<< endl;
            }
            containedIn++;
        }
        if(wigglePointsM1.found(pointI))
        {
            if(verboseWiggle_){
                Info << "Found in 1"<< endl;
            }
            containedIn++;
        }

        if(containedIn == 2)
        {
            scalar lambda_ = wiggleLambda[
                pointToIdLabel[pointI]
            ];
            if(verboseWiggle_){
                Info << "Got 2 infomrations lambda was: "<< lambda_ << endl;
            }
            if(lambda_ > randomLambLimit){
                gradientStep = true;
            }
        }
        if(!gradientStep){
            if(verboseWiggle_){
                Info << "Performing Random Step for Point: "<< pointI << endl; 
            }
            theta_ = rndGen.sample01<scalar>()*2*Foam::constant::mathematical::pi;
            phi_ = rndGen.sample01<scalar>()*2*Foam::constant::mathematical::pi;
        }else{
            if(verboseWiggle_){
                Info << "Performing gradient based step for Point: " << pointI << endl; 
            }
            
            //Pull the last thetas and lambdas
            label id0 = pointToIdLabel[pointI];
            label idM1 = pointToIdLabelM1[pointI];


            scalar lambda0 = wiggleLambda[id0];
            if (complexReward_)
            {
                scalar adjustedSelfFactor_ = selfFactor_;
                //check if we were a branch point, then we dont reward our own movement
                if(branchePoints.found(pointI))
                {
                    Info << "Branch Point, set own weight to 0" << endl;
                    adjustedSelfFactor_ = 0;
                }else{
                    Info << "Native Point" << endl;
                }
                scalar sumNeighborLambda = 0.0;
                label count = 0;

                const labelList& pCells = mesh_.pointCells()[pointI];
                for (label celli : pCells)
                {
                    const labelList& cPts = mesh_.cellPoints()[celli];
                    for (label neighborPt : cPts)
                    {
                        // Check if neighbor is a wiggling point
                        if (neighborPt != pointI && wigglePoints.found(neighborPt))
                        {
                            // Get lambda from last iteration
                            sumNeighborLambda += wiggleLambda[pointToIdLabel[neighborPt]];
                            count++;
                        }
                    }
                }

                if (count > 0)
                {

                    //Info << "Complex Lambda, self value: " << lambda0 << endl;
                    lambda0 = (adjustedSelfFactor_ * lambda0) + sumNeighborLambda;
                    lambda0 = lambda0 / (adjustedSelfFactor_ + scalar(count));
                    ///Info << "Complex Lambda, complex value: " << lambda0 << endl;
                }
            }
            scalar lambdaM1  = wiggleLambdaM1[idM1];

            scalar theta0 = thetaWiggle[id0];
            scalar thetaM1 = thetaWiggleM1[idM1];

            scalar phi0 = phiWiggle[id0];
            scalar phiM1 = phiWiggleM1[idM1];
            
            scalar deltaLambda = lambda0 - lambdaM1;
            scalar deltaTheta = calculateAngleDelta(theta0,thetaM1);
            scalar deltaPhi = calculateAngleDelta(phi0,phiM1);
            //prevent division by 0
            deltaTheta = Foam::sign(deltaTheta)*max(1e-4,mag(deltaTheta));
            deltaPhi = Foam::sign(deltaPhi)*max(1e-4,mag(deltaPhi));
            
            scalar gradTheta = deltaLambda/deltaTheta;
            scalar gradPhi = deltaLambda/deltaPhi;

            theta_ = theta0 + (1-mag(lambda0))*Foam::sign(gradTheta)*gradientDelta*rndGen.sample01<scalar>();
            phi_ = phi0 + (1-mag(lambda0))*Foam::sign(gradPhi)*gradientDelta*rndGen.sample01<scalar>();
            
            theta_ = normalizeAngle(theta_);
            phi_ = normalizeAngle(phi_);

        }

        vector direction(
            Foam::sin(theta_) * Foam::cos(phi_), // x-component
            Foam::sin(theta_) * Foam::sin(phi_), // y-component
            Foam::cos(theta_)                    // z-component
        );

        direction = direction*wiggleFactor*minEdge[i];

        nextDisplacement.append(direction);
        thetaList_.append(theta_);
        phiList_.append(phi_);
        ++i;

        
    }
    //transfer one step back
    pointToIdLabelM1 = pointToIdLabel;
    wigglePointsM1 = wigglePoints;
    wiggleLambdaM1 = wiggleLambda;
    thetaWiggleM1 = thetaWiggle;
    phiWiggleM1 = phiWiggle;

    //transfer new Data
    pointToIdLabel = pointToId;
    wigglePoints = currentBadPoints;
    thetaWiggle = thetaList_;
    phiWiggle = phiList_;
    wiggleDisplacement = nextDisplacement;
    branchePoints = currentBranchPoints;

    //at this step the new displacements musst be stored in wiggleDisplacement
    //as well as pointToIdLable and wigglePoints must be populated
    //calculating the displacement field
    pointVectorField wiggleDisp
    (
        IOobject("wiggleDisp", mesh_.time().timeName(), mesh_, IOobject::NO_READ, IOobject::NO_WRITE),
        pMesh, dimensionedVector(dimLength, vector::zero)
    );

    forAll(badPointIDs,ind)
    {
        label pointI = badPointIDs[ind];
        label index = pointToIdLabel[pointI];
        wiggleDisp[pointI] = wiggleDisplacement[index];
    }

    //sync the wiggle
    syncTools::syncPointList(mesh_, wiggleDisp.primitiveFieldRef(), minMagSqrOp<vector>(), point(GREAT, GREAT, GREAT));

    Info << "Applying wiggle diffusion" << endl;
    applyQualityGuardedMorph(applyRelativeGuardedMorph(wiggleDisp,false));
}

scalar meshFixer::calculateAngleDelta(
    const scalar ang1, const scalar ang2
){
    using namespace Foam::constant::mathematical; 
    
    scalar diff = ang1 - ang2;
    
    // 3*pi is the equivalent of +540 offset. 2*pi is 360. pi is 180.
    return std::fmod((diff + 3.0 * pi), 2.0 * pi) - pi;
}

scalar meshFixer::normalizeAngle(
    const scalar ang
){
    using namespace Foam::constant::mathematical;
    
    scalar result = std::fmod(ang, 2.0 * pi);
    return result + (2.0 * pi * Foam::neg(result));
}


// ------------------------------------------------------------------------- //
// Helper: Calculate Face Pyramid Volume
// ------------------------------------------------------------------------- //
scalar meshFixer::calculateFacePyramidVolume
(
    const pointField& p,
    const label facei,
    const bool isOwner
) const
{
    const face& f = mesh_.faces()[facei];
    vector fc = f.centre(p);
    vector fa = f.areaNormal(p); // Vector pointing from owner to neighbour

    label celli = isOwner ? mesh_.faceOwner()[facei] : mesh_.faceNeighbour()[facei];
    
    scalar vol; 
    vector cc;
    calculateCellVolumeAndCentre(p, celli, vol, cc);

    // Calculate actual volume: 1/3 * (Area . height)
    scalar pyrVol = (1.0/3.0) * (fa & (fc - cc));

    // If evaluating the neighbour, the area normal points inward, so we invert it
    return isOwner ? pyrVol : -pyrVol;
}

// ------------------------------------------------------------------------- //
// Helper: Calculate Face Concavity (Max Angle)
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// Helper: Calculate Face Concavity (Matched OF Float Tolerance)
// ------------------------------------------------------------------------- //
scalar meshFixer::calculateFaceConcavity
(
    const pointField& p,
    const label facei
) const
{
    const face& f = mesh_.faces()[facei];
    
    if (f.size() < 4) 
    {
        return 0.0; // Triangles cannot be concave
    }

    vector faceNormal = f.areaNormal(p);
    faceNormal /= mag(faceNormal) + VSMALL;

    scalar maxConcaveSin = 0.0;

    vector ePrev = p[f.first()] - p[f.last()];
    ePrev /= mag(ePrev) + VSMALL;

    forAll(f, fp)
    {
        vector eNext = p[f.nextLabel(fp)] - p[f[fp]];
        eNext /= mag(eNext) + VSMALL;

        vector edgeNormal = ePrev ^ eNext;
        scalar magEdgeNormal = mag(edgeNormal);

        if (magEdgeNormal > SMALL)
        {
            vector nEdgeNormal = edgeNormal / magEdgeNormal;
            
            // EXACT MATCH: OpenFOAM uses < SMALL (1e-15), not < -SMALL
            if ((nEdgeNormal & faceNormal) < SMALL)
            {
                maxConcaveSin = max(maxConcaveSin, magEdgeNormal);
            }
        }
        ePrev = eNext;
    }

    return radToDeg(asin(min(1.0, maxConcaveSin)));
}

// ------------------------------------------------------------------------- //
// Helper: Calculate Interpolation Weights
// ------------------------------------------------------------------------- //
scalar meshFixer::calculateInterpolationWeight
(
    const pointField& p,
    const label facei
) const
{
    const face& f = mesh_.faces()[facei];
    vector fa = f.areaNormal(p);
    vector fc = f.centre(p);

    label own = mesh_.faceOwner()[facei];
    scalar ownVol; 
    vector ownCc;
    calculateCellVolumeAndCentre(p, own, ownVol, ownCc);

    scalar dOwn = mag(fa & (fc - ownCc));
    scalar dNei = 0.0;

    if (mesh_.isInternalFace(facei))
    {
        label nei = mesh_.faceNeighbour()[facei];
        scalar neiVol; 
        vector neiCc;
        calculateCellVolumeAndCentre(p, nei, neiVol, neiCc);
        
        dNei = mag(fa & (neiCc - fc));
    }
    else
    {
        // For uncoupled boundary faces, weight is implicitly 0.5 in standard FV.
        // For coupled (MPI), a true check requires syncing neighbour cell centers.
        // For local smoothing, substituting dOwn for dNei prevents false positives.
        dNei = dOwn; 
    }

    return min(dNei, dOwn) / (dNei + dOwn + VSMALL);
}

// ------------------------------------------------------------------------- //
// Helper: Calculate Volume Ratio
// ------------------------------------------------------------------------- //
scalar meshFixer::calculateVolumeRatio
(
    const pointField& p,
    const label facei
) const
{
    if (!mesh_.isInternalFace(facei))
    {
        return 1.0; // Cannot easily calculate ratio for boundaries without MPI swap
    }

    label own = mesh_.faceOwner()[facei];
    label nei = mesh_.faceNeighbour()[facei];

    scalar ownVol; vector ownCc;
    calculateCellVolumeAndCentre(p, own, ownVol, ownCc);

    scalar neiVol; vector neiCc;
    calculateCellVolumeAndCentre(p, nei, neiVol, neiCc);

    ownVol = mag(ownVol);
    neiVol = mag(neiVol);

    return min(ownVol, neiVol) / (max(ownVol, neiVol) + VSMALL);
}

// ------------------------------------------------------------------------- //
// Helper: Calculate Face Twist
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// Helper: Calculate Face Twist (OpenFOAM Compliant)
// ------------------------------------------------------------------------- //
scalar meshFixer::calculateFaceTwist
(
    const pointField& p,
    const label facei
) const
{
    const face& f = mesh_.faces()[facei];
    
    if (f.size() <= 3) 
    {
        return 1.0; // Triangles are always planar
    }

    vector fc = f.centre(p);
    
    // 1. Calculate the reference normal (nf) 
    // OpenFOAM evaluates twist against the cell-to-cell vector, NOT the face normal.
    vector nf(Zero);
    
    label own = mesh_.faceOwner()[facei];
    scalar ownVol; vector ownCc;
    calculateCellVolumeAndCentre(p, own, ownVol, ownCc);

    if (mesh_.isInternalFace(facei))
    {
        label nei = mesh_.faceNeighbour()[facei];
        scalar neiVol; vector neiCc;
        calculateCellVolumeAndCentre(p, nei, neiVol, neiCc);
        
        nf = neiCc - ownCc; // Internal: Cell-to-Cell vector
    }
    else
    {
        // Boundary: Face center to Owner cell center
        nf = fc - ownCc;
    }

    nf /= mag(nf) + VSMALL;

    scalar minTwist = GREAT;

    // 2. Check the normal of each face-center decomposed triangle against nf
    forAll(f, fp)
    {
        vector triArea = triPointRef::areaNormal(p[f[fp]], p[f.nextLabel(fp)], fc);
        scalar magTri = mag(triArea);

        if (magTri > VSMALL)
        {
            // (nf & triArea/magTri) is the exact twist metric used in checkFaceTwist
            scalar twist = (nf & triArea) / magTri;
            minTwist = min(minTwist, twist);
        }
    }

    return minTwist;
}

scalarField meshFixer::calculateFaceSkewness() const
{
    scalarField skewness(mesh_.nFaces(), 0.0);
    const surfaceVectorField& Sf = mesh_.Sf();
    const volVectorField& C = mesh_.C();
    const pointField& faceCentres = mesh_.faceCentres();

    for (label facei = 0; facei < mesh_.nInternalFaces(); ++facei)
    {
        label own = mesh_.faceOwner()[facei];
        label nei = mesh_.faceNeighbour()[facei];
        
        vector d = C[nei] - C[own];
        vector n = Sf[facei] / (mag(Sf[facei]) + VSMALL);
        
        scalar t = ((faceCentres[facei] - C[own]) & n) / ((d & n) + VSMALL);
        point intersection = C[own] + t * d;
        
        skewness[facei] = mag(faceCentres[facei] - intersection) / (mag(d) + VSMALL);
    }
    return skewness;
}

} // End namespace Foam