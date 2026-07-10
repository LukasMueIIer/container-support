# Project Structure

```text
.
├── Make
│   ├── files
│   └── options
├── meshFixerApp.C
├── meshFixer.C
├── meshFixer.H
├── meshMetrics.H
└── wrapCode.sh

1 directory, 7 files
```

---

# File Contents
### File: files
**Location:** `Make/files`

```
meshFixer.C
meshFixerApp.C

EXE = $(FOAM_USER_APPBIN)/meshFixer```

---

### File: options
**Location:** `Make/options`

```
EXE_INC = \
    -I$(LIB_SRC)/finiteVolume/lnInclude \
    -I$(LIB_SRC)/meshTools/lnInclude \
    -I$(LIB_SRC)/dynamicMesh/lnInclude 

EXE_LIBS = \
    -lfiniteVolume \
    -lmeshTools \
    -ldynamicMesh ```

---

### File: meshFixerApp.C
**Location:** `meshFixerApp.C`

```
#include "fvCFD.H"
#include "meshFixer.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Foam::meshFixer fixer(mesh);
    fixer.fix();

    return 0;
}```

---

### File: meshFixer.C
**Location:** `meshFixer.C`

```
#include "meshFixer.H"
#include "fvm.H"
#include "fvc.H"
#include "syncTools.H"
#include "polyMesh.H"
#include "volPointInterpolation.H"
#include "mathematicalConstants.H"
#include "unitConversion.H"

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
    writeIntermediate_ = smootherDict.lookupOrDefault<bool>("writeIntermediate", true);  
    writeDisplacement_ = smootherDict.lookupOrDefault<bool>("writeDisplacement", false); 
    localRelaxation_   = smootherDict.lookupOrDefault<scalar>("localRelaxation", 1.0);

    qualityControlDict_ = dict.subDict("meshQualityControls");
    movementRelaxIters_ = qualityControlDict_.lookupOrDefault<label>("nSmoothScale", 1000);
    movementRelaxation_ = qualityControlDict_.lookupOrDefault<scalar>("errorReduction", 0.5);

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

    wiggleFactor = 0.1;
    wiggleInitialized = false;
    gradientDelta = 30/180*Foam::constant::mathematical::pi;
    randomLambLimit = 0.01;
    //rndGen = uniformGeneratorOp();

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
    
    if (mag(s) < VSMALL) return GREAT; 

    const label own = mesh_.faceOwner()[facei];
    scalar ownVol; vector ownCc;
    calculateCellVolumeAndCentre(p, own, ownVol, ownCc);

    vector d;
    if (mesh_.isInternalFace(facei))
    {
        const label nei = mesh_.faceNeighbour()[facei];
        scalar neiVol; vector neiCc;
        calculateCellVolumeAndCentre(p, nei, neiVol, neiCc);
        d = neiCc - ownCc;
    }
    else
    {
        d = f.centre(p) - ownCc;
    }

    scalar dDotS = (d & s) / (mag(d) * mag(s) + VSMALL);
    dDotS = max(-1.0, min(1.0, dDotS)); 
    
    return radToDeg(acos(dDotS));
}

scalar meshFixer::calculateFaceSkewness
(
    const pointField& p,
    const label facei
) const
{
    const face& f = mesh_.faces()[facei];
    vector s = f.areaNormal(p);
    if (mag(s) < VSMALL) return GREAT;

    const label own = mesh_.faceOwner()[facei];
    scalar ownVol; vector ownCc;
    calculateCellVolumeAndCentre(p, own, ownVol, ownCc);
    
    vector fc = f.centre(p);
    vector d;
    vector dOwn = fc - ownCc;

    if (mesh_.isInternalFace(facei))
    {
        const label nei = mesh_.faceNeighbour()[facei];
        scalar neiVol; vector neiCc;
        calculateCellVolumeAndCentre(p, nei, neiVol, neiCc);
        d = neiCc - ownCc;
    }
    else
    {
        d = fc - ownCc;
    }

    scalar t = (s & dOwn) / ((s & d) + VSMALL);
    point intersection = ownCc + t * d;

    return mag(fc - intersection) / (mag(d) + VSMALL);
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

    // 1. Identify Active Topology
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

    // 2. Pre-calculate Old State
    Map<scalar> orthoOld(activeFaces.size());
    Map<scalar> skewOld(activeFaces.size());
    
    for (const label facei : activeFaces)
    {
        orthoOld.insert(facei, calculateFaceNonOrtho(p0, facei));
        skewOld.insert(facei, calculateFaceSkewness(p0, facei));
    }

    // 3. Local Lambda Setup
    scalarField lambda(mesh_.nPoints(), 1.0);
    const label maxIters = movementRelaxIters_; 
    const scalar lambdaReduction = movementRelaxation_;
    
    pointField pTest = p0;
    bool completelyValid = false;

    Info << "   -> Starting Local Guarded Morph on " 
         << returnReduce(activeFaces.size(), sumOp<label>()) << " active faces." << endl;

    // 4. Local Rollback Loop
    for (label iter = 0; iter < maxIters; ++iter)
    {
        
        syncTools::syncPointList(mesh_, lambda, minMagSqrEqOp<scalar>(), 1.0);

        forAll(pTest, ptI)
        {
            pTest[ptI] = p0[ptI] + (lambda[ptI] * disp[ptI]);
        }

        bool iterValid = true;
        
        // --- TELEMETRY COUNTERS ---
        label nFailVol = 0;
        label nFailArea = 0;
        label nFailQual = 0;
        label nPassAbs = 0;
        label nPassImp = 0;

        // --- TIER 1 EVALUATION (Cells) ---
        labelHashSet failedCells;
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

            //check if the cells failed
            label ownerCellI = mesh_.faceOwner()[facei];
            label neighbourCellI = mesh_.faceNeighbour()[facei];

            if(
                failedCells.found(ownerCellI)
                ||
                failedCells.found(neighbourCellI)
            ){
                faceFails = true;
            }

            // Tier 1 (Face Area)
            scalar areaTest = mag(mesh_.faces()[facei].areaNormal(pTest));
            if (areaTest <= VSMALL)
            {
                faceFails = true;
                nFailArea++;
            }
            else
            {
                scalar orthoNew = calculateFaceNonOrtho(pTest, facei);
                scalar skewNew  = calculateFaceSkewness(pTest, facei);

                // Granular Logic for Logging
                bool orthoPassAbs = (orthoNew <= qBad_);
                bool skewPassAbs  = (skewNew <= sBad_);
                
                bool orthoPassImp = (orthoNew < orthoOld[facei]);
                bool skewPassImp  = (skewNew < skewOld[facei]);

                bool orthoPass = orthoPassAbs || orthoPassImp;
                bool skewPass  = skewPassAbs  || skewPassImp;

                if (!orthoPass || !skewPass)
                {
                    faceFails = true;
                    nFailQual++;
                }
                else
                {
                    // If it passed, categorize *how* it passed
                    if (orthoPassAbs && skewPassAbs)
                    {
                        nPassAbs++;
                    }
                    else
                    {
                        // Passed purely because it is better than before
                        nPassImp++;
                    }
                }
            }

            if (faceFails)
            {
                iterValid = false;
                const face& f = mesh_.faces()[facei];
                forAll(f, i)
                {
                    label ptI = f[i];
                    if (activePointsSet.found(ptI))
                    {
                        lambda[ptI] *= lamSign*lambdaReduction;
                    }
                }
            }
        }

        lamFactor *= lambdaReduction;
        if(lamFactor < limit)
        {
            if(pendulumState == 0){
                lamSign = -1;
                pendulumState = 1;
                Info << "Trying displacement flip" << endl;
            }else if(pendulumState == 1){
                pendulumState = 2;
                Info << "Reverting displacement flig" << endl;
            }else if(pendulumState == 3){
                lamSign = 1;
                pendulumState = 4;
                Info << "Back to normal operation" << endl;
            }
        }
        
        // Parallel Reduction for Terminal Printouts
        reduce(nFailVol, sumOp<label>());
        reduce(nFailArea, sumOp<label>());
        reduce(nFailQual, sumOp<label>());
        reduce(nPassAbs, sumOp<label>());
        reduce(nPassImp, sumOp<label>());

        Info << "      Iter " << iter << " | "
             << "Pass(Target): " << nPassAbs << ", "
             << "Pass(Improved): " << nPassImp << " || "
             << "Fail(Vol): " << nFailVol << ", "
             << "Fail(Area): " << nFailArea << ", "
             << "Fail(Qual): " << nFailQual << endl;

        // FIX: Break based on GLOBAL failure counts, so all processors exit together.
        if ((nFailVol + nFailArea + nFailQual) == 0)
        {
            completelyValid = true;
            break;
        }

    }
    
    Info << "Relaxation done." << endl;

    // 5. Final Sync and Commit
    syncTools::syncPointList(mesh_, lambda, minMagSqrEqOp<scalar>(), 1.0);

    pointVectorField adaptedDisplacement(finalPointDisp);
    forAll(pTest, ptI)
    {
        pTest[ptI] = p0[ptI] + (lambda[ptI] * disp[ptI]);
        adaptedDisplacement[ptI] = lambda[ptI] * disp[ptI];
    }

    syncTools::syncPointList(mesh_, pTest, minEqOp<point>(), point(GREAT, GREAT, GREAT));

    // 2. NaN / Inf Check
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
    scalar maxActiveLambda = 0.0;
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

    // 3. Move the mesh
    if(move){
        mesh_.movePoints(pTest);
    }
    
    if(wiggleInitialized){
        Info << "Transfering wiggle lambda";
        DynamicList<scalar> newLambda(wigglePoints.size());
        labelList badPointIDs = wigglePoints.toc();
        forAll(badPointIDs,ind)
        {
            label pointI = badPointIDs[ind];
            scalar localLamb = lambda[pointI];
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
    //Find all cells associated with bad faces
    scalarField nonOrthoAngle = calculateFaceNonOrtho();
    scalarField skewness = calculateFaceSkewness();

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
    
    labelHashSet currentBadPoints;

    for (label facei = 0; facei < mesh_.nInternalFaces(); ++facei)
    {
        bool badFace = (nonOrthoAngle[facei] > qBad_ || skewness[facei] > sBad_);
        if(badFace)
        {
            const face& f = mesh_.faces()[facei];
            forAll(f, i)
            {
                label ptI = f[i];
                if(
                    !boundaryPoints.found(ptI)
                )
                currentBadPoints.insert(ptI);
            }
        }
    }

    labelList badPointIDs = currentBadPoints.toc();

    Info << "Wiggeling Points: " << currentBadPoints.size() << endl;

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

            theta_ = theta0 + (1-mag(lambda0))*Foam::sign(gradTheta)*gradientDelta;
            phi_ = phi0 + (1-mag(lambda0))*Foam::sign(gradPhi)*gradientDelta;
            
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

} // End namespace Foam```

---

### File: meshFixer.H
**Location:** `meshFixer.H`

```
#ifndef meshFixer_H
#define meshFixer_H

#include "fvMesh.H"
#include "volFields.H"
#include "pointFields.H"
#include "IOdictionary.H"
#include "boolList.H"

#include "motionSmoother.H"
#include "pointMesh.H"
#include "indirectPrimitivePatch.H"

#include "Map.H"
#include "DynamicList.H"
#include "Random.H"

namespace Foam
{

class meshFixer
{
    fvMesh& mesh_;
    Time& runTime_;

    // Workflow control
    bool discreteStep_;
    bool laplacianStep_;
    bool wiggleStep_;

    // Dictionary parameters
    label totalIterations_;
    bool writeIntermediate_;
    bool writeDisplacement_;
    scalar localRelaxation_;

    label nSmoothLayers_;
    label laplacianIters_;
    scalar localLaplacianRelaxation_;

    label movementRelaxIters_;
    scalar movementRelaxation_;

    scalar wVol_;
    scalar wLap_;
    scalar wExp_;
    scalar maxExpRatio_;
    scalar maxARShield_;

    scalar qGood_; 
    scalar qBad_;  
    scalar sGood_; 
    scalar sBad_;  
    scalar diffConstant_;
    scalar orthoSpring_;
    scalar skewSpring_;

    Random rndGen;

    dictionary qualityControlDict_;

    // Point mask for locked boundaries
    boolList isBoundaryPoint_;
    boolList isBadPoint_;

    // Private methods
    void identifyBoundaries();
    void calcLocalCorrections(pointVectorField& localDisp);
    void solveGlobalDiffusion(const pointVectorField& localDisp);
    void applyQualityGuardedMorph(const pointVectorField& finalPointDisp);
    pointVectorField applyRelativeGuardedMorph(const pointVectorField& finalPointDisp, bool move); // Fixed scoping
    void solveLocalLaplacian(const pointMesh& pMesh);
    scalar calculateAngleDelta(const scalar ang1, const scalar ang2);
    scalar normalizeAngle(const scalar ang);

    void wiggleSmooth(const pointMesh& pMesh);
    bool wiggleInitialized;
    //data for wiggleSmooth;
    labelHashSet wigglePoints;
    labelHashSet wigglePointsM1; //M1 inplies minus 1 aka. last step
    Map<label> pointToIdLabel;
    Map<label> pointToIdLabelM1;
    DynamicList<vector> wiggleDisplacement;
    DynamicList<scalar> wiggleLambda;
    DynamicList<scalar> wiggleLambdaM1;
    DynamicList<scalar> thetaWiggle;
    DynamicList<scalar> phiWiggle;
    DynamicList<scalar> thetaWiggleM1;
    DynamicList<scalar> phiWiggleM1;
    scalar wiggleFactor;
    scalar gradientDelta;
    scalar randomLambLimit;
    labelHashSet boundaryPoints;
    bool verboseWiggle_;

    scalarField calculateFaceNonOrtho() const;
    scalarField calculateFaceSkewness() const;
    scalar getMinEdgeLength(label ptI) const;

    scalar calculateCellAspectRatio(label celli) const;

    // Dynamic mesh metric functions for the Bisection Engine
    void calculateCellVolumeAndCentre
    (
        const pointField& p,
        const label celli,
        scalar& volume,
        vector& centre 
    ) const; // Fixed scoping

    scalar calculateFaceNonOrtho
    (
        const pointField& p,
        const label facei
    ) const; // Fixed scoping

    scalar calculateFaceSkewness
    (
        const pointField& p,
        const label facei
    ) const; // Fixed scoping

public:

    meshFixer(fvMesh& mesh);
    void fix();
};

} // End namespace Foam

#endif```

---

### File: meshMetrics.H
**Location:** `meshMetrics.H`

```
// ------------------------------------------------------------------------- //
// Helper: Dynamically calculate a cell's center and volume
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
    
    // 1. Estimate center via face centers
    vector cEst = Zero;
    for (const label facei : cFaces)
    {
        cEst += mesh_.faces()[facei].centre(p);
    }
    cEst /= cFaces.size();

    // 2. Sum face-pyramid contributions
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

// ------------------------------------------------------------------------- //
// Helper: Calculate Non-Orthogonality
// ------------------------------------------------------------------------- //
scalar meshFixer::calculateFaceNonOrtho
(
    const pointField& p,
    const label facei
) const
{
    // 1. OPENFOAM COMPLIANCE: Solid walls have perfect orthogonality (0 degrees)
    if (!mesh_.isInternalFace(facei))
    {
        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
        const label patchi = pbm.whichPatch(facei);
        
        // If it is an uncoupled boundary (wall, patch, empty), OpenFOAM ignores it.
        if (!pbm[patchi].coupled())
        {
            return 0.0; 
        }
    }

    const face& f = mesh_.faces()[facei];
    vector s = f.areaNormal(p);
    
    if (mag(s) < VSMALL) return GREAT; 

    const label own = mesh_.faceOwner()[facei];
    scalar ownVol; vector ownCc;
    calculateCellVolumeAndCentre(p, own, ownVol, ownCc);

    vector d;
    if (mesh_.isInternalFace(facei))
    {
        const label nei = mesh_.faceNeighbour()[facei];
        scalar neiVol; vector neiCc;
        calculateCellVolumeAndCentre(p, nei, neiVol, neiCc);
        d = neiCc - ownCc;
    }
    else
    {
        // For coupled patches (MPI/Cyclic), we approximate using the face center.
        // To be 100% compliant, you would need to sync dynamic cell centers 
        // across processors here, but for a local smoother, this approximation 
        // is standard to avoid massive MPI overhead inside the loop.
        d = f.centre(p) - ownCc;
    }

    // 2. OPENFOAM COMPLIANCE: Calculate Cosine and Clamp
    scalar dDotS = (d & s) / (mag(d) * mag(s) + VSMALL);
    dDotS = max(-1.0, min(1.0, dDotS)); 
    
    // Return the Angle (to match your qBad_ logic and checkMesh output)
    return radToDeg(acos(dDotS));
}

// ------------------------------------------------------------------------- //
// Helper: Calculate Skewness
// ------------------------------------------------------------------------- //
scalar meshFixer::calculateFaceSkewness
(
    const pointField& p,
    const label facei
) const
{
    const face& f = mesh_.faces()[facei];
    vector s = f.areaNormal(p);
    if (mag(s) < VSMALL) return GREAT;

    const label own = mesh_.faceOwner()[facei];
    scalar ownVol; vector ownCc;
    calculateCellVolumeAndCentre(p, own, ownVol, ownCc);
    
    vector fc = f.centre(p);
    vector d;
    vector dOwn = fc - ownCc;

    if (mesh_.isInternalFace(facei))
    {
        const label nei = mesh_.faceNeighbour()[facei];
        scalar neiVol; vector neiCc;
        calculateCellVolumeAndCentre(p, nei, neiVol, neiCc);
        d = neiCc - ownCc;
    }
    else
    {
        d = fc - ownCc;
    }

    // Intersect d with the face plane
    scalar t = (s & dOwn) / ((s & d) + VSMALL);
    vector pInt = ownCc + t * d;

    return mag(fc - pInt) / (mag(d) + VSMALL);
}```

---

### File: wrapCode.sh
**Location:** `wrapCode.sh`

```
#!/bin/bash

OUTPUT="codebase_context.md"

# 1. Initialize file and generate the directory tree
echo "# Project Structure" > "$OUTPUT"
echo "" >> "$OUTPUT"
echo '```text' >> "$OUTPUT"
# -a includes hidden files. -I ignores common noisy directories and the output file.
tree -a -I ".git|node_modules|venv|__pycache__|$OUTPUT" >> "$OUTPUT"
echo '```' >> "$OUTPUT"
echo -e "\n---\n" >> "$OUTPUT"

echo "# File Contents" >> "$OUTPUT"

# 2. Iterate through files and append their contents
find . -type f \
    -not -path '*/\.git/*' \
    -not -path '*/node_modules/*' \
    -not -path '*/venv/*' \
    -not -path '*/__pycache__/*' \
    -not -name "$OUTPUT" \
    | sort | while read -r file; do
    
    # Extract just the filename and the clean relative path
    filename=$(basename "$file")
    filepath="${file#./}"
    
    # Append the formatted data to the Markdown file
    echo "### File: $filename" >> "$OUTPUT"
    echo "**Location:** \`$filepath\`" >> "$OUTPUT"
    echo "" >> "$OUTPUT"
    echo '```' >> "$OUTPUT"
    cat "$file" >> "$OUTPUT"
    echo '```' >> "$OUTPUT"
    echo -e "\n---\n" >> "$OUTPUT"
done

echo "Done! Your context file has been generated: $OUTPUT"```

---

