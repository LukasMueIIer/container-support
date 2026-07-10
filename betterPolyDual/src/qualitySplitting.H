//Mark all low quality cells to be devided, merging them would only result in problems
if(
    useQualitySplitting
){
    Info << "Calculating Quality Metrics" << endl;
    labelHashSet lowQualityPoints;
    labelHashSet lowQualityEdges;

    scalarField faceOrthogonality;
    scalarField nonOrthoAngle;
    scalarField faceSkewness;

    volScalarField aspectRatio(
        IOobject
        (
            "aspectRatio",                                                                                                                                                                         
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar(dimless, Zero)
    );
    volScalarField aspectRatioCell(
        IOobject
        (
            "cellAspectRatio",                                                                                                                                                                         
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar(dimless, Zero)
    );
    volScalarField minTetVolume
    (
        IOobject
        (
            "minTetVolume",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar("minTetVolume", dimless, GREAT),
        fvPatchFieldBase::zeroGradientType()
    );
    scalarField cellOpenness;
    faceOrthogonality = polyMeshTools::faceOrthogonality
    (
        mesh,
        mesh.faceAreas(),
        mesh.cellCentres()
    );
    nonOrthoAngle = radToDeg
    (
        Foam::acos(min(scalar(1), max(scalar(-1), faceOrthogonality)))
    );
    faceSkewness = polyMeshTools::faceSkewness
    (
        mesh,
        mesh.points(),
        mesh.faceCentres(),
        mesh.faceAreas(),
        mesh.cellCentres()
    );
    polyMeshTools::cellClosedness
    (                                                                                                                                                                                             
        mesh,
        mesh.geometricD(),
        mesh.faceAreas(),
        mesh.cellVolumes(),
        cellOpenness,
        aspectRatio
    );
    aspectRatio.correctBoundaryConditions(); 
    aspectRatioCell.ref().field() = cellAspectRatio(meshPoly);  
    aspectRatioCell.correctBoundaryConditions(); 
    //calc min tet vol
    {
        const labelList& own = mesh.faceOwner();
        const labelList& nei = mesh.faceNeighbour();
        const pointField& p = mesh.points();
        forAll(own, facei)
        {
            const face& f = mesh.faces()[facei];
            const point& fc = mesh.faceCentres()[facei];
            {
                const point& ownCc = mesh.cellCentres()[own[facei]];
                scalar& ownVol = minTetVolume[own[facei]];
                forAll(f, fp)
                {
                    scalar tetQual = tetPointRef
                    (
                        p[f[fp]],
                        p[f.nextLabel(fp)],
                        ownCc,
                        fc
                    ).quality();
                    ownVol = min(ownVol, tetQual);
                }
            }
            if (mesh.isInternalFace(facei))
            {
                const point& neiCc = mesh.cellCentres()[nei[facei]];
                scalar& neiVol = minTetVolume[nei[facei]];
                forAll(f, fp)
                {
                    scalar tetQual = tetPointRef
                    (
                        p[f[fp]],
                        p[f.nextLabel(fp)],
                        fc,
                        neiCc
                    ).quality();
                    neiVol = min(neiVol, tetQual);
                }
            }
        }
        minTetVolume.correctBoundaryConditions();
    }
    //marking all low quality regions
    for (label faceI = 0; faceI < mesh.faceNeighbour().size(); faceI++){
        
            //Get onwer an neighbour
            label ownerCellI = mesh.faceOwner()[faceI];
            label neighbourCellI = mesh.faceNeighbour()[faceI];
            
            //Check if both arent tet, so we take these faces out of the poly process
        
            //Check if cell is of low quality
            bool cellIsLowQuality = false;
            if(
                aspectRatio[ownerCellI] > aspectRatioLimit
                ||
                aspectRatio[neighbourCellI] > aspectRatioLimit
                ||
                aspectRatioCell[ownerCellI] > aspectRatioLimit
                ||
                aspectRatioCell[neighbourCellI] > aspectRatioLimit
                ||
                minTetVolume[ownerCellI] < minTetVolumeLimit
                ||
                minTetVolume[neighbourCellI] < minTetVolumeLimit
            ){
                cellIsLowQuality = true;
            }
            //here the check
            const face& f = mesh.faces()[faceI];
            const labelList& fEdges = mesh.faceEdges()[faceI];
            forAll(f, fp)
            {
                //check if face is of low quality or parent cell is of low quality
                if(
                    (nonOrthoAngle[faceI] > nonOrthoAngleLimit)
                    ||
                    (faceSkewness[faceI] > skewnessLimit)
                    ||
                    cellIsLowQuality
                ){
                    if(
                        detectionVerbose
                    ){
                        Info << "Splitting low quality face " << faceI << endl;
                    }
                    lowQualityPoints.insert(f[fp]);
                    lowQualityEdges.insert(fEdges[fp]);
                }
            }
    }
    //extend the low quality region if needed
    label levelsLeft = qualityTransitionLevels;
    while (levelsLeft > 0)
    {
        Info<< "Extending low quality region" << endl;
        levelsLeft--;
        labelHashSet extendedLowQualityPoints;
        labelHashSet extendedLowQualityEdges;
          //now if any part of our face is part of low quality we split it
        for (label faceI = 0; faceI < mesh.faceNeighbour().size(); faceI++){

                //Get onwer an neighbour
                label ownerCellI = mesh.faceOwner()[faceI];
                label neighbourCellI = mesh.faceNeighbour()[faceI];

                //Check if both arent tet, so we take these faces out of the poly process

                //Check if cell is of low quality
                bool lowQualityRegion = false;
                const face& f = mesh.faces()[faceI];
                const labelList& fEdges = mesh.faceEdges()[faceI];
                forAll(f, fp)
                {
                    if((
                        lowQualityPoints.found(f[fp]) 
                        ||
                        lowQualityEdges.found(fEdges[fp]) 
                    ))
                    {
                       lowQualityRegion = true;
                    }
                }
                if(
                    lowQualityRegion
                )
                {
                    forAll(f, fp)
                    {
                        extendedLowQualityPoints.insert(f[fp]);
                        extendedLowQualityEdges.insert(fEdges[fp]);
                    }
                }                
        }
        lowQualityPoints = extendedLowQualityPoints;
        lowQualityEdges = extendedLowQualityEdges;
    }
    
    //now if any part of our face is part of low quality we split it
    for (label faceI = 0; faceI < mesh.faceNeighbour().size(); faceI++){
        
            //Get onwer an neighbour
            label ownerCellI = mesh.faceOwner()[faceI];
            label neighbourCellI = mesh.faceNeighbour()[faceI];
            
            //Check if both arent tet, so we take these faces out of the poly process
        
            //Check if cell is of low quality
            bool lowQualityRegion = false;
            const face& f = mesh.faces()[faceI];
            const labelList& fEdges = mesh.faceEdges()[faceI];
            forAll(f, fp)
            {
                if((
                    lowQualityPoints.found(f[fp]) 
                    ||
                    lowQualityEdges.found(fEdges[fp]) 
                ))
                {
                   lowQualityRegion = true;
                }
            }
            if(
                lowQualityRegion
            )
            {
                forAll(f, fp)
                {
                    //check if face is of low quality or parent cell is of low quality
                    singleCellFeaturePointSet.erase(f[fp]);
                    multiCellFeaturePointSet.insert(f[fp]);
                    featureEdgeSet.insert(fEdges[fp]);
                    featureFaceSet.insert(faceI);
                }
            }                
    }
}