/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    polyDualMesh

Group
    grpMeshManipulationUtilities

Description
    Creates the dual of a polyMesh, adhering to all the feature and patch edges.

Usage
    \b polyDualMesh featureAngle

    Detects any boundary edge > angle and creates multiple boundary faces
    for it. Normal behaviour is to have each point become a cell
    (1.5 behaviour)

    Options:
      - \par -concaveMultiCells
        Creates multiple cells for each point on a concave edge. Might limit
        the amount of distortion on some meshes.

      - \par -splitAllFaces
        Normally only constructs a single face between two cells. This single
        face might be too distorted. splitAllFaces will create a single face for
        every original cell the face passes through. The mesh will thus have
        multiple faces in between two cells! (so is not strictly
        upper-triangular anymore - checkMesh will complain)

      - \par -doNotPreserveFaceZones:
        By default all faceZones are preserved by marking all faces, edges and
        points on them as features. The -doNotPreserveFaceZones disables this
        behaviour.

Note
    It is just a driver for meshDualiser. Substitute your own simpleMarkFeatures
    to have different behaviour.

\*---------------------------------------------------------------------------*/

#define dumpFeaturesFlag false


//---------------------------------------------------------------------------

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "unitConversion.H"
#include "polyTopoChange.H"
#include "mapPolyMesh.H"
#include "bitSet.H"
#include "meshTools.H"
#include "OFstream.H"
#include "meshDualiser.H"
#include "ReadFields.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "topoSet.H"
#include "processorMeshes.H"

//New includes
#include "tetMatcher.H"
#include "hexMatcher.H"
#include "polyMeshTools.H"
#include "cellAspectRatio.H"
#include "tetrahedron.H"


using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Naive feature detection. All boundary edges with angle > featureAngle become
// feature edges. All points on feature edges become feature points. All
// boundary faces become feature faces.
void simpleMarkFeatures
(
    const fvMesh& mesh,
    const polyMesh& meshPoly,
    const bitSet& isBoundaryEdge,
    const scalar featureAngle,
    const bool concaveMultiCells,
    const bool doNotPreserveFaceZones,

    labelList& featureFaces,
    labelList& featureEdges,
    labelList& singleCellFeaturePoints,
    labelList& multiCellFeaturePoints
)
{

    #include "readFromDict.H"



    scalar minCos = Foam::cos(degToRad(featureAngle));

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Working sets
    labelHashSet featureEdgeSet;
    labelHashSet singleCellFeaturePointSet;
    labelHashSet multiCellFeaturePointSet;

    // Face centres that need inclusion in the dual mesh
    labelHashSet featureFaceSet(mesh.nFaces());

    // 1. Mark all edges between patches
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelHashSet boundaryPointTracker;

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];
        const labelList& meshEdges = pp.meshEdges();

        // All patch corner edges. These need to be feature points & edges!
        // To prevent the "overlaying" faces, keep all edges and points of the patch
        
        label edgeIStart = pp.nInternalEdges();

        if(flagBoundaryHandlingMethod == 1){ 
            Info << "Retaining all boundary features" << endl;
            edgeIStart = 0;
        }

        for (label edgeI = edgeIStart; edgeI < pp.nEdges(); edgeI++)
        {
            label meshEdgeI = meshEdges[edgeI];
            featureEdgeSet.insert(meshEdgeI);
            multiCellFeaturePointSet.insert(mesh.edges()[meshEdgeI][0]);
            multiCellFeaturePointSet.insert(mesh.edges()[meshEdgeI][1]);
            singleCellFeaturePointSet.erase(mesh.edges()[meshEdgeI][0]);
            singleCellFeaturePointSet.erase(mesh.edges()[meshEdgeI][1]);
            boundaryPointTracker.insert(mesh.edges()[meshEdgeI][0]);
            boundaryPointTracker.insert(mesh.edges()[meshEdgeI][1]);
        }
    }

    if(
        flagBoundaryHandlingMethod == 2
    )
    {
        Info << "Quality based boundary retention, scanning for face quality" << endl;
        Info << "Boundary orthogonality detection limit: " << boundaryOrthoDetectionLimit << endl;

        //calculating the quality measures
        scalarField faceOrthogonality;
        scalarField nonOrthoAngle;
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
        
        //Go through all faces in the mesh (could improve that  ...) and mark all points that are part of low quality faces
        labelHashSet lowQualityPoints;
        for (label faceI = 0; faceI < mesh.faceNeighbour().size(); faceI++){
            
                
                const face& f = mesh.faces()[faceI];

                forAll(f, fp)
                {
                    //check if face is of low quality or parent cell is of low quality
                    if(
                        (nonOrthoAngle[faceI] > boundaryOrthoDetectionLimit)
                    ){
                        if(
                            detectionVerbose
                        ){
                            Info << "Added low quality Point: " << f[fp] << endl;
                        }

                        lowQualityPoints.insert(f[fp]);
                    }
                }
        }


        //now go through boundary faces and retain them if they contain a lowQualityPoint
        forAll(patches, patchi)
        {
            const polyPatch& pp = patches[patchi];
            const labelList& meshEdges = pp.meshEdges();

            label edgeIStart = 0;

            for (label edgeI = edgeIStart; edgeI < pp.nInternalEdges(); edgeI++)
            {
                label meshEdgeI = meshEdges[edgeI];
                if(
                    lowQualityPoints.found(mesh.edges()[meshEdgeI][0]) 
                    ||
                    lowQualityPoints.found(mesh.edges()[meshEdgeI][1])
                )
                {
                    if(
                        detectionVerbose
                    ){
                        Info << "Added edge because it has low quality points: " << edgeI << endl;
                    }

                    featureEdgeSet.insert(meshEdgeI);
                    multiCellFeaturePointSet.insert(mesh.edges()[meshEdgeI][0]);
                    multiCellFeaturePointSet.insert(mesh.edges()[meshEdgeI][1]);
                    singleCellFeaturePointSet.erase(mesh.edges()[meshEdgeI][0]);
                    singleCellFeaturePointSet.erase(mesh.edges()[meshEdgeI][1]);
                    boundaryPointTracker.insert(mesh.edges()[meshEdgeI][0]);
                    boundaryPointTracker.insert(mesh.edges()[meshEdgeI][1]);
                }
                
            }
        }
        

    }

    if(
        flagNonTetCellMethod != 0
    ){
        Info << "Treating non-tet cells" << endl;

        for (label faceI = 0; faceI < mesh.faceNeighbour().size(); faceI++){
            
            //Get onwer an neighbour
            label ownerCellI = mesh.faceOwner()[faceI];
            label neighbourCellI = mesh.faceNeighbour()[faceI];
            
            //Check if the face belongs to at least one non tet cell
            if(
                (
                    (!tetMatcher::test(mesh,ownerCellI))
                    ||
                    (!tetMatcher::test(mesh,neighbourCellI))
                )
                &&  
                (   //check that we dont have a bad precident for later hex cell special treatment
                    (
                        (!hexMatcher::test(mesh,ownerCellI))
                        &&
                        (!hexMatcher::test(mesh,neighbourCellI))
                    )
                    ||
                    (
                        flagHexCellMethod == 0
                    )
                )    
            )
            {
                if(
                    detectionVerbose 
                ){
                    Info << "Treating face: " << faceI << " belonging to non-tet cell" << endl;
                }

                const face& f = mesh.faces()[faceI];
                const labelList& fEdges = mesh.faceEdges()[faceI];

                forAll(f, fp)
                {

                    if(
                       flagNonTetCellMethod == 1 
                    )
                    {   //single cell splitting
                        if( //check we werent multi split before
                            !(multiCellFeaturePointSet.found(f[fp]))
                        ){
                            singleCellFeaturePointSet.insert(f[fp]);
                            //featureEdgeSet.insert(fEdges[fp]);
                            //featureFaceSet.insert(faceI);
                        }
                    }
                    else if(
                        flagNonTetCellMethod == 2
                    )
                    {   //multi cell splitting
                        singleCellFeaturePointSet.erase(f[fp]);
                        multiCellFeaturePointSet.insert(f[fp]);
                        featureEdgeSet.insert(fEdges[fp]);
                        featureFaceSet.insert(faceI);
                    }
                }

            }
        }
    }

    if(
        flagHexCellMethod != 0
    ){
        Info << "Specially treating hex cells" << endl;

        for (label faceI = 0; faceI < mesh.faceNeighbour().size(); faceI++){
            
            //Get onwer an neighbour
            label ownerCellI = mesh.faceOwner()[faceI];
            label neighbourCellI = mesh.faceNeighbour()[faceI];
            
            //Check if the face belongs to at least one non tet cell
            if(
                (hexMatcher::test(mesh,ownerCellI))
                ||
                (hexMatcher::test(mesh,neighbourCellI))
            )
            {
                if(
                    detectionVerbose 
                ){
                    Info << "Treating face: " << faceI << " belonging to hex cell" << endl;
                }

                const face& f = mesh.faces()[faceI];
                const labelList& fEdges = mesh.faceEdges()[faceI];

                forAll(f, fp)
                {

                    if(
                       flagHexCellMethod == 1 
                    )
                    {   
                        if( //check we werent multi split before, if yes its to maintain boundary and takes precedence
                            !(multiCellFeaturePointSet.found(f[fp]))
                        ){
                            singleCellFeaturePointSet.insert(f[fp]);
                            //featureEdgeSet.insert(fEdges[fp]);
                            //featureFaceSet.insert(faceI);
                        }
                    }
                    else if(
                        flagHexCellMethod == 2
                    )
                    {   //multi cell splitting
                        singleCellFeaturePointSet.erase(f[fp]);
                        multiCellFeaturePointSet.insert(f[fp]);
                        featureEdgeSet.insert(fEdges[fp]);
                        featureFaceSet.insert(faceI);
                    }
                }

            }
        }
    }

    //interface face treatment (between tet and non tet cells)
    if(
        flagInterfaceTreatmentMethod != 0
    ){
        Info << "Treating interface faces" << endl;

        for (label faceI = 0; faceI < mesh.faceNeighbour().size(); faceI++){
            
            //Get onwer an neighbour
            label ownerCellI = mesh.faceOwner()[faceI];
            label neighbourCellI = mesh.faceNeighbour()[faceI];
            
            //Check if the face belongs to at least one non tet cell
            if(
                (
                    (!tetMatcher::test(mesh,ownerCellI))
                    &&
                    (tetMatcher::test(mesh,neighbourCellI))
                )
                ||
                (
                    (tetMatcher::test(mesh,ownerCellI))
                    &&
                    (!tetMatcher::test(mesh,neighbourCellI))
                )
            )
            {
                if(
                    detectionVerbose 
                ){
                    Info << "Treating face: " << faceI << " belonging to interface" << endl;
                }

                const face& f = mesh.faces()[faceI];
                const labelList& fEdges = mesh.faceEdges()[faceI];

                forAll(f, fp)
                {

                    if(
                       flagInterfaceTreatmentMethod == 1 
                    )
                    {   //single cell splitting
                        multiCellFeaturePointSet.erase(f[fp]);
                        singleCellFeaturePointSet.insert(f[fp]);
                        featureEdgeSet.erase(fEdges[fp]);
                        featureFaceSet.erase(faceI);
                    }
                    else if(
                        flagInterfaceTreatmentMethod == 2
                    )
                    {   //multi cell splitting
                        singleCellFeaturePointSet.erase(f[fp]);
                        multiCellFeaturePointSet.insert(f[fp]);
                        featureEdgeSet.insert(fEdges[fp]);
                        featureFaceSet.insert(faceI);
                    }
                }

            }
        }
    }


    //interface face treatment (between tet and hex cells)
    if(
        flagHexInterfaceTreatmentMethod != 0
    ){
        Info << "Special treatment for hex-non hex interface faces" << endl;

        labelHashSet hexNonHexCells;
        
        Info << "Marking all non hex cells which are connected to a hex cell" << endl;

        for (label faceI = 0; faceI < mesh.faceNeighbour().size(); faceI++){
            
            //Get onwer an neighbour
            label ownerCellI = mesh.faceOwner()[faceI];
            label neighbourCellI = mesh.faceNeighbour()[faceI];
            
            //Check if the face belongs to a tet and a hex cell
            if(
                (
                    (hexMatcher::test(mesh,ownerCellI))
                    &&
                    (!hexMatcher::test(mesh,neighbourCellI))
                )
            )
            {
                if(
                    detectionVerbose 
                ){
                    Info << "Marking Cell: " << neighbourCellI << " belonging to hex-non interface" << endl;
                }
                hexNonHexCells.insert(neighbourCellI);
            }
            else if(
                (
                    (!hexMatcher::test(mesh,ownerCellI))
                    &&
                    (hexMatcher::test(mesh,neighbourCellI))
                )
            )
            {
                if(
                    detectionVerbose 
                ){
                    Info << "Marking Cell: " << ownerCellI << " belonging to hex-non interface" << endl;
                }
                hexNonHexCells.insert(ownerCellI);
            }
        }

        Info << "Splitting all faces that belong to marked cells" << endl;

        for (label faceI = 0; faceI < mesh.faceNeighbour().size(); faceI++){
            
            //Get onwer an neighbour
            label ownerCellI = mesh.faceOwner()[faceI];
            label neighbourCellI = mesh.faceNeighbour()[faceI];
            
            //Check if the face belongs to a tet and a hex cell
            if(
                (
                    hexNonHexCells.find(ownerCellI)
                    ||
                    hexNonHexCells.find(neighbourCellI)
                )
            )
            {
                if(
                    detectionVerbose 
                ){
                    Info << "Treating face: " << faceI << " belonging to marked cell" << endl;
                }

                const face& f = mesh.faces()[faceI];
                const labelList& fEdges = mesh.faceEdges()[faceI];

                forAll(f, fp)
                {

                    if(
                       flagHexInterfaceTreatmentMethod == 1 
                    )
                    {   //single cell splitting
                        multiCellFeaturePointSet.erase(f[fp]);
                        singleCellFeaturePointSet.insert(f[fp]);
                        featureEdgeSet.erase(fEdges[fp]);
                        featureFaceSet.erase(faceI);
                    }
                    else if(
                        flagHexInterfaceTreatmentMethod == 2
                    )
                    {   //multi cell splitting
                        singleCellFeaturePointSet.erase(f[fp]);
                        multiCellFeaturePointSet.insert(f[fp]);
                        featureEdgeSet.insert(fEdges[fp]);
                        featureFaceSet.insert(faceI);
                    }
                }

            }
        }
    }

    if(
        extrudedClimbBoundarySplitting
    ){
        //run a fake feature edge detection run to add these to boundaryPointTracker
        //run a "fake" feature edge detection to add these cells as well
        // Check for features.
        //build a list of all cells that are boundary cells
        labelHashSet boundaryCellTracker;

        labelHashSet featureEdgePoints;
        primitivePatch allBoundary
        (
            SubList<face>
            (
                mesh.faces(),
                mesh.nBoundaryFaces(),
                mesh.nInternalFaces()
            ),
            mesh.points()
        );

        const labelListList& edgeFaces = allBoundary.edgeFaces();
        const labelList& meshPoints = allBoundary.meshPoints();

        forAll(edgeFaces, edgeI)
        {
            const labelList& eFaces = edgeFaces[edgeI];
        
            if (eFaces.size() == 2)
            {
                label f0 = eFaces[0];
                label f1 = eFaces[1];
            
                // check angle
                const vector& n0 = allBoundary.faceNormals()[f0];
                const vector& n1 = allBoundary.faceNormals()[f1];
            
                if ((n0 & n1) < minCos)
                {
                    const edge& e = allBoundary.edges()[edgeI];
                    label v0 = meshPoints[e[0]];
                    label v1 = meshPoints[e[1]];

                    label meshEdgeI = meshTools::findEdge(mesh, v0, v1);
                    featureEdgePoints.insert(mesh.edges()[meshEdgeI][0]);
                    featureEdgePoints.insert(mesh.edges()[meshEdgeI][1]);
                    if(
                        detectionVerbose
                    ){
                        Info << "Added feature edge points: " << mesh.edges()[meshEdgeI][0] << " " << mesh.edges()[meshEdgeI][1] << endl;
                    }
                }
            }
        }
        //since we only added edges we loop over all faces now and porperly add all associated faces
        for (label faceI = mesh.faceNeighbour().size(); faceI < mesh.faceOwner().size(); faceI++){
            
                const face& f = mesh.faces()[faceI];
                const labelList& fEdges = mesh.faceEdges()[faceI];
                bool isRelated = false;

                forAll(f, fp) //we check that all points of our face are part of the boundary points         
                { //if not we are not a boundary face
                    //check if face is of low quality or parent cell is of low quality
                    if(
                        featureEdgePoints.found(f[fp])
                    ){
                        isRelated = true;
                        break;
                    }
                }
                if(
                    isRelated
                ){
                    boundaryCellTracker.insert(mesh.faceOwner()[faceI]);
                    if(
                        detectionVerbose
                    ){
                        Info << "Added feature edge cell: " << mesh.faceOwner()[faceI] << endl;
                    }
                    forAll(f, fp) //we check that all points of our face are part of the boundary points         
                    { //if not we are not a boundary face
                        //check if face is of low quality or parent cell is of low quality
                        singleCellFeaturePointSet.erase(f[fp]);
                        multiCellFeaturePointSet.insert(f[fp]);
                        featureEdgeSet.insert(fEdges[fp]);
                        featureFaceSet.insert(faceI);
                        boundaryPointTracker.insert(f[fp]);
                    }
                }
        }

        for (label faceI = 0; faceI < mesh.faceOwner().size(); faceI++){
            
                const face& f = mesh.faces()[faceI];
                bool isBoundary = true;

                forAll(f, fp) //we check that all points of our face are part of the boundary points         
                { //if not we are not a boundary face
                    //check if face is of low quality or parent cell is of low quality
                    if(
                        !(boundaryPointTracker.found(f[fp]))
                    ){
                        isBoundary = false;
                        break;
                    }
                }
                if(
                    isBoundary
                ){
                    boundaryCellTracker.insert(mesh.faceOwner()[faceI]);
                    if(
                        detectionVerbose
                    ){
                        Info << "Added boundary cell: " << mesh.faceOwner()[faceI] << endl;
                    }
                    forAll(f, fp) //we check that all points of our face are part of the boundary points         
                    { //if not we are not a boundary face
                        //check if face is of low quality or parent cell is of low quality
                        if(
                            !(boundaryPointTracker.found(f[fp]))
                        ){
                            isBoundary = false;
                            break;
                        }
                    }
                }
        }

        label internalMax = maxClimbingIterations;
        label internalTracker = 0;

        //now loop over the cells and perform the climbing
        while(
            internalTracker < internalMax
        ){
            internalTracker++;
            Info << "Climbing Iteration " << internalTracker << endl;;
            //housekeeping
            label additionTracker = 0;
            labelHashSet localBoundaryCellTracker = boundaryCellTracker;
            labelHashSet localBoundaryPointTracker = boundaryPointTracker;

            //loop over faces
            for (label faceI = 0; faceI < mesh.faceNeighbour().size(); faceI++){
                //check if we are assiciated with marked cell
                if(
                    (boundaryCellTracker.found(mesh.faceOwner()[faceI]))
                    ||
                    (boundaryCellTracker.found(mesh.faceNeighbour()[faceI]))
                ){
                    if(
                        detectionVerbose
                    ){
                        Info << "Found associated cells for face: " << faceI << endl;
                    }  

                    //check that we are not in the boundary cloud
                    bool isNotContained = true;
                    const face& f = mesh.faces()[faceI];
                    const labelList& fEdges = mesh.faceEdges()[faceI];

                    forAll(f, fp) //we check that all points of our face are part of the boundary points         
                    { //if not we are not a boundary face
                        //check if face is of low quality or parent cell is of low quality
                        if(
                            boundaryPointTracker.found(f[fp])
                        ){
                            isNotContained = false;
                            break;
                        }
                    }

                    if(
                        isNotContained
                    ){
                        if(
                        detectionVerbose
                        ){
                            Info << "Face will be added"<< endl;
                        } 
                        additionTracker++;

                        //add to local boundary cells
                        localBoundaryCellTracker.insert(mesh.faceOwner()[faceI]);
                        localBoundaryCellTracker.insert(mesh.faceNeighbour()[faceI]);

                        //loop over points and add details
                        forAll(f, fp)
                        {
                            //multi cell splitting
                            singleCellFeaturePointSet.erase(f[fp]);
                            multiCellFeaturePointSet.insert(f[fp]);
                            featureEdgeSet.insert(fEdges[fp]);
                            featureFaceSet.insert(faceI);
                            localBoundaryPointTracker.insert(f[fp]);
                        }

                    }else{
                        if(
                        detectionVerbose
                        ){
                            Info << "Face disregarded because points were found"<< endl;
                        } 
                    }

                }
            }

            //transfer data
            if(
                additionTracker == 0
            ){
                Info << "Iteration did not add any cells, stopping climbing" << endl;
                internalMax = 0;
            }else{
                Info << "Iteration added " << additionTracker << " faces" << endl;
            }
            boundaryCellTracker = localBoundaryCellTracker;
            boundaryPointTracker = localBoundaryPointTracker;
        }

    }



    if(
        interfaceBufferLevel != 0
    )
    {
        Info << "Creating Interfacfe Buffer Layers" << endl;
        labelHashSet taggedPoints = multiCellFeaturePointSet;
        //copy via loop cause we dont trust
        //forAll(multiCellFeaturePointSet, fp)
        //{
        //    taggedPoints.insert(fp);
        //}
        Info << "Currently Tainted Points: " << taggedPoints.size() << endl;
        labelHashSet taintedFaces;

        label iterationCount = 0;
        while(
            iterationCount < interfaceBufferLevel
        )
        {
            labelHashSet newTaggedPoints = taggedPoints;
            iterationCount++;
            Info << "Tagging Iteration: " << iterationCount << endl;
            //loop over all faces, check if they are purly tet owned
            //if yes check if any of its point belong to the taggedPoints or are multi split
            for (label faceI = 0; faceI < mesh.faceNeighbour().size(); faceI++){
            
                //Get onwer an neighbour
                label ownerCellI = mesh.faceOwner()[faceI];
                label neighbourCellI = mesh.faceNeighbour()[faceI];
                
                //Check if the face belongs to at least one non tet cell
                if(
                    (tetMatcher::test(mesh,ownerCellI))
                    &&
                    (tetMatcher::test(mesh,neighbourCellI))
                )
                {
                    const face& f = mesh.faces()[faceI];

                    bool isTainted = false;

                    forAll(f, fp)
                    {
                        if(
                            taggedPoints.found(f[fp])
                        )
                        {   
                            isTainted = true;
                            if(
                                detectionVerbose
                            ){
                                Info << "Tagged face: " << faceI << endl;
                            }
                        }
                    }

                    if(
                        isTainted
                    ){
                        taintedFaces.insert(faceI);
                        forAll(f, fp)
                        {
                            newTaggedPoints.insert(f[fp]);
                        }
                    }

                }
            }
            taggedPoints = newTaggedPoints;
            Info << "Currently Tainted Points: " << taggedPoints.size() << endl;
            Info << "Currently tagged Faces: " << taintedFaces.size() << endl;
        }
            Info << "Splitting all tagged faces" << endl;

            //now access all faces that have been tainted and add them to multisplit
            forAll(
                taintedFaces, faceI
            )
            {
                const face& f = mesh.faces()[faceI];
                const labelList& fEdges = mesh.faceEdges()[faceI];
                if(
                    detectionVerbose
                ){
                    Info << "Splitting face: " << faceI << endl;
                }

                forAll(f, fp)
                {
                    //multi cell splitting
                    singleCellFeaturePointSet.erase(f[fp]);
                    multiCellFeaturePointSet.insert(f[fp]);
                    featureEdgeSet.insert(fEdges[fp]);
                    featureFaceSet.insert(faceI);
                }

            }

    }

    #include "qualitySplitting.H"
    


    // 2. Mark all geometric feature edges
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Make distinction between convex features where the boundary point becomes
    // a single cell and concave features where the boundary point becomes
    // multiple 'half' cells.

    // Addressing for all outside faces
    primitivePatch allBoundary
    (
        SubList<face>
        (
            mesh.faces(),
            mesh.nBoundaryFaces(),
            mesh.nInternalFaces()
        ),
        mesh.points()
    );

    // Check for non-manifold points (surface pinched at point)
    allBoundary.checkPointManifold(false, &singleCellFeaturePointSet);

    forAllConstIters(multiCellFeaturePointSet, iter)
    {
        if (singleCellFeaturePointSet.found(iter.key()))
        {
            singleCellFeaturePointSet.erase(iter.key());
        }
    }

    // Check for non-manifold edges (surface pinched at edge)
    const labelListList& edgeFaces = allBoundary.edgeFaces();
    const labelList& meshPoints = allBoundary.meshPoints();

    forAll(edgeFaces, edgeI)
    {
        const labelList& eFaces = edgeFaces[edgeI];

        if (eFaces.size() > 2)
        {
            const edge& e = allBoundary.edges()[edgeI];

            //Info<< "Detected non-manifold boundary edge:" << edgeI
            //    << " coords:"
            //    << allBoundary.points()[meshPoints[e[0]]]
            //    << allBoundary.points()[meshPoints[e[1]]] << endl;

            if (!multiCellFeaturePointSet.found(meshPoints[e[0]]))
            {
                singleCellFeaturePointSet.insert(meshPoints[e[0]]);
            }
            if (!multiCellFeaturePointSet.found(meshPoints[e[1]]))
            {
                singleCellFeaturePointSet.insert(meshPoints[e[1]]);
            }

        }
    }

    // Check for features.
    forAll(edgeFaces, edgeI)
    {
        const labelList& eFaces = edgeFaces[edgeI];

        if (eFaces.size() == 2)
        {
            label f0 = eFaces[0];
            label f1 = eFaces[1];

            // check angle
            const vector& n0 = allBoundary.faceNormals()[f0];
            const vector& n1 = allBoundary.faceNormals()[f1];

            if ((n0 & n1) < minCos)
            {
                const edge& e = allBoundary.edges()[edgeI];
                label v0 = meshPoints[e[0]];
                label v1 = meshPoints[e[1]];

                label meshEdgeI = meshTools::findEdge(mesh, v0, v1);
                featureEdgeSet.insert(meshEdgeI);
                singleCellFeaturePointSet.erase(mesh.edges()[meshEdgeI][0]);
                singleCellFeaturePointSet.erase(mesh.edges()[meshEdgeI][1]);
                multiCellFeaturePointSet.insert(mesh.edges()[meshEdgeI][0]);
                multiCellFeaturePointSet.insert(mesh.edges()[meshEdgeI][1]);

                // Check if convex or concave by looking at angle
                // between face centres and normal
                vector c1c0
                (
                    allBoundary[f1].centre(allBoundary.points())
                  - allBoundary[f0].centre(allBoundary.points())
                );

                if (concaveMultiCells && (c1c0 & n0) > SMALL)
                {
                    // Found concave edge. Make into multiCell features
                    Info<< "Detected concave feature edge:" << edgeI
                        << " cos:" << (c1c0 & n0)
                        << " coords:"
                        << allBoundary.points()[v0]
                        << allBoundary.points()[v1]
                        << endl;

                    singleCellFeaturePointSet.erase(v0);
                    multiCellFeaturePointSet.insert(v0);
                    singleCellFeaturePointSet.erase(v1);
                    multiCellFeaturePointSet.insert(v1);
                }
                else
                {
                    // Convex. singleCell feature.
                    if (!multiCellFeaturePointSet.found(v0))
                    {
                        singleCellFeaturePointSet.insert(v0);
                    }
                    if (!multiCellFeaturePointSet.found(v1))
                    {
                        singleCellFeaturePointSet.insert(v1);
                    }
                }
            }
        }
    }


    // 3. Mark all feature faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~


    // A. boundary faces.
    for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); facei++)
    {
        featureFaceSet.insert(facei);
    }

    // B. face zones.
    const faceZoneMesh& faceZones = mesh.faceZones();

    if (doNotPreserveFaceZones)
    {
        if (faceZones.size() > 0)
        {
            WarningInFunction
                << "Detected " << faceZones.size()
                << " faceZones. These will not be preserved."
                << endl;
        }
    }
    else
    {
        if (faceZones.size() > 0)
        {
            Info<< "Detected " << faceZones.size()
                << " faceZones. Preserving these by marking their"
                << " points, edges and faces as features." << endl;
        }

        forAll(faceZones, zoneI)
        {
            const faceZone& fz = faceZones[zoneI];

            Info<< "Inserting all faces in faceZone " << fz.name()
                << " as features." << endl;

            forAll(fz, i)
            {
                label facei = fz[i];
                const face& f = mesh.faces()[facei];
                const labelList& fEdges = mesh.faceEdges()[facei];

                featureFaceSet.insert(facei);
                forAll(f, fp)
                {
                    // Mark point as multi cell point (since both sides of
                    // face should have different cells)
                    singleCellFeaturePointSet.erase(f[fp]);
                    multiCellFeaturePointSet.insert(f[fp]);

                    // Make sure there are points on the edges.
                    featureEdgeSet.insert(fEdges[fp]);
                }
            }
        }
    }

    // Transfer to arguments
    featureFaces = featureFaceSet.toc();
    featureEdges = featureEdgeSet.toc();
    singleCellFeaturePoints = singleCellFeaturePointSet.toc();
    multiCellFeaturePoints = multiCellFeaturePointSet.toc();
}


// Dump features to .obj files
void dumpFeatures
(
    const polyMesh& mesh,
    const labelList& featureFaces,
    const labelList& featureEdges,
    const labelList& singleCellFeaturePoints,
    const labelList& multiCellFeaturePoints
)
{
    {
        OFstream str("featureFaces.obj");
        Info<< "Dumping centres of featureFaces to obj file " << str.name()
            << endl;
        forAll(featureFaces, i)
        {
            meshTools::writeOBJ(str, mesh.faceCentres()[featureFaces[i]]);
        }
    }
    {
        OFstream str("featureEdges.obj");
        Info<< "Dumping featureEdges to obj file " << str.name() << endl;
        label vertI = 0;

        forAll(featureEdges, i)
        {
            const edge& e = mesh.edges()[featureEdges[i]];
            meshTools::writeOBJ(str, mesh.points()[e[0]]);
            vertI++;
            meshTools::writeOBJ(str, mesh.points()[e[1]]);
            vertI++;
            str<< "l " << vertI-1 << ' ' << vertI << nl;
        }
    }
    {
        OFstream str("singleCellFeaturePoints.obj");
        Info<< "Dumping featurePoints that become a single cell to obj file "
            << str.name() << endl;
        forAll(singleCellFeaturePoints, i)
        {
            meshTools::writeOBJ(str, mesh.points()[singleCellFeaturePoints[i]]);
        }
    }
    {
        OFstream str("multiCellFeaturePoints.obj");
        Info<< "Dumping featurePoints that become multiple cells to obj file "
            << str.name() << endl;
        forAll(multiCellFeaturePoints, i)
        {
            meshTools::writeOBJ(str, mesh.points()[multiCellFeaturePoints[i]]);
        }
    }
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Creates the dual of a polyMesh,"
        " adhering to all the feature and patch edges."
    );

    #include "addOverwriteOption.H"
    argList::noParallel();

    argList::addArgument
    (
        "featureAngle",
        "in degrees [0-180]"
    );

    argList::addBoolOption
    (
        "splitAllFaces",
        "Have multiple faces in between cells"
    );
    argList::addBoolOption
    (
        "concaveMultiCells",
        "Split cells on concave boundary edges into multiple cells"
    );
    argList::addBoolOption
    (
        "doNotPreserveFaceZones",
        "Disable the default behaviour of preserving faceZones by having"
        " multiple faces in between cells"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    const word oldInstance = mesh.pointsInstance();

    // Mark boundary edges and points.
    // (Note: in 1.4.2 we can use the built-in mesh point ordering
    //  facility instead)
    bitSet isBoundaryEdge(mesh.nEdges());
    for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); facei++)
    {
        const labelList& fEdges = mesh.faceEdges()[facei];

        forAll(fEdges, i)
        {
            isBoundaryEdge.set(fEdges[i]);
        }
    }

    const scalar featureAngle = args.get<scalar>(1);
    const scalar minCos = Foam::cos(degToRad(featureAngle));

    Info<< "Feature:" << featureAngle << endl
        << "minCos :" << minCos << endl
        << endl;


    const bool splitAllFaces = args.found("splitAllFaces");
    if (splitAllFaces)
    {
        Info<< "Splitting all internal faces to create multiple faces"
            << " between two cells." << nl
            << endl;
    }

    const bool overwrite = args.found("overwrite");
    const bool doNotPreserveFaceZones = args.found("doNotPreserveFaceZones");
    const bool concaveMultiCells = args.found("concaveMultiCells");
    if (concaveMultiCells)
    {
        Info<< "Generating multiple cells for points on concave feature edges."
            << nl << endl;
    }


    // Face(centre)s that need inclusion in the dual mesh
    labelList featureFaces;
    // Edge(centre)s  ,,
    labelList featureEdges;
    // Points (that become a single cell) that need inclusion in the dual mesh
    labelList singleCellFeaturePoints;
    // Points (that become a multiple cells)        ,,
    labelList multiCellFeaturePoints;

    // Sample implementation of feature detection.
    simpleMarkFeatures
    (
        mesh,
        mesh,
        isBoundaryEdge,
        featureAngle,
        concaveMultiCells,
        doNotPreserveFaceZones,

        featureFaces,
        featureEdges,
        singleCellFeaturePoints,
        multiCellFeaturePoints
    );

    // If we want to split all polyMesh faces into one dualface per cell
    // we are passing through we also need a point
    // at the polyMesh facecentre and edgemid of the faces we want to
    // split.
    if (splitAllFaces)
    {
        featureEdges = identity(mesh.nEdges());
        featureFaces = identity(mesh.nFaces());
    }

    // Write obj files for debugging
    if(dumpFeaturesFlag){
        dumpFeatures
        (
            mesh,
            featureFaces,
            featureEdges,
            singleCellFeaturePoints,
            multiCellFeaturePoints
        );
    }
    



    // Read objects in time directory
    IOobjectList objects(mesh, runTime.timeName());

    // Read vol fields.
    PtrList<volScalarField> vsFlds;
    ReadFields(mesh, objects, vsFlds);

    PtrList<volVectorField> vvFlds;
    ReadFields(mesh, objects, vvFlds);

    PtrList<volSphericalTensorField> vstFlds;
    ReadFields(mesh, objects, vstFlds);

    PtrList<volSymmTensorField> vsymtFlds;
    ReadFields(mesh, objects, vsymtFlds);

    PtrList<volTensorField> vtFlds;
    ReadFields(mesh, objects, vtFlds);

    // Read surface fields.
    PtrList<surfaceScalarField> ssFlds;
    ReadFields(mesh, objects, ssFlds);

    PtrList<surfaceVectorField> svFlds;
    ReadFields(mesh, objects, svFlds);

    PtrList<surfaceSphericalTensorField> sstFlds;
    ReadFields(mesh, objects, sstFlds);

    PtrList<surfaceSymmTensorField> ssymtFlds;
    ReadFields(mesh, objects, ssymtFlds);

    PtrList<surfaceTensorField> stFlds;
    ReadFields(mesh, objects, stFlds);


    // Topo change container
    polyTopoChange meshMod(mesh.boundaryMesh().size());

    // Mesh dualiser engine
    meshDualiser dualMaker(mesh);

    // Insert all commands into polyTopoChange to create dual of mesh. This does
    // all the hard work.
    dualMaker.setRefinement
    (
        splitAllFaces,
        featureFaces,
        featureEdges,
        singleCellFeaturePoints,
        multiCellFeaturePoints,
        meshMod
    );

    // Create mesh, return map from old to new mesh.
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);

    // Update fields
    mesh.updateMesh(map());

    // Optionally inflate mesh
    if (map().hasMotionPoints())
    {
        mesh.movePoints(map().preMotionPoints());
    }

    if (!overwrite)
    {
        ++runTime;
    }
    else
    {
        mesh.setInstance(oldInstance);
    }

    Info<< "Writing dual mesh to " << runTime.timeName() << endl;

    mesh.write();
    topoSet::removeFiles(mesh);
    processorMeshes::removeFiles(mesh);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
