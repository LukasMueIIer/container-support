import os
from copy import deepcopy
from collections import defaultdict

from numpy import zeros, unique, where, argsort, searchsorted, allclose, array
import numpy as np
from scipy.spatial import cKDTree

from pyNastran.converters.aflr.ugrid.ugrid_reader import read_ugrid
from pyNastran.converters.aflr.surf.surf_reader import TagReader
import time

def convertMesh(ugrid,boundarySurf, boundaryTagFilePath, outputDirectory,cleanupInternalBoundaries = True, untaggedSinglesToBoundary = True,deepSearchUntaggedBoundaries = True):
    pointsFilename = os.path.join(outputDirectory,"points")
    facesFilename = os.path.join(outputDirectory,"faces")
    neighbourFilename = os.path.join(outputDirectory,"neighbour")
    ownerFilename = os.path.join(outputDirectory,"owner")
    boundaryFilename = os.path.join(outputDirectory,"boundary")

    _write_points(ugrid, pointsFilename)
    _write_faces(ugrid,boundarySurf, boundaryTagFilePath, facesFilename,neighbourFilename,ownerFilename,boundaryFilename, cleanupInternalBoundaries,untaggedSinglesToBoundary=untaggedSinglesToBoundary,deepSearchUntaggedBoundaries=True)

def _write_points(ugrid, points_filename):
    """writes an OpenFOAM points file"""
    with open(points_filename, 'w') as points_file:
        nnodes = ugrid.nodes.shape[0]

        points_file.write(
            '/*--------------------------------*- C++ -*----------------------------------*\\\n'
            '| =========                 |                                                 |\n'
            '| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n'
            '|  \\\\    /   O peration     | Version:  1.7.1                                 |\n'
            '|   \\\\  /    A nd           | Web:      www.OpenFOAM.com                      |\n'
            '|    \\\\/     M anipulation  |                                                 |\n'
            '\\*---------------------------------------------------------------------------*/\n'
            'FoamFile\n'
            '{\n'
            '    version     2.0;\n'
            '    format      ascii;\n'
            '    class       vectorField;\n'
            '    location    "constant/polyMesh";\n'
            '    object      points;\n'
            '}\n'
            '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * /\n'
        )


        points_file.write('\n\n')
        points_file.write('%i\n' % nnodes)
        points_file.write('(\n')
        for node in enumerate(ugrid.nodes):
            points_file.write('    (%-12s %-12s %-12s)\n' % (node[1][0], node[1][1], node[1][2]))
        points_file.write(')\n')


class polyWrap:
    def __init__(self,cellCount):
        #we know we have a maximum of 4 nodes
        self.indexCount = 0
        self.nodes = -2*np.ones((cellCount,4), dtype=int)
        self.sortedNodes = -1*np.ones((cellCount,4), dtype=int)
        self.nNodes = -1*np.ones((cellCount,1), dtype=int)
        self.owner = -1*np.ones((cellCount,1), dtype=int)
        self.neighbour = -1*np.ones((cellCount,1), dtype=int)
        self.boundaryId = -1*np.ones((cellCount,1), dtype=int)
        self.boundaryType = -1*np.ones((cellCount,1), dtype=int)

        self.valid = np.ones((cellCount, 1), dtype=int)

        self.last_print_time = time.time()  # Initialize timer
        self.faceMap = {} #map for linear match search time
    
    def addFace(self,newNodes,ownerCellId,boundaryId = -1, boundaryType = -1):

        if len(newNodes) != len(set(newNodes)):
            print(f"Error: Face has duplicate nodes, skipping: {newNodes}")
        
        #print(self.indexCount)
        width = 4    
        # 2. Create a temporary array filled with -1s
        nodes_padded = -1 * np.ones(width, dtype=int)
        # 3. Copy the input nodes into the temporary array
        # This handles lists like [3, 14, 98783] becoming [3, 14, 98783, -1]
        nodes_padded[:len(newNodes)] = newNodes

        if sum(nodes_padded) == -1:
            print("Error Empty Face")

        self.nodes[self.indexCount] = nodes_padded
        self.sortedNodes[self.indexCount] = np.sort(nodes_padded)
        self.nNodes[self.indexCount] = len(newNodes)
        self.owner[self.indexCount] = ownerCellId
        self.boundaryId[self.indexCount] = boundaryId
        self.boundaryType[self.indexCount] = boundaryType

        key = tuple(np.sort(nodes_padded))
        self.faceMap[key] = self.indexCount
        # -----------------------------------------------

        if True:    # --- Timer Snippet ---
            current_time = time.time()
            # Check if 5 seconds have passed since the last print
            if current_time - self.last_print_time >= 5.0:
                print(f"Current Index Count: {self.indexCount}")
                self.last_print_time = current_time  # Reset the timer
        
        self.indexCount = self.indexCount + 1

    def getBoundaryPatchCount(self):
        """
        Returns the maximum Boundary ID found among VALID faces only.
        Returns -1 if no boundary faces exist.
        """
        # 1. Slice arrays
        valid_b_ids = self.boundaryId[:self.indexCount, 0]
        valid_mask = self.valid[:self.indexCount, 0]

        # 2. Filter: Get boundary IDs where valid == 1
        active_boundaries = valid_b_ids[valid_mask == 1]

        # 3. Return max or -1 if empty
        if active_boundaries.size > 0:
            return np.max(active_boundaries)
        else:
            return -1


    def isBoundary(self,index):
        if(self.boundaryId[index] == -1):
            return False
        else:
            return True
    
    def isExternal(self,index):
        if(self.neighbour[index] == -1):
            return True
        else:
            return False
    
    def checkIfFaceExists(self,nodes):
        # 1. Pad and sort (same as before) to create a unique signature
        width = 4
        nodes_padded = -1 * np.ones(width, dtype=int)
        nodes_padded[:len(nodes)] = nodes
        
        # 2. Convert to a Tuple (dictionaries require immutable keys)
        # We sort it so [3, 1, 2] and [1, 2, 3] generate the same key
        key = tuple(np.sort(nodes_padded))
        
        # 3. Check the dictionary (Instant lookup)
        if key in self.faceMap:
            return True, self.faceMap[key]
        else:
            return False, -1
        
    def setNeighbour(self,nodeID,neighbourCellId):
        self.neighbour[nodeID] = neighbourCellId

    def writeFaceNodesWithBoundaryId(self,fileHandle,targetBoundaryId):
        # 1. Slice arrays
        valid_b_ids = self.boundaryId[:self.indexCount, 0]
        valid_mask = self.valid[:self.indexCount, 0]

        # 2. Create Boolean Mask: Matches Boundary ID AND is Valid
        mask = (valid_b_ids == targetBoundaryId) & (valid_mask == 1)

        # 3. Get indices
        matching_indices = np.where(mask)[0]

        # 4. Write loop
        count = 0
        for idx in matching_indices:
            # Get raw nodes and filter out padding (-1)
            raw_nodes = self.nodes[idx]
            real_nodes = raw_nodes[raw_nodes >= 0]
            
            # Format: count(node1 node2 ...)
            n_count = len(real_nodes)
            nodes_str = " ".join(map(str, real_nodes))
            
            line = f"{n_count}({nodes_str})\n"
            fileHandle.write(line)
            count += 1
            
        return count
    
    def writeFaceOwnerWithBoundaryId(self,fileHandle,targetBoundaryId):
        # 1. Slice arrays
        valid_b_ids = self.boundaryId[:self.indexCount, 0]
        valid_mask = self.valid[:self.indexCount, 0]

        # 2. Create Boolean Mask: Matches Boundary ID AND is Valid
        mask = (valid_b_ids == targetBoundaryId) & (valid_mask == 1)

        # 3. Get indices
        matching_indices = np.where(mask)[0]

        # 4. Write loop
        count = 0
        for idx in matching_indices:
            owner_id = self.owner[idx, 0]
            fileHandle.write(f"{owner_id}\n")
            count += 1
            
        return count
    
    def writeFaceNeighbourWithBoundaryId(self, fileHandle, targetBoundaryId):
        # 1. Slice arrays
        valid_b_ids = self.boundaryId[:self.indexCount, 0]
        valid_mask = self.valid[:self.indexCount, 0]

        # 2. Create Boolean Mask: Matches Boundary ID AND is Valid
        mask = (valid_b_ids == targetBoundaryId) & (valid_mask == 1)

        # 3. Get indices
        matching_indices = np.where(mask)[0]

        # 4. Write loop
        count = 0
        for idx in matching_indices:
            neighbour_id = self.neighbour[idx, 0]
            fileHandle.write(f"{neighbour_id}\n")
            count += 1
            
        return count

    
    def writeOnlyInternalFaces(self,fileHandle):
        return self.writeFaceNodesWithBoundaryId(fileHandle,-1)
    
    def debugMissingNeighbours(self, nodeCoords, targetFilePath):
        """
        Finds faces that are marked as internal (boundaryId == -1) 
        but have no neighbour (neighbour == -1).
        Writes 3D points to CSV including the face node count (3 or 4).
        """
        # 1. Slice valid data ranges
        b_ids = self.boundaryId[:self.indexCount, 0]
        n_ids = self.neighbour[:self.indexCount, 0]

        # 2. Create Mask: Internal faces (-1) with missing neighbour (-1)
        mask = (b_ids == -1) & (n_ids == -1)
        
        problem_indices = np.where(mask)[0]
        count = len(problem_indices)
        
        if count == 0:
            print("Debug: No disconnected internal faces found. Mesh topology looks safe.")
            return 0

        print(f"Debug: Found {count} faces with boundaryId -1 but missing neighbour! Writing to {targetFilePath}...")

        # 3. Write to CSV
        with open(targetFilePath, 'w') as f:
            # Header update: added nodeCount
            f.write("x,y,z,faceIndex,nodeIndex,nodeCount\n")
            
            for face_idx in problem_indices:
                # Get the nodes for this face
                raw_nodes = self.nodes[face_idx]
                
                # Filter out -1 padding
                real_nodes = raw_nodes[raw_nodes >= 0]
                
                # Get the total count (e.g., 3 for Tri, 4 for Quad)
                n_count = len(real_nodes)
                
                # Look up coordinates and write
                for node_idx in real_nodes:
                    x, y, z = nodeCoords[node_idx]
                    
                    # Write line with the new column at the end
                    f.write(f"{x},{y},{z},{face_idx},{node_idx},{n_count}\n")
                    
        return count
    
    def countFacesWithNeighbour(self):
        """
        Counts how many faces have a valid neighbour (neighbour != -1).
        These are typically strictly internal faces.
        """
        # 1. Slice valid data
        valid_neighbours = self.neighbour[:self.indexCount, 0]
        
        # 2. Create Mask: True where neighbour is NOT -1
        mask = (valid_neighbours != -1)
        
        # 3. Sum the True values (True = 1, False = 0)
        count = np.sum(mask)
        
        return count
    
    def cleanUpInternalBoundaries(self):
        """
        Scans all faces. If a face has a valid neighbour (neighbour != -1)
        but is currently marked as a boundary (boundaryId != -1),
        it resets the boundaryId to -1 (internal).
        """
        # 1. Slice the arrays to the current valid count
        # These are views into the main array, so they share memory
        current_neighbours = self.neighbour[:self.indexCount]
        current_boundaries = self.boundaryId[:self.indexCount]

        # 2. Create a Boolean Mask
        # Condition: Has a neighbour AND is marked as a boundary
        mask = (current_neighbours != -1) & (current_boundaries != -1)

        # 3. Count how many we are about to fix (for reporting)
        fixed_count = np.sum(mask)

        # 4. Apply the fix
        # We set boundaryId to -1 wherever the mask is True
        current_boundaries[mask] = -1

        print(f"Cleanup: Fixed {fixed_count} faces that had neighbours but were marked as boundaries.")
        return fixed_count
    
    def untagedBoundaryToPatch(self):
        """
        Finds faces that have no neighbour (-1) and no assigned boundary ID (-1).
        Assigns them a new boundary ID that is incremented from the current maximum.
        """
        # 1. Determine the new ID
        # getBoundaryPatchCount() returns the max ID. 
        # If no boundaries exist, it returns -1. Next ID becomes 0.
        # If max is 5, Next ID becomes 6.
        current_max_id = self.getBoundaryPatchCount()
        new_patch_id = current_max_id + 1
        
        # 2. Slice arrays to valid range
        current_neighbours = self.neighbour[:self.indexCount]
        current_boundaries = self.boundaryId[:self.indexCount]
        
        # 3. Create Mask
        # Condition: External (neighbour == -1) AND Untagged (boundaryId == -1)
        mask = (current_neighbours == -1) & (current_boundaries == -1)
        
        # 4. Apply new ID
        count = np.sum(mask)
        
        if count > 0:
            current_boundaries[mask] = new_patch_id
            print(f"Auto-Patch: Assigned new boundary ID {new_patch_id} to {count} untagged external faces.")
        else:
            print("Auto-Patch: No untagged external faces found.")
            
        return new_patch_id, count
    
    def getUntaggedBoundaryFaceIndeces(self):
        #returns the indeces of all untagged boundary faces
        # 1. Slice valid data to current count
        # shape is (N, 1), so we access column 0 to get 1D arrays
        valid_neighbours = self.neighbour[:self.indexCount, 0]
        valid_boundaries = self.boundaryId[:self.indexCount, 0]

        # 2. Create Boolean Mask
        mask = (valid_neighbours == -1) & (valid_boundaries == -1)

        # 3. Return the indices where the mask is True
        return np.where(mask)[0]
    
    def getValidFaceCount(self):
        """
        Returns the total number of faces marked as valid (1).
        """
        # Slice to valid range and sum the boolean/int column
        # valid is shape (N, 1), so we access [:, 0]
        return np.sum(self.valid[:self.indexCount, 0])
        
    def getSingleBoundaryIndex(self):
        """
        Returns the indices of all faces that:
        1. Do not have a neighbour (neighbour == -1) -> External/Boundary
        2. Do not have a boundary tag (boundaryId == -1) -> Untagged
        """
        # 1. Slice valid data to current count to avoid checking empty memory
        # self.neighbour and self.boundaryId are (N, 1), so we access column 0
        valid_neighbours = self.neighbour[:self.indexCount, 0]
        valid_boundaries = self.boundaryId[:self.indexCount, 0]

        # 2. Create Boolean Mask
        # Condition: External (no neighbour) AND Untagged (no boundary ID)
        mask = (valid_neighbours == -1) & (valid_boundaries == -1)

        # 3. Return the indices where the mask is True
        return np.where(mask)[0]
    
    def getQuadFaceIndices(self):
        """
        Returns an array of indices for all faces that have exactly 4 valid nodes.
        Assumes -1 is used as the padding value for invalid nodes.
        """
        # Count how many nodes are >= 0 in each row (axis 1)
        valid_node_counts = np.sum(self.nodes >= 0, axis=1)
        
        # Return indices where the count is exactly 4
        return np.where(valid_node_counts == 4)[0]
    
    def expandTo8Nodes(self):
        """
        Expands the internal nodes array to have 8 columns.
        Existing data is preserved in the first columns.
        New columns are initialized with -1.
        """
        rows, cols = self.nodes.shape
        
        if cols >= 8:
            return  # Already large enough
            
        # 1. Create a new array of shape (N, 8) filled with -1
        # Preserves the data type of the original array (likely int)
        new_nodes = np.full((rows, 8), -1, dtype=self.nodes.dtype)
        
        # 2. Copy the existing data into the left side
        new_nodes[:, :cols] = self.nodes
        
        # 3. Replace the class attribute
        self.nodes = new_nodes
        print(f"Expanded node storage from {cols} to 8 columns.")

def createSurfaceToVolMapKD(A,B):

    # 1. Build a Tree from List B (The source)
    tree = cKDTree(B)

    # 2. Query the Tree with List A
    # k=1 finds the single closest neighbor.
    # distances: How far apart the match is (should be ~0)
    # map_A_to_B: The indices in B corresponding to items in A
    distances, map_A_to_B = tree.query(A, k=1)

    # 3. Usage
    print("Indices in B corresponding to A:", map_A_to_B)

    # Optional: Verify the match is actually close enough
    # (Filters out 'nearest neighbors' that aren't actually the same point)
    tolerance = 1e-5
    valid_matches = distances < tolerance

    if not np.all(valid_matches):
        print("Warning: Some items in A did not find an exact match in B within tolerance.")

    return map_A_to_B

def createSurfaceToVolMap(A,B):
    import numpy as np


    # Decide on precision (e.g., 5 decimal places)
    PRECISION = 5

    # 1. Create Lookup for B with ROUNDED keys
    # np.round ensures floats are truncated safely before hashing
    b_lookup = {tuple(np.round(arr, PRECISION)): i for i, arr in enumerate(B)}

    # 2. Create the Map
    # We must also round A when looking up
    map_A_to_B = []
    for arr in A:
        key = tuple(np.round(arr, PRECISION))
        if key in b_lookup:
            map_A_to_B.append(b_lookup[key])
        else:
            # Handle case where no match is found
            map_A_to_B.append(-1) 

    print("Mapping List:", map_A_to_B)
    print("Debug")

def save_points_for_paraview(A, point_dump_path):
    """
    Saves a list of 3D points (A) to a CSV file compatible with ParaView.
    
    Args:
        A: List of arrays or (N, 3) numpy array representing points.
        point_dump_path: Full file path (e.g., 'outputs/points.csv')
    """
    # 1. Convert input list to a standardized (N, 3) numpy matrix
    data_matrix = np.array(A)

    # 2. Safety check: Ensure the output directory exists
    directory = os.path.dirname(point_dump_path)
    if directory:
        os.makedirs(directory, exist_ok=True)

    # 3. Write to CSV
    # header="x,y,z" -> Names the columns so ParaView can find them easily.
    # comments=""    -> Important! Removes the '#' that numpy usually adds to headers.
    # fmt="%.8f"     -> Uses 8 decimal places (prevents scientific notation for smaller floats).
    np.savetxt(
        point_dump_path, 
        data_matrix, 
        delimiter=",", 
        header="x,y,z", 
        comments="", 
        fmt="%.8f"
    )
    
    print(f"Successfully wrote {len(data_matrix)} points to {point_dump_path}")

def _write_faces(ugrid,boundarySurf, boundaryTagFilePath, facesFilename,neighbourFilename,ownerFilename,boundaryFilename,removeInternalBoundaries=True,untaggedSinglesToBoundary = True,deepSearchUntaggedBoundaries=True, rematchIceInterface = True):
    # HERE WE SET UP THE MAPPING FROM THE BOUNDARY INFORMATIONS WE HAVE
    print("Importing the boundary")
    #transfering the boundaryInformation
    boundaryTris = boundarySurf.tris
    boundaryTrisProps  = boundarySurf.tri_props
    sortedTrisMap = {}
    i = 0

    nodeMap = createSurfaceToVolMapKD(boundarySurf.nodes,ugrid.nodes)

    for tris in boundaryTris:
        # 1. Get the Surface Node IDs (e.g., [1, 2, 3])
        # 2. Use nodeMap to translate them to Volume Node IDs (e.g., [105, 99, 500])
        #    Note: We perform -1 if boundarySurf is 1-based, but pyNastran usually makes it 0-based.
        #    If your boundarySurf is 1-based, use: nodeMap[tris - 1]
        vol_nodes = nodeMap[tris-1] 
        
        # 3. Sort and create key using VOLUME IDs
        key = tuple(np.sort(vol_nodes))
        sortedTrisMap[key] = i
        i = i + 1

    boundaryQuads = boundarySurf.quads
    boundaryQuadsProps = boundarySurf.quad_props
    sortedQuadsMap = {}

    i = 0
    for quads in boundaryQuads:
        # Translate Surface IDs -> Volume IDs
        vol_nodes = nodeMap[quads-1]
        
        key = tuple(np.sort(vol_nodes))
        sortedQuadsMap[key] = i
        i = i + 1

    #Create the map for the boundary contained within the ugrid
    #the surface ID map is joint from tris and quads
    volSurfaceTrisMap = {}
    i = 0
    for tris in ugrid.tris:
        # 3. Sort and create key using VOLUME IDs
        key = tuple(np.sort(tris-1))
        volSurfaceTrisMap[key] = i
        i = i + 1

    volSurfaceTrisCount = i

    volSurfaceQuadMap = {}
    i = 0
    for quad in ugrid.quads:
        # 3. Sort and create key using VOLUME IDs
        key = tuple(np.sort(quad-1))
        volSurfaceQuadMap[key] = i
        i = i + 1

    

    

    #Create connectivity information for the cells
    nhexas = ugrid.hexas.shape[0]
    npenta6s = ugrid.penta6s.shape[0]
    npenta5s = ugrid.penta5s.shape[0]
    ntets = ugrid.tets.shape[0]

    nquad_faces = nhexas * 6 + npenta5s + npenta6s * 3
    ntri_faces = ntets * 4 + npenta5s * 4 + npenta6s * 2
    nfaces = ntri_faces + nquad_faces
    assert nfaces > 0, nfaces

    polyFaces = polyWrap(nfaces)

    def treatFace(nodes,cellID):

        matchIndex = -1
        match = False

        match, matchIndex = polyFaces.checkIfFaceExists(nodes)

        #if the cell does not exist yet we are the owner, check if we are an internal boundary
        #if we are boundary, set the boundary condition, and our surface id
        if not match:
            
            boundaryId = -1
            boundaryType = -1
            
            #Check if we are a boundary face
            if len(nodes) == 3:

                key = tuple(np.sort(nodes))

                if key in sortedTrisMap:

                    found_index = sortedTrisMap[key]
            
                    boundaryId = boundaryTrisProps[found_index][0]-1
                    if(boundaryId == -1):
                        print("Empty Entry in Boundary Field")
                    boundaryType = boundaryTrisProps[found_index][2]

                #Now check if we are in the volSurf and if so overwrite the boundaryID

                if key in volSurfaceTrisMap:
                    found_index = volSurfaceTrisMap[key]

                    boundaryId = ugrid.pids[found_index]-1
                       
            elif len(nodes) == 4:

                key = tuple(np.sort(nodes))

                if key in sortedQuadsMap:

                    found_index = sortedQuadsMap[key]
            
                    boundaryId = boundaryQuadsProps[found_index][0]
                    boundaryType = boundaryQuadsProps[found_index][2]

                if key in volSurfaceQuadMap:
                    found_index = volSurfaceQuadMap[key]

                    boundaryId = ugrid.pids[found_index+volSurfaceTrisCount]-1
        
            else:
                print("Errorrrrrr face is invalid")

            polyFaces.addFace(nodes,cellID,boundaryId,boundaryType)

            
        if match: #if it does we are neighbour and it is not an external cell, but is a boudary for sure
            #check for double matches
            if(polyFaces.neighbour[matchIndex]==-1):
                polyFaces.setNeighbour(matchIndex,cellID)
            else:
                print("Tripple appreaance")
       

    cellID = 0

    # =========================================================
    # HELPER FUNCTIONS FOR GEOMETRIC CHECKS
    # =========================================================
    
    # Pre-fetch coordinates to speed up lookups
    coords = ugrid.nodes 

    def get_face_normal(p1, p2, p3):
        vecA = p2 - p1
        vecB = p3 - p1
        return np.cross(vecA, vecB)

    def ensure_outward(face_nodes, cell_center, cell_id):
        """
        Checks if face normal points away from cell center.
        Handles both Triangles (3 nodes) and Quads (4 nodes).
        """
        # Get coordinates
        p1 = coords[face_nodes[0]]
        p2 = coords[face_nodes[1]]
        p3 = coords[face_nodes[2]]

        # Calculate Face Center
        if len(face_nodes) == 3:
            face_center = (p1 + p2 + p3) / 3.0
        elif len(face_nodes) == 4:
            p4 = coords[face_nodes[3]]
            face_center = (p1 + p2 + p3 + p4) / 4.0
        
        # Calculate Normal (using first 3 points is sufficient for direction)
        normal = get_face_normal(p1, p2, p3)

        # Vector from Cell Center to Face Center
        center_to_face = face_center - cell_center

        # Check direction
        if np.dot(normal, center_to_face) < 0:
            print(f"Notice: Flipped face {len(face_nodes)}-nodes for Cell {cell_id} to ensure outward normal.")
            
            # Flip winding to reverse normal
            if len(face_nodes) == 3:
                return [face_nodes[0], face_nodes[2], face_nodes[1]]
            else:
                # For Quad [0,1,2,3], swapping 1 and 3 reverses winding -> [0,3,2,1]
                return [face_nodes[0], face_nodes[3], face_nodes[2], face_nodes[1]]
        
        return face_nodes

    # =========================================================s
    # UPDATED CELL LOOPS
    # =========================================================

    cellID = 0

    checkOutward = True

    # --- TETS (Tetrahedra) ---
    print("Processing Tets...")
    for element in ugrid.tets - 1:
        (n1, n2, n3, n4) = element
        
        # Calculate Cell Center
        cell_center = (coords[n1] + coords[n2] + coords[n3] + coords[n4]) / 4.0

        # Define Faces
        face1 = [n1, n3, n2]
        face2 = [n2, n3, n4]
        face3 = [n1, n4, n3]
        face4 = [n1, n2, n4]

        # Check Orientation
        if checkOutward:
            face1 = ensure_outward(face1, cell_center, cellID)
            face2 = ensure_outward(face2, cell_center, cellID)
            face3 = ensure_outward(face3, cell_center, cellID)
            face4 = ensure_outward(face4, cell_center, cellID)

        treatFace(face1, cellID)
        treatFace(face2, cellID)
        treatFace(face3, cellID)
        treatFace(face4, cellID)

        cellID += 1
    
    print("Tet Cells Done, Moving to Penta5; current Index: " + str(cellID))

    # --- PENTA5 (Pyramids) ---
    for element in ugrid.penta5s - 1:
        (n1, n2, n3, n4, n5) = element

        # Calculate Cell Center
        cell_center = (coords[n1] + coords[n2] + coords[n3] + coords[n4] + coords[n5]) / 5.0

        # Define Faces (1 Quad, 4 Tris)
        face1 = [n2, n3, n5]
        face2 = [n4, n5, n3]
        face3 = [n1, n4, n3]
        face4 = [n1, n3, n2]
        face5 = [n1, n2, n5, n4] # Base Quad

        # Check Orientation
        if checkOutward:
            face1 = ensure_outward(face1, cell_center, cellID)
            face2 = ensure_outward(face2, cell_center, cellID)
            face3 = ensure_outward(face3, cell_center, cellID)
            face4 = ensure_outward(face4, cell_center, cellID)
            face5 = ensure_outward(face5, cell_center, cellID)

        treatFace(face1, cellID)
        treatFace(face2, cellID)
        treatFace(face3, cellID)
        treatFace(face4, cellID)
        treatFace(face5, cellID)

        cellID += 1

    print("Penta5 done moving to penta6; current Index: " + str(cellID))

    # --- PENTA6 (Prisms/Wedges) ---
    for element in ugrid.penta6s - 1:
        (n1, n2, n3, n4, n5, n6) = element

        # Calculate Cell Center
        cell_center = (coords[n1] + coords[n2] + coords[n3] + coords[n4] + coords[n5] + coords[n6]) / 6.0

        # Define Faces (2 Tris, 3 Quads)
        face1 = [n1, n4, n5, n2]       # Bottom Tri
        face2 = [n1, n3, n6, n4]       # Top Tri
        face3 = [n3, n2, n5, n6]   # Side Quad
        face4 = [n1, n2, n3]   # Side Quad
        face5 = [n4, n6, n5]   # Side Quad

        # Check Orientation
        if checkOutward:
            face1 = ensure_outward(face1, cell_center, cellID)
            face2 = ensure_outward(face2, cell_center, cellID)
            face3 = ensure_outward(face3, cell_center, cellID)
            face4 = ensure_outward(face4, cell_center, cellID)
            face5 = ensure_outward(face5, cell_center, cellID)

        treatFace(face1, cellID)
        treatFace(face2, cellID)
        treatFace(face3, cellID)
        treatFace(face4, cellID)
        treatFace(face5, cellID)

        cellID += 1

    print("Penta6 done; current Index: " + str(cellID))
    print("Moving to hexas; current Index: " + str(cellID))

    # --- HEXAS (Hexahedra) ---
    for element in ugrid.hexas - 1:
        (n1, n2, n3, n4, n5, n6, n7, n8) = element

        # Calculate Cell Center
        cell_center = (coords[n1] + coords[n2] + coords[n3] + coords[n4] + 
                       coords[n5] + coords[n6] + coords[n7] + coords[n8]) / 8.0

        # Define Faces (6 Quads)
        face1 = [n1, n4, n3, n2]
        face2 = [n1, n2, n6, n5]
        face3 = [n4, n8, n7, n3]
        face4 = [n5, n6, n7, n8]
        face5 = [n1, n5, n8, n4]
        face6 = [n2, n3, n7, n6]

        # Check Orientation
        if checkOutward:
            face1 = ensure_outward(face1, cell_center, cellID)
            face2 = ensure_outward(face2, cell_center, cellID)
            face3 = ensure_outward(face3, cell_center, cellID)
            face4 = ensure_outward(face4, cell_center, cellID)
            face5 = ensure_outward(face5, cell_center, cellID)
            face6 = ensure_outward(face6, cell_center, cellID)

        treatFace(face1, cellID)
        treatFace(face2, cellID)
        treatFace(face3, cellID) 
        treatFace(face4, cellID)
        treatFace(face5, cellID)
        treatFace(face6, cellID)

        cellID += 1

    currentCounter = 0
    boundaryCellArray = []
    boundaryCellCounter = []
    # for each boundary Surface index this holds the Starting cell ID 


    #Perform cleanup of the conversion operations
    if removeInternalBoundaries:
        polyFaces.cleanUpInternalBoundaries()

    if deepSearchUntaggedBoundaries:
        untagged = polyFaces.getUntaggedBoundaryFaceIndeces()

        # Assemble the node list for each boundaryID 
        maxBoundaryID = np.max(ugrid.pids)
        trisPids = ugrid.pids[:volSurfaceTrisCount]
        quadPids = ugrid.pids[volSurfaceTrisCount:]
        
        print(f"Deep Search: Attempting to recover {len(untagged)} untagged faces by checking node containment...")

        

        for i in range(1, maxBoundaryID + 1):
            # 1. Filter: Find surface faces that belong to PID 'i'
            mask_tris = (trisPids == i)
            mask_quads = (quadPids == i)
            
            # 2. Collect Nodes: Get unique nodes from these faces
            # We use a Set for O(1) lookup speed.
            # Note: ugrid.tris/quads are 1-based, polyFaces are 0-based. We apply -1.
            boundaryNodeSet = set()
            
            if np.any(mask_tris):
                # np.unique flattens the array and finds unique elements
                boundaryNodeSet.update(np.unique(ugrid.tris[mask_tris]) - 1)
            
            if np.any(mask_quads):
                boundaryNodeSet.update(np.unique(ugrid.quads[mask_quads]) - 1)
            
            if not boundaryNodeSet:
                continue
            # 3. Match: Check untagged faces against this set
            recovered_count = 0
            for idx in untagged:
                # Skip if we already assigned this face in a previous loop iteration
                if polyFaces.boundaryId[idx] != -1:
                    continue
                
                # Get the specific nodes of the untagged face (ignoring -1 padding)
                raw_nodes = polyFaces.nodes[idx]
                face_nodes = raw_nodes[raw_nodes >= 0]
                
                # Check if ALL nodes of this face are present in the boundary set
                if boundaryNodeSet.issuperset(face_nodes):
                    polyFaces.boundaryId[idx] = i - 1
                    recovered_count += 1
            
            if recovered_count > 0:
                print(f"  -> Recovered {recovered_count} faces into Boundary ID {i-1}")
                


            #!warning we must make sure to do a -1 before assignment
    
    if rematchIceInterface:
        print("Rematching Ice Interface (Spatial Search)...")
        
        distTolerance = 1e-4  # Tolerance for "snapping" virtual nodes
        
        # 1. Gather all unique nodes currently on the interface
        iceFacesIndex = polyFaces.getSingleBoundaryIndex()
        
        unique_interface_node_ids = set()
        for face_idx in iceFacesIndex:
            raw_nodes = polyFaces.nodes[face_idx]
            unique_interface_node_ids.update(raw_nodes[raw_nodes >= 0])
            
        interface_node_list = list(unique_interface_node_ids)
        interface_coords = ugrid.nodes[interface_node_list]
        
        # 2. Build KDTree
        tree = cKDTree(interface_coords)
        print(f"Built search tree with {len(interface_node_list)} interface nodes.")

        # 3. Initialize Connectivity Map & Log
        interfaceNodeMembers = defaultdict(list)
        modification_log = [] 
        
        snapped_count = 0
        failed_count = 0

        # 4. Loop over faces

        class DefaultMap(dict):
            def __missing__(self, key):
                return -1

        edgeSplitingIndex = DefaultMap()
        #Logs which edges must be split 
        #Nodes of the edges to split are stored in sorted order and the return value is the node to be inserted
        #if the edge is not to be split will return -1

        for face_idx in iceFacesIndex:
            
            raw_nodes = polyFaces.nodes[face_idx]
            face_nodes = raw_nodes[raw_nodes >= 0]
            num_nodes = len(face_nodes)

            # --- Step A: Always add known corner nodes ---
            for node_id in face_nodes:
                interfaceNodeMembers[node_id].append(face_idx)

            # --- Step B: Handle Virtual Nodes for Quads ---
            if num_nodes == 4:
                # Get Coordinates
                p1 = ugrid.nodes[face_nodes[0]]
                p2 = ugrid.nodes[face_nodes[1]]
                p3 = ugrid.nodes[face_nodes[2]]
                p4 = ugrid.nodes[face_nodes[3]]
                pArr = [face_nodes[0],face_nodes[1],face_nodes[2],face_nodes[3],face_nodes[0]] #for later easy use
                
                # Calculate Virtual Points
                v5 = (p1 + p2) * 0.5
                v6 = (p2 + p3) * 0.5
                v7 = (p3 + p4) * 0.5
                v8 = (p4 + p1) * 0.5
                v9 = (p1 + p2 + p3 + p4) * 0.25
                
                virtual_points = [
                    (v5, 5), (v6, 6), (v7, 7), (v8, 8), (v9, 9)
                ]
                i = -1
                for v_point, inode_idx in virtual_points:
                    i = i + 1
                    # Query the KDTree for the nearest neighbor
                    dist, tree_idx = tree.query(v_point)
                    
                    if dist < distTolerance:
                        # Success
                        found_node_id = interface_node_list[tree_idx]
                        interfaceNodeMembers[found_node_id].append(face_idx)
                        
                        status = "SNAPPED" if dist > 1e-8 else "EXACT"
                        modification_log.append([face_idx, inode_idx, status, dist])
                        snapped_count += 1
                        
                        #Append information to edge splitting
                        if i<4: #Dont do the center point
                            corneNodes = [pArr[i],pArr[i+1]]
                            cornerKey = tuple(np.sort(corneNodes))
                            edgeSplitingIndex[cornerKey] = found_node_id

                    else:
                        # Failure
                        # print(f"Notification: Face {face_idx} Virtual Node {inode_idx} distance {dist:.2e} > tolerance.")
                        modification_log.append([face_idx, inode_idx, "FAILED", dist])
                        failed_count += 1

        # 5. Print Tabulated Log
        print("\n" + "="*60)
        print(f"{'Face ID':<10} | {'V-Node':<8} | {'Status':<10} | {'Error/Dist':<15}")
        print("-" * 60)
        
        print_count = 0
        for entry in modification_log:
            f_id, v_id, status, dist = entry
            
            # Print failures or the first 20 entries
            if status == "FAILED" or print_count < 20:
                print(f"{f_id:<10} | {v_id:<8} | {status:<10} | {dist:<15.2e}")
                print_count += 1
                
        if len(modification_log) > 20:
            print(f"... and {len(modification_log) - 20} more entries.")
        print("="*60)
        
        # 6. Print Statistics Line
        print(f"STATISTICS: Snapped/Exact: {snapped_count} | Failed/Missing: {failed_count}")
        print("="*60 + "\n")

        # Initialize statistics counters
        count_no_match = 0
        count_one_match = 0
        count_multi_match = 0

        print("Checking face connectivity against interface members...")

        for face_idx in iceFacesIndex:
            
            # Skip if this face was already invalidated (e.g., it was the "partner" of a previous match)
            if polyFaces.valid[face_idx, 0] == 0:
                continue

            # Get the nodes for this face
            raw_nodes = polyFaces.nodes[face_idx]
            face_nodes = raw_nodes[raw_nodes >= 0]
            
            face_sets = []
            possible_match = True
            
            # 1. Pull faceIDs for each node
            for node_id in face_nodes:
                if node_id in interfaceNodeMembers:
                    face_sets.append(set(interfaceNodeMembers[node_id]))
                else:
                    possible_match = False
                    break
            
            if not possible_match or not face_sets:
                polyFaces.valid[face_idx] = 0
                count_no_match += 1
                continue

            # 2. Intersection
            common_faces = face_sets[0]
            for s in face_sets[1:]:
                common_faces = common_faces.intersection(s)
            
            # 3. Remove self
            if face_idx in common_faces:
                common_faces.remove(face_idx)
            
            # 4. Check matches
            match_count = len(common_faces)
            
            if match_count == 0:
                # No match -> Invalid
                polyFaces.valid[face_idx] = 0
                count_no_match += 1
                
            elif match_count == 1:
                # Exactly one match -> Action A
                matched_face_id = list(common_faces)[0]
                
                # Retrieve owners (use [0] to access value from (N,1) array)
                owner_self = polyFaces.owner[face_idx, 0]
                owner_match = polyFaces.owner[matched_face_id, 0]
                
                # Determine new owner (smaller) and neighbour (larger)
                new_owner = min(owner_self, owner_match)
                new_neighbour = max(owner_self, owner_match)

                if new_owner != owner_self:
                    # Reverse the valid nodes
                    flipped_nodes = face_nodes[::-1]
                    
                    # Reconstruct the row with padding
                    # Use -1 * ones ensures we have the correct array shape/type
                    new_node_row = -1 * np.ones_like(raw_nodes)
                    new_node_row[:len(flipped_nodes)] = flipped_nodes
                    
                    # Write back to the main array
                    polyFaces.nodes[face_idx] = new_node_row
                
                # Update the CURRENT face
                polyFaces.owner[face_idx, 0] = new_owner
                polyFaces.neighbour[face_idx, 0] = new_neighbour
                
                # INVALIDATE the MATCHED face (so it isn't written out, preventing duplicates)
                #polyFaces.valid[matched_face_id] = 0
                
                count_one_match += 1
                
            else:
                # > 1 match -> Error
                print(f"Error: Face {face_idx} has multiple ({match_count}) shared face matches: {common_faces}")
                count_multi_match += 1

        #Splitting edges
        quadIndexes = polyFaces.getQuadFaceIndices()
        polyFaces.expandTo8Nodes()
        for index in quadIndexes:
            #pull the local nodes
            tempNodes = polyFaces.nodes[index]
            paddedNodes = [tempNodes[0],
                           tempNodes[1],
                           tempNodes[2],
                           tempNodes[3],
                           tempNodes[0]]
            splittedNodes = []
            for i in range(0,4):
                splittedNodes.append(paddedNodes[i])
                #check if we need to split the face
                key = tuple(np.sort([paddedNodes[i]  ,paddedNodes[i+1]]))
                centerNode = edgeSplitingIndex[key]
                if centerNode != -1:
                    #we need to split the edge
                    print("Splitting the edge: " + str(i) + " of Face: " + str(index))
                    splittedNodes.append(centerNode)
            #Write back the splitted Nodes information
            writeableNode = [-1,-1,-1,-1,
                             -1,-1,-1,-1]
            for j in range(0,len(splittedNodes)):
                writeableNode[j] = splittedNodes[j]

            polyFaces.nodes[index] = writeableNode

        # 5. Print Statistics
        print("\n" + "="*60)
        print("FACE CONNECTIVITY STATISTICS")
        print("-" * 60)
        print(f"No Match (Invalidated)   : {count_no_match}")
        print(f"One Match (Merged)       : {count_one_match}")
        print(f"Multi Match (Error)      : {count_multi_match}")
        print("="*60 + "\n")

    if untaggedSinglesToBoundary:
        polyFaces.untagedBoundaryToPatch()




    totalValidFaces = polyFaces.getValidFaceCount()

    ## No we write the faces file
    with open(facesFilename, 'w') as facesFile:
        facesFile.write(
            '/*--------------------------------*- C++ -*----------------------------------*\\\n'
            '| =========                 |                                                 |\n'
            '| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n'
            '|  \\\\    /   O peration     | Version:  1.7.1                                 |\n'
            '|   \\\\  /    A nd           | Web:      www.OpenFOAM.com                      |\n'
            '|    \\\\/     M anipulation  |                                                 |\n'
            '\\*---------------------------------------------------------------------------*/\n'
            'FoamFile\n'
            '{\n'
            '    version     2.0;\n'
            '    format      ascii;\n'
            '    class       faceList;\n'
            '    location    "constant/polyMesh";\n'
            '    object      faces;\n'
            '}\n'
            '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * /\n'
        )


        facesFile.write('\n\n')
        facesFile.write('%i\n' % totalValidFaces)
        facesFile.write('(\n')
        
        #write the internal faces and keep track of the count
        currentCounter = polyFaces.writeOnlyInternalFaces(facesFile)
        
        #iterate over and print the boundary faces
        for i in range(0,polyFaces.getBoundaryPatchCount()+1):
            boundaryCellArray.append(currentCounter)
            currentCounter = polyFaces.writeFaceNodesWithBoundaryId(facesFile,i)
            boundaryCellCounter.append(currentCounter)
            currentCounter = currentCounter+boundaryCellArray[-1]

        facesFile.write(')\n')

    with open(ownerFilename, 'w') as file:
        file.write(
            '/*--------------------------------*- C++ -*----------------------------------*\\\n'
            '| =========                 |                                                 |\n'
            '| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n'
            '|  \\\\    /   O peration     | Version:  1.7.1                                 |\n'
            '|   \\\\  /    A nd           | Web:      www.OpenFOAM.com                      |\n'
            '|    \\\\/     M anipulation  |                                                 |\n'
            '\\*---------------------------------------------------------------------------*/\n'
            'FoamFile\n'
            '{\n'
            '    version     2.0;\n'
            '    format      ascii;\n'
            '    class       labelList;\n'
            '    location    "constant/polyMesh";\n'
            '    object      owner;\n'
            '}\n'
            '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * /\n'
        )


        file.write('\n\n')
        file.write('%i\n' % totalValidFaces)
        file.write('(\n')
        
        #write the internal faces and keep track of the count
        currentCounter = polyFaces.writeFaceOwnerWithBoundaryId(file,-1)
        
        #iterate over and print the boundary faces
        for i in range(0,polyFaces.getBoundaryPatchCount()+1):
            currentCounter = polyFaces.writeFaceOwnerWithBoundaryId(file,i)

        file.write(')\n')

    with open(neighbourFilename, 'w') as file:
        file.write(
            '/*--------------------------------*- C++ -*----------------------------------*\\\n'
            '| =========                 |                                                 |\n'
            '| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n'
            '|  \\\\    /   O peration     | Version:  1.7.1                                 |\n'
            '|   \\\\  /    A nd           | Web:      www.OpenFOAM.com                      |\n'
            '|    \\\\/     M anipulation  |                                                 |\n'
            '\\*---------------------------------------------------------------------------*/\n'
            'FoamFile\n'
            '{\n'
            '    version     2.0;\n'
            '    format      ascii;\n'
            '    class       labelList;\n'
            '    location    "constant/polyMesh";\n'
            '    object      neighbour;\n'
            '}\n'
            '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * /\n'
        )


        file.write('\n\n')
        file.write('%i\n' % totalValidFaces)
        file.write('(\n')
        
        #write the internal faces and keep track of the count
        currentCounter = polyFaces.writeFaceNeighbourWithBoundaryId(file,-1)
        
        #iterate over and print the boundary faces
        for i in range(0,polyFaces.getBoundaryPatchCount()+1):
            currentCounter = polyFaces.writeFaceNeighbourWithBoundaryId(file,i)

        file.write(')\n')

    #Read the boundary Tag file
    def read_boundary_tags(file_path):
        """
        Reads a boundary tag file.
        Format expected:
        Count
        ID Type Name
        ...

        Returns: dict { int_ID : str_Name }
        """
        boundary_map = {}

        if not os.path.exists(file_path):
            print(f"Warning: Boundary tag file not found at {file_path}. Using generic names.")
            return boundary_map

        with open(file_path, 'r') as f:
            lines = f.readlines()

            # We start at index 1 to skip the count header (index 0)
            # We check length to ensure the file isn't empty
            if len(lines) > 1:
                for line in lines[1:]:
                    parts = line.strip().split()
                    # Ensure we have at least ID, Type, Name (3 parts)
                    if len(parts) >= 3:
                        try:
                            # parts[0] is ID, parts[2] is Name
                            b_id = int(parts[0])
                            b_name = parts[2]
                            boundary_map[b_id] = b_name
                        except ValueError:
                            print(f"Skipping malformed line: {line.strip()}")
                            continue
                        
        print(f"Loaded {len(boundary_map)} boundary names.")
        return boundary_map

    boundaryNameMap = read_boundary_tags(boundaryTagFilePath)

    with open(boundaryFilename, 'w') as boundaryFile:
        boundaryFile.write(
            '/*--------------------------------*- C++ -*----------------------------------*\\\n'
            '| =========                 |                                                 |\n'
            '| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n'
            '|  \\\\    /   O peration     | Version:  1.7.1                                 |\n'
            '|   \\\\  /    A nd           | Web:      www.OpenFOAM.com                      |\n'
            '|    \\\\/     M anipulation  |                                                 |\n'
            '\\*---------------------------------------------------------------------------*/\n'
            'FoamFile\n'
            '{\n'
            '    version     2.0;\n'
            '    format      ascii;\n'
            '    class       polyBoundaryMesh;\n'
            '    location    "constant/polyMesh";\n'
            '    object      boundary;\n'
            '}\n'
            '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * /\n'
        )

        boundaryFile.write('\n\n')
        
        # Write the total number of patches
        numPatches = len(boundaryCellArray)
        boundaryFile.write('%i\n' % numPatches)
        boundaryFile.write('(\n')

        # Iterate through the arrays we populated earlier
        for i in range(numPatches):
            # Define a generic patch name (e.g., patch_1, patch_2)
            # Logic: In treatFace, you did boundaryId = prop[0] - 1
            # So if file had ID 1, you stored 0.
            # Here 'i' is the stored ID (0). So we look up (i + 1) to get the file ID.
            lookup_id = i + 1
            
            if lookup_id in boundaryNameMap:
                patchName = boundaryNameMap[lookup_id]
            else:
                patchName = f"patch_{i+1}" # Fallback
            
            boundaryFile.write(f"    {patchName}\n")
            boundaryFile.write( "    {\n")
            boundaryFile.write( "        type            patch;\n")
            boundaryFile.write(f"        nFaces          {boundaryCellCounter[i]};\n")
            boundaryFile.write(f"        startFace       {boundaryCellArray[i]};\n")
            boundaryFile.write( "    }\n")

        boundaryFile.write(')\n\n')
        boundaryFile.write('// ************************************************************************* //\n')


    ret = polyFaces.debugMissingNeighbours(ugrid.nodes,"testOutput")

    #print("Unconnected Faces: " + str(ret))

    print("Debug")














