import os
import sys
from pyNastran.converters.aflr.ugrid.ugrid_reader import read_ugrid
from pyNastran.converters.aflr.surf.surf_reader import SurfReader

# 1. Get the path of the current script (testScript.py)
current_dir = os.path.dirname(os.path.abspath(__file__))

# 2. construct the path to the root folder (3 levels up)
# polyMesh -> constant -> basicPipeCase -> Root
root_dir = os.path.abspath(os.path.join(current_dir, "../../.."))

# 3. Add that root directory to the system path so Python can "see" the files there
sys.path.append(root_dir)

# 4. Now you can import it simply by its filename (without .py)
import espAflr2Foam

# Get Folder This Script is in
path = os.path.realpath(__file__)
path = os.path.realpath(os.path.join(path , ".."))

uGrid = read_ugrid(os.path.join(path,"mesh.lb8.ugrid"))
boundarySurf = SurfReader()
boundarySurf.read_surf(os.path.join(path,"mesh.surf"))
boundaryTagFile = os.path.join(path,"mesh.mapbc")
espAflr2Foam.convertMesh(uGrid,boundarySurf,boundaryTagFile,path)
