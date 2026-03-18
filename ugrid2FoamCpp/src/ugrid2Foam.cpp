#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <filesystem>
#include <chrono>
#include <cctype>

namespace fs = std::filesystem;

// ============================================================================
// Data Structures & Math Utilities
// ============================================================================

struct Point3D {
    double x, y, z;
    Point3D() : x(0), y(0), z(0) {}
    Point3D(double x, double y, double z) : x(x), y(y), z(z) {}

    Point3D operator-(const Point3D& o) const { return Point3D(x - o.x, y - o.y, z - o.z); }
    Point3D operator+(const Point3D& o) const { return Point3D(x + o.x, y + o.y, z + o.z); }
    Point3D operator*(double s) const { return Point3D(x * s, y * s, z * s); }
    Point3D operator/(double s) const { return Point3D(x / s, y / s, z / s); }
    
    double dot(const Point3D& o) const { return x * o.x + y * o.y + z * o.z; }
    Point3D cross(const Point3D& o) const {
        return Point3D(y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x);
    }
    double distSq(const Point3D& o) const {
        double dx = x - o.x, dy = y - o.y, dz = z - o.z;
        return dx*dx + dy*dy + dz*dz;
    }
};

struct Face {
    int nodes[8];
    int num_nodes;
    int owner;
    int neighbour;
    int boundaryId;
    int boundaryType;
    bool valid;

    Face() : num_nodes(0), owner(-1), neighbour(-1), boundaryId(-1), boundaryType(-1), valid(true) {
        std::fill(nodes, nodes + 8, -1);
    }
};

struct FaceSignature {
    int sorted_nodes[8]; 
    int num_nodes;
    int face_index;

    FaceSignature() : num_nodes(0), face_index(-1) {
        std::fill(sorted_nodes, sorted_nodes + 8, -1);
    }
    
    FaceSignature(const int* n, int num, int idx = -1) : num_nodes(num), face_index(idx) {
        std::fill(sorted_nodes, sorted_nodes + 8, -1);
        for(int i = 0; i < num; ++i) sorted_nodes[i] = n[i];
        std::sort(sorted_nodes, sorted_nodes + num);
    }

    bool operator<(const FaceSignature& o) const {
        if (num_nodes != o.num_nodes) return num_nodes < o.num_nodes;
        for (int i = 0; i < num_nodes; ++i) {
            if (sorted_nodes[i] != o.sorted_nodes[i]) return sorted_nodes[i] < o.sorted_nodes[i];
        }
        return false;
    }
    
    bool operator==(const FaceSignature& o) const {
        if (num_nodes != o.num_nodes) return false;
        for (int i = 0; i < num_nodes; ++i) {
            if (sorted_nodes[i] != o.sorted_nodes[i]) return false;
        }
        return true;
    }
};

// ============================================================================
// KD-Tree Implementation (Array-Backed)
// ============================================================================

class KDTree {
    struct KDNode {
        Point3D pt;
        int original_idx;
        int axis;
    };
    std::vector<KDNode> nodes;

    void build(int start, int end, int depth) {
        if (start >= end) return;
        int axis = depth % 3;
        int mid = start + (end - start) / 2;

        std::nth_element(nodes.begin() + start, nodes.begin() + mid, nodes.begin() + end,
            [axis](const KDNode& a, const KDNode& b) {
                if (axis == 0) return a.pt.x < b.pt.x;
                if (axis == 1) return a.pt.y < b.pt.y;
                return a.pt.z < b.pt.z;
            });

        nodes[mid].axis = axis;
        build(start, mid, depth + 1);
        build(mid + 1, end, depth + 1);
    }

    void nearestNeighborSearch(int start, int end, const Point3D& target, int& best_idx, double& best_distSq) const {
        if (start >= end) return;
        int mid = start + (end - start) / 2;
        const KDNode& node = nodes[mid];

        double dSq = target.distSq(node.pt);
        if (dSq < best_distSq) {
            best_distSq = dSq;
            best_idx = node.original_idx;
        }

        double diff = (node.axis == 0) ? (target.x - node.pt.x) :
                      (node.axis == 1) ? (target.y - node.pt.y) : (target.z - node.pt.z);

        int first = diff < 0 ? start : mid + 1;
        int last = diff < 0 ? mid : end;
        int other_first = diff < 0 ? mid + 1 : start;
        int other_last = diff < 0 ? end : mid;

        nearestNeighborSearch(first, last, target, best_idx, best_distSq);
        if (diff * diff < best_distSq) {
            nearestNeighborSearch(other_first, other_last, target, best_idx, best_distSq);
        }
    }

public:
    KDTree(const std::vector<Point3D>& pts, const std::vector<int>& original_ids = {}) {
        nodes.reserve(pts.size());
        bool use_ids = !original_ids.empty();
        for (size_t i = 0; i < pts.size(); ++i) {
            nodes.push_back({pts[i], use_ids ? original_ids[i] : (int)i, 0});
        }
        build(0, nodes.size(), 0);
    }

    int nearest(const Point3D& target, double& out_dist) const {
        if (nodes.empty()) return -1;
        int best_idx = -1;
        double best_distSq = 1e30;
        nearestNeighborSearch(0, nodes.size(), target, best_idx, best_distSq);
        out_dist = std::sqrt(best_distSq);
        return best_idx;
    }
};

// ============================================================================
// Endianness & I/O Utilities
// ============================================================================

bool isLittleEndian() {
    uint16_t number = 0x1;
    char *numPtr = (char*)&number;
    return (numPtr[0] == 1);
}

template <typename T>
void swapEndian(T& val) {
    char* ptr = reinterpret_cast<char*>(&val);
    std::reverse(ptr, ptr + sizeof(T));
}

// ============================================================================
// UGRID Reader
// ============================================================================

struct UGRID {
    std::vector<Point3D> nodes;
    std::vector<int> tris, quads, pids;
    std::vector<int> tets, penta5s, penta6s, hexas;

    bool load(const std::string& filename) {
        std::ifstream file(filename, std::ios::binary);
        if (!file) throw std::runtime_error("Cannot open UGRID file: " + filename);

        bool fileIsBigEndian = (filename.find(".b8.ugrid") != std::string::npos || filename.find(".b4.ugrid") != std::string::npos);
        if (filename.find(".lb8.") != std::string::npos || filename.find(".lb4.") != std::string::npos) {
            fileIsBigEndian = false;
        }
        bool needsSwap = (isLittleEndian() == fileIsBigEndian);
        bool isDouble = (filename.find("8.ugrid") != std::string::npos);

        int32_t header[7];
        file.read(reinterpret_cast<char*>(header), 7 * sizeof(int32_t));
        if (needsSwap) for (int i = 0; i < 7; ++i) swapEndian(header[i]);

        int nnodes = header[0], ntris = header[1], nquads = header[2];
        int ntets = header[3], npenta5 = header[4], npenta6 = header[5], nhexas = header[6];

        std::cout << "Loading UGRID: nodes=" << nnodes << ", tets=" << ntets << ", penta5s=" << npenta5 
                  << ", penta6s=" << npenta6 << ", hexas=" << nhexas << "\n";

        nodes.resize(nnodes);
        if (isDouble) {
            std::vector<double> buf(nnodes * 3);
            file.read(reinterpret_cast<char*>(buf.data()), buf.size() * sizeof(double));
            if (needsSwap) for (auto& v : buf) swapEndian(v);
            for (int i = 0; i < nnodes; ++i) nodes[i] = Point3D(buf[i*3], buf[i*3+1], buf[i*3+2]);
        } else {
            std::vector<float> buf(nnodes * 3);
            file.read(reinterpret_cast<char*>(buf.data()), buf.size() * sizeof(float));
            if (needsSwap) for (auto& v : buf) swapEndian(v);
            for (int i = 0; i < nnodes; ++i) nodes[i] = Point3D(buf[i*3], buf[i*3+1], buf[i*3+2]);
        }

        auto readInts = [&](std::vector<int>& vec, int count, int multiplier) {
            vec.resize(count * multiplier);
            file.read(reinterpret_cast<char*>(vec.data()), vec.size() * sizeof(int));
            if (needsSwap) for (auto& v : vec) swapEndian(v);
        };

        readInts(tris, ntris, 3);
        readInts(quads, nquads, 4);
        readInts(pids, ntris + nquads, 1);
        readInts(tets, ntets, 4);
        readInts(penta5s, npenta5, 5);
        readInts(penta6s, npenta6, 6);
        readInts(hexas, nhexas, 8);

        return true;
    }
};

// ============================================================================
// SURF and TAGS Readers
// ============================================================================

struct SurfElement {
    int nodes[4];
    int num_nodes;
    int surface_id;
    int grid_bc_flag;
};

struct SURF {
    std::vector<Point3D> nodes;
    std::vector<SurfElement> tris;
    std::vector<SurfElement> quads;

    bool load(const std::string& filename) {
        std::ifstream file(filename);
        if (!file) {
            std::cerr << "Warning: Cannot open SURF file: " << filename << "\n";
            return false;
        }

        std::string line;
        if (!std::getline(file, line)) return false;
        std::stringstream ss(line);
        int ntris = 0, nquads = 0, nnodes = 0;
        ss >> ntris >> nquads >> nnodes;
        
        std::cout << "Loading SURF: nodes=" << nnodes << ", tris=" << ntris << ", quads=" << nquads << "\n";

        nodes.resize(nnodes);
        for (int i = 0; i < nnodes; ++i) {
            std::getline(file, line);
            std::stringstream nss(line);
            nss >> nodes[i].x >> nodes[i].y >> nodes[i].z;
        }

        tris.resize(ntris);
        for (int i = 0; i < ntris; ++i) {
            std::getline(file, line);
            std::stringstream tss(line);
            int recon;
            tss >> tris[i].nodes[0] >> tris[i].nodes[1] >> tris[i].nodes[2] 
                >> tris[i].surface_id >> recon >> tris[i].grid_bc_flag;
            tris[i].num_nodes = 3;
        }

        quads.resize(nquads);
        for (int i = 0; i < nquads; ++i) {
            std::getline(file, line);
            std::stringstream qss(line);
            int recon;
            qss >> quads[i].nodes[0] >> quads[i].nodes[1] >> quads[i].nodes[2] >> quads[i].nodes[3] 
                >> quads[i].surface_id >> recon >> quads[i].grid_bc_flag;
            quads[i].num_nodes = 4;
        }

        return true;
    }
};

struct TagReader {
    std::map<int, std::string> tags;

    bool load(const std::string& filename) {
        std::ifstream file(filename);
        if (!file) {
            std::cerr << "Warning: Cannot open TAGS file: " << filename << ". Using generic patch names.\n";
            return false;
        }

        std::string line;
        int count = 0;
        while (std::getline(file, line)) {
            size_t pos = line.find('#');
            if (pos != std::string::npos) line = line.substr(0, pos);
            if (line.empty()) continue;

            std::stringstream ss(line);
            int id;
            std::string name;
            
            if (ss >> id >> name) {
                tags[id] = name;
                count++;
            }
        }
        std::cout << "Loaded " << count << " boundary tags.\n";
        return true;
    }
};

// ============================================================================
// Converter Logic & Geometry
// ============================================================================

void ensureOutward(int* face_nodes, int num_nodes, const Point3D& cell_center, const std::vector<Point3D>& coords) {
    Point3D p1 = coords[face_nodes[0] - 1]; 
    Point3D p2 = coords[face_nodes[1] - 1];
    Point3D p3 = coords[face_nodes[2] - 1];

    Point3D face_center;
    if (num_nodes == 3) face_center = (p1 + p2 + p3) / 3.0;
    else {
        Point3D p4 = coords[face_nodes[3] - 1];
        face_center = (p1 + p2 + p3 + p4) / 4.0;
    }

    Point3D normal = (p2 - p1).cross(p3 - p1);
    Point3D center_to_face = face_center - cell_center;

    if (normal.dot(center_to_face) < 0) {
        if (num_nodes == 3) {
            std::swap(face_nodes[1], face_nodes[2]);
        } else {
            std::swap(face_nodes[1], face_nodes[3]); // 0,1,2,3 -> 0,3,2,1
        }
    }
}

std::string getOpenFOAMPatchType(const std::string& patchName) {
    std::string lowerName = patchName;
    std::transform(lowerName.begin(), lowerName.end(), lowerName.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    
    if (lowerName.find("wall") != std::string::npos) return "wall";
    if (lowerName.find("symm") != std::string::npos) return "symmetry";
    if (lowerName.find("empty") != std::string::npos) return "empty";
    if (lowerName.find("slip") != std::string::npos) return "symmetryPlane";
    if (lowerName.find("cyclic") != std::string::npos) return "cyclic";
    
    return "patch"; 
}

void writeFoamHeader(std::ofstream& out, const std::string& objectType, const std::string& objectClass) {
    out << "/*--------------------------------*- C++ -*----------------------------------*\\\n"
        << "| =========                 |                                                 |\n"
        << "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n"
        << "|  \\\\    /   O peration     | Version:  2.0.x                                 |\n"
        << "|   \\\\  /    A nd           | Web:      www.OpenFOAM.com                      |\n"
        << "|    \\\\/     M anipulation  |                                                 |\n"
        << "\\*---------------------------------------------------------------------------*/\n"
        << "FoamFile\n{\n    version     2.0;\n    format      ascii;\n    class       " 
        << objectClass << ";\n    location    \"constant/polyMesh\";\n    object      " 
        << objectType << ";\n}\n// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";
}

const std::string FOAM_CLOSE = ");\n\n// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n";

// ============================================================================
// Main Execution
// ============================================================================

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <input.ugrid> <output_dir> [input.surf] [input.tags]\n";
        return 1;
    }

    std::string inputFile = argv[1];
    std::string outputDir = argv[2];
    std::string polyMeshDir = outputDir + "/constant/polyMesh";
    
    std::string surfFile = (argc > 3) ? argv[3] : "";
    std::string tagsFile = (argc > 4) ? argv[4] : "";
    bool rematchIceInterface = true;

    try {
        UGRID ugrid;
        SURF surf;
        TagReader tags;

        auto start = std::chrono::high_resolution_clock::now();
        ugrid.load(inputFile);
        
        bool hasSurf = false;
        if (!surfFile.empty()) hasSurf = surf.load(surfFile);
        if (!tagsFile.empty()) tags.load(tagsFile);

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end - start;
        std::cout << "Data loaded in " << diff.count() << " seconds.\n";

        // Topological Mapping
        std::map<FaceSignature, int> volSurfaceMap;
        for (size_t i = 0; i < ugrid.tris.size(); i += 3) {
            int nodes[3] = { ugrid.tris[i]-1, ugrid.tris[i+1]-1, ugrid.tris[i+2]-1 };
            volSurfaceMap[FaceSignature(nodes, 3)] = ugrid.pids[i/3] - 1; 
        }
        int num_tris = ugrid.tris.size() / 3;
        for (size_t i = 0; i < ugrid.quads.size(); i += 4) {
            int nodes[4] = { ugrid.quads[i]-1, ugrid.quads[i+1]-1, ugrid.quads[i+2]-1, ugrid.quads[i+3]-1 };
            volSurfaceMap[FaceSignature(nodes, 4)] = ugrid.pids[num_tris + i/4] - 1;
        }

        // Spatial Mapping (KDTree)
        std::map<FaceSignature, int> spatialSurfMap;
        if (hasSurf && !surf.nodes.empty()) {
            std::cout << "Building KD-Tree for spatial mapping...\n";
            KDTree tree(ugrid.nodes);
            std::vector<int> nodeMap(surf.nodes.size());
            
            double dist;
            for (size_t i = 0; i < surf.nodes.size(); ++i) {
                nodeMap[i] = tree.nearest(surf.nodes[i], dist);
            }

            for (const auto& tri : surf.tris) {
                int n[3] = { nodeMap[tri.nodes[0]-1], nodeMap[tri.nodes[1]-1], nodeMap[tri.nodes[2]-1] };
                spatialSurfMap[FaceSignature(n, 3)] = tri.surface_id - 1;
            }
            for (const auto& quad : surf.quads) {
                int n[4] = { nodeMap[quad.nodes[0]-1], nodeMap[quad.nodes[1]-1], nodeMap[quad.nodes[2]-1], nodeMap[quad.nodes[3]-1] };
                spatialSurfMap[FaceSignature(n, 4)] = quad.surface_id - 1;
            }
        }

        std::vector<Face> all_faces;
        size_t est_faces = (ugrid.tets.size()/4)*4 + (ugrid.penta5s.size()/5)*5 + 
                           (ugrid.penta6s.size()/6)*5 + (ugrid.hexas.size()/8)*6;
        all_faces.reserve(est_faces);

        std::cout << "Extracting volume faces...\n";
        int cellID = 0;
        
        auto processCellFaces = [&](const std::vector<std::vector<int>>& local_faces, const Point3D& cell_center) {
            for(const auto& lf : local_faces) {
                Face f;
                f.num_nodes = lf.size();
                for(size_t i=0; i<lf.size(); ++i) f.nodes[i] = lf[i];
                ensureOutward(f.nodes, f.num_nodes, cell_center, ugrid.nodes);
                f.owner = cellID;
                for(int i=0; i<f.num_nodes; ++i) f.nodes[i] -= 1; // 0-based
                all_faces.push_back(f);
            }
            cellID++;
        };

        // TETRAHEDRA
        for (size_t i = 0; i < ugrid.tets.size(); i += 4) {
            int n1 = ugrid.tets[i], n2 = ugrid.tets[i+1], n3 = ugrid.tets[i+2], n4 = ugrid.tets[i+3];
            Point3D center = (ugrid.nodes[n1-1] + ugrid.nodes[n2-1] + ugrid.nodes[n3-1] + ugrid.nodes[n4-1]) / 4.0;
            processCellFaces({{n1,n3,n2}, {n2,n3,n4}, {n1,n4,n3}, {n1,n2,n4}}, center);
        }

        // PENTA5
        for (size_t i = 0; i < ugrid.penta5s.size(); i += 5) {
            int n[5]; for(int j=0; j<5; ++j) n[j] = ugrid.penta5s[i+j];
            Point3D center; 
            for(int j=0; j<5; ++j) center = center + ugrid.nodes[n[j]-1];
            center = center / 5.0;
            processCellFaces({
                {n[1], n[2], n[4]}, {n[3], n[4], n[2]}, {n[0], n[3], n[2]}, 
                {n[0], n[2], n[1]}, {n[0], n[1], n[4], n[3]}
            }, center);
        }

        // PENTA6
        for (size_t i = 0; i < ugrid.penta6s.size(); i += 6) {
            int n[6]; for(int j=0; j<6; ++j) n[j] = ugrid.penta6s[i+j];
            Point3D center; 
            for(int j=0; j<6; ++j) center = center + ugrid.nodes[n[j]-1];
            center = center / 6.0;
            processCellFaces({
                {n[0], n[3], n[4], n[1]}, {n[0], n[2], n[5], n[3]}, 
                {n[2], n[1], n[4], n[5]}, {n[0], n[1], n[2]}, {n[3], n[5], n[4]}
            }, center);
        }

        // HEXAHEDRA
        for (size_t i = 0; i < ugrid.hexas.size(); i += 8) {
            int n[8]; for(int j=0; j<8; ++j) n[j] = ugrid.hexas[i+j];
            Point3D center; 
            for(int j=0; j<8; ++j) center = center + ugrid.nodes[n[j]-1];
            center = center / 8.0;
            processCellFaces({
                {n[0], n[3], n[2], n[1]}, {n[0], n[1], n[5], n[4]}, 
                {n[3], n[7], n[6], n[2]}, {n[4], n[5], n[6], n[7]}, 
                {n[0], n[4], n[7], n[3]}, {n[1], n[2], n[6], n[5]}
            }, center);
        }

        std::vector<FaceSignature> sigs(all_faces.size());
        for (size_t i = 0; i < all_faces.size(); ++i) {
            sigs[i] = FaceSignature(all_faces[i].nodes, all_faces[i].num_nodes, i);
        }

        std::sort(sigs.begin(), sigs.end());

        int matched = 0;
        for (size_t i = 0; i < sigs.size();) {
            size_t j = i + 1;
            while (j < sigs.size() && sigs[i] == sigs[j]) j++;
            
            int count = j - i;
            if (count == 2) {
                int f1 = sigs[i].face_index;
                int f2 = sigs[i+1].face_index;
                
                int o1 = all_faces[f1].owner;
                int o2 = all_faces[f2].owner;
                
                int new_owner = std::min(o1, o2);
                int new_neigh = std::max(o1, o2);

                all_faces[f1].owner = new_owner;
                all_faces[f1].neighbour = new_neigh;
                if (new_owner != o1) std::reverse(all_faces[f1].nodes, all_faces[f1].nodes + all_faces[f1].num_nodes);
                
                all_faces[f2].valid = false;
                matched++;
            }
            i = j;
        }

        // ====================================================================
        // Ice Interface Rematching (Edge Splitting)
        // ====================================================================
        if (rematchIceInterface) {
            std::cout << "Rematching Ice Interface (Spatial Edge Splitting)...\n";
            std::set<int> interface_node_set;
            
            // Gather unique nodes from currently unmatched boundary faces
            for (const auto& f : all_faces) {
                if (f.valid && f.neighbour == -1) {
                    for (int i = 0; i < f.num_nodes; ++i) {
                        interface_node_set.insert(f.nodes[i]);
                    }
                }
            }
            
            std::vector<int> interface_node_ids(interface_node_set.begin(), interface_node_set.end());
            std::vector<Point3D> interface_coords;
            interface_coords.reserve(interface_node_ids.size());
            for (int nid : interface_node_ids) {
                interface_coords.push_back(ugrid.nodes[nid]);
            }
            
            KDTree interfaceTree(interface_coords, interface_node_ids);
            std::map<std::pair<int, int>, int> edgeSplitingIndex;
            int snapped_count = 0;
            double distToleranceSq = 1e-8; // 1e-4 squared
            
            // Search for virtual midpoints
            for (const auto& f : all_faces) {
                if (f.valid && f.neighbour == -1 && f.num_nodes == 4) {
                    for (int i = 0; i < 4; ++i) {
                        int nA = f.nodes[i];
                        int nB = f.nodes[(i + 1) % 4];
                        
                        Point3D pA = ugrid.nodes[nA];
                        Point3D pB = ugrid.nodes[nB];
                        Point3D midpoint = (pA + pB) * 0.5;
                        
                        double dist;
                        int found_nid = interfaceTree.nearest(midpoint, dist);
                        if (found_nid != -1 && (dist * dist) < distToleranceSq) {
                            auto key = std::make_pair(std::min(nA, nB), std::max(nA, nB));
                            edgeSplitingIndex[key] = found_nid;
                            snapped_count++;
                        }
                    }
                }
            }
            
            // Apply splits
            if (snapped_count > 0) {
                std::cout << "Applying " << edgeSplitingIndex.size() << " edge splits to quad boundaries...\n";
                for (auto& f : all_faces) {
                    if (f.valid && f.num_nodes == 4) { // Only checking quads
                        std::vector<int> new_nodes;
                        new_nodes.reserve(8);
                        bool split_applied = false;
                        
                        for (int i = 0; i < 4; ++i) {
                            new_nodes.push_back(f.nodes[i]);
                            int nA = f.nodes[i];
                            int nB = f.nodes[(i + 1) % 4];
                            auto key = std::make_pair(std::min(nA, nB), std::max(nA, nB));
                            
                            if (edgeSplitingIndex.count(key)) {
                                new_nodes.push_back(edgeSplitingIndex[key]);
                                split_applied = true;
                            }
                        }
                        
                        if (split_applied) {
                            f.num_nodes = new_nodes.size();
                            for (size_t i = 0; i < new_nodes.size(); ++i) {
                                f.nodes[i] = new_nodes[i];
                            }
                        }
                    }
                }
            }
        }
        // ====================================================================

        // ID Mapping Check
        for (auto& f : all_faces) {
            if (f.valid) {
                if (f.neighbour != -1) {
                    f.boundaryId = -1; 
                } else {
                    FaceSignature sig(f.nodes, f.num_nodes);
                    auto itTop = volSurfaceMap.find(sig);
                    if (itTop != volSurfaceMap.end()) {
                        f.boundaryId = itTop->second;
                    } else {
                        auto itSpa = spatialSurfMap.find(sig);
                        if (itSpa != spatialSurfMap.end()) {
                            f.boundaryId = itSpa->second;
                        } else {
                            f.boundaryId = 0; // Fallback
                        }
                    }
                }
            }
        }

        std::vector<const Face*> valid_faces;
        valid_faces.reserve(all_faces.size() - matched);
        for (const auto& f : all_faces) {
            if (f.valid) valid_faces.push_back(&f);
        }

        std::sort(valid_faces.begin(), valid_faces.end(), [](const Face* a, const Face* b) {
            if (a->boundaryId != b->boundaryId) return a->boundaryId < b->boundaryId;
            if (a->owner != b->owner) return a->owner < b->owner;
            return a->neighbour < b->neighbour;
        });

        fs::create_directories(polyMeshDir);

        std::ofstream pts(polyMeshDir + "/points");
        writeFoamHeader(pts, "points", "vectorField");
        pts << ugrid.nodes.size() << "\n(\n";
        for (const auto& p : ugrid.nodes) pts << "    (" << p.x << " " << p.y << " " << p.z << ")\n";
        pts << FOAM_CLOSE;
        pts.close();

        std::ofstream fcs(polyMeshDir + "/faces");
        writeFoamHeader(fcs, "faces", "faceList");
        fcs << valid_faces.size() << "\n(\n";
        for (const auto* f : valid_faces) {
            fcs << f->num_nodes << "(";
            for (int n = 0; n < f->num_nodes; ++n) fcs << (n>0?" ":"") << f->nodes[n];
            fcs << ")\n";
        }
        fcs << FOAM_CLOSE;
        fcs.close();

        std::ofstream own(polyMeshDir + "/owner");
        writeFoamHeader(own, "owner", "labelList");
        own << valid_faces.size() << "\n(\n";
        for (const auto* f : valid_faces) own << f->owner << "\n";
        own << FOAM_CLOSE;
        own.close();

        int internal_faces = 0;
        for (const auto* f : valid_faces) { if (f->boundaryId == -1) internal_faces++; }
        
        std::ofstream nei(polyMeshDir + "/neighbour");
        writeFoamHeader(nei, "neighbour", "labelList");
        nei << internal_faces << "\n(\n";
        for (const auto* f : valid_faces) {
            if (f->boundaryId == -1) nei << f->neighbour << "\n";
        }
        nei << FOAM_CLOSE;
        nei.close();

        std::ofstream bnd(polyMeshDir + "/boundary");
        writeFoamHeader(bnd, "boundary", "polyBoundaryMesh");
        
        std::map<int, int> patch_counts;
        std::map<int, int> patch_starts;
        int current_idx = 0;
        
        for (const auto* f : valid_faces) {
            if (f->boundaryId != -1) {
                if (patch_counts.find(f->boundaryId) == patch_counts.end()) {
                    patch_starts[f->boundaryId] = current_idx;
                    patch_counts[f->boundaryId] = 0;
                }
                patch_counts[f->boundaryId]++;
            }
            current_idx++;
        }

        bnd << patch_counts.size() << "\n(\n";
        for (const auto& [bId, count] : patch_counts) {
            std::string pName = "patch_" + std::to_string(bId);
            if (tags.tags.find(bId) != tags.tags.end()) {
                pName = tags.tags[bId];
            } else if (tags.tags.find(bId + 1) != tags.tags.end()) {
                pName = tags.tags[bId + 1];
            }

            std::string pType = getOpenFOAMPatchType(pName);

            bnd << "    " << pName << "\n    {\n"
                << "        type            " << pType << ";\n"
                << "        nFaces          " << count << ";\n"
                << "        startFace       " << patch_starts[bId] << ";\n    }\n";
        }
        bnd << FOAM_CLOSE;
        bnd.close();

        std::cout << "Conversion complete! Mesh exported cleanly.\n";

    } catch (const std::exception& e) {
        std::cerr << "Fatal Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}