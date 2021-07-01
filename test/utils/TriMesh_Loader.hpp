#pragma once

#include <map>
#include <cstring>
#include <fstream>

#include "src/util/OBJ_Loader.h"
#include "src/geometry/TriMesh.h"

//#define TRIMESH_VERBOSE

//std::vector<std::string> inline tokenize(std::string s, std::string delimiter = " ") {
//    // this method is taken directly from https://www.geeksforgeeks.org/how-to-split-a-string-in-cc-python-and-java/
//    int start = 0;
//    int end = s.find(delimiter);
//    std::vector<std::string> res;
//    while (end != -1) {
//        res.push_back(s.substr(start, end - start));
//        start = end + delimiter.size();
//        end = s.find(delimiter, start);
//    }
//    res.push_back(s.substr(start, end - start));
//    return res;
//}

// Adapted from https://stackoverflow.com/questions/10058606/splitting-a-string-by-a-character/10058756
__always_inline auto tokenize(const std::string& s, char delimiter) -> std::vector<std::string> {
  std::stringstream stream(s);
  std::string segment;
  std::vector<std::string> seglist;

  while(std::getline(stream, segment, delimiter)) {
    seglist.push_back(segment);
  }

  return seglist;
}

/**
 * Parses .obj file and per object reads vertex positions and faces, then stores them in encountered order in the
 * given arguments
 *
 * pseudocode:
 * 1. read line
 * - if it is a new object, remember or incerment the index or that object
 * - if it is a vertex, store that vertex
 * - if it is a face, store that face's vertices in the faces and the indices in the indices
 * - if it is sth irrelevant, continue
 * - if unexpected: throw a runtime_error after printing what it is.
 * - ???
 * - profit
 *
 * @param objFileName
 * @param vertices
 * @param indices
 * @param faces
 * @return number of found objects or negative on failure
 */
int load_meshes_from_obj(const std::string& objFileName,
                         std::vector<float>& vertices,
                         std::vector<std::vector<int>>& indices,
                         std::vector<std::vector<float>>& faces) {
    std::ifstream objFile(objFileName);
    if (!objFile || !objFile.is_open()) {
        std::cerr << "Failed to open " << objFileName << " for reading!" << std::endl;
        return -1;
    }

    int nth_obj = -1;

    std::string line;

    // handle initial lines until first object
    while (nth_obj < 0 && std::getline(objFile, line)) {
        if (line.starts_with("o")) {
#ifdef TRIMESH_VERBOSE
            std::cout << "\nConverting mesh " << line.substr(2, line.size()) << std::endl;
#endif
            nth_obj++;
            indices.emplace_back();
            faces.emplace_back();
            break;
        }
    }

    while (std::getline(objFile, line)) {
        std::vector<std::string> tokens = tokenize(line, ' ');
        if (tokens[0] == "o") {
#ifdef TRIMESH_VERBOSE
            std::cout << "\nConverting mesh " << tokens[1] << std::endl;
#endif
            nth_obj++;
            indices.emplace_back();
            faces.emplace_back();
        } else if (tokens[0] == "v") { // vertex
//            std::cout << "Extracting vertex " << line << std::endl;
//            std::cout << tokens[1] << std::endl;
//            std::cout << tokens[2] << std::endl;
//            std::cout << tokens[3] << std::endl;
            vertices.push_back(std::stof(tokens[1])); // x
            vertices.push_back(std::stof(tokens[2])); // y
            vertices.push_back(std::stof(tokens[3])); // z
        } else if (tokens[0] == "f") {
//            std::cout << "Extracting face and indices... " << line << std::endl;
            if (tokens.size() == 5) {
#ifdef TRIMESH_VERBOSE
              std::cout << "Found non triangle faces. Converting them to triangles." << std::endl;
#endif
              auto vertex_data_1 = tokenize(tokens[1], '/');
              int idx_1 = std::stoi(vertex_data_1[0]) - 1;
              auto vertex_data_2 = tokenize(tokens[2], '/');
              int idx_2 = std::stoi(vertex_data_2[0]) - 1;
              auto vertex_data_3 = tokenize(tokens[3], '/');
              int idx_3 = std::stoi(vertex_data_3[0]) - 1;
              auto vertex_data_4 = tokenize(tokens[4], '/');
              int idx_4 = std::stoi(vertex_data_4[0]) - 1;

              if (idx_1 * 3 + 2 >= vertices.size()
              || idx_2 * 3 + 2 >= vertices.size()
              || idx_3 * 3 + 2 >= vertices.size()
              || idx_4 * 3 + 2 >= vertices.size()) {
                return -1;
              }

              // One triangle
              indices[nth_obj].push_back(idx_1);
              indices[nth_obj].push_back(idx_2);
              indices[nth_obj].push_back(idx_3);
              // The other triangle
              indices[nth_obj].push_back(idx_3);
              indices[nth_obj].push_back(idx_4);
              indices[nth_obj].push_back(idx_1);

              // One triangle
              faces[nth_obj].push_back(vertices[idx_1 * 3]);
              faces[nth_obj].push_back(vertices[idx_1 * 3 + 1]);
              faces[nth_obj].push_back(vertices[idx_1 * 3 + 2]);
              faces[nth_obj].push_back(vertices[idx_2 * 3]);
              faces[nth_obj].push_back(vertices[idx_2 * 3 + 1]);
              faces[nth_obj].push_back(vertices[idx_2 * 3 + 2]);
              faces[nth_obj].push_back(vertices[idx_3 * 3]);
              faces[nth_obj].push_back(vertices[idx_3 * 3 + 1]);
              faces[nth_obj].push_back(vertices[idx_3 * 3 + 2]);

              // The other triangle
              faces[nth_obj].push_back(vertices[idx_3 * 3]);
              faces[nth_obj].push_back(vertices[idx_3 * 3 + 1]);
              faces[nth_obj].push_back(vertices[idx_3 * 3 + 2]);
              faces[nth_obj].push_back(vertices[idx_4 * 3]);
              faces[nth_obj].push_back(vertices[idx_4 * 3 + 1]);
              faces[nth_obj].push_back(vertices[idx_4 * 3 + 2]);
              faces[nth_obj].push_back(vertices[idx_1 * 3]);
              faces[nth_obj].push_back(vertices[idx_1 * 3 + 1]);
              faces[nth_obj].push_back(vertices[idx_1 * 3 + 2]);

              continue;
            }
            if (tokens.size() > 5 || tokens.size() < 4) {
              std::cout << tokens.size() << std::endl;
              std::cout << line << std::endl;
              throw std::runtime_error("Face with unsupported number of vertices!");
            }

            // iterate over referenced vertices in face
            // todo: assert that this (as given by blender) is actually counterclockwise
            for (size_t i = 1; i < 4; ++i) {
                // get only index and correct blender's indexing to start at zero
                auto vertex_data = tokenize(tokens[i], '/');
                int idx = std::stoi(vertex_data[0]) - 1;
                indices[nth_obj].push_back(idx);

                if (idx * 3 + 2 >= vertices.size()) return -1;

                faces[nth_obj].push_back(vertices[idx * 3]);
                faces[nth_obj].push_back(vertices[idx * 3 + 1]);
                faces[nth_obj].push_back(vertices[idx * 3 + 2]);
            }
        } else if (line.starts_with("vt ") ||
            line.starts_with("vn") ||
            line.starts_with("usemtl ") ||
            line.starts_with("mtl") ||
            line.starts_with("s ") ||
            line.starts_with("l ")) {
            continue;
        } else {
#ifdef TRIMESH_VERBOSE
            std::cout << "Unexpected identifier in line:\n" << line << std::endl;
#endif
          continue;
          //throw std::runtime_error("UNEXPECTED IDENTIFIER");
        }
    }

    return nth_obj + 1;
}

namespace tri_m_direct {

std::vector<TriMesh> load_meshes_to_TriMesh(const std::string& objFile) {
  std::vector<TriMesh> tri_meshes;
  std::vector<float> vertices;
  std::vector<std::vector<int>> indices;
  std::vector<std::vector<float>> faces;

  int n_meshes = load_meshes_from_obj(objFile, vertices, indices, faces);
  if (n_meshes < 0) {
    std::cout << "FAILED TO LOAD MESH!" << std::endl;
  }

  for (int i = 0; i < n_meshes; ++i) {
    struct Vec3f (* ffaces)[3] = static_cast<decltype(ffaces)>(malloc(faces[i].size() * sizeof(float)));
    std::memcpy(ffaces, faces[i].data(), faces[i].size() * sizeof(float));

    int n_faces = (static_cast<int>(indices[i].size())) / 3;
#ifdef TRIMESH_VERBOSE
    printf("Loader: Address of ffaces is %p\n", (void *)ffaces);
    printf("Loader: n_faces=%d\n", n_faces);
    printf("Loader: (faces.size() / 9)=%lu\n", faces[i].size() / 9);
#endif
    struct TriMesh t = {n_faces, ffaces};
    tri_meshes.push_back(t);
    //tri_meshes.push_back(TriMesh({n_faces, ffaces}));
  }
  return tri_meshes;
}

//std::vector<TriMeshI> load_meshes_to_TriMeshI(std::string objFile) {
//    std::vector<TriMeshI> tri_meshes;
//    std::vector<std::vector<float>> vertices;
//    std::vector<std::vector<int>> indices;
//    std::vector<std::vector<float>> faces;
//
//    int n_meshes = load_meshes_from_obj(objFile, vertices, indices, faces);
//
//    for (int i = 0; i < n_meshes; ++i) {
//        float (* vs)[3] = static_cast<decltype(vs)>(malloc(vertices[i].size() * sizeof(float)));
//        std::memcpy(vs, vertices[i].data(), vertices[i].size() * sizeof(float));
//
//        int
//        (* triangle_indices)[3] = static_cast<decltype(triangle_indices)>(malloc(indices[i].size() * sizeof(int)));
//        std::memcpy(triangle_indices, indices[i].data(), indices[i].size() * sizeof(int));
//
//        int n_vertices = vertices[i].size() / 3;
//        int n_faces = indices[i].size() / 3;
//
//        tri_meshes.push_back(TriMeshI({n_vertices, n_faces, vs, triangle_indices}));
//    }
//    return tri_meshes;
//}
}

void print_obj_meshes_info(objl::Loader& objLoader) {
#ifdef TRIMESH_VERBOSE
    std::cout << std::endl;
    std::cout << "OBJ_Loader loaded a total of "
              << objLoader.LoadedMeshes.size()
              << (objLoader.LoadedMeshes.size() == 1 ? " mesh, " : " meshes, ")
              << objLoader.LoadedVertices.size()
              << (objLoader.LoadedVertices.size() == 1 ? " vertex and " : " vertices and ")
              << objLoader.LoadedIndices.size() << (objLoader.LoadedIndices.size() == 1 ? " index:" : " indices:")
              << std::endl;
    for (auto& o : objLoader.LoadedMeshes) {
        std::cout << o.MeshName << " with "
                  << o.Vertices.size() << " vertices and "
                  << o.Indices.size() << " indices." << std::endl;
    }
#endif
}

/**
 * Reads out vertices from given Mesh o and populates remaining paramters
 * @param o objl::Mesh containing the vertex coordinates.
 * @param vertices empty std::vector<float>, gets filled with 3D vertex coordinates in encountered order
 * @param indices empty std::vector<int>, gets filled with new indices for triangles of vertices in encountered order
 * @param faces empty std::vector<float>, gets filled with coords of vertices for each face according to indices
 */
void populate(const objl::Mesh& o,
              std::vector<float>& vertices,
              std::vector<int>& indices,
              std::vector<float>& faces) {

    std::map<std::tuple<float, float, float>, int> index_map;
    std::vector<float> xs;
    std::vector<float> ys;
    std::vector<float> zs;

    int i = 0;
    for (auto& vertex : o.Vertices) {
        float x = vertex.Position.X;
        float y = vertex.Position.Y;
        float z = vertex.Position.Z;
        const auto& coords = std::make_tuple(x, y, z);

        faces.push_back(x);
        faces.push_back(y);
        faces.push_back(z);

        // Known vertex coordinates, reuse index
        const auto it = index_map.find(coords);
        if (it != index_map.end()) {
            indices.push_back(it->second);
            continue;
        }

        // New vertex coordinates
        index_map.insert(std::make_pair(coords, i));
        xs.push_back(x);
        ys.push_back(y);
        zs.push_back(z);
        indices.push_back(i++);
    }

    int n_vertices = xs.size();

    vertices.reserve(n_vertices * 3);
    for (int j = 0; j < n_vertices; ++j) {
        vertices.emplace_back(xs[j]);
        vertices.emplace_back(ys[j]);
        vertices.emplace_back(zs[j]);
    }

    assert(faces.size() == indices.size() * 3);
    assert(indices.size() % 3 == 0);
    assert(faces.size() % 3 == 0);
}

namespace tri_m_I {

void print_trimeshes_info(const std::vector<TriMeshI>& trimeshes,
                          std::vector<std::string> mesh_names,
                          int verbose) {
#ifdef TRIMESH_VERBOSE
    std::cout << std::endl;
    std::cout << "Converted " << trimeshes.size() << (trimeshes.size() == 1 ? " mesh " : " meshes ")
              << "to struct TriMeshI. " << std::endl;

    int nth_mesh = 0;
    for (const auto& trimesh : trimeshes) {
        std::cout << mesh_names[nth_mesh++] << " with "
                  << trimesh.num_vertices
                  << (trimesh.num_vertices == 1 ? " unique vertex and " : " unique vertices and ")
                  << trimesh.num_faces << (trimesh.num_faces == 1 ? " triangle. " : " triangles. ")
                  << std::endl;

        if (verbose <= 0) continue;

        std::cout << "Vertices:" << std::endl;
        for (int i = 0; i < trimesh.num_vertices; ++i)
            std::cout << trimesh.vs[i][0] << " " << trimesh.vs[i][1] << " " << trimesh.vs[i][2] << " "
                      << std::endl;

        std::cout << "Indices: " << std::endl;
        for (int i = 0; i < trimesh.num_faces; ++i)
            std::cout << trimesh.idx[i][0] << " " << trimesh.idx[i][1] << " " << trimesh.idx[i][2] << std::endl;
    }
#endif
}

/**
 * Test object loader utility and conversion of .obj meshes into trimeshes
 * Assumes that each face is a triangle. I.e. when exporting objects from Blender,
 * manually select 'Geometry > Triangulate Faces'.
 */
std::vector<TriMeshI> load_meshes(std::string objFile, int verbosity) {
#ifdef TRIMESH_VERBOSE
    std::cout << std::endl << "Attempting to load meshes..." << std::endl;
#endif

    // load the .obj file
    objl::Loader objLoader;
    if (!objLoader.LoadFile(objFile)) {
#ifdef TRIMESH_VERBOSE
        std::cout << "FAILED to load File " << objFile << std::endl;
#endif
        return std::vector<TriMeshI>();
    }

#ifdef TRIMESH_VERBOSE
    std::cout << "Loaded file " << objFile << std::endl;
#endif
    print_obj_meshes_info(objLoader);

    // populate trimeshes
    std::vector<TriMeshI> tri_meshes;
    std::vector<std::string> mesh_names;
    for (auto& o : objLoader.LoadedMeshes) {
#ifdef TRIMESH_VERBOSE
        std::cout << "Converting mesh " << o.MeshName << std::endl;
#endif
        mesh_names.push_back(o.MeshName);

        std::vector<float> vertices;
        std::vector<int> indices;
        std::vector<float> faces;

        populate(o, vertices, indices, faces);

        float (* vs)[3] = static_cast<decltype(vs)>(malloc(vertices.size() * sizeof(float)));
        std::memcpy(vs, vertices.data(), vertices.size() * sizeof(float));

        int
        (* triangle_indices)[3] = static_cast<decltype(triangle_indices)>(malloc(indices.size() * sizeof(int)));
        std::memcpy(triangle_indices, indices.data(), indices.size() * sizeof(int));

        int n_vertices = vertices.size() / 3;
        int n_faces = indices.size() / 3;

        tri_meshes.push_back(TriMeshI({n_vertices, n_faces, vs, triangle_indices}));
    }

    print_trimeshes_info(tri_meshes, mesh_names, verbosity);
    return tri_meshes;
}

int runTest_load_triangulated_meshes() {
//    std::string objFile = "../../Data/objects/simple_cube_triangulated.obj";
    std::string objFile = "../../Data/objects/fractured_cube_frames/fracturetest_triangulated.obj";
    const std::vector<TriMeshI>& trimeshes = load_meshes(objFile, 0);
    return 0;
}

}

namespace tri_m {

void print_trimeshes_info(const std::vector<TriMesh>& trimeshes,
                          std::vector<std::string> mesh_names,
                          int verbose) {
#ifdef TRIMESH_VERBOSE
    std::cout << std::endl;
    std::cout << "Converted " << trimeshes.size() << (trimeshes.size() == 1 ? " mesh " : " meshes ")
              << "to struct TriMesh. " << std::endl;
    int nth_mesh = 0;
    for (const auto& trimesh : trimeshes) {
        std::cout << mesh_names[nth_mesh++] << " with "
                  << trimesh.num_faces << (trimesh.num_faces == 1 ? " triangle. " : " triangles. ")
                  << std::endl;

        if (verbose <= 0) continue;

        std::cout << "Faces: " << std::endl;
        for (int i = 0; i < trimesh.num_faces; ++i)
            std::cout
                << trimesh.faces[i][0].x << " "
                << trimesh.faces[i][0].y << " "
                << trimesh.faces[i][0].z << std::endl

                << trimesh.faces[i][1].x << " "
                << trimesh.faces[i][1].y << " "
                << trimesh.faces[i][1].z << std::endl

                << trimesh.faces[i][2].x << " "
                << trimesh.faces[i][2].y << " "
                << trimesh.faces[i][2].z << std::endl << std::endl;
    }
#endif
}

/**
 * Test object loader utility and conversion of .obj meshes into trimeshes
 * Assumes that each face is a triangle. I.e. when exporting objects from Blender,
 * manually select 'Geometry > Triangulate Faces'.
 */
std::vector<TriMesh> load_meshes(const std::string& objFile, int verbosity) {
#ifdef TRIMESH_VERBOSE
    std::cout << std::endl << "Attempting to load meshes..." << std::endl;
#endif

    // load the .obj file
    objl::Loader objLoader;
    if (!objLoader.LoadFile(objFile)) {
#ifdef TRIMESH_VERBOSE
        std::cout << "FAILED to load File " << objFile << std::endl;
#endif
        return std::vector<TriMesh>();
    }

#ifdef TRIMESH_VERBOSE
    std::cout << "Loaded file " << objFile << std::endl;
#endif
    print_obj_meshes_info(objLoader);

    // populate trimeshes
    std::vector<TriMesh> tri_meshes;
    std::vector<std::string> mesh_names;
    for (auto& o : objLoader.LoadedMeshes) {
#ifdef TRIMESH_VERBOSE
        std::cout << "Converting mesh " << o.MeshName << std::endl;
#endif
        mesh_names.push_back(o.MeshName);

        std::vector<float> vertices;
        std::vector<int> indices;
        std::vector<float> faces;

        populate(o, vertices, indices, faces);

        struct Vec3f (* ffaces)[3] = static_cast<decltype(ffaces)>(malloc(faces.size() * sizeof(float)));
        std::memcpy(ffaces, faces.data(), faces.size() * sizeof(float));

        int n_faces = indices.size() / 3;

        tri_meshes.push_back(TriMesh({n_faces, ffaces}));
    }

    print_trimeshes_info(tri_meshes, mesh_names, verbosity);
    return tri_meshes;
}

}
