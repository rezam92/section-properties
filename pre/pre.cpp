//
// Created by Reza on 6/10/22.
//

#include "pre.h"
#include "../libdistmesh/include/distmesh/triangulation.h"

#include <utility>

Material::Material(std::string n, float em, float pr, float ys, float den, std::string col) {
    this->name = std::move(n);
    this->elastic_modulus = em;
    this->poissons_ratio = pr;
    this->yield_strength = ys;
    this->density = den;
    this->color = std::move(col);
}

float Material::shear_modulus() {
    return this->elastic_modulus / (2 * (1 + this->poissons_ratio));
}

Eigen::ArrayXXi create_mesh(
        Eigen::ArrayXXd points,
        Eigen::ArrayXXd facets,
        Eigen::ArrayXXd holes,
        Eigen::ArrayXXd control_points,
        std::vector<float> mesh_sizes,
        bool coarse
) {
    std::vector<float> mesh_size;

    if (mesh_sizes.size() > control_points.size())
        throw "Mesh size is bigger than control points";

    if (mesh_sizes.size() == 1) {
        mesh_size = std::vector<float>(control_points.size(), mesh_sizes[0]);
    }
//    std::map<std::string, std::vector<std::vector<float>>> tri;
//    tri["vertices"] = points;
//    tri["segments"] = facets; //set facets
//    if (holes.size() > 0) { tri["holes"] = holes; }  // set holes
//    std::vector<std::vector<float>> regions;
//    for (int i = 0; i < control_points.size(); i++) {
//        std::vector<float> cp = control_points[i];
//        std::vector<float> tmp;
//        tmp.push_back(cp[0]);
//        tmp.push_back(cp[1]);
//        tmp.push_back(i);
//        tmp.push_back(mesh_size[i]);
//        regions.push_back(tmp);
//    }
//    tri["regions"] = regions;

    auto mesh = distmesh::triangulation::delaunay(points);
    return mesh;
}
