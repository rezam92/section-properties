//
// Created by Reza on 6/10/22.
//

#ifndef SECTION_PROPERTIES_CPP_PRE_H
#define SECTION_PROPERTIES_CPP_PRE_H

#include "vector"
#include "string"
#include "variant"
#include "map"
#include "../Eigen/Core"

class Material {
    std::string name;
    float elastic_modulus;
    float poissons_ratio;
    float yield_strength;
    float density;
    std::string color;

public:
    Material(std::string n, float  em, float pr, float ys, float den, std::string col);
    float  shear_modulus();
};


const Material DEFAULT_MATERIAL = Material("default", 1, 0, 1, 1, "w");


Eigen::ArrayXXi create_mesh(
            Eigen::ArrayXXd points,
            Eigen::ArrayXXd facets,
            Eigen::ArrayXXd holes,
            Eigen::ArrayXXd control_points,
            std::variant<std::vector<float>, float> mesh_sizes,
            bool coarse
        );

#endif //SECTION_PROPERTIES_CPP_PRE_H
