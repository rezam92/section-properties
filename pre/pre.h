//
// Created by Reza on 6/10/22.
//

#ifndef SECTION_PROPERTIES_CPP_PRE_H
#define SECTION_PROPERTIES_CPP_PRE_H

#include "vector"
#include "string"
#include "variant"

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


void create_mesh(
            std::vector<std::vector<float>> points,
            std::vector<std::vector<float>> facets,
            std::vector<std::vector<float>> holes,
            std::vector<std::vector<float>> control_points,
            std::variant<std::vector<float>, float> mesh_sizes,
            bool coarse
        );

#endif //SECTION_PROPERTIES_CPP_PRE_H
