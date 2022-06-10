//
// Created by Reza on 6/10/22.
//

#include "pre.h"

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

void create_mesh(
        std::vector<std::vector<float>> points,
        std::vector<std::vector<float>> facets,
        std::vector<std::vector<float>> holes,
        std::vector<std::vector<float>> control_points,
        std::variant<std::vector<float>, float> mesh_sizes,
        bool coarse
) {
    //
}
