//
// Created by Reza on 6/23/22.
//

#include "fea.h"

Tri6::Tri6() {

}

tuple<float> Tri6::geometric_properties() {
    return tuple<float>();
}

tuple<Eigen::ArrayXXd> Tri6::torsion_properties() {
    return tuple<Eigen::ArrayXXd>();
}

tuple<Eigen::ArrayXXd> Tri6::shear_load_vectors(float ixx, float iyy, float ixy, float nu) {
    return tuple<Eigen::ArrayXXd>();
}

tuple<float> Tri6::shear_warping_integrals(float ixx, float iyy, float ixy, Eigen::ArrayXXd omega) {
    return tuple<float>();
}

tuple<float>
Tri6::shear_coefficients(float ixx, float iyy, float ixy, Eigen::ArrayXXd psi_shear, Eigen::ArrayXXd phi_shear,
                         float nu) {
    return tuple<float>();
}

tuple<float> Tri6::monosymmetry_integrals(float phi) {
    return tuple<float>();
}

tuple<Eigen::ArrayXXd> Tri6::element_stress(float N, float Mxx, float Myy, float M11, float M22, float Mzz, float Vx, float Vy, float ea,
                     float cx, float cy, float ixx, float iyy, float ixy, float i11, float i22, float phi, float j,
                     float nu, Eigen::ArrayXXd omega, Eigen::ArrayXXd psi_shear, Eigen::ArrayXXd phi_shear,
                     float Delta_s) {
    return tuple<Eigen::ArrayXXd>();
}

bool Tri6::point_within_element(Eigen::ArrayX<float> pt) {
    return false;
}


// --- define function

Eigen::ArrayXXd gauss_points(int n) {
    return Eigen::ArrayXXd();
}

tuple<Eigen::ArrayXXd, float> shape_function(Eigen::ArrayXXd coords, Eigen::ArrayXXd gauss_point) {
    return tuple<Eigen::ArrayXXd, float>();
}

tuple<Eigen::ArrayXXd> extrapolate_to_nodes(tuple<Eigen::ArrayXXd> w) {
    return tuple<Eigen::ArrayXXd>();
}

tuple<float> principal_coordinate(float phi, float x, float y) {
    return tuple<float>();
}

tuple<float> global_coordinate(float phi, float x11, float y22) {
    return tuple<float>();
}

bool point_above_line(tuple<Eigen::ArrayXXd> u, float px, float py, float x, float y) {
    return false;
}
