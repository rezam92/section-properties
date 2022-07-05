//
// Created by Reza on 6/23/22.
//

#ifndef SECTION_PROPERTIES_CPP_FEA_H
#define SECTION_PROPERTIES_CPP_FEA_H

#include "../Eigen/Core"
#include "../pre/pre.h"
#include "tuple"

using namespace std;

class Tri6 {
public:
//    Tri6();
    tuple<float,float,float,float,float,float,float,float,float> geometric_properties();

    tuple<Eigen::MatrixXf, Eigen::MatrixXf> torsion_properties();

    tuple<Eigen::MatrixXf, Eigen::MatrixXf> shear_load_vectors(float ixx, float iyy, float ixy, float nu);

    tuple<float, float, float, float, float, float>
    shear_warping_integrals(float ixx, float iyy, float ixy, Eigen::MatrixXf omega);

    tuple<float, float, float>
    shear_coefficients(float ixx, float iyy, float ixy, Eigen::MatrixXf psi_shear, Eigen::MatrixXf phi_shear, float nu);

    tuple<float, float, float, float> monosymmetry_integrals(float phi);

    tuple<Eigen::MatrixXf,Eigen::MatrixXf,Eigen::MatrixXf,Eigen::MatrixXf,Eigen::MatrixXf,Eigen::MatrixXf,Eigen::MatrixXf,
            Eigen::MatrixXf,Eigen::MatrixXf,Eigen::MatrixXf,Eigen::MatrixXf,Eigen::MatrixXf> element_stress(
            float N, float Mxx, float Myy, float M11, float M22,
            float Mzz, float Vx, float Vy, float ea, float cx,
            float cy, float ixx, float iyy, float ixy, float i11,
            float i22, float phi, float j, float nu,
            Eigen::MatrixXf omega, Eigen::MatrixXf psi_shear, Eigen::MatrixXf phi_shear,
            float Delta_s
    );

    bool point_within_element(Eigen::ArrayX<float> pt);

    int el_id = 0;
    Eigen::MatrixXf coords;
    Eigen::MatrixXf node_ids;
    Material material = DEFAULT_MATERIAL; //Material("mat", 0, 0, 0, 0, "mat");
private:
};

Eigen::MatrixXf gauss_points(int n);

tuple<Eigen::Matrix<float, 1, 6> , Eigen::Matrix<float, 2, 6>, float>
shape_function(const Eigen::MatrixXf &coords, Eigen::MatrixXf gauss_point);

Eigen::Matrix<float, 1, 6> extrapolate_to_nodes(tuple<Eigen::MatrixXf> w);

tuple<float, float> principal_coordinate(float phi, float x, float y);

tuple<float, float> global_coordinate(float phi, float x11, float y22);

bool point_above_line(tuple<Eigen::MatrixXf> u, float px, float py, float x, float y);

#endif //SECTION_PROPERTIES_CPP_FEA_H
