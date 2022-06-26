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
    tuple<float> geometric_properties();
    tuple<Eigen::ArrayXXd> torsion_properties();
    tuple<Eigen::ArrayXXd> shear_load_vectors(float ixx, float iyy, float ixy, float nu);
    tuple<float> shear_warping_integrals(float ixx, float iyy, float ixy, Eigen::ArrayXXd omega);
    tuple<float> shear_coefficients(float ixx, float iyy, float ixy, Eigen::ArrayXXd psi_shear, Eigen::ArrayXXd phi_shear, float nu);
    tuple<float> monosymmetry_integrals(float phi);
    tuple<Eigen::ArrayXXd> element_stress(
            float N,float Mxx,float Myy,float M11,float M22,
            float Mzz,float Vx,float Vy,float ea,float cx,
            float cy,float ixx,float iyy,float ixy,float i11,
            float i22,float phi,float j,float nu,
            Eigen::ArrayXXd omega, Eigen::ArrayXXd psi_shear,Eigen::ArrayXXd phi_shear,
            float Delta_s
    );
    bool point_within_element(Eigen::ArrayX<float> pt);
private:
    int el_id = 0;
    Eigen::ArrayXXd coords;
    Eigen::ArrayXd node_ids;
    Material material = DEFAULT_MATERIAL; //Material("mat", 0, 0, 0, 0, "mat");
};

Eigen::ArrayXXd gauss_points(int n);

tuple<Eigen::ArrayXXd, float> shape_function(Eigen::ArrayXXd coords, Eigen::ArrayXXd gauss_point);

tuple<Eigen::ArrayXXd> extrapolate_to_nodes(tuple<Eigen::ArrayXXd> w);

tuple<float> principal_coordinate(float phi, float x, float y);

tuple<float> global_coordinate(float phi, float x11, float y22);

bool point_above_line(tuple<Eigen::ArrayXXd> u, float px, float py, float x, float y);

#endif //SECTION_PROPERTIES_CPP_FEA_H
