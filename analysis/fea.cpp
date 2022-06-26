//
// Created by Reza on 6/23/22.
//

#include "fea.h"

//TODO: define
//Tri6::Tri6() = default;

//TODO: define
tuple<float> Tri6::geometric_properties() {
    return tuple<float>();
}

//TODO: define
tuple<Eigen::ArrayXXd> Tri6::torsion_properties() {
    return tuple<Eigen::ArrayXXd>();
}

//TODO: define
tuple<Eigen::ArrayXXd> Tri6::shear_load_vectors(float ixx, float iyy, float ixy, float nu) {
    return tuple<Eigen::ArrayXXd>();
}

//TODO: define
tuple<float> Tri6::shear_warping_integrals(float ixx, float iyy, float ixy, Eigen::ArrayXXd omega) {
    return tuple<float>();
}

//TODO: define
tuple<float>
Tri6::shear_coefficients(float ixx, float iyy, float ixy, Eigen::ArrayXXd psi_shear, Eigen::ArrayXXd phi_shear,
                         float nu) {
    return tuple<float>();
}

//TODO: define
tuple<float> Tri6::monosymmetry_integrals(float phi) {
    return tuple<float>();
}

//TODO: define
tuple<Eigen::ArrayXXd>
Tri6::element_stress(float N, float Mxx, float Myy, float M11, float M22, float Mzz, float Vx, float Vy, float ea,
                     float cx, float cy, float ixx, float iyy, float ixy, float i11, float i22, float phi, float j,
                     float nu, Eigen::ArrayXXd omega, Eigen::ArrayXXd psi_shear, Eigen::ArrayXXd phi_shear,
                     float Delta_s) {
    return tuple<Eigen::ArrayXXd>();
}

//TODO: define
bool Tri6::point_within_element(Eigen::ArrayX<float> pt) {
    return false;
}


// --- define function

Eigen::ArrayXXd gauss_points(int n) {

    if (n == 1) {
        Eigen::ArrayXXd out(1, 4);
        out << 1, 1 / 3.0, 1 / 3.0, 1 / 3.0;
        return out;
    } else if (n == 3) {
        Eigen::ArrayXXd out(3, 4);
        out << 1.0 / 3, 2.0 / 3, 1.0 / 6, 1.0 / 6,
                1.0 / 3, 1.0 / 6, 2.0 / 3, 1.0 / 6,
                1.0 / 3, 1.0 / 6, 1.0 / 6, 2.0 / 3;
        return out;
    } else if (n == 6) {
        Eigen::ArrayXXd out(6, 4);

        double g1 = 1.0 / 18 * (8 - sqrt(10) + sqrt(38 - 44 * sqrt(2.0 / 5)));
        double g2 = 1.0 / 18 * (8 - sqrt(10) - sqrt(38 - 44 * sqrt(2.0 / 5)));
        double w1 = (620 + sqrt(213125 - 53320 * sqrt(10))) / 3720;
        double w2 = (620 - sqrt(213125 - 53320 * sqrt(10))) / 3720;

        out << w2, 1 - 2 * g2, g2, g2,
                w2, g2, 1 - 2 * g2, g2,
                w2, g2, g2, 1 - 2 * g2,
                w1, g1, g1, 1 - 2 * g1,
                w1, 1 - 2 * g1, g1, g1,
                w1, g1, 1 - 2 * g1, g1;
        return out;
    }
    throw "Invalid value of `n`";
}

tuple<Eigen::ArrayXXd, Eigen::ArrayXXd, float> shape_function(Eigen::ArrayXXd coords, Eigen::ArrayXXd gauss_point) {

    // cordes 2x6
    double eta = gauss_point(0);
    double xi = gauss_point(1);
    double zeta = gauss_point(2);

    Eigen::ArrayXXd N(1, 6);
    N << eta * (2 * eta - 1),
            xi * (2 * xi - 1),
            zeta * (2 * zeta - 1),
            4 * eta * xi,
            4 * xi * zeta,
            4 * eta * zeta;

    Eigen::ArrayXXd B_iso(3, 6);
    B_iso << 4 * eta - 1, 0, 0, 4 * xi, 0, 4 * zeta,
            0, 4 * xi - 1, 0, 4 * eta, 4 * zeta, 0,
            0, 0, 4 * zeta - 1, 0, 4 * xi, 4 * eta;

    Eigen::ArrayXXd J_upper(1, 3);
    J_upper << 1, 1, 1;

    Eigen::ArrayXXd J_lower;
    // 2x6 --- 6x3 = 2x3
    J_lower = coords * B_iso.transpose();

    Eigen::ArrayXXd J(3, 3);

    J.row(0) = J_upper;
    J.row(1) = J_lower.row(0);
    J.row(2) = J_lower.row(1);

    double j = 0.5 * (J.matrix()).determinant();
    Eigen::ArrayXXd B;
    if (j != 0) {
        //        # calculate the P matrix
        Eigen::ArrayXXd mul(3, 2);
        mul << 0, 0,
                1, 0,
                0, 1;
        Eigen::MatrixXd P = (J).inverse() * mul;

        //        # calculate the B matrix in terms of cartesian co-ordinates
        B = (B_iso.transpose() * P.array()).transpose();
    } else {
        B = Eigen::ArrayXXd::Zero(2, 6);
    }


    tuple<Eigen::ArrayXXd, Eigen::ArrayXXd, float> ret(N, B, j);

    return ret;
}

//TODO: define
tuple<Eigen::ArrayXXd> extrapolate_to_nodes(tuple<Eigen::ArrayXXd> w) {
    return tuple<Eigen::ArrayXXd>();
}

//TODO: define
tuple<float> principal_coordinate(float phi, float x, float y) {
    return tuple<float>();
}

//TODO: define
tuple<float> global_coordinate(float phi, float x11, float y22) {
    return tuple<float>();
}

//TODO: define
bool point_above_line(tuple<Eigen::ArrayXXd> u, float px, float py, float x, float y) {
    return false;
}
