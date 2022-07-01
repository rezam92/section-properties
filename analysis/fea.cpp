//
// Created by Reza on 6/23/22.
//

#include "fea.h"
#include "../Eigen/Geometry"

//TO-DO: define
//Tri6::Tri6() = default;


///  Calculates the geometric properties for the current finite element.
/// \return Tuple containing the geometric properties and the elastic and shear moduli of the element: *(area, qx, qy, ixx, iyy, ixy, e, g, rho)*
tuple<float, float, float, float, float, float, float, float, float> Tri6::geometric_properties() {

    //# initialise geometric properties

    double area = 0;
    double qx = 0;
    double qy = 0;
    double ixx = 0;
    double iyy = 0;
    double ixy = 0;

    //# Gauss points for 6 point Gaussian integration
    Eigen::MatrixXd gps = gauss_points(6);

    //# loop through each Gauss point
    for (int i = 0; i < gps.rows(); ++i) {
        Eigen::MatrixXd gp = gps.row(i);
        //(N, _, j) = shape_function(self.coords, gp)
        tuple<Eigen::MatrixXd, Eigen::MatrixXd, float> res = shape_function(this->coords, gp);
        float j = get<2>(res);
        area += gp(0) * j;

        qx += gp(0) * (get<0>(res) * this->coords.row(1).transpose())(0) * j;
        qy += gp(0) * (get<0>(res) * this->coords.row(0).transpose())(0) * j;

        ixx += gp(0) * pow((get<0>(res) * this->coords.row(1).transpose())(0), 2) * j;
        iyy += gp(0) * pow((get<0>(res) * this->coords.row(0).transpose())(0), 2) * j;

        ixy += (gp(0) * ((get<0>(res) * this->coords.row(1).transpose())(0)) *
                ((get<0>(res) * this->coords.row(0).transpose())(0)) * j);
    }

    tuple<float, float, float, float, float, float, float, float, float> tup(
            area, qx, qy, ixx, iyy, ixy, this->material.elastic_modulus,
            this->material.shear_modulus(), this->material.density
    );

    return tup;
}

/// Calculates the element stiffness matrix used for warping analysis and the torsion load vector.
/// \return Element stiffness matrix *(k_el)* and element torsion load vector *(f_el)*
tuple<Eigen::MatrixXd, Eigen::MatrixXd> Tri6::torsion_properties() {
    //# initialise stiffness matrix and load vector
    Eigen::MatrixXd k_el;
    Eigen::MatrixXd f_el;

    //# Gauss points for 6 point Gaussian integration
    Eigen::MatrixXd gps = gauss_points(6);

    for (int i = 0; i < gps.rows(); ++i) {
        Eigen::MatrixXd gp = gps.row(i);
        tuple<Eigen::MatrixXd, Eigen::MatrixXd, float> res = shape_function(this->coords, gp);
        float j = get<2>(res);
        Eigen::MatrixXd B = get<1>(res);
        Eigen::MatrixXd N = get<0>(res);

        double Nx = (N * this->coords.row(0).transpose())(0);
        double Ny = (N * this->coords.row(1).transpose())(0);
        Eigen::MatrixXd mat(1, 2);
        mat << Ny, -Nx;
        k_el += (gp(0) * (B.transpose() * B) * j * (this->material.elastic_modulus));
        f_el += (gp(0) * (B.transpose() * mat.transpose()) * j * (this->material.elastic_modulus));
    }

    tuple<Eigen::MatrixXd, Eigen::MatrixXd> tup(k_el, f_el);
    return tup;
}

/// Calculates the element shear load vectors used to evaluate the shear functions.
/// \param ixx : Second moment of area about the centroidal x-axis
/// \param iyy : Second moment of area about the centroidal y-axis
/// \param ixy : Second moment of area about the centroidal xy-axis
/// \param nu : Effective Poisson's ratio for the cross-section
/// \return Element shear load vector psi *(f_psi)* and phi *(f_phi)*
tuple<Eigen::MatrixXd, Eigen::MatrixXd> Tri6::shear_load_vectors(float ixx, float iyy, float ixy, float nu) {

    Eigen::MatrixXd f_psi;
    Eigen::MatrixXd f_phi;


    //# Gauss points for 6 point Gaussian integration
    Eigen::MatrixXd gps = gauss_points(6);

    for (int i = 0; i < gps.rows(); ++i) {
        Eigen::MatrixXd gp = gps.row(i);
        //# determine shape function, shape function derivative and jacobian
        tuple<Eigen::Matrix<double, 1, 6> , Eigen::Matrix<double, 2, 6>, float> res = shape_function(this->coords, gp);
        float j = get<2>(res);
        Eigen::MatrixXd B = get<1>(res);
        Eigen::MatrixXd N = get<0>(res);

        //# determine x and y position at Gauss point
        double Nx = (N * this->coords.row(0).transpose())(0);
        double Ny = (N * this->coords.row(1).transpose())(0);

        //# determine shear parameters
        double r = Nx * Nx - Ny * Ny;
        double q = 2 * Nx * Ny;
        double d1 = ixx * r - ixy * q;
        double d2 = ixy * r + ixx * q;
        double h1 = -ixy * r + iyy * q;
        double h2 = -iyy * r - ixy * q;

        Eigen::Matrix<double, 2, 1> tmp;
        Eigen::Matrix<double, 2, 1> tmp2;
        tmp << d1, d2;
        tmp2 << h1, h2;

        f_psi = gp(0) * ((nu / 2 * (B.transpose() * tmp.transpose()).transpose()) + (2 * (1 + nu) * (N.transpose() * (ixx * Nx - ixy * Ny)))) * j * this->material.elastic_modulus;

        f_phi = gp(0) * ((nu / 2 * (B.transpose() * tmp2.transpose()).transpose()) + (2 * (1 + nu) * (N.transpose() * (iyy * Ny - ixy * Nx)))) * j * this->material.elastic_modulus;

        tuple<Eigen::MatrixXd, Eigen::MatrixXd> tup(f_psi, f_phi);
        return tup;
    }

}

/// """Calculates the element shear centre and warping integrals required for shear analysis of
/// \param ixx Second moment of area about the centroidal x-axis
/// \param iyy Second moment of area about the centroidal y-axis
/// \param ixy Second moment of area about the centroidal xy-axis
/// \param omega Values of the warping function at the element nodes
/// \return hear centre integrals about the x and y-axes *(sc_xint, sc_yint)*, warping integrals *(q_omega, i_omega, i_xomega, i_yomega)*
//TODO: define
tuple<float, float, float, float, float, float>
Tri6::shear_warping_integrals(float ixx, float iyy, float ixy, Eigen::MatrixXd omega) {

}

/// """Calculates the variables used to determine the shear deformation coefficients.
/// \param ixx: Second moment of area about the centroidal x-axis
/// \param iyy: Second moment of area about the centroidal y-axis
/// \param ixy: Second moment of area about the centroidal xy-axis
/// \param psi_shear: Values of the psi shear function at the element nodes
/// \param phi_shear: Values of the phi shear function at the element nodes
/// \param nu Effective Poisson's ratio for the cross-section
/// \return Shear deformation variables *(kappa_x, kappa_y, kappa_xy)*
//TODO: define

tuple<float>
Tri6::shear_coefficients(float ixx, float iyy, float ixy, Eigen::MatrixXd psi_shear, Eigen::MatrixXd phi_shear,
                         float nu) {

}

/// Calculates the integrals used to evaluate the monosymmetry constant about both global axes and both principal axes.
/// \param phi : Principal bending axis angle
/// \return Integrals used to evaluate the monosymmetry constants *(int_x, int_y, int_11, int_22)*
//TODO: define
tuple<float, float, float, float> Tri6::monosymmetry_integrals(float phi) {

}

/// Calculates the stress within an element resulting from a specified loading. Also returns the shape function weights.
/// \param N : Axial force
/// \param Mxx : Bending moment about the centroidal xx-axis
/// \param Myy : Bending moment about the centroidal yy-axis
/// \param M11 : Bending moment about the centroidal 11-axis
/// \param M22 : Bending moment about the centroidal 22-axis
/// \param Mzz : Torsion moment about the centroidal zz-axis
/// \param Vx : Shear force acting in the x-direction
/// \param Vy : Shear force acting in the y-direction
/// \param ea : Modulus weighted area
/// \param cx : x position of the elastic centroid
/// \param cy : y position of the elastic centroid
/// \param ixx : Second moment of area about the centroidal x-axis
/// \param iyy : Second moment of area about the centroidal y-axis
/// \param ixy : Second moment of area about the centroidal xy-axis
/// \param i11 : Second moment of area about the principal 11-axis
/// \param i22 : Second moment of area about the principal 22-axis
/// \param phi : Principal bending axis angle
/// \param j : St. Venant torsion constant
/// \param nu : Effective Poisson's ratio for the cross-section
/// \param omega : Values of the warping function at the element nodes
/// \param psi_shear : Values of the psi shear function at the element nodes
/// \param phi_shear : Values of the phi shear function at the element nodes
/// \param Delta_s : Cross-section shear factor
/// \return Tuple containing element stresses and integration weights
///            (:math:`sigma_{zz,n}`, :math:`sigma_{zz,mxx}`,
///            :math:`sigma_{zz,myy}`, :math:`sigma_{zz,m11}`,
///            :math:`sigma_{zz,m22}`, :math:`sigma_{zx,mzz}`,
///            :math:`sigma_{zy,mzz}`, :math:`sigma_{zx,vx}`,
///            :math:`sigma_{zy,vx}`, :math:`sigma_{zx,vy}`,
///            :math:`sigma_{zy,vy}`, :math:`w_i`)
//TODO: define
tuple<Eigen::MatrixXd>
Tri6::element_stress(float N, float Mxx, float Myy, float M11, float M22, float Mzz, float Vx, float Vy, float ea,
                     float cx, float cy, float ixx, float iyy, float ixy, float i11, float i22, float phi, float j,
                     float nu, Eigen::MatrixXd omega, Eigen::MatrixXd psi_shear, Eigen::MatrixXd phi_shear,
                     float Delta_s) {

}

/// Determines whether a point lies within the current element.
/// \param pt Point to check *(x, y)*
/// \return Whether the point lies within an element
//TODO: define
bool Tri6::point_within_element(Eigen::ArrayX<float> pt) {
    return false;
}


// --- define function

/// Returns the Gaussian weights and locations for *n* point Gaussian integration of a quadratic triangular element.
/// \param n: Number of Gauss points (1, 3 or 6)
/// \return An *n x 4* matrix consisting of the integration weight and the eta, xi and zeta locations for *n* Gauss points
Eigen::MatrixXd gauss_points(int n) {

    if (n == 1) {
        Eigen::MatrixXd out(1, 4);
        out << 1, 1 / 3.0, 1 / 3.0, 1 / 3.0;
        return out;
    } else if (n == 3) {
        Eigen::MatrixXd out(3, 4);
        out << 1.0 / 3, 2.0 / 3, 1.0 / 6, 1.0 / 6,
                1.0 / 3, 1.0 / 6, 2.0 / 3, 1.0 / 6,
                1.0 / 3, 1.0 / 6, 1.0 / 6, 2.0 / 3;
        return out;
    } else if (n == 6) {
        Eigen::MatrixXd out(6, 4);

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

/// Computes shape functions, shape function derivatives and the determinant of the Jacobian matrix for a tri 6 element at a given Gauss point.
/// \param coords Global coordinates of the quadratic triangle vertices [2 x 6]
/// \param gauss_point Gaussian weight and isoparametric location of the Gauss point
/// \return The value of the shape functions *N(i)* at the given Gauss point [1 x 6],
/// the derivative of the shape functions in the j-th global direction *B(i,j)* [2 x 6] and the determinant of the Jacobian matrix *j*
tuple<Eigen::Matrix<double, 1, 6> , Eigen::Matrix<double, 2, 6>, float>
shape_function(const Eigen::Matrix<double, 2, 6> coords, Eigen::MatrixXd gauss_point) {

    // cordes 2x6
    double eta = gauss_point(0);
    double xi = gauss_point(1);
    double zeta = gauss_point(2);

    Eigen::MatrixXd N(1, 6);
    N << eta * (2 * eta - 1),
            xi * (2 * xi - 1),
            zeta * (2 * zeta - 1),
            4 * eta * xi,
            4 * xi * zeta,
            4 * eta * zeta;

    Eigen::MatrixXd B_iso(3, 6);
    B_iso << 4 * eta - 1, 0, 0, 4 * xi, 0, 4 * zeta,
            0, 4 * xi - 1, 0, 4 * eta, 4 * zeta, 0,
            0, 0, 4 * zeta - 1, 0, 4 * xi, 4 * eta;

    Eigen::MatrixXd J_upper(1, 3);
    J_upper << 1, 1, 1;

    Eigen::MatrixXd J_lower;
    // 2x6 --- 6x3 = 2x3
    J_lower = coords * B_iso.transpose();

    Eigen::MatrixXd J(3, 3);

    J.row(0) = J_upper;
    J.row(1) = J_lower.row(0);
    J.row(2) = J_lower.row(1);

    double j = 0.5 * (J.matrix()).determinant();
    Eigen::MatrixXd B;
    if (j != 0) {
        //        # calculate the P matrix
        Eigen::MatrixXd mul(3, 2);
        mul << 0, 0,
                1, 0,
                0, 1;
        Eigen::MatrixXd P = (J).inverse() * mul;

        //        # calculate the B matrix in terms of cartesian co-ordinates
        B = (B_iso.transpose() * P).transpose();
    } else {
        B = Eigen::MatrixXd::Zero(2, 6);
    }


    tuple<Eigen::MatrixXd, Eigen::MatrixXd, float> ret(N, B, j);

    return ret;
}

/// Extrapolates results at six Gauss points to the six nodes of a quadratic triangular element.
/// \param w Result at the six Gauss points [1 x 6]
/// \return Extrapolated nodal values at the six nodes [1 x 6]
Eigen::Matrix<double, 1, 6> extrapolate_to_nodes(Eigen::Matrix<double, 1, 6> w) {

    Eigen::MatrixXd H_inv(6, 6);
    H_inv
            << 1.87365927351160, 0.138559587411935, 0.138559587411935, -0.638559587411936, 0.126340726488397, -0.638559587411935,
            0.138559587411935, 1.87365927351160, 0.138559587411935, -0.638559587411935, -0.638559587411935, 0.126340726488397,
            0.138559587411935, 0.138559587411935, 1.87365927351160, 0.126340726488396, -0.638559587411935, -0.638559587411935,
            0.0749010751157440, 0.0749010751157440, 0.180053080734478, 1.36051633430762, -0.345185782636792, -0.345185782636792,
            0.180053080734478, 0.0749010751157440, 0.0749010751157440, -0.345185782636792, 1.36051633430762, -0.345185782636792,
            0.0749010751157440, 0.180053080734478, 0.0749010751157440, -0.345185782636792, -0.345185782636792, 1.36051633430762;

    return (H_inv * w.transpose()).transpose();
}

/// Determines the coordinates of the cartesian point *(x, y)* in the principal axis system given an axis rotation angle phi.
/// \param phi Principal bending axis angle (degrees)
/// \param x x coordinate in the global axis
/// \param y y coordinate in the global axis
/// \return Principal axis coordinates *(x1, y2)*
tuple<float, float> principal_coordinate(float phi, float x, float y) {
    double phi_rad = phi * M_PI / 180;

    Eigen::ArrayXXd R;

    R << cos(phi_rad), sin(phi_rad), -sin(phi_rad), cos(phi_rad);

    Eigen::ArrayXXd arr;
    arr << x, y;
    Eigen::ArrayXXd x_rotated = R * arr;

    return tuple(x_rotated(0), x_rotated(1));

}


/// Determines the global coordinates of the principal axis point *(x1, y2)* given principal axis rotation angle phi.
/// \param phi Principal bending axis angle (degrees)
/// \param x11 11 coordinate in the principal axis
/// \param y22 22 coordinate in the principal axis
/// \return Global axis coordinates *(x, y)*
tuple<float, float> global_coordinate(float phi, float x11, float y22) {
    double phi_rad = phi * M_PI / 180;

    Eigen::ArrayXXd R;

    R << cos(phi_rad), -sin(phi_rad), sin(phi_rad), cos(phi_rad);

    Eigen::ArrayXXd arr;
    arr << x11, y22;
    Eigen::ArrayXXd x_rotated = R * arr;

    return tuple(x_rotated(0), x_rotated(1));
}

/// Determines whether a point *(x, y)* is a above or below the line defined by the parallel unit vector *u* and the point *(px, py)*.
/// \param u  Unit vector parallel to the line [1 x 2]
/// \param px x coordinate of a point on the line
/// \param py y coordinate of a point on the line
/// \param x x coordinate of the point to be tested
/// \param y y coordinate of the point to be tested
/// \return This method returns *True* if the point is above the line or *False* if the point is below the line
bool point_above_line(Eigen::Matrix<double, 1, 2> u, float px, float py, float x, float y) {
    Eigen::Matrix<double, 1, 3> PQ;
    PQ << px - x, py - y, 0;

    Eigen::Matrix<double, 1, 3> tmp;
    tmp(0) = u(0);
    tmp(1) = u(1);
    tmp(2) = 0;

    Eigen::MatrixXd res = PQ.cross(tmp);
    return res(0) > 0;
}
