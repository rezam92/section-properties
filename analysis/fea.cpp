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

    float area = 0;
    float qx = 0;
    float qy = 0;
    float ixx = 0;
    float iyy = 0;
    float ixy = 0;

    //# Gauss points for 6 point Gaussian integration
    Eigen::MatrixXf gps = gauss_points(6);

    //# loop through each Gauss point
    for (int i = 0; i < gps.rows(); ++i) {
        Eigen::MatrixXf gp = gps.row(i);
        //(N, _, j) = shape_function(self.coords, gp)
        tuple<Eigen::MatrixXf, Eigen::MatrixXf, float> res = shape_function(this->coords, gp);
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
tuple<Eigen::MatrixXf, Eigen::MatrixXf> Tri6::torsion_properties() {
    //# initialise stiffness matrix and load vector
    Eigen::MatrixXf k_el;
    Eigen::MatrixXf f_el;

    //# Gauss points for 6 point Gaussian integration
    Eigen::MatrixXf gps = gauss_points(6);

    for (int i = 0; i < gps.rows(); ++i) {
        Eigen::MatrixXf gp = gps.row(i);
        tuple<Eigen::MatrixXf, Eigen::MatrixXf, float> res = shape_function(this->coords, gp);
        float j = get<2>(res);
        Eigen::MatrixXf B = get<1>(res);
        Eigen::MatrixXf N = get<0>(res);

        float Nx = (N * this->coords.row(0).transpose())(0);
        float Ny = (N * this->coords.row(1).transpose())(0);
        Eigen::MatrixXf mat(2, 1);
        mat << Ny, -Nx;
        k_el += (gp(0) * (B.transpose() * B) * j * (this->material.elastic_modulus));
        f_el += (gp(0) * (B.transpose() * mat.transpose()) * j * (this->material.elastic_modulus));
    }

    tuple<Eigen::MatrixXf, Eigen::MatrixXf> tup(k_el, f_el);
    return tup;
}

/// Calculates the element shear load vectors used to evaluate the shear functions.
/// \param ixx : Second moment of area about the centroidal x-axis
/// \param iyy : Second moment of area about the centroidal y-axis
/// \param ixy : Second moment of area about the centroidal xy-axis
/// \param nu : Effective Poisson's ratio for the cross-section
/// \return Element shear load vector psi *(f_psi)* and phi *(f_phi)*
tuple<Eigen::MatrixXf, Eigen::MatrixXf> Tri6::shear_load_vectors(float ixx, float iyy, float ixy, float nu) {

    Eigen::MatrixXf f_psi;
    Eigen::MatrixXf f_phi;


    //# Gauss points for 6 point Gaussian integration
    Eigen::MatrixXf gps = gauss_points(6);

    for (int i = 0; i < gps.rows(); ++i) {
        Eigen::MatrixXf gp = gps.row(i);
        //# determine shape function, shape function derivative and jacobian
        tuple<Eigen::Matrix<float, 1, 6>, Eigen::Matrix<float, 2, 6>, float> res = shape_function(this->coords, gp);
        float j = get<2>(res);
        Eigen::MatrixXf B = get<1>(res);
        Eigen::MatrixXf N = get<0>(res);

        //# determine x and y position at Gauss point
        float Nx = (N * this->coords.row(0).transpose())(0);
        float Ny = (N * this->coords.row(1).transpose())(0);

        //# determine shear parameters
        float r = Nx * Nx - Ny * Ny;
        float q = 2 * Nx * Ny;
        float d1 = ixx * r - ixy * q;
        float d2 = ixy * r + ixx * q;
        float h1 = -ixy * r + iyy * q;
        float h2 = -iyy * r - ixy * q;

        Eigen::Matrix<float, 1, 2> tmp;
        Eigen::Matrix<float, 1, 2> tmp2;
        tmp << d1, d2;
        tmp2 << h1, h2;

        f_psi = gp(0) * ((nu / 2 * (B.transpose() * tmp.transpose()).transpose()) +
                         (2 * (1 + nu) * (N.transpose() * (ixx * Nx - ixy * Ny)))) * j * this->material.elastic_modulus;

        f_phi = gp(0) * ((nu / 2 * (B.transpose() * tmp2.transpose()).transpose()) +
                         (2 * (1 + nu) * (N.transpose() * (iyy * Ny - ixy * Nx)))) * j * this->material.elastic_modulus;

        tuple<Eigen::MatrixXf, Eigen::MatrixXf> tup(f_psi, f_phi);
        return tup;
    }

}

/// """Calculates the element shear centre and warping integrals required for shear analysis of
/// \param ixx Second moment of area about the centroidal x-axis
/// \param iyy Second moment of area about the centroidal y-axis
/// \param ixy Second moment of area about the centroidal xy-axis
/// \param omega Values of the warping function at the element nodes
/// \return hear centre integrals about the x and y-axes *(sc_xint, sc_yint)*, warping integrals *(q_omega, i_omega, i_xomega, i_yomega)*
tuple<float, float, float, float, float, float>
Tri6::shear_warping_integrals(float ixx, float iyy, float ixy, Eigen::MatrixXf omega) {

    float sc_xint = 0;
    float sc_yint = 0;
    float q_omega = 0;
    float i_omega = 0;
    float i_xomega = 0;
    float i_yomega = 0;

//# Gauss points for 6 point Gaussian integration
    Eigen::MatrixXf gps = gauss_points(6);

    for (int i = 0; i < gps.rows(); ++i) {
        Eigen::MatrixXf gp = gps.row(i);
        //# determine shape function, shape function derivative and jacobian
        tuple<Eigen::Matrix<float, 1, 6>, Eigen::Matrix<float, 2, 6>, float> res = shape_function(this->coords, gp);
        float j = get<2>(res);
        Eigen::MatrixXf B = get<1>(res);
        Eigen::MatrixXf N = get<0>(res);

        //# determine x and y position at Gauss point
        float Nx = (N * this->coords.row(0).transpose())(0);
        float Ny = (N * this->coords.row(1).transpose())(0);
        Eigen::MatrixXf Nomega = N * omega;

        sc_xint = gp(0) * (iyy * Nx + ixy * Ny) * (pow(Nx, 2) + pow(Ny, 2)) * j * this->material.elastic_modulus;
        sc_yint = gp(0) * (ixx * Ny + ixy * Nx) * (pow(Nx, 2) + pow(Ny, 2)) * j * this->material.elastic_modulus;

        q_omega += gp(0) * Nomega(0) * j * this->material.elastic_modulus;
        i_omega += gp(0) * pow(Nomega(0), 2) * j * this->material.elastic_modulus;
        i_xomega += gp(0) * Nx * Nomega(0) * j * this->material.elastic_modulus;
        i_yomega += gp(0) * Ny * Nomega(0) * j * this->material.elastic_modulus;
    }
    return tuple(sc_xint, sc_yint, q_omega, i_omega, i_xomega, i_yomega);
}

/// """Calculates the variables used to determine the shear deformation coefficients.
/// \param ixx: Second moment of area about the centroidal x-axis
/// \param iyy: Second moment of area about the centroidal y-axis
/// \param ixy: Second moment of area about the centroidal xy-axis
/// \param psi_shear: Values of the psi shear function at the element nodes
/// \param phi_shear: Values of the phi shear function at the element nodes
/// \param nu Effective Poisson's ratio for the cross-section
/// \return Shear deformation variables *(kappa_x, kappa_y, kappa_xy)*
tuple<float, float, float>
Tri6::shear_coefficients(float ixx, float iyy, float ixy, Eigen::MatrixXf psi_shear, Eigen::MatrixXf phi_shear,
                         float nu) {
    float kappa_x = 0;
    float kappa_y = 0;
    float kappa_xy = 0;

    Eigen::MatrixXf gps = gauss_points(6);
    for (int i = 0; i < gps.rows(); ++i) {
        Eigen::MatrixXf gp = gps.row(i);
        //# determine shape function, shape function derivative and jacobian
        tuple<Eigen::Matrix<float, 1, 6>, Eigen::Matrix<float, 2, 6>, float> res = shape_function(this->coords, gp);
        float j = get<2>(res);
        Eigen::MatrixXf B = get<1>(res);
        Eigen::MatrixXf N = get<0>(res);

        float Nx = (N * this->coords.row(0).transpose())(0);
        float Ny = (N * this->coords.row(1).transpose())(0);

        float r = Nx * Nx - Ny * Ny;
        float q = 2 * Nx * Ny;
        float d1 = ixx * r - ixy * q;
        float d2 = ixy * r + ixx * q;
        float h1 = -ixy * r + iyy * q;
        float h2 = -iyy * r - ixy * q;

        Eigen::Matrix<float, 1, 2> tmp;
        Eigen::Matrix<float, 1, 2> tmp2;
        tmp << d1, d2;
        tmp2 << h1, h2;

        kappa_x += (
                gp(0)
                * (psi_shear * B.transpose() - nu / 2 * tmp) * (B * (psi_shear) - nu / 2 * tmp)
                * j
                * this->material.elastic_modulus
        )(0);
        kappa_y += (
                gp(0)
                * (phi_shear * B.transpose() - nu / 2 * tmp2) * (B * (phi_shear) - nu / 2 * tmp2)
                * j
                * this->material.elastic_modulus
        )(0);
        kappa_xy += (
                gp(0)
                * (psi_shear * B.transpose() - nu / 2 * tmp) * (B * (phi_shear) - nu / 2 * tmp2)
                * j
                * this->material.elastic_modulus
        )(0);
    }
    return tuple(kappa_x, kappa_y, kappa_xy);
}

/// Calculates the integrals used to evaluate the monosymmetry constant about both global axes and both principal axes.
/// \param phi : Principal bending axis angle
/// \return Integrals used to evaluate the monosymmetry constants *(int_x, int_y, int_11, int_22)*
tuple<float, float, float, float> Tri6::monosymmetry_integrals(float phi) {

    //# initialise properties
    float int_x = 0;
    float int_y = 0;
    float int_11 = 0;
    float int_22 = 0;

    Eigen::MatrixXf gps = gauss_points(6);

    for (int i = 0; i < gps.rows(); ++i) {
        Eigen::MatrixXf gp = gps.row(i);
        //# determine shape function, shape function derivative and jacobian
        tuple<Eigen::Matrix<float, 1, 6>, Eigen::Matrix<float, 2, 6>, float> res = shape_function(this->coords, gp);
        float j = get<2>(res);
        Eigen::MatrixXf B = get<1>(res);
        Eigen::MatrixXf N = get<0>(res);

        float Nx = (N * this->coords.row(0).transpose())(0);
        float Ny = (N * this->coords.row(1).transpose())(0);

        //# determine 11 and 22 position at Gauss point
        tuple<float, float> pc = principal_coordinate(phi, Nx, Ny);
        float Nx_11 = get<0>(pc);
        float Ny_22 = get<1>(pc);

        //# weight the monosymmetry integrals by the section elastic modulus
        int_x += (
                gp(0)
                * (Nx * Nx * Ny + Ny * Ny * Ny)
                * j
                * this->material.elastic_modulus
        );
        int_y += (
                gp(0)
                * (Ny * Ny * Nx + Nx * Nx * Nx)
                * j
                * this->material.elastic_modulus
        );
        int_11 += (
                gp(0)
                * (Nx_11 * Nx_11 * Ny_22 + Ny_22 * Ny_22 * Ny_22)
                * j
                * this->material.elastic_modulus
        );
        int_22 += (
                gp(0)
                * (Ny_22 * Ny_22 * Nx_11 + Nx_11 * Nx_11 * Nx_11)
                * j
                * this->material.elastic_modulus
        );

    }

    return tuple(int_x, int_y, int_11, int_22);
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
tuple<Eigen::MatrixXf, Eigen::MatrixXf, Eigen::MatrixXf, Eigen::MatrixXf, Eigen::MatrixXf, Eigen::MatrixXf, Eigen::MatrixXf,
        Eigen::MatrixXf, Eigen::MatrixXf, Eigen::MatrixXf, Eigen::MatrixXf, Eigen::MatrixXf>
Tri6::element_stress(float N, float Mxx, float Myy, float M11, float M22, float Mzz, float Vx, float Vy, float ea,
                     float cx, float cy, float ixx, float iyy, float ixy, float i11, float i22, float phi, float j,
                     float nu, Eigen::MatrixXf omega, Eigen::MatrixXf psi_shear, Eigen::MatrixXf phi_shear,
                     float Delta_s) {
    Eigen::Matrix<float, 6, 1> un{
            1, 1, 1, 1, 1, 1
    };

    Eigen::MatrixXf sig_zz_n = N * un * this->material.elastic_modulus / ea;
    Eigen::Matrix<float, 6, 1> sig_zz_mxx_gp{0, 0, 0, 0, 0, 0};
    Eigen::Matrix<float, 6, 1> sig_zz_myy_gp{0, 0, 0, 0, 0, 0};
    Eigen::Matrix<float, 6, 1> sig_zz_m11_gp{0, 0, 0, 0, 0, 0};
    Eigen::Matrix<float, 6, 1> sig_zz_m22_gp{0, 0, 0, 0, 0, 0};
    Eigen::Matrix<float, 6, 2> sig_zxy_mzz_gp{{0, 0},
                                              {0, 0},
                                              {0, 0},
                                              {0, 0},
                                              {0, 0},
                                              {0, 0}};
    Eigen::Matrix<float, 6, 2> sig_zxy_vx_gp{{0, 0},
                                             {0, 0},
                                             {0, 0},
                                             {0, 0},
                                             {0, 0},
                                             {0, 0}};
    Eigen::Matrix<float, 6, 2> sig_zxy_vy_gp{{0, 0},
                                             {0, 0},
                                             {0, 0},
                                             {0, 0},
                                             {0, 0},
                                             {0, 0}};

    Eigen::MatrixXf gps = gauss_points(6);
    for (int i = 0; i < gps.rows(); ++i) {
        Eigen::Matrix<float, 2, 6> coords_c{{0, 0, 0, 0, 0, 0},
                                            {0, 0, 0, 0, 0, 0}};

        coords_c.row(0) = (this->coords.row(0).array() - cx).matrix();
        coords_c.row(1) = (this->coords.row(1).array() - cx).matrix();

        Eigen::MatrixXf gp = gps.row(i);

        //# determine shape function, shape function derivative and jacobian
        tuple<Eigen::Matrix<float, 1, 6>, Eigen::Matrix<float, 2, 6>, float> res = shape_function(this->coords, gp);
        Eigen::MatrixXf B = get<1>(res);
        Eigen::MatrixXf N = get<0>(res);

        float Nx = (N * this->coords.row(0).transpose())(0);
        float Ny = (N * this->coords.row(1).transpose())(0);

        tuple<float, float> pc = principal_coordinate(phi, Nx, Ny);
        float Nx_11 = get<0>(pc);
        float Ny_22 = get<1>(pc);

        float r = Nx * Nx - Ny * Ny;
        float q = 2 * Nx * Ny;
        float d1 = ixx * r - ixy * q;
        float d2 = ixy * r + ixx * q;
        float h1 = -ixy * r + iyy * q;
        float h2 = -iyy * r - ixy * q;

        sig_zz_mxx_gp(i, 0) = this->material.elastic_modulus * (
                -(ixy * Mxx) / (ixx * iyy - ixy * ixy) * Nx
                + (iyy * Mxx) / (ixx * iyy - ixy * ixy) * Ny
        );
        sig_zz_myy_gp(i, 0) = this->material.elastic_modulus * (
                -(ixx * Myy) / (ixx * iyy - ixy * ixy) * Nx
                + (ixy * Myy) / (ixx * iyy - ixy * ixy) * Ny
        );
        sig_zz_m11_gp(i, 0) = this->material.elastic_modulus * M11 / i11 * Ny_22;
        sig_zz_m22_gp(i, 0) = this->material.elastic_modulus * -M22 / i22 * Nx_11;


        if (Mzz != 0) {
            Eigen::MatrixXf tmp(2, 1);
            tmp << Ny, -Nx;
            sig_zxy_mzz_gp.row(i) = (
                    this->material.elastic_modulus * Mzz / j * (B * (omega) - tmp)
            );
        }
        if (Vy != 0) {
            Eigen::MatrixXf tmp(2, 1);
            tmp << d1, d2;
            sig_zxy_mzz_gp.row(i) = (
                    this->material.elastic_modulus
                    * Vx
                    / Delta_s
                    * (B * (phi_shear) - nu / 2 * tmp)
            );
        }
        if (Vx != 0) {
            Eigen::MatrixXf tmp(2, 1);
            tmp << h1, h2;
            sig_zxy_mzz_gp.row(i) = (
                    this->material.elastic_modulus
                    * Vy
                    / Delta_s
                    * (B * (phi_shear) - nu / 2 * tmp)
            );
        }
    }

    //# extrapolate results to nodes
    Eigen::Matrix<float, 1, 6> sig_zz_mxx = extrapolate_to_nodes(sig_zz_mxx_gp.col(0));
    Eigen::Matrix<float, 1, 6> sig_zz_myy = extrapolate_to_nodes(sig_zz_myy_gp.col(0));
    Eigen::Matrix<float, 1, 6> sig_zz_m11 = extrapolate_to_nodes(sig_zz_m11_gp.col(0));
    Eigen::Matrix<float, 1, 6> sig_zz_m22 = extrapolate_to_nodes(sig_zz_m22_gp.col(0));
    Eigen::Matrix<float, 1, 6> sig_zx_mzz = extrapolate_to_nodes(sig_zxy_mzz_gp.col(0));
    Eigen::Matrix<float, 1, 6> sig_zy_mzz = extrapolate_to_nodes(sig_zxy_mzz_gp.col(1));
    Eigen::Matrix<float, 1, 6> sig_zx_vx = extrapolate_to_nodes(sig_zxy_vx_gp.col(0));
    Eigen::Matrix<float, 1, 6> sig_zy_vx = extrapolate_to_nodes(sig_zxy_vx_gp.col(1));
    Eigen::Matrix<float, 1, 6> sig_zx_vy = extrapolate_to_nodes(sig_zxy_vy_gp.col(0));
    Eigen::Matrix<float, 1, 6> sig_zy_vy = extrapolate_to_nodes(sig_zxy_vy_gp.col(1));


    return tuple(
            sig_zz_n, // Eigen::MatrixXf
            sig_zz_mxx, //Eigen::Matrix<float, 1, 6>
            sig_zz_myy, // Eigen::Matrix<float, 1, 6>
            sig_zz_m11, // Eigen::Matrix<float, 1, 6>
            sig_zz_m22, // Eigen::Matrix<float, 1, 6>
            sig_zx_mzz, // Eigen::Matrix<float, 1, 6>
            sig_zy_mzz, //Eigen::Matrix<float, 1, 6>
            sig_zx_vx, //Eigen::Matrix<float, 1, 6>
            sig_zy_vx, // Eigen::Matrix<float, 1, 6>
            sig_zx_vy,// Eigen::Matrix<float, 1, 6>
            sig_zy_vy, // Eigen::Matrix<float, 1, 6>
            gps.col(0) // Eigen::MatrixXf
    );
}

/// Determines whether a point lies within the current element.
/// \param pt Point to check *(x, y)*
/// \return Whether the point lies within an element

bool Tri6::point_within_element(Eigen::ArrayX<float> pt) {
    float px = pt(0);
    float py = pt(1);

    //# get coordinates of corner points
    float x1 = this->coords(0, 0);
    float y1 = this->coords(1, 0);
    float x2 = this->coords(0, 1);
    float y2 = this->coords(1, 1);
    float x3 = this->coords(0, 2);
    float y3 = this->coords(1, 2);


    //# compute variables alpha, beta and gamma
    float alpha = ((y2 - y3) * (px - x3) + (x3 - x2) * (py - y3)) / (
            (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)
    );
    float beta = ((y3 - y1) * (px - x3) + (x1 - x3) * (py - y3)) / (
            (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)
    );
    float gamma = 1.0 - alpha - beta;

    //# if the point lies within an element
    if (alpha >= 0 and beta >= 0 and gamma >= 0) { return true; }
    else { return false; }
}


// --- define function

/// Returns the Gaussian weights and locations for *n* point Gaussian integration of a quadratic triangular element.
/// \param n: Number of Gauss points (1, 3 or 6)
/// \return An *n x 4* matrix consisting of the integration weight and the eta, xi and zeta locations for *n* Gauss points
Eigen::MatrixXf gauss_points(int n) {

    if (n == 1) {
        Eigen::MatrixXf out(1, 4);
        out << 1, 1 / 3.0, 1 / 3.0, 1 / 3.0;
        return out;
    } else if (n == 3) {
        Eigen::MatrixXf out(3, 4);
        out << 1.0 / 3, 2.0 / 3, 1.0 / 6, 1.0 / 6,
                1.0 / 3, 1.0 / 6, 2.0 / 3, 1.0 / 6,
                1.0 / 3, 1.0 / 6, 1.0 / 6, 2.0 / 3;
        return out;
    } else if (n == 6) {
        Eigen::MatrixXf out(6, 4);

        float g1 = 1.0 / 18 * (8 - sqrt(10) + sqrt(38 - 44 * sqrt(2.0 / 5)));
        float g2 = 1.0 / 18 * (8 - sqrt(10) - sqrt(38 - 44 * sqrt(2.0 / 5)));
        float w1 = (620 + sqrt(213125 - 53320 * sqrt(10))) / 3720;
        float w2 = (620 - sqrt(213125 - 53320 * sqrt(10))) / 3720;

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
tuple<Eigen::Matrix<float, 1, 6>, Eigen::Matrix<float, 2, 6>, float>
shape_function(const Eigen::Matrix<float, 2, 6> coords, Eigen::MatrixXf gauss_point) {

    // cordes 2x6
    float eta = gauss_point(0);
    float xi = gauss_point(1);
    float zeta = gauss_point(2);

    Eigen::MatrixXf N(1, 6);
    N << eta * (2 * eta - 1),
            xi * (2 * xi - 1),
            zeta * (2 * zeta - 1),
            4 * eta * xi,
            4 * xi * zeta,
            4 * eta * zeta;

    Eigen::MatrixXf B_iso(3, 6);
    B_iso << 4 * eta - 1, 0, 0, 4 * xi, 0, 4 * zeta,
            0, 4 * xi - 1, 0, 4 * eta, 4 * zeta, 0,
            0, 0, 4 * zeta - 1, 0, 4 * xi, 4 * eta;

    Eigen::MatrixXf J_upper(1, 3);
    J_upper << 1, 1, 1;

    Eigen::MatrixXf J_lower;
    // 2x6 --- 6x3 = 2x3
    J_lower = coords * B_iso.transpose();

    Eigen::MatrixXf J(3, 3);

    J.row(0) = J_upper;
    J.row(1) = J_lower.row(0);
    J.row(2) = J_lower.row(1);

    float j = 0.5 * (J.matrix()).determinant();
    Eigen::MatrixXf B;
    if (j != 0) {
        //        # calculate the P matrix
        Eigen::MatrixXf mul(3, 2);
        mul << 0, 0,
                1, 0,
                0, 1;
        Eigen::MatrixXf P = (J).inverse() * mul;

        //        # calculate the B matrix in terms of cartesian co-ordinates
        B = (B_iso.transpose() * P).transpose();
    } else {
        B = Eigen::MatrixXf::Zero(2, 6);
    }


    tuple<Eigen::MatrixXf, Eigen::MatrixXf, float> ret(N, B, j);

    return ret;
}

/// Extrapolates results at six Gauss points to the six nodes of a quadratic triangular element.
/// \param w Result at the six Gauss points [1 x 6]
/// \return Extrapolated nodal values at the six nodes [1 x 6]
Eigen::Matrix<float, 1, 6> extrapolate_to_nodes(Eigen::Matrix<float, 1, 6> w) {

    Eigen::MatrixXf H_inv(6, 6);
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
    float phi_rad = phi * M_PI / 180;

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
    float phi_rad = phi * M_PI / 180;

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
bool point_above_line(Eigen::Matrix<float, 1, 2> u, float px, float py, float x, float y) {
    Eigen::Matrix<float, 1, 3> PQ;
    PQ << px - x, py - y, 0;

    Eigen::Matrix<float, 1, 3> tmp;
    tmp(0) = u(0);
    tmp(1) = u(1);
    tmp(2) = 0;

    Eigen::MatrixXf res = PQ.cross(tmp);
    return res(0) > 0;
}
