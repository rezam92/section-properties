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

    ///  Calculates the geometric properties for the current finite element.
/// \return Tuple containing the geometric properties and the elastic and shear moduli of the element: *(area, qx, qy, ixx, iyy, ixy, e, g, rho)*
    tuple<float, float, float, float, float, float, float, float, float> geometric_properties();

/// Calculates the element stiffness matrix used for warping analysis and the torsion load vector.
/// \return Element stiffness matrix *(k_el)* and element torsion load vector *(f_el)*
    tuple<Eigen::MatrixXf, Eigen::MatrixXf> torsion_properties();

    /// Calculates the element shear load vectors used to evaluate the shear functions.
/// \param ixx : Second moment of area about the centroidal x-axis
/// \param iyy : Second moment of area about the centroidal y-axis
/// \param ixy : Second moment of area about the centroidal xy-axis
/// \param nu : Effective Poisson's ratio for the cross-section
/// \return Element shear load vector psi *(f_psi)* and phi *(f_phi)*
    tuple<Eigen::MatrixXf, Eigen::MatrixXf> shear_load_vectors(float ixx, float iyy, float ixy, float nu);

/// Calculates the element shear centre and warping integrals required for shear analysis of
/// \param ixx Second moment of area about the centroidal x-axis
/// \param iyy Second moment of area about the centroidal y-axis
/// \param ixy Second moment of area about the centroidal xy-axis
/// \param omega Values of the warping function at the element nodes
/// \return hear centre integrals about the x and y-axes *(sc_xint, sc_yint)*, warping integrals *(q_omega, i_omega, i_xomega, i_yomega)*
    tuple<float, float, float, float, float, float>
    shear_warping_integrals(float ixx, float iyy, float ixy, Eigen::MatrixXf omega);

///Calculates the variables used to determine the shear deformation coefficients.
/// \param ixx: Second moment of area about the centroidal x-axis
/// \param iyy: Second moment of area about the centroidal y-axis
/// \param ixy: Second moment of area about the centroidal xy-axis
/// \param psi_shear: Values of the psi shear function at the element nodes
/// \param phi_shear: Values of the phi shear function at the element nodes
/// \param nu Effective Poisson's ratio for the cross-section
/// \return Shear deformation variables *(kappa_x, kappa_y, kappa_xy)*
    tuple<float, float, float>
    shear_coefficients(float ixx, float iyy, float ixy, Eigen::MatrixXf psi_shear, Eigen::MatrixXf phi_shear, float nu);

    /// Calculates the integrals used to evaluate the monosymmetry constant about both global axes and both principal axes.
/// \param phi : Principal bending axis angle
/// \return Integrals used to evaluate the monosymmetry constants *(int_x, int_y, int_11, int_22)*
    tuple<float, float, float, float> monosymmetry_integrals(float phi);

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
            Eigen::MatrixXf, Eigen::MatrixXf, Eigen::MatrixXf, Eigen::MatrixXf, Eigen::MatrixXf> element_stress(
            float N, float Mxx, float Myy, float M11, float M22,
            float Mzz, float Vx, float Vy, float ea, float cx,
            float cy, float ixx, float iyy, float ixy, float i11,
            float i22, float phi, float j, float nu,
            Eigen::MatrixXf omega, Eigen::MatrixXf psi_shear, Eigen::MatrixXf phi_shear,
            float Delta_s
    );

    /// Determines whether a point lies within the current element.
    /// \param pt Point to check *(x, y)*
    /// \return Whether the point lies within an element
    bool point_within_element(Eigen::ArrayX<float> pt);

    int el_id = 0;
    Eigen::MatrixXf coords;
    Eigen::MatrixXf node_ids;
    Material material = DEFAULT_MATERIAL; //Material("mat", 0, 0, 0, 0, "mat");
private:
};

/// Returns the Gaussian weights and locations for *n* point Gaussian integration of a quadratic triangular element.
/// \param n: Number of Gauss points (1, 3 or 6)
/// \return An *n x 4* matrix consisting of the integration weight and the eta, xi and zeta locations for *n* Gauss points
Eigen::MatrixXf gauss_points(int n);

/// Computes shape functions, shape function derivatives and the determinant of the Jacobian matrix for a tri 6 element at a given Gauss point.
/// \param coords Global coordinates of the quadratic triangle vertices [2 x 6]
/// \param gauss_point Gaussian weight and isoparametric location of the Gauss point
/// \return The value of the shape functions *N(i)* at the given Gauss point [1 x 6],
/// the derivative of the shape functions in the j-th global direction *B(i,j)* [2 x 6] and the determinant of the Jacobian matrix *j*
tuple<Eigen::Matrix<float, 1, 6>, Eigen::Matrix<float, 2, 6>, float>
shape_function(const Eigen::MatrixXf &coords, Eigen::MatrixXf gauss_point);

/// Extrapolates results at six Gauss points to the six nodes of a quadratic triangular element.
/// \param w Result at the six Gauss points [1 x 6]
/// \return Extrapolated nodal values at the six nodes [1 x 6]
Eigen::Matrix<float, 1, 6> extrapolate_to_nodes(tuple<Eigen::MatrixXf> w);

/// Determines the coordinates of the cartesian point *(x, y)* in the principal axis system given an axis rotation angle phi.
/// \param phi Principal bending axis angle (degrees)
/// \param x x coordinate in the global axis
/// \param y y coordinate in the global axis
/// \return Principal axis coordinates *(x1, y2)*
tuple<float, float> principal_coordinate(float phi, float x, float y);

/// Determines the global coordinates of the principal axis point *(x1, y2)* given principal axis rotation angle phi.
/// \param phi Principal bending axis angle (degrees)
/// \param x11 11 coordinate in the principal axis
/// \param y22 22 coordinate in the principal axis
/// \return Global axis coordinates *(x, y)*
tuple<float, float> global_coordinate(float phi, float x11, float y22);

/// Determines whether a point *(x, y)* is a above or below the line defined by the parallel unit vector *u* and the point *(px, py)*.
/// \param u  Unit vector parallel to the line [1 x 2]
/// \param px x coordinate of a point on the line
/// \param py y coordinate of a point on the line
/// \param x x coordinate of the point to be tested
/// \param y y coordinate of the point to be tested
/// \return This method returns *True* if the point is above the line or *False* if the point is below the line
bool point_above_line(tuple<Eigen::MatrixXf> u, float px, float py, float x, float y);

#endif //SECTION_PROPERTIES_CPP_FEA_H
