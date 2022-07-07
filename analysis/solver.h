//
// Created by Reza on 6/23/22.
//

#ifndef SECTION_PROPERTIES_CPP_SOLVER_H
#define SECTION_PROPERTIES_CPP_SOLVER_H

#include "../Eigen/Core"

/// Solves a linear system of equations (Ku = f) using the CGS iterative method.
/// \param k : N x N matrix of the linear system
/// \param f : N x 1 right hand side of the linear system
/// \param m : Preconditioner for the linear matrix approximating the inverse of k
/// \param tol : Tolerance for the solver to achieve. The algorithm terminates when either the relative or the absolute residual is below tol.
/// \return The solution vector to the linear system of equations
Eigen::MatrixXf solve_cgs(Eigen::MatrixXf k, Eigen::MatrixXf f, m=None, float tol=1e-5);

/// Solves a linear system of equations (Ku = f) using the CGS iterative method and the Lagrangian multiplier method.
/// \param k_lg (N+1) x (N+1) Lagrangian multiplier matrix of the linear system
/// \param f N x 1 right hand side of the linear system
/// \param m Preconditioner for the linear matrix approximating the inverse of k
/// \param tol Tolerance for the solver to achieve. The algorithm terminates when either the relative or the absolute residual is below tol.
/// \return The solution vector to the linear system of equations
Eigen::MatrixXf solve_cgs_lagrange(Eigen::MatrixXf k_lg, Eigen::MatrixXf f, float tol=1e-5, m=None);

/// Solves a linear system of equations (Ku = f) using the direct solver method.
/// \param k N x N matrix of the linear system
/// \param f N x 1 right hand side of the linear system
/// \return
Eigen::MatrixXf solve_direct(Eigen::MatrixXf k, Eigen::MatrixXf f);

/// Solves a linear system of equations (Ku = f) using the direct solver method and the Lagrangian multiplier method.
/// \param k_lg (N+1) x (N+1) Lagrangian multiplier matrix of the linear system
/// \param f N x 1 right hand side of the linear system
/// \return The solution vector to the linear system of equations
Eigen::MatrixXf solve_direct_lagrange(Eigen::MatrixXf k_lg, Eigen::MatrixXf f);

//def create_progress()
//class CustomTimeElapsedColumn {};

#endif //SECTION_PROPERTIES_CPP_SOLVER_H
