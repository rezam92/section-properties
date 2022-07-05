//
// Created by Reza on 6/23/22.
//

#ifndef SECTION_PROPERTIES_CPP_SOLVER_H
#define SECTION_PROPERTIES_CPP_SOLVER_H

#include "../Eigen/Core"


Eigen::MatrixXf solve_cgs(Eigen::MatrixXf k, Eigen::MatrixXf f, m=None, float tol=1e-5);

Eigen::MatrixXf solve_cgs_lagrange(Eigen::MatrixXf k_lg, Eigen::MatrixXf f, float tol=1e-5, m=None);

Eigen::MatrixXf solve_direct(Eigen::MatrixXf k, Eigen::MatrixXf f);

Eigen::MatrixXf solve_direct_lagrange(Eigen::MatrixXf k_lg, Eigen::MatrixXf f);

//def create_progress()
//class CustomTimeElapsedColumn {};

#endif //SECTION_PROPERTIES_CPP_SOLVER_H
