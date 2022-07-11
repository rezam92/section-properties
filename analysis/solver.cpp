//
// Created by Reza on 6/23/22.
//

#include "solver.h"


Eigen::MatrixXf solve_cgs(Eigen::MatrixXf k, Eigen::MatrixXf f/* m, float tol*/) {
    Eigen::VectorXf x(k.cols()), b(k.cols());
    Eigen::ConjugateGradient<Eigen::SparseMatrix<float>, Eigen::Upper> cg;
    cg.compute(k.sparseView());
    x = cg.solve(b);
    return x.matrix();
}

Eigen::MatrixXf solve_direct(Eigen::MatrixXf k, Eigen::MatrixXf f) {
    return k.householderQr().solve(f).matrix();
}

//Eigen::MatrixXf solve_cgs_lagrange(Eigen::MatrixXf k_lg, Eigen::MatrixXf f/*, float tol, m*/) {
//    return Eigen::MatrixXf();
//}
//


//Eigen::MatrixXf solve_direct_lagrange(Eigen::MatrixXf k_lg, Eigen::MatrixXf f) {
//    return Eigen::MatrixXf();
//}
