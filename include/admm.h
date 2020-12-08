#ifndef ADMM
#define ADMM
#include <Eigen/Core>
#include <iostream>

void admm(const std::function<void(
                         const Eigen::VectorXd &,
                         const Eigen::VectorXd &,
                         const double,
                         Eigen::VectorXd & )> argmin_X,
const std::function<void(
                   const Eigen::VectorXd &,
                   const Eigen::VectorXd &,
                   const double,
                   Eigen::VectorXd & )> argmin_Z, const Eigen::SparseMatrix<double> A, const Eigen::SparseMatrix<double> B, const Eigen::VectorXd c, const Eigen::VectorXd & Z0, Eigen::VectorXd & Z);
#endif
