#include <Eigen/Core>
#include <Eigen/Sparse>
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
                   Eigen::VectorXd & )> argmin_Z, const Eigen::SparseMatrix<double> A, const Eigen::SparseMatrix<double> B, const Eigen::VectorXd c, const Eigen::VectorXd & Z0, Eigen::VectorXd & Z)
{
    Eigen::VectorXd X,U,Xprev,Uprev,Zprev;
    U.resize(B.rows());
    U.setZero();
    X.resize(B.rows());
    X = -B*Z0;
    Z.resize(B.cols());
    Z = Z0;
    // intialized
    double rho = .000001;
    int iterations = 0;
    bool stop = false;
    double residual, dual_residual, tol_abs, tol_rel, eps_pri, eps_dual;
    tol_abs = 0.0001;
    tol_rel = 0.001;
    int check_interval = 10;
    double bmu = 5000;
    double btao_inc = 2;
    double btao_dec = 2;
    
    while (iterations<10000) {
        iterations = iterations + 1;
        Xprev = X;
        //std::cout << "test" << std::endl;
        argmin_X(Z,U,rho,X);
        Zprev = Z;
        //std::cout << "test" << std::endl;
        argmin_Z(X,U,rho,Z);
        Uprev = U;
        
        U = U + A*X + B*Z - c; //CAREFUL SIGN CONVENTION
        //std::cout << "test" << std::endl;
        residual = (A*X + B*Z - c).norm();
        dual_residual = (A.transpose()*(B*(Zprev-Z))).norm();
        int k = c.size();
        if (iterations%check_interval == 0) {
            if (residual > bmu*dual_residual) {
                rho = btao_inc*rho;
                U = U/btao_inc;
            }else if(dual_residual > bmu*residual){
                rho = rho/btao_dec;
                U = U*btao_dec;
            }
        }
        eps_pri = sqrt(k*2)*tol_abs + tol_rel*std::max((A*X).norm(),std::max((B*Z).norm(),c.norm()));
        eps_dual = sqrt(k)*tol_abs + tol_rel*rho*(A.transpose()*U).norm();
        if (residual < eps_pri && dual_residual < eps_dual) {
            std::cout << "converged" << std::endl;
            break;
        }
    }
    
}
