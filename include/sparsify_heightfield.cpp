#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SVD>
#include <iostream>
#include <vector>
#include <math.h>
#include <igl/parallel_for.h>
#include <igl/min_quad_with_fixed.h>
#include "admm.h"


void sparsify_heightfield(const Eigen::MatrixXd & X, const Eigen::MatrixXd & Y, const Eigen::MatrixXd & Z0, Eigen::MatrixXd & Z)
{
    using namespace Eigen;
    std::vector<int> II;
    Eigen::VectorXd Z0col;
    Z0col.resize(X.rows()*X.cols());
    int u_counter = 0;
    for(int i = 0; i<Z0.rows(); i++){
        for(int j = 0; j<Z0.cols(); j++){
            Z0col(u_counter) = Z0(i,j);
            u_counter = u_counter + 1;
        }
    }
    // Step 1 & 2: Build II and operator
    u_counter = 0;
    int ii_counter = 0;
    
    Eigen::MatrixXd B;
    double omega = 10000000;
    double h = X(0,1) - X(0,0);
    double sq = sqrt(3);
    B.resize(4,7);
    B <<    0,    0,    3,   -6,    3,    0,    0,
    2,    2,   -1,   -6,   -1,    2,    2,
    -2*sq, 2*sq,    0,    0,    0, 2*sq,-2*sq,
    -2*sq, 2*sq,    0,    0,    0, 2*sq,-2*sq; // hxy = hyx;
    // Full quadric fit:
    //           -h,    h, -2*h,    0,  2*h,   -h,   -h,
    //         h*sq, h*sq,    0,    0,    0,-h*sq,-h*sq,
    //            0,    0,    0, -6*(h^2),0,    0,    0;
    B = B/(6.0*(pow(h,2.0)));
    Eigen::VectorXi neighbors_odd, neighbors_even, neighbors;
    neighbors_odd.resize(7);
    neighbors_even.resize(7);
    neighbors.resize(7);
    neighbors_even << -1,X.cols()-1,-X.cols(),0,X.cols(),1,X.cols()+1;
    neighbors_odd << -X.cols()-1,-1,-X.cols(),0,X.cols(),-X.cols()+1,1;
    std::vector<Eigen::Triplet<double> > ijv, ijv_speye_II;
    std::vector<int> fixed;
    std::vector<double> fixed_values;
    bool valid = true;
    for(int i = 0; i<X.rows(); i++){
        for(int j = 0; j<X.cols(); j++){
            if (i>0 && i<(X.rows()-1) && j>0 && j<(X.cols()-1) && Z0(i,j)>0 ) {
                if (j%2==0) {
                    neighbors = neighbors_even;
                }else{
                    neighbors = neighbors_odd;
                }
                valid = true;
                for (int nn = 0; nn<7; nn = nn+1) {
                    if (Z0col(u_counter + neighbors(nn))==0.0) {
                        valid = false;
                    }
                }
                if (valid) {
                for (int row = 0; row<4; row = row + 1) {
                    for (int col = 0; col<7; col = col + 1) {
                        ijv.emplace_back(4*ii_counter + row,u_counter + neighbors(col),B(row,col));
                    }
                }
                II.push_back(u_counter); // Interior point
                ijv_speye_II.emplace_back(u_counter,u_counter,1.0);
                ii_counter = ii_counter + 1;
                }
            }
            if (Z0(i,j)==0) {
                fixed.push_back(u_counter);
                fixed_values.push_back(0.0);
            }
            Z0col(u_counter) = Z0(i,j);
            u_counter = u_counter + 1;
        }
    }
    Eigen::SparseMatrix<double> A,speye_II;
    speye_II.resize(Z0col.size(),Z0col.size());
    A.resize(4.*II.size(),Z0col.size());
    A.setFromTriplets(ijv.begin(),ijv.end());
    speye_II.setFromTriplets(ijv_speye_II.begin(),ijv_speye_II.end());
    
    // Step 3: Build two admm functions
    igl::min_quad_with_fixed_data<double>  precomputed_data;
    Eigen::SparseMatrix<double> speye(A.cols(),A.cols());
    speye.setIdentity();
    Eigen::SparseMatrix<double> Q,Aeq;
    double rho = .000001;
    double rho_prev = rho;
    Q = omega*speye_II + 0.001*speye + 0.5*rho*A.transpose()*A;
    Aeq.resize(0,0);
    Eigen::VectorXi known = Eigen::VectorXi::Map(fixed.data(),fixed.size());
    min_quad_with_fixed_precompute(Q,known,Aeq,true,precomputed_data);
    // ARGMIN_Z
    std::function<void(
                       const Eigen::VectorXd &,
                       const Eigen::VectorXd &,
                       const double,
                       Eigen::VectorXd & )> argmin_Z = [&A,&omega,&speye,&speye_II,&Z0col,&fixed,&fixed_values,&precomputed_data,&rho_prev](
                                                                                    const Eigen::VectorXd & X,
                                                                                    const Eigen::VectorXd & U,
                                                                                    const double rho,
                                                                                    Eigen::VectorXd & Z)->void{
                           //std::cout << "test" << std::endl;
                           Eigen::VectorXd C,sol;
                           C = X + U;
                           Eigen::VectorXd L;
                           L = -1.0*((Z0col.transpose()*(omega*speye_II + 0.001*speye )) + (0.5*rho*C.transpose()*A));
                           Eigen::VectorXd Beq;
                           Eigen::VectorXd known_values = Eigen::VectorXd::Map(fixed_values.data(),fixed_values.size());
                           Beq.resize(0);
                           if (rho!=rho_prev) {
                               Eigen::SparseMatrix<double> Q,Aeq;
                               Q = omega*speye_II + 0.001*speye + 0.5*rho*A.transpose()*A;
                               Aeq.resize(0,0);
                               Eigen::VectorXi known = Eigen::VectorXi::Map(fixed.data(),fixed.size());
                               min_quad_with_fixed_precompute(Q,known,Aeq,true,precomputed_data);
                           }
                           igl::min_quad_with_fixed_solve(precomputed_data,L,known_values,Beq,Z,sol);
                           rho_prev = rho;
                       };
    
    
    
    
    // ARGMIN_X (copied from the OG mex implementation)
    std::function<void(
                       const Eigen::VectorXd &,
                       const Eigen::VectorXd &,
                       const double,
                       Eigen::VectorXd & )> argmin_X = [&A](
                                                            const Eigen::VectorXd & Z,
                                                            const Eigen::VectorXd & U,
                                                            const double rho,
                                                            Eigen::VectorXd & X)->void{
                           Eigen::VectorXd C;
                           C = A*Z - U;
                           X.resize(C.size());
                           igl::parallel_for(C.size()/4,[&] (const int i){
                               MatrixXd MM(2,2);
                               MatrixXd Hval(2,2);
                               MatrixXd HH(2,2);
                               VectorXd S(2);
                               double rho_with_weight;
                               rho_with_weight = rho;
                               //rho = rho_with_weight;
                               double s1,h1,s2,h2;
                               HH(1,0) = 0;
                               HH(0,1) = 0;
                               
                               MM(0,0) = C(4*i);
                               MM(0,1) = C(4*i+1);
                               MM(1,0) = C(4*i+2);
                               MM(1,1) = C(4*i+3);
                               //
                               JacobiSVD<MatrixXd> svd( MM, ComputeFullV | ComputeFullU );
                               S = svd.singularValues();
                               
                               if ((S(0)-(1/rho_with_weight))>0) {
                                   h1 = S(0)-(1/rho_with_weight);
                               }else{
                                   if((S(0)+(1/rho_with_weight))<0){
                                       h1 = S(0)+(1/rho_with_weight);
                                   }else{
                                       h1 = 0;
                                   }
                               }
                               if ((S(1)-(1/rho_with_weight))>0) {
                                   h2 = S(1)-(1/rho_with_weight);
                               }else{
                                   if((S(1)+(1/rho_with_weight))<0){
                                       h2 = S(1)+(1/rho_with_weight);
                                   }else{
                                       h2 = 0;
                                   }
                               }
                               HH(0,0) = h1;
                               HH(1,1) = h2;
                               
                               Hval = svd.matrixU()*HH*svd.matrixV().transpose();
                               
                               X(4*i) = Hval(0,0);
                               X(4*i+1) = Hval(0,1);
                               X(4*i+2) = Hval(1,0);
                               X(4*i+3) = Hval(1,1);
                           },0);
                       };
    
    
    
    
    // Step 4: Call admm
    Eigen::VectorXd XX,U,Ztest,Ztest2;
    U.resize(A.rows());
    U.setZero();
    //std::cout << "test" << std::endl;
    XX = A*Z0col;
    Eigen::SparseMatrix<double> BB(A.rows(),A.rows());
    BB.setIdentity();
    Eigen::VectorXd c;
    c.resize(A.rows());
    c.setZero();
    // AX + BZ - c = 0 !!!
    admm(argmin_X,argmin_Z,BB,-A,c,Z0col,Ztest2);
    Z.resize(Z0.rows(),Z0.cols());
    //Z << Ztest2;
    Z = Z0; // Placeholder
    int counter = 0;
    for(int i = 0; i<Z.rows(); i++){
        for(int j = 0; j<Z.cols(); j++){
            Z(i,j) = Ztest2(counter);
            counter = counter + 1;
        }
    }
}
