#include <Eigen/Core>
#include <iostream>
#include <igl/resolve_duplicated_faces.h>

void triangulate_heightfield(const Eigen::MatrixXd & X,const Eigen::MatrixXd & Y,const Eigen::MatrixXd & Z,Eigen::MatrixXd & U, Eigen::MatrixXi & G)
{
    Eigen::MatrixXi G_clean, J;
    int u_counter = 0;
    int g_counter = 0;
    G.resize(2*(X.rows()-2)*(X.cols()-2),3);
    U.resize((X.rows())*(X.cols()),3);
    for(int i = 0; i<X.rows(); i++){
        for(int j = 0; j<X.cols(); j++){
            if (i>0 && i<(X.rows()-1) && j>0 && j<(X.cols()-1) ) {
                if (j%2==0) {
                    G(g_counter,0) = u_counter;
                    G(g_counter,1) = u_counter - 1;
                    G(g_counter,2) = u_counter + X.cols() - 1;
                    G(g_counter + 1,0) = u_counter - 1;
                    G(g_counter + 1,1) = u_counter;
                    G(g_counter + 1,2) = u_counter - X.cols();
                    
                }else{
                    G(g_counter,0) = u_counter;
                    G(g_counter,1) = u_counter - 1;
                    G(g_counter,2) = u_counter + X.cols();
                    G(g_counter + 1,0) = u_counter - 1;
                    G(g_counter + 1,1) = u_counter;
                    G(g_counter + 1,2) = u_counter - X.cols() - 1;
                }
                g_counter = g_counter + 2;
            }
            U(u_counter,0) = X(i,j);
            U(u_counter,1) = Y(i,j);
            U(u_counter,2) = Z(i,j);
            u_counter = u_counter + 1;
        }
    }
    
    // Sure, we could build the mesh in a complicated way and ensure that there
    // are no duplicates, but let's just do this for the moment.
    igl::resolve_duplicated_faces(G,G_clean,J);
    G = G_clean;
}
