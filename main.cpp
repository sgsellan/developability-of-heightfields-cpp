#include <igl/opengl/glfw/Viewer.h>
#include "sparsify_heightfield.h"
#include "triangulate_heightfield.h"
#include "heightfield_from_mesh.h"
#include <iostream>

int main(int argc, char *argv[])
{
  Eigen::MatrixXd V,U,X,Y,Z0,Z,Vh;
  Eigen::MatrixXi F,G,Fh;
    std::string input = "data/range.obj";
    int grid_size = 100;
    int argindex = 0;
    while(argindex+1<argc){
        if(strncmp(argv[argindex+1],"-n",2)==0){
            grid_size = atoi(argv[argindex+2]);
            argindex = argindex+2;
        }else if(strncmp(argv[argindex+1],"-i",2)==0){
            input = argv[argindex+2];
            argindex = argindex+2;
        }
    }
    igl::readOBJ(input,V,F);
    Eigen::MatrixXd V_normalized;
    int n = V.rows();
    V_normalized.resize(n,3);
    V_normalized.col(0) = V.col(0) - Eigen::VectorXd::Constant(n,V.col(0).minCoeff());
    V_normalized.col(1) = V.col(1) - Eigen::VectorXd::Constant(n,V.col(1).minCoeff());
    V_normalized.col(2) = V.col(2) - Eigen::VectorXd::Constant(n,V.col(2).minCoeff());
    V_normalized = V_normalized/(1.3*V_normalized.maxCoeff());
    V_normalized.col(0) = V_normalized.col(0) + Eigen::VectorXd::Constant(n,-0.0);
    V_normalized.col(1) = V_normalized.col(1) + Eigen::VectorXd::Constant(n,-0.0);
    V_normalized.col(2) = V_normalized.col(2) + Eigen::VectorXd::Constant(n,+1.0);
    V = V_normalized;
    
    heightfield_from_mesh(V,F,grid_size,X,Y,Z0);
    triangulate_heightfield(X,Y,Z0,Vh,Fh);
    igl::opengl::glfw::Viewer viewer;
    viewer.callback_key_pressed = [&](decltype(viewer) &,unsigned int key, int mod)
    {
        switch(key)
        {
            case 'S': case 's':
                sparsify_heightfield(X,Y,Z0,Z);
                triangulate_heightfield(X,Y,Z,U,G);
                viewer.data().clear();
                viewer.data().set_mesh(U,G);
                return true;
            case 'H': case 'h':
                viewer.data().clear();
                viewer.data().set_mesh(Vh,Fh);
                return true;
        }
        return false;
    };
    std::cout<<R"(
H,h to turn into heightfield
S,s to make developable
    )";
  viewer.data().set_mesh(V,F);
  viewer.data().set_face_based(true);
  viewer.launch();
}
