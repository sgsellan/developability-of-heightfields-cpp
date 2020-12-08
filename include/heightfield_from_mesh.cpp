#include <Eigen/Core>
#include <igl/Hit.h>
#include <igl/AABB.h>
#include <igl/ray_mesh_intersect.h>



void heightfield_from_mesh(const Eigen::MatrixXd & V,const Eigen::MatrixXi & F, const int grid_size, Eigen::MatrixXd & X,Eigen::MatrixXd & Y,Eigen::MatrixXd & Z0)
{
	Eigen::MatrixXd V_normalized;
	int n = V.rows();
	V_normalized.resize(n,3);
    V_normalized.col(0) = V.col(0) - Eigen::VectorXd::Constant(n,V.col(0).minCoeff());
	V_normalized.col(1) = V.col(1) - Eigen::VectorXd::Constant(n,V.col(1).minCoeff());
	V_normalized.col(2) = V.col(2) - Eigen::VectorXd::Constant(n,V.col(2).minCoeff());
    V_normalized = V_normalized/(1.3*V_normalized.maxCoeff());
    V_normalized.col(0) = V_normalized.col(0) + Eigen::VectorXd::Constant(n,0.1);
    V_normalized.col(1) = V_normalized.col(1) + Eigen::VectorXd::Constant(n,0.1);
    V_normalized.col(2) = V_normalized.col(2) + Eigen::VectorXd::Constant(n,0.1);
	// V normalized done
	// Now, let's create a hexagonal grid of the unit square
	Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(grid_size,0,1);
	Eigen::VectorXd y = Eigen::VectorXd::LinSpaced(grid_size/sin(M_PI/3),0,1);
	double hx = x(1)-x(0);
	double hy = y(1)-y(0);
	X.resize(y.size(),x.size());	
	Y.resize(y.size(),x.size());
	for(int i = 0; i<x.size(); i++){
        if(i%2==0){
            X.col(i) = y;
        }else{
            X.col(i) = y - Eigen::VectorXd::Constant(y.size(),(hy/2));
        }
	}		
	for(int j = 0; j<y.size(); j++){
		Y.row(j) = x;
	}
	// Built the grid, now do ray intersect
//	VectorXd Xvec(Map<VectorXd>(X.data(), X.cols()*X.rows())) ,Yvec(Map<VectorXd>(Y.data(), Y.cols()*Y.rows()));
    Z0.resize(X.rows(),X.cols());
    igl::AABB<Eigen::MatrixXd,3> aabb;
    aabb.init(V_normalized,F);
    
	igl::Hit hit;
	Eigen::Vector3d origin,direction;
	direction << 0,0,-1;
	origin(2) = 1;
	for(int i = 0; i<X.rows(); i++){
		for(int j = 0; j<X.cols(); j++){
			origin(0) = X(i,j);
			origin(1) = Y(i,j);
            if(aabb.intersect_ray(V_normalized,F,origin,direction,hit)){
                Z0(i,j) = 2.0-hit.t;
            }else{
                Z0(i,j) = 0.0;
            }
            //std::cout << hit.t << std::endl;
		}
	}
	
}	
