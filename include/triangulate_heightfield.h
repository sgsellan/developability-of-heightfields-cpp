#ifndef TRIANGULATE_HEIGHTFIELD
#define TRIANGULATE_HEIGHTFIELD
#include <Eigen/Core>

void triangulate_heightfield(const Eigen::MatrixXd & X,const Eigen::MatrixXd & Y,const Eigen::MatrixXd & Z,Eigen::MatrixXd & U,Eigen::MatrixXi & G);
#endif
