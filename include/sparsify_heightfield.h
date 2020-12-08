 #ifndef SPARSIFY_HEIGHTFIELD
#define SPARSIFY_HEIGHTFIELD
#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>
void sparsify_heightfield(const Eigen::MatrixXd & X, const Eigen::MatrixXd & Y, const Eigen::MatrixXd & Z0, Eigen::MatrixXd & Z);
#endif
