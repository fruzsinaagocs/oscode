#include <eigen3/Eigen/Dense>

// complex numbers, vectors, matrices
typedef Eigen::VectorXcd Vector;
typedef Vector::Scalar Scalar;
typedef Eigen::MatrixXcd Matrix;

// complex-valued functions
typedef Scalar (* Scalarfn)(Vector);
typedef Vector (* Vectorfn)(Vector);
typedef Matrix (* Matrixfn)(Vector);
