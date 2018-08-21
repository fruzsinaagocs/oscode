#include <eigen3/Eigen/Dense>

namespace RKWKB{
    // complex numbers, vectors, matrices
    typedef Eigen::VectorXcd Vector;
    typedef Vector::Scalar Scalar;
    typedef Eigen::MatrixXcd Matrix;
    
    // complex-valued functions
    typedef std::function<Scalar(Vector)> Scalarfn;
    typedef std::function<Vector(Vector)> Vectorfn;
    typedef std::function<Matrix(Vector)> Matrixfn;
    typedef std::function<Scalar(Vector, Scalar)> event;
}
