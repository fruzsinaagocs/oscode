#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

namespace RKWKB{
    // complex numbers, vectors, matrices
    typedef Eigen::VectorXcd Vector;
    typedef Vector::Scalar Scalar;
    typedef Eigen::MatrixXcd Matrix;
    typedef Eigen::Tensor<Scalar, 3> Rank3T;
    typedef Eigen::Tensor<Scalar, 4> Rank4T;
    
    // complex-valued functions
    typedef std::function<Scalar(Vector)> Scalarfn;
    typedef std::function<Vector(Vector)> Vectorfn;
    typedef std::function<Matrix(Vector)> Matrixfn;
    typedef std::function<Rank3T(Vector)> Rank3Tfn;
    typedef std::function<Rank4T(Vector)> Rank4Tfn;
    typedef std::function<Scalar(Vector, Scalar)> event;
}
