#include "include/solver.hpp"
#include <boost/math/special_functions/airy.hpp>
#include "include/catch.hpp"

using namespace RKWKB;

Vector F(Vector z){
    Vector result(z.size());
    result << 0.25*std::pow(z(0),-3), 0.25*std::pow(z(1),-3);
    return result;
};

Matrix DF(Vector z){
    Matrix result(z.size(),z.size());
    result << -0.75*std::pow(z(0), -4), 0.0,
               0.0, -0.75*std::pow(z(1), -4);
    return result;

};

Scalar w(Vector z){
    return z(0)*z(1);
};

Vector Dw(Vector z){
    Vector result(z.size());
    result << z(1), z(0);
    return result;
}

Matrix DDw(Vector z){
    Matrix result(z.size(), z.size());
    result << 0.0, 1.0,
              1.0, 0.0;
    return result;
}

Scalar g(Vector z){
    return 0.0;
};

Vector Dg(Vector z){
    Vector result(z.size());
    result << 0.0, 0.0;
    return result;
};

Matrix DDg(Vector z){
    Matrix result(z.size(), z.size());
    result << 0.0, 0.0,
              0.0, 0.0;
    return result;
};

Scalar f_end(Vector y, Scalar t){
    // A function to end integration of the ODE when f_end=0.
    double t_f = 1500.0;
    return (t - t_f);
};
