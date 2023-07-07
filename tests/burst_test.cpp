#include <oscode/solver.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <fstream>
#include <string>
#include <stdlib.h>

/** A variable to control the number of oscillations in the solution of the burst equation  */
double n = 100.0; 

/** Defines the friction term in the ODE to solve
 * @param[in] t Value of the independent variable in the ODE
 * @returns The value of the friction term at \a t
 */
std::complex<double> g(double t){
    return 0.0;
};

/** Defines the frequency term in the ODE to solve
 * @param[in] t Value of the independent variable in the ODE
 * @returns The value of the frequency term at \a t
 */
std::complex<double> w(double t){
    return std::pow(n*n - 1.0,0.5)/(1.0 + t*t);
};

/** Defines the initial value of the solution of the ODE
 * @param[in] t Value of the independent variable in the ODE
 * @returns The value of the solution \a x at \a t
 */
std::complex<double> xburst(double t){
    return 100*std::pow(1.0 + t*t,
    0.5)/n*std::complex<double>(std::cos(n*std::atan(t)),std::sin(n*std::atan(t))); 
};

/** Defines the initial value of the derivative of the solution of the ODE
 * @param[in] t Value of the independent variable in the ODE
 * @returns The value of the derivative of solution \a dx/dt at \a t
 */
std::complex<double> dxburst(double t){
    return 100/std::pow(1.0 + t*t,
    0.5)/n*(std::complex<double>(t,n)*std::cos(n*std::atan(t)) +
    std::complex<double>(-n,t)*std::sin(n*std::atan(t))); 
};

template <typename T>
std::vector<T> linspace(T start, T end, std::size_t points) {
  std::vector<T> res(points);
  float step = (end - start) / (points - 1);
  size_t i = 0;
  for (auto& e : res) {
    e = start + step * i++;
  }
  return res;
}

TEST(SolverTest, SolveBurst) {
    /** File to write solution of ODE to */
    std::string output = "output.txt"; 
   
    /** Define integration range */
    double ti = -2*n;
    double tf = 2*n;

    /** Define initial conditions */
    std::complex<double> x0 = xburst(ti); 
    std::complex<double> dx0 = dxburst(ti); 

    /** Create differential equation "system" */
    /** Method 1: Give frequency and damping term as functions */
    de_system sys(&w, &g);
    std::vector<double> times = linspace(ti, tf, 5000);
    /** Solve the ODE */    
    Solution solution(sys, x0, dx0, ti, tf, times);
    solution.solve();

    EXPECT_TRUE(true);
}

/** \brief Routine to call oscode to solve an example ODE, to illustrate basic
 * functionality.
 *
 * Routine to call oscode to solve the example ODE (the "burst" equation):
 * \f[
 * \ddot{x} + \fract{n^2-1}{(1+t^2)^2}x = 0
 * \f]
 * where \f$ n \f$ is a parameter that controls how many oscillations the solution
 * \f$ x(t) \f$ will go through.
 *
 * The routine defines the differential equation in four different ways:
 * -# Defines the frequency and friction terms (\f$ \omega(t) \f$ and \f$
 *  \gamma(t) \f$) as functions,
 * -# Defines the frequency and friction terms as sequences:
 *      -# as Eigen vectors, 
 *      -# as std::arrays,
 *      -# as std::vectors,
 * and then feeds a pointer to the underlying data array to the class \a de_system representing the ODE.
 *
 * All four methods are included in this routine, with the last three commented
 * out. Test any of the above methods by uncommenting and commenting out all others.
 *
 * After solving the ODE, this routine prints the solution and its derivative to
 * a file called \a output.txt, the contents of which you can plot with the
 * Python script \a plot_burst.py, or any other script of your choice.
 */

