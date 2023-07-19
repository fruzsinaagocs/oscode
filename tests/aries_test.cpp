#include <oscode/solver.hpp>
#include <gtest/gtest.h>
#include <boost/math/special_functions/airy.hpp>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <string>

/** Defines the friction term in the ODE to solve
 * @param[in] t Value of the independent variable in the ODE
 * @returns The value of the friction term at \a t
 */
std::complex<double> g(double t) { return 0.0; }

/** Defines the frequency term in the ODE to solve
 * @param[in] t Value of the independent variable in the ODE
 * @returns The value of the frequency term at \a t
 */
std::complex<double> w(double t) { return std::sqrt(t); }

/** Defines the initial value of the solution of the ODE
 * @param[in] t Value of the independent variable in the ODE
 * @returns The value of the solution \a x at \a t
 */
std::complex<double> xairy(double t) {
  return std::complex<double>(boost::math::airy_ai(-t),
                              boost::math::airy_bi(-t));
}

/** Defines the initial value of the derivative of the solution of the ODE
 * @param[in] t Value of the independent variable in the ODE
 * @returns The value of the derivative of solution \a dx/dt at \a t
 */
std::complex<double> dxairy(double t) {
  return std::complex<double>(-boost::math::airy_ai_prime(-t),
                              -boost::math::airy_bi_prime(-t));
}

/** \brief Routine to call oscode to solve an example ODE, to illustrate basic
 * functionality.
 *
 * Routine to call oscode to solve the example ODE (the "burst" equation):
 * \f[
 * \ddot{x} + \fract{n^2-1}{(1+t^2)^2}x = 0
 * \f]
 * where \f$ n \f$ is a parameter that controls how many oscillations the
 * solution \f$ x(t) \f$ will go through.
 *
 * The routine defines the differential equation in four different ways:
 * -# Defines the frequency and friction terms (\f$ \omega(t) \f$ and \f$
 *  \gamma(t) \f$) as functions,
 * -# Defines the frequency and friction terms as sequences:
 *      -# as Eigen vectors,
 *      -# as std::arrays,
 *      -# as std::vectors,
 * and then feeds a pointer to the underlying data array to the class \a
 * de_system representing the ODE.
 *
 * All four methods are included in this routine, with the last three commented
 * out. Test any of the above methods by uncommenting and commenting out all
 * others.
 *
 * After solving the ODE, this routine prints the solution and its derivative to
 * a file called \a output.txt, the contents of which you can plot with the
 * Python script \a plot_burst.py, or any other script of your choice.
 */
TEST(SolverTest, SolveAriesUnEvenFwd) {
  /** Define integration range */
  double ti = 1.0;
  double tf = 1e6;

  /** Define initial conditions */
  std::complex<double> x0 = xairy(ti);
  std::complex<double> dx0 = dxairy(ti);

  /** We'll ask for dense output at the following points */

  double a = 0.1;
  double t_dense_i = 5e1;
  std::vector<double> t_dense(200);
  std::generate(
      t_dense.begin(), t_dense.end(),
      [n = 0, &a, &t_dense_i]() mutable { return n++ * a + t_dense_i; });

  /** Create differential equation "system" */
  /** Method 1: Give frequency and damping term as functions */
  de_system sys(&w, &g);
  /** Solve the ODE */
  Solution solution(sys, x0, dx0, ti, tf, t_dense, 3, 1e-8);
  solution.solve();
  std::vector<std::complex<double>> sol_vec;
  sol_vec.reserve(t_dense.size());
  for (auto &sol : solution.dosol) {
    sol_vec.push_back(sol);
  }
  std::vector<std::complex<double>> true_sol_vec;
  true_sol_vec.reserve(t_dense.size());
  for (auto &t : t_dense) {
    true_sol_vec.push_back(xairy(t));
  }
  std::vector<std::complex<double>> dsol_vec;
  for (auto &dsol : solution.dodsol) {
    dsol_vec.push_back(dsol);
  }
  std::vector<std::complex<double>> true_dsol_vec;
  true_dsol_vec.reserve(t_dense.size());
  for (auto &t : t_dense) {
    true_dsol_vec.push_back(dxairy(t));
  }
  for (std::size_t i = 0; i < t_dense.size(); ++i) {
    EXPECT_NEAR((std::real(sol_vec[i]) - std::real(true_sol_vec[i])), 0.0f,
                1e-4f);
    EXPECT_NEAR((std::imag(sol_vec[i]) - std::imag(true_sol_vec[i])), 0.0f,
                1e-4f);
    EXPECT_NEAR((std::real(dsol_vec[i]) - std::real(true_dsol_vec[i])), 0.0f,
                1e-4f);
    EXPECT_NEAR((std::imag(dsol_vec[i]) - std::imag(true_dsol_vec[i])), 0.0f,
                1e-4f);

  }
}
/*
TEST(SolverTest, SolveAriesUnEvenBwd) {
  double ti = 1.0;
  double tf = 1e6;
  std::complex<double> x0 = xairy(ti);
  std::complex<double> dx0 = dxairy(ti);

  double a = 0.1;
  double t_dense_i = 5e1;
  std::vector<double> t_dense(200);
  std::generate(
      t_dense.begin(), t_dense.end(),
      [n = 0, &a, &t_dense_i]() mutable { return n++ * a + t_dense_i; });
  std::reverse(t_dense.begin(), t_dense.end());

  de_system sys(&w, &g);
  Solution solution(sys, x0, dx0, tf, ti, t_dense, 3, 1e-8, 0, -0.1);
  solution.solve();
  std::vector<std::complex<double>> sol_vec;
  for (auto &sol : solution.dosol) {
    sol_vec.push_back(sol);
  }
  std::vector<std::complex<double>> true_sol_vec;
  for (auto &t : t_dense) {
    true_sol_vec.push_back(xairy(t));
  }
  std::vector<std::complex<double>> dsol_vec;
  for (auto &dsol : solution.dodsol) {
    dsol_vec.push_back(dsol);
  }

  std::vector<std::complex<double>> true_dsol_vec;
  for (auto &t : t_dense) {
    true_dsol_vec.push_back(dxairy(t));
  }
  for (std::size_t i = 0; i < t_dense.size(); ++i) {
    EXPECT_NEAR((std::real(sol_vec[i]) - std::real(true_sol_vec[i])), 0.0f,
                1e-3f)
        << "i = " << i << " sol = " << sol_vec[i]
        << " true_sol = " << true_sol_vec[i];
    EXPECT_NEAR((std::imag(sol_vec[i]) + std::imag(true_sol_vec[i])), 0.0f,
                1e-4f)
        << "i = " << i << " sol = " << sol_vec[i]
        << " true_sol = " << true_sol_vec[i];
    EXPECT_NEAR((std::real(dsol_vec[i]) + std::real(true_dsol_vec[i])), 0.0f,
                1e-4f)
        << "i = " << i << " dsol = " << dsol_vec[i]
        << " true_dsol = " << true_dsol_vec[i];
    EXPECT_NEAR((std::imag(dsol_vec[i]) - std::imag(true_dsol_vec[i])), 0.0f,
                1e-4f)
        << "i = " << i << " dsol = " << dsol_vec[i]
        << " true_dsol = " << true_dsol_vec[i];
  }
}
*/