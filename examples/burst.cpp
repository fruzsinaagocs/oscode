#include "solver.hpp"
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
int main(){

    std::ofstream f;
    /** File to write solution of ODE to */
    std::string output = "output.txt"; 
    std::complex<double> x0, dx0;
    double ti, tf;
   
    /** Define integration range */
    ti = -2*n;
    tf = 2*n;

    /** Define initial conditions */
    x0 = xburst(ti); 
    dx0 = dxburst(ti); 

    /** Create differential equation "system" */
    /** Method 1: Give frequency and damping term as functions */
    de_system sys(&w, &g);

    /** Method 2: Give frequency and damping term as std::vectors */
    /** 
    int N = 20000;
    std::vector<double> times_arr(N);
    std::vector<std::complex<double>> ws_arr(N), gs_arr(N);
    
    for(int i=0; i<N; i++){
        times_arr[i] = i/500.0;
        ws_arr[i] = w(times_arr[i]);
        gs_arr[i] = 0.0;
    }
    double * times;
    std::complex<double> * ws, * gs;
    times = times_arr.data();
    ws = ws_arr.data();
    gs = gs_arr.data();
    de_system sys(times, ws, gs, times, N);
    */ 

    /** Method 3: Give frequency and damping terms as Eigen::Vectors */
    /**
    int N = 20000;
    Eigen::VectorXd times_arr = Eigen::VectorXd::Zero(N);
    Eigen::VectorXcd ws_arr = Eigen::VectorXcd::Zero(N), gs_arr = Eigen::VectorXcd::Zero(N);

    for(int i=0; i<N; i++){
        times_arr[i] = i/500.0;
        ws_arr[i] = w(times_arr[i]);
        gs_arr[i] = 0.0;
    }
    double * times;
    std::complex<double> * ws, * gs;
    times = times_arr.data();
    ws = ws_arr.data();
    gs = gs_arr.data();
    de_system sys(times, ws, gs, times, N);
    */

    /** Method 4: Define frequency and damping terms as std::arrays */
    /**
    int N = 20000;
    double times_arr[20000];
    std::complex<double> ws_arr[20000], gs_arr[20000];

    for(int j=0; j<N; j++){
        times_arr[j] = j/500.0;
        ws_arr[j] = w(times_arr[j]);
        gs_arr[j] = 0.0;
    }
    double * times;
    std::complex<double> * ws, * gs;
    times = times_arr;
    ws = ws_arr;
    gs = gs_arr;
    de_system sys(times, ws, gs, times, N);
    */

    /** Solve the ODE */    
    Solution solution(sys, x0, dx0, ti, tf);
    solution.solve();

    /** Extract the solution and the types of steps taken by oscode */
    std::list<std::complex<double>> xs = solution.sol;
    std::list<double> ts = solution.times;
    std::list<bool> types = solution.wkbs;
    int steps = solution.ssteps;

    /** Write result in file */
    f.open(output);
    auto it_t = ts.begin();
    auto it_x = xs.begin();
    auto it_ty = types.begin();
    for(int i=0; i<=steps; i++){
        f << *it_t << ", " << std::real(*it_x) << ", " << std::imag(*it_x) << ", " <<  *it_ty << std::endl;
        ++it_t;
        ++it_x;
        ++it_ty;
    };
    f.close();
    std::cout << "Wrote results to " << output << std::endl;
};


