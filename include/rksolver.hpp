#pragma once
#include "system.hpp"

/** Class that defines a 4 and 5th order Runge-Kutta method.
 *
 * This is a "bespoke"
 * Runge-Kutta formula based on the nodes used in 5th and 6th order
 * Gauss-Lobatto integration, as detailed in [1].
 *
 * [1] Agocs, F. J., et al. “Efficient Method for Solving Highly Oscillatory Ordinary Differential Equations with Applications to Physical Systems.” Physical Review Research, vol. 2, no. 1, 2020, doi:10.1103/physrevresearch.2.013030. 
 */
class RKSolver
{
    private: 
    // Frequency and friction term
    //std::function<std::complex<double>(double)> w;
    std::function<std::complex<double>(double)> g;
    

    // Butcher tablaus
    Eigen::Matrix<double,5,5> butcher_a5;
    Eigen::Matrix<double,3,3> butcher_a4;
    Eigen::Matrix<double,6,1> butcher_b5, butcher_c5, dense_b5;
    Eigen::Matrix<double,4,1> butcher_b4, butcher_c4;

    // Current values of w, g
    std::complex<double> wi, gi;
   
    public:

    /** Defines the ODE */
    de_system *de_sys_;

    /** Callable that gives the frequency term in the ODE at a given time */
    std::function<std::complex<double>(double)> w;
    
    // constructors
    RKSolver();
    RKSolver(de_system &de_sys);
    // grid of ws, gs
    /** 6 values of the frequency term per step, evaluated at the nodes of 6th
     * order Gauss-Lobatto quadrature
     */
    Eigen::Matrix<std::complex<double>,6,1> ws;
    /** 6 values of the friction term per step, evaluated at the nodes of 6th
     * order Gauss-Lobatto quadrature
     */
    Eigen::Matrix<std::complex<double>,6,1> gs;
    /** 5 values of the frequency term per step, evaluated at the nodes of 5th
     * order Gauss-Lobatto quadrature
     */
    Eigen::Matrix<std::complex<double>,5,1> ws5;
    /** 5 values of the friction term per step, evaluated at the nodes of 5th
     * order Gauss-Lobatto quadrature
     */
    Eigen::Matrix<std::complex<double>,5,1> gs5;

    // single step function 
    Eigen::Matrix<std::complex<double>,2,2> step(std::complex<double>, std::complex<double>, double, double);
    // ODE to solve
    Eigen::Matrix<std::complex<double>,1,2> f(double t, const Eigen::Matrix<std::complex<double>,1,2> &y);  
    // For dense output
    Eigen::Matrix<std::complex<double>,6,2> k5;
    Eigen::Matrix<std::complex<double>,1,2> dense_point(std::complex<double> x, std::complex<double> dx, const Eigen::Matrix<std::complex<double>,6,2> &k5);
    Eigen::Matrix<std::complex<double>,7,2> k_dense;
    Eigen::Matrix<double,7,4> P_dense;
    void dense_step(double t0, double h0, std::complex<double> y0, std::complex<double> dy0, const std::list<double> &dots, std::list<std::complex<double>> &doxs, std::list<std::complex<double>> &dodxs);

};

/** Default constructor. */
RKSolver::RKSolver(){
};

/** Constructor for the RKSolver class. It sets up the Butcher tableaus for the
 * two Runge-Kutta methods (4th and 5th order) used.
 *
 * @param de_sys[in] the system of first-order equations defining the
 * second-order ODE to solve. 
 */
RKSolver::RKSolver(de_system &de_sys){

   
    de_sys_ = &de_sys;
    // Set Butcher tableaus
    RKSolver::butcher_a5 << 0.1174723380352676535740,0,0,0,0,
                 -0.1862479800651504276304,0.5436322218248278794734,0,0,0,
                 -0.6064303885508280518989,1,0.2490461467911506000559,0,0,
                 2.899356540015731406420,-4.368525611566240669139,2.133806714786316899991,0.2178900187289247091542,0,
                 18.67996349995727204273,-28.85057783973131956546,10.72053408420926869789,1.414741756508049078612,-0.9646615009432702537787;
	RKSolver::butcher_a4 << 0.172673164646011428100,0,0,
                -1.568317088384971429762,2.395643923738960001662,0,
                -8.769507466172720011410,10.97821961869480000808,-1.208712152522079996671;
	RKSolver::butcher_b5 << 0.1127557227351729739820,0,0.5065579732655351768382,0.04830040376995117617928,0.3784749562978469803166,-0.04608905606850630731611;
	RKSolver::dense_b5 << 0.2089555395,0.,0.7699501023,0.009438629906,-0.003746982422,0.01540271068; 
	RKSolver::butcher_c5 << 0,0.117472338035267653574,0.357384241759677451843,0.642615758240322548157,0.882527661964732346426,1;
	RKSolver::butcher_b4 << -0.08333333333333333333558,0.5833333333333333333357,0.5833333333333333333356,-0.08333333333333333333558;
	RKSolver::butcher_c4 << 0,0.172673164646011428100,0.827326835353988571900,1;
    RKSolver::P_dense << 1.        , -2.48711376,  2.42525041, -0.82538093,
                         0.        ,  0.        ,  0.        ,  0.,
                         0.        ,  3.78546138, -5.54469086,  2.26578746,
                         0.        , -0.27734213,  0.74788587, -0.42224334,
                         0.        , -2.94848704,  7.41087391, -4.08391191,
                         0.        ,  0.50817346, -1.20070313,  0.64644062,
                         0.        ,  1.4193081 , -3.8386162 ,  2.4193081;

};

/** Turns the second-order ODE into a system of first-order ODEs as follows:
 * 
 * \f[ y = [x, \dot{x}], \f]
 * \f[ \dot{y[0]} = y[1], \f]
 * \f[ \dot{y[1]} = -\omega^2(t)y[0] -2\gamma(t)y[1]. \f]
 *
 * @param t[in] time \f$ t \f$
 * @param y[in] vector of unknowns \f$ y = [x, \dot{x}] \f$
 * @returns a vector of the derivative of \f$ y \f$
 */
Eigen::Matrix<std::complex<double>,1,2> RKSolver::f(double t, const Eigen::Matrix<std::complex<double>,1,2> &y){
    
    if(de_sys_->islogw_)
        wi = de_sys_->Winterp.expit(t);
    else
        wi = de_sys_->Winterp(t);
    if(de_sys_->islogg_)
        gi = de_sys_->Ginterp.expit(t);
    else
        gi = de_sys_->Ginterp(t);
    Eigen::Matrix<std::complex<double>,1,2> result;
    result << y[1], -wi*wi*y[0]-2.0*gi*y[1];
    return result;
};

/** Gives dense output at a single point during the step "for free", i.e. at no
 * extra evaluations of \f$ \omega(t), \gamma(t) \f$. This solution is roughly
 * mid-way through the step at $\sigma \sim 0.59 \f$, where \f$ \sigma = 0 \f$
 * corresponds to the start, \f$ \sigma = 1 \f$ to the end of the step.
 */
Eigen::Matrix<std::complex<double>,1,2> RKSolver::dense_point(std::complex<double> x, std::complex<double> dx, const Eigen::Matrix<std::complex<double>,6,2> &k5){

    Eigen::Matrix<std::complex<double>,1,2> ydense;
    ydense << x, dx;
    for(int j=0; j<=5; j++)
        ydense += 0.5866586817*dense_b5(j)*k5.row(j);       
   return ydense; 
    
};

/** Calculated dense output at a given set of points during a step after a
 * successful Runge-Kutta type step. 
 */
void RKSolver::dense_step(double t0, double h0, std::complex<double> y0, std::complex<double> dy0, const std::list<double> &dots, std::list<std::complex<double>> &doxs, std::list<std::complex<double>> &dodxs){
    
    int docount = dots.size();
    double h = h0;
    double sig, sig2, sig3, sig4;
    Eigen::Matrix<std::complex<double>,2,4> Q_dense = k_dense.transpose()*P_dense;
    Eigen::MatrixXd  R_dense(4,docount);
    Eigen::MatrixXcd Y_dense(2,docount);
    
    int colcount = 0;
    for(auto it=dots.begin(); it!=dots.end(); it++){
        // Transform intermediate points to be in (-1,1):
        sig = (*it - t0)/h; 
        sig2 = sig*sig;
        sig3 = sig2*sig;
        sig4 = sig3*sig;
        R_dense.col(colcount) << sig, sig2, sig3, sig4;
        colcount += 1;
    }
    Y_dense = Q_dense*R_dense;
    auto it = doxs.begin();
    auto dit = dodxs.begin();
    for(int j=0; j<docount; j++){
        *it = y0 + Y_dense.row(0)(j);
        *dit = dy0 + Y_dense.row(1)(j);
        it++;
        dit++;
    }

};

/** Computes a single Runge-Kutta type step, and returns the solution and its
 * local error estimate. 
 *
 *
 */
Eigen::Matrix<std::complex<double>,2,2> RKSolver::step(std::complex<double> x0, std::complex<double> dx0, double t0, double h){

    Eigen::Matrix<std::complex<double>,1,2> y0, y, y4, y5, delta, k5_i, k4_i;
    y4 = Eigen::Matrix<std::complex<double>,1,2>::Zero();
    y0 << x0, dx0;
    y5 = y4;
    // TODO: resizing of ws5, gs5, insertion
    Eigen::Matrix<std::complex<double>,4,2> k4;
    Eigen::Matrix<std::complex<double>,2,2> result;
    k5.row(0) = h*f(t0, y0);
    ws(0) = wi;
    gs(0) = gi;
    for(int s=1; s<=5; s++){
        y = y0;
        for(int i=0; i<=(s-1); i++)
            y += butcher_a5(s-1,i)*k5.row(i);
        k5_i = h*f(t0 + butcher_c5(s)*h, y);
        k5.row(s) = k5_i;
        ws(s) = wi;
        gs(s) = gi;
    }
    k4.row(0) = k5.row(0);
    ws5(0) = ws(0);
    gs5(0) = gs(0);
    for(int s=1; s<=3; s++){
       y = y0;
       for(int i=0; i<=(s-1); i++)
           y += butcher_a4(s-1,i)*k4.row(i);
        k4_i = h*f(t0 + butcher_c4(s)*h, y);
        k4.row(s) = k4_i;
        ws5(s) = wi;
        gs5(s) = gi;
    }
    for(int j=0; j<=5; j++)
        y5 += butcher_b5(j)*k5.row(j);
    for(int j=0; j<=3; j++)
        y4 += butcher_b4(j)*k4.row(j);
    delta = y5 - y4;
    y5 += y0;
    result << y5, delta;
    // Add in missing w, g at t+h/2
    ws5(4) = ws5(3);
    ws5(3) = ws5(2);
    if(de_sys_->islogw_)
        ws5(2) = de_sys_->Winterp.expit(t0+h/2);
    else
        ws5(2) = de_sys_->Winterp(t0+h/2);
    gs5(4) = gs5(3);
    gs5(3) = gs5(2);
    if(de_sys_->islogg_)
        gs5(2) = de_sys_->Ginterp.expit(t0+h/2);
    else
        gs5(2) = de_sys_->Ginterp(t0+h/2);


    // Fill up k_dense matrix for dense output
    for(int i=0; i<=5; i++)
        k_dense.row(i) = k5.row(i);
    k_dense.row(6) = h*f(t0+h,y5);
    return result;
};

