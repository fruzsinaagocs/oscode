#pragma once
#include <iomanip>
#include "system.hpp"

class WKBSolver
{
    protected: 
    // Derivative terms 
    void d1w1();
    void d1w2();
    void d1w3();
    void d1w4();
    void d1w5();
    void d1w6();
    void d2w1();
    void d2w6();
    void d3w1();
    void d3w6();
    void d4w1();
    void d1g1();
    void d1g6();
    void d2g1();
    void d2g6();
    void d3g1();
    void d1w2_5();
    void d1w3_5();
    void d1w4_5();
    // WKB series and derivatives (order dependent)
    virtual void dds();
    virtual void dsi();
    virtual void dsf();
    virtual void s();
    // WKB solutions and derivatives
    void fp();
    void fm();
    void dfpi();
    void dfmi();
    void dfpf();
    void dfmf();
    void ddfp();
    void ddfm();
    void ap();
    void am();
    void bp();
    void bm();
    // Gauss-Lobatto integration
    Eigen::Matrix<std::complex<double>,2,1> integrate(const
    Eigen::Matrix<std::complex<double>,6,1> &integrand6, const
    Eigen::Matrix<std::complex<double>,5,1> &integrand5);

    // Gauss-Lobatto n=6, 5 weights
    Eigen::Matrix<double,6,1> glws6;
    Eigen::Matrix<double,5,1> glws5;
    // weights for derivatives
    Eigen::Matrix<double,7,1> d4w1_w;
    Eigen::Matrix<double,6,1> d1w1_w, d1w2_w, d1w3_w, d1w4_w, d1w5_w, d1w6_w,
    d2w1_w, d2w6_w, d3w1_w, d3w6_w, d1g1_w, d1g6_w, d2g1_w, d2g6_w, d3g1_w;
    Eigen::Matrix<double,5,1> d1w2_5_w, d1w3_5_w, d1w4_5_w;
    // grid of ws, gs
    Eigen::Matrix<std::complex<double>,7,1> ws7_;
    Eigen::Matrix<std::complex<double>,6,1> ws_, gs_;
    Eigen::Matrix<std::complex<double>,5,1> ws5_, gs5_;
    // derivatives
    std::complex<double> d1w1_, d1w2_, d1w3_, d1w4_, d1w5_, d1w6_, d2w1_, d2w6_,
    d3w1_, d3w6_, d4w1_, d1g1_, d1g6_, d2g1_, d2g6_, d3g1_; 
    std::complex<double> d1w2_5_, d1w3_5_, d1w4_5_;
    Eigen::Matrix<std::complex<double>,6,1> dws_;
    Eigen::Matrix<std::complex<double>,5,1> dws5_;
    // WKB series and their derivatives
    Eigen::Matrix<std::complex<double>,1,4> dds_, dsi_, dsf_, s_; 
    // Error in WKB series
    Eigen::Matrix<std::complex<double>,1,4> s_error;
    // WKB solutions and their derivatives
    std::complex<double> fp_, fm_, dfpi_, dfmi_, dfpf_, dfmf_, ddfp_, ddfm_, ap_, am_, bp_, bm_; 
    // step and stepsize
    double h;
    std::complex<double> x, dx, ddx;
    // order
    int order_;
    // error estimate on step
    std::complex<double> err_fp, err_fm, err_dfp, err_dfm;

    public:
    // constructor
    WKBSolver();
    WKBSolver(de_system &de_sys, int order);
    Eigen::Matrix<std::complex<double>,3,2> step(std::complex<double> x0,
    std::complex<double> dx0, double t0, double h0, const
    Eigen::Matrix<std::complex<double>,6,1> &ws, const
    Eigen::Matrix<std::complex<double>,6,1> &gs, const
    Eigen::Matrix<std::complex<double>,5,1> &ws5, const
    Eigen::Matrix<std::complex<double>,5,1> &gs5); 

};

WKBSolver::WKBSolver(){
};

WKBSolver::WKBSolver(de_system &de_sys, int order){

    // Set order
    order_ = order;
    // Set Gauss-Lobatto weights
    glws6 << 1.0/15.0,0.3784749562978469803166,0.5548583770354863530167,
        0.5548583770354863530167,0.3784749562978469803166,1.0/15.0;
    glws5 << 1.0/10.0,49.0/90.0,32.0/45.0,49.0/90.0,1.0/10.0;
    // Set finite difference weights
    d1w1_w << -15.0000000048537, 20.2828318761850,
        -8.07237453994912, 4.48936929577350, -2.69982662677053, 0.999999999614819;
    d1w2_w << -3.57272991033049, 0.298532922755350e-7  ,
        5.04685352597795, -2.30565629452303, 1.30709499910514, -0.475562350082855;
    d1w3_w << 0.969902096162109, -3.44251390568294,
    -0.781532641131861e-10, 3.50592393061265, -1.57271334190619, 0.539401220892526 ;
    d1w4_w << -0.539401220892533, 1.57271334190621,
        -3.50592393061268, 0.782075077478921e-10, 3.44251390568290, -0.969902096162095;
    d1w5_w << 0.475562350082834, -1.30709499910509,
        2.30565629452296, -5.04685352597787, -0.298533681980831e-7, 3.57272991033053;
    d1w6_w << -0.999999999614890, 2.69982662677075,
        -4.48936929577383, 8.07237453994954, -20.2828318761854, 15.0000000048538;
    d2w1_w << 140.000000016641, -263.163968874741,
        196.996471291466, -120.708905753218, 74.8764032980854, -27.9999999782328;
    d2w6_w << -27.9999999782335, 74.8764032980873,
        -120.708905753221, 196.996471291469, -263.163968874744, 140.000000016642;
    d3w1_w << -840.000000234078, 1798.12714381468,
        -1736.74461287884, 1322.01528240287, -879.397812956524, 335.999999851893;
    d3w6_w << -335.999999851897, 879.397812956534,
        -1322.01528240289, 1736.74461287886, -1798.12714381470, 840.000000234086;
    d4w1_w << //9744.00062637928, -27851.6858893579, 75653.4044616243,
        //-107520.008443354, 61113.3089030151, -15843.0200709916, 4704.00041268436;
            3024.00000383582, -6923.06197480357, 7684.77676018742, 0.0, -6855.31809730784, 5085.60330881706, -2016.00000072890;
    d1g1_w << -15.0000000048537, 20.2828318761850,
    -8.07237453994912, 4.48936929577350, -2.69982662677053, 0.999999999614819 ;
    d1g6_w << -0.999999999614890, 2.69982662677075,
        -4.48936929577383, 8.07237453994954, -20.2828318761854, 15.0000000048538;
    d2g1_w << 140.000000016641, -263.163968874741,
        196.996471291466, -120.708905753218, 74.8764032980854, -27.9999999782328;
    d2g6_w << -27.9999999782335, 74.8764032980873,
        -120.708905753221, 196.996471291469, -263.163968874744, 140.000000016642;
    d3g1_w << -840.000000234078, 1798.12714381468,
        -1736.74461287884, 1322.01528240287, -879.397812956524, 335.999999851893;
    d1w2_5_w << -2.48198050935042, 0.560400997591235e-8,
        3.49148624058567, -1.52752523062733, 0.518019493788063;
    d1w3_5_w << 0.750000000213852, -2.67316915534181,
        0.360673032443906e-10, 2.67316915534181, -0.750000000213853;
    d1w4_5_w << -0.518019493788065, 1.52752523062733,
        -3.49148624058568, -0.560400043118500e-8, 2.48198050935041;
};

Eigen::Matrix<std::complex<double>,3,2> WKBSolver::step(std::complex<double> x0,
std::complex<double> dx0, double t0, double h0, const
Eigen::Matrix<std::complex<double>,6,1> &ws, const
Eigen::Matrix<std::complex<double>,6,1> &gs, const
Eigen::Matrix<std::complex<double>,5,1> &ws5, const
Eigen::Matrix<std::complex<double>,5,1> &gs5){
    
    Eigen::Matrix<std::complex<double>,3,2> result;
    result << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    // Set grid of ws, gs:
    ws_ = ws;
    gs_ = gs;
    ws5_ = ws5;
    gs5_ = gs5;
    ws7_ << ws_(0), ws_(1), ws_(2), ws5_(2), ws_(3), ws_(4), ws_(5);
    // Set i.c.
    x = x0;
    dx = dx0;
    ddx = -std::pow(ws_(0),2)*x - 2.0*gs_(0)*dx;
    // Set derivatives:
    h = h0;
    d1w1(); d1w2(); d1w3(); d1w4(); d1w5(); d1w6(); d2w1(); d2w6(); d3w1();
    d3w6(); d4w1(); d1g1(); d1g6(); d2g1(); d2g6(); d3g1(); d1w2_5(); d1w3_5(); d1w4_5();
    dws_ << d1w1_, d1w2_, d1w3_, d1w4_, d1w5_, d1w6_;
    dws5_ << d1w1_, d1w2_5_, d1w3_5_, d1w4_5_, d1w6_;
    // Higher order step
        // Calculate A, B
        fm_ = 1.0;
        fp_ = 1.0;
        s_ << 0.0, 0.0, 0.0, 0.0; dsi(); dds();
        dfpi(); dfmi();
        ddfp(); ddfm();
        ap(); am(); bp(); bm();
        // Calculate step
        s(); 
        dsf();
        fp(); fm(); 
        dfpf(); dfmf();
        result(0,0) = ap_*fp_ + am_*fm_;
        result(0,1) = bp_*dfpf_ +bm_*dfmf_;
    // Error estimate on this
        err_fp = s_error.cwiseAbs().sum()*fp_;
        err_fm = std::conj(err_fp);
        err_dfp = dfpf_/fp_*err_fp; 
        err_dfm = std::conj(err_dfp);
        result(2,0) = ap_*err_fp + am_*err_fm;
        result(2,1) = bp_*err_dfp + bm_*err_dfm;
    // Lower order step for correction
        // A, B
        dsi_(order_) = 0.0; dds_(order_) = 0.0;
        dfpi(); dfmi();
        ddfp(); ddfm();
        ap(); am(); bp(); bm();
        // Calculate step
        s_(order_) = 0.0;
        dsf_(order_) = 0.0;
        fp(); fm();
        dfpf(); dfmf();
        result(1,0) = result(0,0) - ap_*fp_ - am_*fm_;
        result(1,1) = result(0,1) - bp_*dfpf_ - bm_*dfmf_;

    return result;
};

Eigen::Matrix<std::complex<double>,2,1> WKBSolver::integrate(const
Eigen::Matrix<std::complex<double>,6,1> &integrand6, const
Eigen::Matrix<std::complex<double>,5,1> &integrand5){
    
    std::complex<double> x6 = h/2.0*glws6.dot(integrand6);
    std::complex<double> x5 = h/2.0*glws5.dot(integrand5);
    Eigen::Matrix<std::complex<double>,2,1> result;
    result << x6, x6-x5;
    return result;
};

void WKBSolver::fp(){
    fp_ = std::exp(s_.sum());
};

void WKBSolver::fm(){
    fm_ = std::conj(fp_); 
};

void WKBSolver::dfpi(){
    dfpi_ = dsi_.sum();    
};

void WKBSolver::dfmi(){
    dfmi_ = std::conj(dfpi_);
};

void WKBSolver::dfpf(){
    dfpf_ = dsf_.sum()*fp_;    
};

void WKBSolver::dfmf(){
    dfmf_ = std::conj(dfpf_);
};

void WKBSolver::ddfp(){
    ddfp_ = dds_.sum() + std::pow(dsi_.sum(),2);
};

void WKBSolver::ddfm(){
    ddfm_ = std::conj(ddfp_);
};

void WKBSolver::ap(){
    ap_ = (dx - x*dfmi_)/(dfpi_ - dfmi_);
};

void WKBSolver::am(){
    am_ = (dx - x*dfpi_)/(dfmi_ - dfpi_);
};

void WKBSolver::bp(){
    bp_ = (ddx*dfmi_ - dx*ddfm_)/(ddfp_*dfmi_ - ddfm_*dfpi_);
};

void WKBSolver::bm(){
    bm_ = (ddx*dfpi_ - dx*ddfp_)/(ddfm_*dfpi_ - ddfp_*dfmi_);
};

void WKBSolver::dds(){
    dds_ << std::complex<double>(0.0, 1.0)*d1w1_, 0.0, 0.0, 0.0;
};

void WKBSolver::dsi(){
   dsi_ << std::complex<double>(0.0, 1.0)*ws_(0), 0.0, 0.0, 0.0;  
};

void WKBSolver::dsf(){
   dsf_ << std::complex<double>(0.0, 1.0)*ws_(5), 0.0, 0.0, 0.0;
};

void WKBSolver::s(){
    Eigen::Matrix<std::complex<double>,2,1> s0;
    s0 << std::complex<double>(0,1)*integrate(ws_, ws5_);  
    s_ << s0(0), 0.0, 0.0, 0.0;
    s_error << s0(1), 0.0, 0.0, 0.0; 
};

void WKBSolver::d1w1(){
    d1w1_ = d1w1_w.dot(ws_)/h;
};

void  WKBSolver::d1w2(){
    d1w2_ = d1w2_w.dot(ws_)/h;
};

void WKBSolver::d1w3(){
    d1w3_ = d1w3_w.dot(ws_)/h;
};

void WKBSolver::d1w4(){
    d1w4_ = d1w4_w.dot(ws_)/h;
};

void WKBSolver::d1w5(){
    d1w5_ = d1w5_w.dot(ws_)/h;
};

void WKBSolver::d1w6(){
    d1w6_ = d1w6_w.dot(ws_)/h;
};

void WKBSolver::d2w1(){
    d2w1_ = d2w1_w.dot(ws_)/(h*h);
};

void WKBSolver::d2w6(){
    d2w6_ = d2w6_w.dot(ws_)/(h*h);
};

void WKBSolver::d3w1(){
    d3w1_ = d3w1_w.dot(ws_)/(h*h*h);
};

void WKBSolver::d3w6(){
    d3w6_ = d3w6_w.dot(ws_)/(h*h*h);
};

void WKBSolver::d4w1(){
    d4w1_ =  d4w1_w.dot(ws7_)/(h*h*h*h);
};

void WKBSolver::d1g1(){
    d1g1_ = d1g1_w.dot(gs_)/h;
};

void WKBSolver::d1g6(){
    d1g6_ = d1g6_w.dot(gs_)/h;
};

void WKBSolver::d2g1(){
    d2g1_ = d2g1_w.dot(gs_)/(h*h);
};

void WKBSolver::d2g6(){
    d2g6_ = d2g6_w.dot(gs_)/(h*h);
};

void WKBSolver::d3g1(){
    d3g1_ = d3g1_w.dot(gs_)/(h*h*h);
};

void WKBSolver::d1w2_5(){
    d1w2_5_ = d1w2_5_w.dot(ws5_)/h;
};

void WKBSolver::d1w3_5(){
    d1w3_5_ = d1w3_5_w.dot(ws5_)/h;
};

void WKBSolver::d1w4_5(){
    d1w4_5_ = d1w4_5_w.dot(ws5_)/h;
};

//////////////////////////////////

class WKBSolver1 : public WKBSolver
{
    private: 
        void dds();
        void dsi();
        void dsf();
        void s();       

    public:
        WKBSolver1();
        WKBSolver1(de_system &de_sys, int order);

};

WKBSolver1::WKBSolver1(){
};

WKBSolver1::WKBSolver1(de_system &de_sys, int order) : WKBSolver(de_sys, order){
};

void WKBSolver1::dds(){
    dds_ << std::complex<double>(0,1)*d1w1_,
    1.0/std::pow(ws_(0),2)*std::pow(d1w1_,2)/2.0-1.0/ws_(0)*d2w1_/2.0-d1g1_, 0.0, 0.0;
};

void WKBSolver1::dsi(){
   dsi_ << std::complex<double>(0.0, 1.0)*ws_(0), -0.5*d1w1_/ws_(0)-gs_(0), 0.0, 0.0;
};

void WKBSolver1::dsf(){
   dsf_ << std::complex<double>(0.0, 1.0)*ws_(5), -0.5*d1w6_/ws_(5)-gs_(5), 0.0, 0.0;
};

void WKBSolver1::s(){
    Eigen::Matrix<std::complex<double>,2,1> s0, s1;
    s0 << std::complex<double>(0,1)*integrate(ws_, ws5_);  
    s1 << integrate(gs_, gs5_);
    s1(0) = std::log(std::sqrt(ws_(0)/ws_(5))) - s1(0);
    s_ << s0(0), s1(0), 0.0, 0.0;
    s_error << s0(1), s1(1), 0.0, 0.0;

};

//////////////////////////////////

class WKBSolver2 : public WKBSolver
{
    private: 
        void dds();
        void dsi();
        void dsf();
        void s();       

    public:
        WKBSolver2();
        WKBSolver2(de_system &de_sys, int order);

};

WKBSolver2::WKBSolver2(){
};

WKBSolver2::WKBSolver2(de_system &de_sys, int order) : WKBSolver(de_sys, order){
};

void WKBSolver2::dds(){
    dds_ << std::complex<double>(0,1)*d1w1_,
    1.0/std::pow(ws_(0),2)*std::pow(d1w1_,2)/2.0-1.0/ws_(0)*d2w1_/2.0-d1g1_,
    -std::complex<double>(0,1/8)*(8.0*d1g1_*gs_(0)*std::pow(ws_(0),3)-4.0*d1w1_*std::pow(gs_(0),2)*std::pow(ws_(0),2)+4.0*d2g1_*std::pow(ws_(0),3)-4.0*d1w1_*d1g1_*std::pow(ws_(0),2)+2.0*d3w1_*std::pow(ws_(0),2)-10.0*d1w1_*d2w1_*ws_(0)+9.0*std::pow(d1w1_,3))/std::pow(ws_(0),4), 0.0;
};

void WKBSolver2::dsi(){
    dsi_ << std::complex<double>(0,1)*ws_(0),-1.0/ws_(0)*d1w1_/2.0-gs_(0),
    std::complex<double>(0,1/8)*(-4.0*std::pow(gs_(0),2)*std::pow(ws_(0),2)-4.0*d1g1_*std::pow(ws_(0),2)-2.0*d2w1_
    *ws_(0)+3.0*std::pow(d1w1_,2))/std::pow(ws_(0),3), 0.0;
};

void WKBSolver2::dsf(){
    dsf_ << std::complex<double>(0,1)*ws_(5),-1.0/ws_(5)*d1w6_/2.0-gs_(5),
    std::complex<double>(0,1/8)*(-4.0*std::pow(gs_(5),2)*std::pow(ws_(5),2)-4.0*d1g6_*std::pow(ws_(5),2)-2.0*d2w6_
    *ws_(5)+3.0*std::pow(d1w6_,2))/std::pow(ws_(5),3), 0.0;
};

void WKBSolver2::s(){
    Eigen::Matrix<std::complex<double>,2,1> s0, s1, s2;  
    Eigen::Matrix<std::complex<double>,6,1> integrand6;
    Eigen::Matrix<std::complex<double>,5,1> integrand5;
    integrand6 = 4*gs_.cwiseProduct(gs_).cwiseQuotient(ws_) +
    4*dws_.cwiseProduct(gs_).cwiseQuotient(ws_.cwiseProduct(ws_)) +
    dws_.cwiseProduct(dws_).cwiseQuotient(ws_.cwiseProduct(ws_.cwiseProduct(ws_)));
    integrand5 =  4*gs5_.cwiseProduct(gs5_).cwiseQuotient(ws5_) +
    4*dws5_.cwiseProduct(gs5_).cwiseQuotient(ws5_.cwiseProduct(ws5_)) +
    dws5_.cwiseProduct(dws5_).cwiseQuotient(ws5_.cwiseProduct(ws5_.cwiseProduct(ws5_)));
    s0 << std::complex<double>(0,1)*integrate(ws_, ws5_);  
    s1 << integrate(gs_, gs5_);
    s1(0) = std::log(std::sqrt(ws_(0)/ws_(5))) - s1(0);
    s2 << integrate(integrand6, integrand5);
    s2(0) = -1/4.0*(dws_(5)/std::pow(ws_(5),2)+2.0*gs_(5)/ws_(5)-
        dws_(0)/std::pow(ws_(0),2)-2.0*gs_(0)/ws_(0))-1/8.0*s2(0);
    s_ << s0(0), s1(0), std::complex<double>(0,1)*s2(0), 0.0;
    s_error << s0(1), s1(1), std::complex<double>(0,-1/8)*s2(1), 0.0;

};

//////////////////////////////////

class WKBSolver3 : public WKBSolver
{
    private: 
        void dds();
        void dsi();
        void dsf();
        void s();       

    public:
        WKBSolver3();
        WKBSolver3(de_system &de_sys, int);

};

WKBSolver3::WKBSolver3(){
};

WKBSolver3::WKBSolver3(de_system &de_sys, int order) : WKBSolver(de_sys, order){
};

void WKBSolver3::dds(){
    dds_ << std::complex<double>(0,1)*d1w1_,
    1.0/std::pow(ws_(0),2)*std::pow(d1w1_,2)/2.0-1.0/ws_(0)*d2w1_/2.0-d1g1_,
    -std::complex<double>(0,1.0/8.0)*(8.0*d1g1_*gs_(0)*std::pow(ws_(0),3)-4.0*d1w1_*std::pow(gs_(0),2)*std::pow(ws_(0),2)+4.0*d2g1_*std::pow(ws_(0),3)-4.0*d1w1_*d1g1_*std::pow(ws_(0),2)+2.0*d3w1_*std::pow(ws_(0),2)-10.0*d1w1_*d2w1_*ws_(0)+9.0*std::pow(d1w1_,3))/std::pow(ws_(0),4),
    (d4w1_*std::pow(ws_(0),3)+2.0*d3g1_*std::pow(ws_(0),4)-9.0*d1w1_*d3w1_*std::pow(ws_(0),2)-6.0*std::pow(d2w1_,2)*std::pow(ws_(0),2) + (42.0*ws_(0)*std::pow(d1w1_,2)-4.0*std::pow(ws_(0),3)*(std::pow(gs_(0),2)+d1g1_))*d2w1_+(4.0*gs_(0)*std::pow(ws_(0),4)-8.0*std::pow(ws_(0),3)*d1w1_)*d2g1_-30.0*std::pow(d1w1_,4)+12.0*std::pow(ws_(0),2)*(std::pow(gs_(0),2)+d1g1_)*std::pow(d1w1_,2)-16.0*d1w1_*d1g1_*gs_(0)*std::pow(ws_(0),3)+4.0*std::pow(d1g1_,2)*std::pow(ws_(0),4))/std::pow(ws_(0),6)/8.0;
};

void WKBSolver3::dsi(){
    dsi_ << std::complex<double>(0,1)*ws_(0),-1.0/ws_(0)*d1w1_/2.0-gs_(0),
    std::complex<double>(0,1.0/8.0)*(-4.0*std::pow(gs_(0),2)*std::pow(ws_(0),2)-4.0*d1g1_*std::pow(ws_(0),2)-2.0*d2w1_
    *ws_(0)+3.0*std::pow(d1w1_,2))/std::pow(ws_(0),3) ,
    (d3w1_*std::pow(ws_(0),2)+2.0*d2g1_*std::pow(ws_(0),3)-6.0*d1w1_*d2w1_*ws_(0)+
    6.0*std::pow(d1w1_,3)-4.0*(std::pow(gs_(0),2)+d1g1_)*std::pow(ws_(0),2)*d1w1_+4.0*d1g1_
    *gs_(0)*std::pow(ws_(0),3))/std::pow(ws_(0),5)/8.0;
};

void WKBSolver3::dsf(){
    dsf_ << std::complex<double>(0,1)*ws_(5),-1.0/ws_(5)*d1w6_/2.0-gs_(5),
    std::complex<double>(0,1.0/8.0)*(-4.0*std::pow(gs_(5),2)*std::pow(ws_(5),2)-4.0*d1g6_*std::pow(ws_(5),2)-2.0*d2w6_
    *ws_(5)+3.0*std::pow(d1w6_,2))/std::pow(ws_(5),3) ,
    (d3w6_*std::pow(ws_(5),2)+2.0*d2g6_*std::pow(ws_(5),3)-6.0*d1w6_*d2w6_*ws_(5)+
    6.0*std::pow(d1w6_,3)-4.0*(std::pow(gs_(5),2)+d1g6_)*std::pow(ws_(5),2)*d1w6_+4.0*d1g6_
    *gs_(5)*std::pow(ws_(5),3))/std::pow(ws_(5),5)/8.0;
};

void WKBSolver3::s(){
    Eigen::Matrix<std::complex<double>,2,1> s0, s1, s2;  
    Eigen::Matrix<std::complex<double>,6,1> integrand6;
    Eigen::Matrix<std::complex<double>,5,1> integrand5;
    integrand6 = 4.0*gs_.cwiseProduct(gs_).cwiseQuotient(ws_) +
    4.0*dws_.cwiseProduct(gs_).cwiseQuotient(ws_.cwiseProduct(ws_)) +
    dws_.cwiseProduct(dws_).cwiseQuotient(ws_.cwiseProduct(ws_.cwiseProduct(ws_)));
    integrand5 =  4.0*gs5_.cwiseProduct(gs5_).cwiseQuotient(ws5_) +
    4.0*dws5_.cwiseProduct(gs5_).cwiseQuotient(ws5_.cwiseProduct(ws5_)) +
    dws5_.cwiseProduct(dws5_).cwiseQuotient(ws5_.cwiseProduct(ws5_.cwiseProduct(ws5_)));
    s0 << std::complex<double>(0,1)*integrate(ws_, ws5_);  
    s1 << integrate(gs_, gs5_);
    s1(0) = std::log(std::sqrt(ws_(0)/ws_(5))) - s1(0);
    s2 << integrate(integrand6, integrand5);
    s2(0) = -1/4.0*(dws_(5)/std::pow(ws_(5),2)+2.0*gs_(5)/ws_(5)-
        dws_(0)/std::pow(ws_(0),2)-2.0*gs_(0)/ws_(0))-1/8.0*s2(0);
    std::complex<double> s3 = (1/4.0*(std::pow(gs_(5),2)/std::pow(ws_(5),2) -
    std::pow(gs_(0),2)/std::pow(ws_(0),2)) + 1/4.0*(d1g6_/std::pow(ws_(5),2) -
    d1g1_/std::pow(ws_(0),2))-3/16.0*(std::pow(dws_(5),2)/std::pow(ws_(5),4) -
    std::pow(dws_(0),2)/std::pow(ws_(0),4)) + 1/8.0*(d2w6_/std::pow(ws_(5),3) -
    d2w1_/std::pow(ws_(0),3)));
    s_ << s0(0), s1(0), std::complex<double>(0,1)*s2(0), s3;
    s_error << s0(1), s1(1), std::complex<double>(0,-1.0/8.0)*s2(1), 0.0;
};

