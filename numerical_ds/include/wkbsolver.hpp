#pragma once
#include "system.hpp"

class WKBSolver
{
    private: 
    // Frequency and friction terms 
    std::complex<double> (*w)(double);
    std::complex<double> (*g)(double);   
    // Derivative terms 
    void d1w1(double,double);
    void d1w2(double,double);
    void d1w3(double,double);
    void d1w4(double,double);
    void d1w5(double,double);
    void d1w6(double,double);
    void d2w1(double,double);
    void d2w6(double,double);
    void d3w1(double,double);
    void d3w6(double,double);
    void d4w1(double,double);
    void d1g1(double,double);
    void d1g6(double,double);
    void d2g1(double,double);
    void d2g6(double,double);
    void d3g1(double,double);
    void d1w2_5(double,double);
    void d1w3_5(double,double);
    void d1w4_5(double,double);
    // WKB series and derivatives (order dependent)
    virtual void dds();
    virtual void dsi();
    virtual void dsf();
    virtual void s();
    // WKB solutions and derivatives
    void fp();
    void fm();
    void dfp();
    void dfm();
    void ddfp();
    void ddfm();
    void ap();
    void am();
    void bp();
    void bm();

    // Gauss-Lobatto n=6, 5 weights
    Eigen::Matrix<double,6,1> glws6;
    Eigen::Matrix<double,5,1> glws5;
    // weights for derivatives
    Eigen::Matrix<double,7,1> d4w1_w;
    Eigen::Matrix<double,6,1> d1w1_w, d1w2_w, d1w3_w, d1w4_w, d1w5_w, d1w6_w, d2w1_w, d2w6_w, d3w1_w, d3w6_w, d1g1_w, d1g6_w, d2g1_w, d2g6_w, d3g1_w;
    Eigen::Matrix<double,5,1> d1w2_5_w, d1w3_5_w, d1w4_5_w;
    // grid of ws, gs
    Eigen::Matrix<std::complex<double>,6,1> ws_, gs_;
    Eigen::Matrix<std::complex<double>,5,1> ws5_, gs5_;
    // derivatives
    std::complex<double> d1w1_, d1w2_, d1w3_, d1w4_, d1w5_, d1w6_, d2w1_, d2w6_, d3w1_, d3w6_, d4w1_, d1g1_, d1g6_, d2g1_, d2g6_, d3g1_; 
    std::complex<double> d1w2_5_, d1w3_5_, d1w4_5_;
    // WKB series and their derivatives
    Eigen::Matrix<std::complex<double>,1,4> dds_, dsi_, dsf_, s_; 
    // WKB solutions and their derivatives
    std::complex<double> fp_, fm_, dfp_, dfm_, ddfp_, ddfm_, ap_, am_, bp_, bm_; 

    public:
    // constructor
    WKBSolver();
    WKBSolver(de_system);
    Eigen::Matrix<std::complex<double>,2,2> step(std::complex<double> x0, std::complex<double> dx0, double t0, double h, const Eigen::Matrix<std::complex<double>,6,1> &ws, const Eigen::Matrix<std::complex<double>,6,1> &gs, const Eigen::Matrix<std::complex<double>,5,1> &ws5, const Eigen::Matrix<std::complex<double>,5,1> &gs5); 

};

WKBSolver::WKBSolver(){
};

WKBSolver::WKBSolver(de_system de_sys){
     
    // Set frequency and friction terms
    w = de_sys.w;
    g = de_sys.g;
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
    d4w1_w << 9744.00062637928, -27851.6858893579, 75653.4044616243,
        -107520.008443354, 61113.3089030151, -15843.0200709916, 4704.00041268436;
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

Eigen::Matrix<std::complex<double>,2,2> WKBSolver::step(std::complex<double> x0, std::complex<double> dx0, double t0, double h, const Eigen::Matrix<std::complex<double>,6,1> &ws, const Eigen::Matrix<std::complex<double>,6,1> &gs, const Eigen::Matrix<std::complex<double>,5,1> &ws5, const Eigen::Matrix<std::complex<double>,5,1> &gs5){
    
    // Set grid of ws, gs:
    ws_ = ws;
    gs_ = gs;
    ws5_ = ws5;
    gs5_ = gs5;
    Eigen::Matrix<std::complex<double>,2,2> result;
    std::cout << "ws : " << ws_ << ", gs: " << gs_ << std::endl;
    result << 0.0, 0.0, 0.0, 0.0;

    return result;
};

void WKBSolver::dds(){
};

void WKBSolver::dsi(){
};

void WKBSolver::dsf(){
};

void WKBSolver::s(){
};

void WKBSolver::d1w1(double t, double h){
    d1w1_ = d1w1_w.dot(ws_)/h;
};

void  WKBSolver::d1w2(double t, double h){
    d1w2_ = d1w2_w.dot(ws_)/h;
};

void WKBSolver::d1w3(double t, double h){
    d1w3_ = d1w3_w.dot(ws_)/h;
};

void WKBSolver::d1w4(double t, double h){
    d1w4_ = d1w4_w.dot(ws_)/h;
};

void WKBSolver::d1w5(double t, double h){
    d1w5_ = d1w5_w.dot(ws_)/h;
};

void WKBSolver::d1w6(double t, double h){
    d1w6_ = d1w6_w.dot(ws_)/h;
};

void WKBSolver::d2w1(double t, double h){
    d2w1_ = d2w1_w.dot(ws_)/(h*h);
};

void WKBSolver::d2w6(double t, double h){
    d2w6_ = d2w6_w.dot(ws_)/(h*h);
};

void WKBSolver::d3w1(double t, double h){
    d3w1_ = d3w1_w.dot(ws_)/(h*h*h);
};

void WKBSolver::d3w6(double t, double h){
    d3w6_ = d3w6_w.dot(ws_)/(h*h*h);
};

//void WKBSolver::d4w1(double t, double h){
    // TODO - increase size of ws
    //d4w1_ =  d4w1_w.dot(ws)/(h*h*h*h);
//};

void WKBSolver::d1g1(double t, double h){
    d1g1_ = d1g1_w.dot(gs_)/h;
};

void WKBSolver::d1g6(double t, double h){
    d1g6_ = d1g6_w.dot(gs_)/h;
};

void WKBSolver::d2g1(double t, double h){
    d2g1_ = d2g1_w.dot(gs_)/(h*h);
};

void WKBSolver::d2g6(double t, double h){
    d2g6_ = d2g6_w.dot(gs_)/(h*h);
};

void WKBSolver::d3g1(double t, double h){
    d3g1_ = d3g1_w.dot(gs_)/(h*h*h);
};

void WKBSolver::d1w2_5(double t, double h){
    d1w2_5_ = d1w2_5_w.dot(ws5_)/h;
};

void WKBSolver::d1w3_5(double t, double h){
    d1w3_5_ = d1w3_5_w.dot(ws5_)/h;
};

void WKBSolver::d1w4_5(double t, double h){
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

};

WKBSolver1::WKBSolver1(){
};

void WKBSolver1::dds(){
};

void WKBSolver1::dsi(){
};

void WKBSolver1::dsf(){
};

void WKBSolver1::s(){
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

};

WKBSolver2::WKBSolver2(){
};

void WKBSolver2::dds(){
};

void WKBSolver2::dsi(){
};

void WKBSolver2::dsf(){
};

void WKBSolver2::s(){
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

};

WKBSolver3::WKBSolver3(){
};

void WKBSolver3::dds(){
};

void WKBSolver3::dsi(){
};

void WKBSolver3::dsf(){
};

void WKBSolver3::s(){
};

