#include "test-pps.hpp"
#include <chrono>
// Settings for the Mukhanov--Sasaki equation.
int n=2;
double m=4.5e-6;
double k;
double phi_p=23.293;
double mp=1.0;
double tstart=1e4;

int main(){

    std::string outputfile = "outputs/pps-hightol.txt";
    std::ofstream fout;
    int npts = 100;
    double kmin = 4.0e4;
    double kmax = 1e+5;
    double kinc = std::exp((std::log(kmax) - std::log(kmin))/(double) npts);
    double Rk1, Rk2, power1, power2;
    double t = 1.0;
    double rtol=1e-4;
    double atol=1e-7;
    int order=1;
    
    fout.open(outputfile, std::ios_base::app);
    fout << "# PPS" << std::endl;
    fout << "# parameters: n=" << n << ", m=" << m << ", phi_p=" << phi_p << ", mp=" << mp << ", ic set at tstart=" << tstart << ", rtol=" << rtol << ", atol=" << atol << ", order=" << order << std::endl;
    fout << "# k, R_k1, R_k2, P_R(HD), P_R(RST), dt" << std::endl;
    
    // First solve until t=tstart to see what background is there.
    Vector ic(6), ic1(6), ic2(6);                                      
    k = 0.1;
    ic << 1000.0*k, 0.0, background(t);
    de_system MSSystem(F, DF, w, Dw, DDw, g, Dg, DDg);   
    Solution BGSolution(MSSystem, ic, t, until_start, order, 1e-7, atol, 1.0);
    BGSolution.evolve();
    Vector ybg = BGSolution.y;
    Vector dybg = BGSolution.f_tot(ybg); 
    k = kmin;
    
    // Then solve the Mukhanov--Sasaki equation for each k-mode, obtaining two
    // linearly indepdendent solutions.
    for(int i=0; i<npts; i++){
        ic1 << 100.0*k, 0.0, ybg.segment(2,4);
        ic2 << 0.0, 10*std::pow(k,2), ybg.segment(2,4);
        
        auto start1 = std::chrono::system_clock::now();
        Solution MSSolution1(MSSystem, ic1, tstart, outside_horizon, order, rtol, atol, 1.0);
        MSSolution1.evolve();
        auto end1 = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed1 = (end1-start1);
        
        auto start2 = std::chrono::system_clock::now();
        Solution MSSolution2(MSSystem, ic2, tstart, outside_horizon, order, rtol, atol, 1.0);
        MSSolution2.evolve();
        auto end2 = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed2 = (end2-start2);

        Rk1 = std::real(MSSolution1.y(0));
        Rk2 = std::real(MSSolution2.y(0)); 
        power1 = HD(k, Rk1, Rk2, ybg.segment(2,4), dybg.segment(2,4));
        power2 = RST(k, Rk1, Rk2, ybg.segment(2,4), dybg.segment(2,4));
        std::cout << k << " " << Rk1 << " " << Rk2 << " " << power1 << " " << power2 << " " << elapsed1.count()+elapsed2.count() << std::endl;
        fout << k << " " << Rk1 << " " << Rk2 << " " << power1 << " " << power2  << " " << elapsed1.count()+elapsed2.count() << std::endl;
        k *= kinc;  
    };
    fout.close();
     
};

