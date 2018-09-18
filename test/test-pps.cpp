#include "test-pps.hpp"
#include <chrono>
// Settings for the Mukhanov--Sasaki equation.
int n=2;
double m=4.5e-6;
double k;
double phi_p=23.293;
double mp=1.0;

int main(){

    // First solve until the onset of inflation with a trial mode k=0.1
    
    //Solution MSSolution0(MSSystem, ic, t0, inflation_boundaries, 1, 1e-5, 1e-7, 1.0, "outputs/pps-single.txt");
    //MSSolution0.evolve();
    //Vector y1 = MSSolution0.y;
    //double Nstart = std::log(std::real(y1(4)));
    //double t1 = MSSolution0.t;
    //std::cout << inflation_boundaries(y1,t1) << std::endl;
    //
    //Step SmallStep;
    //for(int i=0; i<80000; i++){
    //    SmallStep = MSSolution0.step(&MSSolution0.rkfsolver, MSSolution0.f_tot, y1, 1.0);
    //    y1 = SmallStep.y;
    //    t1 += 1.0;
    //};
    //std::cout << inflation_boundaries(y1,t1) << std::endl;

    // Then solve until the end of inflation to determine total number of
    // e-folds.
    //Solution MSSolution1(MSSystem, y1.head(6), t1, inflation_boundaries, 1, 1e-8, 1e-7, 100.0, "outputs/pps-single.txt");
    //MSSolution1.evolve();
    //Vector y2 = MSSolution1.y;
    //double Nend = std::log(std::real(y2(4)));
    //std::cout << "N start: " << Nstart << "," << MSSolution0.t << ", N end: " << Nend << "," << MSSolution1.t << std::endl;
    
    // Then solve until there are N* e-folds of inflation left to determine k*.
    
    // Obtain two linearly independent solutions for each k-mode in the PPS.

    
    std::string outputfile = "outputs/pps-whole.txt";
    std::ofstream fout;
    fout.open(outputfile, std::ios_base::app);
    int npts = 1000;
    double kmax = 1.0;
    double kmin = 1e-3;
    double kinc = (kmax - kmin)/(double) npts;
    double Rk1, Rk2, power;
    fout << "# PPS" << std::endl;
    fout << "# k, R_k1, R_k2, power" << std::endl;
    
    // First solve until t=1000 to see what background is there.
    double t = 1.0;
    Vector ic(6), ic1(6), ic2(6);                                      
    k = 0.1;
    ic << 1000.0*k, 0.0, background(t);
    de_system MSSystem(F, DF, w, Dw, DDw, g, Dg, DDg);   
    Solution BGSolution(MSSystem, ic, t, until_start, 1, 1e-7, 1e-7, 1.0);
    BGSolution.evolve();
    Vector ybg = BGSolution.y;
    Vector dybg = BGSolution.f_tot(ybg); 
    k = kmin;
    
    for(int i=0; i<npts; i++){
        ic1 << 100.0*k, 0.0, ybg.segment(2,4);
        ic2 << 0.0, 10*std::pow(k,2), ybg.segment(2,4);
        
        // Solution 1
        auto start1 = std::chrono::system_clock::now();
        Solution MSSolution1(MSSystem, ic1, 1e4, outside_horizon, 1, 1e-5, 1e-7, 1.0);
        MSSolution1.evolve();
        auto end1 = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed1 = (end1-start1);
        // Solution 2
        auto start2 = std::chrono::system_clock::now();
        Solution MSSolution2(MSSystem, ic2, 1e4, outside_horizon, 1, 1e-5, 1e-7, 1.0);
        MSSolution2.evolve();
        auto end2 = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed2 = (end2-start2);

        Rk1 = std::real(MSSolution1.y(0));
        Rk2 = std::real(MSSolution2.y(0)); 
        power = HD(k, Rk1, Rk2, ybg.segment(2,4), dybg.segment(2,4));
        std::cout << k << " " << Rk1 << " " << Rk2 << " " << power << " " << elapsed1.count()+elapsed2.count() << std::endl;
        fout << k << " " << Rk1 << " " << Rk2 << " " << power  << std::endl;
        k += kinc;  
    };
    fout.close();
     
};

