#pragma once
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <boost/math/special_functions/airy.hpp>
#include "system.hpp"
#include "rksolver.hpp"
#include "wkbsolver.hpp"

class Solution
{
    private:
    // Parameters for solver
    double t, tf, rtol, atol, h0;
    std::complex<double> x, dx;
    int order;
    bool fo;
    RKSolver rksolver;
    // TODO: any way not to create all objects here?
    //WKBSolver wkbsolver
    WKBSolver * wkbsolver;
    WKBSolver1 wkbsolver1;
    WKBSolver2 wkbsolver2;
    WKBSolver3 wkbsolver3;

    public:
    // constructor
    Solution(de_system de_sys, std::complex<double> x0, std::complex<double>
    dx0, double t_i, double t_f, int o=3, double r_tol=1e-4, double a_tol=0.0,
    double h_0=1, bool full_output=false, bool interp=true);
    void solve();
    // stats
    int ssteps,totsteps,wkbsteps;
    std::vector<std::complex<double>> sol, dsol;
    std::vector<double> times;
    std::vector<bool> wkbs;
};


Solution::Solution(de_system de_sys, std::complex<double> x0,
std::complex<double> dx0, double t_i, double t_f, int o, double r_tol, double
a_tol, double h_0, bool full_output, bool interp){
    
    // Set parameters for solver
    x = x0;
    dx = dx0;
    t = t_i;
    tf = t_f;
    order = o;
    rtol = r_tol;
    atol = a_tol;
    h0 = h_0;
    fo = full_output;
    rksolver = RKSolver(de_sys);
    switch(order){
        case 1: wkbsolver1 = WKBSolver1(de_sys, order);
                wkbsolver = &wkbsolver1;
                break;
        case 2: wkbsolver2 = WKBSolver2(de_sys, order);
                wkbsolver = &wkbsolver2;
                break;
        case 3: wkbsolver3 = WKBSolver3(de_sys, order);
                wkbsolver = &wkbsolver3;
                break;
    };
};

void Solution::solve(){ 
    
    int nrk, nwkb1, nwkb2;
    // Settings for MS
    nrk = 5;
    nwkb1 = 1;
    nwkb2 = 8;
    Eigen::Matrix<std::complex<double>,2,4> rkstep;
    Eigen::Matrix<std::complex<double>,3,2> wkbstep;
    Eigen::Matrix<std::complex<double>,1,2> rkx, wkbx;
    Eigen::Matrix<std::complex<double>,1,2> rkerr, wkberr, truncerr;
    Eigen::Matrix<double,1,2> rkdeltas; 
    Eigen::Matrix<double,1,4> wkbdeltas;
    double tnext, hnext, h, hrk, hwkb;
    double wkbdelta, rkdelta;
    std::complex<double> xnext, dxnext;
    bool wkb = false;
    Eigen::Index maxindex;
    h = h0;
    tnext = t+h;
    // Initialise stats
    sol.emplace_back(x);
    dsol.emplace_back(dx);
    times.emplace_back(t);
    ssteps = 0;
    totsteps = 0;
    wkbsteps = 0;
    
    while(t < tf){
        // Check if we are reaching the end of integration
        if(tnext>tf){
            h = tf - t;
            tnext = tf;
        }

        // Keep updating stepsize until step is accepted
        while(true){
            // RK step
            rkstep = rksolver.step(x, dx, t, h);
            rkx << rkstep(0,0), rkstep(0,1);
            rkerr << rkstep(1,0), rkstep(1,1);
            // WKB step
            wkbstep = wkbsolver->step(x, dx, t, h, rksolver.ws, rksolver.gs, rksolver.ws5, rksolver.gs5);
            wkbx = wkbstep.row(0);
            wkberr = wkbstep.row(2);
            truncerr = wkbstep.row(1);
            // dominant error calculation
            wkbdeltas << std::abs(truncerr(0))/std::abs(wkbx(0)),
            std::abs(truncerr(1))/std::abs(wkbx(1)),
            std::abs(wkberr(0))/std::abs(wkbx(0)),
            std::abs(wkberr(1))/std::abs(wkbx(1));
            rkdeltas << std::abs(rkerr(0))/std::abs(rkx(0)), std::abs(rkerr(1))/std::abs(rkx(1));
            wkbdelta = std::max(1e-10, wkbdeltas.maxCoeff(&maxindex));
            rkdelta = std::max(1e-10, rkdeltas.maxCoeff()); 
            // predict next stepsize 
            hrk = h*std::pow((rtol/rkdelta),1.0/nrk);
            if(maxindex<=1)
                hwkb = h*std::pow(rtol/wkbdelta,1.0/nwkb1);
            else
                hwkb = h*std::pow(rtol/wkbdelta,1.0/nwkb2);
            // choose step with larger predicted stepsize
            if(hwkb >= hrk)
                wkb = true;
            else
                wkb = false;
            if(wkb){
                xnext = wkbx(0);
                dxnext = wkbx(1);
                // if wkb step chosen, ignore truncation error in
                // stepsize-increase
                wkbdelta = std::max(1e-10, std::abs(wkbdeltas.tail(2).maxCoeff()));
                hnext = h*std::pow(rtol/wkbdelta,1.0/nwkb2);
            }
            else{
                xnext = rkx(0);
                dxnext = rkx(1);
                hnext = hrk;
            };
            totsteps += 1;
            // check if chosen step was successful
            if(hnext>=h){
                //std::cout << "wkb: " << wkb << ", t: "<< tnext << ", x: " << xnext << ", dx: " << dxnext << std::endl;
                sol.emplace_back(xnext);
                dsol.emplace_back(dxnext);
                times.emplace_back(tnext);
                tnext += hnext;
                x = xnext;
                dx = dxnext;
                t += h;
                h = hnext;
                ssteps +=1;
                if(wkb){
                    wkbsteps +=1;
                    wkbs.emplace_back(true);
                }
                else
                    wkbs.emplace_back(false);
                break;
            }
            else{
                if(wkb){
                    if(maxindex<=1){
                        if(nwkb1 > 1)
                            hnext = h*std::pow(rtol/wkbdelta,1.0/(nwkb1-1));
                        else
                            hnext = 0.95*h*rtol/wkbdelta;
                    }
                    else
                        hnext = h*std::pow(rtol/wkbdelta,1.0/(nwkb2-1));
                }
                else
                    hnext = h*std::pow(rtol/rkdelta,1.0/(nrk-1));
                h = hnext;
                tnext = t + hnext;
            };
        };
    };

    // Write output to file if prompted
    if(fo){
        std::string output;
        std::cout << "Please enter filename to print output to: " << std::endl;
        std::cin >> output;
        std::ofstream f;
        f.open(output);
        f << "# Summary:\n# total steps taken: " + std::to_string(totsteps) +
        "\n# of which successful: " + std::to_string(ssteps) + "\n# of which"+
        +"wkb: " + std::to_string(wkbsteps) + "\n# time, x, dx, wkb, Ai(-t)+i*Bi(-t)\n";
        for(int i=0; i<=ssteps; i++)
            f << std::setprecision(20) << times[i] << " " <<
            std::setprecision(20) << sol[i] << " " << std::setprecision(20) <<
            dsol[i] << " " << wkbs[i] 
            //<< " " << std::setprecision(20) << std::complex<double>(boost::math::airy_ai(-times[i]), boost::math::airy_bi(-times[i])) 
            //<< " " << std::setprecision(20) << -std::complex<double>(boost::math::airy_ai_prime(-times[i]), boost::math::airy_bi_prime(-times[i]))
            << "\n"; 
        f.close();
    }
    
};
