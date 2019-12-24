#pragma once
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <iomanip>
#include <limits>
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
    const char* fo;
    WKBSolver * wkbsolver;
    WKBSolver1 wkbsolver1;
    WKBSolver2 wkbsolver2;
    WKBSolver3 wkbsolver3;

    public:
    // constructor
    RKSolver rksolver;
    Solution(de_system &de_sys, std::complex<double> x0, std::complex<double>
    dx0, double t_i, double t_f, int o=3, double r_tol=1e-4, double a_tol=0.0,
    double h_0=1, const char* full_output="");
    void solve();
    // stats
    int ssteps,totsteps,wkbsteps;
    std::list<std::complex<double>> sol, dsol;
    std::list<double> times;
    std::list<bool> wkbs;
};


Solution::Solution(de_system &de_sys, std::complex<double> x0,
std::complex<double> dx0, double t_i, double t_f, int o, double r_tol, double
a_tol, double h_0, const char* full_output){
    
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
    nwkb1 = 2;
    nwkb2 = 4;
    Eigen::Matrix<std::complex<double>,2,2> rkstep;
    Eigen::Matrix<std::complex<double>,3,2> wkbstep;
    Eigen::Matrix<std::complex<double>,1,2> rkx, wkbx;
    Eigen::Matrix<std::complex<double>,1,2> rkerr, wkberr, truncerr;
    Eigen::Matrix<double,1,2> rkdeltas; 
    Eigen::Matrix<double,1,4> wkbdeltas;
    double fend, fnext, tnext, hnext, h, hrk, hwkb;
    double wkbdelta, rkdelta;
    std::complex<double> xnext, dxnext;
    bool wkb = false;
    Eigen::Index maxindex;
    h = h0;
    tnext = t+h;
    // Initialise stats
    sol.push_back(x);
    dsol.push_back(dx);
    times.push_back(t);
    wkbs.push_back(false);
    ssteps = 0;
    totsteps = 0;
    wkbsteps = 0;
    // Determine direction of integration, fend>0 and integration ends when
    // it crosses zero
    if((t>=tf) and h<0){
        // backwards
        fend = t-tf;
        fnext = fend;
    }
    else if((t<=tf) and h>0){
        // forward
        fend = tf-t;
        fnext = fend;
    }
    else{
        throw "Direction of integration in conflict with direction of initial step, terminating. Please check your values for ti, tf, and h. ";
        return;
    }

    while(fend > 0){
        // Check if we are reaching the end of integration
        if(fnext < 0){
            h = tf - t;
            tnext = tf;
        };

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
            // Safety feature for when all wkb steps are 0 (truncer=0), but not
            // necessarily in good WKB regime:
            truncerr(0) = std::max(1e-10,abs(truncerr(0)));
            truncerr(1) = std::max(1e-10,abs(truncerr(1)));
            // dominant error calculation
            wkbdeltas << std::abs(truncerr(0))/std::abs(wkbx(0)),
            std::abs(truncerr(1))/std::abs(wkbx(1)),
            std::abs(wkberr(0))/std::abs(wkbx(0)),
            std::abs(wkberr(1))/std::abs(wkbx(1));
            rkdeltas << std::abs(rkerr(0))/std::abs(rkx(0)), std::abs(rkerr(1))/std::abs(rkx(1));
            rkdelta = std::max(1e-10, rkdeltas.maxCoeff()); 
            if(std::isnan(wkbdeltas.maxCoeff())==false && std::isinf(std::real(wkbx(0)))==false && std::isinf(std::real(wkbx(1)))==false)
                wkbdelta = std::max(1e-10, wkbdeltas.maxCoeff(&maxindex));
            else
                wkbdelta = std::numeric_limits<double>::infinity();

            // predict next stepsize 
            hrk = h*std::pow((rtol/rkdelta),1.0/nrk);
            if(maxindex<=1)
                hwkb = h*std::pow(rtol/wkbdelta,1.0/nwkb1);
            else
                hwkb = h*std::pow(rtol/wkbdelta,1.0/nwkb2);
            // choose step with larger predicted stepsize
            if(std::abs(hwkb) >= std::abs(hrk)){
                wkb = true;
            }
            else{
                wkb = false;
            }
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
            if(std::abs(hnext)>=std::abs(h)){
                sol.push_back(xnext);
                dsol.push_back(dxnext);
                times.push_back(tnext);
                tnext += hnext;
                x = xnext;
                dx = dxnext;
                t += h;
                h = hnext;
                if(h>0){
                    fend=tf-t;
                    fnext=tf-tnext;
                }
                else{
                    fend=t-tf;
                    fnext=tnext-tf;
                };
                ssteps +=1;
                if(wkb){
                    wkbsteps +=1;
                    wkbs.push_back(true);
                }
                else
                    wkbs.push_back(false);
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
                if(h>0){
                    fnext=tf-tnext;
                }
                else{
                    fnext=tnext-tf;
                };
            };
        };
    };

    // Write output to file if prompted
    if(not *fo==0){
        std::string output(fo);
        std::ofstream f;
        f.open(output);
        f << "# Summary:\n# total steps taken: " + std::to_string(totsteps) +
        "\n# of which successful: " + std::to_string(ssteps) + "\n# of which"+
        +" wkb: " + std::to_string(wkbsteps) + "\n# time, sol, dsol, wkb? (type)\n";
        auto it_t = times.begin();
        auto it_w = wkbs.begin();
        auto it_x = sol.begin();
        auto it_dx = dsol.begin();
        for(int i=0; i<=ssteps; i++){
            f << std::setprecision(15) << *it_t << " " <<
            std::setprecision(15) << *it_x << " " << std::setprecision(15) <<
            *it_dx << " " << *it_w << "\n"; 
            ++it_t;
            ++it_x;
            ++it_dx;
            ++it_w;
        };
        f.close();
    }
    
};
