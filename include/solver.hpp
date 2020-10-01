#pragma once
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <iomanip>
#include <limits>
#include <vector>
#include "interpolator.hpp"
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
    // Underlying de_system object
    de_system *de_sys_;
    // Related to direction of integration
    double fend, fnext;

    public:
    // constructor
    RKSolver rksolver;
    Solution(de_system &de_sys, std::complex<double> x0, std::complex<double>
    dx0, double t_i, double t_f, int o=3, double r_tol=1e-4, double a_tol=0.0,
    double h_0=1, const char* full_output="");

    template<typename X = double> Solution(de_system &de_sys,
    std::complex<double> x0, std::complex<double> dx0, double t_i, double t_f,
    const X &do_times, int o=3, double r_tol=1e-4, double a_tol=0.0,
    double h_0=1, const char* full_output="");

    void solve();
    // stats
    int ssteps,totsteps,wkbsteps;
    std::list<std::complex<double>> sol, dsol;
    std::list<double> times;
    std::list<bool> wkbs;
    // dense output
    std::list<double> dotimes, dotimes_rk;
    std::list<std::complex<double>> dosol, dodsol, dosol_rk, dodsol_rk;
    std::list<double>::iterator dotit;

};


Solution::Solution(de_system &de_sys, std::complex<double> x0,
std::complex<double> dx0, double t_i, double t_f, int o, double r_tol, double
a_tol, double h_0, const char* full_output){
    
    // Make underlying equation system accessible
    de_sys_ = &de_sys;

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
    rksolver = RKSolver(*de_sys_);

    // Determine direction of integration, fend>0 and integration ends when
    // it crosses zero
    if((t>=tf) and h0<0){
        // backwards
        fend = t-tf;
        fnext = fend;
        de_sys_->Winterp.sign_ = 0;
        de_sys_->Ginterp.sign_ = 0;

    }
    else if((t<=tf) and h0>0){
        // forward
        fend = tf-t;
        fnext = fend;
        de_sys_->Winterp.sign_ = 1;
        de_sys_->Ginterp.sign_ = 1;
    }
    else{
        throw "Direction of integration in conflict with direction of initial step, terminating. Please check your values for ti, tf, and h. ";
        return;
    }

    // No dense output desired if this constructor was called, so only output
    // answer at t_i and t_f
    dotimes.push_back(t_i);
    dotimes.push_back(t_f);
    dosol.push_back(x0);
    dodsol.push_back(dx0);
    dotit = dotimes.end();

    
    switch(order){
        case 1: wkbsolver1 = WKBSolver1(*de_sys_, order);
                wkbsolver = &wkbsolver1;
                break;
        case 2: wkbsolver2 = WKBSolver2(*de_sys_, order);
                wkbsolver = &wkbsolver2;
                break;
        case 3: wkbsolver3 = WKBSolver3(*de_sys_, order);
                wkbsolver = &wkbsolver3;
                break;
    };
};

template<typename X> Solution::Solution(de_system &de_sys, std::complex<double> x0,
std::complex<double> dx0, double t_i, double t_f, const X &do_times, int o, double r_tol, double
a_tol, double h_0, const char* full_output){

    // Make underlying equation system accessible
    de_sys_ = &de_sys;
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
    rksolver = RKSolver(*de_sys_);

    // Determine direction of integration, fend>0 and integration ends when
    // it crosses zero
    if((t>=tf) and h0<0){
        // backwards
        fend = t-tf;
        fnext = fend;
        de_sys_->Winterp.sign_ = 0;
        de_sys_->Ginterp.sign_ = 0;
    }
    else if((t<=tf) and h0>0){
        // forward
        fend = tf-t;
        fnext = fend;
        de_sys_->Winterp.sign_ = 1;
        de_sys_->Ginterp.sign_ = 1;
    }
    else{
        throw "Direction of integration in conflict with direction of initial step, terminating. Please check your values for ti, tf, and h. ";
        return;
    }

    // Dense output checks: 
    int dosize = do_times.size();
    dotimes.resize(dosize);
    dosol.resize(dosize);
    dodsol.resize(dosize);
    int docount = 0;
    auto doit = do_times.begin();
    if(de_sys_->Winterp.sign_ == 1){
                for(auto it=dotimes.begin(); it!=dotimes.end(); it++){
            *it = *doit;
            docount++; doit++;
        }
    }
    else{
         for(auto it=dotimes.rbegin(); it!=dotimes.rend(); ++it){
            *it = *doit;
            docount++; ++doit;
        }
    }
    dotit = dotimes.begin();
    switch(order){
        case 1: wkbsolver1 = WKBSolver1(*de_sys_, order);
                wkbsolver = &wkbsolver1;
                break;
        case 2: wkbsolver2 = WKBSolver2(*de_sys_, order);
                wkbsolver = &wkbsolver2;
                break;
        case 3: wkbsolver3 = WKBSolver3(*de_sys_, order);
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
    Eigen::Matrix<double,1,2> errmeasure_rk; 
    Eigen::Matrix<double,1,4> errmeasure_wkb;
    double tnext, hnext, h, hrk, hwkb;
    double wkbdelta, rkdelta;
    std::complex<double> xnext, dxnext;
    bool wkb = false;
    Eigen::Index maxindex_wkb, maxindex_rk;
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
    // Dense output
    std::list<double> inner_dotimes;
    std::list<std::complex<double>> inner_dosols;
    auto it_dosol = dosol.begin();
    Eigen::Matrix<std::complex<double>,1,2> y_dense_rk;
    std::complex<double> x_dense_rk, dx_dense_rk;

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
            // Error scale measures
            errmeasure_rk << std::abs(rkerr(0))/(std::abs(rkx(0))*rtol+atol), std::abs(rkerr(1))/(std::abs(rkx(1))*rtol+atol);
            errmeasure_wkb << std::abs(truncerr(0))/(std::abs(wkbx(0))*rtol+atol),
            std::abs(truncerr(1))/(std::abs(wkbx(1))*rtol+atol),
            std::abs(wkberr(0))/(std::abs(wkbx(0))*rtol+atol),
            std::abs(wkberr(1))/(std::abs(wkbx(1))*rtol+atol);
            rkdelta = std::max(1e-10, errmeasure_rk.maxCoeff(&maxindex_rk)); 
            if(std::isnan(errmeasure_wkb.maxCoeff())==false &&
               std::isinf(std::real(wkbx(0)))==false &&
               std::isinf(std::imag(wkbx(0)))==false &&
               std::isinf(std::real(wkbx(1)))==false &&
               std::isinf(std::imag(wkbx(1)))==false &&
               std::isnan(std::real(wkbx(0)))==false &&
               std::isnan(std::imag(wkbx(0)))==false &&
               std::isnan(std::real(wkbx(1)))==false &&
               std::isnan(std::imag(wkbx(1)))==false){
                wkbdelta = std::max(1e-10, errmeasure_wkb.maxCoeff(&maxindex_wkb));
            }
            else{
                wkbdelta = std::numeric_limits<double>::infinity();
            }

            // predict next stepsize 
            hrk = h*std::pow((1.0/rkdelta),1.0/nrk);
            if(maxindex_wkb<=1)
                hwkb = h*std::pow(1.0/wkbdelta,1.0/nwkb1);
            else
                hwkb = h*std::pow(1.0/wkbdelta,1.0/nwkb2);
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
                wkbdelta = std::max(1e-10, errmeasure_wkb.tail(2).maxCoeff());
                hnext = h*std::pow(1.0/wkbdelta,1.0/nwkb2);
            }
            else{
                xnext = rkx(0);
                dxnext = rkx(1);
                hnext = hrk;
            };
            totsteps += 1;
            // Checking for too many steps and low acceptance ratio:
            if(totsteps % 5000 == 0){
                std::cerr << "Warning: the solver took " << totsteps << " steps, and may take a while to converge." << std::endl; 
                if(ssteps/totsteps < 0.05){
                    std::cerr << "Warning: the step acceptance ratio is below 5%, the solver may take a while to converge." << std::endl;
                }
            }

            // check if chosen step was successful
            if(std::abs(hnext)>=std::abs(h)){
                if(dotit!=dotimes.end()){
                    while((*dotit-t>=0 && tnext-*dotit>=0) or (*dotit-t<=0 && tnext-*dotit<=0)){
                        inner_dotimes.push_back(*dotit);
                        dotit++;
                    }
                    if(inner_dotimes.size() > 0){
                        inner_dosols.resize(inner_dotimes.size());
                        if(wkb){
                            // Dense output after successful WKB step
                            wkbsolver->dense_step(t,inner_dotimes,inner_dosols);
                        }
                        else{
                            // Dense output after successful RK step
                            for(auto it=inner_dotimes.begin(); it!=inner_dotimes.end(); it++)
                            rksolver.dense_step(t,h,x,inner_dotimes,inner_dosols);
                        }
                   }
                }
                auto inner_it=inner_dosols.begin();
                while(inner_it!=inner_dosols.end() && it_dosol!=dosol.end()){
                    *it_dosol = *inner_it;
                    it_dosol++;
                    inner_it++;
                }
                inner_dotimes.resize(0);
                inner_dosols.resize(0);
               
                // record type of step
                if(wkb){
                    wkbsteps +=1;
                    wkbs.push_back(true);
                }
                else{
                    wkbs.push_back(false);
                }

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
                // Update interpolation bounds
                de_sys_->Winterp.update_interp_bounds();
                de_sys_->Ginterp.update_interp_bounds();

                break;
            }
            else{
                if(wkb){
                    if(maxindex_wkb<=1){
                        if(nwkb1 > 1)
                            hnext = h*std::pow(1.0/wkbdelta,1.0/(nwkb1-1));
                        else
                            hnext = 0.95*h*1.0/wkbdelta;
                    }
                    else
                        hnext = h*std::pow(1.0/wkbdelta,1.0/(nwkb2-1));
                }
                else
                    hnext = h*std::pow(1.0/rkdelta,1.0/(nrk-1));
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

    // If integrating backwards, reverse dense output (because it will have been
    // reversed at the start)
    if(de_sys_->Winterp.sign_ == 0){
        dosol.reverse();
    }

    // Write output to file if prompted
    if(not (*fo==0)){
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
            f << std::setprecision(15) << *it_t << ";" <<
            std::setprecision(15) << *it_x << ";" << std::setprecision(15) <<
            *it_dx << ";" << *it_w << "\n"; 
            ++it_t;
            ++it_x;
            ++it_dx;
            ++it_w;
        }
        // print all dense output to file
        int dosize = dosol.size();
        auto it_dosol = dosol.begin();
        auto it_dotimes = dotimes.begin();
        for(int i=0; i<dosize; i++){
            f << std::setprecision(20) << *it_dotimes << ";" << std::setprecision(20) << *it_dosol << ";;\n";
            ++it_dosol;
            ++it_dotimes;
        }

        f.close();
    }
    
};
