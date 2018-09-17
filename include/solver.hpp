#include "integrators.hpp"
#include <fstream>
#include <boost/math/special_functions/airy.hpp>

namespace RKWKB{
    
    class Solution
    {
        private:
    
    
        public:
        // default constructor
        Solution();
        // constructor overloads
        Solution(de_system de_sys, Vector ic, double t_ini, event F_end, int o=1, double r_tol=1e-4, double a_tol=0.0, double stepsize=1.0, std::string output="output.txt");
    
        // class data 
        int d, order; // size of solution vector, order of WKB solution
        double t_i, rtol, h; // current stepsize
        double t; // current time
        Vector y; // current solution vector
        Vector error; // current error on solution vector
        bool wkb; // type of step taken, true for wkb and false for rkf
        Vector atol;
        de_system de_sys;
        RKFsolver rkfsolver;
        WKBsolver * wkbsolver;
        WKBsolver1 wkbsolver1;
        WKBsolver2 wkbsolver2;
        WKBsolver3 wkbsolver3;
        event f_end; // function s.t. solution terminates at f(y,t)=0
        std::string outputfile;
        int stepsok, stepsall;
        double waste; // stepsok: no. of successful steps, waste: 1 - stepsok/stepsall
    
        // class functions
        Vector F_tot(Vector); // gives time-derivative of all variables propagated
        Vectorfn f_tot;
        void evolve();
        Step step(Integrator * method, Vectorfn F, Vector Y, double H); 
        void write(std::string output);
        void update_h(double err, bool wkb, bool success);
        
    };
    
    Solution::Solution(){
        // default constructor for a solution object (does nothing)
    };
    
    Solution::Solution(de_system system, Vector ic, double t_ini, event F_end, int o, double r_tol, double a_tol, double stepsize, std::string output){
        // constructor for solution of a system of differential equations
        
        order = o;
        d = ic.size() + (order+2);
        y.resize(d);
        y << ic, Vector::Zero(order+2);
        t_i = t_ini;
        rtol = r_tol;
        atol = a_tol*Vector::Ones(y.size());
        h = stepsize;
        t = t_i;
        f_end = F_end;
        error = Vector::Zero(y.size());
        outputfile = output;
        f_tot = std::bind(&Solution::F_tot, this, std::placeholders::_1);
        stepsok = 0;
        stepsall = stepsok;

        de_sys = system;
        rkfsolver = RKFsolver(de_sys);
        switch(order){
            case 1: wkbsolver1 = WKBsolver1(de_sys,1);
                    wkbsolver = &wkbsolver1;
                    break;
            case 2: wkbsolver2 = WKBsolver2(de_sys,2);
                    wkbsolver = &wkbsolver2;
                    break;
            case 3: wkbsolver3 = WKBsolver3(de_sys,3);
                    wkbsolver = &wkbsolver3;
                    break;
        };
        //write(outputfile);
    };
                
    void Solution::evolve(){
        // function to compute full numerical solution of the de_system.
    
        bool sign = (f_end(y, t_i).real() <= 0.0);
        bool next_sign = sign;
        Vector y_next = y;
        Vector error_next = error;
        Vector scale = rtol*y + atol;
        double t_next = t;
        Step step_rkf, step_wkb;
        double end_error;
        wkb=false;
        double maxerr=0.0;
        Scalar maxerr_c = 0.0;

        //just for testing
        double tend = 100000;
                                
        while(true){
            // keep updating stepsize until step is successful
            while(true){
                y.tail(order+2) = Vector::Zero(order+2);
                error.tail(order+2) = Vector::Zero(order+2);
                wkbsolver->y0 = y;
                wkbsolver->error0 = error;
                step_rkf = step(&rkfsolver, f_tot, y, h);
                wkbsolver->y1 = step_rkf.y;
                wkbsolver->error1 = step_rkf.error;
                step_wkb = step(wkbsolver, de_sys.F, y, h);
                wkb = step_rkf.error.head(2).norm() > wkbsolver->trunc_error.norm();
                if(wkb){
                    y_next = step_wkb.y;
                    error_next = step_wkb.error;
                }
                else{
                    y_next = step_rkf.y;
                    error_next = step_rkf.error;
                }
                // check if step is accepted
                scale = rtol*y_next + atol;
                maxerr_c = error_next.cwiseProduct(scale.cwiseInverse()).cwiseAbs().maxCoeff();
                maxerr = std::real(maxerr_c);
                if(maxerr <= 1.0){
                    t_next += h;   
                    stepsok += 1;
                    stepsall += 1;
                    break;
                }
                else{
                    update_h(maxerr, wkb, false);
                    stepsall += 1;
                };
            };
            
            // check if f_end has changed sign
            next_sign = (f_end(y_next, t_next).real() <= 0.0);
            end_error = std::abs(f_end(y_next, t_next));
            if(next_sign != sign){
                h = (t_next - t)/2.0;
                t_next = t;
            }
            else{
                y = y_next;
                std::string stepc;
                if(wkb)
                    stepc = "wkb";
                else
                    stepc = "rkf";
                t = t_next;
                error = error_next;
                // update stepsize
                update_h(maxerr, wkb, true);
                write(outputfile);
                if(std::abs(end_error) < 1e-4){
                    waste = 1.0 - (double) stepsok/(double) stepsall;
                    break;
                };
            };
        };
    };

    Step Solution::step(Integrator * integrator, Vectorfn F, Vector Y, double H){
        // function to take a single step with a given method in the numerical solution of the
        // de_system, from y with stepsize h.
    
        Step s = integrator->step(F, Y, H);
        return s;
    };
    
    void Solution::update_h(double err, bool wkb, bool success){
        // updates stepsize.

        // TODO: deal with NaNs in the following line
        if(success){
            if(wkb)
                h*=std::pow(err, -1.0/(order+1));
            else
                h*=std::pow(err, -1/5.0);
        }
        else{
            if(wkb)
                h*=std::pow(err, -1.0/(order));
            else
                h*=std::pow(err, -1/4.0);
        };
    };

    void Solution::write(std::string output){
        // write current step to file

        std::ofstream fout;
        fout.open(output, std::ios_base::app); // appends at the end of outputfile
        fout << std::setprecision(10) << t << " ";
        for(int i=0; i<y.size(); i++)
            fout << y(i) << " ";
        for(int i=0; i<error.size(); i++) 
            fout << error(i) << " ";
        fout << wkb;
        fout << std::endl;
        fout.close();
    };

    Vector Solution::F_tot(Vector z){
        // time-derivative of all variables propagated. z = [x, x', y, S_i]
        Vector result(d);
        Vector z_bg = z.segment(2, d-2-(order+2));
        result << z(1), -z(0)*std::pow(de_sys.w(z_bg),2) -2.0*z(1)*de_sys.g(z_bg), de_sys.F(z_bg), wkbsolver->dS(z_bg);
        return result;
    };

} 



