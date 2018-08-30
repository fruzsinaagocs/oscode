# include "system.hpp"
namespace RKWKB{
    
    struct Step
    {
        // data structure to store the result of a step and its error.
        Vector y, error;
        bool wkb;
        // default constructor does nothing
        Step(){};
        // constructor takes in size and resizes y, error
        Step(int d){
            y.resize(d);
            error.resize(d);
        };
    };
    
    class Integrator
    {
        protected:

        public:
        // default constructor
        Integrator();  
        // class functions
        virtual Step step(Vectorfn F, Vector y, double h)=0;
    };
    
    Integrator::Integrator(){
        // Default constructor of a Integrator (does nothing)
    };

    ////////////////////////////////////RKF/////////////////////////////////////
    
    class RKFsolver : public Integrator
    {
        private:
    
        public:
        // default constructor
        RKFsolver(); 
        // constructor overloads
        RKFsolver(de_system);
    
        // class data
        de_system sys;
        Matrix butcher_a;
        Vector butcher_c, butcher_b4, butcher_b5, butcher_r;
        
        // class functions
        Step step(Vectorfn F, Vector y, double h);
    
    };
    
    RKFsolver::RKFsolver(){
        // Default constructor of RKFsolver (does nothing)
    };
    
    RKFsolver::RKFsolver(de_system system){
        // Constructor for an RKFsolver from a system of differential equations
        sys = system;
        
        butcher_a.resize(5,5);
        butcher_c.resize(6);
        butcher_b4.resize(6);
        butcher_b5.resize(6);
        butcher_r.resize(6);
        
        butcher_a << 1.0/4.0, 0.0, 0.0, 0.0, 0.0,
        			3.0/32.0, 9.0/32.0, 0.0, 0.0, 0.0,
        			1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0, 0.0, 0.0,
        			439.0/216.0, -8.0, 3680.0/513.0, -845.0/4104.0, 0.0,
        			-8.0/27.0, 2.0, -3544.0/2565.0, 1859.0/4104.0, -11.0/40.0;
        butcher_c << 0.0, 1.0/4.0, 3.0/8.0, 12.0/13.0, 1.0, 1.0/2.0;
        butcher_b4 << 25.0/216.0 , 0.0, 1408.0/2565.0, 2197.0/4104.0, -1.0/5.0 , 0.0;
        butcher_b5 << 16.0/135.0 , 0.0 , 6656.0/12825.0 , 28561.0/56430.0, -9.0/50.0, 2.0/55.0;
        butcher_r << 1.0/360.0, 0.0, -128/4275.0, -2197/75240.0, 1/50.0, 2/55.0; 
    
    };
    
    Step RKFsolver::step(Vectorfn F, Vector y, double h){
        // Stepper function using the RKF method
        
        int d = y.size();
        Matrix k(6,d);
        Step result(d);
        Vector yarg = y;

        k.row(0) = h*F(y);
        for(int s=1; s<=5; s++){
            yarg = y;
            for(int i=0; i<=(s-1); i++)
                yarg += butcher_a(s-1, i)*k.row(i);
            k.row(s) = h*F(yarg);
        } 
        result.y = y;
        result.error = Vector::Zero(d);
        for(int j=0; j<=5; j++){
            result.y += butcher_b5(j)*k.row(j);
            result.error += (butcher_b4(j) - butcher_b5(j))*k.row(j);
        }   
        result.wkb = false;
        
        return result;
    };
    
    ////////////////////////////////////WKB////////////////////////////////////
    
    class WKBsolver : public Integrator
    {
        private:
    
        public:
        // default constructor
        WKBsolver();    
        // constructor overloads
        WKBsolver(de_system system, int o=1);
    
        // class data
        de_system sys;
        int order;
        Scalar x, dx, ddx;
        Vector ds, ds_all, dds, s, s_all;
        Scalar fp, fm, dfp, dfm, ddfp, ddfm;
        Scalar ap, am, bp, bm;
        Vector y0, y1; // solution vector at beginning and end of current step
        Vector error0, error1; // error on solution vector at beginning and end of current step  
        Vector trunc_error;

        // for error estimation           
        Scalar error_x, error_dx, error_ddx, error_ap, error_am, error_bp, error_bm, error_fp, error_fm, error_dfp, error_dfm;
        Vector error_s;
        Scalar t_error_ap, t_error_am, t_error_bp, t_error_bm, t_error_fp, t_error_fm, t_error_dfp, t_error_dfm, t_error_ddfp, t_error_ddfm; 

        // class methods
        virtual Vector ddS(Vector y); // to be overridden by higher order WKBs
        virtual Vector dS(Vector y);
        void truncation_error(bool);
        void error(bool, Vector);
        Step step(Vectorfn F, Vector y, double h);
        Scalar Fp();
        Scalar Fm();
        Scalar dFp();
        Scalar dFm(); 
        Scalar ddFp();
        Scalar ddFm();
        Scalar Ap();
        Scalar Am();
        Scalar Bp();
        Scalar Bm();
    };
    
    WKBsolver::WKBsolver(){
        // Default constructor for a WKBsolver (does nothing)
    };
    
    WKBsolver::WKBsolver(de_system system, int o){
        // Constructor for a WKBsolver from a system of differential equations
        sys = system;
        order = o;
    };
    
    Step WKBsolver::step(Vectorfn F, Vector y, double h){
        // Stepper function using the RKF method   
    
        int d = y0.size();
        Vector y0_bg = y0.segment(2, d-2-(order+2)); // Background y and its error before step 
        Vector bg_error0 = error0.segment(2, d-2-(order+2));
        Vector y1_bg = y1.segment(2, d-2-(order+2)); // Background and its error after step
        Vector bg_error1 = error1.segment(2, d-2-(order+2));

        ddx = -y0(0)*std::pow(sys.w(y0_bg),2) - 2.0*y0(1)*sys.g(y0_bg);
        dx = y0(1);
        x = y0(0);
    
        s_all = Vector::Zero(order+2);
        s = s_all.head(order+1);
        ds_all = dS(y0_bg);
        ds = ds_all.head(order+1);
        dds = ddS(y0_bg);
        fp = 1.0;
        fm = 1.0;
        dfp = dFp();
        dfm = dFm();
        ddfp = ddFp();
        ddfm = ddFm(); 
        ap = Ap();
        am = Am();
        bp = Bp();
        bm = Bm();

        // Error on quantities at y0
        error(true, y0_bg);
        truncation_error(true);

        for(int i=0; i<(order+2); i++)
            s_all(i) = y1(d-(order+2) + i);
       
        s = s_all.head(order+1);
        ds_all = dS(y1_bg);
        ds = ds_all.head(order+1);
        dds = ddS(y1_bg);
        fp = Fp();
        fm = Fm();
        dfp = dFp();
        dfm = dFm();
        x = ap*fp + am*fm;
        dx = bp*dfp + bm*dfm;

        // Error on quantities at y1
        error(false, y1_bg); 
        truncation_error(false);
        
        Step result(d);
        result.y << x, dx, y1.tail(d-2);
        result.error << error1;
        result.wkb = true; 
        return result;
    };

    Vector WKBsolver::ddS(Vector y){
        Vector result(1);
        result << std::complex<double>(0.0, 1.0)*sys.dw(y);
        return result;
    };
    
    Vector WKBsolver::dS(Vector y){
        Vector result(1);
        result << std::complex<double>(0.0, 1.0)*sys.w(y), -sys.g(y) - 0.5*sys.dw(y)/sys.w(y);
        return result;
    };

    Scalar WKBsolver::Fp(){
        return std::exp(s.sum());
    };
    
    Scalar WKBsolver::Fm(){
        return std::conj(fp);
    };

    Scalar WKBsolver::dFp(){
        return ds.sum()*fp;
    };

    Scalar WKBsolver::dFm(){
        return std::conj(dfp); 
    };

    Scalar WKBsolver::ddFp(){
        return (dds.sum() + std::pow(ds.sum(),2))*fp;
    };

    Scalar WKBsolver::ddFm(){
        return std::conj(ddfp);
    };

    Scalar WKBsolver::Ap(){
        return (dx - x*dfm)/(dfp - dfm);
    };

    Scalar WKBsolver::Am(){
        return (dx - x*dfp)/(dfm - dfp);;
    };

    Scalar WKBsolver::Bp(){
        return (ddx*dfm - dx*ddfm)/(ddfp*dfm - ddfm*dfp);
    };

    Scalar WKBsolver::Bm(){
        return (ddx*dfp - dx*ddfp)/(ddfm*dfp - ddfp*dfm);
    };

    void WKBsolver::error(bool before, Vector y0_bg){
        // Error or on quantities at y0
        if(before){
            error_x = error0(0);
            error_dx = error0(1);
            error_ddx = -2.0*sys.g(y0_bg)*error_dx - std::pow(sys.w(y0_bg),2)*error_x;
            error_ap = (error_dx - dfm*error_x)/(dfp - dfm);
            error_am = (error_dx - dfp*error_x)/(dfm - dfp);
            error_bp = (error_ddx - ddfm*error_dx)/(ddfp*dfm - ddfm*dfp);
            error_bm = (error_ddx - ddfp*error_dx)/(ddfm*dfp - ddfp*dfm);
        }
        // Error on quantities at y1
        else{
            error_s = Vector::Zero(order+1);
            for(int i=0; i<(order+1); i++)
                    error_s(i) = error1(y1.size() - (order+2) + i);
            error_fp = error_s.sum()*fp;
            error_fm = std::conj(error_s.sum())*fm;//
            error_dfp = ds.sum()*error_fp;
            error_dfm = std::conj(ds.sum())*error_fm;//
            error1(0) = error_ap*fp + error_am*fm + ap*error_fp + am*error_fm;
            error1(1) = error_bp*dfp + error_bm*dfm + bp*error_dfp + bm*error_dfm;
        };
    };
    
    void WKBsolver::truncation_error(bool before){
        trunc_error = Vector::Zero(2);
        // Error from truncation at y0
        if(before){
            t_error_dfp = ds_all(order+1);
            t_error_dfm = std::conj(t_error_dfp);
            t_error_ddfp = std::pow(ds_all(order+1),2)+2.0*ds_all(order+1)*dfp;
            t_error_ddfm = std::conj(t_error_ddfp);
            t_error_ap = -(am*t_error_dfm + ap*t_error_dfp)/(dfp - dfm);
            t_error_am = -(ap*t_error_dfp + am*t_error_dfm)/(dfm - dfp);
        }
        // At y1
        else{
            t_error_fp = ds_all(order+1)*fp;
            t_error_fm = std::conj(t_error_fp);
            t_error_dfp = s_all(order+1)*dfp + (1.0 + s_all(order+1))*ds_all(order+1)*fp;
            t_error_dfm = std::conj(t_error_dfp);
            trunc_error(0)  = t_error_ap*fp + t_error_am*fm + ap*t_error_fp + am*t_error_fm;
            trunc_error(1) = t_error_bp*dfp + t_error_bm*dfm + bp*t_error_dfp + bm*t_error_dfm;
        };
    };

    class WKBsolver1 : public WKBsolver
    {
        private:
         
        public:
        // default constructor
        WKBsolver1();
        // constructor overloads
        WKBsolver1(de_system, int o=1);
        // class functions
        Vector ddS(Vector y);
        Vector dS(Vector y);
    };

    WKBsolver1::WKBsolver1(){
    };

    WKBsolver1::WKBsolver1(de_system system, int o) : WKBsolver(system, o){
    };

    Vector WKBsolver1::ddS(Vector y){
        Vector result(2);
        result << std::complex<double>(0.0, 1.0)*sys.dw(y), -sys.g(y)-0.5*sys.ddw(y)/sys.w(y)+0.5*std::pow(sys.dw(y),2)/std::pow(sys.w(y),2);
        return result;
    };
    
    Vector WKBsolver1::dS(Vector y){
        Vector result(3);
        result << std::complex<double>(0.0, 1.0)*sys.w(y), -sys.g(y)-0.5*sys.dw(y)/sys.w(y), std::complex<double>(0.0, 1.0)*1.0/8.0*(-4.0*std::pow(sys.w(y)*sys.g(y),2)-4.0*std::pow(sys.w(y),2)*sys.dg(y)+3.0*std::pow(sys.dw(y),2)-2.0*sys.w(y)*sys.ddw(y))/std::pow(sys.w(y),3);
        return result;
    };
    
    class WKBsolver2 : public WKBsolver
    {
        private:
         
        public:
        // default constructor
        WKBsolver2();
        // constructor overloads
        WKBsolver2(de_system, int o=2);
        // class functions
        Vector ddS(Vector y);
        Vector dS(Vector y);
    };

    WKBsolver2::WKBsolver2(){
    };

    WKBsolver2::WKBsolver2(de_system system, int o) : WKBsolver(system, o){
    };

    Vector WKBsolver2::ddS(Vector y){
        Vector result(2);
        result << std::complex<double>(0.0, 1.0)*sys.dw(y), -sys.g(y)-0.5*sys.ddw(y)/sys.w(y)+0.5*std::pow(sys.dw(y),2)/std::pow(sys.w(y),2);
        return result;
    };
    
    Vector WKBsolver2::dS(Vector y){
        Vector result(2);
        result << std::complex<double>(0.0, 1.0)*sys.w(y), -sys.g(y)-0.5*sys.dw(y)/sys.w(y);
        return result;
    };

    class WKBsolver3 : public WKBsolver
    {
        private:
        
        public:
        // default constructor
        WKBsolver3();
        // constructor overloads
        WKBsolver3(de_system, int o=3);
        // class functions
        Vector ddS(Vector y);
        Vector dS(Vector y);
    };

    WKBsolver3::WKBsolver3(){
    };

    WKBsolver3::WKBsolver3(de_system system, int o) : WKBsolver(system, o){
    };

    Vector WKBsolver3::ddS(Vector y){
        Vector result(2);
        result << std::complex<double>(0.0, 1.0)*sys.dw(y), -sys.g(y)-0.5*sys.ddw(y)/sys.w(y)+0.5*std::pow(sys.dw(y),2)/std::pow(sys.w(y),2);
        return result;
    };
    
    Vector WKBsolver3::dS(Vector y){
        Vector result(2);
        result << std::complex<double>(0.0, 1.0)*sys.w(y), -sys.g(y)-0.5*sys.dw(y)/sys.w(y);
        return result;
    };
    
}
