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
    
    class Method
    {
        protected:

        public:
        // default constructor
        Method();  
        // class functions
        virtual Step step(Vectorfn F, Vector y, double h)=0;
    };
    
    Method::Method(){
        // Default constructor of a Method (does nothing)
    };

    ////////////////////////////////////RKF/////////////////////////////////////
    
    class RKFsolver : public Method
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
    
    class WKBsolver : public Method
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
        //Vector S, dS, ddS;
        Scalar fp, fm, dfp, dfm, ddfp, ddfm;
        Scalar ap, am, bp, bm;
        Vector y0, y1; // solution vector at beginning and end of current step
        Vector error0, error1; // error on solution vector at beginning and end of current step  

        // class methods
        virtual Vector ddS(Vector y);
        virtual Vector dS(Vector y);
        virtual Vector S_odd(Vector y);
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
        // Background y and its error before step 
//        Vector y0_bg = y0.segment(2, d-2-(order+2)/2);
//        Vector error0_bg = error0.segment(2, d-2-(order+2)/2);
        // Background after step
//        Vector y1_bg = y1.segment(2, d-2-(order+2)/2);
        
//        Scalar ddx = -y0(0)*std::pow(sys.w(y0_bg),2) - 2.0*y0(1)*sys.g(y0_bg);
        Vector ds = dS(y);
        std::cout << "dS(y): " << ds << std::endl;




        Step result(d);
        result.wkb = true; 
        return result;
    };

    Vector WKBsolver::ddS(Vector y){
        return y;
    };
    
    Vector WKBsolver::dS(Vector y){
        return y;
    };

    Vector WKBsolver::S_odd(Vector y){
        return y;
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
        Vector S_odd(Vector y);

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
        Vector result(2);
        result << std::complex<double>(0.0, 1.0)*sys.w(y), -sys.g(y)-0.5*sys.dw(y)/sys.w(y);
        return result;
    };

    Vector WKBsolver1::S_odd(Vector y){
        Vector result(1);
        result << -0.5*std::log(sys.w(y));
        return result;
    };
    
    class WKBsolver2 : public WKBsolver
    {
    };
    
    class WKBsolver3 : public WKBsolver
    {
    };



}
