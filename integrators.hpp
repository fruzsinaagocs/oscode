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
        Scalar x, dx, ddx;
        Vector ds, dds, s;
        Scalar fp, fm, dfp, dfm, ddfp, ddfm;
        Scalar ap, am, bp, bm;
        Vector y0, y1; // solution vector at beginning and end of current step
        Vector error0, error1; // error on solution vector at beginning and end of current step  

        // class methods
        virtual Vector ddS(Vector y); // to be overridden by higher order derived classes
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
        Vector y0_bg = y0.segment(2, d-2-(order+2)/2); // Background y and its error before step 
        Vector y1_bg = y1.segment(2, d-2-(order+2)/2); // Background after step
        ddx = -y0(0)*std::pow(sys.w(y0_bg),2) - 2.0*y0(1)*sys.g(y0_bg);
        dx = y0(1);
        x = y0(0);
    
        s = Vector::Zero(order+1);
        ds = dS(y0_bg);
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
        std::cout << "at y0: " << std::endl;
        std::cout << "ddx: " << ddx << ", dx: " << dx << ", x: " << x << std::endl;
        std::cout << "ds: " << ds << std::endl;
        std::cout << "fp " << fp << ", fm: " << fm << std::endl;
        std::cout << "dfp: " << dfp << ", dfm: " << dfm << std::endl;
        std::cout << "ddfp: " << ddfp << ", ddfm: " << ddfm << std::endl;
        std::cout << "ap: " << ap << ", am: " << am << std::endl;
        std::cout << "bp: " << bp << ", bm: " << bm << std::endl;

        Vector s_odd = S_odd(y1_bg) - S_odd(y0_bg);
        for(int i=0; i<(order+1); i++){
            if(i%2==0)
                s(i) = y1(d-(order+2)/2 + i/2);
            else
                s(i) = s_odd(i/2); 
        };
        std::cout << "s at y1: " << s << std::endl;
        ds = dS(y1_bg);
        dds = ddS(y1_bg);
        fp = Fp();
        fm = Fm();
        dfp = dFp();
        dfm = dFm();
        std::cout << "at y1: " << std::endl;
        std::cout << "ds: " << ds << std::endl;
        std::cout << "fp " << fp << ", fm: " << fm << std::endl;
        std::cout << "dfp: " << dfp << ", dfm: " << dfm << std::endl;
        x = ap*fp + am*fm;
        dx = bp*dfp + bm*dfm;

//        Vector error0_bg = error0.segment(2, d-2-(order+2)/2);

        Step result(d);
        result.y << x, dx, y1.tail(d-2);
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
        result << std::complex<double>(0.0, 1.0)*sys.w(y);
        return result;
    };

    Vector WKBsolver::S_odd(Vector y){
        Vector result;
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
