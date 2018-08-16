#include <eigen3/Eigen/Dense>
#include <iostream>

// complex numbers, vectors, matrices
typedef Eigen::VectorXcd V;
typedef V::Scalar S;
typedef Eigen::MatrixXcd M;

// complex-valued functions
typedef S (* Sfn)(V);
typedef V (* Vfn)(V);
typedef M (* Mfn)(V);

class de_system
{
    private:

    public:
        // default constructor
        de_system();        
        // constructor overloads
        de_system(Vfn);//, Mfn, Sfn, Vfn, Mfn, Sfn, Vfn, Mfn);
        
        // class data
        
        // class functions
        Vfn F;
        Mfn DF;
        Sfn w;
        Vfn Dw;
        Mfn DDw;
        Sfn g;
        Vfn Dg;
        Mfn DDg;

};

de_system::de_system(){
    // Default constructor of de_system (does nothing)
}

de_system::de_system(Vfn f){//, Mfn Df, Sfn freq, Vfn Dfreq, Mfn DDfreq, Sfn gam, Vfn Dgam, Mfn DDgam){
    // default constructor for a system of differential equations

    F = f;
    //DF = Df;
    //w = freq;
    //Dw = Dfreq;
    //DDw = DDfreq;
    //g = gam;
    //Dg = Dgam;
    //DDg = DDgam;
    std::cout << "Constructed a de_system object" << std::endl;
}





