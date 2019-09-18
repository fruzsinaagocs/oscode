========================================================================
oscode: Oscillatory ordinary differential equation solver
========================================================================

.. image::
    https://github.com/fruzsinaagocs/oscode/raw/master/pyoscode/images/spectra.gif
    :width: 800

:oscode: oscillatory ordinary differential equation solver
:Author: Fruzsina Agocs, Will Handley, Mike Hobson, and Anthony Lasenby
:Version: 1.0.1
:Homepage: https://github.com/fruzsinaagocs/oscode
:Documentation: https://oscode.readthedocs.io
:Paper: https://arxiv.org/abs/1906.01421

.. image:: https://readthedocs.org/projects/oscode/badge/?version=latest
    :target: https://oscode.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://badges.gitter.im/oscode-help/community.svg
    :target: https://gitter.im/oscode-help/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge
    :alt: Chat on gitter
 

``oscode`` is a C++ tool with a Python interface that solves **osc**\illatory
**o**\rdinary **d**\ifferential **e**\quations efficiently. It is designed to
deal with equations of the form

.. image:: 
    https://github.com/fruzsinaagocs/oscode/raw/master/pyoscode/images/oscillator.png

where |gamma| (friction term) and |omega| (frequency) can be given as

.. |gamma| image:: https://github.com/fruzsinaagocs/oscode/raw/master/pyoscode/images/gamma.png

.. |omega| image:: https://github.com/fruzsinaagocs/oscode/raw/master/pyoscode/images/omega.png

- *In C++*, explicit functions or sequence containers (Eigen::Vectors, arrays,
  std::vectors, lists),
- *In Python*, numpy.arrays.

``oscode`` makes use of an analytic approximation of x(t) embedded in a
stepping procedure to skip over long regions of oscillations, giving a reduction
in computing time. The approximation is valid when the frequency changes slowly
relative to the timescales of integration, it is therefore worth applying when
this condition holds for at least some part of the integration range. 

For the details of the numerical method used by ``oscode``, see the Citations
section.


Installation
------------

Dependencies
~~~~~~~~~~~~

Basic requirements for using the C++ interface:

- C++11 or later
- `Eigen <http://eigen.tuxfamily.org/index.php?title=Main_Page>`__

For using the Python interface, you will additionally need:

- Python 2.7 or 3.5+
- `numpy <https://pypi.org/project/numpy/>`__
- `scipy <https://pypi.org/project/scipy/>`__

for the examples, plotting requires

- `matplotlib <https://pypi.org/project/matplotlib/>`__


Python
~~~~~~

``pyoscode`` can be installed via pip (*not available yet*)

.. code:: bash

   pip install pyoscode

or via the setup.py

.. code:: bash

   git clone https://github.com/fruzsinaagocs/oscode
   cd oscode
   python setup.py install --user

You can then import ``pyoscode`` from anywhere. Omit the ``--user`` option if
you wish to install globally or in a virtual environment. If you have any
difficulties, check out the FAQs_ section below. 

You can check that things are working by running

.. code:: bash

   pytest tests/

C++
~~~

``oscode`` is a header-only C++ package, it requires no installation.

.. code:: bash

   git clone https://github.com/fruzsinaagocs/oscode

and then include the relevant header files in your C++ code:

.. code:: c

    #include "solver.hpp"
    #include "system.hpp"


Quick start
-----------

Try the following quick examples. These and more are available in the `examples
<https://github.com/fruzsinaagocs/oscode/pyoscode/examples/>`__.

Python
~~~~~~


.. code:: python

    # "airy.py" 
    import pyoscode
    import numpy
    from scipy.special import airy
    from matplotlib import pyplot as plt
    
    # Define the frequency and friction term over the range of integration
    ts = numpy.linspace(1,1000,5000)
    ws = numpy.sqrt(ts)
    gs = numpy.zeros_like(ws)
    # Define the range of integration and the initial conditions
    ti = 1.0
    tf = 1000.
    x0 = airy(-ti)[0] + 1j*airy(-ti)[2]
    dx0 = -airy(-ti)[1] - 1j*airy(-ti)[3]
    # Solve the system
    sol = pyoscode.solve(ts, ws, gs, ti, tf, x0, dx0)
    t = numpy.asarray(sol['t'])
    x = numpy.asarray(sol['sol'])
    types = numpy.asarray(sol['types'])
    # Plot the solution
    ana_t = numpy.linspace(1.0,35.0,1000)
    plt.plot(ana_t,[airy(-T)[0] for T in ana_t],label='true solution')
    plt.plot(t[types==0],x[types==0],'.',color='red',label='RK steps')
    plt.plot(t[types==1],x[types==1],'.',color='green',label='WKB steps')
    plt.legend()
    plt.xlim((1.0,35.0))
    ply.ylim((-1.0,1.0))
    plt.xlabel('t')
    plt.ylabel('Ai(-t)')
    plt.savefig('airy-example.png')
    
The above code, stored in ``airy.py``, produces the plot:

.. image::
   https://github.com/fruzsinaagocs/oscode/raw/master/pyoscode/images/airy-example.png
   :width: 800

``cosmology.ipynb`` is a jupyter notebook that demonstrates how ``pyoscode`` can
be used to quickly generate *primordial power spectra*, like these:

.. image::
    https://github.com/fruzsinaagocs/oscode/raw/master/pyoscode/images/spectra.gif
    :width: 800


C++
~~~

Below is an example where the frequency and friction terms are explicit
functions of time, and are defined as functions. The code is found in
``burst.cpp``, the results are plotted with ``plot_burst.py``.

.. code:: c

    // "burst.cpp"
    #include "solver.hpp"
    #include <cmath>
    #include <fstream>
    #include <string>
    #include <stdlib.h>
    
    double n = 40.0;
    
    // Define the gamma term
    std::complex<double> g(double t){
        return 0.0;
    };
    
    // Define the frequency
    std::complex<double> w(double t){
        return std::pow(n*n - 1.0,0.5)/(1.0 + t*t);
    };
    
    // Initial conditions x, dx
    std::complex<double> xburst(double t){
        return 100*std::pow(1.0 + t*t,
        0.5)/n*std::complex<double>(std::cos(n*std::atan(t)),std::sin(n*std::atan(t))); 
    };
    
    std::complex<double> dxburst(double t){
        return 100/std::pow(1.0 + t*t,
        0.5)/n*(std::complex<double>(t,n)*std::cos(n*std::atan(t)) +
        std::complex<double>(-n,t)*std::sin(n*std::atan(t))); 
    };
    
    int main(){
    
        std::ofstream f;
        std::string output = "output.txt";
        std::complex<double> x0, dx0;
        double ti, tf;
        // Create differential equation 'system'
        de_system sys(&w, &g);
        // Define integration range
        ti = -2*n;
        tf = 2*n;
        // Define initial conditions
        x0 = xburst(ti); 
        dx0 = dxburst(ti); 
        // Solve the equation
        Solution solution(sys, x0, dx0, ti, tf); 
        solution.solve();
        // The solution is stored in lists, copy the solution
        std::list<std::complex<double>> xs = solution.sol;
        std::list<double> ts = solution.times;
        std::list<bool> types = solution.wkbs;
        int steps = solution.ssteps;
        // Write result in file
        f.open(output);
        auto it_t = ts.begin();
        auto it_x = xs.begin();
        auto it_ty = types.begin();
        for(int i=0; i<steps; i++){
            f << *it_t << ", " << std::real(*it_x) << ", " << *it_ty << std::endl; 
            ++it_t;
            ++it_x;
            ++it_ty;
        };
        f.close();
    };

To compile and run:

.. code:: bash

    g++ -g -Wall -std=c++11 -c -o burst.o burst.cpp
    g++ -g -Wall -std=c++11 -o burst burst.o
    ./burst

Plotting the results with Python yields

.. image::
   https://github.com/fruzsinaagocs/oscode/raw/master/pyoscode/images/burst-example.png
   :width: 800


Documentation
-------------

Documentation is hosted at `readthedocs <https://oscode.readthedocs.io>`__.

To build your own local copy of the documentation you'll need to install `sphinx
<https://pypi.org/project/Sphinx/>`__. You can then run:

.. code:: bash

   cd pyoscode/docs
   make html

Citation
--------

If the works below are **in prep.**, please email the authors at fa325@cam.ac.uk
for a copy.

If you use ``oscode`` to solve equations for a publication, please cite
as: ::

   Agocs, F., Handley, W., Lasenby, A., and Hobson, M., (2019). An efficient method for solving highly oscillatory
   ordinary differential equations with applications to physical systems. arXiv
   e-prints, arXiv:1906.01421 (2019) [physics.comp-ph].

or using the BibTeX:

.. code:: bibtex

	@ARTICLE{2019arXiv190601421A,
	       author = {{Agocs}, F.~J. and {Handley}, W.~J. and {Lasenby}, A.~N. and
	         {Hobson}, M.~P.},
	        title = "{An efficient method for solving highly oscillatory ordinary differential equations with applications to physical systems}",
	      journal = {arXiv e-prints},
	     keywords = {Physics - Computational Physics, Astrophysics - Instrumentation and Methods for Astrophysics, Mathematics - Numerical Analysis},
	         year = "2019",
	        month = "May",
	          eid = {arXiv:1906.01421},
	        pages = {arXiv:1906.01421},
	archivePrefix = {arXiv},
	       eprint = {1906.01421},
	 primaryClass = {physics.comp-ph},
	       adsurl = {https://ui.adsabs.harvard.edu/abs/2019arXiv190601421A},
	      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
	}


Contributing
------------

Any comments and improvements to this project are welcome. You can contribute
by:

- Opening and `issue <https://www.github.com/fruzsinaagocs/oscode/issues/>`__ to report bugs and propose new features.
- Making a pull request.

FAQs
----

Installation
~~~~~~~~~~~~

1. Eigen import errors:
    .. code:: bash

       pyoscode/_pyoscode.hpp:6:10: fatal error: Eigen/Dense: No such file or directory
        #include <Eigen/Dense>
                  ^~~~~~~~~~~~~

    Try explicitly including the location of your Eigen library via the
    ``CPLUS_INCLUDE_PATH`` environment variable, for example:

    .. code:: bash

       CPLUS_INCLUDE_PATH=/usr/include/eigen3 python setup.py install --user
       # or 
       CPLUS_INCLUDE_PATH=/usr/include/eigen3 pip install pyoscode

    where  ``/usr/include/eigen3`` should be replaced with your system-specific
    eigen location.


Changelog
---------

- 0.1.0:
    - Memory leaks at python interface fixed
    - C++ documentation added 
