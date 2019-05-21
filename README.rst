=========================================================
oscode: Oscillatory ordinary differential equation solver
=========================================================

:oscode: oscillatory ordinary differential equation solver
:Author: Fruzsina Agocs, Will Handley, Mike Hobson, and Anthony Lasenby
:Version: 1.0
:Homepage: https://github.com/agocsfruzsina/oscode
:Documentation:

Installation
------------

``oscode`` can be installed via the setup.py

.. code:: bash

   git clone https://github.com/agocsfruzsina/oscode
   cd oscode
   python setup.py install --user

or 

.. code:: bash

    git clone https://github.com/agocsfruzsina/oscode
    cd oscode
    python setup.py build_ext --inplace


Dependencies
~~~~~~~~~~~~

Basic requirements: 

- Python 2.7 or 3.5+
- `numpy`
- `scipy`
- `Eigen`

Quick start
-----------

Try the following minimal working example:

.. code:: python

    import pyoscode
    import numpy
    from scipy.special import airy
    from matplotlib import pyplot as plt
    
    # Define the frequency and friction term over the range of integration
    ts = numpy.linspace(1,35,5000)
    ws = numpy.sqrt(ts)
    gs = numpy.zeros_like(ws)
    # Define the range of integration and the initial conditions
    ti = 1.0
    tf = 35.0
    x0 = airy(-ti)[0] + 1j*airy(-ti)[2]
    dx0 = -airy(-ti)[1] - 1j*airy(-ti)[3]
    # Solve the system
    sol = pyoscode.solve(ts, ws, gs, ti, tf, x0, dx0)
    t = numpy.asarray(sol['t'])
    x = numpy.asarray(sol['sol'])
    types = numpy.asarray(sol['types'])
    # Plot the solution
    plt.plot(ts,[airy(-T)[0] for T in ts],label='true solution')
    plt.plot(t[types==0],x[types==0],'.',color='red',label='RK steps')
    plt.plot(t[types==1],x[types==1],'.',color='green',label='WKB steps')
    plt.legend()
    plt.xlabel('t')
    plt.ylabel('Ai(-t)')
    plt.show()

This produces

.. image:: images/airy-example.png
   :width: 800


Documentation
-------------

To build your own local copy of the documentation you'll need to install `sphinx
<https://pypi.org/project/Sphinx/>`__. You can then run:

.. code:: bash

   cd pyoscode/docs
   make html

Citation
--------

