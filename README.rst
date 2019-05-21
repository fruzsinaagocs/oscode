========================================================================
oscode: Oscillatory ordinary differential equation solver
========================================================================

:oscode: oscillatory ordinary differential equation solver
:Author: Fruzsina Agocs, Will Handley, Mike Hobson, and Anthony Lasenby
:Version: 1.0
:Homepage: https://github.com/fruzsinaagocs/oscode
:Documentation: https://oscode.readthedocs.io

``oscode`` is a C++ tool with a Python interface that solves **osc**\illatory
**o**\rdinary **d**\ifferential **e**\quations efficiently. It is designed to
deal with equations of the form


.. image:: 
    https://github.com/fruzsinaagocs/oscode/raw/master/pyoscode/images/oscillator.png

where |gamma| and |omega| can be given as

.. |gamma| image:: https://github.com/fruzsinaagocs/oscode/raw/master/pyoscode/images/gamma.png

.. |omega| image:: https://github.com/fruzsinaagocs/oscode/raw/master/pyoscode/images/omega.png

- *In C++*, explicit functions or sequence containers (Eigen::Vectors, arrays,
  std::vectors, lists),
- *In Python*, numpy.arrays.

``oscode`` makes use of an analytic approximation of :math:`x(t)` embedded in a
stepping procedure to skip over long regions of oscillations, giving a reduction
in computing time. The approximation is valid when the frequency
:math:`\omega(t)` changes slowly relative to the timescales of integration, it
is therefore worth applying when this condition holds for at least some part of
the integration range. 

The numerical method used by ``oscode`` is described in detail in `this work.
<https://>`__


Installation
------------

Python
~~~~~~

``oscode`` can be installed via pip

.. code:: bash

   pip install pyoscode

or via the setup.py

.. code:: bash

   git clone https://github.com/fruzsinaagocs/oscode
   cd oscode
   python setup.py install --user

or 

.. code:: bash

    git clone https://github.com/fruzsinaagocs/oscode
    cd oscode
    python setup.py build_ext --inplace


C++
~~~

``oscode`` is a header-only C++ package, it requires no installation.

.. code:: bash

   git clone https://github.com/fruzsinaagocs/oscode

and then ``#include`` them in your C++ code. 


Dependencies
~~~~~~~~~~~~

Basic requirements: 

Python
^^^^^^

- Python 2.7 or 3.5+
- `numpy <https://pypi.org/project/numpy/>`
- `scipy <https://pypi.org/project/scipy/>`

C++ 
^^^

- `Eigen <http://eigen.tuxfamily.org/index.php?title=Main_Page>`

Quick start
-----------

Try the following example(s):

Python
~~~~~~

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

.. image::
   https://github.com/fruzsinaagocs/oscode/raw/master/pyoscode/images/airy-example.png
   :width: 800

C++
~~~



Documentation
-------------

To build your own local copy of the documentation you'll need to install `sphinx
<https://pypi.org/project/Sphinx/>`__. You can then run:

.. code:: bash

   cd pyoscode/docs
   make html

Citation
--------


If you use ``oscode`` to solve equations for a publication, please cite
as: ::

   Agocs et al., (2019). ...

or using the BibTeX:

.. code:: bibtex

   @article{oscode,
       doi = {},
       url = {},
       year  = {},
       month = {},
       publisher = {},
       volume = {},
       number = {},
       author = {},
       title = {},
       journal = {}
   }

