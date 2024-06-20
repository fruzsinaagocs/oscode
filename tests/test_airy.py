import numpy as np
import pyoscode
from scipy.special import airy
import pytest

def test_no_integration():
    # Tests the case when ti=tf
    ts = np.linspace(1,100,5000)
    ws = np.sqrt(ts)
    gs = np.zeros_like(ws)
    ti = 1.0
    tf = ti
    x0 = airy(-ti)[0] + 1j*airy(-ti)[2]
    dx0 = -airy(-ti)[1] - 1j*airy(-ti)[3]
    sol = pyoscode.solve(ts, ws, gs, ti, tf, x0, dx0, even_grid=True)
    t = np.asarray(sol['t'])
    x = np.asarray(sol['sol'])
    dx = np.asarray(sol['dsol'])
    assert (x.shape[0] == 1 and dx.shape[0] == 1 and x[0] == x0 and dx[0] == dx0)


def test_airy_forward_even():
    # Integration forward on even grid
    time_space = np.linspace(1,100,5000)
    weights = np.sqrt(time_space)
    gradients = np.zeros_like(weights)
    initial_timepoint = 1.0
    final_timepoint = 100.0
    initial_value = airy(-initial_timepoint)[0] + 1j*airy(-initial_timepoint)[2]
    initial_deriv_value = -airy(-initial_timepoint)[1] - 1j*airy(-initial_timepoint)[3]
    t_eval = np.linspace(initial_timepoint,final_timepoint,5000)
    sol = pyoscode.solve(time_space, weights, gradients, initial_timepoint, final_timepoint, initial_value, initial_deriv_value, t_eval=t_eval, even_grid=True)
    t = np.asarray(sol['t'])
    x = np.asarray(sol['sol'])
    dense = np.asarray(sol['x_eval'])
    dense_d = np.asarray(sol['dx_eval'])
    ana_t = np.linspace(initial_timepoint,final_timepoint,5000)
    ana_x = np.asarray([airy(-T)[0]+1j*airy(-T)[2] for T in t])
    dense_ana_x = np.asarray([airy(-T)[0]+1j*airy(-T)[2] for T in ana_t])
    dense_ana_dx = np.asarray([-airy(-T)[1]-1j*airy(-T)[3] for T in ana_t])
    dense_error = np.abs((dense-dense_ana_x)/dense_ana_x)/0.01 > 1.0
    dense_error_d = np.abs((dense_d-dense_ana_dx)/dense_ana_dx)/0.01 > 1.0
    error =  np.abs((x-ana_x)/ana_x)/0.01 > 1.0
    assert (np.any(dense_error) == False and np.any(dense_error_d) == False and np.any(error) == False )


def test_airy_forward_uneven():
    # Integration forward on uneven grid
    ts = np.logspace(0,2,5000)
    ws = np.sqrt(ts)
    gs = np.zeros_like(ws)
    ti = 1.0
    tf = 100.0
    x0 = airy(-ti)[0] + 1j*airy(-ti)[2]
    dx0 = -airy(-ti)[1] - 1j*airy(-ti)[3]
    t_eval = np.linspace(ti,tf,5000)
    sol = pyoscode.solve(ts, ws, gs, ti, tf, x0, dx0, t_eval=t_eval)
    t = np.asarray(sol['t'])
    x = np.asarray(sol['sol'])
    dense = np.asarray(sol['x_eval'])
    dense_d = np.asarray(sol['dx_eval'])
    types = np.asarray(sol['types'])
    ana_t = np.linspace(ti,tf,5000)
    ana_x = np.asarray([airy(-T)[0]+1j*airy(-T)[2] for T in t])
    dense_ana_x = np.asarray([airy(-T)[0]+1j*airy(-T)[2] for T in ana_t])
    dense_ana_dx = np.asarray([-airy(-T)[1]-1j*airy(-T)[3] for T in ana_t])
    dense_error = np.abs((dense-dense_ana_x)/dense_ana_x)/0.01 > 1.0
    dense_error_d = np.abs((dense_d-dense_ana_dx)/dense_ana_dx)/0.01 > 1.0
    error =  np.abs((x-ana_x)/ana_x)/0.01 > 1.0
    assert (np.any(dense_error) == False and np.any(dense_error_d) == False and np.any(error) == False )


def test_airy_backward_even():
    # Integration backwards on even grid
    ts = np.linspace(1,100,5000)
    ws = np.sqrt(ts)
    gs = np.zeros_like(ws)
    ti = 1.0
    tf = 100.0
    t_eval = np.linspace(ti,tf,5000)
    x0 = airy(-tf)[0] + 1j*airy(-tf)[2]
    dx0 = -airy(-tf)[1] - 1j*airy(-tf)[3]
    sol = pyoscode.solve(ts, ws, gs, tf, ti, x0, dx0, t_eval=t_eval, even_grid=True)
    t = np.asarray(sol['t'])
    x = np.asarray(sol['sol'])
    dense = np.asarray(sol['x_eval'])
    dense_d = np.asarray(sol['dx_eval'])
    types = np.asarray(sol['types'])
    ana_t = np.linspace(ti,tf,5000)
    ana_x = np.asarray([airy(-T)[0]+1j*airy(-T)[2] for T in t])
    dense_ana_x = np.asarray([airy(-T)[0]+1j*airy(-T)[2] for T in ana_t])
    dense_ana_dx = np.asarray([-airy(-T)[1]-1j*airy(-T)[3] for T in ana_t])
    dense_error = np.abs((dense-dense_ana_x)/dense_ana_x)/0.01 > 1.0
    dense_error_d = np.abs((dense_d-dense_ana_dx)/dense_ana_dx)/0.01 > 1.0
    error =  np.abs((x-ana_x)/ana_x)/0.01 > 1.0
    assert (np.any(dense_error) == False and np.any(dense_error_d) == False and np.any(error) == False )



def test_airy_backward_uneven():
    # Integration backwards on uneven grid
    ts = np.logspace(0,2,5000)
    ws = np.sqrt(ts)
    gs = np.zeros_like(ws)
    ti = 1.0
    tf = 100.0
    t_eval = np.linspace(ti,tf,5000)
    x0 = airy(-tf)[0] + 1j*airy(-tf)[2]
    dx0 = -airy(-tf)[1] - 1j*airy(-tf)[3]
    sol = pyoscode.solve(ts, ws, gs, tf, ti, x0, dx0, t_eval=t_eval)
    t = np.asarray(sol['t'])
    x = np.asarray(sol['sol'])
    dense = np.asarray(sol['x_eval'])
    dense_d = np.asarray(sol['dx_eval'])   
    types = np.asarray(sol['types'])
    ana_t = t_eval
    ana_x = np.asarray([airy(-T)[0]+1j*airy(-T)[2] for T in t])
    dense_ana_x = np.asarray([airy(-T)[0]+1j*airy(-T)[2] for T in ana_t])
    dense_ana_dx = np.asarray([-airy(-T)[1]-1j*airy(-T)[3] for T in ana_t])
    dense_error = np.abs((dense-dense_ana_x)/dense_ana_x)/0.01 > 1.0
    dense_error_d = np.abs((dense_d-dense_ana_dx)/dense_ana_dx)/0.01 > 1.0
    error =  np.abs((x-ana_x)/ana_x)/0.01 > 1.0
    assert (np.any(dense_error) == False and np.any(dense_error_d) == False and np.any(error) == False )

def test_airy_wgfunctions():
    # w, g are supplied as functions instead of arrays
    w = lambda t: np.sqrt(t)
    g = lambda t: 0.0
    ti = 1.0
    tf = 100.0
    t_eval = np.linspace(ti,tf,5000)
    x0 = airy(-ti)[0] + 1j*airy(-ti)[2]
    dx0 = -airy(-ti)[1] - 1j*airy(-ti)[3]
    sol = pyoscode.solve_fn(w, g, ti, tf, x0, dx0, t_eval=t_eval)
    t = np.asarray(sol['t'])
    x = np.asarray(sol['sol'])
    dense = np.asarray(sol['x_eval'])
    dense_d = np.asarray(sol['dx_eval'])   
    types = np.asarray(sol['types'])
    ana_t = t_eval
    ana_x = np.asarray([airy(-T)[0]+1j*airy(-T)[2] for T in t])
    dense_ana_x = np.asarray([airy(-T)[0]+1j*airy(-T)[2] for T in ana_t])
    dense_ana_dx = np.asarray([-airy(-T)[1]-1j*airy(-T)[3] for T in ana_t])
    dense_error = np.abs((dense-dense_ana_x)/dense_ana_x)/0.01 > 1.0
    dense_error_d = np.abs((dense_d-dense_ana_dx)/dense_ana_dx)/0.01 > 1.0
    error =  np.abs((x-ana_x)/ana_x)/0.01 > 1.0
    assert (np.any(dense_error) == False and np.any(dense_error_d) == False and np.any(error) == False )

def test_airy_wgfunctions_back():
    # w, g are supplied as functions instead of arrays
    # Integrating backwards
    w = lambda t: np.sqrt(t)
    g = lambda t: 0.0
    ti = 1.0
    tf = 100.0
    t_eval = np.linspace(ti,tf,5000)
    x0 = airy(-tf)[0] + 1j*airy(-tf)[2]
    dx0 = -airy(-tf)[1] - 1j*airy(-tf)[3]
    sol = pyoscode.solve_fn(w, g, tf, ti, x0, dx0, t_eval=t_eval)
    t = np.asarray(sol['t'])
    x = np.asarray(sol['sol'])
    dense = np.asarray(sol['x_eval'])
    dense_d = np.asarray(sol['dx_eval'])   
    types = np.asarray(sol['types'])
    ana_t = t_eval
    ana_x = np.asarray([airy(-T)[0]+1j*airy(-T)[2] for T in t])
    dense_ana_x = np.asarray([airy(-T)[0]+1j*airy(-T)[2] for T in ana_t])
    dense_ana_dx = np.asarray([-airy(-T)[1]-1j*airy(-T)[3] for T in ana_t])
    dense_error = np.abs((dense-dense_ana_x)/dense_ana_x)/0.01 > 1.0
    dense_error_d = np.abs((dense_d-dense_ana_dx)/dense_ana_dx)/0.01 > 1.0
    error =  np.abs((x-ana_x)/ana_x)/0.01 > 1.0
    assert (np.any(dense_error) == False and np.any(dense_error_d) == False and np.any(error) == False )


