import numpy as np
import pyoscode
from scipy.special import airy
import pytest

def test_airy_forward_even():
    # Integration forward on even grid
    ts = np.linspace(1,100,5000)
    ws = np.sqrt(ts)
    gs = np.zeros_like(ws)
    ti = 1.0
    tf = 100.0
    x0 = airy(-ti)[0] + 1j*airy(-ti)[2]
    dx0 = -airy(-ti)[1] - 1j*airy(-ti)[3]
    t_eval = np.linspace(ti,tf,5000)
    sol = pyoscode.solve(ts, ws, gs, ti, tf, x0, dx0,t_eval=t_eval,even_grid=True)
    t = np.asarray(sol['t'])
    x = np.asarray(sol['sol'])
    dense = np.asarray(sol['x_eval'])
    ana_t = np.linspace(ti,tf,5000)
    ana_x = np.asarray([airy(-T)[0]+1j*airy(-T)[2] for T in t])
    dense_ana_x = np.asarray([airy(-T)[0]+1j*airy(-T)[2] for T in ana_t])
    dense_error = np.abs((dense-dense_ana_x)/dense_ana_x)/0.01 > 1.0
    error =  np.abs((x-ana_x)/ana_x)/0.01 > 1.0
    assert (np.any(dense_error) == False and np.any(error) == False )


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
    sol = pyoscode.solve(ts, ws, gs, ti, tf, x0, dx0,t_eval=t_eval)
    t = np.asarray(sol['t'])
    x = np.asarray(sol['sol'])
    dense = np.asarray(sol['x_eval'])
    types = np.asarray(sol['types'])
    ana_t = np.linspace(ti,tf,5000)
    ana_x = np.asarray([airy(-T)[0]+1j*airy(-T)[2] for T in t])
    dense_ana_x = np.asarray([airy(-T)[0]+1j*airy(-T)[2] for T in ana_t])
    dense_error = np.abs((dense-dense_ana_x)/dense_ana_x)/0.01 > 1.0
    error =  np.abs((x-ana_x)/ana_x)/0.01 > 1.0
    assert (np.any(dense_error) == False and np.any(error) == False )


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
    sol = pyoscode.solve(ts, ws, gs, tf, ti, x0, dx0,t_eval=t_eval,even_grid=True)
    t = np.asarray(sol['t'])
    x = np.asarray(sol['sol'])
    dense = np.asarray(sol['x_eval'])
    types = np.asarray(sol['types'])
    ana_t = np.linspace(ti,tf,5000)
    ana_x = np.asarray([airy(-T)[0]+1j*airy(-T)[2] for T in t])
    dense_ana_x = np.asarray([airy(-T)[0]+1j*airy(-T)[2] for T in ana_t])
    dense_error = np.abs((dense-dense_ana_x)/dense_ana_x)/0.01 > 1.0
    error =  np.abs((x-ana_x)/ana_x)/0.01 > 1.0
    assert (np.any(dense_error) == False and np.any(error) == False )


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
    sol = pyoscode.solve(ts, ws, gs, tf, ti, x0, dx0,t_eval=t_eval)
    t = np.asarray(sol['t'])
    x = np.asarray(sol['sol'])
    dense = np.asarray(sol['x_eval'])
    types = np.asarray(sol['types'])
    ana_t = np.linspace(ti,tf,5000)
    ana_x = np.asarray([airy(-T)[0]+1j*airy(-T)[2] for T in t])
    dense_ana_x = np.asarray([airy(-T)[0]+1j*airy(-T)[2] for T in ana_t])
    dense_error = np.abs((dense-dense_ana_x)/dense_ana_x)/0.01 > 1.0
    error =  np.abs((x-ana_x)/ana_x)/0.01 > 1.0
    assert (np.any(dense_error) == False and np.any(error) == False )



