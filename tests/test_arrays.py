import numpy as np
import pyoscode
from scipy.special import airy
import pytest

def test_different_size():
    # Different size arrays
    ts = np.linspace(1,50,2000)
    ws = np.linspace(1,50,2001)
    gs = np.zeros_like(ws)
    ti = 1.0
    tf = 50.0
    x0 = 1.0
    dx0 = 0.0
    with pytest.raises(TypeError):
        pyoscode.solve(ts,ws,gs,ti,tf,x0,dx0)
    return None;

def test_too_small():
    # Same size arrays but too small
    ts = np.linspace(1,50,1)
    ws = np.linspace(1,50,1)
    gs = np.zeros_like(ws)
    ti = 1.0
    tf = 50.0
    x0 = 1.0
    dx0 = 0.0
    with pytest.raises(TypeError):
        pyoscode.solve(ts,ws,gs,ti,tf,x0,dx0)
    return None;

def test_outside_range_ti():
    # ti outside ts
    ts = np.linspace(1,50,1000)
    ws = np.linspace(1,50,1000)
    gs = np.zeros_like(ws)
    ti = 0.5
    tf = 50.0
    x0 = 1.0
    dx0 = 0.0
    with pytest.raises(TypeError):
        pyoscode.solve(ts,ws,gs,ti,tf,x0,dx0)
    return None;

def test_outside_range_tf():
    # ti outside ts
    ts = np.linspace(1,50,1000)
    ws = np.linspace(1,50,1000)
    gs = np.zeros_like(ws)
    ti = 1.0
    tf = -5.0
    x0 = 1.0
    dx0 = 0.0
    with pytest.raises(TypeError):
        pyoscode.solve(ts,ws,gs,ti,tf,x0,dx0)
    return None;

def test_not_monotonous():
    # ts not strictly monotonous
    ts = np.random.rand(1000)
    ws = ts
    gs = ws
    ti = ts[0]
    tf = ts[-1]
    x0 = 1.0
    dx0 = 0.0
    with pytest.raises(TypeError):
        pyoscode.solve(ts,ws,gs,ti,tf,x0,dx0)
    return None;

def test_airy_forward():
    # Integration forward without dense output, on even grid
    ts = np.linspace(1,100,5000)
    ws = np.sqrt(ts)
    gs = np.zeros_like(ws)
    ti = 1.0
    tf = 100.0
    x0 = airy(-ti)[0] + 1j*airy(-ti)[2]
    print(x0)
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


#def test_airy_forward():
#    # Integration forward without dense output 
#    # Define the frequency and friction term over the range of integration
#    ts = np.logspace(0,2,5000)
#    #ts = np.linspace(1,100,5000)
#    ws = np.sqrt(ts)
#    gs = np.zeros_like(ws)
#    # Define the range of integration and the initial conditions
#    ti = 1.0
#    tf = 100.0
#    x0 = airy(-ti)[0] + 1j*airy(-ti)[2]
#    dx0 = -airy(-ti)[1] - 1j*airy(-ti)[3]
#    t_eval = np.linspace(ti,tf,5000)
#    # Solve the system
#    # Test forward integration
#    #sol = pyoscode.solve(ts, ws, gs, ti, tf, x0, dx0,t_eval=t_eval) # Uneven grid assumed, even grid given
#    #sol = pyoscode.solve(ts, ws, gs, ti, tf, x0, dx0,t_eval=t_eval,even_grid=True) # Even grid expected, even grid given
#    # Test backward integration 
#    x0 = airy(-tf)[0] + 1j*airy(-tf)[2]
#    dx0 = -airy(-tf)[1] - 1j*airy(-tf)[3]
#    #sol = pyoscode.solve(ts, ws, gs, tf, ti, x0, dx0,t_eval=t_eval,even_grid=True) # X
#    sol = pyoscode.solve(ts, ws, gs, tf, ti, x0, dx0,t_eval=t_eval) # X
#    t = np.asarray(sol['t'])
#    x = np.asarray(sol['sol'])
#    dense = np.asarray(sol['x_eval'])
#    types = np.asarray(sol['types'])
#    # Plot the solution
#    ana_t = np.linspace(ti,tf,5000)
#
#def test_airy_forward():
#    # Integration forward without dense output 
#    # Define the frequency and friction term over the range of integration
#    ts = np.logspace(0,2,5000)
#    #ts = np.linspace(1,100,5000)
#    ws = np.sqrt(ts)
#    gs = np.zeros_like(ws)
#    # Define the range of integration and the initial conditions
#    ti = 1.0
#    tf = 100.0
#    x0 = airy(-ti)[0] + 1j*airy(-ti)[2]
#    dx0 = -airy(-ti)[1] - 1j*airy(-ti)[3]
#    t_eval = np.linspace(ti,tf,5000)
#    # Solve the system
#    # Test forward integration
#    #sol = pyoscode.solve(ts, ws, gs, ti, tf, x0, dx0,t_eval=t_eval) # Uneven grid assumed, even grid given
#    #sol = pyoscode.solve(ts, ws, gs, ti, tf, x0, dx0,t_eval=t_eval,even_grid=True) # Even grid expected, even grid given
#    # Test backward integration 
#    x0 = airy(-tf)[0] + 1j*airy(-tf)[2]
#    dx0 = -airy(-tf)[1] - 1j*airy(-tf)[3]
#    #sol = pyoscode.solve(ts, ws, gs, tf, ti, x0, dx0,t_eval=t_eval,even_grid=True) # X
#    sol = pyoscode.solve(ts, ws, gs, tf, ti, x0, dx0,t_eval=t_eval) # X
#    t = np.asarray(sol['t'])
#    x = np.asarray(sol['sol'])
#    dense = np.asarray(sol['x_eval'])
#    types = np.asarray(sol['types'])
#    # Plot the solution
#    ana_t = np.linspace(ti,tf,5000)
#
#def test_airy_forward():
#    # Integration forward without dense output 
#    # Define the frequency and friction term over the range of integration
#    ts = np.logspace(0,2,5000)
#    #ts = np.linspace(1,100,5000)
#    ws = np.sqrt(ts)
#    gs = np.zeros_like(ws)
#    # Define the range of integration and the initial conditions
#    ti = 1.0
#    tf = 100.0
#    x0 = airy(-ti)[0] + 1j*airy(-ti)[2]
#    dx0 = -airy(-ti)[1] - 1j*airy(-ti)[3]
#    t_eval = np.linspace(ti,tf,5000)
#    # Solve the system
#    # Test forward integration
#    #sol = pyoscode.solve(ts, ws, gs, ti, tf, x0, dx0,t_eval=t_eval) # Uneven grid assumed, even grid given
#    #sol = pyoscode.solve(ts, ws, gs, ti, tf, x0, dx0,t_eval=t_eval,even_grid=True) # Even grid expected, even grid given
#    # Test backward integration 
#    x0 = airy(-tf)[0] + 1j*airy(-tf)[2]
#    dx0 = -airy(-tf)[1] - 1j*airy(-tf)[3]
#    #sol = pyoscode.solve(ts, ws, gs, tf, ti, x0, dx0,t_eval=t_eval,even_grid=True) # X
#    sol = pyoscode.solve(ts, ws, gs, tf, ti, x0, dx0,t_eval=t_eval) # X
#    t = np.asarray(sol['t'])
#    x = np.asarray(sol['sol'])
#    dense = np.asarray(sol['x_eval'])
#    types = np.asarray(sol['types'])
#    # Plot the solution
#    ana_t = np.linspace(ti,tf,5000)
#
#def test_airy_forward():
#    # Integration forward without dense output 
#    # Define the frequency and friction term over the range of integration
#    ts = np.logspace(0,2,5000)
#    #ts = np.linspace(1,100,5000)
#    ws = np.sqrt(ts)
#    gs = np.zeros_like(ws)
#    # Define the range of integration and the initial conditions
#    ti = 1.0
#    tf = 100.0
#    x0 = airy(-ti)[0] + 1j*airy(-ti)[2]
#    dx0 = -airy(-ti)[1] - 1j*airy(-ti)[3]
#    t_eval = np.linspace(ti,tf,5000)
#    # Solve the system
#    # Test forward integration
#    #sol = pyoscode.solve(ts, ws, gs, ti, tf, x0, dx0,t_eval=t_eval) # Uneven grid assumed, even grid given
#    #sol = pyoscode.solve(ts, ws, gs, ti, tf, x0, dx0,t_eval=t_eval,even_grid=True) # Even grid expected, even grid given
#    # Test backward integration 
#    x0 = airy(-tf)[0] + 1j*airy(-tf)[2]
#    dx0 = -airy(-tf)[1] - 1j*airy(-tf)[3]
#    #sol = pyoscode.solve(ts, ws, gs, tf, ti, x0, dx0,t_eval=t_eval,even_grid=True) # X
#    sol = pyoscode.solve(ts, ws, gs, tf, ti, x0, dx0,t_eval=t_eval) # X
#    t = np.asarray(sol['t'])
#    x = np.asarray(sol['sol'])
#    dense = np.asarray(sol['x_eval'])
#    types = np.asarray(sol['types'])
#    # Plot the solution
#    ana_t = np.linspace(ti,tf,5000)
#
#def test_airy_forward():
#    # Integration forward without dense output 
#    # Define the frequency and friction term over the range of integration
#    ts = np.logspace(0,2,5000)
#    #ts = np.linspace(1,100,5000)
#    ws = np.sqrt(ts)
#    gs = np.zeros_like(ws)
#    # Define the range of integration and the initial conditions
#    ti = 1.0
#    tf = 100.0
#    x0 = airy(-ti)[0] + 1j*airy(-ti)[2]
#    dx0 = -airy(-ti)[1] - 1j*airy(-ti)[3]
#    t_eval = np.linspace(ti,tf,5000)
#    # Solve the system
#    # Test forward integration
#    #sol = pyoscode.solve(ts, ws, gs, ti, tf, x0, dx0,t_eval=t_eval) # Uneven grid assumed, even grid given
#    #sol = pyoscode.solve(ts, ws, gs, ti, tf, x0, dx0,t_eval=t_eval,even_grid=True) # Even grid expected, even grid given
#    # Test backward integration 
#    x0 = airy(-tf)[0] + 1j*airy(-tf)[2]
#    dx0 = -airy(-tf)[1] - 1j*airy(-tf)[3]
#    #sol = pyoscode.solve(ts, ws, gs, tf, ti, x0, dx0,t_eval=t_eval,even_grid=True) # X
#    sol = pyoscode.solve(ts, ws, gs, tf, ti, x0, dx0,t_eval=t_eval) # X
#    t = np.asarray(sol['t'])
#    x = np.asarray(sol['sol'])
#    dense = np.asarray(sol['x_eval'])
#    types = np.asarray(sol['types'])
#    # Plot the solution
#    ana_t = np.linspace(ti,tf,5000)
#
