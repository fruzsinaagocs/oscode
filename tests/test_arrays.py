import numpy as np
import pyoscode
import pytest


def test_dense_output_sorted():
    # Dense output not sorted
    ts = np.linspace(1,50,2000)
    ws = np.linspace(1,50,2000)
    t_eval = np.linspace(1,50,100)
    t_eval = np.flip(t_eval)
    gs = np.zeros_like(ws)
    ti = 1.0
    tf = 50.0
    x0 = 1.0
    dx0 = 0.0
    with pytest.raises(TypeError):
        pyoscode.solve(ts,ws,gs,ti,tf,x0,dx0,t_eval=t_eval)
    return None;

def test_dense_output_range():
    # Dense output outside integration range
    ts = np.linspace(1,50,2000)
    ws = np.linspace(1,50,2000)
    gs = np.zeros_like(ws)
    t_eval = np.linspace(0,60,1000)
    ti = 1.0
    tf = 50.0
    x0 = 1.0
    dx0 = 0.0
    with pytest.raises(TypeError):
        pyoscode.solve(ts,ws,gs,ti,tf,x0,dx0,t_eval=t_eval)
    return None;



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


