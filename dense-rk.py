import numpy as np

# Special point for 6-stage, 5th order RK formula based on GLO nodes for which a 4th order solution
# is available for free
s = 0.5866586817
#s = 0.5
# Specific Butcher tableau b-coefficients for the RK formula to interpolate

# RK45
#b1 = 35./384 
#b2 = 0
#b3 = 500./1113
#b4 = 125./192
#b5 = -2187./6784
#b6 = 11./84
#b7 = 0.

# RK GLO
b1 = 0.112755722735172
b2 = 0.
b3 = 0.506557973265535
b4 = 0.0483004037699511
b5 = 0.378474956297846
b6 = -0.0460890560685063
b7 = 0.

# Butcher tableau b-coefficients for 4th order solution at sigma=s

# RK GLO
b1s = 0.2089555395
b2s = 0.
b3s = 0.7699501023
b4s = 0.009438629906
b5s = -0.003746982422
b6s = 0.01540271068
b7s = 0.

# RK45
#b1s = 6025192743./30085553152
#b2s = 0.
#b3s = 51252292925./65400821598
#b4s = -2691868925./45128329728
#b5s = 187940372067./1594534317056
#b6s =-1776094331./19743644256
#b7s = 11237099./235043384

# 
M = np.array([[1,0,0,0],[1,1,1,1],[1,2,3,4],[s,s**2,s**3,s**4]])
S = np.array([[1,0,0,0,0,0,0],[b1,b2,b3,b4,b5,b6,b7],[0,0,0,0,0,0,1],[s*b1s,s*b2s,s*b3s,s*b4s,s*b5s,s*b6s,s*b7s]])
M.shape, S.shape

P = (S.T).dot(np.linalg.inv(M.T))
P

