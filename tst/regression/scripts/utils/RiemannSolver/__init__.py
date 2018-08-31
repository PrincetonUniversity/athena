#from numpy import finfo


brent_opt = dict(full_output=True, disp=True)
ode_opt = dict(mxstep=5000000, rtol=1e-15, atol=1e-13)#, rtol=2 * finfo(float).eps, atol=2 * finfo(float).eps)