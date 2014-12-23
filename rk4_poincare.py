from matplotlib.pyplot import *
from scipy import *
from numpy import *

# a simple Runge-Kutta integrator for multiple dependent variables and one independent variable

def rungekutta4(yprime, time, y0):
    # yprime is a list of functions, y0 is a list of initial values of y
    # time is a list of t-values at which solutions are computed
    #
    # Dependency: numpy

    N = len(time)

    y = array([thing*ones(N) for thing in y0]).T
    
    for ii in xrange(N-1):
        dt = time[ii+1] - time[ii]
        k1 = dt*yprime(y[ii], time[ii])
        k2 = dt*yprime(y[ii] + 0.5*k1, time[ii] + 0.5*dt)
        k3 = dt*yprime(y[ii] + 0.5*k2, time[ii] + 0.5*dt)
        k4 = dt*yprime(y[ii] + k3, time[ii+1])
        y[ii+1] = y[ii] + (k1 + 2.0*(k2 + k3) + k4)/6.0

    return y

# Miscellaneous functions

def total_energy(valpair):
    (x, y, px, py) = tuple(valpair)
    return .5*(px**2 + py**2) + .5*(x**2 + y**2 + 2.*x**2.*y - (2.0/3)*y**3)
    
def pqdot(valpair, tval):
    # input: [x, y, px, py], t
    # takes a pair of x and y values and returns \dot{p} according to the Henon-Heiles Hamiltonian
    (x, y, px, py) = tuple(valpair)
    return array([px, py, -x-2*x*y, -y-x*x+y*y]).T

def findcrossings(data):
    # returns indices in 1D data set where the data crossed zero. Useful for generating Poincare map at 0
    prb = list()
    for ii in xrange(len(data)-1):
        if (data[ii] > 0)&(data[ii+1] < 0):
            prb.append(ii)
        if (data[ii] < 0)& (data[ii+1] > 0):
            prb.append(ii)
    return array(prb)

# # for testing:
#
# def test2(valpair, tval):
#     xval = valpair[0]
#     yval = valpair[1]
#     return array([-7*yval, -yval+xval]).T
# kk = rungekutta4(test2, linspace(1,10,100), array([4.0, 2.0]))
# plot(kk)