# Test out the runge kutta library
# William Gilpin 2014

from matplotlib import pyplot
from scipy import *
from numpy import *

from rk4_poincare import *

t = linspace(0, 100.0, 1000)
print ("step size is " + str(t[1]-t[0]))

# Four representative initial conditions for E=1/12 on the Henon-Heiles system
init_cons = [array([0.0, 0.5, 0.0, 0.0]),\
             array([0.408248, 0.0, 0.0, 0.0]),\
             array([0.0, 0.0, 0.408248, 0.0]),\
             array([0.0, 0.0, 0.0, 0.408248])\
             ]

outs = list()
for con in init_cons:
    outs.append( rungekutta4(pqdot, t, con) )


# plot the results
fig1 = figure(1)
for ii in xrange(4):
    subplot(2, 2, ii+1)
    plot(outs[ii][:,1],outs[ii][:,3])
    ylabel("py")
    xlabel("y")
    title("Full trajectory projected onto the plane")

fig1.suptitle('Full trajectories E = 1/12', fontsize=20)



# Plot Poincare sections at x=0
fig2 = figure(2)
for ii in xrange(4):
    subplot(2, 2, ii+1)
    xcrossings = findcrossings(outs[ii][:,0])
    yints = [.5*(outs[ii][cross, 1] + outs[ii][cross+1, 1]) for cross in xcrossings]
    pyints = [.5*(outs[ii][cross, 3] + outs[ii][cross+1, 3]) for cross in xcrossings]
    plot(yints, pyints,'.')
    ylabel("py")
    xlabel("y")
    title("Poincare section x = 0")

fig2.suptitle('Poincare Sections E = 1/12', fontsize=20)


show()


