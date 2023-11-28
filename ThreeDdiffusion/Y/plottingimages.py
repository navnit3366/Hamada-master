import argparse
import math
import scipy
import numpy as np
from scipy import special, optimize
import matplotlib.pyplot as plt
import boutdata
from boututils.datafile import DataFile
from boutdata.collect import collect
from boututils.showdata import showdata
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata
N = collect("N")      #Collecting the density from the dmp file
q = 2.1234556436832658946825               #Setting the q-value. The q-value must be the same as in the grid file.
#showdata(N[:,:,2,:]) #Showing the data for N in the original grid.

#Grid values. They should all be the same as in the grid file
nz = 32 #Number of z points
ny = 64 #Number of y points
nx = 36 #Number of x points
MXG = 1 #Number of guardcells in the x-direction
Ly = 1.0#Length of y
Lx = 1.0#Length of x

# Set the internal points
internal_x_points = nx - 2*MXG
internal_y_points = ny

# Set the internal line segments (one less than points)
internal_x_segments = internal_x_points - 1
internal_y_segments = internal_y_points - 1

dx = Lx / (internal_x_segments + 1.0) #x stepsize
dy = Ly / (internal_y_segments + 1.0) #y stepsize

Pi = math.pi
r1 = 0.8 #The small radius. Approximated to a constant.
R = 2    #The large radius.
epsilon = r1/R #Ratio

# Finding the eta value for each theta using the Newton method
def func(eta, Theta):
    return -Theta+eta/(2*Pi)+np.sin(eta)*epsilon/(2*Pi)

etaarr = np.ndarray([nx,ny])
Thetaarr  = np.ndarray([nx, ny])
for grid_point_ntheta in range(ny):
    Theta = grid_point_ntheta*dy - 0.5 #The y-direction values (Theta=y)
    etause = scipy.optimize.newton(func, x0=Pi, args = (Theta,)) #Finding eta using Newtons method
    Thetaarr[:,grid_point_ntheta] = grid_point_ntheta*dy - 0.5
    etaarr[:,grid_point_ntheta] = etause

r = np.linspace(0.8-0.2,0.8+0.2,nx) #The r-coordinate. (Approximated to a constant earlier)
Thetamesh = Thetaarr[1,:] #Creating a vector from the matrix
s = np.linspace(0,1,nz) #The grid for the used 3rd-axis
etamesh = etaarr[1,:]

time = np.linspace(0,0.5*1e-2,101)

#S,ETA = np.meshgrid(s,etamesh)

Y,X,SMESH = np.meshgrid(Thetamesh,r,s)     #Creating a meshgrid using the grid used for simulations
YETA,XETA,SETA = np.meshgrid(etamesh,r,s)  #Creating a meshgrid with an eta-value corresponding to a specific theta value

#PHI = 2*Pi*(-S-(q/(2*Pi))*(-2*np.arctan(np.sqrt((1-epsilon)/(1+epsilon))*np.tan(ETA/2))))

#Getting Torusoidal coordinates from the used grid for reshaping the results to Torusoidal and cartesian coordinates
phi = 2*Pi*(-SMESH-(q/(2*Pi))*(-2*np.arctan(np.sqrt((1-epsilon)/(1+epsilon))*np.tan(YETA/2))))
rplot = X
etaplot = YETA

for time in range(100):
    

#Making a slice for a fixed phi
    search_phi = np.where(((phi[0,:,:]<-8.9)&(phi[0,:,:]>-9.1))|((phi[0,:,:]< (-8.9 + 2*Pi))&(phi[0,:,:]> (-9.1 + 2*Pi))))
    N_fixed = N[time,:,search_phi[0],search_phi[1]]
    N_fixtrans = N_fixed.T
#points = (rplot,etaplot,phi)
#values = np.array(N[1,:,:,:])
#output = np.array((rplot[:,:,1],etaplot[:,:,1],-6))
#grid_new = griddata(points, values, output, method='linear')


#Checking that it does, indeed become a torus
#X_Tcheck = (R + X*np.cos(2*Pi*SMESH))*np.sin(2*Pi*Y)#
#Y_Tcheck = (R + X*np.cos(2*Pi*SMESH))*np.cos(2*Pi*Y)#
#Z_Tcheck = X*np.sin(2*Pi*SMESH)#

#From torusoidal to cartesian coordinates
    X_RT = (R + X*np.cos(etaplot))*np.sin(phi)
    Y_RT = (R + X*np.cos(etaplot))*np.cos(phi)
    Z_RT = X*np.sin(etaplot)

#The values with a fixed phi
    Const_Phi = -4
    X_const = (0.9 + X)*np.sin(etaplot)
    Y_const = (0.9 + X)*np.cos(etaplot)

#Making slice for fixed phi v. 2.0
#points = X_RT[:,:,:]
#values = N[1,:,:,:]
#output = X_const[:,:,:]
#grid_new = griddata(points, values, output, method='linear')

#Plotting the bad shit
    #levels = np.linspace(1, 1.25, 180)
    #myplot24 = plt.contourf(X_const[:,:,3],Y_const[:,:,3],N_fixtrans[:,0:64],levels=levels)
    #cbar = plt.colorbar(mappable=myplot24)
    #name = '/home/aske/test/ThreeDdiffusion/Y/Images/pardiffusionphi9Zdist' + str(time) + '.png'
    #plt.savefig(name)
    #plt.clf()

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_zlim3d(-4,4)
ax.scatter(X_RT,Y_RT,Z_RT,marker = '.')