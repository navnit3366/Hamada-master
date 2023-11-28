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
import scipy.interpolate as spint

RGI = spint.RegularGridInterpolator

N = collect("N")      #Collecting the density from the dmp file
#q = 1.2               #Setting the q-value. The q-value must be the same as in the grid file.
#showdata(N[:,:,2,:]) #Showing the data for N in the original grid.

#Grid values. They should all be the same as in the grid file
nz = 2 #Number of z points
ny = 128 #Number of y points
nx = 20 #Number of x points
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
#r1 = 0.8 #The small radius. Approximated to a constant.
#R = 2    #The large radius.
#const_phi = -3.5

name = "fieldaligned"

#Reading the values from the grid-generating file
with DataFile(name + ".grd.nc", write=False, create=False) as grid:
    q = grid.read("q")
    r1 = grid.read("r")
    R = grid.read("R")
    s = grid.read("s")
    Eta = grid.read("Eta")
    Thetause = grid.read("Theta")
    
print(q)
epsilon = r1/R #Ratio
z_length = s.size
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
Thetamesh = Thetause[1,:] #Creating a vector from the matrix
#s = np.linspace(0,0.96874999999998668,z_length) #The grid for the used 3rd-axis
etamesh = Eta[1,:]

#s_fixedphi = -(1/(2*Pi))*const_phi + (q/Pi)*np.arctan(np.sqrt((1-epsilon)/(1+epsilon))*np.tan(etaarr[1,:]/2))

time = 0

#S,ETA = np.meshgrid(s,etamesh)

Y,X,Z = np.meshgrid(Thetamesh,r,s)     #Creating a meshgrid using the grid used for simulations
YETA,XETA,ZETA = np.meshgrid(etamesh,r,s)  #Creating a meshgrid with an eta-value corresponding to a specific theta value
#Y_CONSTPHI,X_CONSTPHI,Z_CONSTPHI = np.meshgrid(Thetamesh,r,s_fixedphi)

phi = 2*Pi*(-Z-(q/(2*Pi))*(-2*np.arctan(np.sqrt((1-epsilon)/(1+epsilon))*np.tan(YETA/2))))
rplot = X
etaplot = YETA

X_RT = (R + X*np.cos(etaplot))*np.sin(phi)
Y_RT = (R + X*np.cos(etaplot))*np.cos(phi)
Z_RT = X*np.sin(etaplot)

points = np.ndarray([nx*ny*z_length,3])
values = np.ndarray([nx*ny*z_length,1])

phi_const = -6

h=0
for i in range(nx):
    for j in range(ny):
        for k in range(z_length):
            points[h,0] = X_RT[i,j,k]
            points[h,1] = Y_RT[i,j,k]
            points[h,2] = Z_RT[i,j,k]
            values[h,0] = N[time,i,j,k]
            h+=1

x_grid,y_grid,z_grid = np.mgrid[X_RT.min():X_RT.max():100j , Y_RT.min():Y_RT.max():100j,Z_RT.min():Z_RT.max():100j]

grid_z1 = griddata(points, values, (x_grid, y_grid, z_grid), method='linear')
#PHI = 2*Pi*(-S-(q/(2*Pi))*(-2*np.arctan(np.sqrt((1-epsilon)/(1+epsilon))*np.tan(ETA/2))))

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_zlim3d(-4,4)
ax.scatter(X_RT,Y_RT,Z_RT,marker = '.')

#Getting Torusoidal coordinates from the used grid for reshaping the results to Torusoidal and cartesian coordinates
