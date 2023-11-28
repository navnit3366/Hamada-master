#!/usr/bin/env python
"""Generates a field-aligned grid and backwards mapping grid"""

import argparse
import math
import scipy
import numpy as np
from scipy import special, optimize
import matplotlib.pyplot as plt
import boutdata
from boututils.datafile import DataFile

def generate_grid(nx = 20, ny = 128, nz=22,\
                  Lx = 1.0, Ly = 1.0,\
                  MXG = 1, name="fieldaligned"):
    """
    Function which generates a gridfile for cylindrical coordinates.
    We will here use the following convention:
    x - The radial coordinate (rho)
    y - The height of the cylinder (z)
    z - The azimuthal coordinate (theta)


    Input parameters
    ----------------
    nx      - nx = internal_x_points + 2*MXG
    ny      - ny = internal_y_points
    nz      - nz = Inner points in MZ (NOT containing 1 artificial point)
    Lx      - Full radial length
    Ly      - Height of cylinder
    MXG     - The number of guard cells in x-direction
    name    - Name of the output file (without extention)

    Output
    ------
    A grid file with .nc extension containing the covariant and the
    contravariant metric tensor elements, together with J
    """


    # Set the internal points
    internal_x_points = nx - 2*MXG
    internal_y_points = ny

    # Set the internal line segments (one less than points)
    internal_x_segments = internal_x_points - 1
    internal_y_segments = internal_y_points - 1

    # Position of the separatrix (-1 is non periodic, >ny is periodic)
    # --------Periodic----------
    ixseps1 = nx
    ixseps2 = nx
    # ------------------------------

    # Calculate the grid size
    # In rho, we put the boundary (which is implemented half between the
    # grid points) in r = 0
    # If we plot the innermost rho points on a circle (where the angle
    # between each point is dz), we will have that the distance between
    # two diametrically oposite points is dx
    #{{{Generally about the relationship between d@ and L@
    # In the case where the boundary lies half between grid points:
    # Let:
    # L@    =   distance from boundary to boundary
    # S@    =   number of internal line segments
    # d@    =   distance between line segments
    #
    # We get
    # L@ = S@*d@ + (contributions from one boundary half between grid cells)*2
    #    = S@*d@ + (1/2)*dx*2
    #    = d@*(S@ + 1)
    #
    # =>
    # d@ = L@/(S@ + 1)
    #}}}
    dx = Lx / (internal_x_segments + 1.0)
    dy = Ly / (internal_y_segments + 1.0)
    
    # Making an array for x-values. (Only used to write to output. CAUTION!!! Not the real r-coordinate)
    rtorusarr = np.ndarray([nx])
    for grid_point_nr in range(nx):
        rtorus = dx*grid_point_nr
        rtorusarr[grid_point_nr] = rtorus
        
    r_grid = np.linspace(0.8-0.2,0.8+0.2,nx) #The grid used for the x-coordinate
    
    Pi = math.pi
    r = 0.8       #The small radius. Assumed to be constant for simplicity. (We are looking at a thin slab)
    R = 2         #The large radius. (Usually referred to as R_0)
    epsilon = r/R #The ratio between the large and small ratio.
    
    #Finding the eta values corresponding to a specific theta(y) value
    def func(eta, Theta):
        return -Theta+eta/(2*Pi)+np.sin(eta)*epsilon/(2*Pi) #Defining a function for later use
    etaarr = np.ndarray([nx,ny])
    Thetaarr  = np.ndarray([nx, ny])
    for grid_point_ntheta in range(ny):
        Theta = grid_point_ntheta*dy - 0.5 #Thetavalues (Y=Theta)
        etause = scipy.optimize.newton(func, x0=Pi, args = (Theta,)) #Finding etavalues corresponding to a theta value using the newton method
        Thetaarr[:,grid_point_ntheta] = grid_point_ntheta*dy - 0.5   #Making a matrix varying in 1 dimension for the metric calculations
        etaarr[:,grid_point_ntheta] = etause                         #Same as above
    
    # Define all constants
    q = Pi        #safety factor. It is constant since r is assumed constant.
    qbar = 3.2       #derivative of the safety factor. It is constant since r is assumed constant.
    
    ShiftAngle = np.ndarray([nx])
    ShiftAngle[:] = (2*Pi)*q

    zetazeta1 = 1/(R**2*(1+r*np.cos(etaarr)/R)**2) #a constant used later. Made to avoid too long expressions. Same for the rest.
    zetazeta21 = (q*np.tan(0.5*etaarr)*(-(1/(R*(1+epsilon)))-((1-epsilon)/((1+epsilon)**2*R))))/(np.sqrt((1-epsilon)/(1+epsilon))*(((1-epsilon)*(np.tan(0.5*etaarr))**2)/(1+epsilon)+1))
    zetazeta22 = 2*qbar*np.arctan(np.sqrt((1-epsilon)/(1+epsilon))*np.tan(etaarr/2))
    zetazeta2 = (zetazeta21+zetazeta22)**2
    zetazeta3 = (4*q**2*(1-epsilon)*(0.5+0.5*(np.tan(etaarr/2))**2)**2)/((1+epsilon)*(((1-epsilon)*(np.tan(etaarr/2))**2)/(1+epsilon)+1)**2*r**2)
    thetazeta2 = (2*q*np.sqrt((1-epsilon)/(1+epsilon))*(0.5 + 0.5*(np.tan(0.5*etaarr)**2)))/((((1-epsilon)*(np.tan(0.5*etaarr)**2))/(1+epsilon)+1)*r**2)
    thetazeta3 = (2*q*np.sqrt((1-epsilon)/(1+epsilon))*(0.5 + 0.5*(np.tan(0.5*etaarr)**2))*np.cos(etaarr))/((((1-epsilon)*(np.tan(0.5*etaarr)**2)/(1+epsilon))+1)*r*R)
    
    #Calculating the final Torus coordinate, Phi (Only used to write to output)
    phitorusarr = np.ndarray([ny,nz])
    sarray = np.ndarray([nz])
    for grid_point_nphi in range(nz):
        s = grid_point_nphi*2*Pi/nz
        sarray[grid_point_nphi] = s/(2*Pi)
        for grid_point_nphi2 in range(ny):
            Theta = grid_point_nphi2*dy
            etause = scipy.optimize.newton(func, x0=Pi, args = (Theta,))
            phitorusarr[grid_point_nphi2,grid_point_nphi] = 2*Pi*(q*Theta - s - (q/(2*Pi))*(2*Pi*Theta-2*np.arctan(np.sqrt((1-epsilon)/(1+epsilon))*np.tan(etause/2))))
            
    
    #Hamada metric coefficients
    grr = 1
    gthetatheta = ((1/r**2)+((2*np.cos(etaarr))/(r*R))+((np.cos(etaarr))**2/R**2)+((np.sin(etaarr))**2/R**2))/((2*Pi)**2)
    gzetazeta = (1/((2*Pi)**2))*(zetazeta1+zetazeta2+zetazeta3)
    grtheta = (np.sin(etaarr))/(2*R*Pi)
    grzeta = (1/(2*Pi))*(zetazeta21+zetazeta22)
    gthetazeta = (1/((2*Pi)**2))*((zetazeta21+zetazeta22)*np.sin(etaarr)/R+thetazeta2+thetazeta3)

    # Intialize the covariant metric tensor elements, contravariant
    # metic tensor elements and J
    # NOTE: x are stored as rows, whereas y are stored as columns.
    #       Note that [0, 0] is the top left corner of the matrix.
    #       If one plots this with a imshow plot with matplotlib the
    #       points in x will correspond to the y-axis on the plot.
    #       In matrix[0, 2], the 0 would correspond to the highest
    #       y-axis value, and 2 would correspond to the corresponding
    #       value of index 2 along the x-axis.
    #       As a side note, the same convention is used in the contour
    #       plot, with the difference that the matrix is mirrored around
    #       the horizontal midaxis, so that in matirx[0,2] the 0 would
    #       correspond to the lowest y-axis value in the plot, and the
    #       2 would correspond to the corresponding value of index 2 along
    #       the x-axis
    # NOTE: Ghost points are important:
    #       g33  = 1/rho^2 (no x derivatives used)
    #       g_33 = rho^2   (x derivative used when deriving Christoffel symbols)
    #       J    = rho     (x derivative used when deriving G)
    # NOTE: Ghost cells in ny is not included in the grid file, but is
    #       in the BOUT++ framework set to be equal the first and last
    #       points in y (as above)
    g11  = np.ndarray([nx, ny])
    g22  = np.ndarray([nx, ny])
    g33  = np.ndarray([nx, ny])
    g12  = np.ndarray([nx, ny])
    g13  = np.ndarray([nx, ny])
    g23  = np.ndarray([nx, ny])


    # In order to get the contravariant elements, it could be enough to
    # calculate the covariant elements. In both cases, we get trouble
    # when assigning the innermost ghost point properly.

    # Set the contravariant metric tensor elements
    g11[:,:] = 1.0
    g22[:,:] = gthetatheta
    # g^{\theta \theta} = 1/rho^2
    # Note that \theta corresponds to z, so we need to set g33 rather
    g33[:, :]  = (qbar*Thetaarr)**2*grr + 2*q*qbar*Thetaarr*grtheta-2*qbar*Thetaarr*grzeta+q**2*gthetatheta-2*q*gthetazeta+gzetazeta
    g12[:,:] = grtheta
    g13[:,:] = qbar*Thetaarr*grr+q*grtheta-grzeta
    g23[:,:] = qbar*Thetaarr*grtheta+q*gthetatheta-gthetazeta
    
    det = g11*(g22*g33-g23*g23)-g12*(g12*g33-g23*g13)+g13*(g12*g23-g22*g13)
    


    # Write the grid file
    with DataFile(name + ".grd.nc", write=True, create=True) as grid:
        # Write the number of points
        grid.write("nx", nx)
        grid.write("ny", ny)
        grid.write("nz", nz)

        # Write the position of the separatrix
        grid.write("ixseps1", ixseps1)
        grid.write("ixseps2", ixseps2)

        # Write the mesh spacingtorusarr
        grid.write("dx", dx)
        grid.write("dy", dy)

        # Write the dimensions of the cylinder
        grid.write("Lx", Lx)
        grid.write("Ly", Ly)
        
        # Write the Theta values
        grid.write("Theta",Thetaarr)
        grid.write("Phi",phitorusarr)
        grid.write("r_grid",r_grid)
        
        # Write the Eta values
        grid.write("Eta",etaarr)
        grid.write("gthetazeta",gthetazeta)
        grid.write("det",det)        
        
        #Write the contravarianta Hamada metric
        grid.write("gthetatheta",gthetatheta)
        grid.write("gzetazeta",gzetazeta)
        grid.write("grtheta",grtheta)
        grid.write("grzeta",grzeta)
        grid.write("gthetazeta",gthetazeta)

        # Write the contravariant metric tensor
        grid.write("g11", g11)
        grid.write("g22", g22)
        grid.write("g33", g33)
        grid.write("g12", g12)
        grid.write("g13", g13)
        grid.write("g23", g23)
        
        # Writing the torus coordinates in terms of field-aligned coordinates
        grid.write("rtorus",rtorusarr)
        grid.write("etatorus",etaarr[MXG,:])
        grid.write("phitorus",phitorusarr)

        #Constants and other values
        grid.write("q",q)
        grid.write("qbar",qbar)
        grid.write("r",r)
        grid.write("R",R)
        grid.write("s",sarray)
        grid.write("ShiftAngle",ShiftAngle)
        grid.write("Pi",Pi)
                
        
if __name__ == "__main__":
    generate_grid()