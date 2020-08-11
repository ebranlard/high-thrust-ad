import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
import weio
from VortexSurfaceFlowField import surface_u, determineVorticitySurface, plotVelocityField
# --- Parameters
R               = 1       # Rotor radius
nRings          = 300     # Number of rings used over the vorticity surface region before the semi-inf cylinder
includeCylinder = True    # Include semi-infinite cylinder
includeCylinder = True    # Include semi-infinite cylinder
EddyVisc        = 0.101   # Scaling of vorticity downstream, 0=no decay
U0              = 1       # Freestream
a               = 0.45    # Axial induction
gamma           = -2*a*U0 # intensity of the vorticity surface [m/s]
n_rcp           = 50      # Number of control points in radial direction
n_xcp           = 100     # Number of control points in axial direction
x_max           = 90*R    # Max point after which the surface should be in equilibrium

# --- Control points used for velocity field
# rcp = np.linspace(0,5*R, n_rcp)
# xcp = np.linspace(-3*R,x_max*1.5, n_xcp)

tolerance   = 0.001
alpha_relax = 0.1
iter_max    = 500
dR_diff = 0.1*R
dt = R/U0

vR_surf=[]
x_surf = np.linspace(0, x_max, nRings)
vMethods = ['iteration' , 'optimize']
vMethods = ['iteration']
# vMethods = ['optimize']

for method in vMethods:
    R_surf, resi, geom = determineVorticitySurface(x_surf, R, U0, gamma, EddyVisc=EddyVisc, includeCylinder=includeCylinder, method=method, tolerance=tolerance, alpha_relax=alpha_relax, iter_max=iter_max)
    vR_surf.append(R_surf)


# --- Plotting final surface
if a<0.5:
    R_surf_th = R* np.sqrt( (1-a) / (1-a*(1+x_surf/R/np.sqrt(1+(x_surf/R)**2 )  ))    )

fig,ax = plt.subplots(1,1)
for R_surf,method in zip(vR_surf,vMethods):
    ax.plot(x_surf/R,R_surf/R   ,label = method)
if a<0.5:
    ax.plot(x_surf/R,R_surf_th/R,label = 'Theory'          )
ax.set_xlabel('x/R [-]')
ax.set_ylabel('R/R0 [-]')
ax.legend()

# # --- Plotting velocity field
xcp = np.linspace(-3*R,x_max*1.5, n_xcp)
rcp = np.linspace(0,5*R, n_rcp)
fig, ax, geom = plotVelocityField(xcp, rcp, x_surf, R_surf, R, U0, gamma, EddyVisc, includeCylinder)

# --- Plotting vorticity
x_rings, R_rings, x_cyl, R_cyl, vGamma, gamma_cyl = geom
fig,ax = plt.subplots(1,1)
dx_surf = x_surf[-1]-x_surf[-2] # Distance between the two last surfs
Gamma_surf = gamma*dx_surf        # Circulation of each vortex surf [m^2/s]
ax.plot(x_surf  ,vGamma, label='Ring intensities')
ax.plot(x_surf  , Gamma_surf  * np.exp( - EddyVisc*x_surf), label='computed')
ax.set_xlabel('x')
ax.set_ylabel('Gamma rings')
ax.legend()


plt.show()
