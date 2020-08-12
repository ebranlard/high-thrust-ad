import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
import weio
from VortexSurfaceFlowField import surface_u, determineVorticitySurface, plotVelocityField
from HighThrust import a_Ct

# --- Parameters
R               = 1       # Rotor radius
nRings          = 300     # Number of rings used over the vorticity surface region before the semi-inf cylinder
includeCylinder = True    # Include semi-infinite cylinder
includeCylinder = True    # Include semi-infinite cylinder
U0              = 1       # Freestream
n_rcp           = 160      # Number of control points in radial direction
n_xcp           = 300     # Number of control points in axial direction
x_max           = 70*R    # Max point after which the surface should be in equilibrium

# --- Control points used for velocity field
# rcp = np.linspace(0,5*R, n_rcp)
# xcp = np.linspace(-3*R,x_max*1.5, n_xcp)

tolerance   = 0.01
alpha_relax = 0.05
iter_max    = 500
dR_diff = 0.1*R
dt = R/U0


# --- Parametric

vCT   = [0.1  ,0.5  ,0.8  ,1.1  ,1.5  ,1.5  ,1.7]
vEddy = [0.000,0.000,0.001,0.100,0.100,0.400,0.8]

vCT   = [1.3 , 1.0 , 0.7    , 1.5 ]
vEddy = [0.30, 0.10, 0.0001,  0.5]

# vCT   = [0.1, 0.5]
# vEddy = [0.000, 0.000]

# vCT   = [0.8 , 1.5  ,1.7]
# vEddy = [0.00,0.400,0.8]

va = a_Ct(vCT, method='HAWC2')
print(va)
x_surf = np.linspace(0, x_max, nRings)
method = 'iteration'
# method = 'optimize'
vR_surf=[]
vGeom=[]

for i,(CT,a,Eddy) in enumerate(zip(vCT,va,vEddy)):
    print('>>>>>>>>> C_T={:.2f} a={:.2f} Eddy={:.3f}'.format(CT,a,Eddy))
    gamma  = -2*a*U0 # intensity of the vorticity surface [m/s]
    R_surf, res, geom = determineVorticitySurface(x_surf, R, U0, gamma, EddyVisc=Eddy, includeCylinder=includeCylinder, method=method, tolerance=tolerance, alpha_relax=alpha_relax, iter_max=iter_max)
    vR_surf.append(R_surf)
    vGeom.append(geom)


# --- Plotting final surface
fig,ax = plt.subplots(1, 1, sharey=False, figsize=(7.4,2.8)) # (6.4,4.8)
for i,(CT,a,Eddy,geom) in enumerate(zip(vCT,va,vEddy,vGeom)):

    x_rings, R_rings, x_cyl, R_cyl, vGamma, gamma_cyl = geom
    vGamma = np.asarray(vGamma)
    I=np.where(vGamma/vGamma[1]>0.01)[0]

    ax.plot(x_surf[I]/(2*R),vR_surf[i][I]/(2*R)   ,label = 'C_T={:.2f} a={:.2f} EddyVisc={:.3f}'.format(CT,a,Eddy))
    ax.plot(x_surf[I]/(2*R),-vR_surf[i][I]/(2*R) )
# if a<0.5:
#     R_surf_th = R* np.sqrt( (1-a) / (1-a*(1+x_surf/R/np.sqrt(1+(x_surf/R)**2 )  ))    )
#     ax.plot(x_surf/R,R_surf_th/R,label = 'Theory'          )
ax.set_xlabel('x/D [-]')
ax.set_ylabel('R/R0 [-]')
ax.legend(fontsize=6)
ax.set_ylim([-1.5,1.5])
ax.set_xlim([-2,8])
ax.set_title('Multi Eddy Visc')




xcp = np.linspace(-4*R,16*R, n_xcp)
rcp = np.linspace(-3*R,3*R , n_rcp)
for i,(CT,a,Eddy) in enumerate(zip(vCT,va,vEddy)):
    print('Velocity field...')
    gamma  = -2*a*U0 # intensity of the vorticity surface [m/s]
    nRings=len(x_surf)
    Rcp, Xcp = np.meshgrid(rcp,xcp)
    Rcp_abs = np.abs(Rcp)
    ur, ux, geom = surface_u(x_surf, vR_surf[i], gamma, Xcp, Rcp_abs, nRings = nRings, includeCylinder=includeCylinder, EddyVisc= Eddy)
    x_rings, R_rings, x_cyl, R_cyl, vGamma, gamma_cyl = geom
    ux += U0 # Adding free stream velocity
    ur[Rcp<0] *=-1

    # Save it
    np.savez('CT{:.2f}.dat'.format(CT), x_surf=x_surf, R_surf=vR_surf[i], ux=ux, ur=ur, Rcp=Rcp, Xcp=Xcp, rcp=rcp, xcp=xcp)



# --- Streamlines
for i,(CT,a,Eddy) in enumerate(zip(vCT,va,vEddy)):
    gamma  = -2*a*U0 # intensity of the vorticity surface [m/s]
    xcp = np.linspace(-4*R,16*R, n_xcp)
    rcp = np.linspace(-3*R,3*R , n_rcp)
    fig, ax, geom = plotVelocityField(xcp, rcp, x_surf, R_surf, R, U0, gamma, Eddy, includeCylinder)
    ax.set_ylim([-1.5,1.5])
    ax.set_xlim([-2,8])
    ax.set_aspect('equal')
    ax.set_xlabel('x [D]')
    ax.set_ylabel('z [D]')



plt.show()
