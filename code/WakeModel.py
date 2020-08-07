import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
import weio
from VortexSurfaceFlowField import surface_u
import scipy.integrate as si
import scipy.optimize
try:
    from pybra.curves import streamQuiver
except:
    def streamQuiver(*args, **kwargs):
        pass
# --- Parameters
R               = 1       # Rotor radius
nRings          = 300     # Number of rings used over the vorticity surface region before the semi-inf cylinder
includeCylinder = True    # Include semi-infinite cylinder
scaleIntensity  = False    # Include semi-infinite cylinder
includeCylinder = True    # Include semi-infinite cylinder
scaleIntensity  = True    # Include semi-infinite cylinder
EddyVisc        =0.011
U0              = 1       # Freestream
a               = 0.45    # Axial induction
gamma           = -2*a*U0 # intensity of the vorticity surface [m/s]
n_rcp           = 50      # Number of control points in radial direction
n_xcp           = 100     # Number of control points in axial direction
x_max           = 60*R    # Max point after which the surface should be in equilibrium

# --- Control points used for velocity field
# rcp = np.linspace(0,5*R, n_rcp)
# xcp = np.linspace(-3*R,x_max*1.5, n_xcp)

method = 'iteration'
method = 'optimize'
tolerance   = 0.001
alpha_relax = 0.1
iter_max    = 500
dR_diff = 0.1*R
dt = R/U0

vR_surf=[]
vMethods = ['iteration' , 'optimize']
# vMethods = ['iteration']
vMethods = ['optimize']

for method in vMethods:
    x_surf      = np.linspace(0, x_max, nRings)
    R_surf_init = np.zeros(x_surf.shape) + R
    R_surf      = R_surf_init

    def epsilon(R_surf, outputMore=False):
        # Compute velocity above and below the surface points
#         ur_low, ux_low, _ = surface_u(x_surf, R_surf, gamma, Xcp = x_surf, Rcp = R_surf-dR_diff, nRings=nRings, includeCylinder=includeCylinder, scaleIntensity=scaleIntensity) 
#         ur_hig, ux_hig, _ = surface_u(x_surf, R_surf, gamma, Xcp = x_surf, Rcp = R_surf+dR_diff, nRings=nRings, includeCylinder=includeCylinder, scaleIntensity=scaleIntensity) 
#         ur = (ur_low+ur_hig)/2
#         ux = (ux_low+ux_hig)/2

        # Compute velocity at the surface points
        ur, ux, geom = surface_u(x_surf, R_surf, gamma, Xcp = x_surf, Rcp = R_surf        , nRings=nRings, includeCylinder=includeCylinder, scaleIntensity= scaleIntensity, EddyVisc= EddyVisc) 
#         ux = U0 # Adding free stream velocity
#         ux = 1 # Adding free stream velocity
        ux += U0 # Adding free stream velocity

        dR     = np.concatenate(([0], si.cumtrapz(ur/ux, x_surf) ))
        dGamma = np.concatenate(([0], si.cumtrapz(ur/(ux*R_surf), x_surf) ))

        R_new          = R + dR
        # Constraints
        R_new[0]       = R       # Constraint
        R_new[R_new<0] = 0
        iFW            = np.argmin(np.abs(x_surf-x_max*0.85))
        R_new[iFW:]    = R_new[iFW]

        residual =  (R_surf-R_new)/R

#         ax.plot(x_surf, np.exp(dGamma), label='')
#         ax.plot(x_surf, R_surf, label='')
    #         ax.plot(x_surf, ur, label='')
        
        if outputMore: 
            return residual, R_surf, geom
        else:
            return residual

    # --- Iteration method
    #fig,ax = plt.subplots(1,1)
    if method =='iteration':
        iteration =0
        residual_scalar=np.inf
        while residual_scalar>tolerance and iteration<iter_max:
            residual = epsilon(R_surf)

            R_new = R_surf - residual*R 
            # relazation
            R_surf = alpha_relax * R_new  + (1-alpha_relax) * R_surf
            iteration+=1
            residual_scalar = np.sum(np.abs(residual))
            print('Iteration {:6d} - residual: {:9.1e} - R/R0=[{:.3f} {:.3f} {:.3f}]'.format(iteration, residual_scalar, R_surf[0], R_surf[int(nRings/2)], R_surf[-1] ))

    elif method=='optimize':
        result = scipy.optimize.anderson(epsilon, R_surf_init, f_tol=tolerance, verbose = True)
        R_surf = result/R
        #res, R_surf, _, = epsilon()
    else:
        raise NotImplementedError()


    vR_surf.append(R_surf)
# ax.set_xlabel('')
# ax.set_ylabel('')
# ax.legend()


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
print('Velocity field...')
rcp = np.linspace(0,5*R, n_rcp)
xcp = np.linspace(-3*R,x_max*1.5, n_xcp)
Rcp, Xcp = np.meshgrid(rcp,xcp)
ur, ux, geom = surface_u(x_surf, R_surf, gamma, Xcp, Rcp, nRings = nRings, includeCylinder=includeCylinder, scaleIntensity=scaleIntensity, EddyVisc= EddyVisc)
x_rings, R_rings, x_cyl, R_cyl, vGamma, gamma_cyl = geom
ux += U0 # Adding free stream velocity
speed = ux/U0
speedMin=0
speedMax=1.1
levels=np.linspace(speedMin,speedMax,10)

fig,ax = plt.subplots(1,1, figsize=(9,4))
im=ax.contourf(Xcp/R,Rcp/R,speed, levels=levels) # the easiest way to get contourf to comply with cmap is to give levels

rseed=np.linspace(0/R,1.6/R,7)
start=np.array([rseed*0-0.1*R,rseed])
print(start)
sp=ax.streamplot(xcp/R,rcp/R,ux.T,ur.T,color='k',start_points=start.T,linewidth=0.7,density=30,arrowstyle='-')
# sp=ax.streamplot(xcp/R,rcp/R,ux.T,ur.T,color='k',linewidth=0.7,density=30,arrowstyle='-')
qv=streamQuiver(ax,sp,n=[5,5,5,5,5,5,5],scale=40,angles='xy')
ax.plot(x_surf/R,R_surf/R,'k-', lw=3, label='Input surface')
ax.plot(x_rings/R, R_rings, 'o'  , label='Rings', ms=2)
ax.plot(np.array([x_cyl,np.max(xcp)])/R, np.array([R_cyl, R_cyl])/R, 'k--'  , label='Vortex cylinder', lw=3)
fig.colorbar(im)
ax.set_aspect('equal')
ax.legend()
# 

# --- Plotting vorticity
fig,ax = plt.subplots(1,1)
dx_surf = x_surf[-1]-x_surf[-2] # Distance between the two last surfs
Gamma_surf = gamma*dx_surf        # Circulation of each vortex surf [m^2/s]
ax.plot(x_surf  ,vGamma, label='stupid')
# ax.plot(x_surf  , Gamma_surf  * np.exp( - 0.1*x_surf), label='computed')
ax.set_xlabel('')
ax.set_ylabel('')
ax.legend()
plt.show()

plt.show()
