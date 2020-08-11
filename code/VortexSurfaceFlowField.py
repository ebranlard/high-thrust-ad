import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.integrate as si
import scipy.optimize

from InducedVelocities import *
try:
    from pybra.curves import streamQuiver
except:
    def streamQuiver(*args, **kwargs):
        pass


# --------------------------------------------------------------------------------}
# --- Vorticity surface determination 
# --------------------------------------------------------------------------------{
def determineVorticitySurface(x_surf, R, U0, gamma, EddyVisc=0.0, includeCylinder=True, method='optimize', tolerance=0.001, alpha_relax=0.1, iter_max=500):
    """ 
    Determine the equlibirum vorticity surface based on a vorticity distribution gamma

    INPUTS:    
       x_surf: x-coordinate, where the surface will be determined

       R, U0, gamma: rotor radius,free stream, intensity of sheet 
       EddyVisc: EddyVisc parameter, if 0, ring intensity is constant, otherwise an exponential decay is used

       includeCylinder: if true, a semi-infinite vortex cylinder is used to continue the surface to infinity

       method: 'iteration' (relaxation)  of 'optimize'
       tolerance: tolerance for convergence of surface
       alpha_relax: relaxation  (for `iteration` method)
       iter_max: max number of iterations  (for `iteration` method)

   OUTPUTS:
       R_surf: r-coordinate, where the surface is defined, at the x_surf locations
    """
    x_max = np.max(x_surf)
    nRings= len(x_surf)
    R_surf_init = np.zeros(x_surf.shape) + R
    R_surf      = R_surf_init

    def epsilon(R_surf, outputMore=False):
        # Compute velocity above and below the surface points
        #         dR_diff = 0.1*R
        #         ur_low, ux_low, _ = surface_u(x_surf, R_surf, gamma, Xcp = x_surf, Rcp = R_surf-dR_diff, nRings=nRings, includeCylinder=includeCylinder, scaleIntensity=scaleIntensity) 
        #         ur_hig, ux_hig, _ = surface_u(x_surf, R_surf, gamma, Xcp = x_surf, Rcp = R_surf+dR_diff, nRings=nRings, includeCylinder=includeCylinder, scaleIntensity=scaleIntensity) 
        #         ur = (ur_low+ur_hig)/2
        #         ux = (ux_low+ux_hig)/2

        # Compute velocity at the surface points
        ur, ux, geom = surface_u(x_surf, R_surf, gamma, Xcp = x_surf, Rcp = R_surf , nRings=nRings, includeCylinder=includeCylinder, EddyVisc= EddyVisc) 
        ux += U0 # Adding free stream velocity

        dR     = np.concatenate(([0], si.cumtrapz(ur/ux, x_surf) ))
        #dGamma = np.concatenate(([0], si.cumtrapz(ur/(ux*R_surf), x_surf) ))

        R_new          = R + dR
        # Constraints
        R_new[0]       = R       # Constraint
        R_new[R_new<0] = 0
        iFW            = np.argmin(np.abs(x_surf-x_max*0.90))
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
            # relaxation
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

    residual,_,geom = epsilon(R_surf, outputMore=True)

    return R_surf, residual, geom



# --------------------------------------------------------------------------------}
# --- Vorticity Surface velocity field
# --------------------------------------------------------------------------------{
def surface_u(x_surf, R_surf, gamma, Xcp, Rcp, nRings =50, includeCylinder=True, EddyVisc=0.0):
    """ 
    Computes the induced velocity field by a vorticity surface of constant intensity

    INPUTS:
       x_surf: axial points defining the surface (array)
               NOTE: last point will be continued by a semi-inf cylinder
       R_surf: surface radius, for each axial point (array)
       gamma: intensity of the vorticity surface (typically <0) [m/s]
       Xcp : scalar/array/matrix of Control Points axial coordinates
       Rcp : scalar/array/matrix of Control Points radial coordinates (same dimension as Xcp)

       nRings: number of rings used to discretize the surface up until the last value of x_surf
       includeCylinder: if true, a semi-infinite vortex cylinder is used to continue the surface to infinity
       EddyVisc: EddyVisc parameter, if 0, ring intensity is constant, otherwise an exponential decay is used

    """
    Xcp=np.asarray(Xcp)
    Rcp=np.asarray(Rcp)

    # --- Vortex rings (Quadrature points/discretization of the surface)
    x_rings = np.linspace(np.min(x_surf), np.max(x_surf), nRings) # Equi-spacing the rings 
    R_rings = np.interp(x_rings, x_surf, R_surf) # Interpolating the ring radii based on the input surface distribution 

    dx_rings = x_rings[-1]-x_rings[-2] # Distance between the two last rings
    Gamma_ring = gamma*dx_rings        # Circulation of each vortex ring [m^2/s]
    epsilon_ring  =  dx_rings/2        # Basic epsilon, scaling can be done later (see TODO)

    # --- Induced velocity from all rings
    ur_tot, ux_tot = np.zeros(Xcp.shape), np.zeros(Xcp.shape)
    vGamma=[]
    for i_ring, (x_ring, R_ring) in enumerate(zip(x_rings, R_rings)):
        if i_ring==0:
            Gamma_ring_scaled=Gamma_ring/2 # first ring represent half a ring..
        else:
            Gamma_ring_scaled   = Gamma_ring  * R_rings[0]/R_ring  # TODO insert scaling here
            Gamma_ring_scaled   = Gamma_ring  * np.exp( - EddyVisc*x_ring)
#             Gamma_ring_scaled   = Gamma_ring_scaled  * np.exp( - EddyVisc*x_ring)
            #Gamma_ring_scaled   = Gamma_ring_scaled  /(4*np.pi*EddyVisc*x_ring)
        vGamma.append(Gamma_ring_scaled)

        epsilon_ring_scaled = epsilon_ring # TODO insert epsilon hack here

        ur, ux =  ring_u_polar(Rcp, Xcp, Gamma=Gamma_ring_scaled, r0=R_ring, z0=x_ring, epsilon=epsilon_ring_scaled, reg_method='Saffman')
        ur_tot += ur
        ux_tot += ux

    # --- Cylinder induced velocity 
    x_cyl     = x_rings[-1] + dx_rings/2  # Cylinder starts at dx/2 after the last ring
    R_cyl     = R_rings[-1]               # Cylinder has same radius as last ring
    gamma_cyl = vGamma[-1]/dx_rings       # Cylinder intensity is the same as last ring
    if includeCylinder:
        ur_cyl, ux_cyl = vc_tang_u(Rcp, Rcp*0, Xcp-x_cyl, gamma_t=gamma_cyl, R=R_cyl,polar_out=True)
        ur_tot += ur_cyl
        ux_tot += ux_cyl

    geom = (x_rings, R_rings, x_cyl, R_cyl, vGamma, gamma_cyl)

    # --- Total velocity
    return ur_tot, ux_tot, geom


def plotVelocityField(xcp, rcp, x_surf, R_surf, R, U0, gamma, EddyVisc=0, includeCylinder=True):
    print('Velocity field...')
    nRings=len(x_surf)
    Rcp, Xcp = np.meshgrid(rcp,xcp)
    Rcp_abs = np.abs(Rcp)
    ur, ux, geom = surface_u(x_surf, R_surf, gamma, Xcp, Rcp_abs, nRings = nRings, includeCylinder=includeCylinder, EddyVisc= EddyVisc)
    x_rings, R_rings, x_cyl, R_cyl, vGamma, gamma_cyl = geom
    ux += U0 # Adding free stream velocity

    ur[Rcp<0] *=-1


    speed = ux/U0
    speedMin=0
    speedMax=1.2
    levels=np.linspace(speedMin,speedMax,100)

    fig,ax = plt.subplots(1,1, figsize=(9,4))
    im=ax.contourf(Xcp/(2*R),Rcp/(2*R),speed, levels=levels) # the easiest way to get contourf to comply with cmap is to give levels
    #rseed=np.linspace(0/R,1.6/R,7)
    rseed=np.linspace(-(1.4*R)/(2*R),(1.4*R)/(2*R), 23)
    start=np.array([(rseed*0-2*R)/(2*R),rseed])
    print(start)
    sp=ax.streamplot(xcp/(2*R),rcp/(2*R),ux.T,ur.T,color='k',start_points=start.T,linewidth=0.7,density=30,arrowstyle='-')
    # sp=ax.streamplot(xcp/R,rcp/R,ux.T,ur.T,color='k',linewidth=0.7,density=30,arrowstyle='-')
    qv=streamQuiver(ax,sp,n=[5]*len(rseed),scale=40,angles='xy')

    ax.plot(x_surf /(2*R), R_surf/(2*R),'k-', lw=3, label='Input surface')
    ax.plot(x_surf /(2*R),-R_surf/(2*R),'k-', lw=3                       )

    ax.plot(x_rings/(2*R), R_rings/2   , 'o'  , label='Rings', ms=2)
    ax.plot(x_rings/(2*R),-R_rings/2   , 'o'                 , ms=2)

    ax.plot(np.array([x_cyl,np.max(xcp)])/(2*R), np.array([R_cyl, R_cyl])/(2*R), 'k--'  , label='Vortex cylinder', lw=3)
    fig.colorbar(im)
    # ax.set_aspect('equal')
#     ax.legend()
    # 
    return fig, ax, geom


if __name__=='__main__':
    # --------------------------------------------------------------------------------}
    # --- Main code 
    # --------------------------------------------------------------------------------{
    # --- Parameters
    R               = 1       # Rotor radius
    nRings          = 50     # Number of rings used over the vorticity surface region before the semi-inf cylinder
    includeCylinder = True    # Include semi-infinite cylinder
    U0              = 10      # Freestream
    a               = 0.4     # Axial induction
    gamma           = -2*a*U0 # intensity of the vorticity surface [m/s]
    n_rcp           = 50      # Number of control points in radial direction
    n_xcp           = 300     # Number of control points in axial direction
    x_max           = 10*R    # Max

    # --- Control points used for velocity field
    rcp = np.linspace(0,5*R, n_rcp)
    xcp = np.linspace(-3*R,x_max*1.5, n_xcp)

    # --- Definition of the vorticity surface (NOTE: last point continued by a semi-inf cylinder)

    # Example 1: cylinder
    # x_surf = np.array([0,x_max ])*R
    # R_surf = np.array([1,1] )*R 
    # rcp = rcp[ np.abs(rcp-R)>0.01*R ] # hack to remove points close cylinder

    # Example 2: expanding cylinder
    x_surf = np.linspace(0,x_max, 100)
    R_surf = R* np.sqrt( (1-a) / (1-a*(1+x_surf/R/np.sqrt(1+(x_surf/R)**2 )  ))    )

    # TODO  Insert your own surface here


    # --- Velocity at special points
    ur0, ux0, _ = surface_u(x_surf, R_surf, gamma, [0               ], [0], nRings = nRings, includeCylinder=includeCylinder)
    urm, uxm, _ = surface_u(x_surf, R_surf, gamma, [np.max(x_surf)/2], [0], nRings = nRings, includeCylinder=includeCylinder)
    urw, uxw, _ = surface_u(x_surf, R_surf, gamma, [np.max(x_surf)  ], [0], nRings = nRings, includeCylinder=includeCylinder)
    print('Center velocity     :', ur0/U0, ux0/U0)
    print('Mid surf velocity   :', urm/U0, uxm/U0)
    print('Wake start velocity :', urw/U0, uxw/U0)

    # --- Velocity on a grid of points for plotting
    Rcp, Xcp = np.meshgrid(rcp,xcp)
    ur, ux, geom = surface_u(x_surf, R_surf, gamma, Xcp, Rcp, nRings = nRings, includeCylinder=includeCylinder)
    x_rings, R_rings, x_cyl, R_cyl, vGamma, gamma_cyl = geom
    ux += U0 # Adding free stream velocity


    # --- Plotting
    # speed = np.sqrt(ur**2 + ux**2)/U0
    speed = ux/U0
    speedMin=0
    speedMax=1.1
    levels=np.linspace(speedMin,speedMax,10)

    fig,ax = plt.subplots(1,1, figsize=(9,4))
    im=ax.contourf(Xcp/R,Rcp/R,speed, levels=levels) # the easiest way to get contourf to comply with cmap is to give levels
    fig.colorbar(im)
    ax.plot(x_surf/R,R_surf/R,'k-', lw=3, label='Input surface')
    ax.plot(x_rings/R, R_rings, 'o'  , label='Rings', ms=2)
    ax.plot(np.array([x_cyl,np.max(xcp)])/R, np.array([R_cyl, R_cyl])/R, 'k--'  , label='Vortex cylinder', lw=3)
    ax.set_aspect('equal')
    ax.legend()


    plt.show()
