import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
import weio
def stretchByRegion(Lx, dx_wanted, xStartStretch, stretch_fact):
    """ stretch a region:
         x<xStartStretch : no stretching
         x>xStartStretch : progressively stretch by stretch_fact """
    dxStart=dx_wanted
    rS = np.array([0,xStartStretch,xStartStretch*1.00001,Lx])
    sS = np.array([1,1,stretch_fact, stretch_fact])
    vx       = []
    vStretch = []
    x0       = 0
    dxPrev   = dx_wanted
    while x0<Lx:
        stretch = np.interp(x0,rS,sS)
        vStretch.append(stretch)
        vx.append(x0)
        dx     = dxPrev*stretch
        dxPrev = dx
        x0     = x0 + dx        
    vx = np.array(vx)
    return vx, vStretch

# --- Parameters
D=1
LxHalf1 = 5*D
LxHalf2 = 15*D
LyHalf  = 10*D
LzHalf  = LyHalf
nPerD = 40
dxStart=D/nPerD

# Stretching in y and z
vyHalf, vStretch = stretchByRegion( Lx=LyHalf, dx_wanted=D/nPerD,  xStartStretch=1.0*D, stretch_fact=1.2)

# Upstream
vxHalf1, vStretch = stretchByRegion( Lx=LxHalf1, dx_wanted=D/nPerD,  xStartStretch=D, stretch_fact=1.2)

# Downstream
vxHalf2, vStretch = stretchByRegion( Lx=LxHalf2, dx_wanted=D/nPerD,  xStartStretch=4*D, stretch_fact=1.1)




vy =np.sort(np.unique( np.concatenate((vyHalf, -vyHalf))))
vx =np.sort(np.unique( np.concatenate((-vxHalf1, vxHalf2))))
xRange = np.array([np.min(vx), np.max(vx)])
yRange = np.array([np.min(vy), np.max(vy)])

vx_norm= (vx-xRange[0])/(xRange[1]-xRange[0])
vy_norm= (vy-yRange[0])/(yRange[1]-yRange[0])


# fig,ax = plt.subplots(1,1)
# ax.plot(vStretch)
# 
fig,ax = plt.subplots(1,1)
for xx in vx:
    ax.plot([xx,xx],[-LyHalf,LyHalf],'k')
for yy in vy:
    ax.plot([-LxHalf1,LxHalf2],[yy,yy],'k')
# ax.plot(    , label='')
# ax.set_xlabel('')
ax.set_aspect('equal')
ax.set_ylabel('')
ax.legend()




print('nalu_abl_mesh:')
print('  output_db: mesh_alm.exo')
print('')
print('  spec_type: bounding_box')
print('')
print('  fluid_part_name: fluid')
print('')
print('  vertices:')
print('  - [-630.0,  -1260.0, -1260.0]')
print('  - [ 1890.0,  1260.0,  1260.0]')
print('')
print('')
print('  xmin_boundary_name: west')
print('')
print('  xmax_boundary_name: east')
print('')
print('  ymin_boundary_name: south')
print('')
print('  ymax_boundary_name: north')
print('')
print('  zmin_boundary_name: lower')
print('')
print('  zmax_boundary_name: upper')
print('  mesh_dimensions: [{:d}, {:d}, {:d}]'.format(len(vx)-1,  len(vy)-1, len(vy)-1  ))
print('')
print('  x_spacing:')
print('    spacing_type: user_specified_spacing')
print('    node_spacing_ratios: ',end='')
print('[' + ','.join([str(x) for x in vx_norm]) + ']')
print('')
print('  y_spacing:')
print('    spacing_type: user_specified_spacing')
print('    node_spacing_ratios: ',end='')
print('[' + ','.join([str(y) for y in vy_norm]) + ']')
print('')
print('  z_spacing:')
print('    spacing_type: user_specified_spacing')
print('    node_spacing_ratios: ',end='')
print('[' + ','.join([str(y) for y in vy_norm]) + ']')
print('')




