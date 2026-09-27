import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import rotate
from scipy.spatial.transform import Rotation
import matplotlib as mpl
import sys


plt.rcParams['font.size'] = 16

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / np.sqrt(np.dot(axis, axis))
    a = np.cos(theta / 2.0)
    b, c, d = -axis * np.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])



args = sys.argv
PA = 45
incl = 0

if (len(args) > 1):
    if "-i" in args:
        incl = float(args[args.index("-i") + 1 ])

    if "-pa" in args:
        PA = float(args[args.index("-pa") + 1 ])



# Make the grid
lim = 3
nbpoints = 50
x, y, z = np.meshgrid(np.linspace(-lim, lim, nbpoints),
                      np.linspace(-lim, lim, nbpoints),
                      np.linspace(-lim, lim, nbpoints))

r = (x*x + y*y + z*z)**0.5

Rshell = 2
Rstar = 1

mask = np.logical_and(r<Rshell, r>Rstar)
mask2 = np.logical_or(x>0, (z*z + y*y) > Rstar)
# mask = np.logical_and(mask, mask2)



x = x[mask]
y = y[mask]
z = z[mask]
r = r[mask]



# theta = np.arccos(zprim/rprim)
# phi = np.arctan2(yprim,xprim)

# phi = np.arctan2(yprimprim,xprimprim)
phi = np.arctan2(y,x)

# Make the direction data for the arrows

vexp = 0
vrot = 1


# expressed in the intrisici coordinate system: x'',y'', z'' co-rotating with the star

phi_v = (-np.sin(phi), np.cos(phi),0*phi)
vx = x*vexp/r + vrot*phi_v[0]
vy = y*vexp/r + vrot*phi_v[1]
vz = z*vexp/r + vrot*phi_v[2]
v = [vx,vy,vz]



# inclination
v = [vx,vy,vz]
#axis2 = [0, np.cos(theta1), np.sin(theta1)]
axis2 = [0, 1, 0]
theta2 = -(90- incl)*np.pi/180
rot2 = rotation_matrix(axis2, theta2)
rot2 = np.linalg.inv(rot2)
v = np.dot(rot2, v)
x = np.dot(rot2, [x,y,z])

[x,y,z] = x

#position angle
axis1 = np.array([1, 0, 0])
# axis1 = np.array([1, 0, 0])
theta1 = -PA*np.pi/180
rot1 = rotation_matrix(axis1, theta1)
rot1 = np.linalg.inv(rot1)
vv = np.dot(rot1, v)
xx = np.dot(rot1, [x,y,z])



#rotation vector
a = [0,0,1.5*lim]
a = np.dot(rot2, a)
a = np.dot(rot1, a)

rot_axis = [[0,a[0]],[0,a[1]],[0,a[2]]]


[xx,yy,zz] = xx



vxx = vv[0]
vyy = vv[1]
vzz = vv[2]


mask2 = np.logical_or(xx>0, (zz*zz + yy*yy) > Rstar)

# xx = xx[mask2]
# yy = yy[mask2]
# zz = zz[mask2]
# vxx = vxx[mask2]
# vyy = vyy[mask2]
# vzz = vzz[mask2]
#



# Color by vx value (radial velocity)
c = vxx
print(np.shape(vxx), np.shape(r))
print(np.min(vxx), np.max(vxx) )
# c = vxx*np.cos(theta2)+ vzz*np.sin(theta2)
c = (c.ravel() - vx.min()) / np.ptp(vx)
# c = np.concatenate((c, np.repeat(c, 2)))
c = plt.cm.bwr_r(c)
# c = plt.cm.jet(c)

ax = plt.figure(figsize = (15,30)).add_subplot(projection='3d')
plt.title(r'$i = $'+str(int(incl)) + r'$^\circ$, PA = '+str(int(PA)) + r'$^\circ$' + ' \n' + r'$\vec{v}_{rot}$ = ' + str(vrot) + r'$.\vec{1}_\phi$, $\vec{v}_{exp}$ = '+ str(vexp) + r'$.\vec{1}_r$')
im = ax.quiver(xx, yy, zz, vxx, vyy, vzz, length=0.2, color = c)

ax.plot(rot_axis[0],rot_axis[1],rot_axis[2], '-', lw = 3, color = 'black')

# ax.quiver(x, y, z, phi_v[0], phi_v[1], phi_v[2], length=0.2, color = c)


# Make data for stellar surface
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
X = Rstar * np.outer(np.cos(u), np.sin(v))
Y = Rstar * np.outer(np.sin(u), np.sin(v))
Z = Rstar * np.outer(np.ones(np.size(u)), np.cos(v))
# Plot the surface
ax.plot_surface(X,Y,Z, color = 'gold')


# ax.contourf(xx, yy, zz, zdir='y', offset=-lim, cmap='bwr')

ax.set_xlim([-lim,2*lim])
ax.set_ylim([-lim,lim])
ax.set_zlim([-lim,lim])

ax.set_aspect('equal')
ax.set_zlabel('DEC axis')
ax.set_ylabel('RA axis (<-E)')
ax.set_xlabel('LOS')
ax.view_init(elev=0, azim=0)


# Hide grid lines
ax.grid(False)

# ax.set_axis_off()
plt.savefig('Figures/3D/i='+str(int(incl)) + '_PA='+str(int(PA)) + '_vrot=' + str(vrot) + r'_vexp='+ str(vexp) + '.png', bbox_inches = 'tight')


#ax.set_box_aspect((1, 1, 1))

plt.show()
