import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import rotate
from scipy.spatial.transform import Rotation


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


ax = plt.figure().add_subplot(projection='3d')


# Make the grid
lim = 3
nbpoints = 30
x, y, z = np.meshgrid(np.linspace(-lim, lim, nbpoints),
                      np.linspace(-lim, lim, nbpoints),
                      np.linspace(-lim, lim, nbpoints))

r = (x*x + y*y + z*z)**0.5

Rshell = 2
Rstar = 1

mask = np.logical_and(r<Rshell, r>Rstar)
mask2 = np.logical_or(x>0, (z*z + y*y) > Rstar)
# mask = np.logical_and(mask, mask2)

PA = -90
incl = 0

axis1 = (2,0)
# axis1 = (2,1)
xprim = rotate(x, PA, axes = axis1, reshape=False)
yprim = rotate(y, PA, axes = axis1, reshape=False)
zprim = rotate(z, PA, axes = axis1, reshape=False)

axis2 = (1,0)
xprimprim = rotate(xprim, incl, axes = axis2, reshape=False)
yprimprim = rotate(yprim, incl, axes = axis2, reshape=False)
zprimprim = rotate(zprim, incl, axes = axis2, reshape=False)

rprim = (xprimprim*xprimprim + yprimprim*yprimprim + zprimprim*zprimprim)**0.5

xprim = xprim[mask]
yprim = yprim[mask]
zprim = zprim[mask]

xprimprim = xprimprim[mask]
yprimprim = yprimprim[mask]
zprimprim = zprimprim[mask]


rprim = rprim[mask]
x = x[mask]
y = y[mask]
z = z[mask]
r = r[mask]

#intrinsic angle according to Euler rotation

# gamma = PA*np.pi/180
# beta = (90-incl)*np.pi/180
# mat = np.array([[np.cos(beta), np.sin(beta)*np.sin(gamma), np.sin(beta)*np.cos(gamma)],
#        [0, np.cos(gamma), -np.sin(gamma)],
#        [-np.sin(beta), np.cos(beta)*np.sin(beta), np.cos(beta)*np.cos(gamma)]])



# theta = np.arccos(zprim/rprim)
# phi = np.arctan2(yprim,xprim)

# phi = np.arctan2(yprimprim,xprimprim)
phi = np.arctan2(y,x)

# Make the direction data for the arrows

vexp = 0
vrot = 1
# phi_v = (np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi),np.cos(theta))


# phi_v = rotate(phi_v, PA,  axes = (0,2), reshape=False)
# phi_v = rotate(phi_v, incl,axes = axis2, reshape=False)


# expressed in the intrisici coordinate system: x'',y'', z'' co-rotating with the star

phi_v = (-np.sin(phi), np.cos(phi),0*phi)
vx = x*vexp/r + vrot*phi_v[0]
vy = y*vexp/r + vrot*phi_v[1]
vz = z*vexp/r + vrot*phi_v[2]
v = [vx,vy,vz]





axis1 = np.array([1, 0, 0])
theta1 = PA*np.pi/180

rot = rotation_matrix(axis1, theta1)
print(theta1,axis1)

rot = np.linalg.inv(rot)

vv = np.dot(rot, v)


# rot = Rotation.from_rotvec(theta1*axis1)
# vv = rot.apply(np.transpose(v))
# vv = np.transpose(vv)

vx = v[0]
vy = v[1]
vz = v[2]


v = [vx,vy,vz]
axis2 = [0, np.cos(theta1), np.sin(theta1)]
axis2 = [0, 1, 0]
theta2 = incl*np.pi/180

# vv = np.dot(rotation_matrix(axis2, theta2), v)

vxx = vv[0]
vyy = vv[1]
vzz = vv[2]


# Color by vx value (radial velocity)
c = vx
c = vxx*np.cos(theta2)+ vzz*np.sin(theta2)
c = (c.ravel() - c.min()) / np.ptp(c)
# c = np.concatenate((c, np.repeat(c, 2)))
c = plt.cm.bwr_r(c)
# c = plt.cm.jet(c)

ax.quiver(x, y, z, vxx, vyy, vzz, length=0.2)#, color = c)

# ax.quiver(x, y, z, phi_v[0], phi_v[1], phi_v[2], length=0.2, color = c)


# Make data for stellar surface
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
X = Rstar * np.outer(np.cos(u), np.sin(v))
Y = Rstar * np.outer(np.sin(u), np.sin(v))
Z = Rstar * np.outer(np.ones(np.size(u)), np.cos(v))
# Plot the surface
ax.plot_surface(X,Y,Z, color = 'gold')

ax.set_xlim([-lim,2*lim])
ax.set_ylim([-lim,lim])
ax.set_zlim([-lim,lim])

ax.set_aspect('equal')
ax.set_zlabel('DEC axis (N->)')
ax.set_ylabel('RA axis (<-E)')
ax.set_xlabel('LOS')
ax.view_init(elev=10, azim=90)
#ax.set_box_aspect((1, 1, 1))

plt.show()
