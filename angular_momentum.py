import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import glob

plt.rcParams['font.size'] = 16

G = 6.6743E-11

vrot = 200 #m/s
M = 1
M =  M*1.988E30

d = np.linspace(50,750,10)

Rmas = 5 # limited by ALMA resolution

vrot = 100*(d/50)

R = 1 #Rmas*d
R = R*1.495979E11 #in m

a = 5*R

Jspin = (2/3)*vrot*M*R

Jplanet = Jspin*5

Mplanet = Jplanet/(G*M*a)**0.5

# Mplanet =  Mplanet/5.9722E24           #1.898E27

Mplanet =  Mplanet/1.898E27

plt.figure(figsize = (10,10))
ax = plt.gca()
plt.plot(d, Mplanet)
plt.plot(55, 50, '*', color = 'red', markersize = 10)
plt.annotate('R Dor',xy = (55, 50), color = 'red')
ax.tick_params(direction="in", which = 'both', top = True, right = True)
plt.fill_between(d, Mplanet, y2 = 80, color = 'lightblue', alpha = 0.5)
plt.xlabel('Distance [dc]')
ax.minorticks_on()
plt.xlim(right = np.max(d))
plt.ylim(top = 80)
plt.ylabel(r'Companions mass [$\rm M_\bigoplus$]')
plt.ylabel(r'Companions mass [$\rm M_J$]')
plt.savefig('Sensitivity_ALMA20240.png', bbox_inches = 'tight')

plt.show()
