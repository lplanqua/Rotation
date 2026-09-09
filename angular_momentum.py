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

Mplanet =  Mplanet/5.9722E24           #1.898E27

# Mplanet =  Mplanet/1.898E27


# Orbital period around the Sun, in days
orbital_period_days = [
    87.969,       # Mercury
    224.701,      # Venus
    365.256,      # Earth
    686.980,      # Mars
    4332.59,      # Jupiter
    10759.22,     # Saturn
    30688.5,      # Uranus
    60182.0       # Neptune
]

# Mass, in Earth masses
mass_earth = [
    0.0553,       # Mercury
    0.8150,       # Venus
    1.0000,       # Earth
    0.1074,       # Mars
    317.83,       # Jupiter
    95.16,        # Saturn
    14.54,        # Uranus
    17.15         # Neptune
]

# Planet names, in the same order
planets = [
    "Mercury",
    "Venus",
    "Earth",
    "Mars",
    "Jupiter",
    "Saturn",
    "Uranus",
    "Neptune"
]



plt.figure(figsize = (10,10))
ax = plt.gca()
plt.plot(d, Mplanet)
plt.plot(55, 50*317, '*', color = 'red', markersize = 10)
plt.annotate('R Dor',xy = (55, 50*317), color = 'red')
ax.tick_params(direction="in", which = 'both', top = True, right = True)
plt.plot(orbital_period_days, mass_earth, 'd', color = 'darkgreen')
for i, planet in enumerate(planets):
    plt.annotate(planet, xy = (orbital_period_days[i], mass_earth[i]), color = 'darkgreen')
plt.fill_between(d, Mplanet, y2 = 80*317, color = 'lightblue', alpha = 0.5)

plt.xlabel('Orbital Period [d]')
ax.minorticks_on()
# plt.ylim(top = 80)
plt.yscale('log')
plt.xscale('log')
plt.ylabel(r'Companions mass [$\rm M_\bigoplus$]')
plt.savefig('Sensitivity_ALMA20240_1.png', bbox_inches = 'tight')

plt.figure(figsize = (10,10))
ax = plt.gca()
plt.plot(d, Mplanet)
plt.plot(55, 50*317, '*', color = 'red', markersize = 10)
plt.annotate('R Dor',xy = (55, 50*317), color = 'red')
ax.tick_params(direction="in", which = 'both', top = True, right = True)
plt.fill_between(d, Mplanet, y2 = 80*317, color = 'lightblue', alpha = 0.5)
plt.xlabel('Distance [dc]')
ax.minorticks_on()
plt.xlim(right = np.max(d))
# plt.ylim(top = 80)
plt.yscale('log')
plt.ylabel(r'Companions mass [$\rm M_\bigoplus$]')
plt.savefig('Sensitivity_ALMA20240_2.png', bbox_inches = 'tight')

plt.show()
