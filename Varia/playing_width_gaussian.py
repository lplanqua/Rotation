import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.size'] = 16

def gauss(x, mu, sigma, amp):
    return amp/(sigma * np.sqrt(2 * np.pi)) *np.exp( - (x - mu)**2 / (2 * sigma**2))


#radial velocity
maxi = 20
x = np.linspace(-maxi, maxi, 41)



blue_emission = gauss(x,-4,7, 4)
central_absorption = gauss(x, 0, 5,-5)


plt.figure(figsize =(12,6))
plt.plot(x, blue_emission, c= 'blue', label = 'Blueshifted emission')
print(x[np.argmax(blue_emission)],np.max(blue_emission))
plt.plot(x[np.argmax(blue_emission)],np.max(blue_emission), 'v' , color = 'blue')
plt.plot(x, central_absorption, c= 'black', ls = '--',label = 'Central absortion')
plt.plot(x, central_absorption+ blue_emission, c = 'black', lw = 2, label = 'Sum')

plt.plot(x[np.argmax(np.abs(blue_emission+ central_absorption))],np.max(np.abs(blue_emission+central_absorption)), 'v' , color = 'k')

plt.axvline(x= 0, c = 'red')
plt.legend(loc = 4)
ax = plt.gca()

ax.tick_params(direction="in", which = 'both', top = True, right = True, color = 'black')
ax.minorticks_on()

ax.set_xlabel(r'RV [km/s]')
ax.set_ylabel(r'Flux')
plt.savefig('sum_emission_and_absorption.png', bbox_inches = 'tight')
plt.show()
