import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import matplotlib.colors as colors
import numpy as np
from scipy.optimize import least_squares
from astropy.io import fits
import glob
import emcee, corner


plt.rcParams['font.size'] = 16

starname = 'Betelgeuse' #or R_Dor or Betelgeuse or R_Leo
date = '2023-10-04'


#HARDCODE VARIABLES

if starname == 'Betelgeuse':
    nbpix = 15 #lim*2 - 1#50 # nbpix
    lim = 22.5 # size of the image in mas
    pix_size =2*(lim)/nbpix
    print(pix_size)
    Rstar = 21#22.6
    Rshell = 40
    vsys0 = 0

elif starname == 'R_Dor':
    nbpix = 17 #lim*2 - 1#50 # nbpix
    lim = 25.5 # size of the image in mas
    pix_size =2*(lim)/nbpix
    print(pix_size)
    Rstar = 24#22.6
    Rshell = 40
    vsys0 = 9.5


elif starname == 'R_Leo':
    nbpix = 15 #lim*2 - 1#50 # nbpix
    lim = 22.5 # size of the image in mas
    pix_size =2*(lim)/nbpix
    print(pix_size)
    Rstar = 21#22.6
    Rshell = 30
    vsys0 = 10



#reading an actual file
if starname == 'Betelgeuse':
    transition = '28SiOv2'
    input_dir = 'DATA/' + starname+ '/'
elif starname == 'R_Dor' or starname == 'R_Leo':
    transition = 'SiO_v=2_8-7'
    input_dir = 'DATA/' + starname+ '/' + date + '/'


config = open('config.py', 'w')
config.write('nbpix = '+ str(nbpix) + '\n')
config.write('pix_size = '+str(pix_size) +'\n')
config.write('lim = '+str(lim) +'\n')
config.write('Rstar = '+str(Rstar) +'\n')
config.write('Rshell = '+str(Rshell) +'\n')
config.close()


from functions import *

inputfile = glob.glob(input_dir + '*'+ transition+ '*center.vel*.fits')[0]
inputfile_err = glob.glob(input_dir + '*'+ transition+'*center.err*.fits')[0]
betel_data, core_params = read_fits(inputfile, Rstar)
nbpix, pix_size, lim, Rstar = core_params


print(core_params)


tot_pix = nbpix*nbpix
pix_arr = np.arange(tot_pix)

"""
The model parameters:
"""
axis = 45 #in degrees
axis_radian = axis*np.pi/180

vsini = 15 #km/s, should be positive
vsys = 0#4.9 #km/s
vexp = 2#1.7 #km/s, positive value for outflows
# Rshell = 36 # pixel
params = [axis_radian, vsini, vsys, vexp]

#intial value
x0 = [90*np.pi/180, 0, vsys0, 0]

err_scalar = 0.11
err,_ = read_fits(inputfile_err, Rstar)
err = np.ravel(err) + 1E-11


# plot_surface(np.reshape(betel_data, (nbpix, nbpix)), x0, core_params, cmap = 'seismic',boolShow= True)

#to replace by betelgeuse data
vel_map_1d = betel_data #velocity(pix_arr,params)

#back to 2d-array for the visualization
vel_map = np.reshape(vel_map_1d, (nbpix, nbpix))
# print(vel_map_1d)

"""
plot the result
"""
# plot_surface(vel_map, x0, core_params, boolShow= True, boolSave = True, cmap = 'bwr',filename = 'input.png')

"""
export the result into a fits file
"""

output_dir = 'RESULTS/' + starname+ '/'
if starname == 'R_Dor':
    output_dir += date

filename = output_dir + inputfile[len(input_dir):-5]+'_input.fits'
write_to_fits(vel_map, params, core_params, filename = filename)



"""
read fits file and extract data
"""
# filename = glob.glob(DIR+'*.fits')[0]
print('saved as : ', filename)

with fits.open(filename) as hdul:
   data = hdul[0].data
   hdr = hdul[0].header

#need to extract the parameters for the image

# print(data)
plot_surface(data, x0,core_params, cmap = 'seismic',boolShow= True, boolSave = True, filename = output_dir + inputfile[len(input_dir):-5]+'_input.png')
plot_surface(np.reshape(err, (nbpix, nbpix)), x0, core_params,cmap = 'seismic',boolShow= False, boolSave = True, filename = output_dir + inputfile[len(input_dir):-5]+'_err.png')

data_1d =  np.ravel(data)
print(np.shape(betel_data))

#to work directly on th efile
data_1d = np.ravel(vel_map_1d)

# # x0 = np.ones(np.shape(data_1d))


#least_squares


res_lsq = least_squares(model, x0, args=(pix_arr, data_1d))

print('least squared reulst = ', res_lsq.x)
print('PA rotation axis = ', res_lsq.x[0]*180/np.pi)

res_vel_map_1d = velocity(pix_arr,res_lsq.x)
plot_surface(np.reshape(res_vel_map_1d, (nbpix, nbpix)), res_lsq.x, core_params, cmap = 'seismic',boolShow= True, boolSave = True, filename = output_dir + inputfile[len(input_dir):-5]+'least_squares_fit.png')
#
chi2 = np.sum(((data_1d - res_vel_map_1d)**2)/err**2)
print('chi2 = ', chi2)
chi2r = chi2/len(data_1d)
print('chi2r = ', chi2r)



#MCMC fit
ndim = len(params)
nwalkers = 50
initial_guess = x0
pos = initial_guess + 1e-4*np.random.randn(nwalkers, ndim) #nitial_guess + 1e-4 *

chainname= "Figures/chains/chains.h5"
backend = emcee.backends.HDFBackend(chainname)
backend.reset(nwalkers, ndim)


sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=(pix_arr, data_1d, err), backend=backend)

max_n = 100000
#
# # We'll track how the average autocorrelation time estimate changes
index = 0
autocorr = np.empty(max_n)
#
# # This will be useful to testing convergence
old_tau = np.inf
#
# # Now we'll sample for up to max_n steps
for sample in sampler.sample(pos, iterations=max_n, progress=True):
      # Only check convergence every 100 steps
    if sampler.iteration % 100:
        continue
#
#     # Compute the autocorrelation time so far
#     # Using tol=0 means that we'll always get an estimate even
#     # if it isn't trustworthy
    tau = sampler.get_autocorr_time(tol=0)
    autocorr[index] = np.mean(tau)
    index += 1
#
#     # Check convergence
    converged = np.all(tau * 100 < sampler.iteration) #100
    converged &= np.all(np.abs(old_tau - tau) / tau < 0.01)
    if converged:
        print('converged')
        break
    old_tau = tau



n = 100 * np.arange(1, index + 1)
y = autocorr[:index]
plt.plot(n, n / 100.0, "--k")
plt.plot(n, y)
plt.xlim(0, np.nanmax(n))
plt.ylim(0, np.nanmax(y) + 0.1 * (np.nanmax(y) - np.nanmin(y)))
plt.xlabel("number of steps")
plt.ylabel(r"mean $\hat{\tau}$");
plt.show()





tau = sampler.get_autocorr_time()
burnin = int(2 * np.max(tau))
thin = int(0.5 * np.min(tau))
samples = sampler.get_chain(discard=burnin, flat=True, thin=thin)
log_prob_samples = sampler.get_log_prob(discard=burnin, flat=True, thin=thin)
log_prior_samples = sampler.get_blobs(discard=burnin, flat=True, thin=thin)

print("burn-in: {0}".format(burnin))
print("thin: {0}".format(thin))
print("flat chain shape: {0}".format(samples.shape))
print("flat log prob shape: {0}".format(log_prob_samples.shape))
# print("flat log prior shape: {0}".format(log_prior_samples.shape))










# nstep = 1000000
# print("Running MCMC...")
# sampler.run_mcmc(pos, nstep, progress=True)
#Results
# burn_in = int(nstep*0.1)
samples = sampler.get_chain(discard=burnin, flat=True)


params_MCMC = np.mean(samples, axis=0)
print(f"Estimation MCMC : Phi = {params_MCMC[0]:.3f}, vsini = {params_MCMC[1]:.3f}, vsys = {params_MCMC[2]:.3f}, vexp = {params_MCMC[3]:.3f}")

MCMC_vel_map_1d = velocity(pix_arr,params_MCMC)

chi2 = np.sum(((data_1d - MCMC_vel_map_1d)**2)/err**2)
print('chi2 from MCMC = ', chi2)
chi2r = chi2/len(data_1d)
print('chi2r = ', chi2r)

plot_surface(np.reshape(MCMC_vel_map_1d, (nbpix, nbpix)), params_MCMC, core_params, cmap = 'seismic',boolShow= True, boolSave = True, filename = output_dir + inputfile[len(input_dir):-5]+'_result_mcmc.png')


samples[:,0] = samples[:,0]*180/np.pi

fig = corner.corner(
    samples,
    labels=[r"$\Phi$", r"$|v_{rot}\sin(i)|$", r"$v_{sys}$", r"$v_{exp}$"],
    quantiles=[0.16, 0.5, 0.84],
    show_titles=True,
    title_fmt=".2f"
)


plt.savefig(output_dir + inputfile[len(input_dir):-5]+'_corner.png', bbox_inches = 'tight')
plt.show()

