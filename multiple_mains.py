import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import matplotlib.colors as colors
import numpy as np
from scipy.optimize import least_squares
from astropy.io import fits
import glob, sys
import emcee, corner

plt.rcParams['font.size'] = 16
plt.rcParams["font.family"] = "serif"

starname = 'R_Dor' #or R_Dor or Betelgeuse or R_Leo
dates = ['2023-10-12','2023-10-09','2023-10-06', '2023-10-04']


args = sys.argv

suffix = '.clean.fit2'

if (len(args) > 1):
    if "-i" in args:
        starname = args[args.index("-i") + 1 ]

    if "-date" in args:
        date = args[args.index("-date") + 1 ]
    if "-line" in args:
        transition = args[args.index("-line") + 1 ]
    if "-shell" in args:
        suffix2 = '.clean.large_scale'




# transition += suffix

transitions = [transition+ suffix]#, transition + suffix2]





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
    Rstar = 24#21#22.6
    Rshell = int(24*pix_size) #27
    vsys0 = 6 #5.5


elif starname == 'R_Leo':
    nbpix = 15 #lim*2 - 1#50 # nbpix
    lim = 22.5 # size of the image in mas
    pix_size =2*(lim)/nbpix
    print(pix_size)
    Rstar = 21#22.6
    Rshell = 30
    vsys0 = -0.5

# if boolShell:
#
#     Rshell = int(24*pix_size)
#     lim = Rshell + pix_size/2´
#     nbpix = int(2*lim/pix_size)



config = open('config.py', 'w')
config.write('pix_size = '+str(pix_size) +'\n')
config.write('Rstar = '+str(Rstar) +'\n')
config.write('Rshell = '+str(Rshell) +'\n')
config.write('starname = '+'"'+str(starname)+ '"' +'\n')
config.close()


from functions import *
datasets = []

for date in dates:

    for transition in transitions:

        #reading an actual file
        if starname == 'Betelgeuse':
            transition = '28SiOv2'
            transition = '28SiOv1' # 29SiOv0 28SiOv1  28SiOv2 12CO
            input_dir = 'DATA/' + starname+ '/'
        elif starname == 'R_Dor' or starname == 'R_Leo':
            input_dir = 'DATA/' + starname+ '/' + date + '/'

        output_dir = 'RESULTS/' + starname+ '/'
        if starname == 'R_Dor' and len(dates)==1 :
            output_dir += date

        # READING THE INPUT


        inputfiles = glob.glob(input_dir + '*'+ transition+ '*center.vel*.fits')
        if len(inputfiles) <1:
            continue
        inputfile = inputfiles[0]
        inputfile_err = inputfiles[0].replace('vel_0', 'err_0') #glob.glob(input_dir + '*'+ transition+'*center.err*.fits')[0]

        input_data = read_fits(inputfile, Rstar)
        vel_map_1d = np.copy(input_data)

        plot_surface(input_data, None, cmap = 'seismic',boolShow= True, boolSave = True, filename = output_dir + inputfile[len(input_dir):-5]+'_input.png')

        #PREPARING THE PARAMETERS

        nbpix = np.shape(input_data)[0]
        pix_arr = np.arange(nbpix*nbpix)

        #The error map
        err = error_map(inputfile_err)

        #to work directly on the file
        data_1d = np.ravel(vel_map_1d)

        datasets.append((pix_arr, data_1d, err))



#The intial parameters:
x0 = [300*np.pi/180, 0, vsys0, 0, 0]



#MCMC fit
ndim = len(x0)
nwalkers = 50


lows = np.array([1.5*np.pi, 0.0, vsys0-2, -5, -1.0])
highs = np.array([2.5 * np.pi, 10.0,vsys0+2,5, 1.0])
pos = lows + (highs - lows) * np.random.rand(nwalkers, ndim)

# pos = x0 + 1e-2*np.random.randn(nwalkers, ndim) #nitial_guess + 1e-4 *

chainname= "Figures/chains/chains.h5"
backend = emcee.backends.HDFBackend(chainname)
backend.reset(nwalkers, ndim)


sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior_total, args=[datasets],  backend=backend)

max_n = 100000


index = 0
autocorr = np.empty(max_n)

old_tau = np.inf

for sample in sampler.sample(pos, iterations=max_n, progress=True):
      # Only check convergence every 100 steps
    if sampler.iteration % 100:
        continue

     # Compute the autocorrelation time so far
     # Using tol=0 means that we'll always get an estimate even if it isn't trustworthy
    tau = sampler.get_autocorr_time(tol=0)
    autocorr[index] = np.mean(tau)
    index += 1

    converged = np.all(tau * 100 < sampler.iteration) #100
    converged &= np.all(np.abs(old_tau - tau) / tau < 0.01)
    if converged:
        print('converged')
        break
    old_tau = tau



n = 100 * np.arange(1, index + 1)
y = autocorr[:index]

tau = sampler.get_autocorr_time()
burnin = int(2 * np.max(tau))
thin = int(0.5 * np.min(tau))



"""
Plotting the convergence
"""

fig, axes = plt.subplots(ndim+1, figsize=(10, 3 * ndim), sharex=True)
labels=[r"mean $\hat{\tau}$",r"$\Phi$", r"$|v_{rot}\sin(i)|$", r"$v_{sys}$", r"$v_{exp}$", r"$\alpha$"]
# axes[0].plot(n, n, "--k")
axes[0].plot(n, y)
# axes[0].set_xlim(0, np.nanmax(n))
# axes[0].ylim(0, np.nanmax(y) + 0.1 * (np.nanmax(y) - np.nanmin(y)))
axes[0].tick_params(direction="in", which = 'both', top = True, right = True)
axes[0].set_ylabel(labels[0])

chain = sampler.get_chain()
nsteps, nwalkers, ndim = chain.shape

for i in range(1,ndim+1):
    if i == 1:
        chain[:, :, i-1] *= 180/np.pi

    axes[i].plot(chain[:, :, i-1], color="k", alpha=0.3)

    # Optional: Draw a vertical line to show where your burn-in cutoff is
    axes[i].axvline(burnin, color="r", linestyle="--", label="Burn-in Cutoff")

    axes[i].set_ylabel(labels[i])
    axes[i].yaxis.set_label_coords(-0.1, 0.5)
    axes[i].tick_params(direction="in", which = 'both', top = True, right = True)

axes[-1].set_xlabel("Step number")

plt.tight_layout()

plt.savefig(output_dir +'diagnostic.png', bbox_inches = 'tight')
plt.show()






samples = sampler.get_chain(discard=burnin, flat=True, thin=thin)
log_prob_samples = sampler.get_log_prob(discard=burnin, flat=True, thin=thin)
log_prior_samples = sampler.get_blobs(discard=burnin, flat=True, thin=thin)

print("burn-in: {0}".format(burnin))
print("thin: {0}".format(thin))
print("flat chain shape: {0}".format(samples.shape))
print("flat log prob shape: {0}".format(log_prob_samples.shape))
# print("flat log prior shape: {0}".format(log_prior_samples.shape))






samples = sampler.get_chain(discard=burnin, flat=True)






params_MCMC = np.mean(samples, axis=0)
print(f"Estimation MCMC : Phi = {params_MCMC[0]:.3f}, vsini = {params_MCMC[1]:.3f}, vsys = {params_MCMC[2]:.3f}, vexp = {params_MCMC[3]:.3f} , alpha = {params_MCMC[4]:.3f}")



# loop on all the datasets to compute each chi2r

for date in dates:

    for transition in transitions:

        # READING THE INPUT
        inputfiles = glob.glob(input_dir + '*'+ transition+ '*center.vel*.fits')
        if len(inputfiles) <1:
            continue
        inputfile = inputfiles[0]
        inputfile_err = inputfiles[0].replace('vel_0', 'err_0') #glob.glob(input_dir + '*'+ transition+'*center.err*.fits')[0]

        input_data = read_fits(inputfile, Rstar)
        vel_map_1d = np.copy(input_data)

        #PREPARING THE PARAMETERS

        nbpix = np.shape(input_data)[0]
        pix_arr = np.arange(nbpix*nbpix)

        #The error map
        err = error_map(inputfile_err)

        #to work directly on the file
        data_1d = np.ravel(vel_map_1d)

        #computing the chi2r
        MCMC_vel_map_1d = velocity(pix_arr,params_MCMC)

        chi2 = np.nansum(((data_1d - MCMC_vel_map_1d)**2)/err**2)
        print('chi2 from MCMC = ', chi2)
        chi2r = chi2/np.count_nonzero(data_1d)
        print('chi2r = ', chi2r)

#plot the final results
MCMC_vel_map = np.reshape(MCMC_vel_map_1d, (nbpix, nbpix))
plot_surface(MCMC_vel_map, params_MCMC, cmap = 'seismic',chi2r = chi2r, boolShow= True, boolSave = True, filename = output_dir + 'result_mcmc.png')


"""
save data and results into a  fits file
"""

filename = output_dir + inputfile[len(input_dir):-5]+ '_results.fits'
print('saved as : ', filename)
write_to_fits([input_data, MCMC_vel_map,input_data - MCMC_vel_map], params_MCMC, filename = filename)

"""
Corner plot
"""


samples[:,0] = samples[:,0]*180/np.pi

figure = corner.corner(
    samples,
    labels=[r"$\Phi$", r"$|v_{rot}\sin(i)|$", r"$v_{sys}$", r"$v_{exp}$", r"$\alpha$"],
    quantiles=[0.16, 0.5, 0.84],
    show_titles=True,
    title_fmt=".2f",
    figsize = (20,20)
)

axes = np.array(figure.axes).reshape((ndim, ndim))
print(np.shape(axes))

# Loop over the histograms
for yi in range(ndim): # Loop over the diagonal
    print(yi)
    axes[yi, yi].axvline(params_MCMC[yi], color="r")

    for xi in range(yi):
        ax = axes[yi, xi]
        ax.tick_params(direction="in", which = 'both', top = True, right = True)

        ax.axvline(params_MCMC[xi], color="r")
        # ax.axvline(value2[xi], color="r")
        ax.axhline(params_MCMC[yi], color="r")
        # ax.axhline(value2[yi], color="r")
        ax.plot(params_MCMC[xi], params_MCMC[yi], "sr")
        # ax.plot(value2[xi], value2[yi], "sr")

# plt.subplots_adjust(bottom=0, right=1.0, top=1.0)


figure.subplots_adjust(wspace=0, hspace=0)


plt.savefig(output_dir + 'corner.png', bbox_inches = 'tight')
plt.show()

