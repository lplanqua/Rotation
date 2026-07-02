import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
from scipy.optimize import least_squares
from astropy.io import fits
import glob
import emcee, corner

plt.rcParams['font.size'] = 16


#HARDCODE VARIABLES

nbpix = 15 #lim*2 - 1#50 # nbpix
lim = 22.5 # size of the image in mas
pix_size =2*(lim)/nbpix
print(pix_size)
Rstar = 21#22.6
Rshell = 25

"""
Some useful functions
"""

def model(params, pix_arr, data): #x t y
    model = velocity(pix_arr,params)
    return data - model

def log_prior(params):
    axis_radian, vsini, vsys, vexp = params
    if vsini <0:
        return -np.inf
    if axis_radian< 0:
        return -np.inf
    if axis_radian > 2*np.pi:
        return -np.inf
    # if vexp >100:       #slow wind
    #     return -np.inf
    # if vexp <-100:
    #     return -np.inf
    return 0

# Vraisemblance (log-likelihood)
def log_likelihood(params, pix_arr, data):
    axis_radian, vsini, vsys, vexp = params
    # sigma = err_scalar*np.ones(np.shape(pix_arr))
    sigma = err
    model = velocity(pix_arr,params)
    return -0.5 * np.sum(((data - model) / sigma) ** 2 + np.log(2 * np.pi * sigma ** 2))

# Posterior = prior + likelihood
def log_posterior(params, pix_arr, data):
    lp = log_prior(params)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(params, pix_arr, data)



def plot_surface(vel_map, params, cmap = 'jet', boolShow =True, boolSave = True, filename = None):


    axis_radian, vsini, vsys, vexp = params
    if filename == None:
        filename = 'Figures/' +'axis' + str(int(axis_radian*180/np.pi))+ '_vsini'+"{:.1f}".format(vsini)+'_vexp'+"{:.1f}".format(vexp)+'_vsys'+"{:.1f}".format(vsys)+'_Rshell'+"{:.1f}".format(Rshell)
        filename += '.png'

    plt.figure(figsize = (12,12))
    # h = plt.contourf(X, Y, Rs)
    if params[0] == 0:
        plt.title('INPUT')
    else:
        plt.title(r'$\Phi$ = '+ str(int(axis_radian*180/np.pi))+r'$^\circ$, $|v_{rot}\sin(i)|$ = '+ str(np.round(vsini,1))+ r' km/s'+ '\n'+ r'$v_{exp}$ = '+ str(np.round(vexp,1))+ r' km/s, $v_{sys}$ = '+ str(np.round(vsys,1))+ r' km/s, $R_{shell}$ = '+ str(np.round(Rshell,1))+' pix')
    # vel_map= np.array(vel_map, dtype=float)
    Cmap = plt.get_cmap(cmap)
    Cmap.set_bad(color='gray', alpha=0.25)
    h = plt.imshow(vel_map, extent = [-lim, lim, -lim, lim], origin = 'lower', cmap = Cmap, norm=colors.CenteredNorm(vcenter = 0))

    plt.axis('scaled')
    ax = plt.gca()

    ax.tick_params(direction="in", which = 'both', top = True, right = True)
    if vsini != 0:
        ax.annotate("", xytext=(1.2*Rstar*np.sin(axis_radian), -1.2*Rstar*np.cos(axis_radian)), xy=(-1.2*Rstar*np.sin(axis_radian), 1.2*Rstar*np.cos(axis_radian)),
                    arrowprops=dict(arrowstyle="->",facecolor='black', lw = 2))
    ax.minorticks_on()
    cbar = plt.colorbar(label = 'RV [km/s]', pad = 0.001,  fraction=0.01,aspect=100)
    # cbar.add_lines([vsys])
    cbar.ax.axhline(y=vsys, c='lime', lw = 5)
    # cbar.axvline(y=vsys)
    plt.xlabel(r'$\Delta \alpha$ [mas] (<--E)')
    plt.ylabel(r'$\Delta \delta$ [mas] (N-->)')
    if boolSave:
        plt.savefig(filename, bbox_inches = 'tight')
        print('saved as : ',filename)
    if boolShow:
        plt.show()
    else:
        plt.close()
    return True


def write_to_fits(vel_map, params, filename = None):
    axis_radian, vsini, vsys, vexp = params
    if filename == None:
        filename = 'Figures/' + 'axis' + str(int(axis_radian*180/np.pi))+ '_vsini'+"{:.1f}".format(vsini)+'_vexp'+"{:.1f}".format(vexp)+'_vsys'+"{:.1f}".format(vsys)+'_Rshell'+"{:.1f}".format(Rshell)
        filename += '.fits'
    hdr = fits.Header()
    hdr['BUNIT'] = 'km/s'
    hdr['BTYPE'] = 'velocity'

    hdr['CRVAL1'] = 0
    hdr['CRVAL2'] = 0

    hdr['CRPIX1'] = nbpix/2
    hdr['CRPIX2'] = nbpix/2

    hdr['CRDELT1'] = -pix_size #should be negative for RA axis positive toward left
    hdr['CRDELT2'] = pix_size

    # the units
    hdr['CUNIT1'] = 'mas'
    hdr['CUNIT1'] = 'mas'

    hdu = fits.PrimaryHDU(vel_map, hdr)

    hdu.writeto(filename, overwrite=True)
    return True


def read_fits(filename):
    print('Using file ',filename)
    with fits.open(filename) as hdul:
        data = hdul[0].data
        hdr = hdul[0].header
    #if stokes params
    if data.ndim == 4:
        data = data[0][0]
    elif data.ndim == 3:
        data = data[0]
    data = np.nan_to_num(data) # Replace NaN values with 0
    nbpix = np.shape(data)[0]
    pix_size = 3
    #pix_size =2*(lim)/nbpix
    lim = nbpix*pix_size/2
    Rstar = 21#22.6
    core_params = nbpix, pix_size, lim, Rstar

    return data, core_params


def velocity(pix_arr,params):
    axis_radian, vsini, vsys, vexp = params

    Y, X = np.unravel_index(pix_arr, (nbpix,nbpix)) #X and Y swapped to speak in cartesian coordinates
    X = X*pix_size -lim +pix_size*0.5
    Y = Y*pix_size -lim +pix_size*0.5

    R = (X*X + Y*Y)**0.5
    # print(X)
    #
    mask = R > Rstar
    # X[mask] = np.nan
    # Y[mask] = np.nan

    axis_perp = X*np.cos(axis_radian) + Y*np.sin(axis_radian) #x'
    axis_para =  -X*np.sin(axis_radian) + Y*np.cos(axis_radian)#y'
    vrot = axis_perp*vsini/Rshell
    vtot = vrot + vsys - vexp*(Rshell - R)/Rshell
    vtot[mask] = 0#np.nan
    # h = plt.imshow(np.reshape(R, (nbpix, nbpix)), extent = [-lim, lim, -lim, lim], origin = 'lower', cmap = 'bwr')
    # plt.show()
    # print(stop)

    return vtot.ravel()


#reading an actual file
input_dir = 'DATA/Betelgeuse/'
inputfile = glob.glob(input_dir + '*center.vel*.fits')[2]
inputfile_err = glob.glob(input_dir + '*center.err*.fits')[2]
betel_data, core_params = read_fits(inputfile)
nbpix, pix_size, lim, Rstar = core_params


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
x0 = [0*np.pi/180, 0, 0.0, 0]
err_scalar = 0.11
err,_ = read_fits(inputfile_err)
err = np.ravel(err) + 1E-9


# plot_surface(np.reshape(betel_data, (nbpix, nbpix)), x0, cmap = 'seismic',boolShow= True)

#to replace by betelgeuse data
vel_map_1d = betel_data #velocity(pix_arr,params)

#back to 2d-array for the visualization
vel_map = np.reshape(vel_map_1d, (nbpix, nbpix))
# print(vel_map_1d)

"""
plot the result
"""
# plot_surface(vel_map, x0, boolShow= True, boolSave = True, cmap = 'bwr',filename = 'input.png')

"""
export the result into a fits file
"""

output_dir = 'RESULTS/Betelgeuse/'

filename = output_dir + inputfile[len(input_dir):-5]+'_input.fits'
write_to_fits(vel_map, params, filename = filename)



"""
read fits file and extract data
"""
# DIR = 'Figures/'

# params = [45*np.pi/180, 5, 0.0, 2.0, 36.0]
# axis_radian, vsini, vsys, vexp, Rshell = params


# filename = glob.glob(DIR+'*.fits')[0]
print('saved as : ', filename)

with fits.open(filename) as hdul:
   data = hdul[0].data
   hdr = hdul[0].header

#need to extract the parameters for the image

# print(data)
plot_surface(data, x0, cmap = 'seismic',boolShow= True, boolSave = True, filename = output_dir + inputfile[len(input_dir):-5]+'_input.png')
plot_surface(np.reshape(err, (nbpix, nbpix)), x0, cmap = 'seismic',boolShow= False, boolSave = True, filename = output_dir + inputfile[len(input_dir):-5]+'_err.png')

data_1d =  np.ravel(data)
# print(np.shape(data_1d))

# # x0 = np.ones(np.shape(data_1d))


#least_squares


res_lsq = least_squares(model, x0, args=(pix_arr, data_1d))

print('least squared reulst = ', res_lsq.x)
print('PA rotation axis = ', res_lsq.x[0]*180/np.pi)

res_vel_map_1d = velocity(pix_arr,res_lsq.x)
plot_surface(np.reshape(res_vel_map_1d, (nbpix, nbpix)), res_lsq.x,  cmap = 'seismic',boolShow= False, boolSave = True, filename = output_dir + inputfile[len(input_dir):-5]+'least_squares_fit.png')
#
chi2 = np.sum(((data_1d - res_vel_map_1d)**2)/err**2)
print('chi2 = ', chi2)



#MCMC fit
ndim = len(params)
nwalkers = 50
initial_guess = x0
pos = np.random.randn(nwalkers, ndim) #nitial_guess + 1e-4 *

chainname= "chains.h5"
backend = emcee.backends.HDFBackend(chainname)
backend.reset(nwalkers, ndim)


sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=(pix_arr, data_1d))

# max_n = 100000
#
# # We'll track how the average autocorrelation time estimate changes
# index = 0
# autocorr = np.empty(max_n)
#
# # This will be useful to testing convergence
# old_tau = np.inf
#
# # Now we'll sample for up to max_n steps
# for sample in sampler.sample(pos, iterations=max_n, progress=True):
#     # Only check convergence every 100 steps
#     if sampler.iteration % 100:
#         continue
#
#     # Compute the autocorrelation time so far
#     # Using tol=0 means that we'll always get an estimate even
#     # if it isn't trustworthy
#     tau = sampler.get_autocorr_time(tol=0)
#     autocorr[index] = np.mean(tau)
#     index += 1
#
#     # Check convergence
#     converged = np.all(tau * 100 < sampler.iteration)
#     converged &= np.all(np.abs(old_tau - tau) / tau < 0.01)
#     if converged:
#         break
#     old_tau = tau





















nstep = 1000000
print("Running MCMC...")
sampler.run_mcmc(pos, nstep, progress=True)
#Results
burn_in = int(nstep*0.1)
samples = sampler.get_chain(discard=burn_in, flat=True)


params_MCMC = np.mean(samples, axis=0)
print(f"Estimation MCMC : Phi = {params_MCMC[0]:.3f}, vsini = {params_MCMC[1]:.3f}, vsys = {params_MCMC[2]:.3f}, vexp = {params_MCMC[3]:.3f}")

MCMC_vel_map_1d = velocity(pix_arr,params_MCMC)

chi2 = np.sum(((data_1d - MCMC_vel_map_1d)**2)/err**2)
print('chi2 from MCMC = ', chi2)

plot_surface(np.reshape(MCMC_vel_map_1d, (nbpix, nbpix)), params_MCMC,  cmap = 'seismic',boolShow= True, boolSave = True, filename = output_dir + inputfile[len(input_dir):-5]+'_result_mcmc.png')


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

