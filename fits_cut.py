import argparse
import warnings
import glob, os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from astropy.io import fits
from scipy.ndimage import map_coordinates

plt.rcParams['font.size'] = 16
plt.rcParams["font.family"] = "serif"


def load_image(path):

    with fits.open(path) as hdul:
        data = hdul[0].data
        hdr = hdul[0].header
        if data.ndim == 4:
            data = data[0][0]
        elif data.ndim == 3:
            data = data[0]
    return data


def extract_cut(image, width=1, order=1):

    ny, nx = image.shape

    # Work out the end points

    cx, cy = args.center if args.center else ((nx - 1) / 2, (ny - 1) / 2)
    length = min(nx, ny) - 1
    th = np.deg2rad(90+args.angle)
    dx, dy = 0.5 * length * np.cos(th), 0.5 * length * np.sin(th)
    p0, p1 = (cx - dx, cy - dy), (cx + dx, cy + dy)
    p0, p1 = (cx , cy ), (cx + dx, cy + dy)


    x0, y0 = p0
    x1, y1 = p1
    length = np.hypot(x1 - x0, y1 - y0)*pix_size
    n = max(int(np.ceil(length)) + 1, 2)
    t = np.linspace(0.0, 1.0, n)
    xs = x0 + t * (x1 - x0)
    ys = y0 + t * (y1 - y0)
    dist = t * length

    # Unit vector perpendicular to the cut
    ux, uy = (x1 - x0) / length, (y1 - y0) / length
    px, py = -uy, ux
    offsets = np.linspace(-(width - 1) / 2, (width - 1) / 2, int(width))

    profiles = []
    for off in offsets:
        # map_coordinates expects (row, col) = (y, x); NaNs outside -> NaN
        coords = np.vstack([ys + off * py, xs + off * px])
        profiles.append(map_coordinates(image, coords, order=order,
                                        mode="constant", cval=np.nan))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # all-NaN samples
        return dist, np.nanmean(profiles, axis=0), p0, p1

pix_size = 3

ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("fits_file")
ap.add_argument("--center", type=float, nargs=2, metavar=("X", "Y"),
                    help="Cut centre (default: image centre if --angle is used)")
ap.add_argument("--angle", type=float,
                    help="Angle in degrees, counter-clockwise from the +x axis")
ap.add_argument("--length", type=float,
                    help="Cut length in pixels (default: as long as fits)")
ap.add_argument("--width", type=int, default=1,
                    help="Averaging width perpendicular to the cut, in px")


args = ap.parse_args()


starname = 'R_Dor' #or R_Dor or Betelgeuse or R_Leo
dates = ['2023-10-12','2023-10-09','2023-10-06', '2023-10-04', '2023-11-04']
transitions = ['SiO_v=3_8-7','SiO_v=2_8-7', 'H2O_v=0_16-17', 'H2O_v=0_17-16','SO_3_v=0_7-6',  '29SiO_v=1_8-7', 'SO2_v=0_40-40','SO2_v=0_34-34','H2O_v2_2v2_v1','SO2_v=0_33-32','H2O_v=0_9-8'] #,'SO2_v=0_21-22'
filenames = []
errnames = []
labels = []
for date in dates:
    input_dir = 'DATA/' + starname+ '/' + date + '/'
    for transition in transitions:
        inputfiles = glob.glob(input_dir + '*'+ transition+ '*large_scale*center.vel*.fits')
        # print(inputfiles)

        if len(inputfiles) != 0:
            if len(inputfiles)> 1:
                if 'CO' in inputfiles[1] and 'large_scale' not in inputfiles[1]:
                    continue
            print(inputfiles)
            filenames +=inputfiles
            errnames +=[inputfiles[0].replace('vel_0', 'err_0')]
            labels += [date + ' '+ transition]*len(inputfiles)

filenames += glob.glob('RESULTS/R_Dor/results*large_scale*Kep*.fits')
labels += ['Model']
# filenames = glob.glob('DATA/R_Dor/2023-10-04/*large_scale*center.vel*.fits')
fig, (ax_img, ax_cut) = plt.subplots(1, 2, figsize=(20, 8))

for idx, fits_file in enumerate(filenames):
    # print(errnames[idx])
    image = load_image(fits_file)
    dist, profile, p0, p1 = extract_cut(image, width=args.width)
    if idx < len(filenames)-1:

        err = load_image(errnames[idx])

        err_dist, err_profile, p0, p1 = extract_cut(err, width=args.width)

    # Plot: image with the cut overlaid, and the profile


        ax_cut.plot(dist, profile, ls = '-')
        ax_cut.fill_between(dist, profile- err_profile, profile+ err_profile, alpha = 0.3, label = labels[idx])

    else:
        print(filenames[idx])
        finite = image[np.isfinite(image)]
        vmin, vmax = np.percentile(finite, [1, 99])
        ax_img.imshow(image, origin="lower", cmap="seismic", vmin=vmin, vmax=vmax)
        ax_img.set(xlabel="x [mas]", ylabel="y [mas]", title="Image and cut")

        ax_cut.plot(dist, profile, ls = '-', color = 'black',lw = 4, label = labels[-1])

ax_img.plot([p0[0], p1[0]], [p0[1], p1[1]], "-", color = 'lime', lw=2)
ax_img.plot(*p0, "o", ms=4)
ax_cut.set(xlabel="Distance along cut [mas]",
               ylabel="Radial velocity [km/s]",
               title=f"1D cut (width = {args.width} px)")
ax_cut.grid(alpha=0.3)
# ax_cut.set_ylim([vmin,vmax])
ax_cut.legend(bbox_to_anchor=(1.04, 1), borderaxespad=0, ncol = len(dates))
ax_cut.set_xlim(left = 0)

for ax in (ax_img, ax_cut):
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which="both", direction="in",
                    top=True, bottom=True, left=True, right=True)
    ax.tick_params(which="major", length=6)
    ax.tick_params(which="minor", length=3)

fig.tight_layout()
fig.savefig('test.png', dpi=150, bbox_inches="tight")
plt.show()

