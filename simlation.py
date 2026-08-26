import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits


# ============================================================
# SYNTHETIC ALMA SPECTRAL-LINE CUBE
#
# Stellar surface + rotating + expanding envelope
#
# Coordinate convention:
#
#   x, y : plane of sky [mas]
#   z    : line of sight [mas]
#
#   z < 0 : front side, toward observer
#   z > 0 : back side, away from observer
#
# Therefore radial expansion on the front side has
# negative LOS velocity = blueshift.
# ============================================================


# ============================================================
# USER PARAMETERS
# ============================================================

# ------------------------------------------------------------
# Stellar surface
# ------------------------------------------------------------

stellar_diameter_mas = 20.0

stellar_radius_mas = (
    stellar_diameter_mas / 2.0
)


# ------------------------------------------------------------
# Expanding envelope
# ------------------------------------------------------------

use_expansion = True

# Outer radius of the envelope [mas]

shell_radius_mas = 30.0

# Maximum expansion velocity at the stellar surface [km/s]

vexp_max = 5.0

# Envelope emissivity relative to the stellar surface

envelope_emissivity = 1.0


# ------------------------------------------------------------
# Stellar surface emission
# ------------------------------------------------------------

use_stellar_surface = True

stellar_emissivity = 1.0


# ------------------------------------------------------------
# Solid-body rotation
# ------------------------------------------------------------

use_rotation = True

# Rotation velocity gradient
#
# Units:
#
#       km/s/mas
#
# For example:
#
#       1 km/s/mas
#
# means that the intrinsic rotational velocity is
# 10 km/s at 10 mas from the rotation axis.

rotation_gradient = 1.0


# ------------------------------------------------------------
# Geometry
# ------------------------------------------------------------

# Inclination of rotation axis
#
# 0 deg  = rotation axis along line of sight
#           -> no observed rotation
#
# 90 deg = rotation axis in plane of sky
#           -> maximum observed rotation

inclination_deg = 60.0


# Position angle of the rotation axis
#
# PA = 0:
#   projected rotation axis is along +Y
#   velocity gradient is along X.
#
# Increasing PA rotates the velocity gradient
# counter-clockwise.

position_angle_deg = 0.0


# ------------------------------------------------------------
# Systemic velocity
# ------------------------------------------------------------

v_sys = 0.0       # km/s


# ------------------------------------------------------------
# Intrinsic line profile
# ------------------------------------------------------------

# Gaussian sigma

sigma_v = 2.0     # km/s


# ------------------------------------------------------------
# Spectral axis
# ------------------------------------------------------------

nchan = 201

v_min = -30.0     # km/s
v_max = +30.0     # km/s


# ------------------------------------------------------------
# Spatial grid
# ------------------------------------------------------------

fov_mas = 25.0

nx = 128
ny = 128


# ------------------------------------------------------------
# Line-of-sight grid through envelope
# ------------------------------------------------------------

nz = 256


# ------------------------------------------------------------
# FITS output
# ------------------------------------------------------------

output_fits = (
    "alma_rotating_expanding_envelope.fits"
)


# ============================================================
# VELOCITY AXIS
# ============================================================

velocity = np.linspace(
    v_min,
    v_max,
    nchan
)

channel_width = (
    velocity[1] - velocity[0]
)


# ============================================================
# SPATIAL GRID
# ============================================================

x = np.linspace(
    -fov_mas / 2.0,
    fov_mas / 2.0,
    nx
)

y = np.linspace(
    -fov_mas / 2.0,
    fov_mas / 2.0,
    ny
)

X, Y = np.meshgrid(x, y)


# Projected distance from stellar centre

rho = np.sqrt(
    X**2 + Y**2
)


# Stellar surface projected disc

stellar_disc = (
    rho <= stellar_radius_mas
)


# ============================================================
# ROTATION GEOMETRY
# ============================================================

PA = np.radians(
    position_angle_deg
)

inclination = np.radians(
    inclination_deg
)


# Rotate coordinates so that:
#
# x_rot = direction of maximum rotational LOS velocity
# y_rot = projected rotation axis
#

x_rot = (
    X * np.cos(PA)
    + Y * np.sin(PA)
)

y_rot = (
    -X * np.sin(PA)
    + Y * np.cos(PA)
)


# ============================================================
# SOLID-BODY ROTATION
# ============================================================

if use_rotation:

    #
    # Solid-body rotation:
    #
    # v_phi = Omega * R
    #
    # Projected LOS velocity:
    #
    # v_rot_LOS =
    #       rotation_gradient
    #       * x_rot
    #       * sin(inclination)
    #

    v_rot_surface = (
        rotation_gradient
        * x_rot
        * np.sin(inclination)
    )

else:

    v_rot_surface = np.zeros_like(X)


# ============================================================
# STELLAR SURFACE CUBE
# ============================================================

cube = np.zeros(
    (nchan, ny, nx),
    dtype=np.float64
)


if use_stellar_surface:

    #
    # Every surface pixel has a Gaussian line profile
    # centered on its local rotational velocity.
    #

    for ichan, v in enumerate(velocity):

        profile = np.exp(
            -0.5
            * (
                (
                    v
                    - v_sys
                    - v_rot_surface
                )
                / sigma_v
            )**2
        )

        cube[ichan] += (
            stellar_emissivity
            * profile
            * stellar_disc
        )


# ============================================================
# 3D ENVELOPE GRID
# ============================================================

#
# z is positive AWAY from the observer.
#
# Therefore:
#
# z < 0 -> front hemisphere -> blueshift
# z > 0 -> back hemisphere  -> redshift
#
# We intentionally only integrate z < 0 here.
#

z_values = np.linspace(
    -shell_radius_mas,
    0.0,
    nz
)

dz = (
    z_values[1]
    - z_values[0]
)


# ============================================================
# ENVELOPE INTEGRATION
# ============================================================

if use_expansion:

    print(
        "Integrating expanding envelope..."
    )

    for iz, z in enumerate(z_values):

        # ----------------------------------------------------
        # 3D radius
        # ----------------------------------------------------

        r = np.sqrt(
            X**2
            + Y**2
            + z**2
        )


        # ----------------------------------------------------
        # Select material inside envelope
        # ----------------------------------------------------

        envelope_mask = (
            (r >= stellar_radius_mas)
            &
            (r <= shell_radius_mas)
        )


        if not np.any(envelope_mask):
            continue


        # ----------------------------------------------------
        # Expansion velocity law
        # ----------------------------------------------------
        #
        # v_exp(R_star)  = vexp_max
        #
        # v_exp(R_shell) = 0
        #
        # Linear decrease with radius.
        #

        v_exp = np.zeros_like(r)

        v_exp[envelope_mask] = (
            vexp_max
            * (
                shell_radius_mas
                - r[envelope_mask]
            )
            / (
                shell_radius_mas
                - stellar_radius_mas
            )
        )


        # ----------------------------------------------------
        # LOS projection of radial expansion
        # ----------------------------------------------------
        #
        # Radial velocity vector:
        #
        #     v = v_exp * r_vector / r
        #
        # Therefore:
        #
        #     v_LOS = v_exp * z/r
        #
        # Since z < 0 on the front side,
        # this is automatically negative:
        #
        #     v_LOS < 0  -> blueshift
        #

        v_exp_los = np.zeros_like(r)

        valid = (
            envelope_mask
            & (r > 0)
        )

        v_exp_los[valid] = (
            v_exp[valid]
            * z
            / r[valid]
        )


        # ----------------------------------------------------
        # Rotation
        # ----------------------------------------------------
        #
        # For solid-body rotation, the LOS rotational
        # velocity depends only on projected distance
        # from the rotation axis.
        #
        # It is therefore the same along a given line
        # of sight.
        #

        if use_rotation:

            v_rot_los = (
                rotation_gradient
                * x_rot
                * np.sin(inclination)
            )

        else:

            v_rot_los = np.zeros_like(X)


        # ----------------------------------------------------
        # Total local LOS velocity
        # ----------------------------------------------------

        v_local = (
            v_sys
            + v_rot_los
            + v_exp_los
        )


        # ----------------------------------------------------
        # Add Gaussian line profile
        # ----------------------------------------------------

        for ichan, v in enumerate(velocity):

            profile = np.exp(
                -0.5
                * (
                    (
                        v
                        - v_local
                    )
                    / sigma_v
                )**2
            )

            contribution = (
                envelope_emissivity
                * profile
                * envelope_mask
                * dz
            )

            cube[ichan] += contribution


# ============================================================
# NORMALIZATION
# ============================================================

#
# The envelope is integrated along z and therefore naturally
# has units proportional to path length.
#
# Normalize only if desired.
#

max_cube = cube.max()

if max_cube > 0:

    cube /= max_cube


# ============================================================
# VELOCITY FIELD AT STELLAR SURFACE
# ============================================================

#
# This is useful for visualizing the rotational component.
#

velocity_surface = np.full_like(
    X,
    np.nan
)

velocity_surface[stellar_disc] = (
    v_sys
    + v_rot_surface[stellar_disc]
)


# ============================================================
# EXPANSION VELOCITY MAP
# ============================================================

#
# For visualization, calculate the maximum blueshift
# encountered along each line of sight through the
# envelope.
#

max_expansion_los = np.zeros_like(X)


if use_expansion:

    for z in z_values:

        r = np.sqrt(
            X**2
            + Y**2
            + z**2
        )

        envelope_mask = (
            (r >= stellar_radius_mas)
            &
            (r <= shell_radius_mas)
        )

        v_exp = np.zeros_like(r)

        v_exp[envelope_mask] = (
            vexp_max
            * (
                shell_radius_mas
                - r[envelope_mask]
            )
            / (
                shell_radius_mas
                - stellar_radius_mas
            )
        )

        valid = (
            envelope_mask
            & (r > 0)
        )

        v_exp_los = np.zeros_like(r)

        v_exp_los[valid] = (
            v_exp[valid]
            * z
            / r[valid]
        )

        max_expansion_los = np.minimum(
            max_expansion_los,
            v_exp_los
        )


# ============================================================
# PLOT 1:
# ROTATION MAP
# ============================================================

plt.figure(
    figsize=(7, 6)
)

rotation_plot = np.where(
    stellar_disc,
    v_rot_surface,
    np.nan
)

im = plt.imshow(
    rotation_plot,
    origin="lower",
    extent=[
        x.min(),
        x.max(),
        y.min(),
        y.max()
    ],
    cmap="RdBu_r"
)

plt.colorbar(
    im,
    label="Rotation LOS velocity [km/s]"
)

plt.xlabel(
    "ΔRA [mas]"
)

plt.ylabel(
    "ΔDec [mas]"
)

plt.title(
    "Solid-body rotational velocity"
)

plt.axis("equal")

plt.tight_layout()


# ============================================================
# PLOT 2:
# EXPANSION MAP
# ============================================================

plt.figure(
    figsize=(7, 6)
)

im = plt.imshow(
    max_expansion_los,
    origin="lower",
    extent=[
        x.min(),
        x.max(),
        y.min(),
        y.max()
    ],
    cmap="Blues_r"
)

plt.colorbar(
    im,
    label="Expansion LOS velocity [km/s]"
)

plt.xlabel(
    "ΔRA [mas]"
)

plt.ylabel(
    "ΔDec [mas]"
)

plt.title(
    "Maximum blueshift from expanding envelope"
)

plt.axis("equal")

plt.tight_layout()


# ============================================================
# CHANNEL MAPS
# ============================================================

fig, axes = plt.subplots(
    2,
    3,
    figsize=(12, 7)
)


channels_to_plot = [
    0,
    nchan // 4,
    nchan // 2,
    3 * nchan // 4,
    nchan - 1,
    np.argmax(
        cube.sum(axis=(1, 2))
    )
]


for ax, channel in zip(
    axes.flat,
    channels_to_plot
):

    im = ax.imshow(
        cube[channel],
        origin="lower",
        extent=[
            x.min(),
            x.max(),
            y.min(),
            y.max()
        ],
        cmap="inferno",
        vmin=0,
        vmax=1
    )

    ax.set_title(
        f"v = {velocity[channel]:.1f} km/s"
    )

    ax.set_xlabel(
        "ΔRA [mas]"
    )

    ax.set_ylabel(
        "ΔDec [mas]"
    )

    ax.set_aspect(
        "equal"
    )


fig.colorbar(
    im,
    ax=axes.ravel().tolist(),
    label="Normalized intensity"
)

plt.tight_layout()


# ============================================================
# INTEGRATED SPECTRUM
# ============================================================

integrated_spectrum = cube.sum(
    axis=(1, 2)
)


plt.figure(
    figsize=(8, 4)
)

plt.plot(
    velocity,
    integrated_spectrum,
    color="black",
    linewidth=2
)

plt.axvline(
    v_sys,
    color="gray",
    linestyle="--",
    alpha=0.7
)

plt.xlabel(
    "Velocity [km/s]"
)

plt.ylabel(
    "Integrated intensity [arbitrary units]"
)

plt.title(
    "Integrated spectral line"
)

plt.grid(
    alpha=0.3
)

plt.tight_layout()


# ============================================================
# FITS HEADER
# ============================================================

pixel_scale_deg = (
    (x[1] - x[0])
    / (3600.0 * 1000.0)
)


header = fits.Header()


# ------------------------------------------------------------
# Spatial axes
# ------------------------------------------------------------

header["CTYPE1"] = "RA---TAN"
header["CTYPE2"] = "DEC--TAN"

header["CUNIT1"] = "deg"
header["CUNIT2"] = "deg"

header["CDELT1"] = -pixel_scale_deg
header["CDELT2"] = pixel_scale_deg

header["CRPIX1"] = (
    nx / 2.0 + 1
)

header["CRPIX2"] = (
    ny / 2.0 + 1
)

header["CRVAL1"] = 0.0
header["CRVAL2"] = 0.0


# ------------------------------------------------------------
# Spectral axis
# ------------------------------------------------------------

header["CTYPE3"] = "VELO-LSR"
header["CUNIT3"] = "km/s"

header["CDELT3"] = channel_width

header["CRPIX3"] = 1

header["CRVAL3"] = velocity[0]


# ------------------------------------------------------------
# Model information
# ------------------------------------------------------------

header["BUNIT"] = "normalized"

header["RSTAR"] = (
    stellar_radius_mas,
    "Stellar radius [mas]"
)

header["RSHELL"] = (
    shell_radius_mas,
    "Outer envelope radius [mas]"
)

header["VEXPMAX"] = (
    vexp_max,
    "Expansion velocity at Rstar [km/s]"
)

header["ROTGRAD"] = (
    rotation_gradient,
    "Solid-body rotation gradient [km/s/mas]"
)

header["INCL"] = (
    inclination_deg,
    "Rotation-axis inclination [deg]"
)

header["PA"] = (
    position_angle_deg,
    "Rotation-axis PA [deg]"
)

header["VSYS"] = (
    v_sys,
    "Systemic velocity [km/s]"
)

header["SIGMAV"] = (
    sigma_v,
    "Intrinsic Gaussian sigma [km/s]"
)

header["EXPONLY"] = (
    "FRONT",
    "Envelope: front hemisphere only"
)


# ============================================================
# WRITE FITS
# ============================================================

hdu = fits.PrimaryHDU(
    data=cube.astype(np.float32),
    header=header
)

hdu.writeto(
    output_fits,
    overwrite=True
)


# ============================================================
# SUMMARY
# ============================================================

print()
print("=" * 60)
print(" Synthetic ALMA spectral-line cube")
print(" Rotating stellar surface + expanding envelope")
print("=" * 60)

print(
    f"Stellar diameter       : "
    f"{stellar_diameter_mas:.1f} mas"
)

print(
    f"Envelope radius        : "
    f"{shell_radius_mas:.1f} mas"
)

print(
    f"Expansion velocity     : "
    f"{vexp_max:.1f} km/s at Rstar"
)

print(
    f"Rotation gradient      : "
    f"{rotation_gradient:.2f} km/s/mas"
)

print(
    f"Inclination             : "
    f"{inclination_deg:.1f} deg"
)

print(
    f"Position angle          : "
    f"{position_angle_deg:.1f} deg"
)

print(
    f"Velocity range          : "
    f"{v_min:.1f} to {v_max:.1f} km/s"
)

print(
    f"Channel width           : "
    f"{channel_width:.2f} km/s"
)

print(
    f"Gaussian sigma           : "
    f"{sigma_v:.2f} km/s"
)

print(
    f"Spatial pixels           : "
    f"{nx} x {ny}"
)

print(
    f"LOS integration points  : "
    f"{nz}"
)

print(
    f"FITS output             : "
    f"{output_fits}"
)

print("=" * 60)

plt.show()
