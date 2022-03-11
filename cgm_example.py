# pyXSIM/SOXS LEM Simulations of CGM

# This is a quick script to show how to create a mock LEM X-ray observation 
# of a galaxy from the CGM of simulated galaxies using yt, pyXSIM, and SOXS.

# First we import the necessary modules:

import soxs
soxs.set_mission_config("lem")

import yt
import pyxsim
from regions import RectangleSkyRegion
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy import wcs
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# Set the "sim" variable to one of the following possibilities:
# 
# * TNG
# * FIRE
# 
# This will set the `hot_gas` particle filter, the list of metals, 
# and anything else simulation-specific.
# 
# The `hot_gas` filter is designed to make a new yt particle type "hot_gas" 
# which masks out any particles which should not be X-ray emitting, based 
# on density, temperature, cooling time, star formation rate, etc. This 
# will vary between simulations as they have different fields.

# Set the simulation you want to try
sim = "FIRE"

if sim == "TNG":
    fn = "cutout_65.hdf5"
    metals = ["He_fraction", "C_fraction", "N_fraction", "O_fraction", 
              "Ne_fraction", "Mg_fraction", "Si_fraction", "Fe_fraction"]
    def hot_gas(pfilter, data):
        pfilter1 = data[pfilter.filtered_type, "temperature"] > 3.0e5
        pfilter2 = data["PartType0", "StarFormationRate"] == 0.0
        pfilter3 = data[pfilter.filtered_type, "cooling_time"] > 50*3.15576e+13
        return pfilter1 & pfilter2 & pfilter3
    def _cooling_time(field, data):
        nH = data["gas", "H_nuclei_density"]
        cool_rate = -ds.arr(data["PartType0","GFM_CoolingRate"].d, "erg*cm**3/s")*nH*nH
        e = 1.5*data["PartType0","pressure"]
        return e/cool_rate
elif sim == "FIRE":
    fn = "snapshot_600.0.hdf5"
    metals = ["He_metallicity", "C_metallicity", "N_metallicity", "O_metallicity", 
              "Ne_metallicity", "Mg_metallicity", "Si_metallicity", "S_metallicity", 
              "Ca_metallicity", "Fe_metallicity"]
    def hot_gas(pfilter, data):
        pfilter1 = data[pfilter.filtered_type, "temperature"] > 3.0e5
        pfilter2 = data[pfilter.filtered_type, "density"] < 5e-25
        pfilter3 = data["PartType0", "StarFormationRate"] == 0.0
        return pfilter1 & pfilter2 & pfilter3


# We now tell yt what the particle filter is.

yt.add_particle_filter("hot_gas", function=hot_gas,
                       filtered_type='gas', requires=["temperature","density"])


# Load the dataset and add the particle filter.

ds = yt.load(fn)
if sim in ["TNG"]:  
    ds.add_field(("gas","cooling_time"), _cooling_time, units="s", 
                 sampling_type="local", force_override=True)
ds.add_particle_filter("hot_gas")

# For this example, the center will be set to the potential minimum, but you 
# could set the variable "c" to a different center if you wanted.

_, c = ds.find_min(("PartType0","Potential"))

# Now create a box centered on the galaxy, 1 Mpc in width, which will be 
# used to draw the cells or particles to make the photons.

width = ds.quan(1.0, "Mpc")
le = c - 0.5*width
re = c + 0.5*width
box = ds.box(le, re)

# Optionally, we can make some simple plots from the simulation itself. NOTE that
# for the FIRE simulation this will take a long time, so you may just want to skip
# it. If so, leave make_sim_plots = False. 

make_sim_plots = False

# Just to get a sense of what things look like, make a projection plot 
# of the total gas density:

if make_sim_plots:
    prj = yt.ProjectionPlot(ds, "z", ("gas","density"), width=(0.4, "Mpc"), center=c,
                            data_source=box)
    prj.set_zlim(("gas","density"), 1.0e-5, 1.0)
    prj.save(f"gas_density_{sim}.png")

# And show the same plot of the "hot" gas density:

if make_sim_plots:
    prj = yt.ProjectionPlot(ds, "z", ("hot_gas","density"), width=(0.4, "Mpc"), center=c, 
                            data_source=box)
    prj.set_zlim(("hot_gas","density"), 1.0e-5, 1.0e-3)
    prj.save(f"hot_gas_density_{sim}.png")

# And the "hot" gas temperature, weighted by the emission measure:

if make_sim_plots:
    prj = yt.ProjectionPlot(ds, "z", ("hot_gas","kT"), width=(0.4, "Mpc"), center=c, 
                            weight_field=("hot_gas", "emission_measure"), data_source=box)
    prj.set_log(("hot_gas","kT"), False)
    prj.set_cmap(("hot_gas","kT"), "turbo")
    prj.save(f"hot_gas_temperature_{sim}.png")

# And we also make a phase plot of the hot gas density vs. temperature, 
# showing the mass of gas at different phases.

if make_sim_plots:
    pp = yt.PhasePlot(box, ("hot_gas","density"), ("hot_gas","kT"), ("hot_gas","mass"), 
                    weight_field=None)
    pp.set_unit(("hot_gas","mass"), "Msun")
    pp.save(f"hot_gas_phase_{sim}.png")


# Now, in order to make the mock observation, we have to set up a 
# ThermalSourceModel in pyXSIM, telling it which fields to use from 
# the dataset, and the min, max, and binning of the spectrum:

# Variable elements--other elements not in this list 
# are assumed to have a single metallicity value
var_elem = {elem.split("_")[0]: ("hot_gas", elem) for elem in metals}

emin = 0.05 # The minimum energy to generate in keV
emax = 3.0 # The maximum energy to generate in keV
nbins = 6000 # The number of energy bins between emin and emax
kT_max = 4.0 # The max gas temperature to use
source_model = pyxsim.ThermalSourceModel(
    "apec", emin, emax, nbins, ("hot_gas","metallicity"),
    temperature_field=("hot_gas","temperature"),
    emission_measure_field=("hot_gas", "emission_measure"),
    var_elem=var_elem, kT_max=kT_max
)

# Now we specify fiducial values of the exposure time, collecting area, and redshift 
# of the object. NOTE that the "area" here is not the effective area of the 
# telescope--this is just a parameter that we need to decide how many sample photons 
# to make. This number should be bigger than the peak of the telescope+instrument 
# effective area curve.

exp_time = (1000., "ks") # exposure time
area = (5000.0, "cm**2") # collecting area
redshift = 0.01 # the cosmological redshift of the source
sky_center = (45.0, 30.0) # in degrees, relatively unimportant

# Now we actually make the photons. They are saved to a file which can be used later. 

n_photons, n_cells = pyxsim.make_photons(f"{sim}_photons", box, redshift, area, 
                                         exp_time, source_model)

# Now we project the photons from the file we saved, also including foreground 
# galactic absorption. We project along the "z"-axis of the simulation box, but 
# any on-axis or off-axis direction could be chosen.

n_events = pyxsim.project_photons(f"{sim}_photons", f"{sim}_events", "z", sky_center,
                                  absorb_model="wabs", nH=0.018)

# Open the file containing the projected events, and convert it to SIMPUT format. 
# The SIMPUT format will be used by SOXS's instrument simulation below. 

events = pyxsim.EventList(f"{sim}_events.h5")
events.write_to_simput(sim, overwrite=True)

# We make a separate background events file and spectrum--this will be used to 
# compare to the files without the source. For now, we turn off the CXB point-source 
# component.

soxs.make_background_file("lem_bkgnd_evt.fits", exp_time, "lem", sky_center, 
                          overwrite=True, ptsrc_bkgnd=False)
soxs.write_spectrum("lem_bkgnd_evt.fits", "lem_bkgnd_evt.pi", overwrite=True)

# We can now use this SIMPUT catalog to make a mock LEM observation. We also 
# ignore the point-source background here for now.

soxs.instrument_simulator(f"{sim}_simput.fits", f"{sim}_evt.fits", 
                          exp_time, "lem", sky_center, overwrite=True,
                          ptsrc_bkgnd=False)

# This produces an event file. This can be used to make images, spectra, and radial 
# profiles with FTOOLS, CIAO, etc., in much the same way as real X-ray observations. 
# Below, we'll show some examples using Python, but you can do your own thing. 

# First we select some lines we want to plot:

bands = {
    "O_VII": (0.552, 0.558),
    "O_VIII": (0.643, 0.65),
    "Fe_XVII": (0.715, 0.725)
}

# We can show an image of the source in these narrow bands, where the LEM field of 
# view is shown with the green square:

center_sky = SkyCoord(45, 30, unit='deg', frame='fk5')
region_sky = RectangleSkyRegion(center=center_sky, width=32*u.arcmin, height=32*u.arcmin)
fig, axes = plt.subplots(ncols=3, nrows=1, figsize=(21,7))
width = 0.6 # degrees
for i, (band, ebounds) in enumerate(bands.items()):
    ax = axes[i]
    soxs.write_image(f"{sim}_evt.fits", f"{sim}_{band}_img.fits", emin=ebounds[0], 
                     emax=ebounds[1], overwrite=True)
    with fits.open(f"{sim}_{band}_img.fits") as f:
        w = wcs.WCS(header=f[0].header)
        dx_pix = width / w.wcs.cdelt[0]
        dy_pix = width / w.wcs.cdelt[1]
        im = ax.imshow(f[0].data[::,::-1], norm=LogNorm(), 
                       cmap="afmhot", origin='lower')
        ax.set_xlim(w.wcs.crpix[0] - 0.5*dx_pix, w.wcs.crpix[0] + 0.5*dx_pix)
        ax.set_ylim(w.wcs.crpix[1] - 0.5*dy_pix, w.wcs.crpix[1] + 0.5*dy_pix)
        ax.set_facecolor("k")
        ax.add_artist(region_sky.to_pixel(w).as_artist(color='limegreen', lw=2))
        ax.set_xticks([])
        ax.set_yticks([])
        ax.text(0.1, 0.83, band.replace("_", " "), color='white', fontsize=30, 
                transform=ax.transAxes)
fig.subplots_adjust(wspace=0.05)
fig.savefig(f"{sim}_images.png")
    
# We can also make a spectrum of the whole source:

soxs.write_spectrum(f"{sim}_evt.fits", f"{sim}_evt.pi", overwrite=True)

# and plot the spectrum for the three different redshifted lines:

plt.rc("font", size=18)
plt.rc("axes", linewidth=2)
plt.rc("xtick.major", width=2, size=6)
plt.rc("ytick.major", width=2, size=6)
plt.rc("ytick.minor", width=2, size=3)
fig, axes = plt.subplots(ncols=3, nrows=1, figsize=(21,7))
for i, (band, ebounds) in enumerate(bands.items()):
    ax = axes[i]
    soxs.plot_spectrum("lem_bkgnd_evt.pi",xmin=ebounds[0]-0.02, xmax=ebounds[1]+0.02, 
                       xscale="linear", fig=fig, ax=ax, yscale='log', label="Background Only")
    soxs.plot_spectrum(f"{sim}_evt.pi", xmin=ebounds[0]-0.02, xmax=ebounds[1]+0.02, 
                       xscale="linear", fig=fig, ax=ax, yscale='log', label="All Emission")
    ax.set_ylim(0.5, 350)
    if i > 0:
        ax.set_ylabel('')
    if i == 2:
        ax.legend()
    ax.text(0.1, 0.7, band.replace("_", " "), color='black', fontsize=30, 
            transform=ax.transAxes)
fig.savefig(f"{sim}_spectra.png")

