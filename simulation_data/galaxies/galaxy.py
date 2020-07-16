import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

from simulation_data import get

from io import StringIO
import io

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from astropy.constants import G, h, k_B
h = 0.6774
cosmo = FlatLambdaCDM(H0= (h * 100) * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3)

import scipy
from scipy import stats



#input params: subhalo_id==int(must exist in range, pre-check); z==redshift (num val); populate_dict=boolean, default value ==False
#preconditions: depends on get(), a function, imported from simulation_data
#output: if populate_dict == False: saves file 'cutout_'+str(subhalo_id)+'.hdf5'
#output: if populate_dict == True: saves file 'cutout_'+str(subhalo_id)+'.hdf5'; returns dictionary with data (6 keys)
#output dictionary keys: 'relative_x_coordinates' : #units: physical kpc
                        #'relative_y_coordinates' : #units: physical kpc
                        #'relative_z_coordinates' : #units: physical kpc
                        #'LookbackTime' : #units: Gyr
                        #'stellar_initial_masses' : #units: solar mass
                        #'stellar_metallicities' : #units: solar metallicity
def get_galaxy_particle_data(subhalo_id, z, populate_dict=False):
    stellar_data = {}
    import h5py
    params = {'stars':'Coordinates,GFM_StellarFormationTime,GFM_InitialMass,GFM_Metallicity'}
    #looping through fields makes the code longer for the data manipulation section
    from pathlib import Path

    if Path('cutout_'+str(subhalo_id)+'.hdf5').is_file():
        pass
    else:
        url = "http://www.tng-project.org/api/TNG100-1/snapshots/z=" + str(z) + "/subhalos/" + str(subhalo_id)
        sub = get(url) # get json response of subhalo properties
        saved_filename = get(url + "/cutout.hdf5",params) # get and save HDF5 cutout file
        with h5py.File(saved_filename, mode='r') as f: #read from h5py file
            dx = f['PartType4']['Coordinates'][:,0] - sub['pos_x']
            dy = f['PartType4']['Coordinates'][:,1] - sub['pos_y']
            dz = f['PartType4']['Coordinates'][:,2] - sub['pos_z']
            starFormationTime = f['PartType4']['GFM_StellarFormationTime'][:]
            starInitialMass = f['PartType4']['GFM_InitialMass'][:]
            starMetallicity = f['PartType4']['GFM_Metallicity'][:]

        #selecting star particles only
        dx = dx[starFormationTime>0]
        dy = dy[starFormationTime>0]
        dz = dz[starFormationTime>0]
        starInitialMass = starInitialMass[starFormationTime>0]
        starMetallicity = starMetallicity[starFormationTime>0]
        starFormationTime = starFormationTime[starFormationTime>0]

        scale_factor = a = 1.0 / (1 + z)
        #unit conversions
        dx = dx*a/h #units: physical kpc
        dy = dy*a/h #units: physical kpc
        dz = dz*a/h #units: physical kpc   
        starFormationTime = 1/starFormationTime - 1 #units:scale factor
        starFormationTime = cosmo.age(starFormationTime).value #units:Gyr
        starInitialMass = starInitialMass*1e10/h #units:solar mass
        Gyr_redshift = cosmo.age(z).value #units:Gyr
        LookbackTime = Gyr_redshift - starFormationTime #units:Gyr
        starMetallicity = starMetallicity / 0.0127 #units:primordial sollar metallicity

        #delete pre-existing file since this is faster than replacing each field
        import os
        os.remove(saved_filename)
        #create new file with same filename
        h5f = h5py.File(saved_filename, 'w')
        #writing data
        h5f.create_dataset('relative_x_coordinates', data = dx)
        h5f.create_dataset('relative_y_coordinates', data = dy)
        h5f.create_dataset('relative_z_coordinates', data = dz)
        h5f.create_dataset('LookbackTime', data = LookbackTime)
        h5f.create_dataset('stellar_initial_masses', data = starInitialMass)
        h5f.create_dataset('stellar_metallicities', data = starMetallicity)
        #close file
        h5f.close()

    saved_filename = 'cutout_'+str(subhalo_id)+'.hdf5'
    h5f_open = h5py.File(saved_filename, 'r')
    dx = h5f_open['relative_x_coordinates'][:]
    dy = h5f_open['relative_y_coordinates'][:]
    dz = h5f_open['relative_z_coordinates'][:]
    LookbackTime = h5f_open['LookbackTime'][:]
    starInitialMass = h5f_open['stellar_initial_masses'][:]
    starMetallicity = h5f_open['stellar_metallicities'][:]
    stellar_data = {
                    'relative_x_coordinates' : dx, #units: physical kpc
                    'relative_y_coordinates' : dy, #units: physical kpc
                    'relative_z_coordinates' : dz, #units: physical kpc
                    'LookbackTime' : LookbackTime, #units:Gyr
                    'stellar_initial_masses' : starInitialMass, #units:solar mass
                    'stellar_metallicities' : starMetallicity #units:primordial stellar metallicity
                   }
    if populate_dict==False:
        return
    else:
        return stellar_data

    
#input params: id==int(must exist in range, pre-check); redshift==numerical-val; plot=="True" or "False"
#preconditions: depends on get_galaxy_particle_data(subhalo_id=id , z=redshift, populate_dict=True) output
#output: (plot=False): (SFH, BE): SFH=StellarFormationHistory, units: $M_\odot$/yr;  BE=BinEdges, units: Gyr
#output: (plot=True): plt.hist (SFH, BE)
def get_stellar_formation_history(id, redshift, plot=False, binwidth=0.05): 
    stellar_data = get_galaxy_particle_data(subhalo_id=id , z=redshift, populate_dict=True)
    HistWeights = stellar_data['stellar_initial_masses']/(binwidth*1e9)
    LookbackTime = stellar_data['LookbackTime']
    if plot==False:
        return np.histogram(LookbackTime, bins=np.arange(min(LookbackTime) - binwidth, max(LookbackTime) + binwidth, binwidth), weights=HistWeights)
    else:
        SFH, BE = np.histogram(LookbackTime, bins=np.arange(min(LookbackTime) - binwidth, max(LookbackTime) + binwidth, binwidth), weights=HistWeights)
        bincenters = [(BE[i]+BE[i+1])/2. for i in range(len(BE)-1)]
        plt.figure(figsize=(10,7))
        plt.step(bincenters, SFH, color = 'b')
        plt.title('Histogram for Lookback Times for id = ' + str(id))
        plt.xlim(0, )
        plt.ylim(0, )
        plt.xlabel("Lookback Time (Gyr)")
        plt.ylabel("$M_\odot$/yr")
        plt.show()
        return plt.show()

    
#input params: subhalo_id==int(must exist in range, pre-check); z==redshift (num val)
#preconditions: depends on get_galaxy_particle_data(subhalo_id=subhalo_id , z=z, populate_dict=True) output
#output: mean stellar age, units: Gyr
def mean_stellar_formation_time(subhalo_id, z):
    stellar_data = get_galaxy_particle_data(subhalo_id=subhalo_id , z=z, populate_dict=True)
    LookbackTime = stellar_data['LookbackTime']    
    return np.mean(LookbackTime)
#opt return: print("Mean Stellar Age for subhalo with id=" +str(subhalo_id)+ " at redshift " + str(z) + " is " + str(MeanStellarAge) + " Gyr")


#input params: z=redshift(num val); subhalo_id==int(must exist in range, pre-check); timescale=num val in range(LT) in units of Gyr
#preconditions: depends on get_stellar_formation_history(redshift = z, id = subhalo_id, plot=False) output; first array SFH
#output: average stellar formation rate over a specified timescale, units: $M_\odot$
def timeaverage_stellar_formation_rate(subhalo_id, z, timescale):
    SFH, BE = get_stellar_formation_history(redshift = z, id = subhalo_id, plot=False)
    if timescale<BE[0]:
        TimeAvg_SFR = SFH[0]
    else:
        timescale_indices = np.where(BE<=timescale)
        TimeAvg_SFR = np.sum([SFH[i] for i in timescale_indices]) / len(timescale_indices[0])
        #NOTE: floored bin value, change to SFH[i+1] for ceiling
    return TimeAvg_SFR


#input params: subhalo_id==int(must exist in range, pre-check); z=redshift (num val)
#preconditions: depends on get_galaxy_particle_data(subhalo_id=subhalo_id , z=z, populate_dict=True) output
#output: median stellar age, units: Gyr
def median_stellar_formation_time(subhalo_id, z):
    stellar_data = get_galaxy_particle_data(subhalo_id=subhalo_id , z=z, populate_dict=True)
    LookbackTime = stellar_data['LookbackTime']
    return np.median(LookbackTime) #units: Gyr in Lookback time


#input params: subhalo_id==int(must exist in range, pre-check); z=redshift (num val)
#preconditions: depends on get_galaxy_particle_data(subhalo_id=subhalo_id , z=z, populate_dict=True) output
#output: mean stellar metallicity, units: solar metallicity
def mean_stellar_metallicity(subhalo_id, z):
    stellar_data = get_galaxy_particle_data(subhalo_id=subhalo_id , z=z, populate_dict=True)
    stellar_metallicities = stellar_data['stellar_metallicities']    
    return np.mean(stellar_metallicities)


#input params: id==int(must exist in range, pre-check); z=redshift (num val); n_bins==int(num of bins for percentile-count stellar particle partition, default value = 20)
#preconditions: depends on get_galaxy_particle_data(subhalo_id=id , z=redshift, populate_dict=True) output
#output: #statistic = median ages of binned stellar particles (units: Gyr); 
         #percentile cutoffs for radial distances =  radial_percentiles[1:]/R_e = percentile bin-edges for normalized radial distance of stellar particle from subhalo center (unitless); 
         #R_e = effective (median) radius of stellar particles in subhalo (units: physical kpc)
def age_profile(id, redshift, n_bins=20):
    stellar_data = get_galaxy_particle_data(subhalo_id=id , z=redshift, populate_dict=True)
    LookbackTime = stellar_data['LookbackTime']
    dx = stellar_data['relative_x_coordinates']
    dy = stellar_data['relative_y_coordinates']
    dz = stellar_data['relative_z_coordinates']
    R = (dx**2 + dy**2 + dz**2)**(1/2)#units: physical kpc
    
    radial_percentiles = np.zeros(n_bins + 1) #N+1 for N percentiles 
    for i in range(1, (n_bins+1)):
        radial_percentiles[i] = np.percentile(R, (100/n_bins)*i) 
    R_e = np.nanmedian(R)
    statistic, bin_edges, bin_number = scipy.stats.binned_statistic(R, LookbackTime, 'median', bins=radial_percentiles)

    return statistic, radial_percentiles[1:]/R_e, R_e