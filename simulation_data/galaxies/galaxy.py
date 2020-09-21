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



def get_galaxy_particle_data(id, redshift, populate_dict=False):
    """
    input params: id==int(must exist in range, pre-check); z==redshift (num val); populate_dict=boolean, default value ==False
    preconditions: depends on get(), a function, imported from simulation_data
    output: if populate_dict == False: saves file 'cutout_'+str(id)+'_data.hdf5'
    output: if populate_dict == True: saves file 'cutout_'+str(id)+'_data.hdf5'; returns dictionary with data (6 keys)
    output dictionary keys: 'relative_x_coordinates' : #units: physical kpc
                            'relative_y_coordinates' : #units: physical kpc
                            'relative_z_coordinates' : #units: physical kpc
                            'LookbackTime' : #units: Gyr
                            'stellar_initial_masses' : #units: solar mass
                            'stellar_metallicities' : #units: solar metallicity  
                            'initial_x_coordinates' : #units: physical kpc
                            'initial_y_coordinates' : #units: physical kpc
                            'initial_z_coordinates' : #units: physical kpc
                            'initial_x_velocities' : #units: km/s
                            'initial_y_velocities' : #units: km/s
                            'initial_z_velocities' : #units: km/s
                            'u_band' : #units: Vega magnitudes
                            'v_band' : #units: Vega magnitudes
                            'i_band' : #units: AB magnitudes
    """
    stellar_data = {}
    import h5py
    params = {'stars':'Coordinates,GFM_StellarFormationTime,GFM_InitialMass,GFM_Metallicity,BirthPos,BirthVel,GFM_StellarPhotometrics'}
    #looping through fields makes the code longer for the data manipulation section
    
    import os
    import urllib
    
    from pathlib import Path
    new_saved_filename = os.path.join('redshift_'+str(redshift)+'_data', 'cutout_'+str(id)+'_redshift_'+str(redshift)+'_data.hdf5')

    if Path('redshift_'+str(redshift)+'_data\cutout_'+str(id)+'_redshift_'+str(redshift)+'_data.hdf5').is_file():
        pass
    else:
        url = "http://www.tng-project.org/api/TNG100-1/snapshots/z=" + str(redshift) + "/subhalos/" + str(id)
        sub = get(url) # get json response of subhalo properties
        saved_filename = get(url + "/cutout.hdf5",params) # get and save HDF5 cutout file
        with h5py.File(saved_filename, mode='r') as f: #read from h5py file
            x_init = f['PartType4']['BirthPos'][:,0]
            y_init = f['PartType4']['BirthPos'][:,1]
            z_init = f['PartType4']['BirthPos'][:,2]
            vx_init = f['PartType4']['BirthVel'][:,0]
            vy_init = f['PartType4']['BirthVel'][:,1]
            vz_init = f['PartType4']['BirthVel'][:,2]
            dx = f['PartType4']['Coordinates'][:,0] - sub['pos_x']
            dy = f['PartType4']['Coordinates'][:,1] - sub['pos_y']
            dz = f['PartType4']['Coordinates'][:,2] - sub['pos_z']
            starFormationTime = f['PartType4']['GFM_StellarFormationTime'][:]
            starInitialMass = f['PartType4']['GFM_InitialMass'][:]
            starMetallicity = f['PartType4']['GFM_Metallicity'][:]
            U = f['PartType4']['GFM_StellarPhotometrics'][:,0] #Vega magnitudes
            V = f['PartType4']['GFM_StellarPhotometrics'][:,2] #Vega magnitudes
            I = f['PartType4']['GFM_StellarPhotometrics'][:,6] #AB magnitudes


        #selecting star particles only
        x_init = x_init[starFormationTime>0] #ckpc/h
        y_init = y_init[starFormationTime>0]
        z_init = z_init[starFormationTime>0]
        vx_init = vx_init[starFormationTime>0]
        vy_init = vy_init[starFormationTime>0]
        vz_init = vz_init[starFormationTime>0]
        dx = dx[starFormationTime>0]
        dy = dy[starFormationTime>0]
        dz = dz[starFormationTime>0]
        starInitialMass = starInitialMass[starFormationTime>0]
        starMetallicity = starMetallicity[starFormationTime>0]
        starFormationTime = starFormationTime[starFormationTime>0]
        U = U[starFormationTime>0] #Vega magnitudes
        V = V[starFormationTime>0] #Vega magnitudes
        I = I[starFormationTime>0] #AB magnitudes

        scale_factor = a = 1.0 / (1 + redshift)
        inv_sqrt_a = a**(-1/2)
        #unit conversions
        x_init = x_init*a/h #units: physical kpc
        y_init = y_init*a/h #units: physical kpc
        z_init = z_init*a/h #units: physical kpc
        vx_init = vx_init*inv_sqrt_a #units: km/s
        vy_init = vy_init*inv_sqrt_a #units: km/s
        vz_init = vz_init*inv_sqrt_a #units: km/s
        dx = dx*a/h #units: physical kpc
        dy = dy*a/h #units: physical kpc
        dz = dz*a/h #units: physical kpc   
        starFormationTime = 1/starFormationTime - 1 #units:scale factor
        starFormationTime = cosmo.age(starFormationTime).value #units:Gyr
        starInitialMass = starInitialMass*1e10/h #units:solar mass
        Gyr_redshift = cosmo.age(redshift).value #units:Gyr
        LookbackTime = Gyr_redshift - starFormationTime #units:Gyr
        starMetallicity = starMetallicity / 0.0127 #units: solar metallicity

        #delete pre-existing file since this is faster than replacing each field
        import os
        os.remove('cutout_'+str(id)+'.hdf5')
        #create new file with same filename
        new_saved_filename = os.path.join('redshift_'+str(redshift)+'_data', 'cutout_'+str(id)+'_redshift_'+str(redshift)+'_data.hdf5')
        #new_saved_filename = 'cutout_'+str(id)+'_redshift_'+str(redshift)+'_data.hdf5'
        with h5py.File(new_saved_filename, 'w') as h5f:
            #writing data
            d1 = h5f.create_dataset('relative_x_coordinates', data = dx)
            d2 = h5f.create_dataset('relative_y_coordinates', data = dy)
            d3 = h5f.create_dataset('relative_z_coordinates', data = dz)
            d4 = h5f.create_dataset('LookbackTime', data = LookbackTime)
            d5 = h5f.create_dataset('stellar_initial_masses', data = starInitialMass)
            d6 = h5f.create_dataset('stellar_metallicities', data = starMetallicity)
            d7 = h5f.create_dataset('initial_x_coordinates', data = x_init)
            d8 = h5f.create_dataset('initial_y_coordinates', data = y_init)
            d9 = h5f.create_dataset('initial_z_coordinates', data = z_init)
            d10 = h5f.create_dataset('initial_x_velocities', data = vx_init)
            d11 = h5f.create_dataset('initial_y_velocities', data = vy_init)
            d12 = h5f.create_dataset('initial_z_velocities', data = vz_init)
            d13 = h5f.create_dataset('u_band', data = U) #Vega magnitudes
            d14 = h5f.create_dataset('v_band', data = V) #Vega magnitudes
            d15 = h5f.create_dataset('i_band', data = I) #Vega magnitudes
        #close file
        #h5f.close()
    
    with h5py.File(new_saved_filename, 'r+') as h5f_open:
        dx = h5f_open['relative_x_coordinates'][:]
        dy = h5f_open['relative_y_coordinates'][:]
        dz = h5f_open['relative_z_coordinates'][:]
        LookbackTime = h5f_open['LookbackTime'][:]
        starInitialMass = h5f_open['stellar_initial_masses'][:]
        starMetallicity = h5f_open['stellar_metallicities'][:]
        x_init = h5f_open['initial_x_coordinates'][:]
        y_init = h5f_open['initial_y_coordinates'][:]
        z_init = h5f_open['initial_z_coordinates'][:]
        vx_init = h5f_open['initial_x_velocities'][:]
        vy_init = h5f_open['initial_y_velocities'][:]
        vz_init = h5f_open['initial_z_velocities'][:]
        U = h5f_open['u_band'][:]
        V = h5f_open['v_band'][:]
        I = h5f_open['i_band'][:]
        
    stellar_data = {
                    'relative_x_coordinates' : dx, #units: physical kpc
                    'relative_y_coordinates' : dy, #units: physical kpc
                    'relative_z_coordinates' : dz, #units: physical kpc
                    'LookbackTime' : LookbackTime, #units:Gyr
                    'stellar_initial_masses' : starInitialMass, #units:solar mass
                    'stellar_metallicities' : starMetallicity, #units:solar metallicity
                    'initial_x_coordinates' : x_init, #units: physical kpc
                    'initial_y_coordinates' : y_init, #units: physical kpc
                    'initial_z_coordinates' : z_init, #units: physical kpc
                    'initial_x_velocities' : vx_init, #units: km/s
                    'initial_y_velocities' : vy_init, #units: km/s
                    'initial_z_velocities' : vz_init, #units: km/s
                    'u_band' : U, #units: Vega magnitudes
                    'v_band' : V, #units: Vega magnitudes
                    'i_band' : I, #units: AB magnitudes
                   }
    if populate_dict==False:
        return
    else:
        return stellar_data

    

def get_merger_tree(id, redshift):
    """
    input params: id==int(must exist in range, pre-check); redshift=redshift (num val)
    preconditions: uses get() to access subhalo catalog
    output: saves the merger tree for the subhalo
    """
    import os
    import shutil
    from pathlib import Path
    new_saved_filename = os.path.join('redshift_'+str(redshift)+'_data', 'sublink_'+str(id)+'_redshift_'+str(redshift)+'.hdf5')
    if Path(new_saved_filename).is_file():
        pass
    else:
        url = "http://www.tng-project.org/api/TNG100-1/snapshots/z=" + str(redshift) + "/subhalos/" + str(id)
        sub = get(url) # get json response of subhalo properties
        get(sub['trees']['sublink'])
        shutil.copy2('sublink_'+str(id)+'.hdf5', new_saved_filename) 
        os.remove('sublink_'+str(id)+'.hdf5')
    return 



def get_star_formation_history(id, redshift, plot=False, binwidth=0.05): 
    """
    input params: id==int(must exist in range, pre-check); redshift==numerical-val; plot=="True" or "False"
    preconditions: depends on get_galaxy_particle_data(id=id , redshift=redshift, populate_dict=True) output
    output: (plot=False): (SFH, BE): SFH=StellarFormationHistory, units: $M_\odot$/yr;  BE=BinEdges, units: Gyr
    output: (plot=True): plt.hist (SFH, BE)
    """
    stellar_data = get_galaxy_particle_data(id=id , redshift=redshift, populate_dict=True)
    HistWeights = stellar_data['stellar_initial_masses']
    #HistWeights = stellar_data['stellar_initial_masses']/(binwidth*1e9) #units: logMsol/yr
    LookbackTime = stellar_data['LookbackTime']
    SFH, BE = np.histogram(LookbackTime, bins=np.arange(0, max(LookbackTime), binwidth), weights=HistWeights, density = True)
    #SFH, BE = np.histogram(LookbackTime, bins=np.arange(0, max(LookbackTime), binwidth), weights=HistWeights)
    bincenters = [(BE[i]+BE[i+1])/2. for i in range(len(BE)-1)]
    if plot==False:
        return bincenters, SFH
    else:     
        plt.figure(figsize=(10,7))
        plt.step(bincenters, SFH, color = 'b')
        plt.title('Histogram for Lookback Times for id = ' + str(id))
        plt.xlim(0, )
        plt.ylim(0, )
        plt.xlabel("Lookback Time (Gyr)")
        plt.ylabel("$M_\odot$/yr")
        plt.show()
        return plt.show()

    

def mean_stellar_age(id, redshift):
    """
    input params: id==int(must exist in range, pre-check); redshift==redshift (num val)
    preconditions: depends on get_galaxy_particle_data(id=id , redshift=redshift, populate_dict=True) output
    output: mean stellar age, units: Gyr
    """
    stellar_data = get_galaxy_particle_data(id=id , redshift=redshift, populate_dict=True)
    LookbackTime = stellar_data['LookbackTime']    
    return np.mean(LookbackTime)
#opt return: print("Mean Stellar Age for subhalo with id=" +str(id)+ " at redshift " + str(redshift) + " is " + str(MeanStellarAge) + " Gyr")



def timeaverage_stellar_formation_rate(id, redshift, timescale, start=0, binwidth=0.05):
    """
    input params: redshift=redshift(num val); id==int(must exist in range, pre-check); timescale=num val in range(LT) in units of Gyr; start=num val in range(LT) in units of Gyr, default value == 0
    preconditions: depends on get_stellar_formation_history(redshift = redshift, id = id, plot=False) output; first array BC=bincenters & second array SFH=star formation histories
    output: average stellar formation rate over a specified timescale, units: $M_\odot$ /yr
    """
    BC, SFH = get_star_formation_history(redshift = redshift, id = id, plot=False, binwidth=binwidth)
    timescale_indices = np.where((np.array(BC)<=start+timescale+BC[0])&(np.array(BC)>=start)) 
    TimeAvg_SFR = np.sum([SFH[i] for i in timescale_indices]) / len(timescale_indices[0])
        #NOTE: ceiling bin value by BC[0] to accommodate edge case of timescale=start (including 0)
    return TimeAvg_SFR



def current_star_formation_rate(id, redshift):
    """
    input params: id==int(must exist in range, pre-check); redshift=redshift (num val)
    preconditions: uses get() to access subhalo catalog
    output: current SFR, units: $M_\odot$ /yr
    """
    url = "http://www.tng-project.org/api/TNG100-1/snapshots/z=" + str(redshift) + "/subhalos/" + str(id)
    sub = get(url) # get json response of subhalo properties
    return sub['sfr']




def median_stellar_age(id, redshift):
    """
    input params: id==int(must exist in range, pre-check); redshift=redshift (num val)
    preconditions: depends on get_galaxy_particle_data(id=id , redshift=redshift, populate_dict=True) output
    output: median stellar age, units: Gyr
    """
    stellar_data = get_galaxy_particle_data(id=id , redshift=redshift, populate_dict=True)
    LookbackTime = stellar_data['LookbackTime']
    return np.median(LookbackTime) #units: Gyr in Lookback time



def mean_stellar_metallicity(id, redshift):
    """
    input params: id==int(must exist in range, pre-check); redshift=redshift (num val)
    preconditions: depends on get_galaxy_particle_data(id=id , redshift=redshift, populate_dict=True) output
    output: mean stellar metallicity, units: solar metallicity
    """
    stellar_data = get_galaxy_particle_data(id=id , redshift=redshift, populate_dict=True)
    stellar_metallicities = stellar_data['stellar_metallicities']    
    return np.mean(stellar_metallicities)



def mean_stellar_mass(id, redshift):
    """
    input params: id==int(must exist in range, pre-check); redshift=redshift (num val)
    preconditions: depends on get_galaxy_particle_data(id=id , redshift=redshift, populate_dict=True) output
    output: mean stellar mass, units: solar masses
    """
    stellar_data = get_galaxy_particle_data(id=id , redshift=redshift, populate_dict=True)
    stellar_mass = stellar_data['stellar_initial_masses']    
    return np.mean(stellar_mass)



def total_stellar_mass(id, redshift):
    """
    input params: id==int(must exist in range, pre-check); redshift=redshift (num val)
    preconditions: uses get() to access subhalo catalog
    output: total stellar mass, units: log10 solar masses
    """
    url = "http://www.tng-project.org/api/TNG100-1/snapshots/z=" + str(redshift) + "/subhalos/" + str(id)
    sub = get(url) # get json response of subhalo properties
    return np.log10(sub['mass_stars']*1e10/h)



def halfmass_rad_stars(id, redshift):
    """
    input params: id==int(must exist in range, pre-check); redshift=redshift (num val)
    preconditions: depends on get_galaxy_particle_data(id=id , redshift=redshift, populate_dict=True) output
    output: half-mass radius of all stellar particles, units: physical kpc
    """
    scale_factor = a = 1.0 / (1 + redshift)
    url = "http://www.tng-project.org/api/TNG100-1/snapshots/z=" + str(redshift) + "/subhalos/" + str(id)
    sub = get(url) # get json response of subhalo properties
    return sub['halfmassrad_stars']*a/h



def halflight_rad_stars(id, redshift, band):
    """
    input params: redshift=redshift(num val); id==int(must exist in range, pre-check); band==['U'(Vega magnitude), 
                                                                                              'V'(Vega magnitude), 
                                                                                              'I'(AB magnitude),
                                                                                              'M'(solM) == (test)]
    preconditions: depends on get_galaxy_particle_data(id=id , redshift=redshift, populate_dict=True) output
                   band must be == 'U'/'V'/'I'/'M'==mass-test
    output: half-light radius of a galaxy for a given band (U/V/I) in pkpc or half-mass-rad
    """
    stellar_data = get_galaxy_particle_data(id=id , redshift=redshift, populate_dict=True)
    dx = stellar_data['relative_x_coordinates']
    dy = stellar_data['relative_y_coordinates']
    dz = stellar_data['relative_z_coordinates']
    R = (dx**2 + dy**2 + dz**2)**(1/2)#units: physical kpc
    
    if band=='U':
        mag = stellar_data['u_band']
        flux = 10**(-0.4*mag) #flux: flux = 10**(-0.4*mag)
        zipped_lists = zip(R, flux)
        sorted_pairs = sorted(zipped_lists)

        tuples = zip(*sorted_pairs)
        R_sort, band_sort = [list(tuple) for tuple in  tuples]

        band_indices = np.where(np.cumsum(np.array(band_sort))>=0.5*np.sum(np.array(band_sort)))
        halflight_rad = max(np.array(R_sort)[i] for i in band_indices)
    
    elif band=='V':
        mag = stellar_data['v_band']
        flux = 10**(-0.4*mag) #flux
        zipped_lists = zip(R, flux)
        sorted_pairs = sorted(zipped_lists)

        tuples = zip(*sorted_pairs)
        R_sort, band_sort = [list(tuple) for tuple in  tuples]

        band_indices = np.where(np.cumsum(np.array(band_sort))>=0.5*np.sum(np.array(band_sort)))
        halflight_rad = max(np.array(R_sort)[i] for i in band_indices)
    
    elif band=='I':
        mag = stellar_data['i_band']
        flux = 10**(-0.4*mag) #flux
        zipped_lists = zip(R, flux)
        sorted_pairs = sorted(zipped_lists)

        tuples = zip(*sorted_pairs)
        R_sort, band_sort = [list(tuple) for tuple in  tuples]

        band_indices = np.where(np.cumsum(np.array(band_sort))>=0.5*np.sum(np.array(band_sort)))
        halflight_rad = max(np.array(R_sort)[i] for i in band_indices)
    
    elif band=='M':
        mass = stellar_data['stellar_initial_masses']
        zipped_lists = zip(R, mass)
        sorted_pairs = sorted(zipped_lists)

        tuples = zip(*sorted_pairs)
        R_sort, band_sort = [list(tuple) for tuple in  tuples]

        band_indices = np.where(np.cumsum(np.array(band_sort))>=0.5*np.sum(np.array(band_sort)))
        halflight_rad = max(np.array(R_sort)[i] for i in band_indices)

    return min(halflight_rad)



def age_profile(id, redshift, n_bins=20, scatter=False):
    """
    input params: id==int(must exist in range, pre-check); redshift=redshift (num val); n_bins==int(num of bins for percentile-count stellar particle partition, default value = 20)
    preconditions: depends on get_galaxy_particle_data(id=id , redshift=redshift, populate_dict=True) output
    output: statistic = median ages of binned stellar particles (units: Gyr); 
            percentile cutoffs for radial distances =  radial_percentiles[1:]/R_e = percentile bin-edges for normalized radial distance of stellar particle from subhalo center (unitless); 
            R_e = effective (median) radius of stellar particles in subhalo (units: physical kpc)
    """
    stellar_data = get_galaxy_particle_data(id=id , redshift=redshift, populate_dict=True)
    LookbackTime = stellar_data['LookbackTime']
    dx = stellar_data['relative_x_coordinates']
    dy = stellar_data['relative_y_coordinates']
    dz = stellar_data['relative_z_coordinates']
    metallicity = stellar_data['stellar_metallicities']
    R = (dx**2 + dy**2 + dz**2)**(1/2)#units: physical kpc
    
    radial_percentiles = np.zeros(n_bins + 1) #N+1 for N percentiles 
    for i in range(1, (n_bins+1)):
        radial_percentiles[i] = np.percentile(R, (100/n_bins)*i) 
    R_e = np.nanmedian(R)
    statistic, bin_edges, bin_number = scipy.stats.binned_statistic(R, LookbackTime, 'median', bins=radial_percentiles)
    
    if scatter==False:
        return statistic, radial_percentiles[:-1]/R_e, R_e, R/R_e, LookbackTime, np.log10(metallicity)
    else:
        plt.figure(figsize=(10,7)) # 10 is width, 7 is height
        plt.scatter(R/R_e, LookbackTime, c=np.log10(metallicity), s=0.5, alpha=0.7)#c=np.log10(metallicity)
        plt.plot(np.array(radial_percentiles[1:]/R_e)[4:-4], np.array(statistic)[4:-4], c='black')
        plt.xlim(1e-2, )
        plt.ylim(1e-1, )
        plt.grid()
        plt.colorbar(boundaries=np.linspace(-3.1,1.1,100), label='Metallicities of Stars ($\log_{10}$ $Z_\odot$)')
        plt.title('Radial Distance vs Stellar Ages (log/log scale) with Binned Age Trend for id='+str(id))
        plt.xlabel('Normalized Radial Distance (R/$R_e$)')
        plt.ylabel('Stellar Ages in Lookback Times(Gyr)')
        plt.xscale('log')
        plt.yscale('log')
        plt.show()
        return plt.show()


    
    
def metallicity_profile(id, redshift, n_bins=20, scatter=False):
    """
    input params: id==int(must exist in range, pre-check); redshift=redshift (num val); n_bins==int(num of bins for percentile-count stellar particle partition, default value = 20)
    preconditions: depends on get_galaxy_particle_data(id=id , redshift=redshift, populate_dict=True) output
    output: statistic = median metallicities of binned stellar particles (units: log10zsol); 
            percentile cutoffs for radial distances =  radial_percentiles[1:]/R_e = percentile bin-edges for normalized radial distance of stellar particle from subhalo center (unitless); 
            R_e = effective (median) radius of stellar particles in subhalo (units: physical kpc)
    """
    stellar_data = get_galaxy_particle_data(id=id , redshift=redshift, populate_dict=True)
    LookbackTime = stellar_data['LookbackTime']
    dx = stellar_data['relative_x_coordinates']
    dy = stellar_data['relative_y_coordinates']
    dz = stellar_data['relative_z_coordinates']
    metallicity = stellar_data['stellar_metallicities']
    R = (dx**2 + dy**2 + dz**2)**(1/2)#units: physical kpc
    
    radial_percentiles = np.zeros(n_bins + 1) #N+1 for N percentiles 
    for i in range(1, (n_bins+1)):
        radial_percentiles[i] = np.percentile(R, (100/n_bins)*i) 
    R_e = np.nanmedian(R)
    statistic, bin_edges, bin_number = scipy.stats.binned_statistic(R, metallicity, 'median', bins=radial_percentiles)
    
    if scatter==False:
        return statistic, radial_percentiles[:-1]/R_e, R_e 
    else:
        plt.figure(figsize=(10,7)) # 10 is width, 7 is height
        plt.scatter(R/R_e, metallicity, s=2, alpha=0.05)#c=np.log10(LookbackTime)
        plt.plot(np.array(radial_percentiles[1:]/R_e)[4:-4], np.array(statistic)[4:-4], c='black')
        plt.xlim(1e-2, )
        plt.ylim(1e-1, )
        #plt.colorbar(label='Age of Stars (Gyr)')
        plt.title('Normalized Radial Distance vs Stellar Metallicities (log/log scale) with Binned Metallicity Trend for id='+str(id))
        plt.xlabel('Normalized Radial Distance (R/$R_e$)')
        plt.ylabel('Stellar Metallicity($\log_{10}$ $Z_\odot$)')
        plt.xscale('log')
        plt.yscale('log')
        plt.show()
        return plt.show()

    
    
    
def particle_age_profile(id, redshift, n_bins=70, binwidth = 0.10):
    """
    input params: id==int(must exist in range, pre-check); redshift=redshift (num val); n_bins==int(num of bins for percentile-count stellar particle partition, default value = 70); binwidth==float(width of histy bins)
    preconditions: depends on get_galaxy_particle_data(id=id , redshift=redshift, populate_dict=True) output
    output: plot of scatter ages + median age profile + y-hist of age trends
    """
    
    def colorbar(mappable):
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        import matplotlib.pyplot as plt
        last_axes = plt.gca()
        ax = mappable.axes
        fig = ax.figure
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = fig.colorbar(mappable, cax=cax)
        cbar.set_boundaries=np.linspace(-3.1,1.1,100)
        cbar.set_label('Metallicity ($\log_{10}$ $Z_\odot$)',size=12)
        plt.sca(last_axes)
        return cbar

    # Define colorbar
    colors = [(0.3, 0.76, 1), (0, 0, 0), (1, 0.3, 0.3)]  # Bu ->  Blk  -> R
    #n_bins = [3, 6, 10, 100]  # Discretizes the interpolation into bins
    cmap_name = 'my_list'

    # Create the colormap
    from matplotlib.colors import LinearSegmentedColormap
    cm = LinearSegmentedColormap.from_list(
            cmap_name, colors, N=100) 


    # populating dictionary
    stellar_data = get_galaxy_particle_data(id=id, redshift=redshift, populate_dict=True)

    #getting data from arrays
    LookbackTime = stellar_data['LookbackTime']
    dx = stellar_data['relative_x_coordinates']
    dy = stellar_data['relative_y_coordinates']
    dz = stellar_data['relative_z_coordinates']
    metallicity = stellar_data['stellar_metallicities']
    R = (dx**2 + dy**2 + dz**2)**(1/2)

    
    radial_percentiles = np.zeros(n_bins + 1) #N+1 for N percentiles 
    for i in range(1, (n_bins+1)):
        radial_percentiles[i] = np.percentile(R, (100/n_bins)*i) 
    R_e = np.nanmedian(R)
    statistic, bin_edges, bin_number = scipy.stats.binned_statistic(R, LookbackTime, 'median', bins=radial_percentiles)

    x = R/R_e #np.log10(R/R_e)
    y = LookbackTime #np.log10(LookbackTime)

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.09


    rect_scatter = [left, bottom, width, height]
    #rect_histx = [left, bottom + height + spacing, width, 0.2]
    rect_histy = [left + width + spacing, bottom, 0.2, height]

    # start with a rectangular Figure
    plt.figure(figsize=(10, 10))

    ax_scatter = plt.axes(rect_scatter)
    ax_scatter.tick_params(direction='in', top=True, right=True)

    #ax_histx = plt.axes(rect_histx)
    #ax_histx.tick_params(direction='in', labelbottom=True)
    ax_histy = plt.axes(rect_histy)
    ax_histy.tick_params(direction='in', labelleft=True)


    # the scatter plot:
    im = ax_scatter.scatter(x, y, s=0.5, alpha=0.7, c=np.log10(metallicity), cmap=cm, vmin=-3, vmax=1)
    ax_scatter.plot(np.array(radial_percentiles[1:]/R_e)[4:-4], np.array(statistic)[4:-4], c='k')
    #ax_scatter.colorbar(boundaries=np.linspace(-3.1,1.1,100), label='Metallicities of Stars ($\log_{10}$ $Z_\odot$)')
    ax_scatter.set_xlabel("Radial Distance ($R_e$) | id = " + str(id))
    ax_scatter.set_ylabel("Stellar Age in Lookback Time (Gyr)")
    ax_scatter.set_xscale('log')
    #ax_scatter.set_yscale('log')

    # Colorbar
    colorbar(im)

    # now determine nice limits by hand:
    
    #lim = np.ceil(np.abs([x, y]).max() / binwidth) * binwidth
    ax_scatter.set_xlim((min(x), max(x)))
    ax_scatter.set_ylim((min(y), max(y)))

    x_bins = np.arange(min(x), max(x) + binwidth, binwidth) #np.linspace, #np.logspace
    y_bins = np.arange(min(y), max(y) + binwidth, binwidth)
    #ax_histx.hist(x, bins=x_bins, histtype = 'step') #add xlabels
    ax_histy.hist(y, bins=y_bins, histtype = 'step', orientation='horizontal')

    #ax_histx.set_xlim(ax_scatter.get_xlim())
    ax_histy.set_ylim(ax_scatter.get_ylim())

    #plt.tight_layout()
    return plt.show()    