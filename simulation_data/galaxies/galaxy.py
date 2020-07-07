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

#input params: redshift==numerical-val; id==int(must exist in range, pre-check); plot=="True" or "False"
#output: (plot=False): (SFH, BE): SFH=StellarFormationHistory, BE=BinEdges
#output: (plot=True): plt.hist (SFH, BE)

def get_stellar_formation_history(redshift, id, plot=False, binwidth=0.05): 

    import h5py
    params = {'stars':'GFM_StellarFormationTime,GFM_InitialMass'}
#if exists, load it. If does not exist, download it
    from pathlib import Path

    if Path('cutout_'+str(id)+'.hdf5').is_file():
        saved_filename = 'cutout_'+str(id)+'.hdf5'
        #print ("File exist") #extract if exists
        with h5py.File(saved_filename, mode='r') as f: #store as h5py file
            starFormationTime = f['PartType4']['GFM_StellarFormationTime'][:]
            starInitialMass = f['PartType4']['GFM_InitialMass'][:]
    else:
        #print ("File not exist") #download and extract if does not exist
        url = "http://www.tng-project.org/api/TNG100-1/snapshots/z=" + str(redshift) + "/subhalos/" + str(id)
        sub = get(url) # get json response of subhalo properties
        saved_filename = get(url + "/cutout.hdf5",params) # get and save HDF5 cutout file
        with h5py.File(saved_filename, mode='r') as f: #store as h5py file
            starFormationTime = f['PartType4']['GFM_StellarFormationTime'][:]
            starInitialMass = f['PartType4']['GFM_InitialMass'][:]

    z_starFormationTime = 1/starFormationTime -1
    Gyr_starFormationTime = cosmo.age(z_starFormationTime).value
    M_Odot_starInitialMass = starInitialMass*1e10/h 
    Gyr_redshift = cosmo.age(2.0).value
    HistWeights = M_Odot_starInitialMass/(binwidth*1e9)
    LookbackTime = Gyr_redshift - Gyr_starFormationTime

    if plot==False:
        return np.histogram(LookbackTime, bins=np.arange(min(Gyr_redshift - Gyr_starFormationTime), max(Gyr_redshift - Gyr_starFormationTime) + binwidth, binwidth), weights=(HistWeights))

    else:
        SFH, BE = np.histogram(LookbackTime, bins=np.arange(min(Gyr_redshift - Gyr_starFormationTime), max(Gyr_redshift - Gyr_starFormationTime) + binwidth, binwidth), weights=(HistWeights))
        bincenters = [(BE[i]+BE[i+1])/2. for i in range(len(BE)-1)]
        plt.figure(figsize=(10,7))
        plt.step(bincenters, SFH, color = 'b')
        plt.title('Histogram for Lookback Times for id = ' + str(id))
        plt.xlim(0, Gyr_redshift)
        plt.ylim(0, )
        plt.xlabel("Lookback Time (Gyr)")
        plt.ylabel("$M_\odot$/yr")
        plt.show()
        return plt.show()

#input params: z=redshift (num val); subhalo_id==int(must exist in range, pre-check)
#dependent on SFH_get(plot=False) output
#output: MeanStellarAge, in units of Gyr

def mean_stellar_formation_time(z, subhalo_id):
    SFH, BE = get_stellar_formation_history(redshift = z, id = subhalo_id)
    MeanStellarAge = np.sum([(SFH[i]*1e9)*BE[i] for i in range(len(SFH))])/sum(SFH*1e9)

    return MeanStellarAge
#opt return: print("Mean Stellar Age for subhalo with id=" +str(subhalo_id)+ " at redshift " + str(z) + " is " + str(MeanStellarAge) + " Gyr")

#input params: z=redshift(num val); subhalo_id==int(must exist in range, pre-check); timescale=num val in range(LT) in units of Gyr
#dependent on get_SFH(plot=False) output
#output: TimeAvg_SFR in units of $M_\odot$/yr

def timeaverage_stellar_formation_rate(z, subhalo_id, timescale):
    SFH, BE = get_stellar_formation_history(redshift = z, id = subhalo_id, plot=False)
    if timescale<BE[0]:
        TimeAvg_SFR = SFH[0]
    else:
        timescale_indices = np.where(BE<=timescale)
        TimeAvg_SFR = np.sum([SFH[i] for i in timescale_indices]) / len(timescale_indices[0])
        #NOTE: floored bin value, change to SFH[i+1] for ceiling
    return TimeAvg_SFR

#input params: z=redshift (num val); subhalo_id==int(must exist in range, pre-check)
#dependent on SFH_get(plot=False) output
#output: MedianStellarAge = a Bin Edge in units of Gyr

def median_stellar_formation_time(z, subhalo_id):
    SFH, BE = get_stellar_formation_history(redshift = z, id = subhalo_id)
    total = [(SFH[i]*1e9)*BE[i] for i in range(len(SFH))] #arr
    tot_sum = (0.5 * np.sum(total))
    cumsum_median = np.cumsum(total) < tot_sum #filter_arr
    newarr = np.cumsum(total)[cumsum_median]
    return BE[len(newarr)-1] #units: Gyr in Lookback time