import numpy as np
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

#imported requests
import requests
#import get()
from simulation_data import get

from .galaxy import mean_stellar_formation_time, timeaverage_stellar_formation_rate, median_stellar_formation_time

class GalaxyPopulation():
    
    
    def __init__(self):
        self.ids = []
        self.mass_min = 0
        self.mass_max = 0
        self.redshift = 0
        
        
    #select ids
    def select_galaxies(self, redshift, mass_min, mass_max=12):
        if self.ids == [] or (self.mass_min != mass_min or self.mass_max != mass_max or self.redshift != redshift):
            h = 0.6774
            mass_minimum = 10**mass_min / 1e10 * h
            mass_maximum = 10**mass_max / 1e10 * h
            # form the search_query string by hand for once
            search_query = "?mass_stars__gt=" + str(mass_minimum) + "&mass_stars__lt=" + str(mass_maximum)
            url = "http://www.tng-project.org/api/TNG100-1/snapshots/z=" + str(redshift) + "/subhalos/" + search_query
            subhalos = get(url, {'limit':5000})
            self.mass_min = mass_min
            self.mass_max = mass_max
            self.redshift = redshift
            self.ids = [ subhalos['results'][i]['id'] for i in range(subhalos['count'])]
        return self.ids
    
    
    #mean SFT
    def calc_mean_stellar_formation_time(self ):
        #create and populate array for mean SFT
        ids = self.ids
        means = np.zeros(len(ids))
        for i, id in enumerate(ids): 
            means[i] = mean_stellar_formation_time(z = self.redshift, subhalo_id = id)
        #save file
        np.savetxt('z='+str(self.redshift)+'_Mean_SFT', means)
        mean_SFT = np.loadtxt('z='+str(self.redshift)+'_Mean_SFT', dtype=float)
        return mean_SFT
    
    
    def get_mean_stellar_formation_time(self): #can parameterize for max mass if needed
        import pathlib
        file = pathlib.Path('z='+str(self.redshift)+'_Mean_SFT')
        if file.exists ():
            MeanSFT = np.loadtxt('z='+str(self.redshift)+'_Mean_SFT', dtype=float) #rename pre-existing files before parameterizing further
            return MeanSFT
        else:
            self.calc_mean_stellar_formation_time()
           
    
    #time avg SFR
    def calc_timeaverage_stellar_formation_rate(self, calc_timescale):
        #create and populate array for mean SFT
        ids = self.ids
        time_averages = np.zeros(len(ids))
        for i, id in enumerate(ids): 
            time_averages[i] = timeaverage_stellar_formation_rate(z = self.redshift, subhalo_id = id, timescale = calc_timescale)
        #save file
        np.savetxt( 'z='+str(self.redshift)+ '_TimeAvg_SFR_'+ str(calc_timescale) +'Gyr', time_averages)
        time_avg_SFT = np.loadtxt('z='+ str(self.redshift) +'_TimeAvg_SFR_'+ str(calc_timescale)+ 'Gyr', dtype=float)
        return time_avg_SFT
    
        
    def get_timeaverage_stellar_formation_rate(self, timescale): #can parameterize for max mass if needed
        import pathlib
        file = pathlib.Path('z='+str(self.redshift)+'_TimeAvg_SFR_'+str(timescale)+'Gyr')
        if file.exists ():
            time_avg_SFT = np.loadtxt('z='+str(self.redshift)+ '_TimeAvg_SFR_'+ str(timescale) +'Gyr', dtype=float) #rename pre-existing files before parameterizing further
            return time_avg_SFT
        else:
            self.calc_timeaverage_stellar_formation_rate(calc_timescale=timescale)
    
    
    #current SFR
    def calc_current_stellar_formation_rate(self):
        #create and populate array for mean SFT
        ids = self.ids
        current_SFRs = np.zeros(len(ids))
        for i, id in enumerate(ids): 
            current_SFRs[i] = timeaverage_stellar_formation_rate(z = self.redshift, subhalo_id = id, timescale = 0)
        #save file
        np.savetxt( 'z='+ str(self.redshift) +'_Current_SFR', current_SFRs)
        current_SFR = np.loadtxt('z='+ str(self.redshift) +'_Current_SFR', dtype=float)
        return current_SFR
    
        
    def get_current_stellar_formation_rate(self): #can parameterize for max mass if needed
        import pathlib
        file = pathlib.Path('z='+ str(self.redshift) +'_Current_SFR')
        if file.exists ():
            current_SFR = np.loadtxt('z='+ str(self.redshift) +'_Current_SFR', dtype=float) #rename pre-existing files before parameterizing further
            return current_SFR
        else:
            self.calc_current_stellar_formation_rate()
    
    
    #median SFT
    def calc_median_stellar_formation_time(self):
        #create and populate array for mean SFT
        ids = self.ids
        MedianSFT = np.zeros(len(ids))
        for i, id in enumerate(ids):
            MedianSFT[i] = median_stellar_formation_time(z = self.redshift, subhalo_id = id)
        #save file
        np.savetxt('z='+ str(self.redshift) +'_Mean_SFT', MedianSFT)
        median_SFT = np.loadtxt('z='+ str(self.redshift) +'_Mean_SFT', dtype=float)
        return median_SFT
    
    
    def get_median_stellar_formation_time(self):
        import pathlib
        file = pathlib.Path('z='+ str(self.redshift) +'_Mean_SFT')
        if file.exists ():
            median_SFT = np.loadtxt('z='+ str(self.redshift) +'_Mean_SFT', dtype=float) #rename pre-existing files before parameterizing further
            return median_SFT
        else:
            self.calc_median_stellar_formation_time()
    