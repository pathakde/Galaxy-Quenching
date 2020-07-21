import numpy as np
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

#imported requests
import requests
#import get()
from simulation_data import get

from .galaxy import mean_stellar_age, timeaverage_stellar_formation_rate, median_stellar_age, mean_stellar_metallicity, age_profile

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
    def calc_mean_stellar_age(self ):
        #create and populate array for mean SFT
        ids = self.ids
        means = np.zeros(len(ids))
        for i, id in enumerate(ids): 
            means[i] = mean_stellar_age(redshift = self.redshift, id = id)
        #save file
        np.savetxt('z='+str(self.redshift)+'_Mean_SFT', means)
        mean_SFT = np.loadtxt('z='+str(self.redshift)+'_Mean_SFT', dtype=float)
        return mean_SFT
    
    
    def get_mean_stellar_age(self): #can parameterize for max mass if needed
        import pathlib
        file = pathlib.Path('z='+str(self.redshift)+'_Mean_SFT')
        if file.exists ():
            MeanSFT = np.loadtxt('z='+str(self.redshift)+'_Mean_SFT', dtype=float) #rename pre-existing files before parameterizing further
            return MeanSFT
        else:
            return self.calc_mean_stellar_age()
    
    #time avg SFR
    def calc_timeaverage_stellar_formation_rate(self, calc_timescale):
        #create and populate array for mean SFT
        ids = self.ids
        time_averages = np.zeros(len(ids))
        for i, id in enumerate(ids): 
            time_averages[i] = timeaverage_stellar_formation_rate(redshift = self.redshift, id = id, timescale = calc_timescale)
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
            return self.calc_timeaverage_stellar_formation_rate(calc_timescale=timescale)
            
    
    #current SFR
    def calc_current_stellar_formation_rate(self):
        #create and populate array for mean SFT
        ids = self.ids
        current_SFRs = np.zeros(len(ids))
        for i, id in enumerate(ids): 
            current_SFRs[i] = timeaverage_stellar_formation_rate(redshift = self.redshift, id = id, timescale = 0.01)
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
            return self.calc_current_stellar_formation_rate()
        
        
    #SFR ratio
    def calc_stellar_formation_rate_ratio(self, calc_timescale):
        #create and populate array for mean SFT
        ids = self.ids
        ratios = np.zeros(len(ids))
        for i, id in enumerate(ids): 
            ratios[i] = timeaverage_stellar_formation_rate(redshift = self.redshift, id = id, timescale = 0.01) / timeaverage_stellar_formation_rate(redshift = self.redshift, id = id, timescale = calc_timescale)
        #save file
        np.savetxt( 'z='+str(self.redshift)+ '_SFR_Ratio_'+ str(calc_timescale) +'Gyr', ratios)
        ratios = np.loadtxt('z='+ str(self.redshift) +'_SFR_Ratio_'+ str(calc_timescale)+ 'Gyr', dtype=float)
        return ratios
    
        
    def get_stellar_formation_rate_ratio(self, timescale): #can parameterize for max mass if needed
        import pathlib
        file = pathlib.Path('z='+str(self.redshift)+'_SFR_Ratio_'+str(timescale)+'Gyr')
        if file.exists ():
            ratios = np.loadtxt('z='+str(self.redshift)+ '_SFR_Ratio_' + str(timescale) +'Gyr', dtype=float) #rename pre-existing files before parameterizing further
            return ratios
        else:
            return self.calc_stellar_formation_rate_ratio(calc_timescale=timescale)
    
    
    #median SFT
    def calc_median_stellar_age(self):
        #create and populate array for mean SFT
        ids = self.ids
        MedianSFT = np.zeros(len(ids))
        for i, id in enumerate(ids):
            MedianSFT[i] = median_stellar_age(redshift = self.redshift, id = id)
        #save file
        np.savetxt('z='+ str(self.redshift) +'_Mean_SFT', MedianSFT)
        median_SFT = np.loadtxt('z='+ str(self.redshift) +'_Mean_SFT', dtype=float)
        return median_SFT
    
    
    def get_median_stellar_age(self):
        import pathlib
        file = pathlib.Path('z='+ str(self.redshift) +'_Mean_SFT')
        if file.exists ():
            median_SFT = np.loadtxt('z='+ str(self.redshift) +'_Mean_SFT', dtype=float) #rename pre-existing files before parameterizing further
            return median_SFT
        else:
            return self.calc_median_stellar_age()
    
    
    #median SFT
    def calc_effective_radius(self):
        #create and populate array for mean SFT
        ids = self.ids
        r_effective = np.zeros(len(ids))
        for i, id in enumerate(ids):
            r_effective[i] = age_profile(id=id, redshift=self.redshift, n_bins=20)[2]
        #save file
        np.savetxt('z='+ str(self.redshift) +'_effective_radius', r_effective)
        r_effective = np.loadtxt('z='+ str(self.redshift) +'_effective_radius', dtype=float)
        return r_effective
    
    
    def get_effective_radius(self):
        import pathlib
        file = pathlib.Path('z='+ str(self.redshift) +'_effective_radius')
        if file.exists ():
            r_effective = np.loadtxt('z='+ str(self.redshift) +'_effective_radius', dtype=float) #rename pre-existing files before parameterizing further
            return r_effective
        else:
            return self.calc_effective_radius()
    
    
    #mean stellar metallicity
    def calc_mean_stellar_metallicity(self):
        #create and populate array for mean SFT
        ids = self.ids
        mean_metallicity = np.zeros(len(ids))
        for i, id in enumerate(ids):
            mean_metallicity[i] = mean_stellar_metallicity(id=id, redshift=self.redshift)
        #save file
        np.savetxt('z='+ str(self.redshift) +'_mean_metallicity', mean_metallicity)
        mean_metallicity = np.loadtxt('z='+ str(self.redshift) +'_mean_metallicity', dtype=float)
        return mean_metallicity
    
    
    def get_mean_stellar_metallicity(self):
        import pathlib
        file = pathlib.Path('z='+ str(self.redshift) +'_mean_metallicity')
        if file.exists ():
            mean_metallicity = np.loadtxt('z='+ str(self.redshift) +'_mean_metallicity', dtype=float) #rename pre-existing files before parameterizing further
            return mean_metallicity
        else:
            return self.calc_mean_stellar_metallicity()
    