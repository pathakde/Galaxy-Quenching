import numpy as np
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

#imported requests
import requests
#import get()
from simulation_data import get

from .galaxy import mean_stellar_age, timeaverage_stellar_formation_rate, median_stellar_age, mean_stellar_metallicity, age_profile, mean_stellar_mass, total_stellar_mass, halfmass_rad_stars, halflight_rad_stars

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
            self.ids = np.array(self.ids, dtype=np.int32)
        return self.ids
 


    def download_galaxy_population_data(self, redshift):
        redshift = self.redshift
        galaxy_population_data = {}
        import h5py
        from pathlib import Path
        if Path('galaxy_population_data_'+str(self.redshift)+'.hdf5').is_file():
            pass
        else:
            with h5py.File('galaxy_population_data_'+str(self.redshift)+'.hdf5', 'w') as f:
                #writing data
                d1 = f.create_dataset('ids', data = self.select_galaxies(redshift=redshift, mass_min=10.5, mass_max=12))
                d2 = f.create_dataset('mean_age', data = self.get_mean_stellar_age())
                d3 = f.create_dataset('median_age', data = self.get_median_stellar_age())
                d4 = f.create_dataset('current_SFR', data = self.get_timeaverage_stellar_formation_rate(timescale=0))
                d5 = f.create_dataset('average_SFR_0.5', data = self.get_timeaverage_stellar_formation_rate(timescale=0.5))
                d6 = f.create_dataset('average_SFR_1', data = self.get_timeaverage_stellar_formation_rate(timescale=1))
                d7 = f.create_dataset('average_SFR_2', data = self.get_timeaverage_stellar_formation_rate(timescale=2))
                d8 = f.create_dataset('SFR_ratio_1', data = self.get_stellar_formation_rate_ratio(timescale=1))
                d9 = f.create_dataset('SFR_ratio_0.5', data = self.get_stellar_formation_rate_ratio(timescale=0.5))
                d10 = f.create_dataset('effective_radius', data = self.get_effective_radius())
                d11 = f.create_dataset('halfmass_radius', data = self.get_halfmass_rad_stars())
                d12 = f.create_dataset('mean_metallicity', data = self.get_mean_stellar_metallicity())
                d13 = f.create_dataset('total_mass', data = self.get_total_stellar_mass())
                d14 = f.create_dataset('halflight_radius_U', data = self.get_halflight_rad_stars(band='U'))
                d15 = f.create_dataset('halflight_radius_V', data = self.get_halflight_rad_stars(band='V'))
                d16 = f.create_dataset('halflight_radius_I', data = self.get_halflight_rad_stars(band='I'))
                d17 = f.create_dataset('halfmass_radius_calculated', data = self.get_halflight_rad_stars(band='M'))
                
        with h5py.File('galaxy_population_data_'+str(self.redshift)+'.hdf5', 'r') as f:
            ids = f['ids'][:]
            mean_age = f['mean_age'][:]
            median_age = f['median_age'][:]
            current_SFR = f['current_SFR'][:] 
            average_SFR_0_5 = f['average_SFR_0.5'][:]
            average_SFR_1 = f['average_SFR_1'][:]
            average_SFR_2 = f['average_SFR_2'][:]
            SFR_ratio_1 = f['SFR_ratio_1'][:]
            SFR_ratio_0_5 = f['SFR_ratio_0.5'][:]
            effective_radius = f['effective_radius'][:]
            halfmass_radius = f['halfmass_radius'][:]
            mean_metallicity = f['mean_metallicity'][:]
            total_mass = f['total_mass'][:]
            halflight_radius_U = f['halflight_radius_U'][:]
            halflight_radius_V = f['halflight_radius_V'][:]
            halflight_radius_I = f['halflight_radius_I'][:]
            halfmass_radius_calculated = f['halfmass_radius_calculated'][:]
        galaxy_population_data = {
                                    'ids': ids,
                                    'mean_age': mean_age,
                                    'median_age': median_age,
                                    'current_SFR': current_SFR,
                                    'average_SFR_0_5': average_SFR_0_5,
                                    'average_SFR_1': average_SFR_1,
                                    'average_SFR_2': average_SFR_2,
                                    'SFR_ratio_1': SFR_ratio_1,
                                    'SFR_ratio_0_5': SFR_ratio_0_5,
                                    'effective_radius': effective_radius,
                                    'halfmass_radius': halfmass_radius,
                                    'mean_metallicity': mean_metallicity,
                                    'total_mass': total_mass,
                                    'halflight_radius_U': halflight_radius_U,
                                    'halflight_radius_V': halflight_radius_V,
                                    'halflight_radius_I': halflight_radius_I,
                                    'halfmass_radius_calculated': halfmass_radius_calculated
                                 }
        return galaxy_population_data

    
    
    def post_starburst_selection(self, redshift, current_SFR, average_SFR, timescale, count=True):
        self.redshift = redshift
        galaxy_population_data = self.download_galaxy_population_data(redshift = redshift)
        post_starburst_ids = galaxy_population_data['ids'][(galaxy_population_data['current_SFR'] <= current_SFR) & (galaxy_population_data['average_SFR_'+str(timescale)] >= average_SFR)]
        if count==True:
            return len(post_starburst_ids)
        else:
            return post_starburst_ids
    
    
    
    #mean SFT
    def calc_mean_stellar_age(self ):
        ids = self.ids
        means = np.zeros(len(ids))
        for i, id in enumerate(ids): 
            means[i] = mean_stellar_age(redshift = self.redshift, id = id)
        #save file
        np.savetxt('z='+str(self.redshift)+'_Mean_SFT', means)
        mean_SFT = np.loadtxt('z='+str(self.redshift)+'_Mean_SFT', dtype=float)
        return mean_SFT
    
    
    def get_mean_stellar_age(self): 
        import pathlib
        file = pathlib.Path('z='+str(self.redshift)+'_Mean_SFT')
        if file.exists ():
            MeanSFT = np.loadtxt('z='+str(self.redshift)+'_Mean_SFT', dtype=float) #rename pre-existing files before parameterizing further
            return MeanSFT
        else:
            return self.calc_mean_stellar_age()
    
    #time avg SFR
    def calc_timeaverage_stellar_formation_rate(self, calc_timescale, calc_start=0):
        ids = self.ids
        time_averages = np.zeros(len(ids))
        for i, id in enumerate(ids): 
            time_averages[i] = timeaverage_stellar_formation_rate(redshift = self.redshift, id = id, timescale = calc_timescale, start=calc_start)
        #save file
        np.savetxt( 'z='+str(self.redshift)+ '_TimeAvg_SFR_'+ str(calc_start) + '_' + str(calc_timescale) +'Gyr', time_averages)
        time_avg_SFT = np.loadtxt('z='+str(self.redshift)+ '_TimeAvg_SFR_'+ str(calc_start) + '_' + str(calc_timescale) +'Gyr', dtype=float)
        return time_avg_SFT
    
        
    def get_timeaverage_stellar_formation_rate(self, timescale, start = 0):
        import pathlib
        file = pathlib.Path('z='+str(self.redshift)+ '_TimeAvg_SFR_'+ str(start) + '_' + str(timescale) +'Gyr')
        if file.exists ():
            time_avg_SFT = np.loadtxt('z='+str(self.redshift)+ '_TimeAvg_SFR_'+ str(start) + '_' + str(timescale) +'Gyr', dtype=float) #rename pre-existing files before parameterizing further
            return time_avg_SFT
        else:
            return self.calc_timeaverage_stellar_formation_rate(calc_timescale=timescale, calc_start=start)
            
    
    #current SFR
    def calc_current_stellar_formation_rate(self):
        ids = self.ids
        current_SFRs = np.zeros(len(ids))
        for i, id in enumerate(ids): 
            current_SFRs[i] = timeaverage_stellar_formation_rate(redshift = self.redshift, id = id, timescale = 0, start = 0)
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
    def calc_stellar_formation_rate_ratio(self, calc_timescale, calc_start=0):
        ids = self.ids
        ratios = np.zeros(len(ids))
        for i, id in enumerate(ids): 
            ratios[i] = timeaverage_stellar_formation_rate(redshift = self.redshift, id = id, timescale = 0, start=0) / timeaverage_stellar_formation_rate(redshift = self.redshift, id = id, timescale = calc_timescale, start=calc_start)
        #save file
        np.savetxt( 'z='+str(self.redshift)+ '_SFR_Ratio_'+ str(calc_start) + '_' + str(calc_timescale) +'Gyr', ratios)
        ratios = np.loadtxt('z='+ str(self.redshift) +'_SFR_Ratio_'+ str(calc_start) + '_' + str(calc_timescale) +'Gyr', dtype=float)
        return ratios
    
        
    def get_stellar_formation_rate_ratio(self, timescale, start=0): #can parameterize for max mass if needed
        import pathlib
        file = pathlib.Path('z='+str(self.redshift)+'_SFR_Ratio_'+ str(start) + '_' + str(timescale) +'Gyr')
        if file.exists ():
            ratios = np.loadtxt('z='+str(self.redshift)+ '_SFR_Ratio_' + str(start) + '_' + str(timescale) +'Gyr', dtype=float) 
            return ratios
        else:
            return self.calc_stellar_formation_rate_ratio(calc_timescale=timescale, calc_start=start)
    
    
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
    
    
    #effective radius
    def calc_effective_radius(self):
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
        #create and populate array for mean metallicity
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
    
    
        #mean stellar mass
    def calc_mean_stellar_mass(self):
        #create and populate array for mean mass
        ids = self.ids
        mean_mass = np.zeros(len(ids))
        for i, id in enumerate(ids):
            mean_mass[i] = mean_stellar_mass(id=id, redshift=self.redshift)
        #save file
        np.savetxt('z='+ str(self.redshift) +'_mean_mass', mean_mass)
        mean_mass = np.loadtxt('z='+ str(self.redshift) +'_mean_mass', dtype=float)
        return mean_mass
    
    
    def get_mean_stellar_mass(self):
        import pathlib
        file = pathlib.Path('z='+ str(self.redshift) +'_mean_mass')
        if file.exists ():
            mean_mass = np.loadtxt('z='+ str(self.redshift) +'_mean_mass', dtype=float) #rename pre-existing files before parameterizing further
            return mean_mass
        else:
            return self.calc_mean_stellar_mass()
        
        
        #total stellar mass
    def calc_total_stellar_mass(self):
        ids = self.ids
        total_mass = np.zeros(len(ids))
        for i, id in enumerate(ids):
            total_mass[i] = total_stellar_mass(id=id, redshift=self.redshift)
        #save file
        np.savetxt('z='+ str(self.redshift) +'_total_mass', total_mass)
        total_mass = np.loadtxt('z='+ str(self.redshift) +'_total_mass', dtype=float)
        return total_mass
    
    
    def get_total_stellar_mass(self):
        import pathlib
        file = pathlib.Path('z='+ str(self.redshift) +'_total_mass')
        if file.exists ():
            total_mass = np.loadtxt('z='+ str(self.redshift) +'_total_mass', dtype=float) #rename pre-existing files before parameterizing further
            return total_mass
        else:
            return self.calc_total_stellar_mass()
        
        
        
        #stellar half mass radius
    def calc_halfmass_rad_stars(self):
        ids = self.ids
        halfmass_rad = np.zeros(len(ids))
        for i, id in enumerate(ids):
            halfmass_rad[i] = halfmass_rad_stars(id=id, redshift=self.redshift)
        #save file
        np.savetxt('z='+ str(self.redshift) +'_halfmass_rad', halfmass_rad)
        halfmass_rad = np.loadtxt('z='+ str(self.redshift) +'_halfmass_rad', dtype=float)
        return halfmass_rad
    
    
    def get_halfmass_rad_stars(self):
        import pathlib
        file = pathlib.Path('z='+ str(self.redshift) +'_halfmass_rad')
        if file.exists ():
            halfmass_rad = np.loadtxt('z='+ str(self.redshift) +'_halfmass_rad', dtype=float) #rename pre-existing files before parameterizing further
            return halfmass_rad
        else:
            return self.calc_halfmass_rad_stars()
        
        

         #halflight_rad_stars
    def calc_halflight_rad_stars(self, calc_band):
        ids = self.ids
        halflight_rad = np.zeros(len(ids))
        for i, id in enumerate(ids):
            halflight_rad[i] = halflight_rad_stars(id=id, redshift=self.redshift, band=calc_band)
        #save file
        np.savetxt('z='+ str(self.redshift) +str(calc_band)+'_halflight_rad', halflight_rad)
        halflight_rad = np.loadtxt('z='+ str(self.redshift) +str(calc_band)+'_halflight_rad', dtype=float)
        return halflight_rad
    
    
    
    def get_halflight_rad_stars(self, band):
        import pathlib
        file = pathlib.Path('z='+ str(self.redshift) +str(band)+'_halflight_rad')
        if file.exists ():
            halflight_rad = np.loadtxt('z='+ str(self.redshift) +str(band)+'_halflight_rad', dtype=float) #rename pre-existing files before parameterizing further
            return halflight_rad
        else:
            return self.calc_halflight_rad_stars(calc_band = band)
        
        
        