import numpy as np
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

#imported requests
import requests
#import get()
from simulation_data import get

from .galaxy import mean_stellar_age, timeaverage_stellar_formation_rate, median_stellar_age, mean_stellar_metallicity, age_profile, mean_stellar_mass, total_stellar_mass, halfmass_rad_stars, halflight_rad_stars, accreted_stellar_fraction, age_gradient, age_gradient_Re, current_star_formation_rate, maximum_mass_jump_ratio, max_merger_mass_ratio, median_age_accreted_stellar_fraction, median_age_half_Re

class GalaxyPopulation():
    
    
    def __init__(self):
        self.ids = []
        self.mass_min = 0
        self.mass_max = 0
        self.redshift = 0
        self.snap = 0
        
        
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
                d18 = f.create_dataset('current_sSFR', data = self.get_current_stellar_formation_rate()/10**(self.get_total_stellar_mass()))
                d19 = f.create_dataset('newbin_current_SFR', data = self.get_timeaverage_stellar_formation_rate(timescale=0, binwidth=0.01))
                d20 = f.create_dataset('age_difference', data = self.get_age_difference())
                d21 = f.create_dataset('accreted_stellar_fraction', data = self.get_accreted_stellar_fraction())
                d22 = f.create_dataset('radial_difference', data = self.get_radial_difference())
                d23 = f.create_dataset('age_gradient', data = self.get_age_difference()/self.get_radial_difference())
                d24 = f.create_dataset('median_age_accreted_stellar_fraction', data = self.get_median_age_accreted_stellar_fraction())
                d25 = f.create_dataset('accreted_stellar_mass_fraction', data = self.get_accreted_stellar_fraction())
                d26 = f.create_dataset('age_gradient_3Re', data = self.get_age_gradient_Re(scale=3))
                d27 = f.create_dataset('median_age_half_Re', data = self.get_median_age_half_Re())
                d28 = f.create_dataset('light_percentile85_U', data = self.get_halflight_rad_stars(band='U', bound=0.85))
                d29 = f.create_dataset('light_percentile85_V', data = self.get_halflight_rad_stars(band='V', bound=0.85))
                d30 = f.create_dataset('light_percentile85_I', data = self.get_halflight_rad_stars(band='I', bound=0.85))
                d31 = f.create_dataset('mass_percentile85', data = my_galaxy_population.get_halflight_rad_stars(band='M', bound=0.85))
                d32 = f.create_dataset('light_percentile95_U', data = my_galaxy_population.get_halflight_rad_stars(band='U', bound=0.95))
                d33 = f.create_dataset('light_percentile95_V', data = my_galaxy_population.get_halflight_rad_stars(band='V', bound=0.95))
                d34 = f.create_dataset('light_percentile95_I', data = my_galaxy_population.get_halflight_rad_stars(band='I', bound=0.95))
                d35 = f.create_dataset('mass_percentile95', data = my_galaxy_population.get_halflight_rad_stars(band='M', bound=0.95))



                
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
            current_sSFR = f['current_sSFR'][:]
            newbin_current_SFR = f['newbin_current_SFR'][:]
            age_difference = f['age_difference'][:]
            radial_difference = f['radial_difference'][:]
            age_gradient = f['age_gradient'][:]
            accreted_stellar_fraction = f['accreted_stellar_fraction'][:]
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
                                    'halfmass_radius_calculated': halfmass_radius_calculated,
                                    'current_sSFR': current_sSFR,
                                    'newbin_current_SFR': newbin_current_SFR,
                                    'age_difference': age_difference,
                                    'radial_difference': radial_difference,
                                    'age_gradient': age_gradient,
                                    'accreted_stellar_fraction': accreted_stellar_fraction
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
    def calc_timeaverage_stellar_formation_rate(self, calc_timescale, calc_start=0, calc_binwidth=0.05):
        ids = self.ids
        time_averages = np.zeros(len(ids))
        for i, id in enumerate(ids): 
            time_averages[i] = timeaverage_stellar_formation_rate(redshift = self.redshift, id = id, timescale = calc_timescale, start=calc_start, binwidth=calc_binwidth)
        #save file
        np.savetxt( 'z='+str(self.redshift)+ '_TimeAvg_SFR_'+ str(calc_start) + '_' + str(calc_timescale) +'Gyr', time_averages)
        time_avg_SFT = np.loadtxt('z='+str(self.redshift)+ '_TimeAvg_SFR_'+ str(calc_start) + '_' + str(calc_timescale) +'Gyr', dtype=float)
        return time_avg_SFT
    
        
    def get_timeaverage_stellar_formation_rate(self, timescale, start = 0, binwidth=0.05):
        import pathlib
        file = pathlib.Path('z='+str(self.redshift)+ '_TimeAvg_SFR_'+ str(start) + '_' + str(timescale) +'Gyr')
        if file.exists ():
            time_avg_SFT = np.loadtxt('z='+str(self.redshift)+ '_TimeAvg_SFR_'+ str(start) + '_' + str(timescale) +'Gyr', dtype=float) #rename pre-existing files before parameterizing further
            return time_avg_SFT
        else:
            return self.calc_timeaverage_stellar_formation_rate(calc_timescale=timescale, calc_start=start, calc_binwidth=binwidth)
            
    
    #current SFR
    def calc_current_stellar_formation_rate(self):
        ids = self.ids
        current_SFRs = np.zeros(len(ids))
        for i, id in enumerate(ids): 
            current_SFRs[i] = current_star_formation_rate(redshift = self.redshift, id = id)
        #save file
        np.savetxt( 'z='+ str(self.redshift) +'_Current_SFR', current_SFRs)
        current_SFR = np.loadtxt('z='+ str(self.redshift) +'_Current_SFR', dtype=float)
        return current_SFR
    
        
    def get_current_stellar_formation_rate(self): 
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
        np.savetxt('z='+ str(self.redshift) +'_Median_SFT', MedianSFT)
        median_SFT = np.loadtxt('z='+ str(self.redshift) +'_Median_SFT', dtype=float)
        return median_SFT
    
    
    def get_median_stellar_age(self):
        import pathlib
        file = pathlib.Path('z='+ str(self.redshift) +'_Median_SFT')
        if file.exists ():
            median_SFT = np.loadtxt('z='+ str(self.redshift) +'_Median_SFT', dtype=float) #rename pre-existing files before parameterizing further
            return median_SFT
        else:
            return self.calc_median_stellar_age()
    
    
    
    #median central 1/2 Re SFT
    def calc_median_age_half_Re(self):
        #create and populate array for mean SFT
        ids = self.ids
        MedianSFT = np.zeros(len(ids))
        for i, id in enumerate(ids):
            MedianSFT[i] = median_age_half_Re(redshift = self.redshift, id = id)
        #save file
        np.savetxt('z='+ str(self.redshift) +'_median_age_half_Re', MedianSFT)
        median_SFT = np.loadtxt('z='+ str(self.redshift) +'_median_age_half_Re', dtype=float)
        return median_SFT
    
    
    
    def get_median_age_half_Re(self):
        import pathlib
        file = pathlib.Path('z='+ str(self.redshift) +'_median_age_half_Re')
        if file.exists ():
            median_SFT = np.loadtxt('z='+ str(self.redshift) +'_median_age_half_Re', dtype=float) #rename pre-existing files before parameterizing further
            return median_SFT
        else:
            return self.calc_median_age_half_Re()
    
    
    
    
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
    def calc_halflight_rad_stars(self, calc_band, calc_bound):
        ids = self.ids
        halflight_rad = np.zeros(len(ids))
        for i, id in enumerate(ids):
            halflight_rad[i] = halflight_rad_stars(id=id, redshift=self.redshift, band=calc_band, bound=calc_bound)
        #save file
        np.savetxt('z='+ str(self.redshift) +str(calc_band)+'_halflight_rad'+str(calc_bound), halflight_rad)
        halflight_rad = np.loadtxt('z='+ str(self.redshift) +str(calc_band)+'_halflight_rad'+str(calc_bound), dtype=float)
        return halflight_rad
    
    
    
    def get_halflight_rad_stars(self, band, bound):
        import pathlib
        file = pathlib.Path('z='+ str(self.redshift) +str(band)+'_halflight_rad'+str(bound))
        if file.exists ():
            halflight_rad = np.loadtxt('z='+ str(self.redshift) +str(band)+'_halflight_rad'+str(bound), dtype=float) #rename pre-existing files before parameterizing further
            return halflight_rad
        else:
            return self.calc_halflight_rad_stars(calc_band = band, calc_bound = bound)
        
        

    #accreted stellar fractions    
    def calc_accreted_stellar_fraction(self):
        ids = self.ids
        accreted_frac = np.zeros(len(ids))
        for i, id in enumerate(ids): 
            accreted_frac[i] = accreted_stellar_fraction(redshift = self.redshift, id = id)
        #save file
        np.savetxt('z='+str(self.redshift)+'_accreted_stellar_fraction', accreted_frac)
        accreted_fraction = np.loadtxt('z='+str(self.redshift)+'_accreted_stellar_fraction', dtype=float)
        return accreted_fraction
    
    
    def get_accreted_stellar_fraction(self): 
        import pathlib
        file = pathlib.Path('z='+str(self.redshift)+'_accreted_stellar_fraction')
        if file.exists ():
            accreted_fraction = np.loadtxt('z='+str(self.redshift)+'_accreted_stellar_fraction', dtype=float) #rename pre-existing files before parameterizing further
            return accreted_fraction
        else:
            return self.calc_accreted_stellar_fraction()
        
        
        
    #median age of accreted stellar fractions    
    def calc_median_age_accreted_stellar_fraction(self):
        ids = self.ids
        accreted_frac = np.zeros(len(ids))
        for i, id in enumerate(ids): 
            accreted_frac[i] = median_age_accreted_stellar_fraction(redshift = self.redshift, id = id)
        #save file
        np.savetxt('z='+str(self.redshift)+'_median_age_accreted_stellar_fraction', accreted_frac)
        accreted_fraction = np.loadtxt('z='+str(self.redshift)+'_median_age_accreted_stellar_fraction', dtype=float)
        return accreted_fraction
    
    
    def get_median_age_accreted_stellar_fraction(self): 
        import pathlib
        file = pathlib.Path('z='+str(self.redshift)+'_median_age_accreted_stellar_fraction')
        if file.exists ():
            accreted_fraction = np.loadtxt('z='+str(self.redshift)+'_median_age_accreted_stellar_fraction', dtype=float) #rename pre-existing files before parameterizing further
            return accreted_fraction
        else:
            return self.calc_median_age_accreted_stellar_fraction()
        
        
        
    #age difference: in Gyr 
    def calc_age_difference(self):
        ids = self.ids
        age_grad = np.zeros(len(ids))
        for i, id in enumerate(ids): 
            age_grad[i] = age_gradient(redshift = self.redshift, id = id)[0]
        #save file
        np.savetxt('z='+str(self.redshift)+'_age_difference', age_grad)
        age_grad = np.loadtxt('z='+str(self.redshift)+'_age_difference', dtype=float)
        return age_grad
    
    
    def get_age_difference(self): 
        import pathlib
        file = pathlib.Path('z='+str(self.redshift)+'_age_difference')
        if file.exists ():
            age_grad = np.loadtxt('z='+str(self.redshift)+'_age_difference', dtype=float) #rename pre-existing files before parameterizing further
            return age_grad
        else:
            return self.calc_age_difference()
        
        
        
    #radial profile difference: in kpc
    def calc_radial_difference(self):
        ids = self.ids
        age_grad = np.zeros(len(ids))
        for i, id in enumerate(ids): 
            age_grad[i] = age_gradient(redshift = self.redshift, id = id)[1]
        #save file
        np.savetxt('z='+str(self.redshift)+'_radial_difference', age_grad)
        age_grad = np.loadtxt('z='+str(self.redshift)+'_radial_difference', dtype=float)
        return age_grad
    
    
    def get_radial_difference(self): 
        import pathlib
        file = pathlib.Path('z='+str(self.redshift)+'_radial_difference')
        if file.exists ():
            age_grad = np.loadtxt('z='+str(self.redshift)+'_radial_difference', dtype=float) #rename pre-existing files before parameterizing further
            return age_grad
        else:
            return self.calc_radial_difference()
        

       
    #radial profile difference: in Gyr 
    def calc_age_gradient_Re(self, calc_scale=3):
        ids = self.ids
        age_grad = np.zeros(len(ids))
        for i, id in enumerate(ids): 
            age_grad[i] = age_gradient_Re(redshift = self.redshift, id = id, scale = calc_scale)[0] / age_gradient_Re(redshift = self.redshift, id = id, scale = calc_scale)[1]
        #save file
        np.savetxt('z='+str(self.redshift)+'_age_gradient_Re'+str(calc_scale), age_grad)
        age_grad = np.loadtxt('z='+str(self.redshift)+'_age_gradient_Re'+str(calc_scale), dtype=float)
        return age_grad
    
    
    def get_age_gradient_Re(self, scale=3): 
        import pathlib
        file = pathlib.Path('z='+str(self.redshift)+'_age_gradient_Re'+str(scale))
        if file.exists ():
            age_grad = np.loadtxt('z='+str(self.redshift)+'_age_gradient_Re'+str(scale), dtype=float) #rename pre-existing files before parameterizing further
            return age_grad
        else:
            return self.calc_age_gradient_Re(calc_scale = scale)
        
        
        
        #max_mass_ratio & snapshot # 
    def calc_maximum_mass_ratio(self, calc_snap, calc_earliest_snap, calc_mass_parameter):
        ids = self.ids
        mass_ratio = np.zeros(len(ids))
        for i, id in enumerate(ids): 
            mass_ratio[i] = maximum_mass_jump_ratio(id, snap=calc_snap, earliest_snap=calc_earliest_snap, mass_parameter=calc_mass_parameter)[0]
        #save file
        np.savetxt('z='+str(self.redshift)+'_maximum_mass_ratio'+calc_mass_parameter, mass_ratio)
        mass_ratio = np.loadtxt('z='+str(self.redshift)+'_maximum_mass_ratio'+calc_mass_parameter, dtype=float)
        return mass_ratio
    
    
    def get_maximum_mass_ratio(self, snap, earliest_snap=6, mass_parameter='stars'): 
        import pathlib
        file = pathlib.Path('z='+str(self.redshift)+'_maximum_mass_ratio'+mass_parameter)
        if file.exists ():
            mass_ratio = np.loadtxt('z='+str(self.redshift)+'_maximum_mass_ratio'+mass_parameter, dtype=float) #rename pre-existing files before parameterizing further
            return mass_ratio
        else:
            return self.calc_maximum_mass_ratio(calc_snap=snap, calc_earliest_snap=earliest_snap, calc_mass_parameter=mass_parameter)
        
        
        
    def calc_maximum_mass_ratio_snap(self, calc_snap, calc_earliest_snap, calc_mass_parameter):
        ids = self.ids
        mass_ratio_snap = np.zeros(len(ids))
        for i, id in enumerate(ids): 
            mass_ratio_snap[i] = maximum_mass_jump_ratio(id, snap=calc_snap, earliest_snap=calc_earliest_snap, mass_parameter=calc_mass_parameter)[1]
        #save file
        np.savetxt('z='+str(self.redshift)+'_maximum_mass_ratio_snap'+calc_mass_parameter, mass_ratio_snap)
        mass_ratio_snap = np.loadtxt('z='+str(self.redshift)+'_maximum_mass_ratio_snap'+calc_mass_parameter, dtype=float)
        return mass_ratio_snap
    
    
    def get_maximum_mass_ratio_snap(self, snap, earliest_snap=6, mass_parameter='stars'): 
        import pathlib
        file = pathlib.Path('z='+str(self.redshift)+'_maximum_mass_ratio_snap'+mass_parameter)
        if file.exists ():
            mass_ratio_snap = np.loadtxt('z='+str(self.redshift)+'_maximum_mass_ratio_snap'+mass_parameter, dtype=float) #rename pre-existing files before parameterizing further
            return mass_ratio_snap
        else:
            return self.calc_maximum_mass_ratio_snap(calc_snap=snap, calc_earliest_snap=earliest_snap, calc_mass_parameter=mass_parameter)
        
        
        
        #max_mass_ratio & snapshot # 
    def calc_max_merger_mass_ratio(self):
        ids = self.ids
        mass_ratio = np.zeros(len(ids))
        for i, id in enumerate(ids): 
            mass_ratio[i] = max_merger_mass_ratio(id=id, redshift=self.redshift)
        #save file
        np.savetxt('z='+str(self.redshift)+'_max_merger_mass_ratio', mass_ratio)
        mass_ratio = np.loadtxt('z='+str(self.redshift)+'_max_merger_mass_ratio', dtype=float)
        return mass_ratio
    
    
    def get_max_merger_mass_ratio(self): 
        import pathlib
        file = pathlib.Path('z='+str(self.redshift)+'_max_merger_mass_ratio')
        if file.exists ():
            mass_ratio = np.loadtxt('z='+str(self.redshift)+'_max_merger_mass_ratio', dtype=float) #rename pre-existing files before parameterizing further
            return mass_ratio
        else:
            return self.calc_max_merger_mass_ratio()
        