import numpy as np
from .galaxy import mean_stellar_formation_time
from simulation_data import get

class GalaxyPopulation():
    
    def __init__(self):
        self.ids = []
    
    #load data directly
    get_mean_stellar_age = np.loadtxt('z=2_Mean_SFT', dtype=float)
    get_timeavg_stellar_formation_rate_1gyr = np.loadtxt('z=2_TimeAvg_SFR_1Gyr', dtype=float)
    get_current_stellar_formation_rate = np.loadtxt('z=2_Current_SFR', dtype=float)
    get_median_stellar_age = np.loadtxt('z=2_Median_SFT', dtype=float)
    
    
    #select ids
    def select_galaxies(self, redshift, mass_min, mass_max=12):
        h = 0.6774
        mass_minimum = 10**mass_min / 1e10 * h
        mass_maximum = 10**mass_max / 1e10 * h
        # form the search_query string by hand for once
        search_query = "?mass_stars__gt=" + str(mass_minimum) + "&mass_stars__lt=" + str(mass_maximum)
        url = "http://www.tng-project.org/api/TNG100-1/snapshots/z=" + str(redshift) + "/subhalos/" + search_query
        subhalos = get(url, {'limit':5000})
        ids = [ subhalos['results'][i]['id'] for i in range(subhalos['count'])]
        return ids
    
    
    #mean SFT
    def calc_mean_stellar_formation_time(self, calc_redshift, calc_mass_min):
        #create and populate array for mean SFT
        ids = self.select_galaxies(redshift = calc_redshift, mass_min = calc_mass_min, mass_max=12)
        means = np.zeros(len(ids))
        for i, id in enumerate(ids): 
            means[i] = mean_stellar_formation_time(z = calc_redshift, subhalo_id = id)
        #save file
        np.savetxt('z='+str(calc_redshift)+'_Mean_SFT', means)
        mean_SFT = np.loadtxt('z='+str(calc_redshift)+'_Mean_SFT', dtype=float)
        return mean_SFT
    
    
    def get_mean_stellar_formation_time(self, z, mass_min): #can parameterize for max mass if needed 
        import pathlib
        file = pathlib.Path('z='+str(z)+'_Mean_SFT')
        if file.exists ():
            MeanSFT = np.loadtxt('z='+str(z)+'_Mean_SFT', dtype=float) #rename pre-existing files before parameterizing further
            return MeanSFT
        else:
            self.calc_mean_stellar_formation_time(calc_redshift = z, calc_mass_min = mass_min)
           
    
    #time avg SFR
    def calc_timeaverage_stellar_formation_rate(self, calc_redshift, calc_timescale, calc_mass_min):
        #create and populate array for mean SFT
        ids = self.select_galaxies(redshift = calc_redshift, mass_min = calc_mass_min, mass_max=12)
        time_averages = np.zeros(len(ids))
        for i, id in enumerate(ids): 
            time_averages[i] = timeaverage_stellar_formation_rate(z = calc_redshift, subhalo_id = id, timescale = calc_timescale)
        #save file
        np.savetxt( 'z='+str(calc_redshift)+ '_TimeAvg_SFR_'+ str(calc_timescale) +'Gyr', time_averages)
        time_avg_SFT = np.loadtxt('z='+ str(calc_redshift) +'_TimeAvg_SFR_'+ str(calc_timescale)+ 'Gyr', dtype=float)
        return time_avg_SFT
    
        
    def get_timeaverage_stellar_formation_rate(self, z, timescale, mass_min): #can parameterize for max mass if needed 
        import pathlib
        file = pathlib.Path('z='+str(z)+'_TimeAvg_SFR_'+str(timescale)+'Gyr')
        if file.exists ():
            time_avg_SFT = np.loadtxt('z='+str(z)+ '_TimeAvg_SFR_'+ str(timescale) +'Gyr', dtype=float) #rename pre-existing files before parameterizing further
            return time_avg_SFT
        else:
            self.calc_timeaverage_stellar_formation_rate(calc_redshift=z, calc_timescale=timescale, calc_mass_min=mass_min)
    
    
    #current SFR
    def calc_current_stellar_formation_rate(self, calc_redshift, calc_mass_min):
        #create and populate array for mean SFT
        ids = self.select_galaxies(redshift = calc_redshift, mass_min = calc_mass_min, mass_max=12)
        current_SFRs = np.zeros(len(ids))
        for i, id in enumerate(ids): 
            current_SFRs[i] = timeaverage_stellar_formation_rate(z = calc_redshift, subhalo_id = id, timescale = 0)
        #save file
        np.savetxt( 'z='+ str(calc_redshift) +'_Current_SFR', current_SFRs)
        current_SFR = np.loadtxt('z='+ str(calc_redshift) +'_Current_SFR', dtype=float)
        return current_SFR
    
        
    def get_current_stellar_formation_rate(self, z, mass_min): #can parameterize for max mass if needed 
        import pathlib
        file = pathlib.Path('z='+ str(z) +'_Current_SFR')
        if file.exists ():
            current_SFR = np.loadtxt('z='+ str(z) +'_Current_SFR', dtype=float) #rename pre-existing files before parameterizing further
            return current_SFR
        else:
            self.calc_current_stellar_formation_rate(calc_redshift=z, calc_mass_min = mass_min)    
    
    
    #median SFT
    def calc_median_stellar_formation_time(self, calc_redshift, calc_mass_min):
        #create and populate array for mean SFT
        ids = self.select_galaxies(redshift = calc_redshift, mass_min = calc_mass_min, mass_max=12)
        MedianSFT = np.zeros(len(ids))
        for i, id in enumerate(ids):
            MedianSFT[i] = median_stellar_formation_time(z = calc_redshift, subhalo_id = id)
        #save file
        np.savetxt('z='+ str(calc_redshift) +'_Mean_SFT', MedianSFT)
        median_SFT = np.loadtxt('z='+ str(calc_redshift) +'_Mean_SFT', dtype=float)
        return median_SFT
    
    
    def get_median_stellar_formation_time(self, z, mass_min):
        import pathlib
        file = pathlib.Path('z='+ str(z) +'_Mean_SFT')
        if file.exists ():
            median_SFT = np.loadtxt('z='+ str(z) +'_Mean_SFT', dtype=float) #rename pre-existing files before parameterizing further
            return median_SFT
        else:
            self.calc_median_stellar_formation_time(calc_redshift = z, calc_mass_min = mass_min)  
    