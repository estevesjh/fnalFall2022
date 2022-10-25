from halotools.empirical_models import PrebuiltHodModelFactory
from scipy.integrate import simps
from scipy.interpolate import interp1d
import numpy as np

### Define Cosmology
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
H0 = 70
cosmo = FlatLambdaCDM(H0=H0, Om0=0.283)

class r200SHMR:
    """This class estimates r200c of a cluster
    """
    def __init__(self, rbins, model_name='leauthaud11', z=0, sigma_bg=0) -> None:
        """__init__ start the r200c estimator 

        Define the radial binning scheme, the redshift and the background density

        Parameters
        ----------
        rbins : array
            radial bins used on the computation
        model_name: str, optional
            stellar-to-halo mass relation model name following prebuiltHaloFactory
        z : int, optional
            redshift of the cluster, by default 0
        sigma_bg : int, optional
            bagrkound density, no contamination by defalt
        fmasked : int, optional
            masked area fraction, by default 0
        """
        #print(4*'-----')
        #print('R200c: Stellar to Halo Mass Estimator')
        
        # define radial bins
        self.rbins = rbins
        self.rmed = 0.5*(rbins[1:]+rbins[:-1])
        self.area = np.pi*self.rmed**2
        
        # define other variables
        self.sigma_bg = 10**sigma_bg
        self.z = z
        
        # define stellar-to-halo mass relation model
        self.model_name = model_name
        
        # compute the universe critical density at redshift z
        self._init_critical_density()
        pass
    
    def fit(self, mstar, radii, bias=0.0, fit_der=False):
        """fits the radii with 200 the critical density of the universe

        Parameters
        ----------
        mstar : array
            galaxies stellar masses
        radii : array
            galaxies cluster centric distance
        """        
        # compute cumulative cluster stellar mass, note the background is subtracted.
        self.compute_stellar_mass_density(mstar, radii)
        
        # stellar mass thresholds
        log_ms_low, log_ms_hig = np.nanpercentile(mstar, [0, 100])

        # predict halo mass as a function of radii using the stellar to halo mass relation
        self.compute_halo_mass(self.model_name, log_ms_low, log_ms_hig)
        
        # fit r200c based on the critical density of the unvierse
        self.compute_r200c(delta=200, bias=bias)
        
        if fit_der:
            # use inital guess to recompute the stellar density
            self.compute_r200c_derivative(mstar, radii)
        pass
        
    def compute_r200c_derivative(self, mstar, radii, nobj=3):
        rad_cut = radii[radii <= 1.75*self.r200c]
        nbins = int(rad_cut.size/nobj)+1
        self.rbins = np.percentile(rad_cut, np.linspace(0,100,nbins))
        self.rmed_cut = 0.5*(self.rbins[1:]+self.rbins[:-1])
        self.area = np.pi*self.rmed_cut**2
        self.compute_stellar_mass_density(mstar, radii)
        self.shmr_halo_mass_updated = np.log10(np.interp(10**self.smass_cluster, self.shmr_total_smass, self.bin_halo_mass))
        self.fit_halo_mass_poly_derivative(p=3)
        # output self.r200c_der and self.m200c_der
        pass
    
    def compute_r200c(self, bias=-0.32, delta=200, th=0.015, window=10):
        """compute_r200c
        
        function under construction
        
        # for instance, r200 is the distance associated with the total halo mass
        take the derivative point.

        Parameters
        ----------
        delta : int, optional
            delta times the critical density of the universe, eg. 200 rho_c, by default 200
        """
        self.shmr_halo_mass = self.shmr_halo_mass-bias
        #self.fit_halo_mass_poly_derivative(bias)
        self.fit_critical_density()
        #self.r200c = convertM200toR200(10**self.haloMax, self._rhoc, delta=delta)/(H0/100.)
        pass
        
    def compute_stellar_mass_density(self, mstar, radii):
        """compute stellar mass density a given aperture

        Parameters
        ----------
        mstar : array
            galaxy stellar mass
        radii : array
            cluster centric distance
        """
        # compute the total surface stellar mass density
        smass_density_total = self.compute_density(radii, 10**mstar)
        smass_cum_total = np.cumsum(smass_density_total,axis=0)
        
        # compute the cumulative stellar mass inside R subtracted by the backround
        self.smass_cluster = np.log10(smass_cum_total - self.sigma_bg*self.area)
        self.smass_density_cluster = (smass_density_total - self.sigma_bg)#*self.area/self.volume
    
    def compute_halo_mass(self, model_name, log_ms_low=10., log_ms_hig=12.5, nbins=50):
        """compute_halo_mass based on the stellar-to-halo mass relation (SHMR)
        
        We assume the SHMR parameters from Sunecsh et al. 2022 fited on the COSMOS2020 dataset.
        Following HOD parametrization from Leauthaud et al. 2011 with no redshift evolution.
        """
        # define binning scheme
        self.bin_halo_mass = np.logspace(11.5, 17.0, nbins+25)
        self.bin_log_stellar_mass = np.linspace(log_ms_low, log_ms_hig, nbins)

        # compute cluster stellar mass
        self.shmr_total_cluster_stellar_mass(model_name)
        
        # predict halo mass
        self.shmr_halo_mass = np.log10(np.interp(10**self.smass_cluster, self.shmr_total_smass, self.bin_halo_mass))
    
    def shmr_total_cluster_stellar_mass(self, model_name):
        # compute Ntot = Nsat+Ncen as function of stellar and halo mass
        self.shmr_cen_sat_stellar_mass(model_name)
        
        # compute cluster stellar mass as a function of halo mass
        self.shmr_cumulative_stellar_mass()
        pass
    
    def shmr_cen_sat_stellar_mass(self, model_name):
        mean_ncen_smass = np.ones((self.bin_log_stellar_mass.size, self.bin_halo_mass.size))
        mean_nsat_smass = np.ones((self.bin_log_stellar_mass.size, self.bin_halo_mass.size))
        for i,logMt in enumerate(self.bin_log_stellar_mass):
            model = PrebuiltHodModelFactory(model_name, threshold = logMt,
                                            redshift=self.z, haloprop_key='halo_mvir')
            mean_ncen_smass[i] = model.mean_occupation_centrals(prim_haloprop = self.bin_halo_mass)
            mean_nsat_smass[i] = model.mean_occupation_satellites(prim_haloprop = self.bin_halo_mass)
            
        ## predicted number of galaxies (haloMass, stellarMass)
        self.shmr_ntot = mean_ncen_smass+mean_nsat_smass
        pass
    
    
    def shmr_cumulative_stellar_mass(self):
        # self.shmr_ntot number of central+number of satelites as a function of stellar mass and halo mass
        
        # compute total stellar mass as a function of halo mass
        total_smass = simps(self.shmr_ntot, x=10**self.bin_log_stellar_mass, axis=0) 
        total_smass-= (self.shmr_ntot[-1]*10**self.bin_log_stellar_mass[-1] -self.shmr_ntot[0]*10**self.bin_log_stellar_mass[0])
        
        # total cluster stellar mass as a function of halo mass 
        self.shmr_total_smass = total_smass
        pass
    
    def fit_critical_density(self, rmax=8, delta=200, dX=20):
        self.shmr_mass_density = (10**self.shmr_halo_mass)/self.area
        critdense = self._rhoc*self.rmed

        ratio = self.shmr_mass_density/critdense

        f = interp1d(self.rmed, ratio, fill_value='extrapolate')
        r_new = np.linspace(0.1,rmax,10000)
        
        deltas = np.linspace(delta-dX, delta+dX, 10)
        r200m = np.median(interp1d(np.log10(f(r_new)),r_new)(np.log10(deltas)))
                
        self.m200c = np.interp(r200m, self.rmed, self.shmr_halo_mass)
        self.r200c = convertM200toR200(10**self.m200c, self._rhoc, delta=delta) # in Mpc

    def fit_halo_mass_poly_derivative(self, p=3, bias=0.):
        sm = smoothP(self.rmed_cut, self.shmr_halo_mass_updated, p, deriv=0)
        d1 = smoothP(self.rmed_cut, self.shmr_halo_mass_updated, p, deriv=1)
        d2 = smoothP(self.rmed_cut, self.shmr_halo_mass_updated, p, deriv=2)

        if np.min(d1)>0:
            m200c = interp1d(d1, sm)(np.min(d1))
        elif d2[0]>0:
            m200c = interp1d(d2, sm, fill_value='extrapolate')(0) 
        else:
            m200c = interp1d(d1, sm, fill_value='extrapolate')(0) 

        # r200c = interp1d(sm, self.rmed_cut, fill_value='extrapolate')(m200c) 
        r200c = convertM200toR200(10**m200c, self._rhoc, delta=200) # in Mpc
        self.r200c_der = r200c
        self.m200c_der = m200c
                
    def _init_critical_density(self):
        """compute critical stellar mass of the universe

        $M_{\star,c}$ is based on the critical density. 
        For a given volume, we can associate a critical halo mass.
        For this critical halo mass we assume a stellar to halo mass relation and compute the critical stellar mass of the universe.
        """
        self._rhoc = (cosmo.critical_density(self.z).to(u.Msun/u.Mpc**3)).value # Msun/Mpc^3
            
    def compute_density(self, x, weights):
        """compute radial density

        Parameters
        ----------
        x : array
            x-variable to compute the profile
        weights : array
            cumulative sum of the weights
        """
        return np.histogram(x, weights=weights, bins=self.rbins)[0]

    def compute_density_kde(self, x, weights,bw=0.2):
        """compute radial density

        Parameters
        ----------
        x : array
            x-variable to compute the profile
        weights : array
            cumulative sum of the weights
        """
        p, N = compute_kde(x,weights)
        
        return np.cumsum(p(self.rmed)*N*bw)
    
def smoothP(x,y, window=3, deriv=0):
    idx = np.isfinite(x) & np.isfinite(y)
    coefs = np.polyfit(x[idx], y[idx], deg=window)
    xnew = np.arange(0., 4., 100)
    if deriv==0:
        return np.poly1d(coefs)(x)
    else:
        return np.poly1d(coefs).deriv(deriv)(x)

def convertM200toR200(M200,rho,delta=200):
    ## M200 in solar masses
    ## R200 in Mpc
    R200 = ( M200/(delta*4*np.pi*rho/3) )**(1/3.)
    return R200

#from scipy.optimize import curve_fit
#def a function
# def f(x, a, b, c):
#     return a*x**3 + b*x**2 + c*x

# #do the fitting
# popt, pcov = curve_fit(f, xs, ys)
    
def smoothP(x,y, window=3, deriv=0):
    idx = np.isfinite(x) & np.isfinite(y)
    coefs = np.polyfit(x[idx], y[idx], deg=window)
    xnew = np.arange(0., 4., 100)
    if deriv==0:
        return np.poly1d(coefs)(x)
    else:
        return np.poly1d(coefs).deriv(deriv)(x)

def convertM200toR200(M200,rho,delta=200):
    ## M200 in solar masses
    ## R200 in Mpc
    R200 = ( M200/(delta*4*np.pi*rho/3) )**(1/3.)
    return R200

#from scipy.optimize import curve_fit
#def a function
# def f(x, a, b, c):
#     return a*x**3 + b*x**2 + c*x

# #do the fitting
# popt, pcov = curve_fit(f, xs, ys)