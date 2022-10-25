from halotools.empirical_models import PrebuiltHodModelFactory
from scipy.integrate import simps
import numpy as np

from projector import radec_to_xy

### Define Cosmology
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
H0 = 70
cosmo = FlatLambdaCDM(H0=H0, Om0=0.283)

rad2deg = 180.0/np.pi

class corcovadoShape:
    """This class estimates shape of a given stellar mass cluster distribution
    """
    def __init__(self, coords=[0,0], z=0, sigma_bg=0) -> None:
        """__init__ Corcovado Shape Measurements Tools

        Estimates the shape of a given stellar mass cluster distribution.

        Parameters
        ----------
        coords : list, optional
            ra, dec of the center of the distribution, by default [0,0]
        z : int, optional
            redshift of the cluster, by default 0
        sigma_bg : int, optional
            background density, by default 0
        """
        self.coords = coords
        self.z = z
        self.sigma_bg = sigma_bg
        self.angularDiameter= cosmo.luminosity_distance(z).value/(1+z**2) ## angular diamteter at redshift z
        self.Mpc2theta = self.angularDiameter/rad2deg # conversion factor to degrees
        
        
    def display(self, rmax=3, kind='points', save=False, ax=None):
        """display point density

        Display cutout of the point distribution

        Parameters
        ----------
        rmax : int, optional
            max size of the box, by default 3
        kind : str, optional
            display can be only points or a colored kernel distribtuion, by default 'points'
        """
        if ax is None: ax = plt.gca()
        if kind=='points':
            plot_points(self.dx_weights, self.dy_weights, rmax, save=save, ax=ax)
        if kind=='kde':
            plot_kde(self.dx_weights, self.dy_weights, rmax, save=save, ax=ax)
            
        pass
    
    def load_sky_coord(self, ra, dec, weights=None):
        """load_sky_coord 

        Converts sky coordinates to physical coords system in Mpc

        Parameters
        ----------
        ra : array
            right ascension
        dec : array
            declination
        weights : array, optional
            weighted value for the point distribution, by default None
            
        Returns
        ----------
        dx : array
            x offset from the center
        dy : array
            y offset from the center
        """
        ra_center, dec_center = self.coords
        
        dx,dy,albers = radec_to_xy(ra, dec, ra_center, dec_center,
                                   self.Mpc2theta) ## Mpc
        # dx and dy in Mpc not Mpc/h
        self.dx = dx
        self.dy = dy
        
        if weights is not None:
            self.weight_xy(weights)
        pass
    
    def weight_xy(self, weights):
        """add_weights

        Weights xy quantities.
        Repeats xy values n times the weights.
        Solution proposed at stackOverflow: 
        https://stackoverflow.com/questions/11442962/python-package-that-supports-weighted-covariance-computation
        
        Parameters
        ----------
        weights : array
            weights values to repeat the xy quantities.
        """
        if not isinstance(np.sum(weights), int):
            weights = (weights/np.min(weights)).astype(int)
            
        self.dx_w = np.repeat(self.dx, weights, axis=0)
        self.dy_w = np.repeat(self.dy, weights, axis=0)
    
    
    def load(self, x, y, weights=None):
        """load dataset

        load point distribution

        Parameters
        ----------
        x : array
            x-coordinates
        y : array
            y-coordiantes 
        weights : array, optional
            weighted value for the point distribution, by default None
        """
        pass
    
    def fit(self):
        """fit shape distribution

        Uses the PCA to find the major and minor axis. 
        Converts the results to polar coordinates, eccentricity and position angle.
        """
        pass
    
    def fit_outliers(self, **kwargs):
        """mask_outliers detection algorithm

        Finds the outliers galaxies in the field. 
        Uses the LocalOutlierFactor function from sklearn.
        For the key word arguments (kwargs) take a look at the documentation 
        https://scikit-learn.org/stable/modules/generated/sklearn.neighbors.LocalOutlierFactor.html#sklearn.neighbors.LocalOutlierFactor
        """
    pass