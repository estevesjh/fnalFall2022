import numpy as np


### Define Cosmology
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u

cosmo = FlatLambdaCDM(H0=70, Om0=0.283)


class r200SHMR:
    """This class estimates r200c of a cluster
    """
    def __init__(self, rbins, z=0, sigma_bg=0, fmasked=0) -> None:
        """__init__ start the r200c estimator 

        Define the radial binning scheme, the redshift and the background density

        Parameters
        ----------
        rbins : array
            radial bins used on the computation
        z : int, optional
            redshift of the cluster, by default 0
        sigma_bg : int, optional
            bagrkound density, no contamination by defalt
        fmasked : int, optional
            masked area fraction, by default 0
        """
        print(4*'-----')
        print('R200c: Stellar to Halo Mass Estimator')
        self.rbins = rbins
        self.rmed = 0.5*(rbins[1:]+rbins[:-1])
        
        self.z = z
        self.sigma_bg = sigma_bg
        self.area = np.pi*self._rmed**2
        self.volume = 4*np.pi*self._rmed**3/3
        pass
    
    def fit(self, mstar, radii):
        """fits the radii with 200 the critical density of the universe

        Parameters
        ----------
        mstar : array
            galaxies stellar masses
        radii : array
            galaxies cluster centric distance
        """
        self.compute_stellar_mass_density(mstar, radii)
        self.compute_halo_mass()
        self.compute_critical_stellar_mass()
        pass
    
    def compute_halo_mass(self):
        """compute_halo_mass based on the stellar halo mass relation
        
        We assume the stellar halo mass parameters from bla et al.
        """
        
    def compute_critical_stellar_mass(self):
        """compute critical stellar mass of the universe

        $M_{\star,c}$ is based on the critical density. 
        For a given volume, we can associate a critical halo mass.
        For this critical halo mass we assume a stellar to halo mass relation and compute the critical stellar mass of the universe.
        """
        self._rhoc = (cosmo.critical_density(self.z).to(u.Msun/u.Mpc**3)).value # Msun/Mpc^3
        #self.lo
    
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
        smass_surface_density_total = self.compute_density(radii, mstar)
        
        # compute the cumulative stellar mass inside R subtracted by the backround
        self.smass_cluster = (smass_surface_density_total-self.sigma_bg)*self.area
        
        # compute the halo mass associated with the stellar inside R
        # assume the stellar to halo mass relation
        self.compute_halo_mass()
        pass
        
    def compute_density(self, x, weights):
        """compute radial density

        Parameters
        ----------
        x : array
            x-variable to compute the profile
        weights : array
            cumulative sum of the weights
        """
        return np.histogram(x, weights=weights, rbins=self.rbins)/self.area
        
        

### Compute R200c based on the HOD model
def computeR200(gals, cat, nbkg, rmax=3, defaultMass=1e14,testPz=False,compute=True):
    ## estimate R200
    ncls = len(cat)
    r200m = []

    for idx in range(ncls):
        cls_id, z_cls = cat['CID'][idx], cat['redshift'][idx]
        magLim_i = cat['magLim'][idx,1]

        gal = gals[(gals['CID']==cls_id)&(gals['mag'][:,2]<=magLim_i)]

        if compute:
            r200i = calcR200(gal['R'],np.ones_like(gal['PDFz']),cls_id,z_cls,nbkg[idx],rmax=rmax,testPz=testPz)
        else:
            r200i = 0.1

        r200i = checkR200(r200i,z_cls,M200=defaultMass)
        print('r200:',r200i)
        r200m.append(r200i)
    
    return np.array(r200m)

def calcR200(radii,pz,cls_id,z_cls,nbkg,rmax=3,testPz=False):
    ngals_cls, rbin = doRadialBin(radii,pz,width=0.1,testPz=testPz)

    ngals = ngals_cls-nbkg*np.pi*rbin**2
    ####params=[11.6,12.45,1.0,12.25,-0.69]#parameters for mass conversion - see table 4 in Tinker paper
    params = [11.59,12.94,1.01,12.48,-0.69]#parameters for mass conversion - see table 4 in Tinker paper
    mass = hod_mass_z(ngals,z_cls,params) #calculate mass given ngals (see above functions)

    volume = (4./3)*np.pi*rbin**3
    mass_density=(mass)/volume

    mass_density=np.where(mass_density<0.1,1e11,mass_density)
    pc=200*np.ones_like(radii)
    rho_crit = rhoc(0)
    critdense1 = crit_density(rho_crit,z_cls,0.23,0.77)
    critdense = critdense1*np.ones_like(rbin)

    X=2000 #desired excess over critical density, ex. if X=200, calculates R/M200
    dX=10  #acceptance window around X
    ratio=mass_density/critdense

    f=interp1d(rbin,ratio,fill_value='extrapolate')
    radii_new=np.linspace(0.1,rmax,10000)
    ratio_new=f(radii_new)
    r200m=radii_new[np.where( (ratio_new>=X-dX)&(ratio_new<=X+dX) )] #find possible r200s within acceptance range
    
    if r200m.size > 0:
        r200m=np.median(r200m) #mean of all possible r200s is measured r200

    else:
        ## Try r500
        X=500
        r200m=radii_new[np.where( (ratio_new>=X-dX)&(ratio_new<=X+dX) )] #find possible r200s within acceptance range
        
        if r200m.size > 0:
            r200m=np.median(r200m)/0.65 #mean of all possible r200s is measured r200
        else:
            r200m = 0
            print('bad cluster:',cls_id,'ratio min/max:',min(ratio_new),max(ratio_new))

    return r200m