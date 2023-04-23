import os
import numpy as np
import kmeans_radec
from kmeans_radec import KMeans, kmeans_sample

class JackKniferKmeans(object):
    """
    Define the Jackknife labels using Kmeans
    
    Args:
        Npatches (int): number of patches
        ra (array): right ascension
        dec (array): declination
    """
    def __init__(self, ra, dec, npatches, fname='center.txt'):
        self.pos = np.vstack([ra,dec]).T
        self.npatches = npatches
        
        if os.path.isfile(fname):
            self.load_center(fname)
        else:
            self.fit_centers()
        
        # assign label and cluster size variables
        self.labels = self.km.labels
        self.cluster_sizes = np.bincount(self.km.labels)
    
    def load_center(self, fname):
        # load a saved center file
        centers = np.load(fname)
        # initialize the class
        self.km = KMeans(centers)
        # run the algorithm
        self.km.run(self.pos, maxiter=20)
            
    def fit_centers(self, verbose=0):
        km = kmeans_sample(self.pos, ncen=self.npatches, verbose=verbose, maxiter=100)

        # check convergence
        if not km.converged:
            print('Kmeans did not converge with 100 iterations')
            print('Trying with an undefined number of iterations')
            km.run(pos)
        
        self.km = km
        pass
    
    def show_stats(self):
        print('Number of patches: ', self.npatches)
        print('Cluster Sizes: ', self.cluster_sizes)
        print('Labels: ', self.labels)
        print('')

    def get_mask(self, k_id, labels=None):
        if labels is None: labels = self.labels
        mask = labels != k_id
        return mask
    
    def add_randoms(self, ra, dec):
        pos_ran = np.vstack([np.array(ra),np.array(dec)]).T
        return self.km.find_nearest(pos_ran)
    
    def write(self, fname):
        """Save Kmeans fitted centers
        
        note: used later to predict the labels
        """
        np.save(fname, self.km.centers)
        
    def refit_centers(self, centers=None):
        if centers is None: centers = self.km.centers
        self.km.set_centers()
        self.km.run(self.pos)
        self.update()
        
    def update(self):
        self.labels = self.km.labels
        self.cluster_sizes = np.bincount(self.km.labels)
        
    def plot_groups(self, Npoints=10000):
        plot_proj(self.pos[:,0], self.pos[:,1], self.labels, Npoints=Npoints)
    
    def plot_groups_masked(self, k_id, Npoints=10000):
        plt.clf()
        mask = self.get_mask(k_id)
        plot_proj(self.pos[mask,0], self.pos[mask,1], self.labels[mask], Npoints=Npoints)

#     ############# TBD #############
#     def drop_bad_groups(self, frac_th=0.5):
#         """
#         Groups with less than a threshold fraction of the median number of a values 
#         will be added to next group
        
#         Args:
#             frac_thr (float): minimum fraction of a given jackknife region
#             that must be unmasked for that region to be included in
#             the set of regions.
#         """
#         self.bad_groups = 
#         self.labels = new_labels

###### Plot The Results #######
try:
    import skymapper as skm
except:
    print('Please install skymapper to plot the results')

def plot_proj(ra, dec, labels, Npoints=10000):
    idx = np.random.randint(len(ra), size=Npoints)
    # define the best Albers projection for the footprint
    # minimizing the variation in distortion
    crit = skm.stdDistortion
    proj = skm.WagnerIV.optimize(ra[idx], dec[idx], crit=crit)

    # construct map: will hold figure and projection
    # the outline of the sphere can be styled with kwargs for matplotlib Polygon
    map = skm.Map(proj)

    # add graticules, separated by 15 deg
    # the lines can be styled with kwargs for matplotlib Line2D
    # additional arguments for formatting the graticule labels
    sep=15
    map.grid(sep=sep,)

    # add scatter plot
    map.scatter(ra[idx], dec[idx], c=labels[idx], cmap='tab20', s=25, edgecolor='w', facecolor='None')


    # # focus on relevant region
    map.focus(ra[idx], dec[idx], pad=0.1)

    map.title('DES K-Means')

    # hide x-ticks
    ax = map.ax
    labels = [item.get_text() for item in ax.get_xticklabels()]

    empty_string_labels = ['']*len(labels)
    ax.set_xticklabels(empty_string_labels)
