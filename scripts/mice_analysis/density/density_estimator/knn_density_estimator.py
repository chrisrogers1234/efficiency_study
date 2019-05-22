import math
import numpy as np
import scipy as sp
import scipy.spatial
import multiprocessing.dummy as mp
import ROOT

class kNNDensityEstimator(object):

    # Initialize the kNN density estimator
    #  -data		Array of data points
    #  -rotate		Evaluate the distances in the metric of the covariance matrix
    #  -nthreads	Number of threads available
    #  -norm		Normalisation factor (i.e. transmission, default to 1.)
    def __init__(self, data, rotate=True, nthreads=1, norm=1.):

    # Initialize the base members
    self.data = data
    self.n = len(data)
    self.dim = data.shape[1]
    self.norm = norm
    self.nthreads = nthreads

    # Set the optimal number of neighbors
    self.k = int(pow(self.n/self.norm, 4./(4.+self.dim)))

        # If requested, get the metric of the space and rotate all of the
    # points to fit a simple Euclidean L2 metric (used by the KD tree)
    self.rotate = rotate
    self.scale = 1.
    self.cov = np.cov(data.transpose())
    if rotate:
        # a) Get the scaling factor (generalised radius of the RMS ellipse)
        self.scale = math.sqrt(np.linalg.det(self.cov))

        # b) Find the eigenvalues and eigenvectors of the covariance matrix
        #    Sort them so that each eigenvec is associated with the right eigenval
        eigval, eigvec = np.linalg.eigh(self.cov)
        sort_perm = eigval.argsort()
        eigval.sort()
        eigvec = eigvec[:, sort_perm]

        # c) Find the metric to rotate the points into a Euclidean set
        invLsqrt = np.full([self.dim, self.dim], 0.)
        for i in range(self.dim):
            invLsqrt[i, i] = pow(eigval[i], -0.5)

        Ut = np.matrix(eigvec).transpose()
        self.metric = invLsqrt*Ut

        # d) Rotate all the points
        for i in range(self.n):
            x = self.metric*np.matrix(data[i]).transpose()
            self.data[i] = x.transpose()

    # Finally, initialize the underlying KD tree on the rotated data
    # The tree does not deep copy, do not touch data (!)
    self.kdtree = sp.spatial.cKDTree(self.data)

    # Initialize a container for the contour levels
    self.levels = []

    # Evaluate the kNN density in x
    def eval(self, x):

    # Rotate the vector in the same metric, if necessary
    if self.rotate:
        x = np.array(self.metric*np.matrix(x).transpose()).flatten()

        return self.norm*self.density(x)/self.scale

    # Returns the density profile as two arrays (values and stat errors)
    #  -npoints	Number of points contained in the graph (number of steps)
    #  -bsr	Use bootstrap resampling to evaluate uncertainties
    #  -bsr_n	Number of iterations of the bootstrap resampling
    def profile(self, npoints=1000, bsr=False, bsr_n=10):

        # Set the levels of the training set if they have not been set yet
        if not len(self.levels):
            self.set_levels()

        # Use the requested uncertainty calculation, produce profiles
        values, errors = [0. for i in range(npoints)], [0. for i in range(npoints)]
        if not bsr:
            # Take the quantiles as the estimate, the Gaussian error as the uncertainty
            for i in range(npoints):
                alpha = (float(i+1)/(npoints+1))
            level = 0.
            if alpha < self.norm:
                level = np.quantile(self.levels, 1.-alpha/self.norm)

                values[i] = level
                errors[i] = level*self.level_uncertainty(alpha)
        else:
            # Bootstrap the data, evaluate the profiles for each new sample
            baseline = self.data
            baseline_levels = self.levels
            levels = np.empty((npoints, bsr_n))
            for j in range(bsr_n):
                # Generate a new sample
                sample = np.empty((self.n, self.dim))
                for i in range(self.n):
                    k = np.random.randint(self.n)
                    sample[i] = baseline[k]

                # Reinitialize the KD tree on the new sample, reset the levels
                self.__init__(sample, self.rotate, self.nthreads, self.norm)
                self.set_levels()

                # Evaluate the array of levels of the new sample
                for i in range(npoints):
                    alpha = (float(i+1)/(npoints+1))
                    levels[i][j] = 0.
                    if alpha < self.norm:
                        levels[i][j] = np.quantile(self.levels, 1.-alpha/self.norm)

                # Take the mean as the estimate, the rms as the uncertainty
                for i in range(npoints):
                    values[i] = np.mean(levels[i])
                    errors[i] = np.std(levels[i])

                # Reset the data
                self.__init__(baseline, self.rotate, self.nthreads, self.norm)
                self.levels = baseline_levels

        # Return
        return values, errors

    # Returns the density profile as a graph
    #  -npoints	Number of points contained in the graph (number of steps)
    #  -bsr	Use bootstrap resampling to evaluate uncertainties
    #  -bsr_n	Number of iterations of the bootstrap resampling
    def profile_graph(self, npoints=1000, bsr=False, bsr_n=10):

        # Initialize the graph
        profile = ROOT.TGraphErrors(npoints)
        profile.SetTitle(";Fraction #alpha;#rho_{#alpha}")
        profile.GetYaxis().SetTitleOffset(1)

        # Set the points
        values, errors = self.profile(npoints, bsr, bsr_n)
        for i in range(npoints):
            alpha = (float(i+1)/(npoints+1))
            profile.SetPoint(i, alpha, values[i])
            profile.SetPointError(i, 0., errors[i])

        # Return
        return profile

    # Returns a lower dimensional section of the space at x
    #  -axes	List of axes ids to produce a Poincare section of
    #  -x	d-point that the Poincare section contains
    #  -mins	Minimums in each of the axes
    #  -maxs	Maximums in each of the axes
    def section(self, axes, x=[], mins=[], maxs=[]):

        # Check that the dimension of the section is supported and sensible
        secdim = len(axes)
        if not secdim:
            raise ValueError('No axis(es) specified')
        if secdim > 2:
            raise ValueError('Sections of dimension larger than 2 not supported')
        for axis in axes:
            if axis >= self.dim:
                raise ValueError('Axis id %d is not compatible with space dimension %d'\
                                     %(axis, self.dim))

        # If an intersection point is provided, check that it is of the right dimension
        # If not provided, take the point at the origin of the space
        if len(x) and len(x) != self.dim:
            raise ValueError('The intersection point must be of dimension %d'%self.dim)
        if not len(x):
            x = np.full(self.dim, 0.)

        # Checks if the boundaries are of the right dimension
        # If not provided, take the 4 sigma interval around the mean
        if len(mins) and len(mins) != secdim:
            raise ValueError('The minimum point must be of dimension %d'%self.dim)
        if len(maxs) and len(maxs) != secdim:
            raise ValueError('The maximum point must be of dimension %d'%self.dim)
        if not len(mins) or not len(maxs):
            mins = np.empty(secdim)
            maxs = np.empty(secdim)
            for i in range(secdim):
                sample = self.data.transpose()[axes[i]]
                mean = np.mean(sample)
                std = math.sqrt(self.cov[axes[i]][axes[i]])
                mins[i] = mean-4*std
                maxs[i] = mean+4*std

        # Create a histogram
        section = None
        if secdim == 1:
            nbins = 100
            section = ROOT.TH1F("section_%d"%(axes[0]),\
                ";x_{%d};#rho(#bf{x})"%(axes[0]),\
                nbins, mins[0], maxs[0])
            for i in range(nbins):
                pos = x
                pos[axes[0]] = section.GetXaxis().GetBinCenter(i+1)
                val = self.eval(pos)
                section.SetBinContent(i+1, val)

        if secdim == 2:
            nbins = 50
            section = ROOT.TH2F("section_%d%d"%(axes[0], axes[1]),\
                ";x_{%d};x_{%d};#rho(#bf{x})"%(axes[0], axes[1]),\
                nbins, mins[0], maxs[0], nbins, mins[1], maxs[1])
            for i in range(nbins):
                for j in range(nbins):
                    pos = x
                    pos[axes[0]] = section.GetXaxis().GetBinCenter(i+1)
                    pos[axes[1]] = section.GetYaxis().GetBinCenter(j+1)
                    val = self.eval(pos)
                    section.SetBinContent(i+1, j+1, val)

        return section

    # Returns the density in the Euclidean metric in x
    def density(self, x):
        vol = self.volume(x)
        return float(self.k)/self.n/vol

    # Returns the volume of the d-ball centred in x of radius R_k
    def volume(self, x):
        dist = self.distance(x)
        return pow(math.sqrt(math.pi)*dist, self.dim)/sp.special.gamma(self.dim/2+1)

    # Returns the distance to the k^th closest point, R_k
    def distance(self, x):
        dists, ids = self.kdtree.query(x, self.k)
        return dists[-1]

    # Sets the contour levels (evaluate the density in all of the points)
    def set_levels(self):
        self.levels = np.empty(self.n)
        # Single-threaded approach
        if self.nthreads == 1:
            for i in range(self.n):
                set_level(i)

        # Multi-threaded approach
        else:
            pool = mp.Pool(processes=self.nthreads)
            pool.map(self.set_level, range(self.n))
            pool.close()
            pool.join()

    # Set a specific contour level (in one of the points)
    def set_level(self, i):
        self.levels[i] = self.norm*self.density(self.data[i])/self.scale

    # Returns the Gaussian uncertainty associated with a contour level
    def level_uncertainty(self, frac):
        return pow(4./(frac*(1-frac)), 1./4)/math.sqrt(self.n)
