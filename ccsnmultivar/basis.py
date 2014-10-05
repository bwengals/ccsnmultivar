# basis classes for ccsnmultivar
from sklearn.decomposition import KernelPCA
import csv as csv
import numpy as np



class PCA(object):
    """
    Performs SVD: Y = USV^\dagger

    PCA has 4 methods:
       - fit(waveforms)
       update class instance with pca fit, done via svd

       - fit_transform()
       do what fit() does, but additionally return the projection onto PC space
       i.e. returns A = U*S

       - inverse_transform(A)
       inverses the decomposition, returns waveforms for an input A, using Z
       returns Y = A*Z^\dagger

       - get_params()
       returns metadata used for fits.
            Number of components (pcs), name of fit ('PCA'), waveform catalog used (if provided)
    """

    def __init__(self,num_components=10,catalog_name='unknown'):
        self.decomposition  = 'PCA'
        self.num_components = num_components
        self.catalog_name   = catalog_name

    def fit(self,waveforms):
        # TODO make sure there are more columns than rows (transpose if not)
        self.waveforms = waveforms
        self.__do_pca()

    def fit_transform(self,waveforms):
        # TODO make sure there are more columns than rows (transpose if not)
        self.waveforms = waveforms
        self.__do_pca()
        return self.A

    def inverse_transform(self,A):
        # convert basis back to waveforms using fit
        Z = self.Z
        # for PCA, (linear) transform just matrix multiplication
        new_waveforms = np.dot(A,(self.Z).T)
        return new_waveforms

    def get_params(self):
        # return a dictionary with the pca instance metadata
        # TODO know what catalog was used! (include waveform metadata)
        params = {i: self.__dict__.get(i,None) for i in ('num_components')}
        return params

    def __do_pca(self):
        U,s,V = np.linalg.svd(self.waveforms,full_matrices=False)
        Z = V[0:self.num_components,:]
        Z = np.array(Z.T)
        A = np.array(np.matrix(U)*np.matrix(np.diag(s)))
        A = A[:,0:self.num_components]
        self.A = A
        self.Z = Z



class KPCA(object):
    """
    Wrapper for sklearn package.  Performs nonlinear kernel PCA

    PCA has 4 methods:
       - fit(waveforms)
       update class instance with pca fit

       - fit_transform()
       do what fit() does, but additionally return the projection onto PC space

       - inverse_transform(A)
       inverses the decomposition, returns waveforms for an input A, using Z

       - get_params()
       returns metadata used for fits.
            Number of components (pcs), name of fit ('PCA'), waveform catalog used (if provided)
    """
    def __init__(self,num_components=10,catalog_name='unknown',kernel='polynomial',
                 alpha = 0.1,
                 gamma = 1.0,
                 degree = 3):
        self.decomposition  = 'Kernel PCA'
        self.num_components = num_components
        self.catalog_name   = catalog_name
        self.gamma          = gamma
        self.kernel         = kernel
        self.alpha          = alpha
        self.KPCA = KernelPCA(kernel=self.kernel,
                         n_components=num_components,
                         fit_inverse_transform=True,
                         gamma = self.gamma)
        self.norm_constant = 1e20

    def fit(self,waveforms):
        # TODO make sure there are more columns than rows (transpose if not)
        # normalize waveforms
        self.waveforms = waveforms*self.norm_constant
        self.KPCA.fit(self.waveforms)

    def fit_transform(self,waveforms):
        # TODO make sure there are more columns than rows (transpose if not)
        # normalize waveforms
        self.waveforms = waveforms*self.norm_constant
        self.A = self.KPCA.fit_transform(self.waveforms)
        return self.A

    def inverse_transform(self,A):
        # convert basis back to waveforms using fit
        new_waveforms = self.KPCA.inverse_transform(A)*(1./self.norm_constant)
        return new_waveforms

    def get_params(self):
        # TODO know what catalog was used! (include waveform metadata)
        params = self.KPCA.get_params()
        params['num_components'] = params.pop('n_components')
        return params

def load_waveforms(path_to_waveforms):
    """
    Can load in either parameter file or waveforms

    Parameter files must be:
    - in .csv or .dat format
    - first row is the header
    - each column is a value of physical parameter

     Waveform matrices must be:
    - each row is a waveform, corresponding to the same row in parameter file
    - each waveform is preprocessed
        -- aligned to core bounce
        -- all have same sampling frequencies
        -- all have same number of time samples
    """
    cr = csv.reader(open(path_to_waveforms,"rb"))
    P = list(cr)
    Y = np.array(P).astype('float')
    return Y
