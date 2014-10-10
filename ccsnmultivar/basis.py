# basis classes for ccsnmultivar
from sklearn.decomposition import KernelPCA, FastICA
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

    def __init__(self, num_components=10,catalog_name='unknown'):
        self._decomposition  = 'PCA'
        self._num_components = num_components
        self._catalog_name   = catalog_name

    def fit(self, waveforms):
        # TODO make sure there are more columns than rows (transpose if not)
        self._waveforms = waveforms
        self._do_pca()

    def fit_transform(self, waveforms):
        # TODO make sure there are more columns than rows (transpose if not)
        self._waveforms = waveforms
        self._do_pca()
        return self._A

    def inverse_transform(self,A):
        # convert basis back to waveforms using fit
        Z = self._Z
        # for PCA, (linear) transform just matrix multiplication
        new_waveforms = np.dot(A,(self._Z).T)
        return new_waveforms

    def get_params(self):
        # return a dictionary with the pca instance metadata
        params = {}
        params['Decomposition'] = self._decomposition
        params['num_components'] = self._num_components
        return params

    def _do_pca(self):
        U,s,V = np.linalg.svd(self._waveforms,full_matrices=False)
        Z = V[0:self._num_components,:]
        Z = np.array(Z.T)
        A = np.array(np.matrix(U)*np.matrix(np.diag(s)))
        A = A[:,0:self._num_components]
        self._A = A
        self._Z = Z



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
    def __init__(self, num_components=10,catalog_name='unknown',kernel='polynomial',
                 alpha = 0.1,
                 gamma = 10.0,
                 degree = 2):

        self._decomposition  = 'Kernel PCA'
        self._num_components = num_components
        self._catalog_name   = catalog_name
        self._gamma          = gamma
        self._kernel         = kernel
        self._alpha          = alpha
        self._KPCA = KernelPCA(kernel=self._kernel,
                         n_components=num_components,
                         fit_inverse_transform=True,
                         gamma = self._gamma)

    def fit(self,waveforms):
        # TODO make sure there are more columns than rows (transpose if not)
        # normalize waveforms
        self._waveforms = waveforms
        self._KPCA.fit(self._waveforms)

    def fit_transform(self,waveforms):
        # TODO make sure there are more columns than rows (transpose if not)
        # normalize waveforms
        self._waveforms = waveforms
        self._A = self._KPCA.fit_transform(self._waveforms)
        return self._A

    def inverse_transform(self,A):
        # convert basis back to waveforms using fit
        new_waveforms = self._KPCA.inverse_transform(A)
        return new_waveforms

    def get_params(self):
        # TODO know what catalog was used! (include waveform metadata)
        params = self._KPCA.get_params()
        params['num_components'] = params.pop('n_components')
        params['Decompositon'] = self._decomposition
        return params


class ICA(object):
    """
    Wrapper for sklearn package.  Performs fast ICA (Independent Component Analysis)

    ICA has 4 methods:
       - fit(waveforms)
       update class instance with ICA fit

       - fit_transform()
       do what fit() does, but additionally return the projection onto ICA space

       - inverse_transform(A)
       inverses the decomposition, returns waveforms for an input A, using Z

       - get_params()
       returns metadata used for fits.
    """
    def __init__(self, num_components=10,
                 catalog_name='unknown',
                 whiten=True,
                 fun = 'logcosh',
                 fun_args = None,
                 max_iter = 400,
                 tol = .00001,
                 w_init = None,
                 random_state = None,
                 algorithm = 'parallel'):

        self._decomposition  = 'Fast ICA'
        self._num_components = num_components
        self._catalog_name   = catalog_name
        self._whiten         = whiten
        self._fun            = fun
        self._fun_args       = fun_args
        self._max_iter       = max_iter
        self._tol            = tol
        self._w_init         = w_init
        self._random_state   = random_state
        self._algorithm      = algorithm

        self._ICA = FastICA(n_components=self._num_components,
                             whiten       = self._whiten,
                             fun          = self._fun,
                             fun_args     = self._fun_args,
                             max_iter     = self._max_iter,
                             tol          = self._tol,
                             w_init       = self._w_init,
                             random_state = self._random_state,
                             algorithm    = self._algorithm)


    def fit(self,waveforms):
        # TODO make sure there are more columns than rows (transpose if not)
        # normalize waveforms
        self._waveforms = waveforms
        self._ICA.fit(self._waveforms)

    def fit_transform(self,waveforms):
        # TODO make sure there are more columns than rows (transpose if not)
        # normalize waveforms
        self._waveforms = waveforms
        self._A = self._ICA.fit_transform(self._waveforms)
        return self._A

    def inverse_transform(self,A):
        # convert basis back to waveforms using fit
        new_waveforms = self._ICA.inverse_transform(A)
        return new_waveforms

    def get_params(self):
        # TODO know what catalog was used! (include waveform metadata)
        params = self._ICA.get_params()
        params['num_components'] = params.pop('n_components')
        params['Decompositon'] = self._decomposition
        return params




