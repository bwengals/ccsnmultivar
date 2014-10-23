# basis classes for ccsnmultivar
from sklearn.decomposition import FastICA, SparsePCA, DictionaryLearning
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

    def get_basis(self):
        """ Return the PCA basis vectors (Z^\dagger)"""
        return self._Z

    def _do_pca(self):
        U,s,V = np.linalg.svd(self._waveforms,full_matrices=False)
        Z = V[0:self._num_components,:]
        Z = np.array(Z.T)
        A = np.array(np.matrix(U)*np.matrix(np.diag(s)))
        A = A[:,0:self._num_components]
        self._A = A
        self._Z = Z



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
                 max_iter = 600,
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

    def get_basis(self):
        """ Return the ICA basis vectors (Z^\dagger)"""
        return self._ICA.get_mixing_matrix()




class SPCA(object):
    """
    Wrapper for sklearn package.  Performs sparse PCA

    SPCA has 5 methods:
       - fit(waveforms)
       update class instance with ICA fit

       - fit_transform()
       do what fit() does, but additionally return the projection onto ICA space

       - inverse_transform(A)
       inverses the decomposition, returns waveforms for an input A, using Z

       - get_basis()
       returns the basis vectors Z^\dagger

       - get_params()
       returns metadata used for fits.
    """
    def __init__(self, num_components=10,
                 catalog_name='unknown',
                 alpha = 0.1,
                 ridge_alpha = 0.01,
                 max_iter = 2000,
                 tol = 1e-9,
                 n_jobs = 1,
                 random_state = None):

        self._decomposition  = 'Sparse PCA'
        self._num_components = num_components
        self._catalog_name   = catalog_name
        self._alpha          = alpha
        self._ridge_alpha    = ridge_alpha
        self._n_jobs         = n_jobs
        self._max_iter       = max_iter
        self._tol            = tol
        self._random_state   = random_state

        self._SPCA = SparsePCA(n_components=self._num_components,
                              alpha        = self._alpha,
                              ridge_alpha  = self._ridge_alpha,
                              n_jobs       = self._n_jobs,
                              max_iter     = self._max_iter,
                              tol          = self._tol,
                              random_state = self._random_state)

    def fit(self,waveforms):
        # TODO make sure there are more columns than rows (transpose if not)
        # normalize waveforms
        self._waveforms = waveforms
        self._SPCA.fit(self._waveforms)

    def fit_transform(self,waveforms):
        # TODO make sure there are more columns than rows (transpose if not)
        # normalize waveforms
        self._waveforms = waveforms
        self._A = self._SPCA.fit_transform(self._waveforms)
        return self._A

    def inverse_transform(self,A):
        # convert basis back to waveforms using fit
        new_waveforms = self._SPCA.inverse_transform(A)
        return new_waveforms

    def get_params(self):
        # TODO know what catalog was used! (include waveform metadata)
        params = self._SPCA.get_params()
        params['num_components'] = params.pop('n_components')
        params['Decompositon'] = self._decomposition
        return params

    def get_basis(self):
        """ Return the SPCA basis vectors (Z^\dagger)"""
        Zt = self._SPCA.components_
        return Zt







class SC(object):
    """
    Wrapper for sklearn package.  Performs sparse coding

    Sparse Coding, or Dictionary Learning has 5 methods:
       - fit(waveforms)
       update class instance with Sparse Coding fit

       - fit_transform()
       do what fit() does, but additionally return the projection onto new basis space

       - inverse_transform(A)
       inverses the decomposition, returns waveforms for an input A, using Z^\dagger

       - get_basis()
       returns the basis vectors Z^\dagger

       - get_params()
       returns metadata used for fits.
    """
    def __init__(self, num_components=10,
                 catalog_name='unknown',
                 alpha = 0.001,
                 transform_alpha = 0.01,
                 max_iter = 2000,
                 tol = 1e-9,
                 n_jobs = 1,
                 verbose = True,
                 random_state = None):

        self._decomposition   = 'Sparse Coding'
        self._num_components  = num_components
        self._catalog_name    = catalog_name
        self._alpha           = alpha
        self._transform_alpha = 0.001
        self._n_jobs          = n_jobs
        self._random_state    = random_state

        self._DL = DictionaryLearning(n_components=self._num_components,
                              alpha           = self._alpha,
                              transform_alpha = self._transform_alpha,
                              n_jobs          = self._n_jobs,
                              verbose         = verbose,
                              random_state    = self._random_state)

    def fit(self,waveforms):
        # TODO make sure there are more columns than rows (transpose if not)
        # normalize waveforms
        self._waveforms = waveforms
        self._DL.fit(self._waveforms)

    def fit_transform(self,waveforms):
        # TODO make sure there are more columns than rows (transpose if not)
        # normalize waveforms
        self._waveforms = waveforms
        self._A = self._DL.fit_transform(self._waveforms)
        return self._A

    def inverse_transform(self,A):
        # convert basis back to waveforms using fit
        new_waveforms = self._DL.inverse_transform(A)
        return new_waveforms

    def get_params(self):
        # TODO know what catalog was used! (include waveform metadata)
        params = self._DL.get_params()
        params['num_components'] = params.pop('n_components')
        params['Decompositon'] = self._decomposition
        return params

    def get_basis(self):
        """ Return the SPCA basis vectors (Z^\dagger)"""
        return self._DL.components_
