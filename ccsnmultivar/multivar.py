import numpy as np
import scipy  as sp
import scipy.stats as stats
import re as re
import csv as csv
from tabulate import tabulate
import GWUtils


# TODO detectorsetup class.  instead of simulating waveforms in one det, expand to det network

class Multivar(object):
    """
    Multivar objects combine the settings in a Catalog, DesignMatrix and Basis object.
    It performs the regression, and has methods 

    .predict()
        Predicts waveforms given a new set of physical parameters

    .reconstruct()
        Reconstructs the waveforms in the original Catalog object

    .summary()
        Prints a summary of hypothesis tests encoded into the design matrix
        Also shows diagnostics of the fit, and metadata about the settings used to 
        produce it

    .get_waveforms()
        Returns the original core-collapse waveforms as a numpy array

    """
    def __init__(self, catalog_object, designmatrix_object, basis_object):
        # save input objects
        self._catalog_object      = catalog_object
        self._basis_object        = basis_object
        self._designmatrix_object = designmatrix_object
        # first fit objects, match waveforms to design matrix rows
        Y,X,A,col_names,row_names = self._run_arguement_object_fits()
        # save fit results
        self._X                   = X
        self._Y                   = Y
        self._A                   = A
        self._col_names           = col_names
        self._row_names           = row_names


    def deleteme_psd(self):
        psd = _resample_GetDetectorPSD('H1')
        return psd

    def _run_arguement_object_fits(self):
        """
        This function fits objects passed to Multivar, guarantees wave ordering in
        Catalog object and DesignMatrix object matches up
        """
        # run fits
        # fit catalog object
        Y_dict = self._catalog_object.get_transformed_Y()
        # fit designmatrix object
        X = self._designmatrix_object.get_matrix()
        Xcol_names = self._designmatrix_object.get_columnnames()
        row_names = self._designmatrix_object.get_rownames()
        # make sure waveforms are matched
        # loop through Xrow_names, use as key for Y_dict, populate Y_matrix
        Y = np.empty((len(row_names),len(Y_dict[Y_dict.keys()[0]])))

        # check if waveforms are complex valued.  if so, instantiate Y as complex type
        if sum(np.iscomplex(Y_dict[Y_dict.keys()[0]])):
            # then it is complex
            Y = np.empty((len(row_names),len(Y_dict[Y_dict.keys()[0]]))).astype(np.complex)

        for i in np.arange(0,len(row_names)):
            Y[i,:] = Y_dict[row_names[i]]
        # fit basis object
        A = self._basis_object.fit_transform(Y)
        return Y, X, A, Xcol_names, row_names

    def get_waveforms(self):
        """ returns waveforms as a numpy array"""
        return self._Y

    def fit(self):
        """ fit waveforms in any domain"""
        # solve for estimator of B
        n,p          = np.shape(self._X)
        self._df     = float(n - p)
        self._Cx     = np.linalg.pinv(np.dot(self._X.T,self._X))
        self._Bhat   = np.dot(np.dot(self._Cx, self._X.T), self._A)
        self._Y_rec  = self._compute_prediction(self._X)


    def _compute_prediction(self,X):
        """ compute predictions given a new X """
        A_pred = np.dot(X,self._Bhat)
        Y_pred = self._basis_object.inverse_transform(A_pred)
        return Y_pred


    def summary(self):
        """ prints results of hotellings T2 """
        transform = self._catalog_object._transform
        if   transform == 'time':
            self._hotellings_time()
        elif transform == 'fourier':
            self._hotellings_fourier()
        elif transform == 'spectrogram':
            self._hotellings_spectrogram()
        elif transform == 'amplitude':
            self._hotellings_amplitude()
        elif transform == 'phase':
            self._hotellings_phase()
        else:
            raise ValueError("Unknown catalog transformation")
        # print out to terminal
        self._make_summary_tables()


    def _hotellings_time(self):
        """ hotelling's T2 tests for time domain waveforms"""
        # get residuals
        df = self._df
        R = self._A - np.dot(self._X, self._Bhat)
        Sigma_Z = np.dot(R.T,R)*(1./df)
        # compute p-values
        T_2_list = []
        p_value_list = []
        for i in np.arange(0,self._Bhat.shape[0]):
            Bstar       = self._Bhat[i,np.arange(0,np.shape(self._A)[1])]
            lstar       = float(np.shape(Bstar)[0])
            cx          = self._Cx[i,i]
            Einv        = np.linalg.pinv(Sigma_Z)
            Zs          = Bstar/np.sqrt(cx)
            T_2         = ((df - lstar + 1.)/(df*lstar)*np.dot(np.dot(Zs,Einv),Zs.T))
            p_value     = 1. - stats.f.cdf(T_2, lstar, df - lstar + 1.)

            p_value_list.append(p_value)
            T_2_list.append(T_2)
        # save pvalue results in a table
        self._results = [['Comparison','Hotellings T^2', "p-value"]]
        for i in np.arange(0,len(self._col_names)):
            self._results.append([self._col_names[i], T_2_list[i], p_value_list[i]])


    # TODO
    def _hotellings_fourier(self):
        """ hotelling's T2 tests for fourier domain waveforms"""
        # first load in PSD and fit results, compute residual
        psd = _resample_GetDetectorPSD('H1')
        # convert psd to noise standard deviations
        sigma2 = psd/(2.*((1./16384.)**2))
        # compute residual
        df = self._df #degrees of freedom
        R = self._A - np.dot(self._X, self._Bhat)
        R = np.matrix(R)
        # residual covariance matrix
        Sigma_R = (R.H*R) * (1./df)
        # noise covariance matrix
        #raise TypeError("fourier domain hypothesis tests not implemented yet")
        Zft = np.matrix(self._basis_object._Z)
        print np.shape(Zft), np.shape(sigma2)
        Sigma_S = Zft.H*np.matrix(np.diag(sigma2))*Zft
        # add covariances (sum of multivariate normals is normal)
        Sigma_Z = np.array(Sigma_R + Sigma_S)

        T_2_list = []
        p_value_list = []
        for i in np.arange(0,self._Bhat.shape[0]):

            Bstar       = self._Bhat[i,np.arange(0,np.shape(self._A)[1])]
            lstar       = float(np.shape(Bstar)[0])
            cx          = self._Cx[i,i]
            Einv        = np.linalg.pinv(Sigma_Z)
            Zs          = Bstar/np.sqrt(cx)
            last_part   = np.array(np.matrix(Zs)*np.matrix(Einv)*np.matrix(Zs).H).real
            T_2         = ((df - lstar + 1.)/(df*lstar))*last_part
            p_value     = 1. - stats.f.cdf(T_2, 2.*lstar, 2.*(df - lstar + 1.))
            p_value_list.append(p_value)
            T_2_list.append(T_2)
        # save pvalue results in a table
        results = [['Comparison','Hotellings T^2', "p-value"]]
        for i in np.arange(0,len(self._col_names)):
            results.append([self._col_names[i], T_2_list[i], p_value_list[i]])
        self._results = results

    # # TODO
    # def _hotellings_spectrogram(self):
    #     """ hotelling's T2 tests for time domain waveforms"""

    #         p_value_list.append(p_value)
    #         T_2_list.append(T_2)
    #     # save pvalue results in a table
    #     results = [['Comparison','Hotellings T^2', "p-value"]]
    #     for i in np.arange(0,len(self._col_names)):
    #         results.append([self._col_names[i], T_2_list[i], p_value_list[i]])
    #     self._results = results

    # # TODO
    # def _hotellings_amplitudephase(self):
    #     """ hotelling's T2 tests for amplitude/phase domain waveforms"""

    #         p_value_list.append(p_value)
    #         T_2_list.append(T_2)
    #     # save pvalue results in a table
    #     results = [['Comparison','Hotellings T^2', "p-value"]]
    #     for i in np.arange(0,len(self._col_names)):
    #         results.append([self._col_names[i], T_2_list[i], p_value_list[i]])
    #     self._results = results


    def _make_summary_tables(self):
        """
        prints the summary of the regression.  It shows
        the waveform metadata, diagnostics of the fit, and results of the
        hypothesis tests for each comparison encoded in the design matrix
        """
        try:
            self._Bhat
        except:
            raise Exception("Regression hasn't been fit yet.  run .fit()")
        else:
            # check degrees of freedom
            num_pcs = self._basis_object.get_params()['num_components']
            total_dof = self._X.shape[0] - self._X.shape[1] - num_pcs
            if total_dof <= 0.0:
                raise ValueError("degrees of freedom <= 0, Hotellings T2 not defined")

            # print catalog and basis info
            cat_table = self._catalog_object.get_params().items()
            bas_table = self._basis_object.get_params().items()
            print tabulate(cat_table+bas_table,tablefmt='plain')

            # then print pvalues
            # make T^2 & pvalue table
            headers = self._results[0]
            table   = self._results[1:]
            print tabulate(table, headers, tablefmt="rst")
            print "Formula Used: %s" % self._designmatrix_object._formula
            print "Degrees of Freedom (n - p - k): %s" % str(total_dof)
            print "Condition Number of X^T*X: %.2f" % np.linalg.cond(np.dot(self._X.T, self._X))

    def reconstruct(self):
        """ return reconstructed waveforms fit by the multivar object"""
        return self._Y_rec

    # TODO print overlap table with [waveform name, parameters, overlap]
    def overlap_summary(self):
        """ print summary of reconstruction overlaps """
        olaps = self.compute_overlaps()

        # compute min, 25% 50% (median), mean, 75%, max
        table =  [["5%: ",np.percentile(olaps,5)],
                 ["25%: ",np.percentile(olaps,25)],
                 ["50%: ",np.percentile(olaps,50)],
                 ["75%: ",np.percentile(olaps,75)],
                 ["95%: ",np.percentile(olaps,95)],
                 ["  " , "  "],
                 ["Min: ",np.min(olaps)],
                 ["Mean: ",np.mean(olaps)],
                 ["Max: ",np.max(olaps)]]

        header = ["Percentile","Overlap"]
        print tabulate(table,header,tablefmt="rst")

    def compute_overlaps(self, *args):
        """ compute overlaps """
        # add Ymean back in
        Y_rec  = self._Y_rec + self._catalog_object.mean_subtract
        Y      = self.get_waveforms() + self._catalog_object.mean_subtract

        olaps = []
        psd = _resample_GetDetectorPSD('H1')
        for i in np.arange(0,Y.shape[0]):
            olaps.append(_overlap(Y[i,:], Y_rec[i,:], psd))
        return olaps


    def predict(self,param_dict):
        """#TODO predict new waveforms using multivar fit """
        encoder_dict = self._designmatrix_object.encoder
        X, col_names = self._designmatrix_object.run_encoder(param_dict, encoder_dict)
        # compute predictions
        Y_rec = self._compute_prediction(X)
        return Y_rec

def _resample_GetDetectorPSD(detector):
    """ Resamples noise PSD to 16384 Hz """
    # Load in detector noise curve for zero_det_high_p
    psd = GWUtils.GetDetectorPSD(detector)
    # psd has freq resolution = 1/3 with 6145 samples
    dF = 1./3.
    N_fd = 6145.
    # interpolate to get resolution=1 and 8192 samples
    # make the vector of frequencies for the PSD
    psd_freqs = dF*np.ones(np.shape(psd))*np.arange(0,N_fd)
    psd_interp = sp.interpolate.interp1d(psd_freqs,psd)
    psd = psd_interp(np.arange(1,2048))

    # concatenate with large frequencies out to frequency bin 8192
    large_f = psd[2046]
    psd = np.concatenate((psd,np.ones(8191-2046)*large_f))

    # make unruley amplitudes at low and high frequencies ineffective
    psd[0:11] = 1.
    psd[2000:8192] = 1.
    psd = np.concatenate((psd,psd[::-1]))
    return psd


def _inner_product(y,yr,psd):
    """
    Compute inner product between two time domain waveforms, weighted by noisecurve.
    """
    fmin = 40.
    fmax = 2000.
    fs = 16384.
    # fourier transform y and yr
    y = sp.fft(y,n=None)
    yr = sp.fft(yr,n=None)
    # compute product
    y = (1./fs)*y
    yr = (1./fs)*yr
    p = np.multiply(y,np.conjugate(yr))/psd
    p = fs*sp.ifft(p)
    product = max(abs(fs*sp.ifft(np.multiply(y,np.conjugate(yr))/psd)))
    return product

def _overlap(y,yr,psd):
    """ returns the detector noise weighted inner product """
    yyr  = _inner_product(y,yr,psd)
    yy   = _inner_product(y,y,psd)
    yryr = _inner_product(yr,yr,psd)
    olap = yyr/np.sqrt(yy*yryr)
    return olap



