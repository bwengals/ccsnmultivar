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


    def _run_arguement_object_fits(self):
        """
        This function fits objects passed to Multivar, guarantees wave ordering in
        Catalog object and DesignMatrix object matches up
        """
        # run fits
        # fit catalog object
        Y_dict = self._catalog_object.fit_transform()
        # fit designmatrix object
        X,Xcol_names,row_names = self._designmatrix_object.fit_transform()
        # make sure waveforms are matched
        # loop through Xrow_names, use as key for Y_dict, populate Y_matrix
        Y = np.empty((len(row_names),len(Y_dict[Y_dict.keys()[0]])))

        # check if waveforms are complex valued.  if so, instantiate Y as complex type
        if sum(np.iscomplex(Y_dict[Y_dict.keys()[0]])):
            # then it is complex
            Y.astype(np.complex,copy=False)

        for i in np.arange(0,len(row_names)):
            Y[i,:] = Y_dict[row_names[i]]
        # fit basis object
        A = self._basis_object.fit_transform(Y)
        return Y, X, A, Xcol_names, row_names


    def _get_waveforms(self):
        """ return the (unnormalized) waveforms as a numpy array """
        return self._Y

    def get_waveforms(self):
        """ returns renormalized waveforms as a numpy array"""
        Y = self._Y
        Y = Y*(1./self._catalog_object.get_normalization())
        return Y

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
        elif transform == 'amplitudephase':
            self._hotellings_amplitudephase()
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
        raise TypeError("fourier domain hypothesis tests not implemented yet")
        # TODO NEED Z_ft, basis Z of PCs in fourier domain.  I dont think sklearn will
        # be happy about this.  Can you even get PCs from nonlinear types of PCA?
        Sigma_S = Zft.H*np.matrix(diag(sigma2))*Zft
        # add covariances (sum of multivariate normals is normal)
        Sigma_Z = Sigma_R + Sigma_S

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
            print "Residual Sum-of-Squares: %.2f" % np.sum(np.sum(np.square(self._Y - self._Y_rec)))


    def _renormalize_catalog(self,Y):
        """ undo normalization constant and mean subtraction applied by Catalog object"""
        Y = Y + self._catalog_object._Ymean
        Y = Y*(1./self._catalog_object.get_normalization())
        return Y


    def reconstruct(self):
        """ return reconstructed waveforms fit by the multivar object"""
        Y_rec = self._renormalize_catalog(self._Y_rec)
        return Y_rec

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

    def leave_one_out_crossvalidation(self):
        """ print summary of a leave-one-out crossvalidation run """ 
        print "not implemented yet"


    def compute_overlaps(self, *args):
        """ compute overlaps """
        if len(args) != 2:
            # renormalize the catalog and predictions and add mean back in
            Y_rec  = self._renormalize_catalog(self._Y_rec)
            Y      = self._renormalize_catalog(self._get_waveforms())
        else:
            # then args[0] and args[1] should be numpy catalog arrays
            Y_rec  = self._renormalize_catalog(args[0])
            Y      = self._renormalize_catalog(args[1])

        olaps = []
        psd = _resample_GetDetectorPSD('H1')
        for i in np.arange(0,Y.shape[0]):
            olaps.append(_overlap(Y[i,:], Y_rec[i,:], psd))
        return olaps


    # TODO code works, idea doesnt...
    def compute_effect_size(self,col_name_list):
        """ compute effect size summary
        Define effect size of X column called (a) as the average catalog
            overlap difference between the True Waveforms and the
            average when column a is removed.  Similarly for multiple columns
        """
        # set behavior is weird if just one column is given as a string
        if type(col_name_list) == str:
            # make it into a list of one element
            col_name_list = [col_name_list]

        # renormalize the catalog and predictions and add mean back in
        Y_rec  = self._renormalize_catalog(self._Y_rec)
        Y      = self._renormalize_catalog(self._get_waveforms())

        # overlaps for full X reconstructions
        fullX_olaps = self.compute_overlaps()

        # make X_reduced by zeroing out the columns of col_name_list
        col_names = self._designmatrix_object.get_columnnames()
        zero_idx = []
        X_reduced = self._designmatrix_object.get_X().copy()
        for item in col_name_list:
            idx = col_names.index(item)
            X_reduced[:,idx] = 0.0

        # compute overlaps for the X_reduced fit
        Y_reduced = self._compute_prediction(X_reduced)
        reducedX_olaps = self.compute_overlaps(Y,Y_reduced)

        efx_size = np.mean(fullX_olaps) - np.mean(reducedX_olaps)
        return efx_size


    def predict(self,param_dict):
        """#TODO predict new waveforms using multivar fit """
        X, col_names = _parse_with_encoder_dict(param_dict, self._designmatrix_object)
        # compute predictions
        Y_rec = self._compute_prediction(X)
        Y_rec = self._renormalize_catalog(Y_rec)
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
    psd[0:11] = 200.
    psd[2000:8192] = 200.
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

def _logL(y_actual, y_template,psd):
    """ compute the noise weighted log likelihood of a template
    waveform, given the true waveform (ssq residual assumed gaussian
    detector noise distributed)
    CODE BY SARAH GOSSAN FROM GWUtils.py
    """
    # frequency spacing
    dF = 1.
    # fourier transform y_actual and y_template
    y_actual = sp.fft(y_actual,n=None)
    y_template = sp.fft(y_template,n=None)
    # calculate log-likelihood
    chi2 = np.sum(pow(abs(y_actual - y_template), 2.)/psd)
    logL = -2.*dF*chi2
    return logL


def _parse_with_encoder_dict(param_dict,designmatrix_object):
    """ parse new param_dict using the encoder_dict from 
    design matrix object """
    encoder_dict = designmatrix_object._encoder_dict
    inter_list = designmatrix_object._inter_list
    # predict first with non 'twoway' or 'threeway' functions
    X_dict = {}
    Xcol_dict = {}
    # put each column of X in Xbycol_dict
    Xbycol_dict = {}
    for key in encoder_dict:
        if (key != 'twoway') and (key != 'threeway'):
            encoder = encoder_dict[key]
            param_values = param_dict[key]
            Xsub,names = encoder(key,param_values)
            X_dict[key] = Xsub
            Xcol_dict[key] = names
            for i in np.arange(0,len(names)):
                Xbycol_dict[names[i]] = Xsub[:,i]

    # now do interactions
    for interaction in inter_list:
        if len(interaction) == 2:
            encoder = encoder_dict['twoway']
            param_name1 = interaction[0]
            param_name2 = interaction[1]
            col_names1 = Xcol_dict[param_name1]
            col_names2 = Xcol_dict[param_name2]
            X_int, names = encoder(param_name1,param_name2, \
                                   col_names1,col_names2, X_dict)

            # put columns into Xbycol_dict
            for i in np.arange(0,len(names)):
                Xbycol_dict[names[i]] = X_int[:,i]

        elif len(interaction) == 3:
            encoder = encoder_dict['threeway']
            param_name1 = interaction[0]
            param_name2 = interaction[1]
            param_name3 = interaction[2]
            col_names1 = Xcol_dict[param_name1]
            col_names2 = Xcol_dict[param_name2]
            col_names3 = Xcol_dict[param_name3]
            X_int, names = encoder(param_name1,param_name2,param_name3, \
                                   col_names1,col_names2,col_names3, X_dict)

            # put columns into Xbycol_dict
            for i in np.arange(0,len(names)):
                Xbycol_dict[names[i]] = X_int[:,i]

        else:
            print "this shouldnt have happenend"
    # get original design matrix column ordering
    col_names = designmatrix_object.get_columnnames()
    X = []
    # col_names[1:] to exclude "Intercept"
    for name in col_names[1:]:
        X.append(Xbycol_dict[name])
    # always add intercept column last
    X.insert(0,np.ones(np.shape(X[0])))
    # final design matrix
    X = np.vstack(X).T
    return X, col_names
