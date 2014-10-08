# hidden imports of dependency modules
import numpy as np
import pandas as pd
import patsy  as pt
import scipy  as sp
import scipy.stats as stats
import re as re
import warnings as warnings
import csv as csv
from tabulate import tabulate

# TODO
#waveform class.  each waveform is named (matching design matrix so no screwups can happen)
#   each waveform set can be given metadata (what catalog)
#   different transforms: fft, amp/phase, spectrograms
#detectorsetup class.  instead of simulating waveforms in one det, expand to det network

class Multivar(object):
    """
    When multivar objects are instantiated, the waveforms are loaded and 
    """
    def __init__(self, waveforms, designmatrix_object, basis_object):
        self._Y             = waveforms
        self._A             = basis_object.fit_transform(self._Y)
        self._X             = designmatrix_object.get_X()
        self._detector      = 'H1'
        # other useful parameters
        self._formula       = designmatrix_object.get_formula()
        self._column_names  = designmatrix_object.get_columnnames()
        self._dof           = designmatrix_object.get_dof()
        # load basis_object
        self._num_components = basis_object.num_components
        # give basis object, and desigm matrix object into DM class self
        self._basis_object = basis_object
        self._designmatrix_object = designmatrix_object

    def fit(self, type_of_fit):
        # TODO: if waveforms complex valued, fit freq domain, otherwise do time domain
        if 'time' in type_of_fit:
            self._fit_time_domain()
        elif 'freq' in type_of_fit:
            self._fit_freq_domain()
        elif 'no' in type_of_fit:
            self._fit_no_basis()
        else:
            raise Exception("Unknown fit-type key word")

    def _fit_time_domain(self):
        # solve for estimators of B and Sigma_Z
        X      = np.matrix(self._X)
        Y      = np.matrix(self._Y)
        A      = np.matrix(self._A)
        n,p    = np.shape(X)
        df     = float(self._dof)
        Cx     = np.linalg.pinv(X.T*X)

        Bhat = Cx*X.T*A
        R = A - X*Bhat
        Sigma_Z = R.T*R*(1./df)
        self._Bhat = np.array(Bhat)
        self._sumofsquares = np.sum(np.sum(np.square(R)))

        columns = np.arange(0,self._num_components)
        T_2_list = []
        p_value_list = []
        z_score_list = []
        for row in range(0,p):
            Bstar       = Bhat[row,columns]
            pstar,lstar = np.shape(Bstar)
            cx          = float(Cx[row,row])
            Einv        = np.linalg.pinv(Sigma_Z)
            Zs          = Bstar/np.sqrt(cx)
            T_2         = float(((df - lstar + 1.)/(df*lstar))*Zs*Einv*Zs.T)
            p_value     = 1. - stats.f.cdf(T_2,float(lstar), df - float(lstar) + 1.)
            z_score     = stats.f.ppf(q=1.-p_value, dfn=float(lstar) , \
                                         dfd =df-float(lstar) + 1.)
            z_score     = stats.norm.ppf(1-p_value)
            p_value_list.append(p_value)
            z_score_list.append(z_score)

            T_2_list.append(T_2)
        results = [['Comparison','Hotellings T^2', "p-value", "Sigma"]]
        for i in np.arange(0,len(self._column_names)):
            results.append([self._column_names[i], T_2_list[i], 
                            p_value_list[i],z_score_list[i]] )
        self._results = results

    def _fit_freq_domain(self):
        print "stuff"

    def _fit_no_basis(self):
        print "stuff"


    def summary(self):
        if self._dof < 1:
            raise Exception("Number of waveforms < number of X columns, Hotellings T2 undefined")
        try:
            self._results
        except:
            raise Exception("Regression hasn't been fit yet.  run .fit()")
            
        else:
            basis_info = [["Catalog Used:",str(self._basis_object.catalog_name)]]
            basis_info.append(["Basis Type Used:",str(self._basis_object.decomposition)])
            # collect table of fit metadata
            basis_params = self._basis_object.get_params()
            # convert dictionary to list
            for key,value in basis_params.iteritems():
                basis_info.append([key,value])
            # print basis information first, then print pvalue information
            print tabulate(basis_info,tablefmt='plain')

            # make T^2 & pvalue table
            headers = self._results[0]
            table   = self._results[1:]
            print tabulate(table, headers, tablefmt="rst")
            X = np.matrix(self._X)
            # print null Hypothesis
            print "Null Hypothesis: B_{row i} = B_{row i+1} = ... = B_{row N} = 0"
            print "Probablitity of observing our results, given the Null is true = p-value"
            # print condition number of X.T*X
            cond_num = np.linalg.cond(X.T*X)
            print "Condition Number of X^T*X: " + str(cond_num)
            print "Sum-of-Squares of Fit Residual: %s" % self._sumofsquares

    def predict(self,*arg):
        if len(arg) == 0:
            # then predict self._Y
            Xpred = self._X
            # make sure to just save the phys param names we want
            formula_dict = parse_formula(self._formula)[0]
            column_names = formula_dict.keys()
            # save prediction parameters
            self._prediction_params = self._parameter_df[column_names]
        elif len(arg) == 1:
            # then predict the extra dataframe parameters
            new_df = arg[0]
            #concat new_df underneath old df
            # keys are column names of new_df and old_df that we want to keep
            formula_dict = parse_formula(self._formula)[0]
            column_names = formula_dict.keys()
            new_df = new_df[column_names]
            old_df = self._parameter_df[column_names]
            # save prediction parameters
            self._prediction_params = new_df[column_names]
            # add column with label for new or old
            new_df['Predict'] = True
            old_df['Predict'] = False
            # concat vertically
            big_df = pd.concat([old_df,new_df])
            # convert to design matrix
            Xdf = formula_to_design_matrix(self._formula,big_df)
            # which rows of Xdf to keep
            keep = np.array(big_df[big_df['Predict'] == True].index)
            Xpred = np.array(Xdf)[keep,:]
        else:
            raise Exception("Only one or two arguements allowed, not more")
        # make predictions
        if self._Bhat == 'Null':
            raise Exception("Regression coefficients haven't been fit yet")
        else:
            Xpred = np.matrix(Xpred)
            Bhat = np.matrix(self._Bhat)
            Z = np.matrix(self._Z)
            Y_new = np.array(Xpred*Bhat*Z.T)
        self._Y_predicted = Y_new


    def load_detector(self, detector='H1'):
        self._detector = detector
        self._PSD = set_pst()

def _set_psd(self):
    detector = self._detector
    # Load in detector noise curve for zero_det_high_p
    psd = _GetDetectorPSD(detector)
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

    # fix unruly values
    psd[0:11] = 200.
    psd[2000:8192] = 200.
    psd = np.concatenate((psd,psd[::-1])) 
    self._psd = psd


def _GetDetectorPSD(detectorName, LIGO3FLAG=0):
    """ GetDetectorPSD - function to return the PSD of the detector described by
    detectorName.

    detectorName - Name of required GW IFO.  Can be 'H1', 'L1', 'V1', 'I1',
                       or 'K1'.
    LIGO3FLAG - Set to 1 to use LIGO3 PSD instead of aLIGO PSD for H1, L1
        and I1 detectors.

    Returns PSD - the noise PSD for the detector described by detectorName.

        Sarah Gossan 2012. Last updated 02/18/14. """

    # H1, L1 or I1 - use aLIGO ZERODET_HIGHP configuration (fLow = 9Hz)
    # Noise curve is super complicated so just load in PSD file for now
    if detectorName == 'H1' or detectorName == 'L1' or detectorName == 'I1':
        if LIGO3FLAG:
            # Read in PSD
            PSD = np.loadtxt('LIGO3_PSD.txt')
        else:
            # Read in PSD
            PSD = np.loadtxt('ZERO_DET_high_P_PSD.txt')
        # Only want second column of file
        PSD = PSD[:,1]
        # V1 - use analytical expression for AdVirgo (fLow = 10Hz)
    elif detectorName == 'V1':
        # Use analytical expression from arXiv:1202.4031v2
        x = np.linspace(0,N_fd-1,num=N_fd)*dF/300.
        x[0] = x[1] # Not going to use f=0Hz component anyway, but this stops 
        # the log fn complaining
        x = np.log(x)
        xSq = x*x
        asd = 1.259e-24*(0.07*np.exp(-0.142 - 1.437*x + 0.407*xSq) + \
                         3.1*np.exp(-0.466 - 1.043*x - 0.548*xSq) + \
                         0.4*np.exp(-0.304 + 2.896*x - 0.293*xSq) + \
                         0.09*np.exp(1.466 + 3.722*x - 0.984*xSq))
        PSD = asd**2 
        # K1 - use analytical expression for KAGRA (fLow = 10Hz) 
    elif detectorName == 'K1':
        # Use analytical expression from arXiv:1202.4031v2
        x = np.linspace(0,N_fd-1,num=N_fd)*dF/300.
        x[0] = x[1] # Not going to use f=0Hz component anyway, but this stops 
                    # the log fn complaining
        x = np.log(x)
        xSq = x*x
        asd = 6.499e-25*(9.72e-9*np.exp(-1.43 - 9.88*x - 0.23*xSq) + \
                 1.17*np.exp(0.14 - 3.10*x - 0.26*xSq) + \
             1.70*np.exp(0.14 + 1.09*x - 0.013*xSq) + \
             1.25*np.exp(0.071 + 2.83*x - 4.91*xSq))
        PSD = asd**2
    return PSD




