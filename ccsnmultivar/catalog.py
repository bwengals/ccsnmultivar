import csv as csv
import numpy as np
import warnings as warnings
import scipy as sp
from scipy import signal

class Catalog(object):
    """
    Catalog Object
    Holds waveforms, catalog meta data.  Does transformations.  When Catalog object is
         instantiated with a path, it loads waveforms at that directory.
    Catalog objects have 4 methods:

       - fit_transform(typeof)
         typeof equals one of ['time','fourier','spectrogram','amplitudephase']
         transforms raw waveforms (in case of 'time' does no transformation)

       - get_params()
         returns dictionary of catalog metadata

       - get_normalization()
         returns the normalization constant used

       - get_catalog()
         return a dictionary of waveforms, where the name is the key

       - set_name(name)
         give catalog name (default is the waveform file name)
    """
    def __init__(self, Y_dict, metadata, catalog_name="Unspecified",
                 transform_type='time'):
        self._catalog_name           = catalog_name
        self._transform              = transform_type
        self.Y_dict                  = Y_dict
        metadata['Catalog Name']     = self._catalog_name
        metadata['Waveform Domain']  = self._transform
        self._metadata               = metadata

    def get_catalog(self):
        """
        return dictionary of named waveforms
        """
        return self.Y_dict

    def get_transformed_Y(self):
        """ return dictionary of named and transformed waveforms """
        return self.Y_transformed

    def get_params(self):
        return self._metadata

    def set_name(self,name):
        self._catalog_name = name

    def fit_transform(self):
        if self._transform.lower() == 'time':
            self._time()
        elif self._transform.lower() == 'fourier':
            self._fourier()
        elif self._transform.lower() == 'spectrogram':
            self._spectrogram()
        elif self._transform.lower() == 'amplitude':
            self._amplitudephase()
        else:
            raise ValueError("Bad input %s, transformation not defined" % typeof)

    def _time(self):
        """ time domain (no transform) """
        self.Y_transformed = Y_dict

    def _fourier(self):
        """ 1 side Fourier transform and scale by dt all waveforms in catalog """

        freq_bin_upper = 2000
        freq_bin_lower = 40
        fs = self._metadata['fs']
        Y_transformed = {}
        for key in self.Y_dict.keys():
            # normalize by fs, bins have units strain/Hz
            fs = self._metadata["fs"]
            Y_transformed[key] = (1./fs)*np.fft.fft(self.Y_dict[key])[freq_bin_lower:freq_bin_upper]
        self.Y_transformed = Y_transformed
        # add in transform metadata
        self._metadata["dF"] = 1./self._metadata["T"]

        # because we are in the fourier domain, we will need the psd
        self.psd = load_psd()[freq_bin_lower:freq_bin_upper]
        dF = self._metadata['dF']
        self.sigma = convert_psd_to_sigma(self.psd, dF)

    def _spectogram(self):
        print "not implemented yet"

    def _amplitude(self):
        print "not implemented yet" 

    def _phase(self):
        print "not implemented yet" 

    def mcmc_hook(self):
        dF = self._metadata["dF"]
        psd = self.psd
        noise = sample_of_fd_noise(psd, dF)
        return noise, psd

def sample_of_fd_noise(psd, dF):
    N = len(psd)
    n = np.zeros((N),complex)
    n_real = np.random.randn(N)*np.sqrt(psd/(4.*dF))
    n_imag = np.random.randn(N)*np.sqrt(psd/(4.*dF))
    n = n_real + 1j*n_imag
    return n

def load_psd():
    """ Resamples advLIGO noise PSD to 4096 Hz """
    # psd has freq resolution = 1/3 with 6145 samples
    psd = np.loadtxt("ZERO_DET_high_P_PSD.txt")[:,1]
    down_factor = 3
    pad_size = int(np.ceil(float(psd.size)/down_factor)*down_factor - psd.size)
    psd_padded = np.append(psd, np.zeros(pad_size)*np.NaN)
    psd = sp.nanmean(psd_padded.reshape(-1,down_factor), axis=1)
    # now dF = 1
    # length of psd = 2048
    return psd

def convert_psd_to_sigma(psd, dF):
    sigma = np.sqrt(psd/(4.*dF))
    return sigma

def resample_waveforms(Y_dict):
    """
        INPUT: - Dictionary of waveforms loaded from text file
               - ALSO, Dictionary of timeseries indexed by name 

        OUTPUT:
               - Y_dict: resampled, normalized waveforms, ready for FFT
                    - new parameters:
                        + fs = 4096
                        + time duration: 1 second
               - metadata: dictionary of sampling information

    """
    # right now, Y_dict, using my ccsnmultivar github waveform sets
    for key in Y_dict.keys():
        Y_dict[key] = signal.decimate(Y_dict[key], 4, ftype='fir')
    metadata = {}
    metadata["fs"] = 4096. # in Hz
    metadata["T"] = 1. # in seconds
    metadata["source distance"] = 10. # in kpc
    return Y_dict, metadata



def load_waveforms_from_file(path_to_waveforms):
    """
     Waveforms must be:
    - in a .csv or .dat file
    - each row is a waveform, values comma separated
    - each waveform is preprocessed
        -- aligned to core bounce
        -- all have same sampling frequencies
        -- all have same number of time samples
    - the name of the waveform is the first column
    - subsequent columns are the waveform time samples
    """
    wave_list = list(csv.reader(open(path_to_waveforms,"rb")))
    # delete empty elements (if any)
    wave_list = [x for x in wave_list if x != []]
    Y_dict = {}
    for i in np.arange(0,len(wave_list)):
        Y_dict[wave_list[i][0]] = map(np.float, wave_list[i][1:])
    return Y_dict







