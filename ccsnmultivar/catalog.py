import csv as csv
import numpy as np
import warnings as warnings
import scipy as sp

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
    def __init__(self, path_to_waveforms,catalog_name=None,
                 transform_type='time',mean_subtract=False):
        self._catalog_name      = catalog_name
        self.mean_subtract      = mean_subtract
        if self._catalog_name == None:
            self._catalog_name = path_to_waveforms.split('/')[-1]

        self._path_to_waveforms = path_to_waveforms
        self.Y_dict             = self._load_waveforms(path_to_waveforms)
        self._transform         = transform_type
        self._set_metadata()

    def _set_metadata(self):
        # make metadata dictionary
        metadata = {}
        metadata['Catalog Name']            = self._catalog_name
        metadata['Number of Waveforms']     = self._n_waves
        metadata['Waveform Domain']         = self._transform
        metadata['Normalization Factor']    = self._norm_factor
        metadata['Catalog Mean Subtracted?']= self.mean_subtract
        self._metadata                      = metadata

    def _load_waveforms(self,path_to_waveforms):
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
        # how many waveforms in catalog
        self._n_waves = len(wave_list)
        # form wave_list into Y and wave_names
        self._Y_array = np.empty((len(wave_list),np.shape(wave_list)[1]-1))
        wave_names = []
        for i in np.arange(0,len(wave_list)):
            self._Y_array[i,:] = map(np.float, wave_list[i][1:])
            wave_names.append(wave_list[i][0])

        # mean subtract
        if self.mean_subtract == True:
            self._Ymean = np.mean(self._Y_array,0) # needed for predict()
            self._Y_array = self._Y_array - self._Ymean
            warnings.warn("Catalog has mean waveform subtracted")
        else:
            # just use the zero vector so Multivar.predict is happy
            self._Ymean = np.zeros(np.shape(self._Y_array)[1])

        # normalizing this way is *sort of* like normalizing each WF to have ssq = 1)
        self._norm_factor = (1./np.linalg.norm(self._Y_array,ord='fro'))*len(self._Y_array)
        self._Y_array = self._Y_array*self._norm_factor
        # now, replace Y_dict waveforms with the normalized waveforms
        Y_dict = {}
        for i in np.arange(0,self._n_waves):
            Y_dict[wave_names[i]] = self._Y_array[i,:]
        return Y_dict

    def get_catalog(self):
        """
        return dictionary of named waveforms
        """
        return self.Y_dict

    def get_transformed_Y(self):
        return self.Y_transformed

    def get_normalization(self):
        return self._norm_factor

    def get_params(self):
        return self._metadata

    def set_name(self,name):
        self._catalog_name = name

    def fit_transform(self):
        if self._transform.lower() == 'time':
            self.Y_transformed = self.Y_dict
        elif self._transform.lower() == 'fourier':
            self._fourier()
        elif self._transform.lower() == 'spectrogram':
            self._spectrogram()
        elif self._transform.lower() == 'amplitudephase':
            self._amplitudephase()
        else:
            raise ValueError("Bad input %s, transformation not defined" % typeof)
        return self.Y_transformed

    def _fourier(self):
        """ Fourier transform all waveforms in catalog, assumes fs = 16384 """
        Y_transformed = self.Y_dict.copy()
        for key in self.Y_dict.keys():
            Y_transformed[key] = sp.fft(self.Y_dict[key],n=None)
        self.Y_transformed = Y_transformed

    def _spectogram(self):
        print "not implemented yet"

    def _amplitudephase(self):
        print "not implemented yet" 






