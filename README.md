# ccsnmultivar companion code


This [Python](http://www.python.org/) module aids the analysis of core-collapse supernova gravitational waves.  It is the companion code for [this paper](http://arxiv.org/abs/1406.1164).

* **Multivariate Regression** of Fourier Transformed or Time Domain waveforms
* Hypothesis testing for measuring the influence of physical parameters
* Optionally incorporate additional uncertainty due to detector noise
* **Waveform reconstruction** and prediction
* Fit and analyze waveforms in Principal Component space
* Includes the [Abdikamalov et. al.](http://arxiv.org/abs/1311.3678) for use examples 

## Details
* Uses pandas for data handling
* A simplified formula language specific to this domain
* [Documentation](http://ccsnmultivar.readthedocs.org/en/latest/)


## Installation
Make sure that the python packages numpy, scipy, pandas, and patsy are already installed.
pip installer will install patsy, pandas and tabular if they aren't installed already.

    cd /path/to/ccsnmultivar

1. Download github zip file here
2. Unzip

```python
# 
cd /CCSNMultivar-master

python setup.py install
```
or

    pip install -U ccsnmultivar

## Basic Walkthrough
Using the code happens in four steps:

1. Making a Basis object.
2. Making a DesignMatrix object.
3. Wrapping them in a Multivar object.
4. Analysis using the Multivar object's methods.

```python
# import code
import ccsnmultivar as cm

# load waveforms
path_to_waveforms = "/path/to/waveforms.dat"
Y = cm.load_data(path_to_waveforms)
```

Note that Abdikamalov et al's 2014 waveform catalog and parameter file are included
in the ProcessedWaveforms directory as an example of how to format the raw files for input.  

Now we need to make two objects, a Basis object and a design matrix object

First we instantiate a Basis object.  Currently, there are two Basis objects.
 
1. PCA - using the Singular Value Decompostion (SVD)
2. KCPA - Kernel PCA.  A wrapper for sklearns [Kernel PCA](http://scikit-learn.org/stable/modules/generated/sklearn.decomposition.KernelPCA.html#sklearn.decomposition.KernelPCA) 

```python
# use a PCA basis keeping the first 10 Principal Components
pca = cm.PCA(num_components=10)
```    
Next we instantiate a DesignMatrix object.

```python
# first, define a formula string describing how the physical parameters
#    need to be translated to the design matrix.  Say we only want to use
#    encodings of the parameters A and B (A is discrete, B is continuous)

formula = "A + B + A*B | Dev(A,omit=2), Poly(B,degree=4)"
```

The formula contains 5 peices of information that determine how the design matrix is 
encoded.  Reading the formula from left to right:

1. Include columns for the physical parameter named "A".
2. Include columns for the physical parameter named "B".
3. Include columns for interaction terms between parameters "A" and "B"
The "|" character seperates instructions for *what* goes into the design matrix from 
*how* it goes in.
4. Use a deviation encoding on parameter "A".  One value of "A" needs to be omitted 
from the design matrix in a deviation encoding, this value is "2".
4. Use a chebyshev polynomial encoding on parameter "B".  Fit "B" with a 4th degree polynomial.

Now we instantiate the DesignMatrix object with two arguements, the formula, and the
path to the parameter file.
```python

path_to_parameterfile = "/path/to/parameterfile.dat"
params = cm.load_data(path_to_parameterfile)

X_obj = cm.DesignMatrix(path_to_parameterfile, formula)
```

Now with the waveforms in numpy array Y, the Basis object, and design matrix object on hand,
we instantiate a multivar object with these three arguements.

```python
# instantiate Multivar object
M = cm.Multivar(Y,X_obj, pca)

```

This makes it easy to create many different DesignMatrix, Basis, and Multivar objects to test 
different fits and parameter influences quickly.

```python
# now we do a fit to time domain waveforms (solve for B)
M.fit('time')

# print/save summary of the hypothesis tests, fit metadata, and other
# facts defined by the particular formula and basis used to make M.

M.summary()
=================  ================  =========
Comparison           Hotellings T^2    p-value
=================  ================  ========= 
[A1 - A2]                3.59           0.1
[A3 - A2]                6.88           0.52
[A4 - A2]                0.4            0.99
B^1                     10.6            0.012
B^2                     50.9            0.000
B^3                      4.3            0.3
[A1 - A2]*B^1            1.1            0.78
[A1 - A2]*B^2           12.99           0.044
    .                     .              .
    .                     .              .
    .                     .              .
    
# reconstruct waveforms from the original set of physical parameters (param_df)
M.predict()

# extract waveform reconstructions from M
Y_rec = M.Y_predicted

# predict waveforms given a new set of physical parameters
M.predict(new_df)

# extract waveform predictions from M
Y_pred = M.Y_predicted
```

## Dependencies
* numpy
* scipy
* pandas
* patsy
* scikits-learn
* tabulate

## Planned
* More than one detector (the GW detector network)
* waveform objects
  - amplitude/phase decomposition, spectrograms, metadata
* other PC basis methods 
  - ICA, kmeans, fix KernelPCA, etc.
* other design matrix fitting methods 
  - splines, rbfs, etc.
* crossvalidation with printed summary()
* Gaussian Process (or other interpolation/machine learning method) reconstructions




