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

## Installation
Make sure that the python packages numpy, scipy, pandas, and patsy are already installed.
pip installer will install patsy, pandas and tabular if they aren't installed already.

    cd /path/for/code
1. Download github zip file here
2. Unzip
```
# 
cd /CCSNMultivar-master

python setup.py install
```
or

    pip install -U ccsnmultivar

# Basic Walkthrough

```
# import code
import ccsnmultivar as cm

# load waveforms
path_to_waveforms = "/path/to/waveforms.dat"
Y = cm.load_data(path_to_waveforms)

path_to_parameterfile = "/path/to/parameterfile.dat"
params = cm.load_data(path_to_parameterfile)
```
Now we need to make two objects, a Basis object and a design matrix object

- Instantiate a basis object
    Currently, I've implemented PCA with SVD, but also sklearns Kernel PCA, a 
    nonlinear basis decomposition.
```
# use a PCA basis keeping the first 10 Principal Components
pca = cm.PCA(num_components=10)
```    
- Instantiate a DesignMatrix object
```
# first, define a formula string describing how the physical parameters
#    need to be translated to the design matrix.  Say we only want to use
#    encodings of the parameters A and B (A is discrete, B is continuous)

formula = "A + B + A*B | Dev(A,omit=2), Poly(B,degree=4)"
```
The formula reads: designmatrix includes parameters A and B, and interaction 
    terms between A and B (A*B).  A is deviation encoded, A=2 is left out.
    B is encoded by a 4th degree **Chebyshev polynomial**.

Now we wrap the waveforms, basis, and design matrix (everything we need to work
with the equation (Y = XBZ^\dagger).

```
# instantiate Multivar object
M = cm.Multivar(Y,X_obj, pca)

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
B^2                      4.3            0.3
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
* crossvalidation with printed summary


* Gaussian Process reconstructions




