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
Using the code happens in five steps:

1. Instantiate a Catalog object
2. Instantiate a Basis object.
3. Instantiate a DesignMatrix object.
4. Wrapping them in a Multivar object.
5. Analysis using the Multivar object's methods.

```python
# import code
import ccsnmultivar as cc

# load waveforms
path_to_waveforms = "/path/to/waveforms.dat"
# the Abdikamalov waveform file is called "Abdika13_waveforms.csv"
# instantiate Catalog object
Y = cc.Catalog(path_to_waveforms)
```

Note that Abdikamalov et al's 2013 waveform catalog and parameter file are included
in the Example_Waveforms directory of the GitHub repo as an example of how to format
the raw files for input.  To access these for the walkthrough, look at the right side
of the GitHub page, there is a toolbar with a Download button.  Download, then unzip.  

Now we need to make two objects, a Basis object and a DesignMatrix object

First we instantiate a Basis object.  Currently, there are three available types of 
Basis objects.
 
1. PCA - using the Singular Value Decompostion (SVD)
2. KCPA - Kernel PCA.  A wrapper for sklearns [Kernel PCA](http://scikit-learn.org/stable/modules/generated/sklearn.decomposition.KernelPCA.html#sklearn.decomposition.KernelPCA) 
3. ICA - Independent Component Ananlysis.  A wrapper for skearns [FastICA](http://scikit-learn.org/stable/modules/generated/sklearn.decomposition.FastICA.html)

```python
# use a PCA basis keeping the first 10 Principal Components
pca = cc.PCA(num_components=10)
```    
Next we instantiate a DesignMatrix object.

```python
# first, define a formula string describing how the physical parameters
#    need to be translated to the design matrix.  Say we only want to use
#    encodings of the parameters A and B (A is discrete, B is continuous)

formula = "A + B + A*B | Dum(A,omit=2), Poly(B,degree=4)"
```

The formula contains 5 peices of information that determine how the design matrix is 
encoded.  Reading the formula from left to right:

1. Include columns for the physical parameter named "A".
2. Include columns for the physical parameter named "B".
3. Include columns for interaction terms between parameters "A" and "B".  
The "|" character seperates instructions for *what* goes into the design matrix from 
*how* it goes in.
4. Use a dummy variable encoding on parameter "A".  One value of "A" needs to be used as a
reference in a dummy variable encoding, we chose value "2".
4. Use a chebyshev polynomial encoding on parameter "B".  Fit "B" with a 4th degree polynomial.

Now we instantiate the DesignMatrix object with two arguments: the formula, and the
path to the parameter file.
```python

# note that the provided Abdikamalov+ parameterfile is called "Abdika13_params.csv"
path_to_parameterfile = "/path/to/parameterfile.csv"

# note that we dont need to load the paramfile, just supply the path.
X = cc.DesignMatrix(path_to_parameterfile, formula)
```

Now with the waveforms in the Catalog object Y, the Basis object pca, and DesignMatrix object 
X on hand, we instantiate a Multivar object with these three arguements.

```python
# instantiate Multivar object
M = cc.Multivar(Y,X, pca)

```

This makes it easy to create many different Catalog, DesignMatrix, Basis, and Multivar
objects to test different fits and parameter influences very quickly.

```python
# now we do a fit to time domain waveforms (solve for B)
M.fit('time')

# print/save summary of the hypothesis tests, metadata, and other
# facts defined by the particular formula and basis used to make M.

M.summary()

Waveform Domain      time
Number of Waveforms  92
Catalog Name         Example_Catalogs/Abdika13_waveforms.csv
Decomposition        PCA
num_components       30
=================  ================  =========
Comparison           Hotellings T^2    p-value
=================  ================  ========= 
A:[1 - 2]                3.59           0.1
A:[3 - 2]                6.88           0.52
A:[4 - 2]                0.4            0.99
B^1                     10.6            0.012
B^2                     50.9            0.000
B^3                      4.3            0.3
A:[1 - 2]*B^1            1.1            0.78
A:[1 - 2]*B^2           12.99           0.044
    .                     .              .
    .                     .              .
    .                     .              .



# we can view the  waveform reconstructions with the Multivar method .reconstruct()
Y_reconstructed = M.reconstruct()
```

## Dependencies
* numpy
* scipy
* scikits-learn
* tabulate

## Planned
* More than one detector (the GW detector network)
* Catalog objects
  - amplitude/phase decomposition, spectrograms, metadata
* other PC basis methods 
  - kmeans, fix KernelPCA, etc.
* other design matrix fitting methods 
  - splines, rbfs, etc.
* cross validation with printed summary()
* Gaussian Process (or other interpolation/machine learning method) reconstructions




