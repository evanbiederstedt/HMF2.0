#
# Python interface to run .R and .Rdata files from Tsalmantza, Hogg 2011
#
import math
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import pyfits as pf
import astropy as ap
import os

os.chdir('/Users/evanbiederstedt/desktop/HMF2.0')
os.getcwd()

# Files in HMF2.0 subdirectory:
# functions.R
# iterations20f.R
# spectra.Rdata
# noise.Rdata
# mask.Rdata
# lamda_tel.Rdata
#


###
### Tsalmantza, Hogg 2011
### http://arxiv.org/pdf/1201.3370.pdf
### http://arxiv.org/pdf/1106.1180.pdf
###
###
### CODE FILES:
###
### functions.R
### interations20f.R
###
###
### DATA FILES:
###
### noise.Rdata
### mask.Rdata
### lamda_tel.Rdata  NOTE: not lamBda_tel.Rdata
### spectra.Rdata
###
###
###
### noise.Rdata:     length(noise)=20280000,   median(nosie)=1.217133e-17
### lamda_tel.Rdata: length(lamda_tel)=4056,   from 3580.964 to 9109.615
### spectra.Rdata:   length(spectra)=20280000, median(spectra)=1.490924e-16
### mask.Rdata:      length(mask)=20280000,
### mask 1.000, length()=1052522
### mask 0.000, length()=18350960
### mask other, length()=876518
###
### noise:   dim()= 4056 4056
### spectra: dim()= 4056 4056
### mask:    dim()= 4056 4056
###
### outputs:
### e.g. a, at, g = dim() =  15 by something
###
### p. 17, All the spectra were moved to the rest-frame by keeping the energy
### constant in each spectral bin while relabeling the wavelength axis.
###
### FINAL WAVELENGTH COVERAGE:  3580.964 < λ < 9109.615 ˚A for Main Galaxies,
### corresponding to 4056 pixels
###

##
##  K components g_{jk} and coefficients a_{ik}, over all the N spectra and the M pixels of the training set
##
##  [N,K] coefficients a_{ik}
##  [K,M] parameters g_{kj} in basis spectra
##  data:  [N,M]
##  innvar: either a scalar or else [N,M]
##  \epsilon: scalar that sets strength of smoothing; regularization equals smoothness prior
##  Adjustments: Impose (1) regularization and (2) non-negativity a_{ik} > 0, g_{kj} > 0
##  All the negative fluxes and their corresponding errors in observed spectra are set to zero
##
##  N: number of M dimensional data points
##  M: number of wavelengths
##  K: number of components permitted
##
##  HMF optimizes probabilistic objective function, i.e. minimize chi**2,
##  returns K components whose linear combination best reproduces the entire training set of the observations,
##  given the reported observational uncertainty variances, \sigma^2_ij
##
##  CONTRAST TO PCA RESULTS
##  (1) HMF produces basis functions that fit the real data much better
##  (2) HMF is also able to extract more information from the training set (more data points, more data dimensions)
##
##  For SDSS double redshift objects, there are N observed spectra i, \vec{f}_i such that
##  f_{ij}=f_{\lambda, i}(\lambda_j)=\Sum^K_{k=1} a_{ik} g_{k}(\lambda_j)+n_{ij}
##  that is, we model spectrum of each object i with sum of K linear components
##  Note: use SDSS pixel flux densities f_{ij} interpolated by cubic-spline interpolation
##  onto a common rest-frame wavelength grid, logarithmically spaced in wavelength
##
##  training set =  subset of data at rest-frame wavelengths
##
##  We use the training set to define a set of basis functions that minimize scalar \Chi_{\epsilon}**2
##  Training points f_{ij} with errors \sigma_{ij} by set of K components g_{kj} and coefficients a_{ik}
##  over all N spectra nd M pixels of a training set.
##
##  In order to minimize \Chi_{\epsilon}**2, do step-by-step procedure:
##  a-step: fix g_kj ---> estimate optimal a_ij
##  g-step: fix a_ik ---> estimate optimal g_kj
##
##  Start at an all-positive guess and then iterate!
##  Iteration proceeds to convergence
##
##  PROCEDURE
##  Use HMF to define a small number of components that is sufficient for modeling SDSS spectra.
##  Using these components, fit each observed spectrum at the redshift provided by SDSS.
##  Then repeat the fitting, but this time using one set of components at the SDSS redshift
##  and one set of components at values of redshift that lie on a nominal grid. If a second
##  object is present we expect the fit to be significantly improved when we use two sets of
##  components, one at the redshift of SDSS and one at that of the second object.
##
##  This second redshift is scanning a regular grid of values selected to be uniform in a logarithmic scale,
##  i.e. in the same way as the wavelengths in SDSS. In this way moving to the next value of the redshift grid is
##  equivalent to shifting the spectrum by one pixel.
##
##
##  INITIALIZATION
##  Use the mean spectrum and the top [K − 1] eigenvectors from the PCA in the orthogonal subspace
##  i.e. use PCA to produce a mean spectrum and then the [K − 1] top eigenvectors orthogonal to it
##  i.e. before performing PCA, we project the spectra into a subspace orthogonal to the mean spectrum,
##  by scaling them and subtracting the mean spectrum from each one of them
##
##  PCA results used as an initialization. The training data were first
##  projected into a hyperplane orthogonal to the mean spectrum that is passing from the
##  zero point. i.e.
##  Each spectrum was scaled appropriately and the mean spectrum was
##  subtracted from it. The flux in each spectral bin was divided with the RMS
##  of the error in that pixel for all the non-masked pixels in the training sample.
##
##
##
##

###
### e.g.
### http://arxiv.org/pdf/1106.1180.pdf
### BHB candidates in SDSS QSO spectra based on large velocity shifts between NLs and BLs
### NL = narrow line
### BL = broad line
### QSO spectra 0.1 < z < 1.5
###
### Assuming one BH is active, BLs associated to the accreting BH expected to show systematic velocity shifts
### w/r.t. NL, which trace the rest-frame of the galaxy
###
### GOAL: Search for two sets of emission lines NL and BL withsmall separation between them,
### caused by the Keplerian rotation of one component of the binary system
###
### Training set = subset of data at rest-frame wavelengths
###
### STEP 1: training set defines set of basis functions that minimize \chi squared
### STEP 2: use resulting components a, g to fit each observed spectrum at redshift z_SDSS, redshift provided by SDSS
### STEP 3: repeat the fitting using two sets of components at different redshifts, i.e.
###         one set z_SDSS and the other set z_2nd free to vary over broad range of redshifts s.t.
###         z corresponding to velocity differences up to 30 000 km/s
###         If second redshift system present, fit significantly improves when
###         second set of components added
###

##
## GOAL: Find double redshifts in SDSS LRG spectroscopy
##
## http://hoggideas.blogspot.com/2012/10/find-lrg-lrg-double-redshifts.html
##
## Why?
## DWH: "Rarely have double redshifts been found without emission lines
## (LRG spectra are almost purely stellar with no nebular lines),
## and because the LRGs sometimes host radio sources,
## you might even get a Hubble-constant-measuring golden lens"
##
## **Hubble-constant-measuring golden lens** It's worth a try!
##
## Recall LRG: http://arxiv.org/abs/astro-ph/0212087
##             http://iopscience.iop.org/article/10.1086/323717/pdf
## LRGs show little rotation, have smooth radial profiles,
## are massive and kinematically hot, have a narrow range of stellarmass-to-light
## ratios, show high metallicities, and reside
## preferentially in the Universe’s denser environments
##
## LRGs are old elliptical galaxies, well-relaxed systems with very uniformed properties,
## i.e. have very regular colors.
## Uniformity means clear relationship between matter and galaxies.
## LRGs are extremely luminous, so perfect for LSS to probe volume of Universe via spatial clustering.
##

from rpy2 import *
import rpy2.robjects as robjects

print("spectra.Rdata")
robjects.r['load']("spectra.Rdata")
spectraRdata = robjects.r['spectra']
print type(spectraRdata)
print len(spectraRdata)
spectra = np.asarray(spectraRdata)  # convert to numpy array
print spectra.shape

print("noise.Rdata")
robjects.r['load']("noise.RData")
noiseRdata = robjects.r['noise']
print type(noiseRdata)
print len(noiseRdata)
noise = np.asarray(noiseRdata)  # convert to numpy array
print noise.shape

print("mask.Rdata")
robjects.r['load']("mask.RData")
maskRdata = robjects.r['mask']
print type(maskRdata)
print len(maskRdata)
mask = np.asarray(maskRdata)  # convert to numpy array
print mask.shape

print("lamda_tel.Rdata")
print("wavelength coverage 3580.964 < λ < 9109.615 ˚A for Main Galaxies")
print("")
robjects.r['load']("lamda_tel.RData")  # wavelengths coverage of Massive Galaxies
lamda_telRdata = robjects.r['lamda_tel']
print type(lamda_telRdata)
print len(lamda_telRdata)
lamda_tel = np.asarray(lamda_telRdata)  # convert to numpy array
print lamda_tel.shape

# In order to run the R code from Tsalmantza/Hogg, simply run this
# It will run for hours
#
r_source = robjects.r['source'] # R function source()
r_source('functions.R')
#
r_source = robjects.r['source'] # R function source()
r_source('iterations20f.R')
#
# SAMPLE OUTPUT
#
#[1] "loading..."
#[1] "initializing..."
#[1] "iterating..."
#[1] 1
#[1]  1.052918e+07  9.388070e+06 -1.848939e-15  1.116262e-15 -2.799753e+00
#[6]  4.560790e-01
#[1] 2
#[1]  9.238428e+06  9.188643e+06 -1.759880e-15  1.770310e-15 -2.528090e+00
#[6]  1.320058e+00
#[1] 3
#[1]  9.169830e+06  9.163552e+06 -1.771468e-15  1.782187e-15 -1.321773e+00
#[6]  2.442542e+00
#[1] 4
#[1]  9.161093e+06  9.160487e+06 -1.812610e-15  1.772111e-15 -2.414206e+00
#[6]  1.322311e+00
# ....
#

