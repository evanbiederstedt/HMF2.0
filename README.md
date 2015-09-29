# HMF2.0
HMF to find LRG LRG double z

Currently data files too large. Waiting for Git Large File Storage to kick in. 

=====

HMF outperforms PCA by modeling *measured uncertainties* within the data, providing a more accurate set of derived basis functions. Particularly useful for missing and low S/N data.

Details: HMF finds basis set and coefficients that minimize chi-squared, takes into account individual spectral pixel uncertainty variances and missing data. 

Formalism: there are K components g_{jk} and coefficients a_{ik}, over all the N spectra and the M pixels of the training set
[N,K] coefficients a_{ik}
[K,M] parameters g_{kj} in basis spectra
[N,M] data
where K is the manually set complexity parameter. 

For succinct summary of HMF applied to search of massive BHBs in galxies and quasars, see Section 2 in the paper: http://arxiv.org/pdf/1106.1180.pdf

=======

##### Steps:
1) Use a subset of data at rest-frame wavelengths as a training set to define the set of basis functions that minimize chi-squared

2) Use the resulting set of components to fit each observed spectrum at the redshift provided by SDSS, z_SDSS

3) Repeat the fitting using two sets of components at different redshifts, the first z_SDSS, the second free to vary over a broad range of z, corresponding to velocity differences up to 30 000 km s^−1

========
#####Why does PCA fall short? 

From Vaderplas, Connolly 2009; http://arxiv.org/pdf/0907.2238.pdf

"[PCA] cannot easily (nor compactly) express inherently non-linear relations within the data (such as dust
obscuration or the variation in spectral line widths). Spectra that have a broad range in
line-widths require a large number of eignespectra to capture their intrinsic variance which
often results in the continuum emission and emission lines being treated independently
(Gy¨ory et al. 2008)." 

