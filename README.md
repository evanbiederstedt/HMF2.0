# HMF2.0
HMF to find LRG LRG double z

Currently data files too large. Waiting for Git Large File Storage to kick in. 

HMF outperforms PCA by modeling *measured uncertainties* within the data, providing a more accurate set of derived basis functions. Particularly useful for missing and low S/N data.

Details: HMF finds basis set and coefficients that minimize chi-squared, takes into account individual spectral pixel uncertainty variances and missing data. 

Formalism: there are K components g_{jk} and coefficients a_{ik}, over all the N spectra and the M pixels of the training set
[N,K] coefficients a_{ik}
[K,M] parameters g_{kj} in basis spectra
[N,M] data
where K is the manually set complexity parameter. 

For succinct summary of HMF applied to search of massive BHBs in galxies and quasars, see Section 2 in the paper: http://arxiv.org/pdf/1106.1180.pdf
