#
# Roughly using http://gastonsanchez.com/plsdepot_nipals.pdf and
# Sergey Koposov (U. Cambridge) https://github.com/segasai/pyhmf/blob/master/nipals.py
#
# data = input data
# Kcomps = number of components to be calculated, K dim
# conv  = convergence value
#
# Matrix X: the rows of X are the observations, the columns of X are the variables
#
# REMEMBER WHEN SETTING Kcomps, K = number of components permitted
# If K is too small, the model doesn’t have enough freedom to fit the training set well,
# and if K is too large, the model over-fits the training set and does worse on the test set.
#
#
# Recall PCA decomposition of matrix X such that X = \Sum^{q}_{h=1} t_h (p_h).T
# t denotes principal components: t_1, t_2, ..., t_h, rows
# p denotes principal axes:       p_1, p_2, ..., p_h, columns
#
#
import numpy as np
#
# Directions: input data; manually set Kcomps, conv, max_it
#
def nipals_pca(data, Kcomps = None, conv=1e-8, max_it=100000):
    X = np.matrix(data) # (N, M) = (rows, columns) = (N observed spectra, M pixels)
    Mpix = X.shape[1]   # M pixels
    # X.shape[1] is the number of columns, X.shape = (rows, columns)
    
    # Initialize eigenvec to zero matrix, dim [Kcomps, Mpix]
    eigenvec = np.zeros((Kcomps, Mpix)) # creates zero matrix sized [Kcomps, Mpix]
    for i in range(Kcomps): # range from 1 to Kcomps
        iterat = 0 # iterat = iteration counter, begin at iter = 0
        # t = principal components
        # STEP 1: set t_h to first column of X
        t = X[0, :] # X[0,:] lists all row entries, dimension [1, MPix]
        diff = conv + 1 # initialize the difference
        
        # STEP 2: repeat the following until convergence of p_h
        while diff > conv:
            iterat += 1 # increate iteration counter
            # Project X onto t to find the corresponding loading p
            # and normalize loading vector p to 1
            p = (X * t.T) / (t * t.T) # dim [N, 1]
            p = p/(np.sqrt(p.T * p))
            
            # project X onto p to find corresponding score vector t_new
            t_new = p.T * X
            # difference between new and old score vector
            tdiff = t_new - t
            diff = (tdiff * tdiff.T)
            t = t_new
            if iterat > max_it:
                msg = ('PC#%d: no convergence after'
                       ' %d iterations.'% (i, max_it))
                raise Exception(msg)
    
        eigenvec[i, :] = t # store ith eigenvector in result matrix
        # remove the estimated principal component from X
        D = np.matrix(np.outer(p, t), copy=False)
    X = X - D
return eigenvec


# Also include funtion to calculate PCA via numpy SVD
def svd_pca(data):
    # Factors the matrix M as U*np.diag(S)*V,
    # where U, V are unitary and S is a 1-D array of M‘s singular values
    U, S, eigenvec = np.linalg.svd(data)
    eigenval = S ** 2. / (data.shape[0] - 1.)
    return eigenvec, eigenval


