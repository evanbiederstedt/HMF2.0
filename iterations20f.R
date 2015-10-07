## Copyright 2010 Paraskevi Tsalmantza & David W. Hogg.
## All rights reserved.

setwd("/Users/evanbiederstedt/desktop/Tsalmantza_Routines")

nonnegative<-FALSE
print("loading...")
load("spectra.Rdata")
load("noise.Rdata")
load("mask.Rdata")
load("lamda_tel.Rdata")
N<-nrow(spectra)      # nrow(spectra) = 5000
M<-length(lamda_tel)  # length(lamda_tel) = 4056

source("functions.R")

# length(mask) = 20280000
# str(mask) = [1:5000, 1:4056]
# nrow(mask) = 5000
# ncol(mask) = 4056
#
# Each of the N observe spectra i can be considered a column vector
# on a grid of M observer-frame wavelengths j
# HUGE MATRIX!
# f_{ij} i = N spectra, j = M wavelengths
# i.e. n_{ij} is the noise in pixel j of spectrum i
#
## Check which spectra have the whole spectrum at the beggining and at the end
all_blue<-which(mask[,1]==0)          # "which first column of mask equals 0" GIVES INDICES of MASK == 0 IN FIRST COLUMN
all_red<-which(mask[,ncol(mask)]==0)  # "which columns of mask equal 0" GIVES INDICES of MASK == 0 in LAST COLUMN
all_br<-union(all_blue,all_red)       # str(all_br) = int [1:24] 337 1127 1203 1500 1803 2790 3922 ...

#
# all_br =  263  909  973 1213 1451 2262 3188 3446 3740 4015 4048  312  579  704 1168
#                1403 1485 1550 2306 2423 3838 3906 3946 4009
#
set.seed(8)
# N is nrow(spectra) = 5000
# M is length(lamda_tel) = 4056
# N-M is 944
index<-sample(setdiff(c(1:N),all_br),(N-M)) # length(index) = 944, str(index) = int [1:944] 2333 1038 3994 3257 1607 3589 1453

#
# setdiff(c(1:N), all_br)) --- all all_br subtracted , "set difference" from union
#
# sample takes a sample of the specified size from the elements of x using either with or without replacement.
# sample(x, size, replace = FALSE, prob = NULL)
#
# sample size = 944
#
#
#
# ARRAY1  = 'SPECTRUM'                      / units of (10^-17 erg/s/cm^2/A
# ARRAY2  = 'CONTINUUM-SUBTRACTED SPECTRUM' / units of 10^-17 erg/s/cm^2/A
# ARRAY3  = 'ERROR   '                      / units of 10^-17 erg/s/cm^2/A
# ARRAY4  = 'MASK    '                      / mask bit array
#

spectra1<-spectra[-index,]   # str(spectra1) = num [1:4056, 1:4056] 1.37e-15 1.70e-16 1.06e-16 1.94e-17 ...
spectrat<-spectra[index,]    # str(spectrat) = num [1:944, 1:4056] 3.27e-17 7.19e-17 1.03e-16 5.47e-17 ...
spectra<-spectra1            # str(spectra) =  num [1:4056, 1:4056] 1.37e-15 1.70e-16 1.06e-16 1.94e-17 ...
noise1<-noise[-index,]       # str(noise1) = num [1:4056, 1:4056] 1e-12 1e-12 1e-12 1e-12 ...
noiset<-noise[index,]        # str(noiset) = num [1:944, 1:4056] 1e-12 1e-12 1e-12 1e-12 ...
noise<-noise1                # str(noise) = num [1:4056, 1:4056] 1e-12 1e-12 1e-12 1e-12 ...
noise<-abs(noise) # should be a no-op
mask1<-mask[-index,]         # str(mask1) = num [1:4056, 1:4056] 1 1 1 1 1 ...
maskt<-mask[index,]          # str(maskt) = num [1:944, 1:4056] 1 1 1 1 1 ...
mask<-mask1                  # str(mask) = num [1:4056, 1:4056] 1 1 1 1 1 ...

all_blue<-which(mask[,1]==0)
all_red<-which(mask[,ncol(mask)]==0)
all_br<-union(all_blue,all_red)

## Keep the the above spectra and the 200 first ones (because this is a test)
#M<-100
spectra200<-as.data.frame(spectra[union(all_br,1:1000),1:M]) # num [1:1018, 1:4056]
noise200<-as.data.frame(noise[union(all_br,1:1000),1:M])     # 'data.frame':	1018 obs. of  4056 variables:
N<-nrow(spectra200)                                          # nrow(spectra200) = 1018

if (nonnegative){
print("making the data non-negative...")
## make the data non-negative, as weakly as possible
negativedata<-(spectra200<0)
spectra200[negativedata]<-0
}

## get everything into matrix form for iterations
spectra200<-as.matrix(spectra200)
invvar200<-as.matrix(1.0/noise200^2)

Kminus1list<-c(1:15)              # 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
epsilonlist<-c(30.,10.,3.,1.)     # 30 10  3  1
nK<-length(Kminus1list)           # 15
nepsilon<-length(epsilonlist)     # 4
xit<-matrix(0.,nepsilon,nK)       # str(xit) = num [1:4, 1:15] 0 0 0 0 0 0 0 0 0 0 ...
Kindex<-0
for (Kminus1 in Kminus1list) {
Kindex<-Kindex+1
K<-Kminus1+1

epsilonindex<-0
for (epsilon in epsilonlist) {
epsilonindex<-epsilonindex+1

## initialize
print("initializing...")
#initialize with Kmeans
g<-kmeans(spectra200,K)$center
g<-g/normbase(g)
a<-outer(sqrt(rowMeans(spectra200^2)),rep(1./K,K))
if (nonnegative){
  for (aiter in 1:128){
    a<-astepnn(a,spectra200,invvar200,g)
  }
}
if (nonnegative){
if (min(spectra200)<0){ stop("spectra200 has negative elements!") }
if (min(a)<0){ stop("a has negative elements!") }
if (min(g)<0){ stop("g has negative elements!") }
epsilon<-0
onomadir<-as.character(paste("K",K,"M",M,"N",N,"_nn",sep=""))
}else{
onomadir<-as.character(paste("K",K,"M",M,"N",N,sep=""))
}
system(paste("mkdir -p",onomadir,sep=" "))
onomabase0<-as.character(paste(onomadir,"/g_0.R",sep=""))
save(g,file=onomabase0)

print("iterating...")
if(nonnegative){
  n_iter<-2048
}else{
  n_iter<-16
}
xi1wn<-vector(length=n_iter)
xk1wn<-vector(length=n_iter)
for (m in 1:n_iter){
print(m)

oldg<-g

if(nonnegative){
a<-astepnn(a,spectra200,invvar200,g)
}else{
a<-astep(spectra200,invvar200,g)
}
xi1wn[m]<-badness(a,g,spectra200,invvar200,0.0)

if(nonnegative){
g<-gstepnn(g,spectra200,invvar200,a)
}else{
g<-gstep(spectra200,invvar200,a,g,epsilon)
}
xk1wn[m]<-badness(a,g,spectra200,invvar200,epsilon)

## reorder and normalize
if(!nonnegative){
  foo<-reorder(a,g)
  a<-foo$a
  g<-foo$g
}
norm<-normbase(g)
g<-g/norm
a<-t(t(a)*norm)
print(c(xi1wn[m],xk1wn[m],range(a),range(g)))

## Save coefficients from the 1st fitting (a), the new base from the 2nd f (base), the sum of the x2 of all the spectra (xi, xi3 - 1st f), the sum of the the x2 for all the lamda (xk, xk3 - 2nd f) and the x2 for every spectrum (xi5 -1st f)
if (m == 2^round(log2(m))){
epssuffix<-as.character(paste("_",log10(epsilon),".R",sep=""))
suffix<-as.character(paste("_",m,epssuffix,sep=""))
onoma<-as.character(paste(onomadir,"/a",suffix,sep=""))
save(a,file=onoma)
onoma1<-as.character(paste(onomadir,"/g",suffix,sep=""))
save(g,file=onoma1)
}
}

## How well does this converged model do on new data?
spectrat<-as.matrix(spectrat[,1:M])
invvart<-as.matrix(1.0/noiset[,1:M]^2)
at<-astep(spectrat,invvart,g)
xit[epsilonindex,Kindex]<-badness(at,g,spectrat,invvart,0.)

onoma11<-as.character(paste(onomadir,"/xi1wn",epssuffix,sep=""))
save(xi1wn,file=onoma11)
onoma12<-as.character(paste(onomadir,"/xk1wn",epssuffix,sep=""))
save(xk1wn,file=onoma12)
}
}
if (nonnegative){
onoma13<-as.character(paste("xit_",(min(Kminus1list)+1),"_",(max(Kminus1list)+1),"_nn.R",sep=""))
  save(xit,file=onoma13)
}else{
onoma14<-as.character(paste("xit_",(min(Kminus1list)+1),"_",(max(Kminus1list)+1),".R",sep=""))
  save(xit,file=onoma14)
}

