"""
Magnolya: de novo copy number variation detection

Usage: magnolya [OPTIONS] -i [FILE]

Input:
  -i    Tab delimited text file with read counts (required)
  -l    Tab delimited text file with contig locations on reference
  
Output:
  -o    CNV calls [stdout]
  -f    Model parameters
  
Parameters:
  -m    Number of Poisson components [Estimate using Baysian Information Criterion]
  -g    Gamma settings {haploid, diploid, None} [None]
  -s    Minimum contig size to use for EM [500]
  -e    Log-likelihood change EM stopping criterion [10e-6]
  -r    Maximum iterations of EM [200]
  -t    CNV threshold: only contigs with CNV probability higher than this value are reported [report all]
  -S    Report Size threshold: don't report contigs smaller than this threshold [500]
  -M    Multiplier for zero copy number threshold. CN = 0 if count/length < mu - M*stdev [5]

Date:        July 25, 2012
Contact:     Jurgen F. Nijkamp, Dick de Ridder
Publication: De novo detection of copy number variation by co-assembly (2012) Bioinformatics

"""

import getopt
from numpy import array, mean,median, std, sum, float64, exp, int32, log10, argmax, argmin,zeros, sqrt
import scipy.stats as stats
from math import isnan
from sys import stdout, stderr
import sys
from copy import deepcopy
import numpy as np
import scipy as sp
import scipy.stats as stats
import pdb
import math
import time
import random
from sys import stderr
import itertools
import pdb

def poisson_geometric(x,t,A,mix,mu,tries,N,alpha):

    x_delay  = x - np.int32(mu[-1]*t)
    minxt    = mu[-1]

    # P(data|model i) and P(data)
    pxi = np.array([np.zeros(N)]*int(A)); px = np.zeros(N);
    for i in range(int(A)):
        if i < int(A)-1:
            # Fill in poisson
            pxi[i,:]          = stats.poisson.pmf(x,mu[i]*t) # Add machine precision
            px                = px + mix[i] * pxi[i,:]

        else:
            # Fill in geometric
            pxi[i,:] = stats.distributions.geom.pmf(x_delay,alpha/np.float64(t))
            px       = px + mix[i] * pxi[i,:]
    return px, pxi, np.sum([x_delay > 0])


def poisson(x,t,A,mix,mu,tries,N):

    # P(data|model i) and P(data)
    pxi = np.array([np.zeros(N)]*A); px = np.zeros(N);
    for i in range(int(A)):
      # Fill in poisson
      pxi[i,:]          = stats.poisson.pmf(x,mu[i]*t)
      px                = px + mix[i] * pxi[i,:]
    return px,pxi

def dirichlet(gamma,n,I):
    """Initialize Dirichlet parameter
       INPUT:
       gamma:    'haploid', 'diploid', 'uniform', or an list of gamma's of length n
       n:        Number of models
       I:        Hyperparameter that controls the prior's impact during model selection (I>0)

       OUTPUT:
       gamma array
    """
    if n<2:
        raise Exception, "Poisson mixture needs at leat two Poisson models"
    if I <= 0:
        raise Exception, "Hyperparameter I should be larger then 0"
    if gamma == 'haploid':
        res = (np.ones(n))
        res[0] += I
    elif gamma == 'diploid':
        res = (np.ones(n))
        res[1] += I
    elif gamma == 'allo-triploid':
        res = (np.ones(n))
        res[0] += np.float64(I)/2
        res[1] += np.float64(I)/2
    elif gamma == 'di-triploid':
        res = (np.ones(n))
        res[1] += np.float64(I)/2
        res[2] += np.float64(I)/2
    elif (type(gamma) == type(np.array([])) or type(gamma) == type([])) and len(gamma)==n:
        res = np.array(gamma)
    elif gamma is None or gamma == 'uniform':
        res = (np.ones(n)+np.float64(I)/n)
    else:
        raise Exception, "Provide a valid gamma, see the docstring of em.Dirichlet"
    return res




def emPoissonBic(count,interval,k=None,epsilon=2e-13,\
                 max_iters=100,gamma='haploid',figname='poisson.pdf', addUniform=1,\
                 initPercList=None,bw=True,plotContigSize=1000,maxval=None,ymax=None):

    print "****Starting model selection using BIC****"
    print "Number of contigs used to train EM: " , len(count)
    print "Median contig length: " , np.median(interval)

    plot_mixture=0

    if k is None:
        k=[3,4,6,8,10,12,15,20]
    elif type(k) is not type([]):
        raise Exception, "Please provide a range of cluster numbers in a List"

    if len(count) == 0:
      print >>stderr, "Warning: Zero reads detected"
      return [0]*len(k),[0]*len(k),[0]*len(k),[0]*len(k),k

    likelihoods = []; mus = []; mixs = []; alphas = [];
    for ki in k:
        likelihood, mu, mix, alpha = emPoissonInit(count,interval,ki,False,epsilon,\
                                            max_iters,gamma,figname, addUniform,initPercList,\
                                            plotContigSize,maxval)
        likelihoods.append(likelihood)
        mus.append(mu)
        mixs.append(mix)
        alphas.append(alpha)

    # Calculate Bayesian Information Criterion
    n    = len(count)
    BICS = []
    for i in range(len(likelihoods)):
        BICS.append(-2 * likelihoods[i] + (k[i]+1+2*addUniform)*np.log(n)) # pi's + lambda + geom (alpha and pi)
        print 'k = %d; BIC = %2.2f; L = %2.2f; lambda = %2.2e;' % \
            (k[i], BICS[i], likelihoods[i], mus[i][0])

    best = np.argmin(BICS)

    return likelihoods, mus, mixs, BICS, k, alphas


def emPoissonInit(count,interval,k=10,plot_mixture=0,epsilon=2e-13,\
              max_iters=100,gamma='haploid',figname='poisson.pdf', addUniform=1,initPercList=None,\
              plotContigSize=1000,maxval=None):
    """

    INITPERCLIST    List with percentiles to initialize the EM
    """


    if initPercList is None and (gamma is not 'haploid' and gamma is not 'diploid'):
      initPercList = [1, 5, 12.5, 25, 50, 70]
      #initPercList = [1, 50, 70]
      print "Performing initializations at 1, 5, 12.5, 25 and 50th percentile"
    elif initPercList is None and (gamma == 'haploid' or gamma == 'diploid'):
      initPercList = [50]
      print "Performing initialization at 50th percentile"
    elif type(initPercList) is not type([]):
      raise Exception, "Invalid initPercList. Provide a valid list of initialisation percentages or \'None\' for default initialisation"

    likelihoods = []; mus = []; mixs = [];alphas = [];
    for perc in initPercList:
        lambda_init = stats.scoreatpercentile(count/np.float64(interval),perc)
        likelihood, mu, mix, alpha = emPoisson(count,interval,k,False,epsilon,max_iters,\
                                        gamma,figname,addUniform,lambda_init,plotContigSize,maxval)
        likelihoods.append(likelihood)
        mus.append(mu)
        mixs.append(mix)
        alphas.append(alpha)
        print "Lambda init: %2.2e for %2.2f-th percentile -> L = %2.2f  " % (lambda_init, perc,likelihood)
    best = np.argmax(likelihoods)

    return likelihoods[best], mus[best], mixs[best], alphas[best]


def emPoisson(count,interval,k=10,plot_mixture=0,epsilon=2e-13,\
              max_iters=100,gamma='haploid',figname='poisson.pdf', addUniform=1, lambda_init=None,\
              plotContigSize=1000,maxval=None):
    """
    [L, MU, S2, MIX] = EM (COUNT, INTERVAL, K, STEPK, FIX_S2, REG, PLOT, EPS, MAXITER)

    Uses the EM algorithm to fit a mixture of K 1D Poissons, where the means
    and the variances of the models are spaced as (1 + k/STEPK), k = 0..K-1.
    Optional parameters are:
    FIX_S2	set the variance to the mean (1, dflt.) or estimate (0)
    REG     regularisation constant (dflt. 0)
    PLOT    plot progress (dflt. 0)
    EPS     epsilon, the change in likelihood to stop training (dflt. 1e-5)
    MAXITER maximum number of iterations (dlft. 100)
    Returns the likelihood L, the means MU(), the variances S2() and priors MIX.
    """

    start = time.clock()

    if addUniform not in [0,1]: raise Exception, "addUniform should be 0 or 1"
    k    = k+addUniform  # One additional k for the Uniform distribution
    x    = np.array(count)
    t    = np.array(interval)
    N    = np.float64(x.size)
    A    = np.float64(k)
    mu   = np.array([0.0]*k)
    mix  = np.array([0.0]*k)
    sumR = np.array([0.0]*k)
    alpha    = 0.0001
    mix_new  = np.array([0.0]*k)
    mu_new   = np.array([0.0]*k)
    min_data = max(min(x/np.float64(t)),epsilon); max_data = max(x/np.float64(t)) # Normalized for interval
    tries = 0; retry = 1;

    while ((retry == 1) and (tries <= 5)):

        # Initialisation
        if gamma == 'haploid':
            gamma = dirichlet(gamma,k,N)
            if lambda_init is None:
                lambda_init = np.median(x/np.float64(t))
                print ("Haploid initialization, lamda: " + str(lambda_init))
            for i in range(int(A)):
                mu[i]  = lambda_init * (i+1)
                mix[i] = 0.4/(k-1)
            mix[0] = 0.6
        elif gamma == 'diploid':
            gamma = dirichlet(gamma,k,N)
            if lambda_init is None:
                lambda_init = np.median(x/np.float64(t)) / 2
                print ("Diploid initialization, lamda: " + str(lambda_init))
            for i in range(int(A)):
                mu[i]  = lambda_init * (i+1)
                mix[i] = 0.4/(k-1)
            mix[0] = 0.6
        elif gamma == 'allo-triploid':
            gamma = dirichlet(gamma,k,N)
            if lambda_init is None:
                lambda_init = np.median(x/np.float64(t)) / 3
            for i in range(int(A)):
                mu[i]  = lambda_init * (i+1)
                mix[i] = 0.2/(k-2)
            mix[1] = 0.4
            mix[2] = 0.4
            print ("haplo-diploid initialization, lamda: " + str(lambda_init))
        elif gamma == 'di-triploid':
            gamma = dirichlet(gamma,k,N)
            if lambda_init is None:
                lambda_init = np.median(x/np.float64(t)) / 3
                print ("Di-triploid initialization, lamda: " + str(lambda_init))
            for i in range(int(A)):
                mu[i]  = lambda_init * (i+1)
                mix[i] = 0.2/(k-2)
            mix[1] = 0.4
            mix[2] = 0.4
        elif gamma == 'uniform':
            gamma = dirichlet(gamma,k,N)
            if lambda_init is None:
                lambda_init = np.median(x/np.float64(t))
                print ("Default initialization, lamda: " + str(lambda_init))
            for i in range(int(A)):
                mu[i]  = lambda_init * (i+1)
                mix[i] = 1./(k)
            #mix[-1] = 0.001
        elif gamma is None:
            if lambda_init is None:
                lambda_init = np.median(x/np.float64(t)) / 2
            print ("Default gamma initialization, lamda: " + str(lambda_init))
            for i in range(int(A)):
                mu[i]  = lambda_init * (i+1)
                mix[i] = 1./(k)
        else:
            gamma = dirichlet(gamma,k,N)
            if lambda_init is None:
                lambda_init = np.median(x/np.float64(t))
                print ("Default gamma initialization, lamda: " + str(lambda_init))
            for i in range(int(A)):
                mu[i]  = lambda_init * (i+1)
                mix[i] = 0.4/(k-1)
            mix[0] = 0.6

        R     = np.random.rand(N,A)      # N points     x A models
        done = 0; retry = 0;
        iter = 0; prev_likelihood = 1.0e20;

        while ((not done) and (not retry)):
            iter = iter + 1;
            done = 1;

            if addUniform:
                px,pxi,outliers = poisson_geometric(x,t,A,mix,mu,tries,N,alpha)
            else:
                px,pxi = poisson(x,t,A,mix,mu,tries,N)


            if (not retry):

                for i in range(int(A)):
                    for contig in range(len(pxi[i,:])):
                        if px[contig] > 0:
                            # Aandeel van model i in px
                            R[contig,i] = np.transpose((pxi[i,contig] * mix[i]) / px[contig])
                        else:
                            R[contig,i] = 0
                    sumR[i]     = sum(R[:,i])
                    if gamma is not None:
                        mix_new[i]  = ((1/N)*sumR[i] + (1/N)*(gamma[i]-1)) / (1 + (sum(gamma)-k)/N);
                    else:
                        mix_new[i] = sumR[i] / N;  # Gemiddelde aandeel alle datapunten

                likelihood = sum(np.log(px+epsilon))

                for i in range(int(A)):
                    likelihood += sum(R[:,i]*np.log(mix[i]+epsilon))
                    if gamma is not None:
                      likelihood += (gamma[i]-1)*np.log(mix[i]+epsilon)

                denomgeom                  = A*t*R[:,-1]*np.log(1-alpha/np.float64(t))
                denomgeom[x-A*t*mu[0] < 0] = 0
                denomgeom_sum              = sum(denomgeom)

                denompois                     = 0
                for i in range(0,int(A)-addUniform):
                    denompois += sum(R[:,i]*(i+1)*t)
                denom = denomgeom_sum + denompois

                numer = 0
                for i in range(0,int(A)-addUniform):
                    numer += sum(R[:,i]*x)

                mu_new[0] = numer/denom
                for i in range(1,int(A)):
                    mu_new[i] = mu_new[0]*(i+1)

                alpha = sumR[-1] / sum(R[:,-1]*(x/np.float64(t)-mu[-1]+1))  # alpha geometric

                if math.isnan(alpha):
                  alpha = 10e-6

                if (iter%100 == 0):
                    print '[%3d] L: %2.2f (change: %2.2e); sum (P(j|x)) = ; alpha = %2.2e; lambda = %2.2e' % \
                                 (iter, likelihood, abs((likelihood - prev_likelihood)/likelihood),alpha,mu_new[0]),
                    for i in range(int(A)):
                        print '%2.2f ' % sumR[i],
                    print '; P(j) = ',
                    for i in range(int(A)-1):
                        print '%2.2f ' % mix[i],
                    print '\n',
                    print 'Outliers:',outliers,
                    print '\n'


                done = (abs ((likelihood - prev_likelihood)/likelihood) < 1e-5); # Controle voor verbetering
                if iter >= max_iters:
                    print "Maximum number of iterations reached"
                    done = True

                if done:
                  print "Number of iterations needed to converge: ", iter

                prev_likelihood = likelihood;

                # Update

                for i in range(int(A)):
                    mix[i]    = mix_new[i];
                    mu[i]     = mu_new[i];

    elapsed = (time.clock() - start)
    print "EM took ", elapsed, " to finish for " , len(count), " contigs"

    return likelihood, mu[:k-addUniform], mix[:k-addUniform], alpha


class Contiglocs(dict):
  
  def __init__(self, tabfile):
    """Read contig locations from tab delimited file. Columns in file should be:
       chromosome <tab> start <tab> end <tab> contigID
       (show-tiling out/.delta | tiling2tab.py)
    """
        
    H  = open(tabfile)
    for line in H:
        chr,start,end,id,orient,clen = line.split("\t")
        id = id.rstrip()
        self[id] = {"chr":chr,"start":int32(start),"end":int32(end),"orient":orient,"clen":int32(clen)}
    H.close()


class Cnv():
    """Calculate CNV based on mixture model"""
    def __init__(self, data, nrsamples, k, epsilon, max_iters, gammas, minSize,minCP,maxCP):
        
        if type(k) == type(2):
          k = [k]
        self.mixtures   = []
        self.posteriors = []
        first           = True
        figname = 'poisson' + str(0) + '.pdf'
        [counts,intervals] = self.getCountsAndlenghts(data, 0)
        [counts,intervals] = self._filterContigs(counts,intervals,minSize,minCP,maxCP)
        print "Fitting mixture model sample: " + str(0)
        self.mixtures.append(PoissonMixture(counts, intervals,k, epsilon,max_iters,gammas[0]))

        # Calculate optimal k base on a combined BIC
        BICall = zeros(len(self.mixtures[0].BICs))
        if self.mixtures[0].BICs is not None:
           BICall += self.mixtures[0].BICs
        best = argmin(BICall)
        print "Combined BIC optimal is: ", self.mixtures[0].ks[best]
        
        # Set new k
        self.mixtures[0].setModelNumber(best)
        
        print "Calculating posteriors"
        self.posteriors.append(self.posteriorAllContigs(self.mixtures[0],data,0))
          
    def printParameters(self,prefix):
      if prefix == "":
        modelfile = "model.txt"
      else:
        modelfile = prefix + ".model.txt"
      h = open(modelfile,"w")
      for mixture in self.mixtures:
        str = '%d\t%.6f\t%.6f' % (len(mixture.mu), mixture.mu[0], mixture.alpha)
        for i in range(len(mixture.mix)):
          str += '\t%.6f' % (mixture.mix[i])
        print >>h, str
      h.close()
      return
      

    def getCountsAndlenghts(self,data,sample):
      counts     = []; lengths    = []
      contigs = data.keys()
      for contig in contigs:
        counts.append(data[contig]['counts'][sample])
        lengths.append(data[contig]['clen'])
      return counts, lengths
    
    def cnv(self, sample1, sample2):
      pCnv = {}
      for node in self.posteriors[sample1]:
          pCnv[node] = self.probabilityCNV(self.posteriors[sample1][node],self.posteriors[sample2][node])
      return pCnv
  
    def cn(self,sample):
      res = {}
      for node in self.posteriors[sample]:
          res[node] = argmax(self.posteriors[sample][node]) + 1
      return res
        
    def isZeroCN(self, contigid, data, sample, multiplier):
      """Set CN to zero if it is lower than 'multiplier x lambda'"""
      mu     = self.mixtures[sample].mu[0]
      counts = data[contigid]['counts'][sample]
      if (counts/float64(data[contigid]['clen']) < (mu - multiplier*sqrt(mu))):
        return True
      else:
        return False
        
    #def printCN(self,am,sample,ref=None,file=stdout):
    #  if file != stdout:
    #      file = open(file,"w")
    #  cns = self.cn(sample)
    #  for node in am.ass.contigGraph.contigStats:
    #      count  = am.ass.contigGraph.contigStats[node]['counts'][0]
    #      length = am.ass.contigGraph.contigStats[node]['length']
    #      cn     = cns[node]
    #      if ref and ref.contigLocation.has_key(node):
    #          chr    = ref.contigLocation[node]['chr']
    #          start  = ref.contigLocation[node]['feat'].location.nofuzzy_start+1
    #          end    = ref.contigLocation[node]['feat'].location.nofuzzy_end
    #      else:
    #          chr    = ''
    #          start  = 0
    #          end    = 0
    #      # contig ID, length, count, CN, CHR, start, stop
    #      print >>file,'%s\t%d\t%d\t%d\t%s\t%d\t%d' % (node,length,count,cn+1,chr,start,end)
    #  if file != stdout:
    #      file.close()    
    
    def posteriorAllContigs(self,mixture, data, sample):
      """Calculate posteriar probabilities for each contig"""
      prob = {}
      for node in data.keys():
          # TODO fix this in absence of am
          #if am.ass.contigGraph.contigStats[node]['sampleMems'][sample]:
          if data[node]['counts'][sample] >= 1:
            prob[node] = mixture._posterior(data[node]['counts'][sample],data[node]['clen'])
          else:
            prob[node] = [0]*len(mixture.mix)
      return prob
    
    def probabilityCNVAllContigs(self, s1, s2):
      res   = {}
      post1 = self.posteriors[s1]
      post2 = self.posteriors[s2]
      for node in post1.keys():
          res[node] = self.probabilityCNV(post1[node], post2[node])
      return res
        
    def probabilityCNV(self, post1, post2):
      """Probability a contig has a different copy number in the two samples
      Input: Two lists with posterior probabilities for all copy numbers (k)"""
      p_equal = 0
      for k in range(len(post1)):
          p_equal += post1[k] * post2[k]
      return 1-p_equal
        
    def _filterContigs(self, counts, intervals, minSize, minCP, maxCP):
      """Filter contigs with too little reads or too small
      minSize: minimum contig size
      minCP  : minimum number of reads per 100 nucleotides [1]
      maxCP  : maximum number of reads per 100 nucleotides [mean(count/100)*5]"""
      counts    = array(counts) 
      intervals = array(intervals)
      countsCP  = (counts / float64(intervals)) * 100
      if maxCP == "Default":
          maxCP = mean(countsCP) * 5
          
      print "Contigs < ", minSize, " are not used to fit the model"
      
      ind  =  (intervals >= minSize) & (countsCP  >= minCP) & (countsCP <= maxCP)                 
      return counts[ind], intervals[ind]


 
class PoissonMixture():
    """Poisson mixture model to extimate copy number for each contig"""
    def __init__(self,counts,intervals,k=None, epsilon=1e-15,max_iters=100,\
                      gamma=None):
                      
        print "Gamma = " , gamma
        
        self.counts    = counts
        self.intervals = intervals
        [self.L,self.mu,self.mix] = \
                self._emPoisson(counts, intervals,k, epsilon,max_iters,gamma)
        print "L: ",  self.L
        print "lambda: ", self.mu
        print "p(k): ", self.mix
    
    def setModelNumber(self, KsIndex):
        self.L = self.likelihoods[KsIndex]; self.mu = self.mus[KsIndex]; 
        self.mix = self.mixs[KsIndex]; self.alpha = self.alphas[KsIndex];
        print "Changed model paramters to:"
        print "L: ",  self.L
        print "lambda: ", self.mu
        print "p(k): ", self.mix
        
    def _priorDepth(self,x,t):
        px = float64(0)
        for k in range(len(self.mix)):
            px += self.mix[k]*self._pmf(x, self.mu[k],t)
        return px
        
    def _posterior(self,x,t):
        """Posterior probability:
        p(k|x) = p(k). (px|k).p(x))"""
        posterior = []
        px = self._priorDepth(x,t)
        for k in range(len(self.mix)):
            z = self.mix[k] * self._pmf(x, self.mu[k],t) / px
            if(isnan(z)):
                z = 10e-6;
            posterior.append(z)
        return posterior
        
    def _pmf(self,x,mu,t):
        """evaluate gaussian at x"""
        return stats.poisson.pmf(x,mu*t)
    
    def _emPoisson(self,counts, intervals,k, epsilon,max_iters,gamma):
        
        # Remove contigs with lenght zero
        counts    = array([counts[i]    for i in range(len(intervals)) if intervals[i]!=0])
        intervals = array([intervals[i] for i in range(len(intervals)) if intervals[i]!=0])
        
        # filter repeats with too high copy number
        x         = counts/float64(intervals)
        
        print "k: ", k                      
        [self.likelihoods, self.mus, self.mixs, self.BICs, self.ks, self.alphas] = \
              emPoissonBic(counts, intervals,k,epsilon,max_iters,gamma,\
                              addUniform=1,initPercList=None)
                              
                              
        best = argmin(self.BICs)
        L = self.likelihoods[best]; mu = self.mus[best]; mix = self.mixs[best]
            
        return L,mu,mix 

class Mutation():
    def __init__(self,type="mutation", location=None):
        self.type     = type
        self.location = location

    def setLocation(self,chr,start,end):
        """provide Python range for start to end, i.e. a deletion of
        length 5 at 2nd position of a sequence has coordinated start=1, end=6"""
        self.location = {"chr":chr,"loc":FeatureLocation(start,end)}


class CNVcall(Mutation):
    def __init__(self, cid, probCNV, posteriors,type="cnv"):
        Mutation.__init__(self, type=type)
        self.cid        = cid
        self.probCNV    = probCNV
        self.posteriors = posteriors
    

class Run():
    """Fit mixture models and find copy number variation"""
    
    def __init__(self,infile,gammas,minSize,locs,k,epsilon,max_iters,minCP,maxCP,cnvThrs,minReportSize,prefix,multiplier):
      
      data      = {}
      nrsamples = 0
      init      = False
      
      if locs is not None:
        locs = Contiglocs(locs)
      
      # Read counts from file
      h = open(infile)
      for line in h.readlines():
        vals          = line.split("\t")
        data[vals[0]] = {"clen":int32(vals[1]), "counts":int32(vals[2:])}
        if not init:
          nrsamples = 1
      
      if len(gammas) == 1:
        gammas *= nrsamples  # Each sample gets the sample gamma if 1 provided
      else:
        if len(gammas) != nrsamples:
          raise Exception, "Number of provided gamma's should be equal to number of samples"
      
      if nrsamples > 1:
        # Train model on samples to find CNVs
        self.cnv,self.probCNV,self.cnvC = self.getCnvContigs(data, k, epsilon, max_iters, gammas, minSize,minCP,maxCP,nrsamples,cnvThrs)
        # Print CNVs to file
        self.printCnv(data,locs,minReportSize,prefix,multiplier)

      elif nrsamples == 1:
        # Train model on single sample to find CNs
        self.cnv = self.getCNContigs(data, k, epsilon, max_iters, gammas, minSize,minCP,maxCP,nrsamples,cnvThrs)
        self.printCN(data,locs,minReportSize,prefix,multiplier)
                    
      # Print model parameters to file (for plotting etc)
      self.cnv.printParameters(prefix)

    def getCNContigs(self, data, k, epsilon, max_iters, gammas, minSize,minCP,maxCP,nrsamples,cnvThrs):
        res  = []
                
        cnv = Cnv(data, nrsamples, k, epsilon, max_iters, gammas, minSize, minCP, maxCP)
        return cnv
        
    def getCnvContigs(self, data, k, epsilon, max_iters, gammas, minSize,minCP,maxCP,nrsamples,cnvThrs):
        res  = []
        
        # TODO Change this to allow for multiple samples?
        sid1 = 0
        sid2 = 1
        
        cnv = Cnv(data, nrsamples, k, epsilon, max_iters, gammas, minSize, minCP, maxCP)
        probCNV = cnv.probabilityCNVAllContigs(sid1,sid2)
        for k in probCNV:
            if probCNV[k] > cnvThrs:
                posteriors = []
                posteriors.append(cnv.posteriors[sid1][k])
                posteriors.append(cnv.posteriors[sid2][k])
                res.append(CNVcall(k,probCNV[k],posteriors,type="cnv"))
        return cnv,probCNV,res
        

    def printCN(self,data,locs,minReportSize,prefix,multiplier):
        
        if prefix == "":
          file = "copynumbers.txt"
        else:
          file = prefix + ".copynumbers.txt"
        file = open(file,"w")
        
        cns = self.cnv.cn(0)
        
        for cid,cn in cns.iteritems():
          
          if data[cid]["clen"] < minReportSize:
            continue
          
          # Set copy number to zero for very low read counts
          if self.cnv.isZeroCN(cid, data, 0, multiplier): cn = 0
          
          if locs is not None and locs.has_key(cid):
            chr   = locs[cid]['chr']
            start = locs[cid]["start"]
            end   = locs[cid]["end"]
          else:
            chr   = ""
            start = -1
            end   = -1
          
          print >>file,'%s\t%d' % (cid, cn)
        if file != stdout:
            file.close()

    def printCnv(self,data,locs,minReportSize,prefix,multiplier):
        
        if prefix == "":
          file = "cnvcalls.txt"
        else:
          file = prefix + ".cnvcalls.txt"
        file = open(file,"w")
        
        for contig in sorted(self.cnvC, key=lambda cnvcall:(1-cnvcall.probCNV)):
          
          if data[contig.cid]["clen"] < minReportSize:
            continue
          
          cn1 = array(self.cnv.posteriors[0][contig.cid]).argmax() + 1
          cn2 = array(self.cnv.posteriors[1][contig.cid]).argmax() + 1
          
          # Set copy number to zero for very low read counts
          if self.cnv.isZeroCN(contig.cid, data, 0, multiplier): cn1 = 0
          if self.cnv.isZeroCN(contig.cid, data, 1, multiplier): cn2 = 0
          
          if locs is not None and locs.has_key(contig.cid):
            chr   = locs[contig.cid]['chr']
            start = locs[contig.cid]["start"]
            end   = locs[contig.cid]["end"]
          else:
            chr   = ""
            start = -1
            end   = -1
          
          print >>file,'%s\t%.6f\t%s\t%d\t%d\t%d\t%d\t%d' % (contig.cid,\
                                                  contig.probCNV,\
                                                  chr,\
                                                  start,\
                                                  end,\
                                                  data[contig.cid]["clen"],\
                                                  cn1,\
                                                  cn2)
        if file != stdout:
            file.close()
    
        
def main():
    
    # Set default
    infile    = None   # Input file with counts
    #outfile   = stdout # Output file name
    prefix    = ""
    gammas    = [None] # Gamma settings (optional)
    minSize   = 500    # Minimum contig size (optional)
    locs      = None   # Contig locations (optional)
    m         = None   # Number of models (optional)
    epsilon   = 10e-6
    max_iters = 500
    minCP     = 1
    maxCP     = 10e5
    cnvThrs   = -1
    minReportSize = 500
    multiplier    = 5

    
    # parse command line options
    try:
        opts, args = getopt.getopt(sys.argv[1:], "g:s:i:l:e:r:m:t:p:S:hM:", ["help"])
    except getopt.error, msg:
        print msg
        print "for help use --help"
        sys.exit(2)
    # process options
    for o, a in opts:
        if o in ("-h", "--help"):
            print __doc__
            sys.exit(0)
        if o == "-g":
            if a == "None":
                gamma = [None]
            else:
                gammas = list(a.split(","))
                #gamma = a
        if o == "-s":
            minSize = int(a)
        if o =="-i":
            infile = a
        if o =="-l":
            locs = a
        #if o =="-o":
        #    outfile = a
        if o =="-e":
            epsilon = a
        if o =="-r":
            max_iters = int32(a)
        if o =="-m":
            m = [int32(a)]
        if o =="-t":
            cnvThrs = float64(a)
        if o =="-S":
            minReportSize = float64(a)
        if o =="-p":
            prefix = a
        if o == "-M":
            multiplier = float64(a)
            
    if infile is None:
        print "Specifiy input file with counts"
        print "for help use --help"
        sys.exit(2)
    
    
    session = Run(infile,gammas,minSize,locs,m,epsilon,max_iters,minCP,maxCP,cnvThrs,minReportSize,prefix,multiplier)
    
    
if __name__ == "__main__":
    main()
        
