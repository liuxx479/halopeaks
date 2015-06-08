#!python
# Jia Liu 2015/05/24
# This code reads in AHF maps and Gadget maps
# and return the peaks/ps for each set

#import WLanalysis
import numpy as np
from scipy import *
import sys, glob, os
#sys.modules["mpi4py"] = None
from emcee.utils import MPIPool
#from lenstools import Ensemble
from lenstools import ConvergenceMap 
from lenstools.defaults import load_fits_default_convergence

kmin = -0.08
kmax = 0.12
thresholds = linspace(kmin, kmax, 26)

peaksGen = lambda fn: ConvergenceMap.load(fn).peakCount(thresholds)

home = '/work/02977/jialiu/lenstools_home/'

storage = '/scratch/02977/jialiu/lenstools_storage/'

glob.glob()

peaks_fn_amiga = os.path.join(home,'Om0.300_Ol0.700/512b240/peaks_amiga.npy')

peaks_fn_gadget = os.path.join(home,'Om0.300_Ol0.700/512b240/peaks_gadget.npy')
