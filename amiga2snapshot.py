#!python
# Jia Liu 2015/05/19
# This code reads in (1) .AHF_particles and (2) gadget2 snapshots
# Then match the IDs, and keep only particles in gadget2
# Create new snapshots that contains only halo particles

import WLanalysis
import numpy as np
from scipy import *
import sys, glob, os
sys.modules["mpi4py"] = None
#from scipy.integrate import quad
#import scipy.optimize as op
#from scipy import interpolate
#import matplotlib.pyplot as plt
#from pylab import *
#import matplotlib.gridspec as gridspec
#from scipy import ndimage as snd
from lenstools.simulations import Gadget2Snapshot

home = '/work/02977/jialiu/lenstools_home/'
storage = '/scratch/02977/jialiu/lenstools_storage/'

#snapshots_dir = '/scratch/02977/jialiu/lenstools_storage/Om0.300_Ol0.700/512b240/ic1/snapshots/'
#amiga_dir = '/scratch/02977/jialiu/lenstools_storage/Om0.300_Ol0.700/512b240/ic1/amiga/'
#snap_amiga_dir = '/scratch/02977/jialiu/lenstools_storage/Om0.300_Ol0.700/512b240/ic1/snapshots_amiga/'
ID = 'Om0.300_Ol0.700|512b240|ic1'
######### reads in .AHF_particles ###########

def halo_particles(ID, snap_id):
	'''
	input: 
	ID = e.g. 'Om0.300_Ol0.700|512b240|ic1'
	snap_id = N, where N runs from 0 to total number of snapshots
	i = n, where n is one of the 16 split files for each snapshots
	'''
	cosmo_id,geometry_id, ic_id = ID.split("|")
	
	#### file names for gadget snap and AHF particles
	#isnap_fn = lambda i: os.path.join(storage, cosmo_id, geometry_id, ic_id, 'snapshots/snapshot_%03d.%i'%(snap_id, i))#####
	#ihalo_fn = lambda i: os.path.join(storage, cosmo_id, geometry_id, ic_id, 'amiga/snap%i.%04d.z0.000.AHF_particles'%(snap_id, i))
	halo_fn_arr = glob.glob(os.path.join(storage, cosmo_id, geometry_id, ic_id, 'amiga'))
	snap_fn_arr = glob.glob(os.path.join(storage, cosmo_id, geometry_id, ic_id, 'snapshots/snapshot_%03d.*'%(snap_id)))
	
	halo_fn_arr = glob.glob("*particles")
	snap_fn_arr = glob.glob('snapshot_060.*')
	
	txt_amiga = concatenate([genfromtxt(halo_fn, skiprows=2) for halo_fn in halo_fn_arr], axis = 0).T
	ID_amiga = txt_amiga[0][txt_amiga[1]==1]
	
	def ihalo_ID_position(snap_fn):
		print snap_fn
		snaps_gadget = Gadget2Snapshot.open(snap_fn)
		ID_gadget = snaps_gadget.getID() 
		idx = where(in1d (ID_gadget, ID_amiga, assume_unique=1) == True)[0]
		ID_HaloParticles = ID_gadget[idx]
		Positions_HaloParticles = snaps_gadget.getPositions()[idx] 
		return ID_HaloParticles, Positions_HaloParticles
	
	halo_ID_position = map(ihalo_ID_position, snap_fn_arr)
	halo_ID = concatenate([halo_ID_position[i][0] for i in range(len(halo_ID_position))])
	halo_position = concatenate([halo_ID_position[i][1] for i in range(len(halo_ID_position))], axis=0)
	return halo_ID, halo_position