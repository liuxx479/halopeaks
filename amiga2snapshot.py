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
from lenstools.simulations import Gadget2Snapshot
from astropy.units import Mpc,m,s
from emcee.utils import MPIPool

home = '/work/02977/jialiu/lenstools_home/'
storage = '/scratch/02977/jialiu/lenstools_storage/'
ID_arr = genfromtxt(os.path.join(home, 'realizations.txt'), dtype=str)
#ID = 'Om0.300_Ol0.700|512b240|ic1'

def halo_particles(IDsnap_id):
	'''
	input: 
	ID = e.g. 'Om0.300_Ol0.700|512b240|ic1'
	snap_id = N, where N runs from 0 to total number of snapshots
	i = n, where n is one of the 16 split files for each snapshots
	'''
	ID, snap_id = IDsnap_id
	cosmo_id,geometry_id, ic_id = ID.split("|")
	
	#### file names for gadget snap and AHF particles
	amiga_dir = os.path.join(storage, cosmo_id, geometry_id, ic_id, 'amiga')
	os.system('mkdir -p %s'%(amiga_dir))
	halo_fn_arr = glob.glob(amiga_dir+'/*particle*')
	snap_fn_arr = glob.glob(os.path.join(storage, cosmo_id, geometry_id, ic_id, 'snapshots/snapshot_%03d.*'%(snap_id)))
	
	##### test on laptop #####
	#halo_fn_arr = glob.glob("*particles")
	#snap_fn_arr = glob.glob('snapshot_060.*')
	
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
	
	halo_snap = Gadget2Snapshot()
	hg = Gadget2Snapshot.open(snap_fn_arr[0]).header #header_gadget
	halo_snap.setHeaderInfo(Om0=hg['Om0'], Ode0=hg['Ode0'], w0=hg['w0'], wa=hg['wa'], h=hg['h'], redshift=hg['redshift'], box_size=hg['box_size'])
	halo_snap.setPositions(array(halo_position)*Mpc)
	
	halo_snap.write(os.path.join(storage, cosmo_id, geometry_id, ic_id, 'snapshots_amiga/snapshot_%03d.*'%(snap_id)), files = len(snap_fn_arr))

	### test on laptop
	#halo_snap.write('snapshots_amiga/snapshot_%03d'%(snap_id), files = len(snap_fn_arr))

pool = MPIPool()
pool.map(halo_particles, [[ID, snap_id] for ID in ID_arr, snap_id in range(61)])

