# SQUAD-CTHYB 
# CT-HYB QMC for a quantum dot connected to two superconducting and one normal leads
# uses TRIQS cthyb as the solver, tested with TRIQS 1.4
# Vladislav Pokorny; 2016-2017; pokornyv@fzu.cz

import scipy as sp
from os import listdir
from ConfigParser import SafeConfigParser
from libqmc14 import OmegaN

config = SafeConfigParser()

if 'cthyb.in' not in listdir('.'): 
	print('Main input file cthyb.in is missing. Exit.')
	exit(1)

config.read('cthyb.in')

P = {}

## defaults values ########################################
## floats in params
P['beta']     = 40.0 # inverse temperature
P['U']        = 10.0 # Coulomb interaction
P['Delta']    = 1.0	# superconduncting gap, usually an energy unit
P['GammaL']   = 0.5	# coupling to left electrode
P['GammaR']   = 0.5	# coupling to right electrode
P['GammaN']   = 0.0 # coupling to normal electrode
P['eps']      = 0.0	# local energy w.r.t. half-filled dot
P['P']        = 0.0	# phase difference
P['B']        = 0.0	# magnetic field
P['BW']       = 50.0 # half band-width
## integers in params
P['NFit']          = 8	# order of the tail fitting polynomial
P['NBins']         = 1
P['PrevBins']      = 0
## strings in params
P['InParam']       = 'U' # independent variable for output
## strings in ct-hyb
P['partition']     = 'autopartition' # Hamiltonian partitioning mathod
## bools in ct-hyb
P['measure_Gtau']  = True
P['measure_Gleg']  = False
P['measure_Gat']   = False
P['measure_pert']  = True
## integers in ct-hyb
P['NWarmup']       = int(1e5)
P['NCycles']       = int(1e7)
P['LengthCycle']   = 100
P['NMats']         = 400
P['NLeg']          = 40
P['NTau']          = 1001
## floats in ct-hyb
P['iw_cut']       = 100.0

## params part ############################################
for par in ['beta','U','Delta','GammaL','GammaR','GammaN','eps','P','B','BW']:
	if config.has_option('params',par):
	     P[par]  = float(config.get('params',par))
	else: print('# Warning: parameter '+str(par)+' not set, using default value: '+str(P[par]))
for par in ['NFit','NBins','PrevBins']:
	if config.has_option('params',par):
		P[par]  = int(config.get('params',par))
	else: print('# Warning: parameter '+str(par)+' not set, using default value: '+str(P[par]))
if config.has_option('params','InputParam'):
	P['InParam'] = str(config.get('params','InputParam'))
else: print('# Warning: parameter InParam not set, using default value: '+str(P['InParam']))

## ct-hyb part ############################################
if config.has_option('ct-hyb','partition'):
	P['partition'] = str(config.get('ct-hyb','partition'))
else: print('# Warning: parameter partition not set, using default value: '+str(P['partition']))
for par in ['measure_Gtau','measure_Gleg','measure_Gat','measure_pert']:
	if config.has_option('ct-hyb',par):
		P[par]    = bool(int(config.get('ct-hyb',par)))
	else: print('# Warning: parameter '+str(par)+' not set, using default value: '+str(P[par]))
for par in ['NWarmup','NCycles','LengthCycle','NMats','NLeg','NTau']:
	if config.has_option('ct-hyb',par):
		P[par]    = int(float(config.get('ct-hyb',par)))
	else: print('# Warning: parameter '+str(par)+' not set, using default value: '+str(P[par]))
if config.has_option('ct-hyb','iw_cut'):
	P['iw_cut']    = float(config.get('ct-hyb','iw_cut'))
else: print('# Warning: parameter iw_cut not set, using default value: '+str(P['iw_cut']))

## process the input ######################################
## just in case so we don't lose the statistics on G(tau):
if P['NBins'] > 1: P['measure_Gtau'] = True

## iw_cut is the cutoff in the Matsubara frequencies.
## If non-zero, the setting of NMats will be overriden and modified.
if P['iw_cut'] > 0.0:	# override of NMats
	NM = int(0.5*(P['beta']*P['iw_cut']/sp.pi)+1.0)
	if NM > P['NMats']: P['NMats'] = NM

## Fourier transform in TRIQS requires this:
if P['NTau'] <= 2*P['NMats']:
	print '# Warning: NTau too small, setting it to {0: 3d}'.format(2*P['NMats']+1)
	P['NTau'] = 2*P['NMats']+1

