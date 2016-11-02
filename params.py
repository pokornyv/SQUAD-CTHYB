import scipy as sp
from ConfigParser import SafeConfigParser

config = SafeConfigParser()
config.read('cthyb.in')

## params part ############################################
beta     = float(config.get('params','beta'))
U        = float(config.get('params','U'))
Delta    = float(config.get('params','Delta'))
GammaL   = float(config.get('params','GammaL'))
GammaR   = float(config.get('params','GammaR'))
GammaN   = float(config.get('params','GammaN'))
eps      = float(config.get('params','eps'))
P        = float(config.get('params','P'))
BW       = float(config.get('params','BW'))
NFit     = int(config.get('params','NFit'))
SymmS    = bool(int(config.get('params','SymmS')))
SymmG    = bool(int(config.get('params','SymmG')))
NBins    = int(config.get('params','NBins'))
PrevBins = int(config.get('params','PrevBins'))


## ct-hyb part ############################################
partition     = str(config.get('ct-hyb','partition'))

measure_Gtau  = bool(int(config.get('ct-hyb','measure_Gtau')))
measure_Gleg  = bool(int(config.get('ct-hyb','measure_Gleg')))
measure_Gat   = bool(int(config.get('ct-hyb','measure_Gat')))
measure_pert  = bool(int(config.get('ct-hyb','measure_pert')))

NWarmup       = int(float(config.get('ct-hyb','NWarmup')))
NCycles       = int(float(config.get('ct-hyb','NCycles')))
LengthCycle   = int(float(config.get('ct-hyb','LengthCycle')))

NMats         = int(config.get('ct-hyb','NMats'))
NLeg          = int(config.get('ct-hyb','NLeg'))
NTau          = int(config.get('ct-hyb','NTau'))

## process the input ######################################

if NBins > 1: measure_Gtau = True

if NTau <= 2*NMats:
	print '# Warning: NTau too small, setting it to {0: 3d}'.format(2*NMats+1)
	NTau = 2*NMats+1

