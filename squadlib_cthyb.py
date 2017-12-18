# SQUAD-CTHYB
# CT-HYB QMC for a quantum dot connected to two superconducting and one normal leads 
# uses TRIQS cthyb as the solver, tested with TRIQS 1.4
# Vladislav Pokorny; 2016-2017; pokornyv@fzu.cz

import scipy as sp
from time import ctime
from itertools import product
from pytriqs.gf.local import *
from pytriqs.gf.local.descriptors import Function
from libqmc14 import *

###########################################################
## non-interacting Green function #########################
"""
def SqrtFiniteBW(Delta,W,iw):
	''' square root factor enering hybridizations for finite bandwidth 2W '''
	sq = lambda x: sp.sqrt(Delta**2+sp.imag(x)**2)
	slog = lambda x: sp.log((1.0j*sq(iw)-x)/(1.0j*sq(iw)+x))
	return -1.0j/(2.0*sp.pi*sq(iw))*(slog(W)-slog(-W))
"""

def IntFiniteBW(Delta,W,iw):
	''' integral enering hybridizations for finite bandwidth 2W '''
	sq = lambda x: sp.sqrt(Delta**2+sp.imag(x)**2)
	return 2.0/(sp.pi*sq(iw))*sp.arctan2(W,sq(iw))


def HybDiagFiniteW(GammaL,GammaR,Delta,W,iw):
	''' diagonal part of the finite-bandwidth sc lead hybridization, real 
	    caution, does not contain the iw_n factor from the matrix!!! '''
	return (GammaL+GammaR)*IntFiniteBW(Delta,W,iw)


def HybOffDiagFiniteW(GammaL,GammaR,Delta,P,W,iw):
	''' off-diagonal part of the finite-bandwidth sc lead hybridization 
	PhiS angle keeps the hybridization real '''
	Phi = P*sp.pi
	PhiS = sp.arctan((GammaL-GammaR)/(GammaL+GammaR+1e-12)*sp.tan(Phi/2.0))
	return Delta*sp.exp(1.0j*PhiS)*\
	(GammaL*sp.exp(-1.0j*Phi/2.0)+GammaR*sp.exp(1.0j*Phi/2.0))*IntFiniteBW(Delta,W,iw)


def GFzero(params_F,W,NMats,FitMin,FitMax,NTerms):
	''' constructs the non-interacting Green function as the input '''
	[beta,U,Delta,GammaL,GammaR,GammaN,eps,P,B] = params_F
	Phi = P*sp.pi
	epsUp =  eps + U/2.0 + B
	epsDn = -eps + U/2.0 + B
	V2 = W*GammaN/sp.pi  # hybridization with normal lead
	## define lambdas (hybridizations are real)
	GF11 = lambda x: x*(1.0 + HybDiagFiniteW(GammaL,GammaR,Delta,W,x)) - epsUp
	GF12 = lambda x: HybOffDiagFiniteW(GammaL,GammaR,Delta,P,W,x)
	GF21 = lambda x: sp.conj(HybOffDiagFiniteW(GammaL,GammaR,Delta,P,W,x))
	GF22 = lambda x: x*(1.0 + HybDiagFiniteW(GammaL,GammaR,Delta,W,x)) - epsDn
	## define GF objects
	GFinv = GfImFreq(indices = ['up','dn'],beta = beta,n_points = NMats)
	GF    = GfImFreq(indices = ['up','dn'],beta = beta,n_points = NMats)
	## fill GF objects with lambdas
	## W is half-bandwidth in Flat() descriptor from TRIQS
	GFinv['up','up'] << Function(GF11) - V2*Flat(W)
	GFinv['up','dn'] << Function(GF12)
	GFinv['dn','up'] << Function(GF21)
	GFinv['dn','dn'] << Function(GF22) - V2*Flat(W)
	## fit the tail
	GFinv.tail.zero()
	fixed_tail = TailGf(2,2,1,-1)
	fixed_tail[-1] = sp.eye(2)
	GFinv.fit_tail(fixed_tail,8,FitMin,FitMax)
	## calculate inverse
	GF << inverse(GFinv)
	## refit the tail
	GF.tail.zero()
	fixed_tail = TailGf(2,2,3,-1)
	fixed_tail[-1] = sp.zeros([2,2])
	fixed_tail[ 0] = sp.zeros([2,2])
	fixed_tail[ 1] = sp.eye(2)
	GF.fit_tail(fixed_tail,8,FitMin,FitMax)
	return GF

"""
def GFzeroRenorm(params_F,W,NMats,FitMin,FitMax,NTerms):
	''' constructs the non-interacting Green function as the input 
	uses the atomic limit renormalization procedure 
	normal lead not implemented!!! '''
	[beta,U,Delta,GammaL,GammaR,GammaN,eps,P,B] = params_F
	if GammaN!=0.0: print('Warning: Normal lead not implemented in renormalized atomic function!')
	Phi = P*sp.pi
	RF = 1.0+IntFiniteBW(Delta,W,iw)
	epsUp =  eps + U/2.0 + B
	epsDn = -eps + U/2.0 + B
	V2 = W*GammaN/sp.pi  # hybridization with normal lead
	## define lambdas (hybridizations are real)
	GF11_0 = lambda x: x*(1.0 + HybDiagFiniteW(GammaL,GammaR,Delta,W,x)) - epsUp
	GF12_0 = lambda x: HybOffDiagFiniteW(GammaL,GammaR,Delta,P,W,x)
	GF21_0 = lambda x: sp.conj(HybOffDiagFiniteW(GammaL,GammaR,Delta,P,W,x))
	GF22_0 = lambda x: x*(1.0 + HybDiagFiniteW(GammaL,GammaR,Delta,W,x)) - epsDn

	GF11_1 = lambda x: 
	GF12_1 = lambda x: 
	GF21_1 = lambda x: 
	GF22_1 = lambda x: 
	## define GF objects
	GFinv = GfImFreq(indices = ['up','dn'],beta = beta,n_points = NMats)
	GF    = GfImFreq(indices = ['up','dn'],beta = beta,n_points = NMats)
	## fill GF objects with lambdas
	## W is half-bandwidth in Flat() descriptor from TRIQS
	GFinv['up','up'] << Function(GF11_0)+Function(GF11_1)
	GFinv['up','dn'] << Function(GF12_0)+Function(GF12_1)
	GFinv['dn','up'] << Function(GF21_0)+Function(GF21_1)
	GFinv['dn','dn'] << Function(GF22_0)+Function(GF22_1)
	## fit the tail
	GFinv.tail.zero()
	fixed_tail = TailGf(2,2,1,-1)
	fixed_tail[-1] = sp.eye(2)
	GFinv.fit_tail(fixed_tail,8,FitMin,FitMax)
	## calculate inverse
	GF << inverse(GFinv)
	## refit the tail
	GF.tail.zero()
	fixed_tail = TailGf(2,2,3,-1)
	fixed_tail[-1] = sp.zeros([2,2])
	fixed_tail[ 0] = sp.zeros([2,2])
	fixed_tail[ 1] = sp.eye(2)
	GF.fit_tail(fixed_tail,8,FitMin,FitMax)
	return GF
"""

###########################################################
## processing the output Green functions ##################

def TotalDensity(G,Gzero,bands_T,gtype):
	''' calculates the density from Green function '''
	N = len(bands_T)
	G_iw = Gzero.copy()
	G_iw.zero()
	if gtype == 'leg':   G_iw['0'].set_from_legendre(G['0'])
	elif gtype == 'tau': G_iw['0'] = Fourier(G['0'])
	else:                G_iw['0'] = Gzero['0']  # gtype=='mats'
	N_F = sp.zeros([N,N])
	for i,j in product(range(N), repeat = 2):
		if gtype == 'leg': N_F[i][j] = sp.real(G['0'][bands_T[i],bands_T[j]].total_density())
		else:              N_F[i][j] = sp.real(G_iw['0'][bands_T[i],bands_T[j]].total_density())
	return N_F


def TotalDensity2(G,bands_T,beta,NMats,gtype):
	''' calculates the density from Green function '''
	N = len(bands_T)
	N_F = sp.zeros([N,N])
	for i,j in product(range(N), repeat = 2):
		if gtype == 'leg': N_F[i][j] = sp.real(G[bands_T[i],bands_T[j]].total_density())
		elif gtype == 'tau':
			G_iw = GfImFreq(indices = bands_T,beta = beta,n_points = NMats)
			G_iw << Fourier(G)
			N_F[i][j] = sp.real(G_iw[bands_T[i],bands_T[j]].total_density())
		else:              N_F[i][j] = sp.real(G[bands_T[i],bands_T[j]].total_density())
	return N_F


def JosephsonCurrent(G,params_F,W,NMats,direction,gtype):
	''' Josephson current calculated from the Nambu Green function '''
	[beta,U,Delta,GammaL,GammaR,GammaN,eps,P,B] = params_F
	Phi = P*sp.pi
	## select direction, in theory JC_R = -JC_L
	Gamma = GammaR  if direction == 'R' else -GammaL
	Phi_s = Phi/2.0 if direction == 'R' else -Phi/2.0
	Hyb = lambda x: Gamma*IntFiniteBW(Delta,W,x)
	## define empty GFs
	J    = GfImFreq(indices = ['up','dn'],beta = beta,n_points = NMats)
	H    = GfImFreq(indices = ['up','dn'],beta = beta,n_points = NMats)
	G_iw = GfImFreq(indices = ['up','dn'],beta = beta,n_points = NMats)
	## calculate the Matsubara function from input
	if   gtype == 'leg': G_iw.set_from_legendre(G)
	elif gtype == 'tau': G_iw.set_from_fourier(G)
	else:                G_iw = G.copy()
	## define lambdas and fill the hybridizations
	H11 = lambda x:  x*Hyb(x)
	H12 = lambda x:  Delta*sp.exp( 1.0j*Phi_s)*Hyb(x)
	H21 = lambda x: -Delta*sp.exp(-1.0j*Phi_s)*Hyb(x)
	H22 = lambda x: -x*Hyb(x)
	H['up','up'] << Function(H11)
	H['up','dn'] << Function(H12)
	H['dn','up'] << Function(H21)
	H['dn','dn'] << Function(H22)
	## fill the GF i<c{dag}d-d{dag}c>
	for band1,band2 in product(['up','dn'], repeat=2):
		J[band1,band2] << H[band1,band2]*G_iw[band1,band2]
	J.tail.zero()	## tail is small and hard to fit
	JC = 2.0*sp.imag(J['up','dn'].density()[0][0]+J['dn','up'].density()[0][0])
	return JC

###########################################################
## functions for writing data files #######################

def AppendBinFile(nbin,Data_F,JC,filename):
	''' appends file with data from current bin '''
	NBand = len(Data_F)
	f = open(filename,'a')
	if nbin == 1: f.write('# new calculation started on '+str(ctime())+'\n')
	f.write(str(nbin)+'\t')
	for i,j in product(range(NBand), repeat = 2):
		f.write('{0: .8f}\t'.format(float(sp.real(Data_F[i][j]))))
	f.write('{0: .8f}'.format(float(JC)))
	f.write('\n')
	f.close()

## end squadlib_cthyb.py

