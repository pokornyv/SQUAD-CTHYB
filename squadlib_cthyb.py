# SQUAD-CTHYB 
# CT-HYB QMC for a quantum dot connected to two superconducting and one normal leads 
# uses TRIQS cthyb as the solver, tested with TRIQS 1.4
# Vladislav Pokorny; 2016; pokornyv@fzu.cz

import scipy as sp
from time import ctime
from itertools import product
from pytriqs.gf.local import *
from pytriqs.gf.local.descriptors import Function
import params as p

###########################################################
## general functions ######################################

def OmegaN(n,beta):
	''' calculates the n-th fermionic Matsubara frequency '''
	return (2.0*n+1.0)*sp.pi/beta


def PrintAndWrite(line,fname):
	'''	print the same line to stdout and to file fname '''
	print(line)
	f = open(fname,'a')
	f.write(line+'\n')
	f.close()


###########################################################
## non-interacting Green function #########################

def SqrtFiniteBW(params_F,W,iw):
	''' square root factor enering hybridizations for finite bandwidth 2W '''
	[beta,U,Delta,GammaL,GammaR,GammaN,eps,P] = params_F
	sq = lambda x: sp.sqrt(Delta**2+sp.imag(x)**2)
	slog = lambda x: sp.log((1.0j*sq(iw)-x)/(1.0j*sq(iw)+x))
	return -1.0j/(2.0*sp.pi*sq(iw))*(slog(W)-slog(-W))
	

def HybDiagFiniteW(params_F,W,iw):
	''' diagonal part of the finite-bandwidth sc lead hybridization, real '''
	[beta,U,Delta,GammaL,GammaR,GammaN,eps,P] = params_F
	return (GammaL+GammaR)*SqrtFiniteBW(params_F,W,iw)


def HybOffDiagFiniteW(params_F,W,iw):
	''' off-diagonal part of the finite-bandwidth sc lead hybridization 
	PhiS angle keeps the hybridization real '''
	[beta,U,Delta,GammaL,GammaR,GammaN,eps,P] = params_F
	Phi = P*sp.pi
	PhiS = sp.arctan((GammaL-GammaR)/(GammaL+GammaR+1e-12)*sp.tan(Phi/2.0))
	return Delta*sp.exp(1.0j*PhiS)*\
	(GammaL*sp.exp(-1.0j*Phi/2.0)+GammaR*sp.exp(1.0j*Phi/2.0))*SqrtFiniteBW(params_F,W,iw)


def HybDiag(params_F,iw):
	''' diagonal part of the sc lead hybridization
	for infinite bandwidth, serves as a test for finite-bandwidth version '''
	[beta,U,Delta,GammaL,GammaR,GammaN,eps,P] = params_F
	return (GammaL+GammaR)/sp.sqrt(Delta**2+sp.imag(iw)**2)


def HybOffDiag(params_F,iw):
	''' off-diagonal part of the sc lead hybridization 
	for infinite bandwidth, serves as a test for finite-bandwidth version
	PhiS angle keeps the hybridization real '''
	[beta,U,Delta,GammaL,GammaR,GammaN,eps,P] = params_F
	Phi = P*sp.pi
	PhiS = sp.arctan((GammaL-GammaR)/(GammaL+GammaR+1e-12)*sp.tan(Phi/2.0))
	return Delta*sp.exp(1.0j*PhiS)*(GammaL*sp.exp(-1.0j*Phi/2.0)+GammaR*sp.exp(1.0j*Phi/2.0))\
	/sp.sqrt(Delta**2+sp.imag(iw)**2)

	
def GFzero(params_F,W,NMats,FitMin,FitMax,NTerms):
	''' constructs the non-interacting Green function as the input '''
	[beta,U,Delta,GammaL,GammaR,GammaN,eps,P] = params_F
	Phi = P*sp.pi
	V2 = 2.0*W*GammaN/sp.pi  # hybridization with normal lead
	## define lambdas (hybridizations are real)
	GF11 = lambda x: x*(1.0 + HybDiagFiniteW(params_F,W,x)) - eps
	GF12 = lambda x: HybOffDiagFiniteW(params_F,W,x)
	GF21 = lambda x: sp.conj(HybOffDiagFiniteW(params_F,W,x))
	GF22 = lambda x: x*(1.0 + HybDiagFiniteW(params_F,W,x)) + eps
	## define GF objects
	GFinv = GfImFreq(indices = ['up','dn'],beta = beta,n_points = NMats)
	GF    = GfImFreq(indices = ['up','dn'],beta = beta,n_points = NMats)
	## fill GF objects with lambdas
	GFinv['up','up'] << Function(GF11) - V2*Wilson(W)
	GFinv['up','dn'] << Function(GF12)
	GFinv['dn','up'] << Function(GF21)
	GFinv['dn','dn'] << Function(GF22) - V2*Wilson(W)
	## fit the tail
	GFinv.tail.zero()
	fixed_tail = TailGf(2,2,1,-1)
	fixed_tail[-1] = sp.eye(2)
	GFinv.fit_tail(fixed_tail,8,FitMin,FitMax)
	## calculate inverse and refit tail
	GF << inverse(GFinv)
	GF.tail.zero()
	fixed_tail = TailGf(2,2,3,-1)
	fixed_tail[-1] = sp.zeros([2,2])
	fixed_tail[ 0] = sp.zeros([2,2])
	fixed_tail[ 1] = sp.eye(2)
	GF.fit_tail(fixed_tail,8,FitMin,FitMax)
	return GF


###########################################################
## processing the output Green functions ##################

def TailCoeffs(G,bands_T):
	''' reads the tail of GF '''
	Gtail1_D = {}
	Gtail2_D = {}
	Gtail3_D = {}
	Gtail4_D = {}
	for band1,band2 in product(bands_T, repeat=2):
		Gtail1_D[band1,band2] = sp.real(G['0'][band1,band2].tail[1][0][0])
		Gtail2_D[band1,band2] = sp.real(G['0'][band1,band2].tail[2][0][0])
		Gtail3_D[band1,band2] = sp.real(G['0'][band1,band2].tail[3][0][0])
		Gtail4_D[band1,band2] = sp.real(G['0'][band1,band2].tail[4][0][0])
	return [Gtail1_D,Gtail2_D,Gtail3_D,Gtail4_D]


def TotalDensity(G,Gzero,bands_T,gtype):
	''' calculates the density from Green function '''
	G_iw = Gzero.copy()
	G_iw.zero()
	if gtype == 'leg':   G_iw['0'].set_from_legendre(G['0'])
	elif gtype == 'tau': G_iw['0'] = Fourier(G['0'])
	else:                G_iw['0'] = Gzero['0']  # gtype=='mats'
	N_F = sp.zeros([2,2])
	for i,j in product([0,1], repeat = 2):
		if gtype == 'leg': N_F[i][j] = sp.real(G['0'][bands_T[i],bands_T[j]].total_density())
		else:              N_F[i][j] = sp.real(G_iw['0'][bands_T[i],bands_T[j]].total_density())
	return N_F


def TotalDensity2(G,bands_T,beta,NMats,gtype):
	''' calculates the density from Green function '''
	N_F = sp.zeros([2,2])
	for i,j in product([0,1], repeat = 2):
		if gtype == 'leg': N_F[i][j] = sp.real(G[bands_T[i],bands_T[j]].total_density())
		elif gtype == 'tau':
			G_iw = GfImFreq(indices = bands_T,beta = beta,n_points = NMats)
			G_iw << Fourier(G)
			N_F[i][j] = sp.real(G_iw[bands_T[i],bands_T[j]].total_density())
		else:              N_F[i][j] = sp.real(G[bands_T[i],bands_T[j]].total_density())
	return N_F


def JosephsonCurrent(G,params_F,W,NMats,direction,gtype):
	''' Josephson current calculated from the Nambu Green function '''
	[beta,U,Delta,GammaL,GammaR,GammaN,eps,P] = params_F
	Phi = P*sp.pi
	## select direction, in theory JC_R = -JC_L
	Gamma = GammaR  if direction == 'R' else -GammaL
	Phi_s = Phi/2.0 if direction == 'R' else -Phi/2.0
	Hyb = lambda x: Gamma*SqrtFiniteBW(params_F,W,x)
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
	#for band1,band2 in product(['up','dn'], repeat=2):
	#	print J[band1,band2].density()
	JC = 2.0*sp.imag(J['up','dn'].density()[0][0]+J['dn','up'].density()[0][0])
	#print sp.imag(J['up','dn'].density()[0][0]),sp.imag(J['dn','up'].density()[0][0])
	return JC

###########################################################
## functions for writing data files #######################

def WriteG_iw(GF,bands_T,beta,NMats,fname,logfname):
	''' writes a Matsubara function to a file '''
	fout = open(fname,'w')
	G  = GfImFreq(indices = [0],beta = beta,n_points = NMats)
	MatsFreq_F = sp.zeros(2*NMats)
	k = 0
	for iw in G.mesh:		# filling the array of Matsubara frequencies
		MatsFreq_F[k] = sp.imag(iw)
		k += 1
	for i,j in product([0,1], repeat = 2):
		G << GF['0'][bands_T[i],bands_T[j]]
		for nw in range(2*NMats):
			fout.write('{0:.12f}\t{1:.12f}\t{2:.12f}\n'\
			.format(float(MatsFreq_F[nw]),float(sp.real(G.data[nw][0][0])),float(sp.imag(G.data[nw][0][0]))))
		fout.write('\n\n')
	fout.close()
	PrintAndWrite('File '+fname+' written.',logfname)


def WriteG_tau(GF,bands_T,beta,NTau,fname,logfname):
	''' writes an imaginary-time function to a file '''
	fout = open(fname,'w')
	G  = GfImTime(indices = [0],beta = beta,n_points = NTau)
	Times_F = sp.empty(NTau)
	k = 0
	for tau in G.mesh:	# filling the array of imaginary times
		Times_F[k] = sp.real(tau)
		k += 1	
	for i,j in product([0,1], repeat = 2):
		G << GF['0'][bands_T[i],bands_T[j]]
		for tau in range(NTau):
			fout.write('{0:.12f}\t{1:.12f}\t{2:.12f}\n'\
			.format(float(Times_F[tau]),float(sp.real(G.data[tau][0][0])),float(sp.imag(G.data[tau][0][0]))))
		fout.write('\n\n')
	fout.close()
	PrintAndWrite('File '+fname+' written.',logfname)


def WriteGleg(GF,bands_T,beta,NLeg,fname,logfname):
	''' writes GF in Legendre prepresentation to file '''
	fout = open(fname,'w')
	G = GfLegendre(indices = [0],beta = beta,n_points = NLeg)
	Leg_F = sp.array(range(NLeg))
	for i,j in product([0,1], repeat = 2):
		G << GF['0'][bands_T[i],bands_T[j]]
		for leg in Leg_F:
			fout.write('{0: 4d}\t{1:.12f}\n'.format(int(leg),float(sp.real(G.data[leg][0][0]))))
		fout.write('\n\n')
	fout.close()
	PrintAndWrite('File '+fname+' written.',logfname)


def WriteMatrix(X_D,bands_T,Xtype,logfname):
	''' writes a dict or matrix with two indices to output in a matrix form 
	input: Xtype = "D" = dict, "M" = matrix '''
	out_text = ''
	N = len(bands_T)
	for i in range(N):
		for j in range(N):
			X = X_D[bands_T[i],bands_T[j]] if Xtype == 'D' else X_D[i][j]
			if sp.imag(X) == 0: out_text += '{0: .6f}\t'.format(float(sp.real(X)))
			else: out_text += '{0: .6f}{1:+0.6f}i\t'.format(float(sp.real(X)),float(sp.imag(X)))
		out_text += '\n'
	if logfname != '': PrintAndWrite(out_text,logfname)
	else:              print out_text


def WriteTail(G,bands_T,logfname):
	''' writes the tail of the non-interacting GF to check for inconsistencies '''
	[Gtail1_D,Gtail2_D,Gtail3_D,Gtail4_D] = TailCoeffs(G,bands_T)
	PrintAndWrite('\n',logfname)	
	PrintAndWrite('G0(iw) tail fit: 1 / iw (-)',logfname)
	WriteMatrix(Gtail1_D,bands_T,'D',logfname)
	PrintAndWrite('G0(iw) tail fit: 1 / iw^2 (-) (local impurity levels)',logfname)
	WriteMatrix(Gtail2_D,bands_T,'D',logfname)
	PrintAndWrite('G0(iw) tail fit: 1 / iw^3 (+)',logfname)
	WriteMatrix(Gtail3_D,bands_T,'D',logfname)
	PrintAndWrite('G0(iw) tail fit: 1 / iw^4 (+)',logfname)
	WriteMatrix(Gtail4_D,bands_T,'D',logfname)


def AppendBinFile(nbin,Data_F,JC,filename):
	''' appends file with data from current bin '''
	f = open(filename,'a')
	if nbin == 1: f.write('# new calculation started on '+str(ctime())+'\n')
	f.write(str(nbin)+'\t')
	for i,j in product([0,1], repeat = 2):
		f.write('{0: .8f}\t'.format(float(sp.real(Data_F[i][j]))))
	f.write('{0: .8f}'.format(float(JC)))
	f.write('\n')
	f.close()


def WriteEig(eig,logfname):
	''' writes the eigenspectrum of the local Hamiltonian '''
	from string import zfill
	PrintAndWrite('Hamiltonian structure:',logfname)
	PrintAndWrite('  Hilbert space dimension:  {0: 3d}'.format(int(eig.full_hilbert_space_dim)),logfname)
	PrintAndWrite('  Number of blocks: {0: 3d}'.format(int(eig.n_blocks)),logfname)
	PrintAndWrite('  Energies:',logfname)
	for i in range(eig.n_blocks):
		for j in range(len(eig.energies[i])):
			PrintAndWrite('    {0: 3d}\t{1: 3d}\t{2: .8f}\t'\
			.format(i+1,j+1,eig.energies[i][j])+zfill(bin(eig.fock_states[i][j])[2:],2),logfname)
	PrintAndWrite('  Ground state energy: {0: .8f}'.format(float(eig.gs_energy)),logfname)


def WriteHisto(po):
	''' writes the histogram of perturbation order into file '''
	[low,high] = po.limits
	N = len(po.data)
	f = open('histo_total.dat','w')
	f.write('# Histogram of total perturbation order.\n# Written on '+ctime()+'\n')
	f.write('# limits: {0: 3d} - {1: 3d},\t'.format(int(po.limits[0]),int(po.limits[1])))
	f.write('number of data points: {0: 3d}\n'.format(int(po.n_data_pts)))	
	orders_F = sp.array(range(N))
	avg_pert = sp.sum(orders_F*po.data)/sum(po.data)
	sd_pert  = sp.sqrt(sp.sum((orders_F-avg_pert)**2*po.data)/sum(po.data))
	skew_pert = sp.sum((orders_F-avg_pert)**3/sd_pert**3*po.data)/sum(po.data)
	kurt_pert = sp.sum((orders_F-avg_pert)**4*po.data)/sum(po.data)/sd_pert**4-3
	f.write('# average pert. order: {0: .5f}, standard deviation: {1: .5f}\n'\
	.format(float(avg_pert),float(sd_pert)))
	f.write('# skew: {0: .5f}, excess kurtosis: {1: .5f}\n'\
	.format(float(skew_pert),float(kurt_pert)))
	for i in orders_F:
		f.write('{0: 3d}\t{1: 3d}\n'.format(int(i),int(po.data[i])))
	f.close()	

## end squadlib_cthyb.py

