# SQUAD-CTHYB 
# CT-HYB QMC for a quantum dot connected to two superconducting and one normal leads 
# uses TRIQS cthyb as the solver, tested with TRIQS 1.4
# Vladislav Pokorny; 2016; pokornyv@fzu.cz

# mpirun -np 2 /storage/praha1/home/pokornyv/bin/triqs1.4/bin/pytriqs ../squad_cthyb.py

import scipy as sp
from numpy.random import randint
from time import time,ctime,sleep
from sys import argv,exit
from os import getcwd,environ,uname,listdir,mkdir

from pytriqs.gf.local import BlockGf,InverseFourier
import pytriqs.utility.mpi as mpi
from pytriqs.operators import *
from pytriqs.archive import *
from pytriqs.applications.impurity_solvers.cthyb import Solver

from squadlib_cthyb import *
try: 	# all threads are trying to read the same file at same time
	import params as p
except ConfigParser.NoSectionError:
	sleep(0.5)
	import params as p

mpi.barrier()

stars  = '**********************************************************'
hashes = '##########################################################'

if p.NBins == 1:		# single run
	logfname = 'cthyb.log'
else:					# binning mode
	logfname = 'bin.log'
	p.measure_Gtau = True	# just in case so we don't lose statistics on G(tau)

if mpi.is_master_node() and p.NBins > 1:
	if 'bin_files' not in listdir('.'): 
		try: mkdir('bin_files',0744)
		except OSError: exit()

MaxMats = OmegaN(p.NMats-1,p.beta) # Matsubaras begin with n = 0
FitMin  = int(0.9*p.NMats)         # index of frequency, not frequency itself
FitMax  = p.NMats
Phi     = p.P*sp.pi
#mu      = p.U/2.0

params_F = [p.beta,p.U,p.Delta,p.GammaL,p.GammaR,p.GammaN,p.eps,p.P]
bands_T  = ['up','dn']

if mpi.is_master_node() and p.NBins == 1:
	PrintAndWrite('\n'+stars+'\n* SQUAD CT-HYB START at '+ctime()+'\n'+stars+'\n',logfname)
	PrintAndWrite('Environment parameters:\n  Master node name: '\
	+str(uname()[1])+', OS '+str(uname()[0])+'-'+str(uname()[4]),logfname)
	if environ.has_key('PBS_NUM_NODES') and environ.has_key('PBS_NUM_PPN'):
		PrintAndWrite('  Number CPUs: '+str(environ['PBS_NUM_NODES'])\
		+' nodes @ '+str(environ['PBS_NUM_PPN'])+' cores/node\n',logfname)
	PrintAndWrite('  Working directory: '+str(getcwd())+'\n',logfname)
	PrintAndWrite('Model parameters:\n  U ={0: .4f}, eps ={1: .4f}, beta = {2: .4f}'
	.format(p.U,p.eps,p.beta),logfname)
	PrintAndWrite('  GammaL ={0: .4f}, GammaR ={1: .4f}, GammaN ={2: .4f}'\
	.format(p.GammaL,p.GammaR,p.GammaN),logfname)
	PrintAndWrite('  Delta ={0: .4f}, Phi/pi ={1: .4f}'.format(p.Delta,p.P),logfname)
	PrintAndWrite('Using {0: 3d} Matsubara frequencies, cutoff: {1: .4f}'.format(p.NMats,float(MaxMats)),logfname)
	PrintAndWrite('Half-bandwidth: {0: .4f}'.format(p.BW),logfname)

## construct the Green function and the solver ############
gf_struct = { '0': bands_T }

S = Solver(beta = p.beta, gf_struct = gf_struct, n_iw = p.NMats, n_l = p.NLeg, n_tau = p.NTau)

GF = GFzero(params_F,p.BW,p.NMats,FitMin,FitMax,p.NFit)
for block, gw0 in S.G0_iw: gw0 << GF

'''
if mpi.is_master_node():
	h1 = lambda iw: HybOffDiagFiniteW(params_F,p.BW,iw)
	h2 = lambda iw: HybDiagFiniteW(params_F,p.BW,iw)
	for iw in GF.mesh:
		H1 = h1(iw)
		H2 = h2(iw)
		print sp.real(iw),sp.imag(iw),sp.real(H1),sp.imag(H1),sp.real(H2),sp.imag(H2)
exit()
'''

if mpi.is_master_node(): 
	## write the non-interacting GF to file
	WriteTail(S.G0_iw,bands_T,logfname)
	WriteG_iw(S.G0_iw,bands_T,p.beta,p.NMats,'gw0',logfname)
	NtauZero_F = TotalDensity(S.G0_iw,S.G0_iw,bands_T,'mats')
	PrintAndWrite('Occupation matrix from G0(iw):',logfname)
	WriteMatrix(NtauZero_F,bands_T,'M',logfname)
	## transform non-interacting GF to imaginary times and write to file
	GTzero = S.G_tau.copy()
	GTzero.zero()
	for band1,band2 in product(bands_T, repeat=2):
		GTzero['0'][band1,band2] << InverseFourier(S.G0_iw['0'][band1,band2])
	WriteG_tau(GTzero,bands_T,p.beta,p.NTau,'gtau0',logfname)

## model Hamiltonian and other operators ##################
h_int = Operator()
N_tot = Operator()

h_int = -p.U*(n('0','up')-0.5)*(n('0','dn')-0.5)
N_tot = n('0','up') + n('0','dn')

## solver parameters ######################################
p_D = {}
p_D['h_int']                       = h_int
p_D['partition_method']            = p.partition
p_D['quantum_numbers']             = [N_tot]
p_D['n_cycles']                    = p.NCycles
p_D['length_cycle']                = p.LengthCycle
p_D['n_warmup_cycles']             = p.NWarmup
p_D['random_name']                 = ''
p_D['random_seed']                 = 123*mpi.rank + 567*int(time())
p_D['max_time']                    = -1
p_D['measure_g_l']                 = p.measure_Gleg
p_D['measure_g_tau']               = p.measure_Gtau
p_D['measure_pert_order']          = p.measure_pert
p_D['measure_density_matrix']      = False
p_D['use_norm_as_weight']          = True
p_D['move_shift']                  = True
p_D['move_double']                 = True
p_D['use_trace_estimator']         = False
p_D['imag_threshold']              = 1e-15

mpi.barrier()

if mpi.is_master_node():
	PrintAndWrite('\n'+hashes+'\nCT-HYB solver START',logfname)
	if p.measure_Gleg: PrintAndWrite('Measuring G(L).',logfname)
	if p.measure_Gtau: PrintAndWrite('Measuring G(tau).',logfname)

## run the solver #########################################
t = time()
for num_bin in range(p.PrevBins,p.PrevBins+p.NBins):
	if mpi.is_master_node() and p.NBins>1:
		PrintAndWrite('\n'+hashes+'\nRunning bin '+str(num_bin+1)+'/'+str(p.NBins+p.PrevBins),logfname)
	p_D['random_name'] = ''
	p_D['random_seed'] = randint(0,128)*mpi.rank + 567*int(time())
	S.solve(**p_D)
	mpi.barrier()
	if p.NBins>1:
		if mpi.is_master_node():
			if p.measure_Gtau:
				Ntau_F = TotalDensity2(S.G_tau['0'],bands_T,p.beta,p.NMats,'tau')
				JC_tau = JosephsonCurrent(S.G_tau['0'],params_F,p.BW,p.NMats,'R','tau')
				AppendBinFile(num_bin+1,Ntau_F,JC_tau,'bin_ntau.dat')
				WriteG_tau(S.G_tau,bands_T,p.beta,p.NTau,'bin_files/gtau_bin'+str(num_bin+1),logfname)
			if p.measure_Gleg:
				Nleg_F = TotalDensity2(S.G_l['0'],bands_T,p.beta,p.NMats,'leg')
				JC_leg = JosephsonCurrent(S.G_l['0'],params_F,p.BW,p.NMats,'R','leg')
				AppendBinFile(num_bin+1,Nleg_F,JC_tau,'bin_nleg.dat')
				WriteGleg(S.G_l,bands_T,p.beta,p.NLeg,'bin_files/gl_bin'+str(num_bin+1),logfname)
run_time = sp.around((time()-t)/60.0,2)

mpi.barrier()

if mpi.is_master_node():
	PrintAndWrite('CT-HYB solver STOP after {0: .2f} minutes = {1: .2f} hours.'\
	.format(run_time,run_time/60.0),logfname)
	PrintAndWrite('Solver status: {0: 3d}'.format(int(S.solve_status)),logfname)
	if p.NBins==1:
		PrintAndWrite('Number of cycles @ cycle length: {0: .2e} @{1: 3d}'\
		.format(p.NCycles,p.LengthCycle),logfname)
		PrintAndWrite('Average sign: {0: .3f}\n'.format(float(S.average_sign)),logfname)

		if p.measure_Gtau:
			WriteG_tau(S.G_tau,bands_T,p.beta,p.NTau,'gtau',logfname)
			Ntau_F = TotalDensity2(S.G_tau['0'],bands_T,p.beta,p.NMats,'tau')
			JC_tau = JosephsonCurrent(S.G_tau['0'],params_F,p.BW,p.NMats,'R','tau')
			PrintAndWrite('Total density from G(tau): {0: .6f}'.format(float(sp.trace(Ntau_F))),logfname)
			PrintAndWrite('Occupation matrix from G(tau):',logfname)
			WriteMatrix(Ntau_F,bands_T,'M',logfname)
			PrintAndWrite(':OUT_TAU\t{0: .6f}\t{1: .6f}\t{2: .6f}\t{3: .6f}\t{4: .6f}\t{5: .6f}\t{6: .6f}'\
			.format(float(Ntau_F[0][0]),float(Ntau_F[0][1]),float(Ntau_F[1][0]),float(Ntau_F[1][1])
			,float(Ntau_F[0][0]-Ntau_F[1][1]),float(Ntau_F[0][1]-Ntau_F[1][0]),float(JC_tau)),logfname)
			PrintAndWrite('\nJosephson current from G(tau): {0: .6f}\n'.format(float(JC_tau)),logfname)

		if p.measure_Gleg:
			WriteGleg(S.G_l,bands_T,p.beta,p.NLeg,'gl',logfname)
			Nleg_F = TotalDensity2(S.G_l['0'],bands_T,p.beta,p.NMats,'leg')
			JC_leg = JosephsonCurrent(S.G_l['0'],params_F,p.BW,p.NMats,'R','leg')
			PrintAndWrite('Total density from G(L): {0: .6f}'.format(float(sp.trace(Nleg_F))),logfname)
			PrintAndWrite('Occupation matrix from G(L):',logfname)
			WriteMatrix(Nleg_F,bands_T,'M',logfname)
			PrintAndWrite(':OUT_LEG\t{0: .6f}\t{1: .6f}\t{2: .6f}\t{3: .6f}\t{4: .6f}\t{5: .6f}\t{6: .6f}'\
			.format(float(Nleg_F[0][0]),float(Nleg_F[0][1]),float(Nleg_F[1][0]),float(Nleg_F[1][1])
			,float(Nleg_F[0][0]-Nleg_F[1][1]),float(Nleg_F[0][1]-Nleg_F[1][0]),float(JC_leg)),logfname)
			PrintAndWrite('\nJosephson current from G(L): {0: .6f}\n'.format(float(JC_leg)),logfname)

		if p.measure_pert:
			PO = S.perturbation_order_total
			WriteHisto(S.perturbation_order_total)
		WriteEig(S.h_loc_diagonalization,logfname)
	## write the results to a HDF5 archive for further processing
	R = HDFArchive('squad_cthyb.h5','w')
	R['g0_iw']  = S.G0_iw
	R['g0_tau'] = GTzero
	R['gtau']   = S.G_tau
	R['gl']     = S.G_l
	R['sgn']    = S.average_sign
	del R

if mpi.is_master_node():
	PrintAndWrite(argv[0]+' done, '+ctime(),logfname)

## end squad_cthyb.py

