Description of the *cthyb.in* file
==================================

The file *cthyb.in* defines various parameters used by the code.
Default values can be found in *params.py*, code that reads the input file.  
  
## Description

### [params] section

- beta       : inverse temperature, β = 1/kT  
- U          : Coulomb interaction strength  
- Delta      : superconducting gap  
- GammaL     : coupling to the left lead  
- GammaR     : coupling to the right lead  
- GammaN     : coupling to the normal lead  
- eps        : local energy level w.r.t e-h symmetric point (ε=0 is a half-filled dot)  
- P          : phase difference, in units of π  
- B          : external magnetic field  
- BW         : half-bandwidth of the non-interacting band  
- NFit       : number of terms fitted in tail fitting of Matsubara functions  
- InputParam : the independent parameter, used for output only  
- NBins      : number of bins in the binning procedure, NBins=1 means a single run  
- PrevBins   : number of bins in a previous run, useful if we continue binning  

### [ct-hyb] section

- partition     : Hilbert space partition method (autopartition is OK)  
- measure_Gtau  : measure G(tau) ? (0/1)  
- measure_Gleg  : measure G(L) ? (0/1)  
- measure_Gat   : measure atomic Green function ? (0/1)  
- measure_pert  : measure perturbation order histograms ? (0/1)  
- NWarmup       : number of warmup cycles  
- NCycles       : number of cycles  
- LengthCycle   : length of one cycle (longer cycle removes autocorrelation errors)  
- NMats         : number of Matsubara frequencies  
- NLeg          : number of Legendre polynomials  
- NTau          : number of imaginary time slices, TRIQS requires NTau > 2*NMats  
- iw_cut        : cutoff in Matsubara frequencies. It can override NMats is NMats is too small.  

