SQUAD-CTHYB
===========
#### Description:
**SQUAD-CTHYB** is a set of codes to solve the single-impurity Anderson model connected to two BCS superconducting leads
and optionally one normal lead using the continuous-time, hybridization-expansion quantum Monte Carlo solver included in the TRIQS package [Seth2016, Parcollet2015]. Developed and tested using TRIQS 1.4. The implementation is motivated by the CT-INT method described in [Luitz2010].  

Some results obtained using this code were published in [Domański2017,Pokorny2017].  

These codes are subject to big changes and require some optimization and cleanup. Please do not use for other than educational purposes unless you really know what you are doing.  

The code uses functions from *limqmc14.py* library that is included in different repository (https://github.com/pokornyv/libQMC-TRIQS). 

Home: https://github.com/pokornyv/SQUAD-CTHYB  

#### List of files:
- *squad_cthyb.py* - main code  
- *squadlib_cthyb.py* - library of functions  
- *params.py* - script to read the configuration file *cthyb.in*  
- *cthyb.in* - configuration file  
- *infile.md* - description of the *cthyb.in* file  

#### References:
- T. Domański, M. Žonda, V. Pokorný, G. Górski, V. Janiš, and T. Novotný, *Phys. Rev. B* **95**, 045104 (2017).
- V. Pokorný and M. Žonda, *arXiv:1706.08783*, (2017).
- D. J. Luitz, and F. F. Assaad, *Phys. Rev. B* **81**, 024509 (2010).  
- P. Seth, I. Krivenko, M. Ferrero, and O. Parcollet, *Comp. Phys. Commun.* **200**, 274 (2016).  
- O. Parcollet, M. Ferrero, T. Ayral, H. Hafermann, I. Krivenko, L. Messio, and P. Seth, *Comput. Phys. Commun.* **196**, 398 (2015).  

