import os
import sys
import time
import numpy as np
import psi4

psi4.core.set_output_file('output.dat', False)

numpy_memory = 14

mol = psi4.geometry("""
O
H 1 1.1
H 1 1.1 2 104
symmetry c1
# """)
                    # H-He 	                Li-Ne                   	Na-Ar 
# cc-pVDZ   	[2s1p] → 5 func. 	    [3s2p1d] →   14 func. 	    [4s3p1d] → 18 func.
# cc-pVTZ   	[3s2p1d] → 14 func. 	[4s3p2d1f] → 30 func. 	    [5s4p2d1f] → 34 func.
# cc-pVQZ   	[4s3p2d1f] → 30 func. 	[5s4p3d2f1g] → 55 func. 	[6s5p3d2f1g] → 59 func.
# aug-cc-pVDZ 	[3s2p] → 9 func. 	    [4s3p2d] → 23 func.     	[5s4p2d] → 27 func.
# aug-cc-pVTZ 	[4s3p2d] → 23 func. 	[5s4p3d2f] → 46 func.   	[6s5p3d2f] → 50 func.
# aug-cc-pVQZ 	[5s4p3d2f] → 46 func. 	[6s5p4d3f2g] → 80 func. 	[7s6p4d3f2g] → 84 func. 
psi4.set_options({'basis': 'STO-6G', 'scf_type': 'pk', 'e_convergence': 1e-8, 'd_convergence': 1e-8})

print('\nStarting SCF and integral build...')

scf_e, wfn = psi4.energy('SCF', return_wfn=True)

C = wfn.Ca()

ndocc = wfn.doccpi()[0]
nmo = wfn.nmo()
nvirt = nmo - ndocc
print(ndocc,nmo,nvirt)
t = time.time()
print("may be energy :",psi4.energy('DETCI'))
print('\nPSI4 cisd cal: %.3f seconds.\n' % (time.time() - t))
# Compute size of Hamiltonian in GB
from scipy.special import comb
nDet_S = ndocc * nvirt * 2
nDet_D = 2 * comb(ndocc, 2) * comb(nvirt, 2) + ndocc**2 * nvirt**2
nDet = 1 + nDet_S + nDet_D
print("CISD nDet :" ,nDet)


nDet = comb(nmo, ndocc)**2

print("FCI nDet :" ,nDet)
