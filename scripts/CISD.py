
import time
import numpy as np
np.set_printoptions(precision=5, linewidth=200, suppress=True)
import psi4

psi4.core.set_output_file('output.dat', False)

numpy_memory = 2

mol = psi4.geometry("""
O
symmetry c1
""")

psi4.set_options({'basis': 'cc-pVDZ', 'scf_type': 'pk', 'e_convergence': 1e-8, 'd_convergence': 1e-8})

print('\nStarting SCF and integral build...')

scf_e, wfn = psi4.energy('SCF', return_wfn=True)

C = wfn.Ca()

ndocc = wfn.doccpi()[0]
print(wfn.nalpha())
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
H_Size = nDet**2 * 8e-9
print("nDet :" ,nDet)
if H_Size > numpy_memory:
    clean()
    raise Exception("Estimated memory utilization (%4.2f GB) exceeds numpy_memory \
                    limit of %4.2f GB." % (H_Size, numpy_memory))

mints = psi4.core.MintsHelper(wfn.basisset())
H = np.asarray(mints.ao_kinetic()) + np.asarray(mints.ao_potential())


MO = np.asarray(mints.mo_spin_eri(C, C))
Vee = np.asarray(mints.mo_eri(C,C,C,C))

H = np.einsum('uj,vi,uv', C, C, H)
Hone = H.copy()
H = np.repeat(H, 2, axis=0)
H = np.repeat(H, 2, axis=1)

spin_ind = np.arange(H.shape[0], dtype=np.int) % 2
H *= (spin_ind.reshape(-1, 1) == spin_ind)


from helper_CI import Determinant, HamiltonianGenerator
from itertools import combinations

print('Generating %d CISD Determinants...' % (nDet))

occList = [i for i in range(ndocc)]

det_ref = Determinant(alphaObtList=occList, betaObtList=occList)

detList = det_ref.generateSingleAndDoubleExcitationsOfDet(nmo)
detList.append(det_ref)


# for i in detList :
#     print("detList",type(i))



print('Generating Hamiltonian Matrix...')

t = time.time()
Hamiltonian_generator = HamiltonianGenerator(H, MO)
Hamiltonian_matrix = Hamiltonian_generator.generateMatrix(detList)

print('..finished generating Matrix in %.3f seconds.\n' % (time.time() - t))



e_cisd, wavefunctions = np.linalg.eigh(Hamiltonian_matrix)

cisd_mol_e = e_cisd[0] + mol.nuclear_repulsion_energy()

print('# Determinants:      % 16d' % (len(detList)))

print('SCF energy:          % 16.10f' % (scf_e))
print('CISD correlation:    % 16.10f' % (cisd_mol_e - scf_e))
print('Total CISD energy:   % 16.10f' % (cisd_mol_e))



print('Exporting integrals to C.txt , Hone.txt , Vee.txt')
f = open("C.txt", "w")
cmat= np.array(C)
for i in cmat:
    for j in i:
        f.write(str(j)+"\n")
f.close()
f = open("Hone.txt", "w")
for i in Hone:
    for j in i:
        f.write(str(j)+",0.0"+"\n")
f.close()
f = open("Vee.txt", "w")
for i in Vee:
    for j in i:
        for k in j:
            for l in k:
                f.write(str(l)+",0.0"+"\n")
f.close()


print('Sending Integrals to Our Ci Code and Building Hamiltonian')
import subprocess
import sys
t = time.time()
subprocess.run(["./../target/release/TDCI","8","14","Singlet","Hone.txt","Vee.txt","2"])
print('..finished generating Matrix in %.3f seconds.\n' % (time.time() - t))
print('Importing Hamiltonian from Our Ci Code back to Python for digonalization.')
f = open("ham.txt", "r")
Lines = f.readlines()
hman =[]
for i in Lines:
    # hman.append(float(i))
    hman.append(float(i.split('+')[0]))

Hmat=np.zeros((int(nDet),int(nDet)))
line=0
for i in range(int(nDet)):
    for j in range(int(nDet)):
        Hmat[i,j]+=hman[line]
        line+=1
# Hmat = Hmat + Hmat.T - Hmat.diagonal(1)
# w,v = np.linalg.eigh(Hmat,UPLO='L')
b_symm = (Hmat + Hmat.T)
np.fill_diagonal(b_symm, b_symm.diagonal()/2)
w,v = np.linalg.eigh(b_symm)
# print(Hmat)
print("Energy By Our Ci code :",w[0])
print("Psi4 Eigen Value      :",e_cisd[0] )
print(" cisd_mol_e : ",cisd_mol_e)
for i in range(0,10):
    print("Energies ",i,"Our Ci :% 16.15f"% (w[i]),"Psi4   :% 16.15f"% (e_cisd[i]))

# import scipy.sparse as sparse
# import matplotlib.pyplot as plt
# plt.figure(figsize=(15, 15))
# plt.spy(Hmat, markersize=1)

# # plt.matshow(Hmat)

# plt.savefig('temp.png')

