 
def run(x):
    import psi4
    import numpy as np
    # psi4.core.set_output_file('output.dat', False)
    psi4.core.be_quiet()
    psi4.set_memory('2 GB')
    numpy_memory = 14
    mol = psi4.geometry("""
    O
    symmetry c1
    """)
    # feb-cc-pcv_6pd_z
    psi4.set_options({'basis': 'cc-pvdz', 'scf_type': 'pk', 'e_convergence': 1e-8, 'd_convergence': 1e-8})
    scf_e, wfn = psi4.energy('SCF', return_wfn=True)
    C = wfn.Ca()
    # Number of molecular orbitals
    nmo = wfn.nmo()
    elec = wfn.nalpha() + wfn.nbeta()
    mints = psi4.core.MintsHelper(wfn.basisset())
    H = np.asarray(mints.ao_kinetic()) + np.asarray(mints.ao_potential())
    MO = np.asarray(mints.mo_spin_eri(C, C))
    # Vee = np.asarray(mints.mo_eri(C,C,C,C))
    H = np.einsum('uj,vi,uv', C, C, H)
    # Hone = H.copy()
    H = np.repeat(H, 2, axis=0)
    H = np.repeat(H, 2, axis=1)
    spin_ind = np.arange(H.shape[0], dtype=np.int64) % 2
    H *= (spin_ind.reshape(-1, 1) == spin_ind)
    psi4.core.clean()
    return H,MO,nmo,elec
