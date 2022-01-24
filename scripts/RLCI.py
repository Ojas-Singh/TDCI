import numpy as np
from tqdm import tqdm


def RL(A, k, mode='rl', learning_rate=0.5, discount=0.99, max_episode=30,
       max_pick=None,silent=False):
    """Returns the approximate eigenvalue and row indices for the k-sparse
       approximate solution to the lowest eigenpair of a matrix

    Parameters
    ----------
    A : numpy.ndarray
        The matrix to be approximately diagonalized
    k : int
        The maximum number of rows to retain in k-sparse apprx
    mode: str
        What type of diagonalization
        'rl': use reinforcement learning
        'apsci': a posteriori selected CI (top-k rows from full eigensolution)
        'greedy': greedy selection of top-k rows to include
    learning_rate: float
        Rate of weight update during learning. Choose between (0.0,1.0]
        Only makes sense with 'mode' == 'rl'
    discount : float
        How much to discount future gain in RL.
        Only makes sense with 'mode' == 'rl'
        Close to 0.0 : near-sighted learning
        Close to 1.0 : far-sighted learning
    max_episode: int
        Maximum number of training episodes
        Only makes sense with 'mode' == 'rl'
    max_pick: int
        Maximum number of external rows to consider during exploration
        If None, defaults to maximal space possible (not efficient!)
        Only makes sense with 'mode' == 'rl'

    Returns
    -------
    E_best: float
        Lowest eigenvalue from k-sparse approximation
    s: array_like
        Indices corresponding to rows for k-sparse approx solution of eigen

    """

    NDIM = len(A)

    if max_pick is None:
        max_pick = NDIM - k

    if mode == 'rl':
        # initialize with greedy; random init can be good choice too
        E, idx = RL(A, k, mode='greedy')
        #idx = np.random.choice(range(NDIM),size=k,replace=False) 

        # state vector is 1 if "active" SD, 0 otherwise
        state = np.zeros((NDIM), dtype=np.int)
        state[idx] = 1  # set topk rows from initialization to "active"

        # initial guess; note argwhere returns a 1D array, but we want 0D
        active = np.argwhere(state == 1).reshape(-1,)
        inactive = np.argwhere(state == 0).reshape(-1,)

        # check we conserve row indices
        assert len(active) == k
        assert len(inactive) == NDIM - k

        E, C = np.linalg.eigh(A[np.ix_(active, active)])
        E0, C0 = E[0], C[:, 0]

        # keep initial guess as our "best" so far
        E_best, state_best = E0, active

        # now populate the weights, renormalizing according to size
        w = np.zeros((NDIM))
        #w = np.random.randn(NDIM)

        w[active] = np.abs(C0)  # from initial guess
        w[active] /= np.linalg.norm(w[active])  # normalize
        #w[active] *= len(active)/NDIM  # scale

        w[inactive] = perturb(A, state)  # PT guess at weights
        w[inactive] /= np.linalg.norm(w[inactive])  # normalize
        #w[inactive] *= len(inactive)/NDIM  # scale

        v = np.zeros_like(w)  # auxilliary weights


        if not silent:
            progress_bar = tqdm(range(max_episode))  # keep track of progress
        else:
            progress_bar = range(max_episode)  

        # outer loop is training episode
        for episode in progress_bar:
            # current status printed out to command line

            if not silent:
                progress_bar.set_description("Best energy: %.6f" % E_best)
                #progress_bar.set_description("Current energy: %.6f" % E0)

            explore_rate = np.exp(-learning_rate*(episode+1))  # can tweak as desired

            # set state to be top-k weights from Q-learning
            state *= 0  # reset state vector

            # sometimes top k, sometimes best state reached so far
            if np.random.rand(1) < 0.75:
                #print("top")
                active = np.argsort(w)[::-1][:k].reshape(-1,)
            else:
                #print("best")
                active = state_best

            state[active] = 1
            inactive = np.argwhere(state == 0).reshape(-1,)

            # get energy from top-k
            E, C = np.linalg.eigh(A[np.ix_(active, active)])
            E0, C0 = E[0], C[:, 0]
            #print(E0)

            # check if top-k is globally optimal, save if it is
            if E0 < E_best:
                E_best = E0
                state_best = active

            # inner loop, partition into active search space
            # I'm using their jargon -- "selected" and "expanded"

            # "selected" are the lowest weighted rows in current selection
            active_idx_sorted = np.argsort(w[active])  # sort asc
            selected = np.arange(len(w))[active][active_idx_sorted]

            # "expanded" are the highest estimated rows outside of curr sele
            cj = perturb(A, state)  # PT guess at weights
            inact_idx_sort = np.argsort(np.abs(cj))[::-1]  # desc
            expanded = np.arange(len(w))[inactive][inact_idx_sort][:max_pick]

            # if we don't replace anything, we can't improve, so we exit
            total_replaced = 0
            for j in range(max_pick):
                assert sum(state) == k  # make sure we have consistent # row
                is_replaced = False
                state[expanded[j]] = 1
                for i in range(k):
                    state[selected[i]] = 0
                    active = np.argwhere(state == 1).reshape(-1,)
                    inactive = np.argwhere(state == 0).reshape(-1,)

                    # check we conserve row indices
                    assert len(active) == k
                    assert len(inactive) == NDIM - k

                    E, C = np.linalg.eigh(A[np.ix_(active, active)])
                    Enew, Cnew = E[0], C[:, 0]
                    if Enew < E0 * (1 - explore_rate*np.random.rand(1)):
                        #print("Update",Enew, E0)
                        is_replaced = True

                        total_replaced += 1

                        # note initial state
                        state1 = np.zeros_like(state)
                        state1[selected] = 1

                        # delete selected[i] from selected, add expanded[j]
                        # put [expanded[j]] to make 0-D array
                        p = selected[i]
                        q = expanded[j]

                        sele1, sele2 = np.split(selected, [i])
                        selected = np.concatenate((sele1,
                                                   sele2[1:],
                                                   [expanded[j]]))

                        # update state
                        state2 = np.zeros_like(state)
                        state2[selected] = 1
                        active = np.argwhere(state2 == 1).reshape(-1,)
                        inactive = np.argwhere(state2 == 0).reshape(-1,)

                        # uncomment to test that selection worked as expected
                        #E,C = np.linalg.eigh(A[np.ix_(selected,selected)])
                        #np.testing.assert_allclose(E[0],Enew)

                        # insert update here
                        # calculate local reward
                        R = E0 - Enew

                        # note if swap was globally optimal
                        if Enew < E_best:
                            E_best = Enew
                            state_best = active

                        E0 = Enew
                        C0 = Cnew

                        # get best possible next move (max Q(s',a'))
                        # a' = (p',q')
                        pp = np.arange(len(w))[active][np.argmin(w[active])]
                        qp = np.arange(len(w))[inactive][np.argmax(w[inactive])]

                        #print(p,q,pp,qp)

                        assert w[pp] == min(w[active])
                        assert w[qp] == max(w[inactive])

                        assert p not in active

                        # normal weight update; delta is TD error
                        delta = R + discount*np.dot(w,f(state2,(pp,qp))) - np.dot(w,f(state1,(p,q)))
                        aux = np.dot(f(state1,(p,q)),v) 

                        # update w
                        w += learning_rate *\
                             (delta*f(state1,(p,q)) - discount*aux*f(state2,(pp,qp)))

                        # update v 
                        beta = np.sqrt(learning_rate) 
                        v += beta*(delta - aux)*f(state1,(p,q))

                        break

                    # "p" from "selected" not selected: reset it
                    state[selected[i]] = 1

                # "q" from "expanded" not selected: reset it
                if is_replaced is False:
                    state[expanded[j]] = 0

            if total_replaced == 0:
                break

        return E_best, state_best

    elif mode == 'apsci':

        ## a posteriori selected CI
        ## full matrix diagonalization first, take top-k rows/columns
        E, C = np.linalg.eigh(A)
        idx = np.argpartition(np.abs(C[:, 0]), -k)[-k:]  # top-k eigenvec comps
        #
        ## form top-k submatrix and diagonalize
        A1 = A[np.ix_(idx, idx)]
        E_apsci = np.linalg.eigvalsh(A1)[0]

        return E_apsci, idx

    elif mode == 'greedy':
        # naive greedy selected CI
        state = np.zeros((NDIM))
        for i in range(k):
            if i == 0:
                # to begin, just take arbitrary element (not optimal)
                state[0] = 1
            # form top-k submatrix and diagonalize
            active = np.argwhere(state == 1).reshape(-1,)
            inactive = np.argwhere(state == 0).reshape(-1,)
            A1 = A[np.ix_(active, active)]
            E, C = np.linalg.eigh(A1)
            E0, C0 = E[0], C[:, 0]

            if i < k - 1:
                cj = perturb(A, state)
                cj_idx = np.argmax(np.abs(cj))
                state[np.arange(NDIM)[inactive][cj_idx]] = 1
            else:
                break

        E_greedy = E0
        idx = np.argwhere(state == 1)
        assert len(idx) == k
        return E_greedy, idx

    elif mode == 'full':
        # does full matrix diagonalization ("exact" result)
        return np.linalg.eigvalsh(A)[0], None

def f(state,a):
    active = np.argwhere(state == 1).reshape(-1,) 
    vector = np.zeros_like(state)
    p,q = a
    assert state[p] == 1
    assert state[q] == 0
    vector[active] = 1
    vector[q] =  1
    vector[p] = -1
    return vector/np.linalg.norm((vector))

def perturb(A, state):
    """

    Parameters
    ----------
    A : numpy.ndarray
        The matrix under consideration
    state : array_like
        Vector or list with '1' if index is active or '0' if inactive

    Returns
    -------
    c: array_like (dimension is number of "0's" in "state")
        Approximate magnitude of coefficients in the inactive space as obtained
        from (Epstein-Nesbet) perturbation theory

    """

    active = np.argwhere(state == 1).reshape(-1,)  # active row indices
    inactive = np.argwhere(state == 0).reshape(-1,)  # inactive row indices

    # diagonalize submatrix of active x active dimension
    E, C = np.linalg.eigh(A[np.ix_(active, active)])
    E0, C0 = E[0], C[:, 0]

    # obtain PT approximation to importance of inactive rows (determinants)
    a = np.minimum(1e5, np.abs(np.dot(A[np.ix_(inactive, active)], C0)))
    b = np.maximum(1e-5, np.abs(E0-np.diagonal(A[np.ix_(inactive, inactive)])))
    c = np.true_divide(a, b)
    return c


if __name__ == '__main__':
    np.random.seed(1)

    # generate random symmetric matrix
    # NDIM = 1000
    # A = np.random.rand(NDIM, NDIM)
    # A = A + A.T


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

    A= Hamiltonian_matrix

    k = 400  # sparsity
    print("k:    ", k)
    print("NDet: ", len(A))

    E_exact = RL(A, k, mode='full') 
    print("Exact:  ", E_exact)

    E_apsci, s = RL(A, k, mode='apsci')
    print("ap-sCI: ", E_apsci)

    E_greedy, s = RL(A, k, mode='greedy')
    print("greedy: ", E_greedy)

    E_rl, s = RL(A, k, mode='rl', max_pick=600)
    # E_rl, s = RL(A, k, mode='rl', max_pick=None)
    print("RL2:    ", E_rl,s)
    # s=[  16 , 17 ,  20   ,23   ,56  , 57  , 60  , 63,  492,  504 , 505 , 508 , 509,  518 , 519,  547 , 570  ,586,  624 , 625  ,663  ,700 , 703 , 820,  830,  860 , 870,  895,  898 , 899 , 902 , 940 , 950 , 953 , 980,  983 ,1014 ,1015 ,1022,1031, 1032 ,1050 ,1053, 1090 ,1093, 1106 ,1107, 1116, 1117, 1145 ,1149 ,1158 ,1171, 1185 ,1188 ,1198 ,1199 ,1210 ,1227 ,1250 ,1266 ,1298, 1299, 1301 ,1312 ,1343, 1383, 1421, 1545, 1585, 1610, 1626, 1627 ,1650 ,1666 ,1667 ,1819, 1822,1823 ,1839 ,1840, 1867 ,1868 ,1915, 1940 ,1941 ,1945 ,1948 ,2089, 2092 ,2093 ,2109 ,2110 ,2137, 2138 ,2185, 2210 ,2211 ,2215 ,2218]
    E, C = np.linalg.eigh(A[np.ix_(s, s)])
    for i in range(0,10):
        print("Energies ",i,"Our Ci :% 16.15f"% (E[i]),"Psi4   :% 16.15f"% (e_cisd[i]))
    
