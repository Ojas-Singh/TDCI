use bit_array::BitArray;
use num::complex::Complex;
use rayon::prelude::*;
use std::sync::Mutex;
use typenum::U64;

pub fn sign(n: usize, binstate: &BitArray<u64, U64>) -> f64 {
    let mut s = 1.0;
    for i in binstate.iter().take(n) {
        if i {
            s *= -1.0;
        }
    }
    s
}

pub fn addparticle(n: usize, binstate: &mut BitArray<u64, U64>) {
    if binstate.get(binstate.len() - 1) != Some(true) {
        let mut a = BitArray::<u64, U64>::from_elem(false);
        a.set(n, true);
        let comp = a.intersect(&binstate);
        if comp {
            binstate.set(n, true);
        } else {
            binstate.set(binstate.len() - 1, true);
        }
    }
}

pub fn removeparticle(n: usize, binstate: &mut BitArray<u64, U64>) {
    if binstate.get(binstate.len() - 1) != Some(true) {
        let mut a = BitArray::<u64, U64>::from_elem(false);
        a.set(n, true);
        let comp = a.intersect(&binstate);
        if !comp {
            binstate.set(n, false);
        } else {
            binstate.set(binstate.len() - 1, true);
        }
    }
}

pub fn secondQuantizationOneBodyOperator(
    p: usize,
    q: usize,
    state: &BitArray<u64, U64>,
    state2: &BitArray<u64, U64>,
) -> f64 {
    let mut phase = 1.0;
    let mut state1 = state.clone();
    let k = state1.len() - 1;
    removeparticle(q, &mut state1);
    if state1.get(k) == Some(true) {
        return 0.0;
    }
    phase *= sign(q, &state1);
    addparticle(p, &mut state1);
    if state1.get(k) == Some(true) {
        return 0.0;
    }
    phase *= sign(p, &state1);
    if &state1 != state2 {
        return 0.0;
    }
    phase
}

pub fn secondQuantizationTwoBodyOperator(
    p: usize,
    q: usize,
    r: usize,
    s: usize,
    state: &BitArray<u64, U64>,
    state2: &BitArray<u64, U64>,
) -> f64 {
    let mut phase = 1.0;
    let mut state1 = state.clone();
    let k = state1.len() - 1;

    removeparticle(r, &mut state1);
    if state1.get(k) == Some(true) {
        return 0.0;
    }
    phase *= sign(r, &state1);
    removeparticle(s, &mut state1);
    if state1.get(k) == Some(true) {
        return 0.0;
    }
    phase *= sign(s, &state1);
    addparticle(q, &mut state1);
    if state1.get(k) == Some(true) {
        return 0.0;
    }
    phase *= sign(q, &state1);
    addparticle(p, &mut state1);
    if state1.get(k) == Some(true) {
        return 0.0;
    }
    phase *= sign(p, &state1);
    if &state1 != state2 {
        return 0.0;
    }
    phase
}

pub fn computeHamiltonianMatrix(
    binstates: Vec<BitArray<u64, U64>>,
    v: Vec<Vec<Vec<Vec<Complex<f64>>>>>,
    h: Vec<Vec<Complex<f64>>>,
    M: usize,
) -> Mutex<Vec<Vec<Complex<f64>>>> {
    let nslater = binstates.len();
    let hamiltonian = Mutex::new(vec![vec![Complex::new(0.0, 0.0); nslater]; nslater]);
    (0..nslater).into_par_iter().for_each(|m| {
        for n in m..nslater {
            let mut hmn = Complex::new(0.0, 0.0);
            for p in 0..M {
                for q in 0..M {
                    let mut phase = 0.0;
                    phase += secondQuantizationOneBodyOperator(
                        2 * p,
                        2 * q,
                        &binstates[n],
                        &binstates[m],
                    );
                    phase += secondQuantizationOneBodyOperator(
                        2 * p + 1,
                        2 * q + 1,
                        &binstates[n],
                        &binstates[m],
                    );
                    hmn += h[p][q] * phase;

                    for r in 0..M {
                        for s in 0..M {
                            let V = v[p][r][q][s];
                            if V != Complex::new(0.0, 0.0) {
                                phase = 0.0;
                                phase += secondQuantizationTwoBodyOperator(
                                    2 * p,
                                    2 * q,
                                    2 * r,
                                    2 * s,
                                    &binstates[n],
                                    &binstates[m],
                                );
                                phase += secondQuantizationTwoBodyOperator(
                                    2 * p + 1,
                                    2 * q + 1,
                                    2 * r + 1,
                                    2 * s + 1,
                                    &binstates[n],
                                    &binstates[m],
                                );
                                phase += secondQuantizationTwoBodyOperator(
                                    2 * p,
                                    2 * q + 1,
                                    2 * r,
                                    2 * s + 1,
                                    &binstates[n],
                                    &binstates[m],
                                );
                                phase += secondQuantizationTwoBodyOperator(
                                    2 * p + 1,
                                    2 * q,
                                    2 * r + 1,
                                    2 * s,
                                    &binstates[n],
                                    &binstates[m],
                                );
                                hmn += 0.5 * V * phase;
                            }
                        }
                    }
                }
            }
            let mut hamiltonian = hamiltonian.lock().unwrap();
            hamiltonian[m][n] = hmn;
            if m != n {
                hamiltonian[n][m] = hamiltonian[m][n];
            }
        }
    });
    hamiltonian
}
