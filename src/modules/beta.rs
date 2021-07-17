// Psi4Numpy COnfiguration Interaction Algorithm Implementation
use bit_array::BitArray;
use num::complex::Complex;
use rayon::prelude::*;
use std::sync::Mutex;
use typenum::U64;

#[path = "transform.rs"]
mod transform;

pub fn computeHamiltonianMatrix(
    binstates: Vec<BitArray<u64, U64>>,
    v: Vec<Vec<Vec<Vec<Complex<f64>>>>>,
    h: Vec<Vec<Complex<f64>>>,
    M: usize,
) -> Mutex<Vec<Vec<Complex<f64>>>> {
    let nslater = binstates.len();
    let hamiltonian = Mutex::new(vec![Vec::new(); nslater]);
    let hspin = transform::honetoSpin(M, &h);
    let vspin = transform::VeetoSpin(M, &v);
    (0..nslater).into_par_iter().for_each(|m| {
        let mut hm = vec![Complex::new(0.0, 0.0); nslater];
        for n in m..nslater {
            let hmn;
            let numuniqueorbitals: usize;
            numuniqueorbitals = totalDiff(&binstates[n], &binstates[m]);
            match numuniqueorbitals {
                0 => {
                    hmn = calcMatrixElementIdentialDet(&binstates[n], &vspin, &hspin);
                }
                1 => {
                    hmn = calcMatrixElementDiffIn1(&binstates[n], &binstates[m], &vspin, &hspin);
                }
                2 => {
                    hmn = calcMatrixElementDiffIn2(&binstates[n], &binstates[m], &vspin);
                }
                _ => {
                    hmn = Complex::new(0.0, 0.0);
                }
            }
            hm[n] = hmn;
        }
        let mut hamiltonian = hamiltonian.lock().unwrap();
        hamiltonian[m] = hm;
        
    });
    {
        let mut hamiltonian = hamiltonian.lock().unwrap();
        for i in 0..nslater {
            for j in i..nslater {
                if i != j {
                    hamiltonian[j][i] = hamiltonian[i][j]
                }
            }
        }
    }
    hamiltonian
}

pub fn signfun(n: usize, binstate: &BitArray<u64, U64>) -> f64 {
    let mut s = 1.0;
    for i in binstate.iter().take(n) {
        if i {
            s *= -1.0;
        }
    }
    s
}

fn getsign(state: &BitArray<u64, U64>, statelist: &Vec<usize>) -> f64 {
    let mut sign = 1.0;
    for i in statelist {
        if state[*i] {
            sign *= signfun(*i, state)
        }
    }
    sign
}

fn getUniqueOrbitalsInMixIndexListsPlusSign(
    state1: &BitArray<u64, U64>,
    state2: &BitArray<u64, U64>,
) -> (Vec<usize>, Vec<usize>, f64) {
    let mut unique1 = Vec::new();
    let mut unique2 = Vec::new();
    let mut state = state1.clone();
    let _common = BitArray::intersect(&mut state, &state2);
    for i in num_iter::range(0, state.len()) {
        if state[i] ^ state1[i] {
            unique1.push(i)
        }
        if state[i] ^ state2[i] {
            unique2.push(i)
        }
    }
    let sign = getsign(state1, &unique1) * getsign(state2, &unique2);

    return (unique1, unique2, sign);
}
fn getOrbitalMixedIndexList(state: &BitArray<u64, U64>) -> Vec<usize> {
    let mut list = Vec::new();
    for (j, i) in state.iter().enumerate() {
        if i {
            list.push(j)
        }
    }
    list
}

fn getCommonOrbitalsInMixedSpinIndexList(
    state1: &BitArray<u64, U64>,
    state2: &BitArray<u64, U64>,
) -> Vec<usize> {
    let mut list = Vec::new();
    let mut state = state1.clone();
    let _common = BitArray::intersect(&mut state, &state2);
    for i in 0..state.len() {
        if state[i] {
            list.push(i)
        }
    }
    list
}

fn totalDiff(state: &BitArray<u64, U64>, state2: &BitArray<u64, U64>) -> usize {
    let mut diff: usize = 0;
    for i in num_iter::range(0, state.len()) {
        if state[i] ^ state2[i] {
            diff += 1
        }
    }
    diff / 2
}

fn calcMatrixElementIdentialDet(
    state1: &BitArray<u64, U64>,
    v: &Vec<Vec<Vec<Vec<Complex<f64>>>>>,
    h: &Vec<Vec<Complex<f64>>>
) -> Complex<f64> {
    let mut helem = Complex::new(0.0, 0.0);
    let mut relem = Complex::new(0.0, 0.0);
    let list = getOrbitalMixedIndexList(&state1);
    for m in &list {
        helem += h[*m][*m]
    }

    for m in num_iter::range(0, list.len() - 1) {
        for n in num_iter::range(m + 1, list.len()) {
            relem += v[list[m]][list[n]][list[m]][list[n]]
        }
    }
    helem + relem
}

fn calcMatrixElementDiffIn1(
    state1: &BitArray<u64, U64>,
    state2: &BitArray<u64, U64>,
    v: &Vec<Vec<Vec<Vec<Complex<f64>>>>>,
    h: &Vec<Vec<Complex<f64>>>
) -> Complex<f64> {
    let (unique1, unique2, sign) = getUniqueOrbitalsInMixIndexListsPlusSign(state1, state2);
    let m = unique1[0];
    let p = unique2[0];
    let helem = h[m][p];
    let common = getCommonOrbitalsInMixedSpinIndexList(state1, state2);
    let mut relem = Complex::new(0.0, 0.0);
    for n in common {
        relem += v[m][n][p][n]
    }
    sign * (helem + relem)
}

fn calcMatrixElementDiffIn2(
    state1: &BitArray<u64, U64>,
    state2: &BitArray<u64, U64>,
    v: &Vec<Vec<Vec<Vec<Complex<f64>>>>>,
) -> Complex<f64> {
    let (unique1, unique2, sign) = getUniqueOrbitalsInMixIndexListsPlusSign(state1, state2);
    sign * v[unique1[0]][unique1[1]][unique2[0]][unique2[1]]
}
