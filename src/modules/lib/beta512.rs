// Psi4Numpy COnfiguration Interaction Algorithm Implementation
use bit_array::BitArray;
use num::complex::Complex;
use rayon::prelude::*;
// use std::sync::Mutex;
// use std::time::Instant;
use typenum::U512;
// use arrayfire::*;

pub fn computeHamiltonianMatrix(
    binstates: Vec<BitArray<u64, U512>>,
    vspin: Vec<Vec<Vec<Vec<Complex<f64>>>>>,
    hspin: Vec<Vec<Complex<f64>>>,
    M: usize,
) -> Vec<Vec<Complex<f64>>> {
    let nslater = binstates.len();
    // let (hspin, vspin) = transform::tospin(M, &h, &v);

    let hamiltonian = (0..nslater)
        .into_par_iter()
        .map(|m| {
            let mut hm = vec![Complex::new(0.0, 0.0); nslater];
            for n in m..nslater {
                let mut hmn = Complex::new(0.0, 0.0);
                let numuniqueorbitals: usize = totalDiff(&binstates[n], &binstates[m], M);
                match numuniqueorbitals {
                    0 => {
                        hmn = calcMatrixElementIdentialDet(&binstates[n], &vspin, &hspin, M);
                    }
                    1 => {
                        hmn = calcMatrixElementDiffIn1(
                            &binstates[n],
                            &binstates[m],
                            &vspin,
                            &hspin,
                            M,
                        );
                    }
                    2 => {
                        hmn = calcMatrixElementDiffIn2(&binstates[n], &binstates[m], &vspin, M);
                    }
                    _ => {}
                }
                hm[n] = hmn;
            }
            hm
        })
        .collect();

    // {
    //     let mut hamiltonian = hamiltonian.lock().unwrap();
    //     for i in 0..nslater {
    //         for j in i..nslater {
    //             if i != j {
    //                 hamiltonian[j][i] = hamiltonian[i][j]
    //             }
    //         }
    //     }
    // }
    hamiltonian
}
#[inline(always)]
pub fn signfun(n: usize, binstate: &BitArray<u64, U512>) -> f64 {
    let s = binstate.iter().take(n).filter(|x| *x).count();
    if s % 2 == 0 {
        return 1 as f64;
    } else {
        return -1 as f64;
    }
}
#[inline(always)]
fn getsign(state: &BitArray<u64, U512>, statelist: &Vec<usize>) -> f64 {
    let mut sign = 1.0;
    for i in statelist {
        if state[*i] {
            sign *= signfun(*i, state)
        }
    }
    sign
}

#[inline(always)]
fn xor(state1: &BitArray<u64, U512>, state2: &BitArray<u64, U512>) -> BitArray<u64, U512> {
    let mut state11 = state1.clone();
    let mut state22 = state2.clone();
    state11.difference(&state2);
    state22.difference(&state1);
    state11.union(&state22);
    state11
}

#[inline(always)]
fn getUniqueOrbitalsInMixIndexListsPlusSign(
    state1: &BitArray<u64, U512>,
    state2: &BitArray<u64, U512>,
    M: usize,
) -> (Vec<usize>, Vec<usize>, f64) {
    let mut unique1 = Vec::with_capacity(2);
    let mut unique2 = Vec::with_capacity(2);
    let mut state = state1.clone();
    let _common = BitArray::intersect(&mut state, &state2);
    let state11 = xor(&state, state1);
    let state22 = xor(&state, state2);

    for (j, i) in state11.iter().take(2 * M + 1).enumerate() {
        if i {
            unique1.push(j)
        }
    }
    for (j, i) in state22.iter().take(2 * M + 1).enumerate() {
        if i {
            unique2.push(j)
        }
    }

    let sign = getsign(state1, &unique1) * getsign(state2, &unique2);

    return (unique1, unique2, sign);
}
#[inline(always)]
fn getOrbitalMixedIndexList(state: &BitArray<u64, U512>, M: usize) -> Vec<usize> {
    let mut list = Vec::new();
    for (j, i) in state.iter().take(2 * M + 1).enumerate() {
        if i {
            list.push(j)
        }
    }
    list
}
#[inline(always)]
fn getCommonOrbitalsInMixedSpinIndexList(
    state1: &BitArray<u64, U512>,
    state2: &BitArray<u64, U512>,
    M: usize,
) -> Vec<usize> {
    let mut list = Vec::new();
    let mut state = state1.clone();
    let _common = BitArray::intersect(&mut state, &state2);
    for i in 0..(2 * M) {
        if state[i] {
            list.push(i)
        }
    }
    list
}
#[inline(always)]
fn totalDiff(state: &BitArray<u64, U512>, state2: &BitArray<u64, U512>, M: usize) -> usize {
    let mut state11 = state.clone();
    let mut state22 = state2.clone();
    state11.difference(&state2);
    state22.difference(state);
    state11.union(&state22);
    state11.iter().take(2 * M + 1).filter(|x| *x).count() / 2
}
// #[inline(always)]
fn calcMatrixElementIdentialDet(
    state1: &BitArray<u64, U512>,
    v: &Vec<Vec<Vec<Vec<Complex<f64>>>>>,
    h: &Vec<Vec<Complex<f64>>>,
    M: usize,
) -> Complex<f64> {
    let mut helem = Complex::new(0.0, 0.0);
    let mut relem = Complex::new(0.0, 0.0);
    let list = getOrbitalMixedIndexList(&state1, M);
    for m in &list {
        helem += h[*m][*m]
    }

    for m in num_iter::range(0, list.len() - 1) {
        for n in num_iter::range(m + 1, list.len()) {
            relem += unsafe {
                v.get_unchecked(list[m])
                    .get_unchecked(list[n])
                    .get_unchecked(list[m])
                    .get_unchecked(list[n])
            }
        }
    }
    helem + relem
}
#[inline(always)]
fn calcMatrixElementDiffIn1(
    state1: &BitArray<u64, U512>,
    state2: &BitArray<u64, U512>,
    v: &Vec<Vec<Vec<Vec<Complex<f64>>>>>,
    h: &Vec<Vec<Complex<f64>>>,
    M: usize,
) -> Complex<f64> {
    let (unique1, unique2, sign) = getUniqueOrbitalsInMixIndexListsPlusSign(state1, state2, M);
    let m = unique1[0];
    let p = unique2[0];
    let helem = h[m][p];
    let common = getCommonOrbitalsInMixedSpinIndexList(state1, state2, M);
    let mut relem = Complex::new(0.0, 0.0);
    for n in common {
        relem += unsafe {
            v.get_unchecked(m)
                .get_unchecked(n)
                .get_unchecked(p)
                .get_unchecked(n)
        }
    }
    sign * (helem + relem)
}
#[inline(always)]
fn calcMatrixElementDiffIn2(
    state1: &BitArray<u64, U512>,
    state2: &BitArray<u64, U512>,
    v: &Vec<Vec<Vec<Vec<Complex<f64>>>>>,
    M: usize,
) -> Complex<f64> {
    let (unique1, unique2, sign) = getUniqueOrbitalsInMixIndexListsPlusSign(state1, state2, M);
    // println!("{}",unique1.len());
    sign * unsafe {
        v.get_unchecked(unique1[0])
            .get_unchecked(unique1[1])
            .get_unchecked(unique2[0])
            .get_unchecked(unique2[1])
    }
}
// v.get_unchecked(unique1.get_unchecked(0)).get_unchecked(unique1.get_unchecked(1)).get_unchecked(unique2.get_unchecked(0)).get_unchecked(unique2.get_unchecked(1))
