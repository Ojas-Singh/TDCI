// Psi4Numpy COnfiguration Interaction Algorithm Implementation
use bit_array::BitArray;
use num::complex::Complex;
use rayon::prelude::*;
use std::sync::Mutex;
use typenum::U64;
#[path = "transform.rs"]
mod transform;

#[inline(always)]
pub fn sign(n: usize, binstate: &BitArray<u64, U64>) -> f64 {
    let mut s = 1.0;
    for i in binstate.iter().take(n) {
        if i {
            s *= -1.0;
        }
    }
    s
}
#[inline(always)]
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
#[inline(always)]
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
#[inline(always)]
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
#[inline(always)]
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

pub fn calcMatrixElement(state1: &BitArray<u64, U64>,state2: &BitArray<u64, U64>,v: &Vec<Vec<Vec<Vec<Complex<f64>>>>>,
    h: &Vec<Vec<Complex<f64>>>,
    M: usize) -> Complex<f64> {
    let mut hmn = Complex::new(0.0, 0.0);
    
    for p in 0..M {
        for q in 0..M {
            let mut phase = 0.0;
            phase += secondQuantizationOneBodyOperator(
                2 * p,
                2 * q,
                &state1,
                &state2,
            );
            phase += secondQuantizationOneBodyOperator(
                2 * p + 1,
                2 * q + 1,
                &state1,
                &state2,
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
                            &state1,
                            &state2,
                        );
                        phase += secondQuantizationTwoBodyOperator(
                            2 * p + 1,
                            2 * q + 1,
                            2 * r + 1,
                            2 * s + 1,
                            &state1,
                            &state2,
                        );
                        phase += secondQuantizationTwoBodyOperator(
                            2 * p,
                            2 * q + 1,
                            2 * r,
                            2 * s + 1,
                            &state1,
                            &state2,
                        );
                        phase += secondQuantizationTwoBodyOperator(
                            2 * p + 1,
                            2 * q,
                            2 * r + 1,
                            2 * s,
                            &state1,
                            &state2,
                        );
                        hmn += 0.5 * V * phase;
                    }
                }
            }
        }
    }
    hmn
}
pub fn computeHamiltonianMatrix(binstates: Vec<BitArray<u64, U64>>,
    v: Vec<Vec<Vec<Vec<Complex<f64>>>>>,
    h: Vec<Vec<Complex<f64>>>,
    M: usize,) -> Mutex<Vec<Vec<Complex<f64>>>> {
        
    let nslater = binstates.len();
    let hamiltonian = Mutex::new(vec![vec![Complex::new(0.0, 0.0); nslater]; nslater]);
    let hspin = transform::honetoSpin(M, &h);
    let vspin = transform::VeetoSpin(M, &v);
    (0..nslater).into_par_iter().for_each(|m| {
        for n in m..nslater {
            let mut hmn = Complex::new(0.0, 0.0);
            let numuniqueorbitals :usize;
            numuniqueorbitals = totalDiff(&binstates[n],&binstates[m]);
            match numuniqueorbitals {
                0 => {
                    hmn = calcMatrixElementIdentialDet(&binstates[n],&binstates[m],&vspin,&hspin,M);
                }
                1 => {
                    hmn = calcMatrixElement(&binstates[n],&binstates[m],&v,&h,M);
                }
                2 => {
                    hmn = calcMatrixElement(&binstates[n],&binstates[m],&v,&h,M);
                }
                _ => {
                    hmn = Complex::new(0.0, 0.0);
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

fn getUniqueOrbitalsInMixIndexListsPlusSign() {

}
fn getOrbitalMixedIndexList(M:usize,state: &BitArray<u64, U64>)->Vec<usize> {
    let mut list = Vec::new();
    for i in num_iter::range(0, M) {
        if state[i] { list.push(i)}
    }
    list
}

fn getCommonOrbitalsInMixedSpinIndexList() {

}
fn totalDiff(state: &BitArray<u64, U64>,
    state2: &BitArray<u64, U64>) -> usize {
    let mut diff:usize = 0;
    for i in num_iter::range(0,state.len()) {
        if state[i] ^ state2[i] { diff +=1}
    }
    diff/2
}

fn calcMatrixElementIdentialDet(state1: &BitArray<u64, U64>,state2: &BitArray<u64, U64>,v: &Vec<Vec<Vec<Vec<Complex<f64>>>>>,
    h: &Vec<Vec<Complex<f64>>>,
    M: usize) -> Complex<f64> {
    let mut helem = Complex::new(0.0, 0.0);
    let mut relem = Complex::new(0.0, 0.0);
    let list = getOrbitalMixedIndexList(M,&state1);
    println!("{}",list.len());
    for m in &list {
        helem += h[*m][*m]
    }
    for m in num_iter::range(0, list.len()-2) {
        for n in num_iter::range(m+1, list.len()-1) {
            relem += v[list[m]][list[n]][list[m]][list[n]]
        }
    }
    helem + relem
}
// fn calcMatrixElementDiffIn1(state1: &BitArray<u64, U64>,
//     state2: &BitArray<u64, U64>) -> Complex<f64> {
//     let hmn = Complex::new(0.0, 0.0);
//     hmn
// }
// fn calcMatrixElementDiffIn2(state1: &BitArray<u64, U64>,
//     state2: &BitArray<u64, U64>) -> Complex<f64> {
//     let hmn = Complex::new(0.0, 0.0);
//     hmn
// }

// pub fn computeHamiltonianMatrix(binstates: Vec<BitArray<u64, U64>>,
//     vspin: Vec<Vec<Vec<Vec<Complex<f64>>>>>,
//     hspin: Vec<Vec<Complex<f64>>>,
//     M: usize,) -> Mutex<Vec<Vec<Complex<f64>>>> {
        
//     let nslater = binstates.len();
//     let hamiltonian = Mutex::new(vec![vec![Complex::new(0.0, 0.0); nslater]; nslater]);
//     (0..nslater).into_par_iter().for_each(|m| {
//         for n in m..nslater {
//             let mut hmn = Complex::new(0.0, 0.0);
//             let numuniqueorbitals :usize;
//             numuniqueorbitals = totalDiff(&binstates[n],&binstates[m]);
//             match numuniqueorbitals {
//                 0 => {
//                     hmn = calcMatrixElementIdentialDet(&binstates[n]);
//                 }
//                 1 => {
//                     hmn = calcMatrixElementDiffIn1(&binstates[n],&binstates[m]);
//                 }
//                 2 => {
//                     hmn = calcMatrixElementDiffIn2(&binstates[n],&binstates[m]);
//                 }
//                 _ => {
//                     hmn = Complex::new(0.0, 0.0);
//                 }
//             }
            
//             let mut hamiltonian = hamiltonian.lock().unwrap();
//             hamiltonian[m][n] = hmn;
//             if m != n {
//                 hamiltonian[n][m] = hamiltonian[m][n];
//             }
//         }
//     });
//     hamiltonian
// }