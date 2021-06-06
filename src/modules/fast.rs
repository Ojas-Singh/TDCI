// Psi4Numpy COnfiguration Interaction Algorithm Implementation
use bit_array::BitArray;
use num::complex::Complex;
use rayon::prelude::*;
use std::sync::Mutex;
use typenum::U64;

#[path = "second_quantization.rs"]
mod second_quantization;

pub fn calcMatrixElement(
    state1: &BitArray<u64, U64>,
    state2: &BitArray<u64, U64>,
    v: &Vec<Vec<Vec<Vec<Complex<f64>>>>>,
    h: &Vec<Vec<Complex<f64>>>,
    M: usize,
) -> Complex<f64> {
    let mut hmn = Complex::new(0.0, 0.0);

    for p in 0..M {
        for q in 0..M {
            let mut phase = 0.0;
            phase += second_quantization::secondQuantizationOneBodyOperator(
                2 * p,
                2 * q,
                &state1,
                &state2,
            );
            phase += second_quantization::secondQuantizationOneBodyOperator(
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
                        phase += second_quantization::secondQuantizationTwoBodyOperator(
                            2 * p,
                            2 * q,
                            2 * r,
                            2 * s,
                            &state1,
                            &state2,
                        );
                        phase += second_quantization::secondQuantizationTwoBodyOperator(
                            2 * p + 1,
                            2 * q + 1,
                            2 * r + 1,
                            2 * s + 1,
                            &state1,
                            &state2,
                        );
                        phase += second_quantization::secondQuantizationTwoBodyOperator(
                            2 * p,
                            2 * q + 1,
                            2 * r,
                            2 * s + 1,
                            &state1,
                            &state2,
                        );
                        phase += second_quantization::secondQuantizationTwoBodyOperator(
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
            let numuniqueorbitals: usize;
            numuniqueorbitals = totalDiff(&binstates[n], &binstates[m]);
            match numuniqueorbitals {
                0 => {
                    hmn = calcMatrixElement(&binstates[n], &binstates[m], &v, &h, M);
                }
                1 => {
                    hmn = calcMatrixElement(&binstates[n], &binstates[m], &v, &h, M);
                }
                2 => {
                    hmn = calcMatrixElement(&binstates[n], &binstates[m], &v, &h, M);
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

fn totalDiff(state: &BitArray<u64, U64>, state2: &BitArray<u64, U64>) -> usize {
    let mut diff: usize = 0;
    for i in num_iter::range(0, state.len()) {
        if state[i] ^ state2[i] {
            diff += 1
        }
    }
    diff / 2
}
