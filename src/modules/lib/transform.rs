// MO to Spin MO Integrals
// https://github.com/erikkjellgren/SlowQuant/blob/master/slowquant/integraltransformation/IntegralTransform.py
use num::complex::Complex;
use rayon::prelude::*;
use std::thread;
fn honetospin(M: usize, h: &Vec<Vec<Complex<f64>>>) -> Vec<Vec<Complex<f64>>> {
    let mut spin = vec![vec![Complex::new(0.0, 0.0); h.len() * 2]; h.len() * 2];
    for p in num_iter::range(1, h.len() * 2 + 1) {
        for q in num_iter::range(1, h.len() * 2 + 1) {
            if p % 2 == q % 2 {
                let i = ((p as f64 + 1.0) / 2.0).floor() - 1.0;
                let j = ((q as f64 + 1.0) / 2.0).floor() - 1.0;
                spin[p - 1][q - 1] = h[i as usize][j as usize];
            }
        }
    }

    spin
}

fn veetospin(M: usize, v: &Vec<Vec<Vec<Vec<Complex<f64>>>>>) -> Vec<Vec<Vec<Vec<Complex<f64>>>>> {
    let mut spin =
        vec![
            vec![vec![vec![Complex::new(0.0, 0.0); v.len() * 2]; v.len() * 2]; v.len() * 2];
            v.len() * 2
        ];
    for p in num_iter::range(1, v.len() * 2 + 1) {
        for r in num_iter::range(1, v.len() * 2 + 1) {
            for q in num_iter::range(1, v.len() * 2 + 1) {
                for s in num_iter::range(1, v.len() * 2 + 1) {
                    let i = ((p as f64 + 1.0) / 2.0).floor() - 1.0;
                    let j = ((r as f64 + 1.0) / 2.0).floor() - 1.0;
                    let k = ((q as f64 + 1.0) / 2.0).floor() - 1.0;
                    let l = ((s as f64 + 1.0) / 2.0).floor() - 1.0;
                    let bool1: f64;
                    let bool2: f64;
                    if (p % 2 == r % 2) & (q % 2 == s % 2) {
                        bool1 = 1.0
                    } else {
                        bool1 = 0.0
                    }
                    if (p % 2 == s % 2) & (q % 2 == r % 2) {
                        bool2 = 1.0
                    } else {
                        bool2 = 0.0
                    }
                    let value1 = v[i as usize][j as usize][k as usize][l as usize] * bool1;
                    let value2 = v[i as usize][l as usize][k as usize][j as usize] * bool2;
                    spin[p - 1][q - 1][r - 1][s - 1] = value1 - value2;
                }
            }
        }
    }

    spin
}

pub fn tospin(
    M: usize,
    h: &Vec<Vec<Complex<f64>>>,
    v: &Vec<Vec<Vec<Vec<Complex<f64>>>>>,
) -> (Vec<Vec<Complex<f64>>>, Vec<Vec<Vec<Vec<Complex<f64>>>>>) {
    // let  hspin = honetospin(M, h);
    let vspin = veetospin(M, v);
    // (hspin, vspin)

    let mut hspin: Option<Vec<Vec<Complex<f64>>>> = None;
    // let mut vspin: Option<Vec<Vec<Vec<Vec<Complex<f64>>>>>> = None;
    rayon::scope(|s| {
        s.spawn(|_| hspin = Some(honetospin(M, h)));
        // s.spawn(|_| vspin = Some(veetospin(M, v)));
    });
    (hspin.unwrap(), vspin)
    // let hspin:Vec<Vec<Complex<f64>>> = ;
    // let handle = thread::spawn(move || {
    //     hspinn = honetospin(M, h);
    // });
    // let  vspin = veetospin(M, v);
    // handle.join().unwrap();
    // (hspin, vspin)
}
