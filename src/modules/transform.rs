// MO to Spin MO Integrals 
// https://github.com/erikkjellgren/SlowQuant/blob/master/slowquant/integraltransformation/IntegralTransform.py

use num::complex::Complex;
pub fn honetoSpin(M: usize,h: &Vec<Vec<Complex<f64>>>) -> Vec<Vec<Complex<f64>>>{
    let mut spin = vec![vec![Complex::new(0.0, 0.0); M*2]; M*2];
    for p in num_iter::range(1, M*2) {
        for q in num_iter::range(1, M*2) {
            if p%2 == q%2 {
                let i = (p as f64 )/2.0  ;
                let j = (q as f64 )/2.0 ;
                spin[p-1][q-1] = h[i.floor() as usize][j.floor() as usize];
            }
        }
    }


    spin
}

pub fn VeetoSpin(M: usize,v :&Vec<Vec<Vec<Vec<Complex<f64>>>>>) -> Vec<Vec<Vec<Vec<Complex<f64>>>>>{
    let mut spin = vec![vec![vec![vec![Complex::new(0.0, 0.0); M*2]; M*2]; M*2]; M*2];
    for p in num_iter::range(1, M*2) {
        for r in num_iter::range(1, M*2) {
            if p%2==r%2 {
                for q in num_iter::range(1, M*2) {
                    for s in num_iter::range(1, M*2) {
                        if q%2 == s%2 {
                            let i = (p as f64 )/2.0 ;
                            let j = (r as f64 )/2.0 ;
                            let k = (q as f64 )/2.0 ;
                            let l = (s as f64 )/2.0 ;
                            spin[p-1][q-1][r-1][s-1] = v[i.floor() as usize][j.floor() as usize][k.floor() as usize][l.floor() as usize]
                        }
                    }
                }
            }
        }
    }
    
    spin
}