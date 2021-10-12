use num::complex::Complex;
use pyo3::{prelude::*, types::PyModule};
#[path = "lib/read_write.rs"]
mod read_write;

#[path = "lib/beta.rs"]
mod beta;
#[path = "lib/beta128.rs"]
mod beta128;
#[path = "lib/beta512.rs"]
mod beta512;
#[path = "lib/configurations.rs"]
mod configurations;
#[path = "lib/configurations128.rs"]
mod configurations128;
#[path = "lib/configurations512.rs"]
mod configurations512;

fn getvariable() -> PyResult<(
    Vec<Vec<Complex<f64>>>,
    Vec<Vec<Vec<Vec<Complex<f64>>>>>,
    usize,
    usize,
)> {
    Python::with_gil(|py| {
        let activators = PyModule::from_code(py, &read_write::pycode(), "code.py", "activators")?;
        let (H, V, M, N): (Vec<Vec<f64>>, Vec<Vec<Vec<Vec<f64>>>>, usize, usize) =
            activators.getattr("run")?.call1((1.0,))?.extract()?;
        let mut Hcom: Vec<Vec<Complex<f64>>> = Vec::new();
        let mut Vcom: Vec<Vec<Vec<Vec<Complex<f64>>>>> = Vec::new();
        for i in H {
            let mut t: Vec<Complex<f64>> = Vec::new();
            for j in i {
                t.push(Complex::new(j, 0.0));
            }
            Hcom.push(t)
        }
        for i in V {
            let mut t2: Vec<Vec<Vec<Complex<f64>>>> = Vec::new();
            for j in i {
                let mut t3: Vec<Vec<Complex<f64>>> = Vec::new();
                for k in j {
                    let mut t4: Vec<Complex<f64>> = Vec::new();
                    for l in k {
                        t4.push(Complex::new(l, 0.0))
                    }
                    t3.push(t4)
                }
                t2.push(t3)
            }
            Vcom.push(t2)
        }

        Ok((Hcom, Vcom, M, N))
    })
}

pub fn gethamiltonian() -> Option<Vec<Vec<Complex<f64>>>> {
    let (H, V, M, N) = getvariable().unwrap();
    println!("{:?} and {:?}", M, N);
    if M * 2 <= 64 {
        let binstates = configurations::bit_slaterdeterminants("Singlet".to_string(), N, M, 2);
        println!("Total Generated States :{}", binstates.len());
        let ham = beta::computeHamiltonianMatrix(binstates, V, H, M);

        return Some(ham);
    }
    if 64 < M * 2 && M * 2 <= 128 {
        let binstates = configurations128::bit_slaterdeterminants("Singlet".to_string(), N, M, 2);
        println!("Total Generated States :{}", binstates.len());
        let ham = beta128::computeHamiltonianMatrix(binstates, V, H, M);

        return Some(ham);
    }
    if 128 < M * 2 && M * 2 < 512 {
        let binstates = configurations512::bit_slaterdeterminants("Singlet".to_string(), N, M, 2);
        println!("Total Generated States :{}", binstates.len());
        let ham = beta512::computeHamiltonianMatrix(binstates, V, H, M);

        return Some(ham);
    }
    return None;
}
