use num::complex::Complex;
use std::fs::File;
use std::io::prelude::*;
use std::io::{self, BufRead, BufWriter, Write};
use std::path::Path;
use std::sync::Mutex;

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where
    P: AsRef<Path>,
{
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

pub fn Hone(f1: String, M: usize) -> Vec<Vec<Complex<f64>>> {
    let mut Hvec = vec![vec![Complex::new(0.0, 0.0); M]; M];
    let mut line = 0;
    let mut Hstring = vec![];
    let mut Hab = Vec::new();
    if let Ok(lines) = read_lines(f1) {
        for line in lines {
            if let Ok(ip) = line {
                Hstring.push(ip);
            }
        }
    }
    for i in Hstring {
        let k = i.replace("(", "").replace(")", "").replace(" ", "");
        let vec: Vec<&str> = k.split(',').collect();
        let a: f64 = vec[0].parse().unwrap();
        let b: f64 = vec[1].parse().unwrap();
        let x = Complex::new(a, b);
        Hab.push(x);
    }
    for i in 0..M {
        for j in 0..M {
            Hvec[i][j] = Hab[line];
            line += 1;
        }
    }
    Hvec
}

pub fn Vpqrs(f2: String, M: usize) -> Vec<Vec<Vec<Vec<Complex<f64>>>>> {
    let mut Vvec = vec![vec![vec![vec![Complex::new(0.0, 0.0); M]; M]; M]; M];
    let mut line = 0;
    let mut Vstring = vec![];
    let mut Vab = Vec::new();
    if let Ok(lines) = read_lines(f2) {
        for line in lines {
            if let Ok(ip) = line {
                Vstring.push(ip);
            }
        }
    }
    for i in Vstring {
        let k = i.replace("(", "").replace(")", "").replace(" ", "");
        let vec: Vec<&str> = k.split(',').collect();
        let a: f64 = vec[0].parse().unwrap();
        let b: f64 = vec[1].parse().unwrap();
        let x = Complex::new(a, b);
        Vab.push(x);
    }
    for i in 0..M {
        for j in 0..M {
            for k in 0..M {
                for l in 0..M {
                    Vvec[i][j][k][l] = Vab[line];
                    line += 1;
                }
            }
        }
    }
    Vvec
}

pub fn save_hamiltonian_txt(hamiltonian: Vec<Vec<Complex<f64>>>, file: String) {
    // let hamiltonian = hamiltonian.lock().unwrap();
    let buffer = File::create(file).expect("Unable to create file");
    let mut f = BufWriter::new(buffer);
    for i in 0..hamiltonian.len() {
        for j in 0..hamiltonian[i].len() {
            let data = hamiltonian[i][j].to_string() + " \n";
            f.write_all(data.as_bytes()).expect("Unable to write data");
        }
    }
}

pub fn save_hamiltonian_txt_real(hamiltonian: Vec<Vec<Complex<f64>>>, file: String) {
    // let hamiltonian = hamiltonian.lock().unwrap();
    let buffer = File::create(file).expect("Unable to create file");
    let mut f = BufWriter::new(buffer);
    for i in 0..hamiltonian.len() {
        for j in 0..hamiltonian[i].len() {
            let data = (hamiltonian[i][j].re).to_string() + " \n";
            f.write_all(data.as_bytes()).expect("Unable to write data");
        }
    }
}
pub fn pycode() -> String {
    let mut content = String::new();
    let filename = "code.py";
    match File::open(filename) {
        Ok(mut file) => {
            file.read_to_string(&mut content).unwrap();
        }
        Err(error) => {
            println!("Error opening file {}: {}", filename, error);
        }
    }
    content
}
