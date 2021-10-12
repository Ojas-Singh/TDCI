#![allow(non_snake_case)]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;
use colored::*;
use ndarray::*;
use ndarray_linalg::*;

use std::env;

#[path = "modules/h_build.rs"]
mod h_build;
#[path = "modules/psi4.rs"]
mod psi4;
#[path = "modules/lib/read_write.rs"]
mod read_write;


fn main() {
    println!("{}", r#" _______ _____   _____ _____ "#.green().bold());
    println!("{}", r#"|__   __|  __ \ / ____|_   _|"#.green().bold());
    println!("{}", r#"   | |  | |  | | |      | |  "#.green().bold());
    println!("{}", r#"   | |  | |  | | |      | |  "#.green().bold());
    println!("{}", r#"   | |  | |__| | |____ _| |_ "#.green().bold());
    println!("{}", r#"   |_|  |_____/ \_____|_____|"#.green().bold());

    println!("{}", "Configuration Interaction Method".blue().bold());
    println!("{}", "https://github.com/Ojas-Singh/TDCI".blue().italic());
    let args: Vec<String> = env::args().collect();

    match args.len() {
        1 => {
            println!("Pass args");
        }
        2 => {
            let result = psi4::gethamiltonian();
            match result {
                Some(x) => {
                    println!("{:?}", x.len());
                    let mut data = Vec::new();

                    let ncols = x.first().map_or(0, |row| row.len());
                    let mut nrows = 0;

                    for i in 0..x.len() {
                        data.extend_from_slice(&x[i]);
                        nrows += 1;
                    }

                    // let arr = Array2::from_shape_vec((nrows, ncols), data).unwrap();
                    // println!("eigenvalues starting ...");
                    // let e = arr.eigvalsh(UPLO::Lower).unwrap();
                    // println!("eigenvalues = \n{:?}", e);

                    println!("Writing to file ...");
                    read_write::save_hamiltonian_txt_real(x, "ham.txt".to_string());
                }
                None => println!("Size too large!"),
            }
        }
        7 => {
            let setting: h_build::Config = arg2cfg(args);
            let result = h_build::h_builder(setting);
            match result {
                Some(x) => {
                    let mut data = Vec::new();

                    let ncols = x.first().map_or(0, |row| row.len());
                    let mut nrows = 0;

                    for i in 0..x.len() {
                        data.extend_from_slice(&x[i]);
                        nrows += 1;
                    }

                    let arr = Array2::from_shape_vec((nrows, ncols), data).unwrap();
                    println!("eigenvalues starting ...");
                    let e = arr.eigvalsh(UPLO::Lower).unwrap();
                    println!("eigenvalues = \n{:?}", e);

                    println!("Writing to file ...");
                    read_write::save_hamiltonian_txt_real(x, "ham.txt".to_string());
                }

                None => println!("Size too large!"),
            }
        }
        _ => {
            help();
        }
    }
}

pub fn arg2cfg(args: Vec<String>) -> h_build::Config {
    let n0 = &args[1];
    let m0 = &args[2];
    let excite = &args[3];
    let f1 = &args[4];
    let f2 = &args[5];
    let truncation = &args[6];
    let nf = n0.trim().parse().unwrap();
    let mf = m0.trim().parse().unwrap();
    let t = truncation.trim().parse().unwrap();
    h_build::Config {
        n: nf,
        m: mf,
        excitation: excite.to_string(),
        oneelectronfilename: f1.to_string(),
        twoelectronfilename: f2.to_string(),
        truncation: t,
    }
}

fn help() {
    println!(
        "usage:: pass system args as :
n , m , excite[Singlet, Triplet] , oneElectronFilename , twoElectronFilename, Truncation[0:FCI, 1:CIS, 2:CISD, i:CI upto ith exitaction] "
    );
}
