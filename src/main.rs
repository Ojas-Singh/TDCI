#![allow(non_snake_case)]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use std::env;
use std::time::Instant;

#[path = "modules/configurations.rs"]
mod configurations;
#[path = "modules/read_write.rs"]
mod read_write;
#[path = "modules/second_quantization.rs"]
mod second_quantization;

pub struct Config {
    n: usize,
    m: usize,
    excitation: String,
    oneelectronfilename: String,
    twoelectronfilename: String,
    truncation: usize,
}

fn main() {
    println!("Configuration Interaction Method");
    let args: Vec<String> = env::args().collect();

    match args.len() {
        1 => {
            println!("Pass args");
        }
        7 => {
            println!("Reading Files ...");
            let setting: Config = arg2cfg(args);
            println!(
                "n : {}, m : {}, excite : {}, oneElectron : {}, twoElectron : {}, truncation: {}",
                setting.n,
                setting.m,
                setting.excitation,
                setting.oneelectronfilename,
                setting.twoelectronfilename,
                setting.truncation
            );
            let binstates = configurations::bit_slaterdeterminants(
                setting.excitation,
                setting.n,
                setting.m,
                setting.truncation,
            );
            println!("Total Generated States :{}", binstates.len());
            let Honemat = read_write::Hone(setting.oneelectronfilename, setting.m);
            let Vmat = read_write::Vpqrs(setting.twoelectronfilename, setting.m);
            let start = Instant::now();
            let ham =
                second_quantization::computeHamiltonianMatrix(binstates, Vmat, Honemat, setting.m);
            let duration = start.elapsed();
            println!(
                "Time elapsed in computeHamiltonianMatrix is: {:?}",
                duration
            );
            read_write::save_hamiltonian_txt(ham, "ham.txt".to_string());
        }

        _ => {
            help();
        }
    }
}

pub fn arg2cfg(args: Vec<String>) -> Config {
    let n0 = &args[1];
    let m0 = &args[2];
    let excite = &args[3];
    let f1 = &args[4];
    let f2 = &args[5];
    let truncation = &args[6];
    let nf = n0.trim().parse().unwrap();
    let mf = m0.trim().parse().unwrap();
    let t = truncation.trim().parse().unwrap();
    Config {
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
n , m , excite , oneElectronFilename , twoElectronFilename "
    );
}
