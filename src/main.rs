#![allow(non_snake_case)]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use colored::*;
use std::env;
use std::process;
use std::time::Instant;
use sysinfo::SystemExt;

#[path = "modules/beta.rs"]
mod beta;
#[path = "modules/configurations.rs"]
mod configurations;
#[path = "modules/fast.rs"]
mod fast;
#[path = "modules/read_write.rs"]
mod read_write;
#[path = "modules/second_quantization.rs"]
mod second_quantization;
#[path = "modules/transform.rs"]
mod transform;
pub struct Config {
    n: usize,
    m: usize,
    excitation: String,
    oneelectronfilename: String,
    twoelectronfilename: String,
    truncation: usize,
}

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

    let mut system = sysinfo::System::new();
    system.refresh_all();
    match args.len() {
        1 => {
            println!("Pass args");
        }
        7 => {
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
            println!("Reading Files ...");
            let start0 = Instant::now();
            let Honemat = read_write::Hone(setting.oneelectronfilename, setting.m);
            let Vmat = read_write::Vpqrs(setting.twoelectronfilename, setting.m);
            let duration0 = start0.elapsed();
            println!(
                "**[Timing] Time elapsed in Reading files is: {:?}",
                duration0
            );
            println!("Generating States ...");
            let start1 = Instant::now();
            let binstates = configurations::bit_slaterdeterminants(
                setting.excitation,
                setting.n,
                setting.m,
                setting.truncation,
            );
            println!("Total Generated States :{}", binstates.len());
            let duration1 = start1.elapsed();
            println!(
                "**[Timing] Time elapsed in Generating States is: {:?}",
                duration1
            );
            let memory = ((binstates.len() as isize).pow(2) * 16 / 1000000000) as u64;
            println!(
                "Estimated Memory utilization : {} GB and System has {} GB",
                memory,
                system.total_memory() / 1000000
            );
            if memory > (system.total_memory() / 1000000) {
                println!("{}", "insufficient memory !!!".red().bold());
                process::abort();
            }
            let start = Instant::now();
            let ham = beta::computeHamiltonianMatrix(binstates, Vmat, Honemat, setting.m);
            let duration = start.elapsed();
            println!(
                "**[Timing] Time elapsed in computeHamiltonianMatrix is: {:?}",
                duration
            );
            println!("Writing to file ...");
            let start2 = Instant::now();
            read_write::save_hamiltonian_txt(ham, "ham.txt".to_string());
            let duration2 = start2.elapsed();
            println!("**[Timing] Time elapsed in Writing  is: {:?}", duration2);
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
