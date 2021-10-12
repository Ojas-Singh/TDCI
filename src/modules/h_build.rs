use colored::*;
use num::complex::Complex;
use std::process;
use std::time::Instant;
use sysinfo::SystemExt;

pub struct Config {
    pub n: usize,
    pub m: usize,
    pub excitation: String,
    pub oneelectronfilename: String,
    pub twoelectronfilename: String,
    pub truncation: usize,
}

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
#[path = "lib/read_write.rs"]
mod read_write;
#[path = "lib/transform.rs"]
mod transform;

pub fn h_builder(setting: Config) -> Option<Vec<Vec<Complex<f64>>>> {
    let mut system = sysinfo::System::new();
    system.refresh_all();
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
    let m = &setting.m;
    let (hspin, vspin) = transform::tospin(*m, &Honemat, &Vmat);
    if m * 2 <= 64 {
        println!("Generating States ...");
        let start1 = Instant::now();
        let binstates = configurations::bit_slaterdeterminants(
            setting.excitation.clone(),
            setting.n.clone(),
            setting.m.clone(),
            setting.truncation.clone(),
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
            // process::abort();
        }
        let start = Instant::now();
        let ham = beta::computeHamiltonianMatrix(
            binstates,
            vspin.clone(),
            hspin.clone(),
            setting.m.clone(),
        );

        let duration = start.elapsed();
        println!(
            "**[Timing] Time elapsed in computeHamiltonianMatrix is: {:?}",
            duration
        );
        // let hamiltonian = ham.into_inner().unwrap();
        let hamiltonian = ham;
        return Some(hamiltonian);
        // println!("Writing to file ...");
        // let start2 = Instant::now();
        // read_write::save_hamiltonian_txt(ham, "ham.txt".to_string());
        // let duration2 = start2.elapsed();
        // println!("**[Timing] Time elapsed in Writing  is: {:?}", duration2);
    }
    if 64 < m * 2 && m * 2 <= 128 {
        println!("Generating States ...");
        let start1 = Instant::now();
        let binstates = configurations128::bit_slaterdeterminants(
            setting.excitation.clone(),
            setting.n.clone(),
            setting.m.clone(),
            setting.truncation.clone(),
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
        let ham = beta128::computeHamiltonianMatrix(
            binstates,
            vspin.clone(),
            hspin.clone(),
            setting.m.clone(),
        );
        let duration = start.elapsed();
        println!(
            "**[Timing] Time elapsed in computeHamiltonianMatrix is: {:?}",
            duration
        );
        // println!("Writing to file ...");
        // let start2 = Instant::now();
        // read_write::save_hamiltonian_txt(ham, "ham.txt".to_string());
        // let duration2 = start2.elapsed();
        // println!("**[Timing] Time elapsed in Writing  is: {:?}", duration2);
        // }
        // let hamiltonian = ham.into_inner().unwrap();
        let hamiltonian = ham;
        return Some(hamiltonian);
    }
    if 128 < m * 2 && m * 2 < 512 {
        println!("Generating States ...");
        let start1 = Instant::now();
        let binstates = configurations512::bit_slaterdeterminants(
            setting.excitation.clone(),
            setting.n.clone(),
            setting.m.clone(),
            setting.truncation.clone(),
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
        let ham = beta512::computeHamiltonianMatrix(
            binstates,
            vspin.clone(),
            hspin.clone(),
            setting.m.clone(),
        );
        let duration = start.elapsed();
        println!(
            "**[Timing] Time elapsed in computeHamiltonianMatrix is: {:?}",
            duration
        );
        // println!("Writing to file ...");
        // let start2 = Instant::now();
        // read_write::save_hamiltonian_txt(ham, "ham.txt".to_string());
        // let duration2 = start2.elapsed();
        // println!("**[Timing] Time elapsed in Writing  is: {:?}", duration2);
        // let hamiltonian = ham.into_inner().unwrap();
        let hamiltonian = ham;
        return Some(hamiltonian);
    }
    return None;
}
