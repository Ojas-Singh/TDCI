#![allow(non_snake_case)]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use std::env;
extern crate num_iter;
extern crate typenum;
use bit_array::BitArray;
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufRead, BufWriter, Write};
use std::path::Path;
use std::sync::Mutex;
use typenum::U64;
use std::time::{ Instant};


fn help() {
    println!(
        "usage:: pass system args as :
n , m , excite , oneElectronFilename , twoElectronFilename "
    );
}

fn main() {
    // rayon::ThreadPoolBuilder::new().num_threads(12).build_global().unwrap();
    println!("Reading Files ...");
    let args: Vec<String> = env::args().collect();
    let n0 = &args[1];
    let m0 = &args[2];
    let excite = &args[3];
    let f1 = &args[4];
    let f2 = &args[5];
    let n: isize = match n0.parse() {
        Ok(n0) => n0,
        Err(_) => {
            eprintln!("error: second argument not an integer");
            help();
            return;
        }
    };
    let m: isize = match m0.parse() {
        Ok(m0) => m0,
        Err(_) => {
            eprintln!("error: second argument not an integer");
            help();
            return;
        }
    };
    println!(
        "n : {}, m : {}, excite : {}, oneElectron : {}, twoElectron : {}",
        n, m, excite, f1, f2
    );
    let Honemat = Hone(f1.to_string(), m as usize);
    let Vmat = Vpqrs(f2.to_string(), m as usize);
    let binstates = createslaterdeterminants_SD(n as usize, m as usize, excite.to_string());
    println!("Total Generated States :{}", binstates.len());
    let start = Instant::now();
    let ham = computeHamiltonianMatrix(binstates, Vmat, Honemat, m as usize);
    let duration = start.elapsed();
    println!("Time elapsed in expensive_function() is: {:?}", duration);

    // print!("Done! {:?}",ham);
    save_hamiltonian_txt(ham, "ham.txt".to_string());
}

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where
    P: AsRef<Path>,
{
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

pub fn Hone(f1: String, M: usize) -> Vec<Vec<f64>> {
    let mut Hvec = vec![vec![0.0; M]; M];
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
        let vec: Vec<&str> = k.split(",").collect();
        let a: f64 = vec[0].parse().unwrap();
        // let b: f64 = vec[1].parse().unwrap();
        Hab.push([a]);
    }
    for i in 0..M {
        for j in 0..M {
            Hvec[i][j] = Hab[line][0];
            line += 1;
        }
    }
    return Hvec;
}

pub fn Vpqrs(f2: String, M: usize) -> Vec<Vec<Vec<Vec<f64>>>> {
    let mut Vvec = vec![vec![vec![vec![0.0; M]; M]; M]; M];
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
        let vec: Vec<&str> = k.split(",").collect();
        let a: f64 = vec[0].parse().unwrap();
        // let b: f64 = vec[1].parse().unwrap();
        Vab.push([a]);
    }
    for i in 0..M {
        for j in 0..M {
            for k in 0..M {
                for l in 0..M {
                    Vvec[i][j][k][l] = Vab[line][0];
                    line += 1;
                }
            }
        }
    }
    return Vvec;
}

pub fn creatinitialstate(excite: String, n: usize) -> Vec<Vec<isize>> {
    let mut stateout = Vec::new();
    if excite == "Singlet" {
        if n % 2 == 0 {
            let mut stateup = Vec::new();
            let mut statedown = Vec::new();
            for i in num_iter::range(0, n / 2) {
                stateup.push((i + 1) as isize);
                statedown.push((i + 1) as isize);
            }
            stateout.push(stateup);
            stateout.push(statedown);
            return stateout;
        } else {
            let mut stateup = Vec::new();
            let mut statedown = Vec::new();
            for i in num_iter::range(0, n / 2) {
                stateup.push((i + 1) as isize);
                statedown.push((i + 1) as isize);
            }
            stateup.push((n / 2 + 1) as isize);
            stateout.push(stateup);
            stateout.push(statedown);
            return stateout;
        }
    }
    if excite == "Triplet" {
        //Left for Future!
        if n % 2 == 0 {
            let mut stateup = Vec::new();
            let mut statedown = Vec::new();
            for i in num_iter::range(0, n / 2) {
                stateup.push((i + 1) as isize);
                statedown.push((i + 1) as isize);
            }
            stateout.push(stateup);
            stateout.push(statedown);
            return stateout;
        } else {
            let mut stateup = Vec::new();
            let mut statedown = Vec::new();
            for i in num_iter::range(0, n / 2) {
                stateup.push((i + 1) as isize);
                statedown.push((i + 1) as isize);
            }
            stateup.push((n / 2 + 1) as isize);
            stateout.push(stateup);
            stateout.push(statedown);
            return stateout;
        }
    }
    return stateout;
}

pub fn odometer(state: Vec<isize>, n: isize, m: isize) -> Vec<isize> {
    let mut newstate = state;
    for j in num_iter::range_step(n - 1, -1, -1) {
        if newstate[j as usize] < m + 1 - n + j {
            let l = newstate[j as usize];
            for k in num_iter::range(j, n) {
                newstate[k as usize] = l + 1 + k - j;
            }
            if newstate[j as usize] != l {
                return newstate;
            }
        }
    }
    newstate.iter_mut().for_each(|x| *x = 0);
    return newstate;
}

pub fn compare(state: Vec<isize>, ground: Vec<isize>) -> usize {
    let mut numberofexited = 0;
    for i in &state {
        if !ground.contains(&i) {
            numberofexited += 1;
        }
    }
    return numberofexited;
}

pub fn createbinarystatearray(state: Vec<isize>) -> BitArray<u64, U64> {
    let mut binstate = BitArray::<u64, U64>::from_elem(false);
    for i in state {
        let k: usize = (i - 1) as usize;
        binstate.set(k, true);
    }
    return binstate;
}

pub fn createslaterdeterminants(n: usize, m: usize, excite: String) -> Vec<BitArray<u64, U64>> {
    let mut binstates = Vec::new();
    let N: usize;
    if n % 2 == 0 {
        N = n / 2;
    } else {
        N = n / 2 + 1;
    }
    let mut stateup = creatinitialstate(excite.to_string(), n as usize)[0].clone();
    let mut statedown = creatinitialstate(excite.to_string(), n as usize)[1].clone();
    let mut statesup = Vec::new();
    statesup.push(stateup.clone());
    let mut statesdown = Vec::new();
    statesdown.push(statedown.clone());
    let mut up = true;
    let mut down = true;
    while up {
        stateup = odometer(stateup, N as isize, m as isize);
        let sm: isize = stateup.iter().sum();
        if sm == 0 {
            up = false;
        } else {
            statesup.push(stateup.clone());
        }
    }
    while down {
        statedown = odometer(statedown, (n / 2) as isize, m as isize);
        let sm: isize = statedown.iter().sum();
        if sm == 0 {
            down = false;
        } else {
            statesdown.push(statedown.clone());
        }
    }
    for i in statesup {
        for j in &statesdown {
            let binstate = createbinarystatearray(mix(i.to_vec(), j.to_vec()));
            binstates.push(binstate);
        }
    }

    return binstates;
}

pub fn createslaterdeterminants_SD(n: usize, m: usize, excite: String) -> Vec<BitArray<u64, U64>> {
    let mut binstates = Vec::new();
    let N: usize;
    if n % 2 == 0 {
        N = n / 2;
    } else {
        N = n / 2 + 1;
    }
    let mut stateup = creatinitialstate(excite.to_string(), n as usize)[0].clone();
    let mut statedown = creatinitialstate(excite.to_string(), n as usize)[1].clone();
    let mut statesup = Vec::new();
    statesup.push(stateup.clone());
    let mut statesdown = Vec::new();
    statesdown.push(statedown.clone());
    let mut up = true;
    let mut down = true;
    let ground = mix(statesup[0].to_vec(), statesdown[0].to_vec());
    while up {
        stateup = odometer(stateup, N as isize, m as isize);
        let sm: isize = stateup.iter().sum();
        if sm == 0 {
            up = false;
        } else {
            if compare(stateup.clone(), ground.clone()) < 3 {
                statesup.push(stateup.clone());
            }
        }
    }
    while down {
        statedown = odometer(statedown, (n / 2) as isize, m as isize);
        let sm: isize = statedown.iter().sum();
        if sm == 0 {
            down = false;
        } else {
            if compare(statedown.clone(), ground.clone()) < 3 {
                statesdown.push(statedown.clone());
            }
        }
    }
    for i in statesup {
        for j in &statesdown {
            let state = mix(i.to_vec(), j.to_vec());
            if compare(state.clone(), ground.clone()) < 3 {
                let binstate = createbinarystatearray(state);
                binstates.push(binstate);
            }
        }
    }
    return binstates;
}

pub fn mix(state1: Vec<isize>, state2: Vec<isize>) -> Vec<isize> {
    let mut state = Vec::new();
    for i in state1 {
        state.push(2 * i - 1);
    }
    for i in state2 {
        state.push(2 * i);
    }
    return state;
}
#[inline(always)]
pub fn sign(n: usize, binstate: &BitArray<u64, U64>) -> f64 {
    let mut s = 1.0;
    for i in binstate.iter().take(n) {
        if i {
            s *= -1.0;
        }
    }
    return s; 
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
    return phase;
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
    return phase;
}

pub fn computeHamiltonianMatrix(
    binstates: Vec<BitArray<u64, U64>>,
    v: Vec<Vec<Vec<Vec<f64>>>>,
    h: Vec<Vec<f64>>,
    M: usize,
) -> Mutex<Vec<Vec<f64>>> {
    let nslater = binstates.len();
    let hamiltonian = Mutex::new(vec![vec![0.0; nslater]; nslater]);
    (0..nslater).into_par_iter().for_each(|m| {
        for n in m..nslater {
            let mut hmn = 0.0;
            for p in 0..M {
                for q in 0..M {
                    let mut phase = 0.0;
                    phase += secondQuantizationOneBodyOperator(
                        2 * p,
                        2 * q,
                        &binstates[n],
                        &binstates[m],
                    );
                    phase += secondQuantizationOneBodyOperator(
                        2 * p + 1,
                        2 * q + 1,
                        &binstates[n],
                        &binstates[m],
                    );
                    hmn += h[p][q] * phase;
                    
                    for r in 0..M {
                        for s in 0..M {
                            let V = v[p][r][q][s];
                            if  V != 0.0 {
                                phase = 0.0;
                                phase += secondQuantizationTwoBodyOperator(
                                    2 * p,
                                    2 * q,
                                    2 * r,
                                    2 * s,
                                    &binstates[n],
                                    &binstates[m],
                                );
                                phase += secondQuantizationTwoBodyOperator(
                                    2 * p + 1,
                                    2 * q + 1,
                                    2 * r + 1,
                                    2 * s + 1,
                                    &binstates[n],
                                    &binstates[m],
                                );
                                phase += secondQuantizationTwoBodyOperator(
                                    2 * p,
                                    2 * q + 1,
                                    2 * r,
                                    2 * s + 1,
                                    &binstates[n],
                                    &binstates[m],
                                );
                                phase += secondQuantizationTwoBodyOperator(
                                    2 * p + 1,
                                    2 * q,
                                    2 * r + 1,
                                    2 * s,
                                    &binstates[n],
                                    &binstates[m],
                                );
                                hmn += 0.5 * V * phase;
                            }
                            
                            
                        }
                    }
                }
            }
            let mut hamiltonian = hamiltonian.lock().unwrap();
            hamiltonian[m][n] = hmn;
            if m != n {
                hamiltonian[n][m] = hamiltonian[m][n];
            }
        }
    });
    return hamiltonian;
}


pub fn save_hamiltonian_txt(hamiltonian: Mutex<Vec<Vec<f64>>>, file: String) {
    let hamiltonian = hamiltonian.lock().unwrap();
    let buffer = File::create(file).expect("Unable to create file");
    let mut f = BufWriter::new(buffer);
    for i in 0..hamiltonian.len() {
        for j in 0..hamiltonian[i].len() {
            let data = hamiltonian[i][j].to_string() + " \n";
            f.write_all(data.as_bytes()).expect("Unable to write data");
        }
    }
}
