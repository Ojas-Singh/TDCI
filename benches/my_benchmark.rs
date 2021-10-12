use bit_array::BitArray;
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use rand::Rng;
use typenum::U128;
use typenum::U512;
use typenum::U64;

fn XOR_BITarray64(state: &BitArray<u64, U64>, state2: &BitArray<u64, U64>, n: &u64) -> usize {
    let mut state11 = state.clone();
    let mut state22 = state2.clone();
    state11.difference(&state2);
    state22.difference(state);
    state11.union(&state22);
    state11.into_iter().take(*n as usize).filter(|x| *x).count()
}
fn XOR_BITarray128(state: &BitArray<u64, U128>, state2: &BitArray<u64, U128>, n: &u64) -> usize {
    let mut state11 = state.clone();
    let mut state22 = state2.clone();
    state11.difference(&state2);
    state22.difference(state);
    state11.union(&state22);
    state11.into_iter().take(*n as usize).filter(|x| *x).count()
}
fn XOR_BITarray512(state: &BitArray<u64, U512>, state2: &BitArray<u64, U512>, n: &u64) -> usize {
    let mut state11 = state.clone();
    let mut state22 = state2.clone();
    state11.difference(&state2);
    state22.difference(state);
    state11.union(&state22);
    state11.into_iter().take(*n as usize).filter(|x| *x).count()
}

fn XOR_BOOLvec(state: &Vec<bool>, state2: &Vec<bool>, n: &u64) -> usize {
    let mut diff: usize = 0;
    for i in 0..(*n as usize) {
        if state[i] ^ state2[i] {
            diff += 1
        }
    }
    diff
}

fn XOR_BOOLarray(state: &[bool], state2: &[bool], n: &u64) -> usize {
    let mut diff: usize = 0;
    for i in 0..(*n as usize) {
        if state[i] ^ state2[i] {
            diff += 1
        }
    }
    // for (d, s) in state.iter().take(*n as usize).zip(state2.iter()) {
    //     if d^s {diff +=1};
    // }
    diff
}

fn bench_fibs(c: &mut Criterion) {
    let mut group = c.benchmark_group("XOR BIT vs BOOL");
    for i in [8u64, 16u64, 24u64, 32u64, 40u64, 48u64, 56u64, 64u64].iter() {
        let mut rng = rand::thread_rng();

        let mut bool1 = Vec::new();
        let mut bool2 = Vec::new();
        let mut boola1 = [false; 64];
        let mut boola2 = [false; 64];
        let mut bit1 = BitArray::<u64, U64>::from_elem(false);
        let mut bit2 = BitArray::<u64, U64>::from_elem(false);
        for j in 0..(*i as isize) {
            let a: bool = rng.gen();
            let b: bool = rng.gen();
            bool1.push(a);
            bool2.push(b);
            if a {
                boola1[j as usize] = true;
            }
            if b {
                boola2[j as usize] = true;
            }
            bit1.set(j as usize, a);
            bit2.set(j as usize, b);
        }
        group.bench_with_input(BenchmarkId::new("BitArray", i), i, |b, i| {
            b.iter(|| XOR_BITarray64(black_box(&bit1), black_box(&bit2), i))
        });
        group.bench_with_input(BenchmarkId::new("BooleanVec", i), i, |b, i| {
            b.iter(|| XOR_BOOLvec(black_box(&bool1), black_box(&bool2), i))
        });
        group.bench_with_input(BenchmarkId::new("BooleanArray", i), i, |b, i| {
            b.iter(|| XOR_BOOLarray(black_box(&boola1), black_box(&boola2), i))
        });
    }
    for i in [72u64, 80u64, 88u64, 96u64, 104u64, 112u64, 120u64, 128u64].iter() {
        let mut rng = rand::thread_rng();

        let mut bool1 = Vec::new();
        let mut bool2 = Vec::new();
        let mut boola1 = [false; 128];
        let mut boola2 = [false; 128];
        let mut bit1 = BitArray::<u64, U128>::from_elem(false);
        let mut bit2 = BitArray::<u64, U128>::from_elem(false);
        for j in 0..(*i as isize) {
            let a: bool = rng.gen();
            let b: bool = rng.gen();
            bool1.push(a);
            bool2.push(b);
            if a {
                boola1[j as usize] = true;
            }
            if b {
                boola2[j as usize] = true;
            }
            bit1.set(j as usize, a);
            bit2.set(j as usize, b);
        }
        group.bench_with_input(BenchmarkId::new("BitArray", i), i, |b, i| {
            b.iter(|| XOR_BITarray128(black_box(&bit1), black_box(&bit2), i))
        });
        group.bench_with_input(BenchmarkId::new("BooleanVec", i), i, |b, i| {
            b.iter(|| XOR_BOOLvec(black_box(&bool1), black_box(&bool2), i))
        });
        group.bench_with_input(BenchmarkId::new("BooleanArray", i), i, |b, i| {
            b.iter(|| XOR_BOOLarray(black_box(&boola1), black_box(&boola2), i))
        });
    }
    for i in [
        136u64, 144u64, 152u64, 160u64, 168u64, 176u64, 184u64, 192u64, 200u64, 208u64, 216u64,
        232u64, 240u64, 248u64, 256u64,
    ]
    .iter()
    {
        let mut rng = rand::thread_rng();

        let mut bool1 = Vec::new();
        let mut bool2 = Vec::new();
        let mut boola1 = [false; 512];
        let mut boola2 = [false; 512];
        let mut bit1 = BitArray::<u64, U512>::from_elem(false);
        let mut bit2 = BitArray::<u64, U512>::from_elem(false);
        for j in 0..(*i as isize) {
            let a: bool = rng.gen();
            let b: bool = rng.gen();
            bool1.push(a);
            bool2.push(b);
            if a {
                boola1[j as usize] = true;
            }
            if b {
                boola2[j as usize] = true;
            }
            bit1.set(j as usize, a);
            bit2.set(j as usize, b);
        }
        group.bench_with_input(BenchmarkId::new("BitArray", i), i, |b, i| {
            b.iter(|| XOR_BITarray512(black_box(&bit1), black_box(&bit2), i))
        });
        group.bench_with_input(BenchmarkId::new("BooleanVec", i), i, |b, i| {
            b.iter(|| XOR_BOOLvec(black_box(&bool1), black_box(&bool2), i))
        });
        group.bench_with_input(BenchmarkId::new("BooleanArray", i), i, |b, i| {
            b.iter(|| XOR_BOOLarray(black_box(&boola1), black_box(&boola2), i))
        });
    }
    group.finish();
}

criterion_group!(benches, bench_fibs);
criterion_main!(benches);
