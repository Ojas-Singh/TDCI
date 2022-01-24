# Configuration Interaction


### Overview
TDCI is a High-performance Configuration Interaction program written in Rust.


Will implement time-propagation in future!

### Getting Started
##### Dependencies
[Psi4](http://psicode.org/psi4manual/1.1/build_obtaining.html) [Download installer](http://vergil.chemistry.gatech.edu/psicode-download/1.1.html) and install according to [instructions](http://psicode.org/psi4manual/1.1/conda.html#how-to-install-a-psi4-binary-with-the-psi4conda-installer-command-line).

##### Building Procedure 
```
  curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
  rustup default nightly
  export RUSTFLAGS="-C target-cpu=native" or export RUSTFLAGS="-C target-feature=+avx2"
  cargo build --release
```

##### To set Threads limit 
```
export RAYON_NUM_THREADS="6"
```

##### Environment Variables
```
export LD_LIBRARY_PATH=$HOME/psi4conda/lib/
export KMP_DUPLICATE_LIB_OK=TRUE
```

### Benchmarks

**Sulphur**  aug-cc-pVDZ  [5s4p2d] → 27 function [Old Benchmark , Current Version is 2 times faster]

```
CPU(s):  48 
Intel(R) Xeon(R) CPU E5-2650 v4 @ 2.20GHz 
[Timing] Time elapsed in Reading files is: 319.593169ms
Total Generated States : 32985
[Timing] Time elapsed in Generating States is: 5.16356ms
Memory usage ~ 18 GB
```
| Cores | *Compute Hamiltonian Matrix* (*s*) |
| ----: | :--------------------------------: |
|     1 |                33.9                |
|     2 |                17.6                |
|     4 |                9.9                 |
|     6 |                6.8                 |
|     8 |                5.2                 |
|    12 |                3.6                 |
|    18 |                2.8                 |
|    24 |                2.5                 |
|    30 |                2.11                |
|    36 |                1.84                |
|    42 |                1.79                |
|    48 |                1.69                |

![](https://raw.githubusercontent.com/Ojas-Singh/TDCI/master/docs/images/speedup.png)

<!-- ![](https://raw.githubusercontent.com/Ojas-Singh/TDCI/master/docs/1.PNG) -->


> Against Psi4Numpy

**Oxygen**  cc-pVDZ  [3s2p1d] → 14 function. **CISD**

```
CPU(s):  6 
AMD Ryzen 5 1600 (6) @ 3.775GHz 
Total Generated States :2221
```

|                               | Our Code | Psi4Numpy |
| ----------------------------- | -------- | --------- |
| Configuration Generation Time | 250 µs   | 4 ms      |
| Matrix Building Time          | 35 ms    | 15 s      |



```
CPU(s):  48 
Intel(R) Xeon(R) CPU E5-2650 v4 @ 2.20GHz 
```
##### H2O

| Basis Set   | Determinants | Transformation Time | Matrix Build Time |
| ----------- | ------------ | ------------------- | ----------------- |
| cc-pVDZ     | 12636        | 103.6 ms            | 371.3 ms          |
| aug-cc-pVDZ | 45361        | 1.153 s             | 9.635 s           |
| cc-pVTZ     | 98316        | 4.747 s             | 48.83 s           |

##### NH3

| Basis Set   | Determinants | Transformation Time | Matrix Build Time |
| ----------- | ------------ | ------------------- | ----------------- |
| cc-pVDZ     | 20161        | 233.8 ms            | 911.93 ms         |
| aug-cc-pVDZ | 70876        | 2.46 s              | 26.33 s           |
| cc-pVTZ     | 157116       | 10.99 s             | 137.85 s          |
