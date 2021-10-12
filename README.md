# Time Dependent Configuration Interaction (TDCI)


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

##### Psi4 API 
```
export LD_LIBRARY_PATH=$HOME/psi4conda/lib/
```

##### Openmpi 
```
export KMP_DUPLICATE_LIB_OK=TRUE
```