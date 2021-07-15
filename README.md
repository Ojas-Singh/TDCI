# Time Dependent Configuration Interaction (TDCI)


##### Building Procedure 
```
  curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
  rustup default nightly
  export RUSTFLAGS="-C target-cpu=native"
  cargo build --release
```
