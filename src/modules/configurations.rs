use bit_array::BitArray;
use typenum::U64;

pub fn bit_slaterdeterminants(
    excitation: String,
    n: usize,
    m: usize,
    truncation: usize,
) -> Vec<BitArray<u64, U64>> {
    fn creatinitialstate(excite: String, n: usize) -> Vec<Vec<isize>> {
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

    fn odometer(state: Vec<isize>, n: isize, m: isize) -> Vec<isize> {
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

    fn compare(state: Vec<isize>, ground: Vec<isize>) -> usize {
        let mut numberofexited = 0;
        for i in &state {
            if !ground.contains(&i) {
                numberofexited += 1;
            }
        }
        return numberofexited;
    }

    fn createbinarystatearray(state: Vec<isize>) -> BitArray<u64, U64> {
        let mut binstate = BitArray::<u64, U64>::from_elem(false);
        for i in state {
            let k: usize = (i - 1) as usize;
            binstate.set(k, true);
        }
        return binstate;
    }

    fn mix(state1: Vec<isize>, state2: Vec<isize>) -> Vec<isize> {
        let mut state = Vec::new();
        for i in state1 {
            state.push(2 * i - 1);
        }
        for i in state2 {
            state.push(2 * i);
        }
        return state;
    }

    fn createslaterdeterminants(
        n: usize,
        m: usize,
        excite: String,
        t: usize,
    ) -> Vec<BitArray<u64, U64>> {
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
                if compare(stateup.clone(), ground.clone()) < t + 1 {
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
                if compare(statedown.clone(), ground.clone()) < t + 1 {
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
    let excite = excitation.clone();
    let binstates = createslaterdeterminants(n, m, excite, truncation);
    return binstates;
}
