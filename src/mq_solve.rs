struct MQSystem;

impl MQSystem {
    fn fix_variables(&self, vars: &[(usize, bool)]) -> MQSystem {}
}

struct F2Poly

fn find_solution(system: MQSystem) -> Vec<bool> {
    let found_variables = Vec::new();
    for i in 0..system.n {
        let mut vars = self.found_variables.clone();
        vars.push((i, false));
        if has_solution(system.fix_vars(vars)) {
            found_variables.push((i, false));
        }
        else {
            found_variables.push((i, true));
        }
    }
    found_variables
}

fn has_solution(system: MQSystem) -> bool {
    for k in 0..=system.n {
        let system2 = system.add_linear_eqs(k);
        if parity_count(system2) {
            return true
        }
    }
    false
}

fn parity_count(system: MQSystem) -> bool {
    let n1 = (k0*n).floor();
    let vs = multi_parity_count(system, n1, system.n - n1);
    let parity = vs.fold(false, |acc, x| acc ^ x);
    parity
}

fn multi_parity_count(system: MQSystem, n1: usize, w: usize) -> Vec<bool> {
    let n2 = (n1 - lambda*n).floor();
    let l = n2 + 2;
    let t = 48*n + 1;
    if n2 <= 0 {
        return bruteforce_multi_parity(system, n1, w)
    }
    let mut sb = vec![0i32; ncr(n - n1, w)*(1 << (n1 - n2))];
    
    for k in 1..=t {
        let system2 = system.new_random();
        let v1 = multi_parity_count(system2, n2, 2*l - n2);
        // interpolate G_k(y, u): Möbius transform
        let g = fmt(v1, n - n1);
        
    }

    let mut vs = vec![false; ncr(n - n2, w)];
}


// Fast möbius transform
// Implementation from https://codeforces.com/blog/entry/72488
// Potential new reference: https://bitbucket.org/fes/fes/src/master/src/moebius_transform.c
fn fmt(f: &[bool], n: usize) -> Vec<bool> {
    let mut g = f.to_owned();
    for i in 0..n {
        for mask in 0..(1 << n) {
            if (mask & (1 << i)) != 0 {
                g[mask] ^= g[mask ^ (1 << i)]
            }
        }
    }
    g
}
