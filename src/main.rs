use std::io::prelude::*;
use std::collections::BTreeMap;
use rand;

mod bitset;
use bitset::Bitset;

#[derive(Debug, Clone)]
struct Poly {
    coefs: Vec<bool>,
    c: bool,
    n: usize,
}

impl Poly {
    fn parse(n: usize, s: &str) -> Poly {
        let mut coefs = vec![false; n*(n+1) / 2];
        let mut c = false;

        for (i, v) in s.split_whitespace().enumerate() {
            if i < n*(n+1) / 2 {
                coefs[i] = v.trim() == "1";
            }
            else if i < n*(n+1) / 2 + n {
                let j = i - n*(n+1) / 2;
                coefs[j*(j+1) / 2 + j] ^= v.trim() == "1";
            }
            else {
                c = v.trim() == "1";
            }
        }

        Poly { coefs, c, n }
    }

    fn zero(n: usize) -> Poly {
        Poly {
            coefs: vec![false; n*(n+1) / 2],
            c: false,
            n
        }
    }

    fn fix_vars(&self, vars: &[(usize, bool)]) -> Poly {
        let mut poly = self.clone();
        poly.n -= vars.len();
        poly.coefs = vec![false; poly.n*(poly.n+1) / 2];
        let old_coefs = self.coefs.clone();
        let mut has_val = vec![false; self.n];
        let mut new_val = vec![false; self.n];
        for (var, val) in vars {
            has_val[*var] = true;
            new_val[*var] = *val;
        }
        let mut var_map = vec![0; self.n];
        {
            let mut count = 0;
            for i in 0..self.n {
                if !has_val[i] {
                    var_map[i] = count;
                    count += 1;
                }
            }
        }

        //dbg!(&has_val, &new_val, &var_map);

        for i in 0..self.n {
            let new_i = var_map[i];
            if has_val[i] {
                for j in 0..=i {
                    let new_j = var_map[j];
                    if has_val[j] {
                        poly.c ^= new_val[i] & new_val[j] 
                            & old_coefs[i*(i+1)/2 + j];
                    }
                    else {
                        poly.coefs[new_j*(new_j+1) / 2 + new_j] ^= new_val[i]
                            & old_coefs[i*(i+1)/2 + j];
                    }
                }
            }
            else {
                for j in 0..=i {
                    let new_j = var_map[j];
                    if has_val[j] {
                        //dbg!((new_i, new_val[j]));
                        poly.coefs[new_i*(new_i+1) / 2 + new_i] ^= new_val[j]
                            & old_coefs[i*(i+1)/2 + j];
                    }
                    else {
                        //dbg!((new_i, new_j, old_coefs[i*(i+1)/2 + j]));
                        poly.coefs[new_i*(new_i+1) / 2 + new_j] ^= 
                            old_coefs[i*(i+1)/2 + j];
                    }
                }
            }
        }
        poly
    }

    fn random_affine(n: usize) -> Poly {
        let mut coefs = vec![false; n*(n+1) / 2];
        for i in 0..n {
            coefs[i*(i + 1) / 2 + i] = rand::random();
        }
        let c = rand::random();
        Poly { coefs, c, n }
    }

    fn evaluate(&self, params: &[bool]) -> bool {
         let mut res = self.c;
         for i in 0..self.n {
            for j in 0..=i {
                res ^= params[i] & params[j] & self.coefs[i*(i+1)/2 + j];
            }
         }
         res
    }
}

impl std::ops::AddAssign for Poly {
    fn add_assign(&mut self, other: Self) {
        self.c ^= other.c;
        for (i, v) in self.coefs.iter_mut().enumerate() {
            *v ^= other.coefs[i];
        }
    }
}

impl std::fmt::Display for Poly {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
         let mut start = true;
         for i in 0..self.n {
            for j in 0..=i {
                if self.coefs[i*(i+1)/2 + j] {
                    if !start {
                        write!(f, " + ")?;
                    }
                    else {
                        start = false;
                    }
                    if i == j {
                        write!(f, "x_{}", j + 1)?;
                    }
                    else {
                        write!(f, "x_{}x_{}", j + 1, i + 1)?;
                    }
                }
            }
         }
         if self.c {
            if !start {
                write!(f, " + 1")?;
            }
            else {
                write!(f, "1")?;
                start = false;
            }
         }
         if start {
            write!(f, "0")?;
         }

         Ok(())
    }
}

#[derive(Debug, Clone)]
struct System {
    eqs: Vec<Poly>,
    n: usize,
    m: usize,
    seed: u32
}

impl System {
    fn read<R: BufRead>(mut rdr: R) -> System {
        let mut lines = rdr.lines();
        assert_eq!(lines.next().unwrap().unwrap(), "Galois Field : GF(2)");
        let n = lines.next().unwrap().unwrap()
            .split(':').next_back().unwrap().trim()
            .parse::<usize>().unwrap();
        let m = lines.next().unwrap().unwrap()
            .split(':').next_back().unwrap().trim()
            .parse::<usize>().unwrap();
        let seed = lines.next().unwrap().unwrap()
            .split(':').next_back().unwrap().trim()
            .parse::<u32>().unwrap();
        while "*********************" != lines.next().unwrap().unwrap() {}

        let mut eqs = Vec::new();
        for _ in 0..m {
            let l = lines.next().unwrap().unwrap();
            eqs.push(Poly::parse(n, l.trim().strip_suffix(";").unwrap()));
        }

        System { eqs, n, m, seed }
    }

    fn check(&self, params: &[bool]) -> bool {
        for eq in &self.eqs {
            if eq.evaluate(params) {
                return false
            }
        }
        true
    }

    fn solve_rec(&self, params: &mut [bool], i: usize, skip: bool) -> bool {
        if i == self.n {
            return self.check(params)
        }

        if params[i] && i + 1 == self.n && skip {
            params[i] = false;
            return false
        }
        else if i + 1 == self.n && skip {
            params[i] = true;
            if self.solve_rec(params, i+1, skip) {
                return true
            }
            params[i] = false;
            return false;
        }
        else if params[i] {
            if self.solve_rec(params, i+1, skip) {
                return true
            }
            params[i] = false;
            return false;
        }
        
        for v in &[false, true] {
            params[i] = *v;
            if self.solve_rec(params, i+1, skip) {
                return true
            }
        }
        params[i] = false;

        false
    }

    fn solve(&self) -> Option<Vec<bool>> {
        let mut params = vec![false; self.n];
        if self.solve_rec(&mut params, 0, false) {
            Some(params)
        }
        else {
            None
        }
    }

    fn solve_all(&self) -> Vec<Vec<bool>> {
        if (self.n == 0) {
            return vec![Vec::new()]
        }
        let mut params = vec![false; self.n];
        let mut solutions = Vec::new();
        let mut skip = false;
        while self.solve_rec(&mut params, 0, skip) {
            solutions.push(params.clone());
            skip = true;
        }
        solutions
    }

    fn fix_vars(&self, vars: &[(usize, bool)]) -> System {
        let mut system = self.clone();
        system.n -= vars.len();
        system.eqs = system.eqs.into_iter()
            .map(|p| p.fix_vars(vars))
            .collect();
        system
    }

    fn add_affine_eqs(&self, k: usize) -> System {
        let mut system = self.clone();
        system.m += k;
        for _ in 0..k {
            system.eqs.push(Poly::random_affine(system.n));
        }
        system
    }

    fn new_random(&self, l: usize) -> System {
        let mut eqs = vec![Poly::zero(self.n); l];
        for eq in eqs.iter_mut() {
            for p in self.eqs.iter() {
                if rand::random() {
                    *eq += p.to_owned();
                }
            }
        }

        System  {
            eqs,
            n: self.n,
            m: l,
            seed: self.seed
        }
        
    }
}

impl std::fmt::Display for System {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for eq in &self.eqs {
            write!(f, "0 = {}\n", eq)?;
        }
        Ok(())
    }
}

fn find_solution(system: System) -> Vec<bool> {
    let mut found_variables = Vec::new();
    for i in 0..system.n {
        for x in std::iter::repeat(&[false, true]).flatten() {
            dbg!(i, x);
            let mut vars = found_variables.clone();
            vars.push((i, *x));
            if has_solution(system.fix_vars(&vars)) {
                found_variables.push((i, *x));
                break;
            }
        }
    }
    let mut res = vec![false; system.n];
    for (i, val) in found_variables {
        res[i] = val;
    }
    res
}

fn has_solution(system: System) -> bool {
    for k in 0..=system.n {
        let system2 = system.add_affine_eqs(k);
        if parity_count(system2) {
            return true
        }
    }
    false
}


fn parity_count_bf(system: System) -> bool {
    system.solve_all().len() % 2 == 1
}

fn parity_count(system: System) -> bool {
    // let k0 = 0.7;
    let n = system.n;
    let n1 = system.n*3 / 10;
    let v = multi_parity_count(system, n1, n - n1);
    v.into_iter().fold(false, |acc, v| acc ^ v)
}

fn multi_parity_count(system: System, n1: usize, w: usize) -> Vec<bool> {
    //let lambda = 0.1;
    //let n2 = (n1 - lambda*n).floor();
    let n = system.n;
    let n2 = (n1*10 - std::cmp::min(n1*10, n)) / 10;
    let l = n2 + 2;
    let t = 48*n + 1;
    dbg!(n, n1, w, n2, l);
    if n2 <= 0 {
        return multi_parity_count_bf(system, n1, w)
    }
    let mut sb = vec![0i32; ncre(n - n1, w)*(1 << (n1 - n2))];
    
    for k in 1..=t {
        let system2 = system.new_random(l);
        let v1 = multi_parity_count(system2, n2, 2*l - n2);
        dbg!(n, n1, n2, v1.len(), 1 << (n - n2), ncre(n-n2, 2*l - n2));
        // interpolate G_k(y, u): Möbius transform
        let g = fmt(&v1, n - n2, 2*l - n2);
        //dbg!(g);
    }

    //dbg!(sb);

    let mut vs = vec![false; ncre(n - n1, w)];
    vs
}

fn multi_parity_count_bf(system: System, n1: usize, w: usize) -> Vec<bool> {
    dbg!(system.n, n1, w, ncre(system.n, w));
    let mut evals = vec![true; ncre(system.n - n1, w)*(1 << n1)];
    for j in 0..system.m {
        Bitset::hamming_weight_iter(system.n - n1, w as u32)
            .flat_map(|s1| Bitset::full_iter(n1).map(move |s2| s1.append(s2, n1)))
            .map(|s| system.eqs[j].evaluate(&s.to_boolvec(system.n)))
            .enumerate()
            .for_each(|(i, x)| evals[i] &= !x);
    }
    let mut v = vec![false; ncre(system.n - n1, w)];
    for (i, vi) in v.iter_mut().enumerate() {
        for j in 0..(1 << n1) {
            *vi ^= evals[i*(1<< n1) + j];
        }
    }
    v
}

fn ncr(n: usize, r: usize) -> usize {
    if r == 0 {
        1
    }
    else {
        n * ncr(n-1, r-1) / r
    }
}

fn ncre(n: usize, r: usize) -> usize {
    (0..=r).map(|r| ncr(n, r)).sum()
}

// Fast möbius transform
// Implementation inspired by https://codeforces.com/blog/entry/72488
// Potential new reference: https://bitbucket.org/fes/fes/src/master/src/moebius_transform.c
fn fmt(f: &[bool], n: usize, w: usize) -> Vec<bool> {
    let mapping = Bitset::hamming_weight_iter(n, w as u32)
        .enumerate()
        .map(|(i, s)| (s, i))
        .collect::<BTreeMap<_, _>>();
    dbg!(mapping.len());
    let mut g = f.to_owned();
    for i in 0..n {
        for mask in Bitset::hamming_weight_iter(n, w as u32) {
            if mask.get_bit(i) {
                let mut mask2 = mask.clone();
                mask2.set_bit(i, false);
                g[mapping[&mask]] ^= g[mapping[&mask2]];
            }
        }
    }
    g
}


fn main() {
    let file = std::fs::File::open("example.txt").unwrap();
    //let file = std::fs::File::open("toy_example/ToyExample-type4-n15-seed0").unwrap();
    let system = System::read(std::io::BufReader::new(file));
    eprintln!("{}", &system);
    eprintln!("{}", &system.fix_vars(&[(0, true), (1, true)]));
    dbg!(system.solve_all());
    
    /*
    */
    dbg!(has_solution(system.clone()));
    let solution = dbg!(find_solution(system.clone()));
    dbg!(system.check(&solution));


    /*
    let mut evals = vec![true; 1 << system.n];
    for j in 0..system.m {
        for (i, s) in Bitset::full_iter(system.n).enumerate() {
            let v = s.to_boolvec(system.n);
            let r = system.eqs[j].evaluate(&v);
            dbg!(&v, r);
            evals[i] &= !r;
        }
    }
    dbg!(&evals);
    dbg!(system.solve());
    */
}
