#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd)]
pub struct Bitset {
    inner: u64
}

impl Bitset {
    pub fn from_u64(n: u64) -> Self {
        Bitset { inner: n }
    }

    pub fn full_iter(n: usize) -> FullIter {
        FullIter {
            v: Bitset::from_u64(0),
            n,
            first: true
        }
    }

    pub fn hamming_weight_iter(n: usize, w: u32) -> HammingWeightIter {
        HammingWeightIter {
            v: Bitset::from_u64(0),
            w,
            n,
            first: true
        }
    } 

    pub fn append(&self, other: Bitset, n: usize) -> Bitset {
        Bitset {
            inner: self.inner | (other.inner << n)
        }
    }

    pub fn set_bit(&mut self, n: usize, v: bool) {
        if v {
            self.inner |= 1 << n;
        }
        else {
            self.inner &= !(1 << n);
        }
    }

    pub fn get_bit(&self, n: usize) -> bool {
        self.inner & (1 << n) != 0
    }

    pub fn hamming_weight(&self) -> u32 {
        self.inner.count_ones()
    }

    pub fn to_boolvec(&self, n: usize) -> Vec<bool> {
        (0..n).map(|i| self.get_bit(i)).collect()
    }
}

#[derive(Clone)]
pub struct HammingWeightIter {
    v: Bitset,
    w: u32,
    n: usize,
    first: bool
}

impl Iterator for HammingWeightIter {
    type Item = Bitset;

    fn next(&mut self) -> Option<Bitset> {
        if self.v.inner == 0 {
            if !self.first {
                return None
            }
            self.first = false;
        }
        let res = self.v.clone();
        if self.v.hamming_weight() == self.w {
            let x = self.v.inner & !(self.v.inner.wrapping_sub(1));
            self.v.inner = (self.v.inner + x) & ((1 << self.n) - 1);
        }
        else {
            self.v.inner = (self.v.inner + 1) & ((1 << self.n) - 1);
        }
        Some(res)
    }
}

#[derive(Clone)]
pub struct FullIter {
    v: Bitset,
    n: usize,
    first: bool
}

impl Iterator for FullIter {
    type Item = Bitset;

    fn next(&mut self) -> Option<Bitset> {
        if self.v.inner == 0 {
            if !self.first {
                return None
            }
            self.first = false;
        }
        let res = self.v.clone();
        self.v.inner = (self.v.inner + 1) & ((1 << self.n) - 1);
        Some(res)
    }
}
