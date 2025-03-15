use std::cmp:: {
    Ordering,
    min,
    max,
};

use crate::primitives::*;

use crate::basis::{
    Basis,
};

// 2^17
const INITIAL_HASH_TABLE_SIZE: usize = 131072;

pub struct Monomial {
    pub degree: Degree,
    divisor_mask: DivisorMask,
    pub exponents: ExpVec,
    last_known_divisor: BasisLength,
}

pub struct HashTable {
    pub monomials    : Vec<Monomial>,
    random_seed      : Vec<HashValue>,
    values           : Vec<HashValue>,
    map              : Vec<HashTableLength>,
    divisor_bounds   : ExpVec,
    pub indices      : Vec<HashTableLength>,
    pub nr_variables : usize,
    length           : usize,
}


impl HashTable {
    pub fn new(initial_exponents: &Vec<Vec<ExpVec>>) -> HashTable {
        debug_assert!(initial_exponents.len() > 0);
        debug_assert!(initial_exponents[0].len() > 0);
        debug_assert!(initial_exponents[0][0].len() > 0);
        let mut ht = HashTable {
            monomials      : Vec::new(),
            random_seed    : Vec::new(),
            values         : vec![0; INITIAL_HASH_TABLE_SIZE],
            map            : vec![HashTableLength::MAX; INITIAL_HASH_TABLE_SIZE],
            divisor_bounds : Vec::new(),
            indices        : vec![0; INITIAL_HASH_TABLE_SIZE],
            nr_variables   : initial_exponents[0][0].len(),
            length         : INITIAL_HASH_TABLE_SIZE,
        };
        ht.generate_random_seed(ht.nr_variables);
        ht.generate_divisor_bounds(initial_exponents);

        return ht;
    }

    fn update_seed(&self, mut seed: HashValue) -> HashValue {
        seed ^= seed << 13;
        seed ^= seed >> 17;
        seed ^= seed << 5;
        return seed;
    }

    // We assume all exponents from a are greater or equal than
    // the corresponding exponents from b. Used to get multiples
    // during symbolic preprocessing, usually a represents the lcm
    // of b with some other monommial.
    pub fn get_difference(
        &mut self, mon1: HashTableLength, mon2:HashTableLength)
        -> HashTableLength {
        let lm_e = &self.monomials[mon1 as usize].exponents;
        let lm_f = &self.monomials[mon2 as usize].exponents;
        debug_assert!(lm_e.into_iter()
            .zip(lm_f)
            .all(|(a,b)| *a>=*b));
        return self.insert(lm_e.into_iter()
            .zip(lm_f)
            .map(|(a,b)| *a-*b)
            .collect());
    }

    pub fn get_lcm(
        &mut self, mon1: HashTableLength, mon2:HashTableLength)
        -> HashTableLength {
        let lm_e = &self.monomials[mon1 as usize].exponents;
        let lm_f = &self.monomials[mon2 as usize].exponents;
        return self.insert(lm_e.into_iter()
            .zip(lm_f)
            .map(|(a,b)| max(*a,*b))
            .collect());
    }

    pub fn are_monomials_coprime(
        &self, mon1: HashTableLength, mon2:HashTableLength)
        -> bool {
        let lm_e = &self.monomials[mon1 as usize].exponents;
        let lm_f = &self.monomials[mon2 as usize].exponents;
        if lm_e.into_iter().zip(lm_f).any(|(a,b)| *a!=0 && *b!=0) {
            return false;
        }  else {
            return true;
        }
    }

    fn generate_random_seed(&mut self, nr_variables: usize) {
        let mut seed: HashValue = 2463534242;
        for _i in 0..nr_variables {
            seed = self.update_seed(seed);
            self.random_seed.push(seed);
        }
    }

    // The divisor mask template is generated once the
    // input polynomials are read in.
    // TODO: Opptimize the divisor mask: If #variables < usize::BITS we leave
    // a part of the divisor mask 0 and do not use it.
    fn generate_divisor_bounds(&mut self, initial_exponents: &Vec<Vec<ExpVec>>) {
        let length_divmask = min(
            self.nr_variables, usize::BITS.try_into().unwrap());

        for i in 0..length_divmask {
            let mut max = initial_exponents[0][0][i];
            let mut min = initial_exponents[0][0][i];
            for e in initial_exponents.into_iter().flatten() {
                if e[i] > max {
                    max = e[i];
                    continue;
                } else if e[i] < min {
                    min = e[i];
                }
            }
            self.divisor_bounds.push(max-min);
        }
    }

    fn get_divisor_mask(&mut self, exp: &ExpVec) -> DivisorMask {
        let divisor_bounds = &self.divisor_bounds;
        let mut divisor_mask = 0usize;
        let min_len = min(exp.len(), divisor_bounds.len());
        for i in 0..min_len {
            if exp[i] >= divisor_bounds[i] {
                divisor_mask |= 1 << i;
            }
        }
        return divisor_mask;
    }

    // Tests if monomial ma divides monomial mb
    pub fn divides(&self, ma: HashTableLength, mb: HashTableLength) -> bool {
        if self.monomials[ma as usize].divisor_mask & !self.monomials[mb as usize].divisor_mask != 0 {
            return false;
        }
        let ea = &self.monomials[ma as usize].exponents;
        let eb = &self.monomials[mb as usize].exponents;
        if ea.into_iter().zip(eb).any(|(a,b)| a > b) {
            return false;
        }  else {
            return true;
        }
    }

    pub fn find_divisor(&mut self, mon: HashTableLength, basis: &Basis)
        -> Option<(BasisLength, HashTableLength)> {

        let start_idx = self.monomials[mon as usize].last_known_divisor;
        for (bi, be) in basis.elements[start_idx..].iter().enumerate() {
            if !be.is_redundant && self.divides(be.monomials[0], mon) {
                self.monomials[mon as usize].last_known_divisor = bi + start_idx;
                self.indices[mon as usize] = 2;
                return Some((bi+start_idx, self.get_difference(mon, be.monomials[0])));
            }
        }
        return None;
    }

    pub fn generate_multiplied_monomials(&mut self, divisor_idx: BasisLength,
        mult_idx: HashTableLength, basis: &Basis) -> MonomVec {

        let vec_len = self.nr_variables;
        let mons = &basis.elements[divisor_idx].monomials;
        let mut mult_mons: MonomVec = vec!(0; mons.len());
        for (idx, m) in mons.iter().enumerate() {
            let mut exps: ExpVec = vec!(0; vec_len);
            for i in 0..vec_len {
                exps[i] = self.monomials[mult_idx as usize].exponents[i]
                    + self.monomials[*m as usize].exponents[i];
            }
            mult_mons[idx] = self.insert(exps);
        }
        return mult_mons;
    }

    fn get_hash(&self, exp: &ExpVec) -> HashValue {
        let mut h: HashValue = 0;
        for i in 0..exp.len() {
            h = h.wrapping_add((exp[i] as HashValue).wrapping_mul(self.random_seed[i]));
        }
        return h;
    }

    pub fn cmp_monomials_by_degree(&self, a: HashTableLength, b:HashTableLength) -> Ordering {
        debug_assert!(
            self.monomials[a as usize].exponents.len() == self.monomials[b as usize].exponents.len());
        let ma = &self.monomials[a as usize];
        let mb = &self.monomials[b as usize];
        if ma.degree != mb.degree {
            if ma.degree > mb.degree { return Ordering::Greater; }
            else { return Ordering::Less; }
        }
        return Ordering::Equal;
    }

    pub fn cmp_monomials_by_drl(&self, a: HashTableLength, b:HashTableLength) -> Ordering {
        debug_assert!(
            self.monomials[a as usize].exponents.len() == self.monomials[b as usize].exponents.len());
        let ma = &self.monomials[a as usize];
        let mb = &self.monomials[b as usize];
        if ma.degree != mb.degree {
            if ma.degree > mb.degree { return Ordering::Greater; }
            else { return Ordering::Less; }
        }
        // From here on, we know that the degrees are the same
        for i in (0..ma.exponents.len()).rev() {
            if ma.exponents[i] < mb.exponents[i] { return Ordering::Greater; }
            else if ma.exponents[i] > mb.exponents[i] { return Ordering::Less; }
        }
        return Ordering::Equal;
    }

    fn enlarge(&mut self) {
        self.length *= 2;
        self.map = vec![HashTableLength::MAX; self.length];
        self.values.resize(self.length, 0);
        self.indices.resize(self.length, 0);

        // reinsert elements
        let div = (self.length - 1) as HashTableLength;

        for i in 0..self.monomials.len() {
            let h = self.values[i];
            let mut k = h;
            for j in 0..self.length {
                k = (k+j as HashTableLength) & div;
                if self.map[k as usize] == HashTableLength::MAX {
                    self.map[k as usize] = i as HashTableLength;
                    break;
                }
            }
        }
    }

    pub fn insert(&mut self, exp: ExpVec) -> HashTableLength {
        let div = (self.length - 1) as HashTableLength;
        let h = self.get_hash(&exp);
        let mut k = h;
        let mut i = 0;
        while i < self.map.len() {
            k = (k+i as HashTableLength) & div;
            let hm = self.map[k as usize];
            if hm == HashTableLength::MAX {
                break;
            }
            if self.values[hm as usize] != h {
                i += 1;
                continue;
            }
            let eh = &self.monomials[hm as usize].exponents;
            if *eh != exp {
                i += 1;
                continue;
            }
            return hm;
        }
        let pos = self.monomials.len() as HashTableLength;
        self.map[k as usize] = pos;
        let monomial = Monomial {
            degree: get_degree(&exp),
            divisor_mask: self.get_divisor_mask(&exp),
            exponents: exp,
            last_known_divisor: 0,
        };
        self.monomials.push(monomial);
        self.values[pos as usize] = h;

        // enlarge AFTER insertion to correctly adjust
        // the just inserted element
        if self.monomials.len() == self.length {
            self.enlarge();
        }

        return pos;
    }
}

fn get_degree(exp: &ExpVec) -> Degree {
    let mut deg = 0;
    for e in exp {
        deg += *e as Degree;
    }
    return deg;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_enlarge() {
        let exp: Vec<Vec<ExpVec>> = vec!(vec!(vec!(1,1,1,1,1)));
        let mut ht = HashTable::new(&exp);
        assert_eq!(ht.monomials.len(), 0);
        assert_eq!(ht.indices.len(), INITIAL_HASH_TABLE_SIZE);
        assert_eq!(ht.map.len(), INITIAL_HASH_TABLE_SIZE);
        assert_eq!(ht.values.len(), INITIAL_HASH_TABLE_SIZE);
        ht.enlarge();
        assert_eq!(ht.indices.len(), 2*INITIAL_HASH_TABLE_SIZE);
        assert_eq!(ht.map.len(), 2*INITIAL_HASH_TABLE_SIZE);
        assert_eq!(ht.values.len(), 2*INITIAL_HASH_TABLE_SIZE);
        assert_eq!(ht.indices.iter().all(|a| *a == 0), true);
        assert_eq!(ht.values.iter().all(|a| *a == 0), true);

    }
    #[test]
    fn test_random_seed() {
        let exp: Vec<Vec<ExpVec>> = vec!(vec!(vec!(1,1,1,1,1)));
        let ht = HashTable::new(&exp);
        assert_eq!(ht.random_seed,
            [723471715, 2497366906, 2064144800, 2008045182, 3532304609]);
    }
    #[test]
    fn test_cmp_monomials_by_degree() {
        let exps: Vec<Vec<ExpVec>> = vec!(vec!(
            vec![1,1,3],
            vec![2,1,3],
            vec![2,0,3]));
        let mut ht = HashTable::new(&exps);
        for e in exps.into_iter().flatten() {
            ht.insert(e);
        }
        assert_eq!(ht.cmp_monomials_by_degree(0,1), Ordering::Less);
        assert_eq!(ht.cmp_monomials_by_degree(1,2), Ordering::Greater);
        assert_eq!(ht.cmp_monomials_by_degree(2,0), Ordering::Equal);
    }
    #[test]
    fn test_cmp_monomials_by_drl() {
        let exps: Vec<Vec<ExpVec>> = vec!(vec!(
            vec![1,1,3],
            vec![2,1,3],
            vec![2,0,3]));
        let mut ht = HashTable::new(&exps);
        for e in exps.into_iter().flatten() {
            ht.insert(e);
        }
        assert_eq!(ht.cmp_monomials_by_drl(0,1), Ordering::Less);
        assert_eq!(ht.cmp_monomials_by_drl(0,0), Ordering::Equal);
        assert_eq!(ht.cmp_monomials_by_drl(2,0), Ordering::Greater);
    }
    #[test]
    fn test_insert() {
        let exp: Vec<Vec<ExpVec>> = vec!(vec!(vec![1,1,1]));
        let mut ht = HashTable::new(&exp);
        for e in exp.into_iter().flatten() {
            ht.insert(e);
        }
        let pos = ht.insert(vec![1 as Exponent,1,1]);
        assert_eq!(pos, 0);
    }
    #[test]
    fn test_find_divisor() {
        let fc : Characteristic = 65521;
        let cfs : Vec<CoeffVec> = vec![vec![-2,65523], vec![1, -3]];
        let exps : Vec<Vec<ExpVec>> = vec![vec![vec![0,3], vec![1,1]], vec![vec![0,2], vec![1,1]]];
        let mut hash_table = HashTable::new(&exps);
        let basis = Basis::new::<i32>(&mut hash_table, fc, cfs, exps);
        let mon1 = hash_table.insert(vec![0 as Exponent,4]);
        assert_eq!(hash_table.find_divisor(mon1, &basis), Some((1, 4)));
        assert_eq!(hash_table.monomials[4].exponents, [0,1]);
        let mon2 = hash_table.insert(vec![5 as Exponent,0]);
        assert_eq!(hash_table.find_divisor(mon2, &basis), None);
    }
    #[test]
    fn test_generate_divisor_bounds() {
        let exps: Vec<Vec<ExpVec>> = vec!(vec!(
            vec![1,1,3],
            vec![2,0,3]));
        let ht = HashTable::new(&exps);
        assert_eq!(ht.divisor_bounds, [1,1,0]);
    }
    #[test]
    fn test_get_divisor_mask() {
        let exps: Vec<Vec<ExpVec>> = vec!(vec!(
            vec![1;65]));
        let mut ht = HashTable::new(&exps);
        println!("{:?}", ht.get_divisor_mask(&exps[0][0]));
        // assert_eq!(ht.divisor_bounds, [1,1,0]);
    }
    #[test]
    fn test_init_hash_table() {
        let exps: Vec<Vec<ExpVec>> = vec!(vec!(
            vec![1,1,3],
            vec![2,0,4]));
        let mut ht = HashTable::new(&exps);
        for e in exps.into_iter().flatten() {
            ht.insert(e);
        }
        assert_eq!(*ht.monomials[0].exponents, [1,1,3]);
        assert_eq!(*ht.monomials[1].exponents, [2,0,4]);
    }
    #[test]
    fn test_divides() {
        let exps: Vec<Vec<ExpVec>> = vec!(vec!(vec![1,1,3], vec![2,0,3]));
        let mut ht = HashTable::new(&exps);
        let ma = ht.insert(exps[0][0].clone());
        let mb = ht.insert(vec![1,2,3]);
        let mc = ht.insert(vec![3,0,4]);
        assert_eq!(ht.divides(ma, mb), true);
        assert_eq!(ht.divides(ma, mc), false);
        assert_eq!(ht.divides(ma, ma), true);
    }
    #[test]
    #[should_panic]
    fn test_get_difference_panic() {
        let exps: Vec<Vec<ExpVec>> = vec!(vec!(
            vec![1,1,3],
            vec![2,1,2]));
        let mut ht = HashTable::new(&exps);
        ht.insert(exps[0][0].clone());
        ht.insert(exps[0][1].clone());
        let _diff = ht.get_difference(1,0);
    }
    #[test]
    fn test_get_difference() {
        let exps: Vec<Vec<ExpVec>> = vec!(vec!(
            vec![1,1,3],
            vec![2,1,4]));
        let mut ht = HashTable::new(&exps);
        ht.insert(exps[0][0].clone());
        ht.insert(exps[0][1].clone());
        let diff = ht.get_difference(1,0);
        assert_eq!(ht.monomials[diff as usize].exponents, [1,0,1]);
    }
    #[test]
    fn test_get_lcm() {
        let exps: Vec<Vec<ExpVec>> = vec!(vec!(
            vec![1,1,3],
            vec![2,0,3]));
        let mut ht = HashTable::new(&exps);
        ht.insert(exps[0][0].clone());
        ht.insert(exps[0][1].clone());
        let lcm = ht.get_lcm(0,1);
        assert_eq!(ht.monomials[lcm as usize].exponents, [2,1,3]);
    }
    #[test]
    fn test_are_monomials_prime() {
        let exps: Vec<Vec<ExpVec>> = vec!(vec!(
            vec![0,1,0],
            vec![2,0,3],
            vec![2,1,0]));
        let mut ht = HashTable::new(&exps);
        // let te : Vec<_> = exps[0][0].iter().zip(exps[0][1].clone()).map(|(a,b)| a+b).collect();
        let m1 = ht.insert(exps[0][0].clone());
        let m2 = ht.insert(exps[0][1].clone());
        let m3 = ht.insert(exps[0][1].clone());
        assert_eq!(ht.are_monomials_coprime(m1, m2), true);
        assert_eq!(ht.are_monomials_coprime(m3, m2), false);
    }
}
