use std::cmp:: {
    Ordering,
    max,
};

use crate::primitives::*;

use crate::basis:: {
    Basis,
};

// 2^17
const INITIAL_HASH_TABLE_SIZE: usize = 131072;

pub struct HashTable {
    pub degrees      : Vec<Degree>,
    pub divisor_masks: Vec<DivisorMask>,
    last_known_divisors: Vec<BasisLength>,
    pub exponents    : Vec<ExpVec>,
    random_seed      : Vec<HashValue>,
    values           : Vec<HashValue>,
    map              : Vec<HashTableLength>,
    divisor_bounds   : ExpVec,
    pub indices      : Vec<HashTableLength>,
    pub nr_variables : usize,
    length           : usize,
    divisor_mask_bits_per_variable: usize,
    divisor_mask_variable_range: usize,
}


impl HashTable {
    pub fn new(initial_exponents: &Vec<Vec<ExpVec>>) -> HashTable {
        debug_assert!(initial_exponents.len() > 0);
        debug_assert!(initial_exponents[0].len() > 0);
        debug_assert!(initial_exponents[0][0].len() > 0);
        let mut ht = HashTable {
            exponents      : Vec::new(),
            random_seed    : Vec::new(),
            degrees        : vec![0; INITIAL_HASH_TABLE_SIZE],
            divisor_masks  : vec![0; INITIAL_HASH_TABLE_SIZE],
            last_known_divisors : vec![0; INITIAL_HASH_TABLE_SIZE],
            values         : vec![0; INITIAL_HASH_TABLE_SIZE],
            map            : vec![0; 2*INITIAL_HASH_TABLE_SIZE],
            divisor_bounds : Vec::new(),
            indices        : vec![0; INITIAL_HASH_TABLE_SIZE],
            nr_variables   : initial_exponents[0][0].len(),
            length         : INITIAL_HASH_TABLE_SIZE,
            divisor_mask_bits_per_variable: 0,
            divisor_mask_variable_range: 0,
        };
        ht.exponents.push(vec!(0; ht.nr_variables));
        // let ev: ExpVec = vec!(0; INITIAL_HASH_TABLE_SIZE * ht.nr_variables);
        // for i in 0..INITIAL_HASH_TABLE_SIZE {
        //     println!("i {}", i);
        //     ht.exponents.push(ev[i*ht.nr_variables..(i+1)*ht.nr_variables-1].to_vec());
        // }
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
        let lm_e = &self.exponents[mon1 as usize];
        let lm_f = &self.exponents[mon2 as usize];
        debug_assert!(lm_e.into_iter().zip(lm_f).all(|(a,b)| *a>=*b));
        return self.insert(lm_e.into_iter()
            .zip(lm_f)
            .map(|(a,b)| *a-*b)
            .collect());
    }

    pub fn get_lcm(
        &mut self, mon1: HashTableLength, mon2:HashTableLength)
        -> HashTableLength {
        let lm_e = &self.exponents[mon1 as usize];
        let lm_f = &self.exponents[mon2 as usize];
        return self.insert(lm_e.into_iter()
            .zip(lm_f)
            .map(|(a,b)| max(*a,*b))
            .collect());
    }

    pub fn are_monomials_coprime(
        &self, mon1: HashTableLength, mon2:HashTableLength)
        -> bool {
        let lm_e = &self.exponents[mon1 as usize];
        let lm_f = &self.exponents[mon2 as usize];
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
        // self.random_seed.reverse();
    }

    // The divisor mask template is generated once the
    // input polynomials are read in.
    // TODO: Opptimize the divisor mask: If #variables < usize::BITS we leave
    // a part of the divisor mask 0 and do not use it.
    fn generate_divisor_bounds(&mut self, initial_exponents: &Vec<Vec<ExpVec>>) {
        let bits_for_divisor_mask:usize =  DivisorMask::BITS.try_into().unwrap();
        self.divisor_mask_bits_per_variable = bits_for_divisor_mask / self.nr_variables;
        if self.divisor_mask_bits_per_variable == 0 {
            self.divisor_mask_bits_per_variable = 1;
        }
        if self.nr_variables < bits_for_divisor_mask {
            self.divisor_mask_variable_range = self.nr_variables;
        } else {
            self.divisor_mask_variable_range = bits_for_divisor_mask;
        }
        println!("bpv {}", self.divisor_mask_bits_per_variable);
        println!("dvr {}", self.divisor_mask_variable_range);
        for i in 0..self.divisor_mask_variable_range {
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
            let mut ctr = (max - min) / self.divisor_mask_bits_per_variable as u16;
            if ctr == 0 {
                ctr = 1;
            }
            for _j in 0..self.divisor_mask_bits_per_variable {
                self.divisor_bounds.push(ctr);
                ctr += 1;
            }
        }
    }

    fn get_divisor_mask(&mut self, exp: &ExpVec) -> DivisorMask {
        let divisor_bounds = &self.divisor_bounds;
        let mut divisor_mask: DivisorMask = 0;
        let mut ctr = 0;
        for i in 0..self.divisor_mask_variable_range {
            for _j in 0..self.divisor_mask_bits_per_variable {
                if exp[i] >= divisor_bounds[ctr] {
                    divisor_mask |= 1 << ctr;
                }
                ctr += 1;
            }
        }
        return divisor_mask;
    }

    // Tests if monomial ma divides monomial mb
    pub fn divides_pairs(&self, ma: HashTableLength, mb: HashTableLength) -> bool {
        if (self.divisor_masks[ma as usize] & !self.divisor_masks[mb as usize]) != 0 {
            return false;
        }
        let ea = &self.exponents[ma as usize];
        let eb = &self.exponents[mb as usize];
        if ea.into_iter().zip(eb).all(|(a,b)| *a <= *b) {
            return true;
        }  else {
            return false;
        }
    }

    pub fn divides(&self, ma: HashTableLength, dma: DivisorMask, mb: HashTableLength, ndmb: DivisorMask) -> bool {
        if (dma & ndmb) != 0 {
            return false;
        }
        let ea = &self.exponents[ma as usize];
        let eb = &self.exponents[mb as usize];
        if ea.into_iter().zip(eb).all(|(a,b)| *a <= *b) {
            return true;
        }  else {
            return false;
        }
    }

    pub fn find_divisor(&mut self, mon: HashTableLength,
            divsv: &Vec<(DivisorMask,HashTableLength,BasisLength)>)
        -> Option<(BasisLength, HashTableLength)> {

        let divs = divsv.as_slice();
        let start_idx = self.last_known_divisors[mon as usize];
        let ndmon = !self.divisor_masks[mon as usize];
        let mut j = 0;
        while divs[j].2 < start_idx {
            j += 1;
        }
        for i in j..divs.len() {
            let d = divs[i];
            if self.divides(d.1, d.0, mon, ndmon) {
                self.last_known_divisors[mon as usize] = d.2;
                return Some((d.2, self.get_difference(mon, d.1)));
            }
        }
        return None;
    }

    pub fn generate_multiplied_monomials(&mut self, divisor_idx: BasisLength,
        mult_idx: HashTableLength, basis: &Basis) -> MonomVec {

        let mons = &basis.elements[divisor_idx as usize].monomials;
        let mut mult_mons: MonomVec = vec!(0; mons.len());
        for (idx, m) in mons.iter().enumerate() {
            let e_mon = &self.exponents[*m as usize];
            let e_mult = &self.exponents[mult_idx as usize];
            let mut exps: ExpVec = vec!(0; self.nr_variables);
            for ((e, mu), mo) in exps.iter_mut().zip(e_mult).zip(e_mon) {
                *e = mu +mo;
            }
            mult_mons[idx] = self.insert(exps);
        }
        return mult_mons;
    }

    #[inline(always)]
    fn get_hash(&self, exp: &ExpVec) -> HashValue {
        return exp.iter().zip(&self.random_seed)
            .map(|(e,r) | (*e as HashTableLength).wrapping_mul(*r))
            .fold(0, |acc, x| acc.wrapping_add(x));
    }

    pub fn cmp_monomials_by_degree(&self, a: HashTableLength, b:HashTableLength) -> Ordering {
        debug_assert!(
            self.exponents[a as usize].len() == self.exponents[b as usize].len());
        let da = &self.degrees[a as usize];
        let db = &self.degrees[b as usize];
        if da != db {
            if da > db { return Ordering::Greater; }
            else { return Ordering::Less; }
        }
        return Ordering::Equal;
    }

    pub fn cmp_monomials_by_drl(&self, a: HashTableLength, b:HashTableLength) -> Ordering {
        debug_assert!(
            self.exponents[a as usize].len() == self.exponents[b as usize].len());
        let da = &self.degrees[a as usize];
        let db = &self.degrees[b as usize];
        if da != db {
            if da > db { return Ordering::Greater; }
            else { return Ordering::Less; }
        }
        // From here on, we know that the degrees are the same
        let ea = &self.exponents[a as usize];
        let eb = &self.exponents[b as usize];
        for i in (0..ea.len()).rev() {
            if ea[i] < eb[i] { return Ordering::Greater; }
            else if ea[i] > eb[i] { return Ordering::Less; }
        }
        return Ordering::Equal;
    }

    fn enlarge(&mut self) {
        let previous_length = self.length;
        self.length *= 2;
        self.map = vec!(0; 2*self.length);
        self.degrees.resize(self.length, 0);
        self.divisor_masks.resize(self.length, 0);
        self.last_known_divisors.resize(self.length, 0);
        self.values.resize(self.length, 0);
        self.indices.resize(self.length, 0);

        let map_len = self.map.len();
        // reinsert elements
        let div = (map_len - 1) as HashTableLength;
        for i in 1..previous_length {
            let mut k = self.values[i];
            for j in 0..map_len {
                k = (k+j as HashTableLength) & div;
                if self.map[k as usize] == 0 {
                    self.map[k as usize] = i as HashTableLength;
                    break;
                }
            }
        }
    }

    #[inline(always)]
    pub fn insert(&mut self, exp: ExpVec) -> HashTableLength {
        let div = (self.map.len() - 1) as HashTableLength;
        let h = self.get_hash(&exp);
        let mut k = h;
        let map_len = self.map.len();
        for  i in 0..map_len {
            k = (k+i as HashTableLength) & div;
            let hm = self.map[k as usize];
            if hm == 0 {
                break;
            }
            if self.values[hm as usize]!= h || self.exponents[hm as usize] != exp {
                continue;
            }
            return hm;
        }
        let pos = self.exponents.len();
        self.map[k as usize] = pos as HashTableLength;
        self.degrees[pos] = get_degree(&exp);
        self.divisor_masks[pos] = self.get_divisor_mask(&exp);
        // self.last_known_divisors[pos] = 0;
        self.exponents.push(exp);
        self.values[pos] = h;

        // enlarge AFTER insertion to correctly adjust
        // the just inserted element
        if self.exponents.len() == self.length {
            self.enlarge();
        }

        return pos as HashTableLength;
    }
}

fn get_degree(exp: &ExpVec) -> Degree {
    return exp.iter().fold(0, |acc, x| acc.wrapping_add(*x as Degree));
}

#[cfg(test)]
mod tests {
    use rand::Rng;
    use super::*;

    #[test]
    fn test_enlarge() {
        let exp: Vec<Vec<ExpVec>> = vec!(vec!(vec!(1,1,1,1,1)));
        let mut ht = HashTable::new(&exp);
        assert_eq!(ht.exponents.len(), 1);
        assert_eq!(ht.indices.len(), INITIAL_HASH_TABLE_SIZE);
        assert_eq!(ht.map.len(), 2*INITIAL_HASH_TABLE_SIZE);
        assert_eq!(ht.values.len(), INITIAL_HASH_TABLE_SIZE);

        // add random data to hash table, otherwise enlargement
        // is not useful and slow
        let mut rng = rand::rng();
        for _i in 0..INITIAL_HASH_TABLE_SIZE-1 {
            let e = vec!(rng.random::<u16>(),
                            rng.random::<u16>(),
                            rng.random::<u16>(),
                            rng.random::<u16>(),
                            rng.random::<u16>());
            ht.insert(e);
        }
        assert_eq!(ht.indices.len(), 2*INITIAL_HASH_TABLE_SIZE);
        assert_eq!(ht.map.len(), 4*INITIAL_HASH_TABLE_SIZE);
        assert_eq!(ht.values.len(), 2*INITIAL_HASH_TABLE_SIZE);
        assert_eq!(ht.indices.iter().all(|a| *a == 0), true);
        assert_eq!(ht.values[INITIAL_HASH_TABLE_SIZE..].iter().all(|a| *a == 0), true);

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
        assert_eq!(ht.cmp_monomials_by_degree(1,2), Ordering::Less);
        assert_eq!(ht.cmp_monomials_by_degree(2,3), Ordering::Greater);
        assert_eq!(ht.cmp_monomials_by_degree(3,1), Ordering::Equal);
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
        assert_eq!(ht.cmp_monomials_by_drl(1,2), Ordering::Less);
        assert_eq!(ht.cmp_monomials_by_drl(1,1), Ordering::Equal);
        assert_eq!(ht.cmp_monomials_by_drl(3,1), Ordering::Greater);
    }
    #[test]
    fn test_insert() {
        let exp: Vec<Vec<ExpVec>> = vec!(vec!(vec![1,1,1]));
        let mut ht = HashTable::new(&exp);
        for e in exp.into_iter().flatten() {
            ht.insert(e);
        }
        let pos = ht.insert(vec![1 as Exponent,1,1]);
        assert_eq!(pos, 1);
    }
    #[test]
    fn test_find_divisor() {
        let fc : Characteristic = 65521;
        let cfs : Vec<CoeffVec> = vec![vec![-2,65523], vec![1, -3]];
        let exps : Vec<Vec<ExpVec>> = vec![vec![vec![0,3], vec![1,1]], vec![vec![0,2], vec![1,1]]];
        let mut hash_table = HashTable::new(&exps);
        let basis = Basis::new::<i32>(&mut hash_table, fc, cfs, exps);
        let mon1 = hash_table.insert(vec![0 as Exponent,4]);
        let mut divs: Vec<(DivisorMask,HashTableLength,BasisLength)> = Vec::new();
        for i in 0..basis.elements.len() {
            if basis.elements[i].is_redundant == false {
                divs.push((
                    hash_table.divisor_masks[basis.elements[i].monomials[0] as usize],
                    basis.elements[i].monomials[0],
                    i as BasisLength));
            }
        }
        assert_eq!(hash_table.find_divisor(mon1, &divs), Some((1, 5)));
        assert_eq!(hash_table.exponents[5], [0,1]);
        let mon2 = hash_table.insert(vec![5 as Exponent,0]);
        assert_eq!(hash_table.find_divisor(mon2, &divs), None);
    }
    #[test]
    fn test_generate_divisor_bounds() {
        let exps: Vec<Vec<ExpVec>> = vec!(vec!(
            vec![1,1,3],
            vec![2,0,3]));
        let ht = HashTable::new(&exps);
        assert_eq!(ht.divisor_bounds,
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
             1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
             1, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
    }
    #[test]
    fn test_get_divisor_mask() {
        let exps: Vec<Vec<ExpVec>> = vec!(vec!(
            vec![1;65]));
        let mut ht = HashTable::new(&exps);
        ht.insert(vec![1;65]);
        assert_eq!(ht.divisor_masks[1], 4294967295);
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
        assert_eq!(*ht.exponents[1], [1,1,3]);
        assert_eq!(*ht.exponents[2], [2,0,4]);
    }
    #[test]
    fn test_divides() {
        let exps: Vec<Vec<ExpVec>> = vec!(vec!(vec![1,1,3], vec![2,0,3]));
        let mut ht = HashTable::new(&exps);
        let ma  = ht.insert(exps[0][0].clone());
        let dma = ht.divisor_masks[ma as usize];
        let mb  = ht.insert(vec![1,2,3]);
        let dmb = ht.divisor_masks[mb as usize];
        let mc  = ht.insert(vec![3,0,4]);
        let dmc = ht.divisor_masks[mc as usize];
        assert_eq!(ht.divides(ma, dma, mb, !dmb), true);
        assert_eq!(ht.divides(ma, dma, mc, !dmc), false);
        assert_eq!(ht.divides(ma, dma, ma, !dma), true);
    }
    #[test]
    fn test_divides_pairs() {
        let exps: Vec<Vec<ExpVec>> = vec!(vec!(vec![1,1,3], vec![2,0,3]));
        let mut ht = HashTable::new(&exps);
        let ma = ht.insert(exps[0][0].clone());
        let mb = ht.insert(vec![1,2,3]);
        let mc = ht.insert(vec![3,0,4]);
        assert_eq!(ht.divides_pairs(ma, mb), true);
        assert_eq!(ht.divides_pairs(ma, mc), false);
        assert_eq!(ht.divides_pairs(ma, ma), true);
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
        let _diff = ht.get_difference(2,1);
    }
    #[test]
    fn test_get_difference() {
        let exps: Vec<Vec<ExpVec>> = vec!(vec!(
            vec![1,1,3],
            vec![2,1,4]));
        let mut ht = HashTable::new(&exps);
        ht.insert(exps[0][0].clone());
        ht.insert(exps[0][1].clone());
        let diff = ht.get_difference(2,1);
        assert_eq!(ht.exponents[diff as usize], [1,0,1]);
    }
    #[test]
    fn test_get_lcm() {
        let exps: Vec<Vec<ExpVec>> = vec!(vec!(
            vec![1,1,3],
            vec![2,0,3]));
        let mut ht = HashTable::new(&exps);
        ht.insert(exps[0][0].clone());
        ht.insert(exps[0][1].clone());
        let lcm = ht.get_lcm(1,2);
        assert_eq!(ht.exponents[lcm as usize], [2,1,3]);
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
