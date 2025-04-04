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
    pub length           : usize,
    pub nr_in       : usize,
    pub nr_in_ex       : usize,
    pub nr_in_new       : usize,
    divisor_mask_bits_per_variable: usize,
    divisor_mask_variable_range: usize,
    pub exponent_buffer : ExpVec,
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
            nr_in : 0,
            nr_in_ex : 0,
            nr_in_new : 0,
            exponent_buffer : vec![0; initial_exponents[0][0].len()],
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
        -> ExpVec {
        let lm_e = &self.exponents[mon1 as usize];
        let lm_f = &self.exponents[mon2 as usize];
        debug_assert!(lm_e.into_iter().zip(lm_f).all(|(a,b)| *a>=*b));
        return lm_e.into_iter()
            .zip(lm_f)
            .map(|(a,b)| *a-*b)
            .collect();
    }

    pub fn get_lcm(
        &mut self, mon1: HashTableLength, mon2:HashTableLength)
        -> HashTableLength {
        let lm_e = &self.exponents[mon1 as usize];
        let lm_f = &self.exponents[mon2 as usize];
        for i in 0..self.nr_variables {
            self.exponent_buffer[i] = max(lm_e[i], lm_f[i]);
        }
        return self.insert();
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
        self.divisor_bounds = vec!(0; self.divisor_mask_bits_per_variable*self.divisor_mask_variable_range);
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
            for j in 0..self.divisor_mask_bits_per_variable {
                self.divisor_bounds[i+self.divisor_mask_variable_range*j] = ctr;
                ctr += 1;
            }
        }
    }

    // Tests if monomial ma divides monomial mb
    #[inline(always)]
    pub fn divides_pairs(&self, ma: HashTableLength, mb: HashTableLength) -> bool {
        if (self.divisor_masks[ma as usize] & !self.divisor_masks[mb as usize]) != 0 {
            return false;
        }
        if self.degrees[ma as usize] > self.degrees[mb as usize] {
            return false;
        }
        let ea = &self.exponents[ma as usize];
        let eb = &self.exponents[mb as usize];
        return ea.into_iter().zip(eb).all(|(&a,&b)| a <= b);
    }

    #[inline(always)]
    pub fn divides(&self, ma: HashTableLength, dma: DivisorMask, mb: HashTableLength, neg_dmb: DivisorMask) -> bool {
        if (dma & neg_dmb) != 0 {
            return false;
        }
        let ea = &self.exponents[ma as usize];
        let eb = &self.exponents[mb as usize];
        return ea.into_iter().zip(eb).all(|(a,b)| *a <= *b);
    }

    pub fn find_divisor(&mut self, mon: HashTableLength,
            divisor_data_vec: &Vec<(DivisorMask,HashTableLength,BasisLength)>,basis: &Basis)
        -> Option<(BasisLength, ExpVec)> {

        let divisor_data = divisor_data_vec.as_slice();
        let start_idx = self.last_known_divisors[mon as usize];
        let ndmon = !self.divisor_masks[mon as usize];
        if start_idx != 0 {
            if !basis.elements[start_idx as usize].is_redundant {
                return Some((start_idx, self.get_difference(mon, basis.elements[start_idx as usize].monomials[0])))
            }
        }
        for d in divisor_data {
            if self.divides(d.1, d.0, mon, ndmon) {
                self.last_known_divisors[mon as usize] = d.2;
                return Some((d.2, self.get_difference(mon, d.1)));
            }
        }
        return None;
    }

    pub fn generate_multiplied_monomials(&mut self, divisor_idx: BasisLength,
        multiplier: ExpVec, basis: &Basis) -> MonomVec {

        let mons = &basis.elements[divisor_idx as usize].monomials;
        let mut mult_mons: MonomVec = vec!(0; mons.len());
        // let mut exps: ExpVec = vec!(0; self.nr_variables);
        for (idx, m) in mons.iter().enumerate() {
            let e_mon = &self.exponents[*m as usize];
            for i in 0..self.nr_variables {
                self.exponent_buffer[i] = multiplier[i] + e_mon[i];
            }
            // let mut exps: ExpVec = vec!(0; self.nr_variables);
            // for ((e, mu), mo) in exps.iter_mut().zip(e_mult).zip(e_mon) {
            //     *e = mu +mo;
            // }
            mult_mons[idx] = self.insert();
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

    pub fn cmp_monomials_by_index_then_drl(&self, a: HashTableLength, b:HashTableLength) -> Ordering {
        debug_assert!(
            self.exponents[a as usize].len() == self.exponents[b as usize].len());
        
        let ia = &self.indices[a as usize];
        let ib = &self.indices[b as usize];
        if ia != ib {
            if ia > ib { return Ordering::Greater; }
            else { return Ordering::Less; }
        }
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
    pub fn insert(&mut self) -> HashTableLength {
        self.nr_in += 1;
        let div = self.map.len() - 1;
        let h = self.get_hash(&self.exponent_buffer);
        let mut k = h as usize;
        let map_len = self.map.len();
        for  i in 0..map_len {
            k = (k+i) & div;
            let hm = self.map[k];
            if hm == 0 {
                break;
            }
            self.nr_in_ex += 1;
            if self.values[hm as usize]!= h
                || self.exponents[hm as usize] != self.exponent_buffer {
                continue;
            }
            return hm;
        }
        self.nr_in_new += 1;
        let pos = self.exponents.len();
        self.map[k as usize] = pos as HashTableLength;
        self.degrees[pos] = get_degree(&self.exponent_buffer);
        self.divisor_masks[pos] = get_divisor_mask(
            &self.exponent_buffer, &self.divisor_bounds, self.divisor_mask_variable_range);
        // self.last_known_divisors[pos] = 0;
        self.values[pos] = h;
        self.exponents.push(self.exponent_buffer.clone());

        // enlarge AFTER insertion to correctly adjust
        // the just inserted element
        if self.exponents.len() == self.length {
            self.enlarge();
        }

        return pos as HashTableLength;
    }
}

#[inline(always)]
fn get_divisor_mask(exp: &ExpVec, divisor_bounds: &ExpVec, range: usize) -> DivisorMask {
    let mut divisor_mask: DivisorMask = 0;
    let e = &exp[0..range];
    e.into_iter().cycle()
        .zip(divisor_bounds.into_iter().enumerate())
        .for_each(|(e,(i,d))| { if *e >= *d { divisor_mask |= 1 << i;}} );
    return divisor_mask;
}

#[inline(always)]
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
            ht.exponent_buffer = vec!(rng.random::<u16>(),
                            rng.random::<u16>(),
                            rng.random::<u16>(),
                            rng.random::<u16>(),
                            rng.random::<u16>());
            ht.insert();
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
            ht.exponent_buffer = e;
            ht.insert();
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
            ht.exponent_buffer = e;
            ht.insert();
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
            ht.exponent_buffer = e;
            ht.insert();
        }
        ht.exponent_buffer = vec![1 as Exponent,1,1];
        let pos = ht.insert();
        assert_eq!(pos, 1);
    }
    #[test]
    fn test_find_divisor() {
        let fc : Characteristic = 65521;
        let cfs : Vec<CoeffVec> = vec![vec![-2,65523], vec![1, -3]];
        let exps : Vec<Vec<ExpVec>> = vec![vec![vec![0,3], vec![1,1]], vec![vec![0,2], vec![1,1]]];
        let mut hash_table = HashTable::new(&exps);
        let basis = Basis::new::<i32>(&mut hash_table, fc, cfs, exps);
        hash_table.exponent_buffer = vec![0 as Exponent,4];
        let mon1 = hash_table.insert();
        let mut divs: Vec<(DivisorMask,HashTableLength,BasisLength)> = Vec::new();
        for i in 0..basis.elements.len() {
            if basis.elements[i].is_redundant == false {
                divs.push((
                    hash_table.divisor_masks[basis.elements[i].monomials[0] as usize],
                    basis.elements[i].monomials[0],
                    i as BasisLength));
            }
        }
        let multiplier: ExpVec = vec![0, 1];
        assert_eq!(hash_table.find_divisor(mon1, &divs, &basis), Some((1, multiplier)));
        assert_eq!(hash_table.exponents.len(), 5);
        hash_table.exponent_buffer = vec![5 as Exponent,0];
        let mon2 = hash_table.insert();
        assert_eq!(hash_table.find_divisor(mon2, &divs, &basis), None);
    }
    #[test]
    fn test_generate_divisor_bounds() {
        let exps: Vec<Vec<ExpVec>> = vec!(vec!(
            vec![1,1,3],
            vec![2,0,3]));
        let ht = HashTable::new(&exps);
        assert_eq!(ht.divisor_bounds,
            [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5,
             6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 10,
             11, 11, 11, 12, 12, 12, 13, 13, 13, 14, 14, 14,
             15, 15, 15, 16, 16, 16, 17, 17, 17, 18, 18, 18,
             19, 19, 19, 20, 20, 20, 21, 21, 21]);
    }
    #[test]
    fn test_get_divisor_mask() {
        let exps: Vec<Vec<ExpVec>> = vec!(vec!(
            vec![1;65]));
        let mut ht = HashTable::new(&exps);
        ht.exponent_buffer = vec![1;65];
        ht.insert();
        assert_eq!(ht.divisor_masks[1],18446744073709551615);
    }
    #[test]
    fn test_init_hash_table() {
        let exps: Vec<Vec<ExpVec>> = vec!(vec!(
            vec![1,1,3],
            vec![2,0,4]));
        let mut ht = HashTable::new(&exps);
        for e in exps.into_iter().flatten() {
            ht.exponent_buffer = e;
            ht.insert();
        }
        assert_eq!(*ht.exponents[1], [1,1,3]);
        assert_eq!(*ht.exponents[2], [2,0,4]);
    }
    #[test]
    fn test_divides() {
        let exps: Vec<Vec<ExpVec>> = vec!(vec!(vec![1,1,3], vec![2,0,3]));
        let mut ht = HashTable::new(&exps);
        ht.exponent_buffer = exps[0][0].clone();
        let ma  = ht.insert();
        let dma = ht.divisor_masks[ma as usize];
        ht.exponent_buffer = vec![1,2,3];
        let mb  = ht.insert();
        let dmb = ht.divisor_masks[mb as usize];
        ht.exponent_buffer = vec![3,0,4];
        let mc  = ht.insert();
        let dmc = ht.divisor_masks[mc as usize];
        assert_eq!(ht.divides(ma, dma, mb, !dmb), true);
        assert_eq!(ht.divides(ma, dma, mc, !dmc), false);
        assert_eq!(ht.divides(ma, dma, ma, !dma), true);
    }
    #[test]
    fn test_divides_pairs() {
        let exps: Vec<Vec<ExpVec>> = vec!(vec!(vec![1,1,3], vec![2,0,3]));
        let mut ht = HashTable::new(&exps);
        ht.exponent_buffer = exps[0][0].clone();
        let ma  = ht.insert();
        ht.exponent_buffer = vec![1,2,3];
        let mb  = ht.insert();
        ht.exponent_buffer = vec![3,0,4];
        let mc  = ht.insert();
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
        ht.exponent_buffer = exps[0][0].clone();
        ht.insert();
        ht.exponent_buffer = exps[0][1].clone();
        ht.insert();
        let _diff = ht.get_difference(2,1);
    }
    #[test]
    fn test_get_difference() {
        let exps: Vec<Vec<ExpVec>> = vec!(vec!(
            vec![1,1,3],
            vec![2,1,4]));
        let mut ht = HashTable::new(&exps);
        ht.exponent_buffer = exps[0][0].clone();
        ht.insert();
        ht.exponent_buffer = exps[0][1].clone();
        ht.insert();
        assert_eq!(ht.get_difference(2,1), [1,0,1]);
    }
    #[test]
    fn test_get_lcm() {
        let exps: Vec<Vec<ExpVec>> = vec!(vec!(
            vec![1,1,3],
            vec![2,0,3]));
        let mut ht = HashTable::new(&exps);
        ht.exponent_buffer = exps[0][0].clone();
        ht.insert();
        ht.exponent_buffer = exps[0][1].clone();
        ht.insert();
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
        ht.exponent_buffer = exps[0][0].clone();
        let m1 = ht.insert();
        ht.exponent_buffer = exps[0][1].clone();
        let m2 = ht.insert();
        ht.exponent_buffer = exps[0][2].clone();
        let m3 = ht.insert();
        assert_eq!(ht.are_monomials_coprime(m1, m2), true);
        assert_eq!(ht.are_monomials_coprime(m3, m2), false);
    }
}
