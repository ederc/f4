use crate::primitives::*;
use std::cmp:: {
    Ordering,
    min,
    max,
};
// 2^17
const INITIAL_HASH_TABLE_SIZE: usize = 131072;

struct Monomial {
    degree: Degree,
    divisor_mask: DivisorMask,
    exponents: ExpVec,
    last_known_divisor: BasisLength,
}

pub struct HashTable {
    monomials     : Vec<Monomial>,
    random_seed   : Vec<HashValue>,
    values        : Vec<HashValue>,
    map           : Vec<HashTableLength>,
    divisor_bounds: ExpVec,
    indices       : Vec<HashValue>,
    nr_variables  : usize,
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
            map            : vec![usize::MAX; INITIAL_HASH_TABLE_SIZE],
            divisor_bounds : Vec::new(),
            indices        : vec![0; INITIAL_HASH_TABLE_SIZE],
            nr_variables   : initial_exponents[0][0].len(),
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

    pub fn get_lcm(
        &mut self, mon1: HashTableLength, mon2:HashTableLength)
        -> HashTableLength {
        let lm_e = &self.monomials[mon1].exponents;
        let lm_f = &self.monomials[mon2].exponents;
        return self.insert(lm_e.into_iter()
            .zip(lm_f)
            .map(|(a,b)| max(*a,*b))
            .collect());
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
        for i in 0..exp.len() {
            if exp[i] >= divisor_bounds[i] {
                divisor_mask |= 1 << i;
            }
        }
        return divisor_mask;
    }

    fn get_hash(&self, exp: &ExpVec) -> HashValue {
        let mut h: HashValue = 0;
        for i in 0..exp.len() {
            h = h.wrapping_add((exp[i] as HashValue).wrapping_mul(self.random_seed[i]));
        }
        return h;
    }

    pub fn cmp_monomials_by_drl(&self, a: HashTableLength, b:HashTableLength) -> Ordering {
        debug_assert!(
            self.monomials[a].exponents.len() == self.monomials[b].exponents.len());
        let ma = &self.monomials[a];
        let mb = &self.monomials[b];
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

    pub fn insert(&mut self, exp: ExpVec) -> HashTableLength {
        let div = self.map.len() - 1;
        let h = self.get_hash(&exp);
        let mut k = h;
        let mut i = 0;
        while i < self.map.len() {
            k = (k+i) & div;
            let hm = self.map[k];
            if hm == usize::MAX {
                break;
            }
            if self.values[hm] != h {
                i += 1;
                continue;
            }
            let eh = &self.monomials[hm].exponents;
            if *eh != exp {
                i += 1;
                continue;
            }
            return hm;
        }
        let pos = self.monomials.len();
        self.map[k] = pos;
        let monomial = Monomial {
            degree: exp.iter().sum(),
            divisor_mask: self.get_divisor_mask(&exp),
            exponents: exp,
            last_known_divisor: 0,
        };
        self.monomials.push(monomial);
        self.values[pos] = h;

        return pos;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_random_seed() {
        let exp: Vec<Vec<ExpVec>> = vec!(vec!(vec!(1,1,1,1,1)));
        let ht = HashTable::new(&exp);
        if std::mem::size_of::<usize>() == 4 {
        assert_eq!(ht.random_seed,
            [723471715, 2497366906, 2064144800, 2008045182, 3532304609]);
        }
        if std::mem::size_of::<usize>() == 8 {
        assert_eq!(ht.random_seed,
            [660888219700579, 3396719463693796860,
            17326311066685913516, 2586175631380707723, 16544630075375549064]);
        }
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
    fn test_generate_divisor_bounds() {
        let exps: Vec<Vec<ExpVec>> = vec!(vec!(
            vec![1,1,3],
            vec![2,0,3]));
        let ht = HashTable::new(&exps);
        assert_eq!(ht.divisor_bounds, [1,1,0]);
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
    fn test_get_lcm() {
        let exps: Vec<Vec<ExpVec>> = vec!(vec!(
            vec![1,1,3],
            vec![2,0,3]));
        let mut ht = HashTable::new(&exps);
        ht.insert(exps[0][0].clone());
        ht.insert(exps[0][1].clone());
        let lcm = ht.get_lcm(0,1);
        assert_eq!(ht.monomials[lcm].exponents, [2,1,3]);
    }
}
