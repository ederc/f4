use crate::types::*;
use std::cmp;
// 2^17
const INITIAL_HASH_TABLE_SIZE: usize = 131072;

struct Monomial<'a> {
    degree: Degree,
    divisor_mask: DivisorMask,
    exponents: &'a Vec<Exponent>,
    last_known_divisor: BasisLength,
}

pub struct HashTable<'a> {
    monomials     : Vec<Monomial<'a>>,
    random_seed   : Vec<HashValue>,
    values        : Vec<HashValue>,
    map           : Vec<HashTableLength>,
    divisor_bounds: Vec<Exponent>,
    indices       : Vec<HashValue>,
    nr_variables  : usize,
}

impl<'a> HashTable<'a> {
    pub fn new(initial_exponents: &'a Vec<Vec<Vec<Exponent>>>) -> HashTable<'a> {
        debug_assert!(initial_exponents.len() > 0);
        debug_assert!(initial_exponents[0].len() > 0);
        debug_assert!(initial_exponents[0][0].len() > 0);
        let mut ht = HashTable {
            monomials      : Vec::new(),
            random_seed    : Vec::new(),
            values         : Vec::new(),
            map            : vec![0; INITIAL_HASH_TABLE_SIZE],
            divisor_bounds : Vec::new(),
            indices        : vec![0; INITIAL_HASH_TABLE_SIZE],
            nr_variables   : initial_exponents[0][0].len(),
        };
        ht.generate_random_seed(ht.nr_variables);
        ht.generate_divisor_bounds(&initial_exponents);
        for e in initial_exponents.into_iter().flatten() {
            ht.insert(e);
        }

        return ht;
    }

    fn update_seed(&self, mut seed: HashValue) -> HashValue {
        seed ^= seed << 13;
        seed ^= seed >> 17;
        seed ^= seed << 5;
        return seed;
    }

    fn generate_random_seed(& mut self, nr_variables: usize) {
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
    fn generate_divisor_bounds(& mut self, initial_exponents: &Vec<Vec<Vec<Exponent>>>) {
        let length_divmask = cmp::min(
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

    fn get_divisor_mask(& mut self, exp: &Vec<Exponent>) -> DivisorMask {
        let divisor_bounds = &self.divisor_bounds;
        let mut divisor_mask = 0usize;
        for i in 0..exp.len() {
            if exp[i] >= divisor_bounds[i] {
                divisor_mask |= 1 << i;
            }
        }
        return divisor_mask;
    }

    fn get_hash(&self, exp: &Vec<Exponent>) -> HashValue {
        let mut h: HashValue = 0;
        for i in 0..exp.len() {
            h = h.wrapping_add((exp[i] as HashValue).wrapping_mul(self.random_seed[i]));
        }
        return h;
    }

    fn insert(& mut self, exp: &'a Vec<Exponent>) -> HashValue {
        let div = self.map.len() - 1;
        let h = self.get_hash(exp);
        let mut k = h;
        let mut i = 0;
        while i < self.map.len() {
            k = (k+i) & div;
            let hm = self.map[k];
            if hm != 0 {
                break;
            }
            if self.values[hm] != h {
                i += 1;
                continue;
            }
            let eh = self.monomials[hm].exponents;
            if *eh != *exp {
                i += 1;
                continue;
            }
            return hm;
        }
        let pos = self.monomials.len();
        self.map[k] = pos;
        println!("exponent {:?}", exp);
        let monomial = Monomial {
            degree: exp.iter().sum(),
            divisor_mask: self.get_divisor_mask(&exp),
            exponents: exp,
            last_known_divisor: 0,
        };
        self.monomials.push(monomial);
        self.values.push(h);

        return pos;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_random_seed() {
        let exp: Vec<Vec<Vec<Exponent>>> = vec!(vec!(vec!(1,1,1,1,1)));
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
    fn test_insert() {
        let exp: Vec<Vec<Vec<Exponent>>> = vec!(vec!(vec![1,1,1]));
        let mut ht = HashTable::new(&exp);
        let pos = ht.insert(&vec![1 as Exponent,1,1]);
        assert_eq!(pos, 0);
    }
    #[test]
    fn test_generate_divisor_bounds() {
        let exps: Vec<Vec<Vec<Exponent>>> = vec!(vec!(
            vec![1,1,3],
            vec![2,0,3]));
        let ht = HashTable::new(&exps);
        assert_eq!(ht.divisor_bounds, [1,1,0]);
    }
    #[test]
    fn test_init_hash_table() {
        let exps: Vec<Vec<Vec<Exponent>>> = vec!(vec!(
            vec![1,1,3],
            vec![2,0,4]));
        let mut ht = HashTable::new(&exps);
        assert_eq!(*ht.monomials[0].exponents, [1,1,3]);
        assert_eq!(*ht.monomials[1].exponents, [2,0,4]);
    }

        




}
