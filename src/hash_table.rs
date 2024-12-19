use crate::types::*;
// 2^17
const INITIAL_HASH_TABLE_SIZE: usize = 131072;

pub struct HashTable {
    exponents     : Vec<Vec<Exponent>>,
    values        : Vec<HashValue>,
    random_seed   : Vec<HashValue>,
    map           : Vec<HashTableLength>,
    indices       : Vec<HashValue>,
    divisor_masks : Vec<DivisorMask>,
}

impl HashTable {
    pub fn new(nr_variables: usize) -> HashTable {
        let mut ht = HashTable {
            exponents     : Vec::new(),
            values        : Vec::new(),
            random_seed   : Vec::new(),
            map           : vec![0; INITIAL_HASH_TABLE_SIZE],
            indices       : vec![0; INITIAL_HASH_TABLE_SIZE],
            divisor_masks : vec![0; INITIAL_HASH_TABLE_SIZE],
        };
        // initialize entry at index 0 with useless data
        ht.exponents.push(vec![u32::MAX]);
        ht.values.push(0);
        ht.generate_random_seed(nr_variables);

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

    fn get_hash(&self, exp: &Vec<Exponent>) -> HashValue {
        let mut h: HashValue = 0;
        for i in 0..exp.len() {
            h = h.wrapping_add(exp[i] as HashValue * self.random_seed[i]);
        }
        return h;
    }

    fn insert(& mut self, exp: Vec<Exponent>) -> HashValue {
        let div = self.map.len() - 1;
        let h = self.get_hash(&exp);
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
            let eh = &self.exponents[hm];
            if *eh != exp {
                i += 1;
                continue;
            }
            return hm;
        }
        let pos = self.exponents.len();
        self.map[k] = pos;
        self.exponents.push(exp);
        self.values.push(h);

        return pos;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_random_seed() {
        let ht = HashTable::new(5);
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
        let mut ht = HashTable::new(3);
        let pos = ht.insert(vec![1,1,1]);
        assert_eq!(pos, 1);
    }
}
