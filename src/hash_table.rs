type Exponent = i32;
type HashValue = u32;

// 2^17
const INITIAL_HASH_TABLE_SIZE: usize = 131072;

pub struct HashTable {
    exponents     : Vec<Vec<Exponent>>,
    values        : Vec<HashValue>,
    random_seed   : Vec<HashValue>,
    map           : Vec<HashValue>,
    indices       : Vec<HashValue>,
    divisor_masks : Vec<HashValue>,
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
        ht.exponents.push(vec![-1,1]);
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
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_random_seed() {
        let ht = HashTable::new(5);
        assert_eq!(ht.random_seed,
            [723471715, 2497366906, 2064144800, 2008045182, 3532304609]);
    }
}
