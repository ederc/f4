use std::collections::{
    BTreeSet,
};

use crate::primitives::*;

use crate::basis::{
    Basis,
    Element,
};

use crate::hash_table::{
    HashTable,
};

#[derive(Debug)]
pub struct LcmSet {
    lcm: HashTableLength,
    generators: BTreeSet<BasisLength>,
}

pub type LcmVec = Vec<LcmSet>;

pub struct LcmData {
    cur: LcmVec, // lcm data to be handled
    old: LcmVec, // lcm data already handled
}

impl LcmData {
    pub fn new() -> LcmData {
        let lcms = LcmData {
            cur: Vec::new(),
            old: Vec::new(),
        };
        return lcms;
    }

    pub fn update(&mut self, basis: &Basis, hash_table: &mut HashTable) {

        let mut e_idx = basis.previous_length as usize;
        while e_idx < basis.elements.len() {
            let e = &basis.elements[e_idx];
            for (f_idx, f) in basis.elements[0..e_idx].iter().enumerate() {
                let lcm = hash_table.get_lcm(e.monomials[0], f.monomials[0]);
                let pos = self.cur.iter()
                    .position(|x| x.lcm == lcm)
                    .unwrap_or(self.cur.len());
                println!("pos {}", pos);

                if pos == self.cur.len() {
                    self.cur.push(LcmSet { lcm: lcm, generators: BTreeSet::new(), });
                }
                self.cur[pos].generators.insert(e_idx as BasisLength);
                self.cur[pos].generators.insert(f_idx as BasisLength);
                println!("{:?}", self.cur[pos].generators);
                // lcm_vec[idx].generators.insert(&e);
                // lcm_vec[idx].generators.insert(&f);
            }
            e_idx += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_update() {
        let fc : Characteristic = 65521;
        let cfs : Vec<CoeffVec> = vec![vec![-2,65523], vec![1, -3]];
        let exps : Vec<Vec<ExpVec>> = vec![vec![vec![0,3], vec![1,1]], vec![vec![0,2], vec![1,1]]];
        let mut hash_table = HashTable::new(&exps);
        let basis = Basis::new::<i32>(&mut hash_table, fc, cfs, exps);
        let mut lcms = LcmData::new();

        lcms.update(&basis, &mut hash_table);

        assert_eq!(lcms.old.len(), 0);
        assert_eq!(lcms.cur.len(), 1);
        assert_eq!(lcms.cur[0].generators, BTreeSet::from([0,1]));
        assert_eq!(lcms.cur[0].lcm, 3);
        assert_eq!(hash_table.exponents[lcms.cur[0].lcm as usize], [1,3]);
    }
}
