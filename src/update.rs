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
pub struct LcmData {
    lcm: HashTableLength,
    generators: BTreeSet<BasisLength>,
}

pub type LcmVec = Vec<LcmData>;

pub fn update_lcms(lcm_vec: &mut LcmVec, basis: &Basis, hash_table: &mut HashTable) {

    let mut e_idx = basis.previous_length as usize;
    while e_idx < basis.elements.len() {
        let e = &basis.elements[e_idx];
        for (f_idx, f) in basis.elements[0..e_idx].iter().enumerate() {
            let lcm = hash_table.get_lcm(e.monomials[0], f.monomials[0]);
            let pos = lcm_vec.iter().position(|x| x.lcm == lcm).unwrap_or(lcm_vec.len());
            println!("pos {}", pos);

            if pos == lcm_vec.len() {
                lcm_vec.push(LcmData { lcm: lcm, generators: BTreeSet::new(), });
            }
            lcm_vec[pos].generators.insert(e_idx as BasisLength);
            lcm_vec[pos].generators.insert(f_idx as BasisLength);
            println!("{:?}", lcm_vec[pos].generators);
            // lcm_vec[idx].generators.insert(&e);
            // lcm_vec[idx].generators.insert(&f);
        }
        e_idx += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_update_lcms() {
        let fc : Characteristic = 65521;
        let cfs : Vec<CoeffVec> = vec![vec![-2,65523], vec![1, -3]];
        let exps : Vec<Vec<ExpVec>> = vec![vec![vec![0,3], vec![1,1]], vec![vec![0,2], vec![1,1]]];
        let mut hash_table = HashTable::new(&exps);
        let basis = Basis::new::<i32>(&mut hash_table, fc, cfs, exps);
        let mut lcms = LcmVec::new();

        update_lcms(&mut lcms, &basis, &mut hash_table);

        assert_eq!(lcms.len(), 1);
        assert_eq!(lcms[0].generators, BTreeSet::from([0,1]));
        assert_eq!(lcms[0].lcm, 3);
        assert_eq!(hash_table.exponents[lcms[0].lcm as usize], [1,3]);
    }
}
