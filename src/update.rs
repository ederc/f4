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

pub fn update_lcms(lcm_vec: &mut LcmVec, basis: &mut Basis, hash_table: &mut HashTable) {

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
