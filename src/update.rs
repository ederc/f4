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
pub struct LcmGenerators<'a> {
    lcm: HashTableLength,
    generators: BTreeSet<&'a Element>,
}

pub type LcmVec<'a> = Vec<LcmGenerators<'a>>;

pub fn update_lcms(lcm_vec: &mut LcmVec, basis: &mut Basis, hash_table: &mut HashTable) {

    while (basis.previous_length as usize) < basis.elements.len() {
        let e = &basis.elements[basis.previous_length as usize];
        for f in &basis.elements[0..(basis.previous_length as usize)] {
            let lcm = hash_table.get_lcm(e.monomials[0], f.monomials[0]);
            match lcm_vec.iter().find(|x| x.lcm == lcm) {
                Some(olcm) => {
                    olcm.generators.insert(&e);
                    olcm.generators.insert(&f);
                },
                None => {
                    lcm_vec.push(LcmGenerators {
                        lcm: lcm,
                        generators: BTreeSet::new(),
                    }),
                    lcm_vec[lcm_vec.len()-1].generators.insert(&e),
                    lcm_vec[lcm_vec.len()-1].generators.insert(&f),
                },
            }
            println!("{}", lcm_vec.len());
            // lcm_vec[idx].generators.insert(&e);
            // lcm_vec[idx].generators.insert(&f);
        }
        basis.previous_length += 1;
    }
}
