use std::collections::HashSet;

use crate::primitives::*;

use crate::pairs::{
    PairSet,
};

use crate::hash_table::{
    HashTable,
};

use crate::basis::{
    Basis,
};

struct Row {
    basis_index: BasisLength,
    columns: HashTableLength,
}

pub struct Matrix {
    known: Vec<Row>,
    unknown: Vec<Row>,
}

impl Matrix {
    pub fn new() -> Matrix {
        let mat = Matrix {
            known : Vec::new(),
            unknown : Vec::new(),
        };

        return mat;
    }

    pub fn get_next_bunch_of_pairs(&mut self, basis: &Basis, pairs: &mut PairSet,
        hash_table: &mut HashTable) {

        let mut next_pairs = pairs.select_pairs_by_minimal_degree(hash_table);
        debug_assert!(next_pairs.len() > 0);

        next_pairs.sort_by(|a,b| hash_table.cmp_monomials_by_drl(a.lcm, b.lcm));
        let mut start = 0;
        let mut gens = HashSet::new();
        while start < next_pairs.len()  {
            let lcm = next_pairs[start].lcm;
            let stop = next_pairs[start..]
                .iter()
                .position(|p| p.lcm != lcm)
                .unwrap_or(next_pairs.len());

            for i in start..stop {
                gens.insert(next_pairs[i].generators.0);
                gens.insert(next_pairs[i].generators.1);
            }
            for g in &gens {
                let mons = &basis.elements[*g].monomials;
                let multiplier = hash_table.get_difference(lcm, mons[0]);
                for m in mons.iter() {

                // known.appened(Row { });
                }
            }
            start = stop;
            gens.clear();
        }

    }
}
