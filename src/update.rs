use crate::primitives::*;

use crate::basis::{
    Basis,
};

use crate::hash_table::{
    HashTable,
};

enum Criterion {
    Keep,
    Chain,
    Product,
}

pub struct Pair {
    lcm: HashTableLength,
    generators: (BasisLength, BasisLength),
    criterion: Criterion,
}

type PairVec = Vec<Pair>;

pub struct PairSet {
    list: PairVec,
    nr_pairs_reduced: usize,
    nr_criteria_applied: usize,
}

impl PairSet {
    pub fn new() -> PairSet {
        return PairSet {
            list               : Vec::new(),
            nr_pairs_reduced   : 0,
            nr_criteria_applied: 0,
        };
    }

    pub fn is_empty(&self) -> bool {
        return self.list.len() == 0;
    }

    pub fn update(
        &mut self,
        basis: &Basis,
        hash_table: &mut HashTable
    ) {
        for (i,e) in basis.elements[basis.previous_length..].iter().enumerate() {
            let mut new_pairs = basis.elements[..basis.previous_length].iter().enumerate().map(|(j,f)|
                Pair {
                    lcm: hash_table.get_lcm(e.monomials[0], f.monomials[0]),
                    generators: (i, j),
                    criterion: if hash_table.are_monomials_coprime(
                        e.monomials[0], f.monomials[0])
                        { Criterion::Product } else { Criterion::Keep },
                }).collect();
            self.list.append(&mut new_pairs);
        }
    }
}
