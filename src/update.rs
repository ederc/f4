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
    pub fn update(
        &mut self,
        basis: &Basis,
        hash_table: &mut HashTable
    ) {
    }
}
