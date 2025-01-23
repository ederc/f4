use crate::primitives::*;

use crate::basis::{
    Basis,
};

use crate::hash_table::{
    HashTable,
};

use crate::meta_data::{
    MetaData,
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

pub type PairVec = Vec<Pair>;

impl PairVec {
    pub fn update(
        &mut self,
        basis: &Basis,
        hash_table: &mut HashTable,
        meta_data: &MetaData
    ) {
    }
}
