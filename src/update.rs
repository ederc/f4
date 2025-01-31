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
        for (i,e) in basis.elements[basis.previous_length..]
        .iter().enumerate() {
            // generate new pairs with basis element e
            let mut new_pairs: PairVec = basis.elements[..basis.previous_length]
                .iter().enumerate().map(|(j,f)|
                    Pair {
                        lcm: hash_table.get_lcm(e.monomials[0], f.monomials[0]),
                        generators: (i, j),
                        criterion: if hash_table.are_monomials_coprime(
                            e.monomials[0], f.monomials[0])
                        { Criterion::Product } else if f.is_redundant { Criterion::Chain}
                        else { Criterion::Keep },
                    }).collect();

            // Gebauer-Möller: testing old pairs
            for p in self.list {
                let deg_p = hash_table.monomials[p.lcm].degree;
                if deg_p > hash_table.monomials[new_pairs[p.generators.0].lcm].degree
                    && deg_p > hash_table.monomials[new_pairs[p.generators.1].lcm].degree
                    && hash_table.divides(e.monomials[0], p.lcm) {
                        p.criterion = Criterion::Chain;
                }
            }
            // no sorting here, we sort just before extracting
            // the pairs in symbolic preprocessing
            self.list.append(&mut new_pairs);
        }
    }
}
