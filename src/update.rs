use crate::primitives::*;

use crate::basis::{
    Basis,
};

use crate::hash_table::{
    HashTable,
};

#[derive(PartialEq)]
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
            for p in &mut self.list {
                let deg_p = hash_table.monomials[p.lcm].degree;
                if deg_p > hash_table.monomials[new_pairs[p.generators.0].lcm].degree
                    && deg_p > hash_table.monomials[new_pairs[p.generators.1].lcm].degree
                    && hash_table.divides(e.monomials[0], p.lcm) {
                        p.criterion = Criterion::Chain;
                }
            }
            // sort new pairs for following Gebauer-Möller steps
            new_pairs.sort_by(|a,b| hash_table.cmp_monomials_by_drl(a.lcm, b.lcm));

            // Gebauer-Möller: remove real multiples from new pairs
            for i in 0..new_pairs.len() {
                if new_pairs[i].criterion == Criterion::Keep {
                    for j in 0..i {
                        if new_pairs[j].criterion != Criterion::Chain
                            && new_pairs[j].lcm != new_pairs[i].lcm
                            && hash_table.divides(new_pairs[j].lcm, new_pairs[i].lcm) {
                                new_pairs[i].criterion = Criterion::Chain;
                                break;
                        }
                    }
                }
            }

            // Gebauer-Möller: remove same lcm pairs from new pairs
            for i in 0..new_pairs.len() {
                if new_pairs[i].criterion == Criterion::Product {
                    for j in 0..new_pairs.len() {
                        if i != j && new_pairs[i].lcm == new_pairs[j].lcm {
                            new_pairs[j].criterion = Criterion::Chain;
                        }
                    }
                } else if new_pairs[i].criterion == Criterion::Keep {
                    for j in i-1..=0 {
                        if new_pairs[j].lcm != new_pairs[i].lcm {
                            break;
                        } else if new_pairs[j].criterion == Criterion::Keep {
                            new_pairs[i].criterion = Criterion::Chain;
                            break;
                        }
                    }
                }
            }

            // remove useless new pairs
            new_pairs.retain(|p| p.criterion == Criterion::Keep);

            // no sorting here, we sort just before extracting
            // the pairs in symbolic preprocessing
            self.list.append(&mut new_pairs);

        }
    }
}
