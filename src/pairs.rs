use crate::primitives::*;

use crate::basis::{
    Basis,
};

use crate::hash_table::{
    HashTable,
};

#[derive(PartialEq, Debug)]
enum Criterion {
    Keep,
    Chain,
    Product,
}

#[derive(PartialEq, Debug)]
pub struct Pair {
    lcm: HashTableLength,
    generators: (BasisLength, BasisLength),
    criterion: Criterion,
}

type PairVec = Vec<Pair>;

#[derive(Debug)]
pub struct PairSet {
    pub list: PairVec,
    pub nr_pairs_reduced: usize,
    nr_product_criterion_applied: usize,
    nr_chain_criterion_applied: usize,
}

impl PairSet {
    pub fn new() -> PairSet {
        return PairSet {
            list                         : Vec::new(),
            nr_pairs_reduced             : 0,
            nr_product_criterion_applied : 0,
            nr_chain_criterion_applied   : 0,
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
            let mut new_pairs: PairVec = basis.elements[..i]
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
                } else if new_pairs[i].criterion == Criterion::Keep && i > 0 {
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

            // no sorting here, we sort just before extracting
            // the pairs in symbolic preprocessing
            self.list.append(&mut new_pairs);

            // bookkeeping of applied criteria
            self.nr_product_criterion_applied += self.list.iter().filter(
                |p| p.criterion == Criterion::Product).count();
            self.nr_chain_criterion_applied += self.list.iter().filter(
                |p| p.criterion == Criterion::Chain).count();
            // remove useless pairs
            self.list.retain(|p| p.criterion == Criterion::Keep);
        }
    }

    fn select_pairs_by_minimal_degree(
        &mut self,
        basis: &Basis,
        hash_table: &HashTable,
    ) -> PairVec {

        // sort pair list by degree
        self.list.sort_by(|a,b| hash_table.cmp_monomials_by_degree(a.lcm, b.lcm));
        // get minimal degree pairs
        let min_degree = hash_table.monomials[self.list[0].lcm].degree;
        let idx = self.list.iter()
            .position(|p| hash_table.monomials[p.lcm].degree > min_degree)
            .unwrap();
        let min_degree_pairs = self.list.split_off(idx);
        // bookkeeping
        self.nr_pairs_reduced += min_degree_pairs.len();

        return min_degree_pairs;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_update() {
        let fc : Characteristic = 65521;
        let cfs : Vec<CoeffVec> = vec![vec![-2,65523], vec![1, -3], vec![1, -1], vec![1, 1]];
        let exps : Vec<Vec<ExpVec>> = vec![vec![vec![0,3,1], vec![1,1,0]], vec![vec![0,2,0], vec![1,1,0]], vec![vec![0,0,2], vec![1,0,0]], vec![vec![0,0,1], vec![0,0,0]]];
        let mut hash_table = HashTable::new(&exps);
        let basis = Basis::new::<i32>(&mut hash_table, fc, cfs, exps);

        let mut pairs = PairSet::new();
        pairs.update(&basis, &mut hash_table);

        assert_eq!(pairs.list.len(), 2);
        assert_eq!(pairs.nr_product_criterion_applied, 2);
        assert_eq!(pairs.nr_chain_criterion_applied, 2);
        assert_eq!(pairs.list[0], Pair { lcm: 3, generators: (1, 0), criterion: Criterion::Keep } );
        assert_eq!(pairs.list[1], Pair { lcm: 0, generators: (3, 0), criterion: Criterion::Keep } );
    }
    #[test]
    fn test_select_pairs_by_minimal_degree() {
        let fc : Characteristic = 65521;
        let cfs : Vec<CoeffVec> = vec![vec![-2,65523], vec![1, -3], vec![1, -1], vec![1, 1]];
        let exps : Vec<Vec<ExpVec>> = vec![vec![vec![0,3,1], vec![1,1,0]], vec![vec![0,2,0], vec![1,1,0]], vec![vec![0,0,2], vec![1,0,0]], vec![vec![0,0,1], vec![0,0,0]]];
        let mut hash_table = HashTable::new(&exps);
        let basis = Basis::new::<i32>(&mut hash_table, fc, cfs, exps);

        let mut pairs = PairSet::new();
        pairs.update(&basis, &mut hash_table);
        let min_degree_pairs = pairs.select_pairs_by_minimal_degree(&basis, &hash_table);

        assert_eq!(min_degree_pairs.len(), 1);
        assert_eq!(min_degree_pairs[0], Pair { lcm: 0, generators: (3, 0), criterion: Criterion::Keep } );
    }
}
