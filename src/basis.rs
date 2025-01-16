use crate::types::*;
use crate::hash_table::{
    HashTable,
};

pub struct Element<'a> {
    coefficients: &'a CoeffVec,
    monomials: &'a MonomVec,
    is_redundant: bool,
}

pub struct Basis<'a> {
    previous_length: BasisLength,
    is_constant: bool,
    maximum_total_degree: Degree,
    elements: Vec<Element<'a>>,
    hash_table: HashTable<'a>,
}

impl<'a> Basis<'a> {
    pub fn new<T>(
        coefficients: &'a Vec<CoeffVec>,
        exponents: &'a Vec<Vec<ExpVec>>) -> Basis<'a> {
        debug_assert!(coefficients.len() > 0);
        debug_assert!(coefficients.len() == exponents.len());
        let mut basis = Basis {
            previous_length      : 0,
            is_constant          : false,
            maximum_total_degree : 0,
            hash_table           : HashTable::new(exponents),
            elements             : Vec::new(),
        };
        basis.add_initial_elements(coefficients, exponents);
        basis.sort_elements_by_drl();
        return basis;
    }

    fn add_initial_elements(&mut self, cfs: & Vec<CoeffVec>, exps: & Vec<Vec<ExpVec>>) {
        let mut hash_table = self.hash_table;
        let prime = self.meta_data.prime;
        for (c,e) in cfs.iter().zip(exps) {
            c.iter().for_each(|a| *a =* a % prime);
            let mons = e.iter().map(|a| hash_table.insert(a)).collect();
            self.elements.push(Element {
                coefficients: &c,
                monomials: &mons,
                is_redundant: false});
        }
    }

    fn sort_elements_by_drl(&mut self) {
        return self.elements.sort_by(
            |a,b| self.hash_table.cmp_monomials_by_drl(
                a.monomials[0], b.monomials[0]));
    }
}
