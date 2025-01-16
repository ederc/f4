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
        let mut basis = Basis {
            previous_length      : 0,
            is_constant          : false,
            maximum_total_degree : 0,
            hash_table           : HashTable::new(exponents),
            elements             : Vec::new(),
        };
        basis.sort_elements_by_drl();
        return basis;
    }
    
    fn sort_elements_by_drl(&mut self) {
        return self.elements.sort_by(
            |a,b| self.hash_table.cmp_monomials_by_drl(
                a.monomials[0], b.monomials[0]));
    }
}
