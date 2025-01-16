use crate::types::*;
use crate::meta_data::{
    MetaData,
};
use crate::hash_table::{
    HashTable,
};

pub struct Element {
    coefficients: CoeffVec,
    monomials: MonomVec,
    is_redundant: bool,
}

pub struct Basis<'a> {
    previous_length: BasisLength,
    is_constant: bool,
    maximum_total_degree: Degree,
    elements: Vec<Element>,
    hash_table: HashTable<'a>,
    meta_data: MetaData,
}

impl<'a> Basis<'a> {
    pub fn new<T>(
        field_characteristic: Coefficient,
        coefficients: &'a Vec<CoeffVec>,
        exponents: &'a Vec<Vec<ExpVec>>) -> Basis<'a> {
        debug_assert!(coefficients.len() > 0);
        debug_assert!(coefficients.len() == exponents.len());
        let mut basis = Basis {
            previous_length      : 0,
            is_constant          : false,
            maximum_total_degree : 0,
            hash_table           : HashTable::new(exponents),
            meta_data            : MetaData {
                characteristic : field_characteristic,
                nr_pairs_reduced : 0,
                nr_redundant_elements : 0,
                nr_input_generators : coefficients.len(),
            },
            elements             : Vec::new(),
        };
        basis.add_initial_elements(coefficients, exponents);
        basis.sort_elements_by_drl();
        return basis;
    }

    fn add_initial_elements(&mut self, cfs: &Vec<CoeffVec>, exps: &'a Vec<Vec<ExpVec>>) {
        let hash_table = &mut self.hash_table;
        let characteristic = self.meta_data.characteristic;
        for (c,e) in cfs.iter().zip(exps) {
            // c.iter_mut().for_each(|a| *a =* a % characteristic);
            // let mons = e.iter().map(|a| hash_table.insert(a)).collect();
            self.elements.push(Element {
                coefficients: c.iter().map(|a| a % characteristic).collect(),
                monomials: e.iter().map(|a| hash_table.insert(a)).collect(),
                is_redundant: false});
        }
    }

    fn sort_elements_by_drl(&mut self) {
        return self.elements.sort_by(
            |a,b| self.hash_table.cmp_monomials_by_drl(
                a.monomials[0], b.monomials[0]));
    }
}
