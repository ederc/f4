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
        basis.sort_terms_by_drl();
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
    fn sort_terms_by_drl(&mut self) {
        for e in &mut self.elements {
            let mut zipped: Vec<_> = e.monomials.clone().into_iter().zip(e.coefficients.clone()).collect();
            zipped.sort_unstable_by(|(a, _), (b, _)| self.hash_table
                 .cmp_monomials_by_drl(*a, *b));
            (e.monomials, e.coefficients) = zipped.into_iter().unzip();
        }
    }

    fn sort_elements_by_drl(&mut self) {
        return self.elements.sort_by(
            |a,b| self.hash_table.cmp_monomials_by_drl(
                a.monomials[0], b.monomials[0]));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_initialize_basis() {
        let fc : Characteristic = 65521;
        let cfs : Vec<CoeffVec> = vec![vec![-2,65523], vec![1, -3]];
        let exps : Vec<Vec<ExpVec>> = vec![vec![vec![0,3], vec![1,1]], vec![vec![0,2], vec![1,1]]];
        let basis = Basis::new::<i32>(fc, &cfs, &exps);
        assert_eq!(basis.elements[0].coefficients, [1,-3]);
        assert_eq!(basis.elements[0].monomials, [2,1]);
        assert_eq!(basis.elements[1].coefficients, [2,-2]);
        assert_eq!(basis.elements[1].monomials, [1,0]);
    }
}
