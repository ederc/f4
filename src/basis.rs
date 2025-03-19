use crate::primitives::*;

use crate::arithmetic::i32::modular_inverse;

use std::cmp:: {
    max,
};

use crate::hash_table::{
    HashTable,
};

#[derive(Clone)]
pub struct Element {
    pub coefficients: CoeffVec,
    pub monomials: MonomVec,
    pub is_redundant: bool,
}

pub struct Basis {
    pub characteristic: Characteristic,
    pub previous_length: BasisLength,
    is_constant: bool,
    maximum_total_degree: Degree,
    pub elements: Vec<Element>,
    nr_redundant_elements: usize,
    nr_input_generators: usize,
}

impl Basis {
    pub fn new<T>(
        hash_table: &mut HashTable,
        characteristic: Characteristic,
        coefficients: Vec<CoeffVec>,
        exponents: Vec<Vec<ExpVec>>) -> Basis {
        debug_assert!(coefficients.len() > 0);
        debug_assert!(coefficients.len() == exponents.len());
        let mut basis = Basis {
            characteristic       : characteristic,
            previous_length      : 0,
            is_constant          : false,
            maximum_total_degree : 0,
            elements             : Vec::new(),
            nr_redundant_elements: 0,
            nr_input_generators  : coefficients.len(),
        };
        basis.add_initial_elements(hash_table, coefficients, exponents);
        basis.sort_terms_by_drl(hash_table);
        basis.sort_elements_by_drl(hash_table);
        basis.normalize_elements();
        return basis;
    }
    
    fn normalize_elements(&mut self) {
        let characteristic:DenseRowCoefficient = self.characteristic.into();
        for e in &mut self.elements {
            // make all coefficients >= 0
            e.coefficients.iter_mut()
                .for_each(|c| { while *c < 0 { *c += self.characteristic; } });
            // normalize if needed
            if e.coefficients[0] != 1 {
            let inv = modular_inverse(e.coefficients[0],
                characteristic as Characteristic) as DenseRowCoefficient;
            e.coefficients.iter_mut()
                .for_each(|c|
                    *c = ((inv * *c as DenseRowCoefficient) % characteristic)
                    as Coefficient);
            }
        }
    }

    fn add_initial_elements(
        &mut self,
        hash_table: &mut HashTable,
        cfs: Vec<CoeffVec>, exps: Vec<Vec<ExpVec>>) {
        for (c,e) in cfs.iter().zip(exps) {
            self.elements.push(Element {
                coefficients: c.iter().map(|a| *a % self.characteristic as Coefficient).collect(),
                monomials: e.iter().map(|a| hash_table.insert(a.to_vec())).collect(),
                is_redundant: false});
        }
    }
    fn sort_terms_by_drl(&mut self, hash_table: &HashTable) {
        for e in &mut self.elements {
            let mut zipped: Vec<_> = e.monomials.clone().into_iter()
                .zip(e.coefficients.clone()).collect();
            zipped.sort_unstable_by(|(a, _), (b, _)| hash_table
                .cmp_monomials_by_drl(*b, *a));
                (e.monomials, e.coefficients) = zipped.into_iter().unzip();
        }
    }

    fn sort_elements_by_drl(&mut self, hash_table: &HashTable) {
        return self.elements.sort_by(
            |a,b| hash_table.cmp_monomials_by_drl(
                a.monomials[0], b.monomials[0]));
    }

    pub fn update_data(&mut self, hash_table: &HashTable) {
        // assumes that the last new element added to the
        // basis is the one of largest degree
        debug_assert!((self.previous_length..self.elements.len()).collect::<Vec<BasisLength>>()
                .windows(2)
                .all(|x| hash_table.degrees[self.elements[x[0]].monomials[0] as usize]
                    <= hash_table.degrees[self.elements[x[1]].monomials[0] as usize]));

        // check for constant
        if hash_table.degrees[
            self.elements[self.previous_length].monomials[0] as usize] == 0 {
            self.is_constant = true;
        }

        // update maximal total degree
        self.maximum_total_degree = max(
            self.maximum_total_degree,
            hash_table.degrees[self.elements[
                self.elements.len()-1].monomials[0] as usize]);

        // check redundancy due to resp. of new basis elements
        for i in self.previous_length..self.elements.len() {
            for j in 0..i {
                if !self.elements[j].is_redundant
                    && hash_table.divides(
                        self.elements[i].monomials[0],
                        self.elements[j].monomials[0]) {
                        self.elements[j].is_redundant = true;
                }
            }
            for j in self.previous_length..i {
                if !self.elements[j].is_redundant
                    && hash_table.divides(
                        self.elements[j].monomials[0],
                        self.elements[i].monomials[0]) {
                        self.elements[i].is_redundant = true;
                        break;
                }
            }
        }

        // update range for newly added elements
        self.previous_length = self.elements.len();
    }

    pub fn is_constant(&self) -> bool {
        return self.is_constant;
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
        let mut hash_table = HashTable::new(&exps);
        let basis = Basis::new::<i32>(&mut hash_table, fc, cfs, exps);
        assert_eq!(basis.elements[0].coefficients, [1,21840]);
        assert_eq!(basis.elements[0].monomials, [1,2]);
        assert_eq!(basis.elements[1].coefficients, [1,65520]);
        assert_eq!(basis.elements[1].monomials, [0,1]);
    }

    #[test]
    fn test_update_data() {
        let fc : Characteristic = 65521;
        let cfs : Vec<CoeffVec> = vec![vec![-2,65523], vec![1]];
        let exps : Vec<Vec<ExpVec>> = vec![vec![vec![0,3], vec![1,1]], vec![vec![0,0]]];
        let mut hash_table = HashTable::new(&exps);
        let mut basis = Basis::new::<i32>(&mut hash_table, fc, cfs, exps);
        assert_eq!(basis.previous_length, 0);
        basis.update_data(&hash_table);
        assert_eq!(basis.previous_length, basis.elements.len());
        assert_eq!(basis.maximum_total_degree, 3);
        assert_eq!(basis.is_constant, true);
    }

}
