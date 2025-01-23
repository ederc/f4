use crate::primitives::*;

use crate::hash_table::{
    HashTable,
};

use crate::meta_data::{
    MetaData,
};

pub struct Element {
    pub coefficients: CoeffVec,
    pub monomials: MonomVec,
    pub is_redundant: bool,
}

pub struct Basis {
    previous_length: BasisLength,
    is_constant: bool,
    maximum_total_degree: Degree,
    elements: Vec<Element>,
}

impl Basis {
    pub fn new<T>(
        hash_table: &mut HashTable,
        meta_data: &MetaData,
        coefficients: Vec<CoeffVec>,
        exponents: Vec<Vec<ExpVec>>) -> Basis {
        debug_assert!(coefficients.len() > 0);
        debug_assert!(coefficients.len() == exponents.len());
        let mut basis = Basis {
            previous_length      : 0,
            is_constant          : false,
            maximum_total_degree : 0,
            elements             : Vec::new(),
        };
        basis.add_initial_elements(hash_table, meta_data, coefficients, exponents);
        basis.sort_terms_by_drl(hash_table);
        basis.sort_elements_by_drl(hash_table);
        return basis;
    }

    fn add_initial_elements(
        &mut self,
        hash_table: &mut HashTable,
        meta_data: &MetaData,
        cfs: Vec<CoeffVec>, exps: Vec<Vec<ExpVec>>) {
        let characteristic = meta_data.characteristic;
        for (c,e) in cfs.iter().zip(exps) {
            self.elements.push(Element {
                coefficients: c.iter().map(|a| a % characteristic).collect(),
                monomials: e.iter().map(|a| hash_table.insert(a.to_vec())).collect(),
                is_redundant: false});
        }

    }
    fn sort_terms_by_drl(&mut self, hash_table: &HashTable) {
        for e in &mut self.elements {
            let mut zipped: Vec<_> = e.monomials.clone().into_iter()
                .zip(e.coefficients.clone()).collect();
            zipped.sort_unstable_by(|(a, _), (b, _)| hash_table
                .cmp_monomials_by_drl(*a, *b));
                (e.monomials, e.coefficients) = zipped.into_iter().unzip();
        }
    }

    fn sort_elements_by_drl(&mut self, hash_table: &HashTable) {
        return self.elements.sort_by(
            |a,b| hash_table.cmp_monomials_by_drl(
                a.monomials[0], b.monomials[0]));
    }

    // pub fn update_pairs(&mut self, hash_table: &mut HashTable) {
    //     for (i,e) in self.elements[self.previous_length..].iter().enumerate() {
    //         let new_pairs = self.elements[..self.previous_length].iter().enumerate().map(|(j,f)|
    //             Pair {
    //                 lcm: hash_table.get_lcm(e.monomials[0], f.monomials[0]),
    //                 generators: (i, j),
    //                 // criterion: are_monomials_prime(e, f),
    //                 criterion: Criterion::Product,
    //             });
    //     }
    //     self.previous_length = self.elements.len();
    // }
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
        let meta_data = MetaData::new(fc, cfs.len());
        let basis = Basis::new::<i32>(&mut hash_table, &meta_data, cfs, exps);
        assert_eq!(basis.elements[0].coefficients, [1,-3]);
        assert_eq!(basis.elements[0].monomials, [2,1]);
        assert_eq!(basis.elements[1].coefficients, [2,-2]);
        assert_eq!(basis.elements[1].monomials, [1,0]);
    }
}
