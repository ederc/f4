use std::collections::HashSet;

use crate::primitives::*;

use crate::pairs::{
    PairSet,
};

use crate::hash_table::{
    HashTable,
};

use crate::basis::{
    Basis,
};

struct Row {
    basis_index: BasisLength,
    columns: MonomVec,
}

pub struct Matrix {
    rows: Vec<Row>,
}

impl Matrix {
    pub fn new() -> Matrix {
        let mat = Matrix {
            rows : Vec::new(),
        };

        return mat;
    }

    pub fn get_next_bunch_of_pairs(&mut self, basis: &Basis, pairs: &mut PairSet,
        hash_table: &mut HashTable) {

        let mut next_pairs = pairs.select_pairs_by_minimal_degree(hash_table);
        debug_assert!(next_pairs.len() > 0);

        next_pairs.sort_by(|a,b| hash_table.cmp_monomials_by_drl(a.lcm, b.lcm));
        let mut start = 0;
        let mut gens = HashSet::new();
        while start < next_pairs.len()  {
            let lcm = next_pairs[start].lcm;
            let stop = next_pairs[start..]
                .iter()
                .position(|p| p.lcm != lcm)
                .unwrap_or(next_pairs.len());

            for i in start..stop {
                gens.insert(next_pairs[i].generators.0);
                gens.insert(next_pairs[i].generators.1);
            }
            let vec_len = hash_table.nr_variables;
            for g in &gens {
                let mons = &basis.elements[*g].monomials;
                let mult_idx = hash_table.get_difference(lcm, mons[0]);
                let mut mult_mons: MonomVec = vec!(0; mons.len());
                for (idx, m) in mons.iter().enumerate() {
                    let mut exps: ExpVec = vec!(0; vec_len);
                    for i in 0..vec_len {
                        exps[i] = hash_table.monomials[mult_idx].exponents[i]
                            + hash_table.monomials[*m].exponents[i];
                    }
                    mult_mons[idx] = hash_table.insert(exps);
                }
                self.rows.push(Row { basis_index : *g, columns : mult_mons} );
            }
            start = stop;
            gens.clear();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_next_bunch_of_pairs() {
        let fc : Characteristic = 65521;
        let cfs : Vec<CoeffVec> = vec![vec![-2,65523], vec![1, -3], vec![1, -1], vec![1, 1]];
        let exps : Vec<Vec<ExpVec>> = vec![vec![vec![0,3,1], vec![1,1,0]],
            vec![vec![0,2,0], vec![1,1,0]], vec![vec![0,0,2], vec![1,0,0]],
            vec![vec![0,0,1], vec![0,0,0]]];
        let mut hash_table = HashTable::new(&exps);
        let basis = Basis::new::<i32>(&mut hash_table, fc, cfs, exps);

        let mut pairs = PairSet::new();
        pairs.update(&basis, &mut hash_table);

        let mut matrix = Matrix::new();

        matrix.get_next_bunch_of_pairs(&basis, &mut pairs, &mut hash_table);

        assert_eq!(matrix.rows.len(), 2);
        assert_eq!(matrix.rows[0].basis_index, 3);
        assert_eq!(matrix.rows[0].columns, [0,1]);
        assert_eq!(matrix.rows[1].basis_index, 0);
        assert_eq!(matrix.rows[1].columns, [0,11]);
    }
}
