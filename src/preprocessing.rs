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
    pivots: Vec<Row>,
    todo: Vec<Row>,
    columns: Vec<HashTableLength>,
    nr_known_pivots: usize,
}

impl Matrix {
    pub fn new() -> Matrix {
        let mat = Matrix {
            pivots  : Vec::new(),
            todo    : Vec::new(),
            columns : Vec::new(),
            nr_known_pivots: 0,
        };

        return mat;
    }

    fn get_next_bunch_of_pairs(&mut self, basis: &Basis, pairs: &mut PairSet,
        hash_table: &mut HashTable) {

        let mut next_pairs = pairs.select_pairs_by_minimal_degree(hash_table);
        debug_assert!(next_pairs.len() > 0);

        next_pairs.sort_by(|a,b| hash_table.cmp_monomials_by_drl(a.lcm, b.lcm));
        let mut start = 0;
        let mut gens = HashSet::new();
        while start < next_pairs.len()  {
            let first_generator = next_pairs[start].generators.0;
            let lcm = next_pairs[start].lcm;
            // set index of lcm as done since we have at least a second generator
            // which plays the role as reducer of this monomial
            hash_table.indices[lcm] = 1;
            self.columns.push(lcm);
            let stop = next_pairs[start..]
                .iter()
                .position(|p| p.lcm != lcm)
                .unwrap_or(next_pairs.len());

            gens.insert(next_pairs[start].generators.1);
            for i in start+1..stop {
                gens.insert(next_pairs[i].generators.0);
                gens.insert(next_pairs[i].generators.1);
            }
            debug_assert!(gens.len() > 0);
            let mult_idx = hash_table.get_difference(
                lcm, basis.elements[first_generator].monomials[0]);
            self.columns.push(lcm);
            self.add_pivot(first_generator, mult_idx, basis, hash_table);

            for g in &gens {
                let mult_idx = hash_table.get_difference(
                    lcm, basis.elements[*g].monomials[0]);
                self.add_todo(*g, mult_idx, basis, hash_table);
            }
            start = stop;
            gens.clear();
        }
    }

    fn add_pivot(&mut self,
        divisor_idx: BasisLength, mult_idx: HashTableLength,
        basis: &Basis, hash_table: &mut HashTable) {

        let mult_mons = hash_table.generate_multiplied_monomials(
            divisor_idx, mult_idx, basis);
        self.pivots.push(
            Row { basis_index : divisor_idx, columns : mult_mons} );
    }

    fn add_todo(&mut self,
        divisor_idx: BasisLength, mult_idx: HashTableLength,
        basis: &Basis, hash_table: &mut HashTable) {

        let mult_mons = hash_table.generate_multiplied_monomials(
            divisor_idx, mult_idx, basis);
        self.todo.push(
            Row { basis_index : divisor_idx, columns : mult_mons} );
    }

    fn get_reducers(&mut self, basis: &Basis, hash_table: &mut HashTable) {

        let mut i = 0;
        while i < self.todo.len() {
            for j in 0..self.todo[i].columns.len() {
                if hash_table.indices[self.todo[i].columns[j]] == 0 {
                    hash_table.indices[self.todo[i].columns[j]] = 1;
                    self.columns.push(self.todo[i].columns[j]);
                    match hash_table.find_divisor(self.todo[i].columns[j], basis) {
                        Some((divisor_idx, multiplier)) =>
                            self.add_pivot(divisor_idx, multiplier, basis, hash_table),
                        None => continue,
                    }
                }
            }
            i += 1;
        }
        i = 0;
        while i < self.pivots.len() {
            for j in 0..self.pivots[i].columns.len() {
                if hash_table.indices[self.pivots[i].columns[j]] == 0 {
                    hash_table.indices[self.pivots[i].columns[j]] = 1;
                    self.columns.push(self.pivots[i].columns[j]);
                    match hash_table.find_divisor(self.pivots[i].columns[j], basis) {
                        Some((divisor_idx, multiplier)) =>
                            self.add_pivot(divisor_idx, multiplier, basis, hash_table),
                        None => continue,
                    }
                }
            }
            i += 1;
        }
    }

    fn convert_hashes_to_columns(&mut self, hash_table: &mut HashTable) {

        self.nr_known_pivots = self.pivots.len();
        // set colum index for corresponding monomial hash in hash table
        self.columns.sort_by(|a,b| hash_table.cmp_monomials_by_drl(*b, *a));
        for i in 0..self.columns.len() {
            hash_table.indices[self.columns[i]] = i;
        }
        // map hashes to columns in matrix
        for i in 0..self.todo.len() {
            for j in 0..self.todo[i].columns.len() {
                self.todo[i].columns[j] =
                    hash_table.indices[self.todo[i].columns[j]];
            }
        }
        for i in 0..self.pivots.len() {
            for j in 0..self.pivots[i].columns.len() {
                self.pivots[i].columns[j] =
                    hash_table.indices[self.pivots[i].columns[j]];
            }
        }
        // reset indices
        for i in 0..self.columns.len() {
            hash_table.indices[self.columns[i]] = 0;
        }
    }

    fn link_pivots_to_columns(&mut self) {
        self.columns.iter_mut().for_each(|a| *a = HashTableLength::MAX);

        for i in 0..self.pivots.len() {
            self.columns[self.pivots[i].columns[0]] = i;
        }
    }

    pub fn preprocessing(&mut self, basis: &Basis,
        pairs: &mut PairSet, hash_table: &mut HashTable) {

        self.get_next_bunch_of_pairs(basis, pairs, hash_table);
        self.get_reducers(basis, hash_table);
        self.convert_hashes_to_columns(hash_table);
        self.pivots.sort_by(|a,b| b.columns[0].cmp(&a.columns[0]));
        self.link_pivots_to_columns();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_convert_hashes_to_columns() {
        let fc : Characteristic = 65521;
        let cfs : Vec<CoeffVec> = vec![vec![-2,65523], vec![1, -3],
                    vec![1, -1], vec![1, 1]];
        let exps : Vec<Vec<ExpVec>> = vec![vec![vec![0,3,1], vec![1,1,0]],
            vec![vec![0,2,0], vec![1,1,0]], vec![vec![0,0,2], vec![1,0,0]],
            vec![vec![0,0,1], vec![0,0,0]]];
        let mut hash_table = HashTable::new(&exps);
        let basis = Basis::new::<i32>(&mut hash_table, fc, cfs, exps);
        let mult: ExpVec = vec![1,1,0];
        let mult_idx = hash_table.insert(mult);

        let mut matrix = Matrix::new();

        matrix.add_todo(0, mult_idx, &basis, &mut hash_table);
        matrix.get_reducers(&basis, &mut hash_table);
        matrix.convert_hashes_to_columns(&mut hash_table);

        assert_eq!(matrix.columns.len(), 3);
        assert_eq!(matrix.todo.len(), 1);
        assert_eq!(matrix.todo[0].columns[0], 0);
        assert_eq!(matrix.todo[0].columns[1], 1);
        assert_eq!(matrix.pivots.len(), 2);
        assert_eq!(matrix.pivots[0].columns[0], 0);
        assert_eq!(matrix.pivots[0].columns[1], 1);
        assert_eq!(matrix.pivots[1].columns[0], 1);
        assert_eq!(matrix.pivots[1].columns[1], 2);
        assert_eq!(matrix.nr_known_pivots, 2);
    }

    fn test_get_reducers() {
        let fc : Characteristic = 65521;
        let cfs : Vec<CoeffVec> = vec![vec![-2,65523], vec![1, -3],
                    vec![1, -1], vec![1, 1]];
        let exps : Vec<Vec<ExpVec>> = vec![vec![vec![0,3,1], vec![1,1,0]],
            vec![vec![0,2,0], vec![1,1,0]], vec![vec![0,0,2], vec![1,0,0]],
            vec![vec![0,0,1], vec![0,0,0]]];
        let mut hash_table = HashTable::new(&exps);
        let basis = Basis::new::<i32>(&mut hash_table, fc, cfs, exps);
        let mult: ExpVec = vec![1,1,0];
        let mult_idx = hash_table.insert(mult);

        let mut matrix = Matrix::new();

        matrix.add_todo(0, mult_idx, &basis, &mut hash_table);
        matrix.get_reducers(&basis, &mut hash_table);

        assert_eq!(matrix.todo.len(), 1);
        assert_eq!(matrix.todo[0].basis_index, 0);
        assert_eq!(hash_table.monomials[matrix.todo[0].columns[0]].exponents, [1,1,1]);
        assert_eq!(hash_table.monomials[matrix.todo[0].columns[1]].exponents, [1,1,0]);
        assert_eq!(matrix.pivots.len(), 2);
        assert_eq!(matrix.pivots[0].basis_index, 0);
        assert_eq!(hash_table.monomials[matrix.pivots[0].columns[0]].exponents, [1,1,1]);
        assert_eq!(hash_table.monomials[matrix.pivots[0].columns[1]].exponents, [1,1,0]);
        assert_eq!(matrix.pivots[1].basis_index, 2);
        assert_eq!(hash_table.monomials[matrix.pivots[1].columns[0]].exponents, [1,1,0]);
        assert_eq!(hash_table.monomials[matrix.pivots[1].columns[1]].exponents, [0,2,0]);
    }
    #[test]
    fn test_add_todo() {
        let fc : Characteristic = 65521;
        let cfs : Vec<CoeffVec> = vec![vec![-2,65523], vec![1, -3],
                    vec![1, -1], vec![1, 1]];
        let exps : Vec<Vec<ExpVec>> = vec![vec![vec![0,3,1], vec![1,1,0]],
            vec![vec![0,2,0], vec![1,1,0]], vec![vec![0,0,2], vec![1,0,0]],
            vec![vec![0,0,1], vec![0,0,0]]];
        let mut hash_table = HashTable::new(&exps);
        let basis = Basis::new::<i32>(&mut hash_table, fc, cfs, exps);
        let mult: ExpVec = vec![0,3,0];
        let mult_idx = hash_table.insert(mult);

        let mut matrix = Matrix::new();

        matrix.add_todo(0, mult_idx, &basis, &mut hash_table);
        assert_eq!(matrix.todo[0].basis_index, 0);
        assert_eq!(matrix.todo[0].columns, [0,7]);
    }

    #[test]
    fn test_add_pivot() {
        let fc : Characteristic = 65521;
        let cfs : Vec<CoeffVec> = vec![vec![-2,65523], vec![1, -3],
                    vec![1, -1], vec![1, 1]];
        let exps : Vec<Vec<ExpVec>> = vec![vec![vec![0,3,1], vec![1,1,0]],
            vec![vec![0,2,0], vec![1,1,0]], vec![vec![0,0,2], vec![1,0,0]],
            vec![vec![0,0,1], vec![0,0,0]]];
        let mut hash_table = HashTable::new(&exps);
        let basis = Basis::new::<i32>(&mut hash_table, fc, cfs, exps);
        let mult: ExpVec = vec![0,3,0];
        let mult_idx = hash_table.insert(mult);

        let mut matrix = Matrix::new();

        matrix.add_pivot(0, mult_idx, &basis, &mut hash_table);
        assert_eq!(matrix.pivots[0].basis_index, 0);
        assert_eq!(matrix.pivots[0].columns, [0,7]);
    }

    #[test]
    fn test_get_next_bunch_of_pairs() {
        let fc : Characteristic = 65521;
        let cfs : Vec<CoeffVec> = vec![vec![-2,65523], vec![1, -3],
                    vec![1, -1], vec![1, 1]];
        let exps : Vec<Vec<ExpVec>> = vec![vec![vec![0,3,1], vec![1,1,0]],
            vec![vec![0,2,0], vec![1,1,0]], vec![vec![0,0,2], vec![1,0,0]],
            vec![vec![0,0,1], vec![0,0,0]]];
        let mut hash_table = HashTable::new(&exps);
        let basis = Basis::new::<i32>(&mut hash_table, fc, cfs, exps);

        let mut pairs = PairSet::new();
        pairs.update(&basis, &mut hash_table);

        let mut matrix = Matrix::new();

        matrix.get_next_bunch_of_pairs(&basis, &mut pairs, &mut hash_table);

        assert_eq!(matrix.todo.len(), 1);
        assert_eq!(matrix.todo[0].basis_index, 0);
        assert_eq!(matrix.todo[0].columns, [0,11]);
        assert_eq!(matrix.pivots.len(), 1);
        assert_eq!(matrix.pivots[0].basis_index, 3);
        assert_eq!(matrix.pivots[0].columns, [0,1]);
    }
}
