use std::collections::HashSet;

use crate::arithmetic::i32::{
    modular_inverse,
};

use crate::primitives::*;

use crate::pairs::{
    PairSet,
};

use crate::hash_table::{
    HashTable,
};

use crate::basis::{
    Basis,
    Element,
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

    fn apply_reducer(&self, dense_row: &mut DenseRow, col_idx: usize, basis: &Basis) {
        let characteristic_2 = (basis.characteristic as DenseRowCoefficient).pow(2);
        let reducer = &self.pivots[self.columns[col_idx]];
        let reducer_coefficients =
            &basis.elements[reducer.basis_index].coefficients;
        let reducer_columns = &reducer.columns;
        debug_assert!(
            reducer_columns.len() == reducer_coefficients.len());
        let multiplier = dense_row[col_idx];

        // update dense row applying multiplied reducer
        for (col, cf) in reducer_columns
            .iter().zip(reducer_coefficients) {

            dense_row[*col] -= multiplier * *cf as DenseRowCoefficient;
            if dense_row[*col].is_negative() {
                dense_row[*col] += characteristic_2;
            }
        }
    }

    fn update_interreduced_pivot(&mut self, dense_row: DenseRow, col_idx: usize, basis: &mut Basis) {

        let (cols, cfs) = generate_sparse_row_from_dense_row(
            dense_row, col_idx, basis.characteristic as DenseRowCoefficient);

        let pivot_idx = self.columns[col_idx];
        self.pivots[pivot_idx].columns = cols;
        basis.elements[self.pivots[pivot_idx].basis_index].coefficients = cfs;
    }

    fn add_new_pivot(&mut self, dense_row: DenseRow, col_idx: usize, basis: &mut Basis) {

        let (cols, cfs) = generate_sparse_row_from_dense_row(
            dense_row, col_idx, basis.characteristic as DenseRowCoefficient);

        self.pivots.push(
            Row {
                basis_index: basis.elements.len(),
                columns: cols,});
        self.columns[col_idx] = self.pivots.len();
        basis.elements.push(
            Element {
                coefficients: cfs,
                monomials: Vec::new(),
                is_redundant: false,});
    }

    fn reduce_row(&mut self, idx: usize, basis: &mut Basis) {

        let row = &self.todo[idx];
        let mut dense_row: DenseRow = vec!(0; self.columns.len());

        let cfs = &basis.elements[row.basis_index].coefficients;
        debug_assert!(cfs.len() == row.columns.len());

        let characteristic = basis.characteristic as DenseRowCoefficient;

        let start_column = row.columns[0];
        let last_column  = self.columns.len();

        for (i,c) in row.columns.iter().enumerate() {
            dense_row[*c] = cfs[i] as DenseRowCoefficient;
        }

        let mut new_pivot_index = 0;
        for i in start_column..last_column {
            if dense_row[i] != 0 {
                dense_row[i] %= characteristic;
                if dense_row[i] != 0 {
                    if self.columns[i] != HashTableLength::MAX {
                        self.apply_reducer(&mut dense_row, i, basis);
                    } else {
                        if new_pivot_index == 0 {
                            new_pivot_index = i;
                        }
                        continue;
                    }
                }
            }
        }
        if new_pivot_index != 0 {
            self.add_new_pivot(dense_row, new_pivot_index, basis);
        }
    }

    fn interreduce_row(&mut self, idx: usize, basis: &mut Basis) {

        let row = &self.pivots[idx];
        let mut dense_row: DenseRow = vec!(0; self.columns.len());

        let cfs = &basis.elements[row.basis_index].coefficients;
        debug_assert!(cfs.len() == row.columns.len());

        let characteristic = basis.characteristic as DenseRowCoefficient;

        let start_column = row.columns[0+1];
        let last_column  = self.columns.len();

        for (i,c) in row.columns.iter().enumerate() {
            dense_row[*c] = cfs[i] as DenseRowCoefficient;
        }

        let mut new_pivot_index = 0;
        for i in start_column..last_column {
            if dense_row[i] != 0 {
                dense_row[i] %= characteristic;
                if dense_row[i] != 0 {
                    if self.columns[i] != HashTableLength::MAX {
                        self.apply_reducer(&mut dense_row, i, basis);
                    } else {
                        if new_pivot_index == 0 {
                            new_pivot_index = i;
                        }
                        continue;
                    }
                }
            }
        }
        self.update_interreduced_pivot(dense_row, new_pivot_index, basis);
    }

    pub fn reduce(&mut self, basis: &mut Basis) {

        // find new pivots, reduce todo rows correspondingly
        for i in 0..self.todo.len() {
            self.reduce_row(i, basis);
        }

        let nr_known_pivots = self.nr_known_pivots;
        // sort newly found pivots by decreasing column index
        // to prepare interreduction process
        self.pivots[nr_known_pivots..].sort_by(|a,b| a.columns[0].cmp(&b.columns[0]));
        for (i,r) in self.pivots[nr_known_pivots..].iter().enumerate() {
            self.columns[r.columns[0]] = i + nr_known_pivots;
        }

        // interreduce newly found pivots
        for i in nr_known_pivots..self.pivots.len() {
            self.reduce_row(i, basis);
        }
    }
}

fn generate_sparse_row_from_dense_row(
    dense_row: DenseRow, col_idx: usize, characteristic: DenseRowCoefficient)
    -> (MonomVec, CoeffVec) {

    let mut cols: MonomVec = Vec::new();
    let mut cfs: CoeffVec = Vec::new();

    let lc = dense_row[col_idx] as Coefficient;

    if lc != 1 {
        let inv = modular_inverse(lc, characteristic as Characteristic);

        for (i, c) in dense_row[col_idx..].iter().enumerate() {
            if *c != 0 {
                cfs.push(((inv as DenseRowCoefficient * *c) % characteristic) as Coefficient);
                cols.push(i+col_idx);
            }
        }
    } else {
        for (i, c) in dense_row[col_idx..].iter().enumerate() {
            if *c != 0 {
                cfs.push(*c as Coefficient);
                cols.push(i+col_idx);
            }
        }
    }
    return (cols, cfs);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reduce_row() {
        let fc : Characteristic = 65521;
        let cfs : Vec<CoeffVec> = vec![vec![-2,65523], vec![1, -3],
        vec![1, -1], vec![1, 1]];
        let exps : Vec<Vec<ExpVec>> = vec![vec![vec![0,3,1], vec![1,1,0]],
        vec![vec![0,2,0], vec![1,1,0]], vec![vec![0,0,2], vec![1,0,0]],
        vec![vec![0,0,1], vec![0,0,0]]];
        let mut hash_table = HashTable::new(&exps);
        let mut basis = Basis::new::<i32>(&mut hash_table, fc, cfs, exps);

        let mut matrix = Matrix::new();

        let mult: ExpVec = vec![0,0,2];
        let mult_idx = hash_table.insert(mult);

        let mut matrix = Matrix::new();

        matrix.add_todo(3, mult_idx, &basis, &mut hash_table);
        matrix.get_reducers(&basis, &mut hash_table);
        matrix.convert_hashes_to_columns(&mut hash_table);
        matrix.pivots.sort_by(|a,b| a.columns[0].cmp(&b.columns[0]));
        matrix.link_pivots_to_columns();
        matrix.reduce_row(0, &mut basis);
        assert_eq!(matrix.pivots.len(), 7);
        assert_eq!(matrix.pivots[6].columns, [4,6,7]);
        assert_eq!(matrix.pivots[6].basis_index, 4);
        assert_eq!(basis.elements[matrix.pivots[6].basis_index].coefficients, [1,4,65520]);
    }
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
    #[test]
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
