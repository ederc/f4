use std::collections::HashSet;
use std::io::Write;
use std::io::stdout;
use rayon::prelude::*;

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
    pivot_lookup: Vec<usize>,
    nr_known_pivots: usize,
    density: f64,
}

impl Matrix {
    pub fn new() -> Matrix {
        let mat = Matrix {
            pivots  : Vec::new(),
            todo    : Vec::new(),
            columns : Vec::new(),
            pivot_lookup : Vec::new(),
            nr_known_pivots: 0,
            density: 0.0,
        };

        return mat;
    }

    fn get_next_bunch_of_pairs(&mut self, basis: &Basis, pairs: &mut PairSet,
        hash_table: &mut HashTable) {

        let mut next_pairs = pairs.select_pairs_by_minimal_degree(hash_table);
        debug_assert!(next_pairs.len() > 0);

        print!("{:3} {:7} {:7} ", hash_table.degrees[next_pairs[0].lcm as usize],
            next_pairs.len(), next_pairs.len() + pairs.list.len());
            next_pairs.sort_by(|a,b| hash_table.cmp_monomials_by_drl(a.lcm, b.lcm));
        stdout().flush().unwrap();
        let mut start = 0;
        let mut gens = HashSet::new();
        while start < next_pairs.len()  {
            let first_generator = next_pairs[start].generators.1;
            let lcm = next_pairs[start].lcm;
            // set index of lcm as done since we have at least a second generator
            // which plays the role as reducer of this monomial
            hash_table.indices[lcm as usize] = 2;
            self.columns.push(lcm);
            let stop = next_pairs[start..]
                .iter()
                .position(|p| p.lcm != lcm)
                .unwrap_or(next_pairs.len()-start) + start;

            gens.insert(next_pairs[start].generators.0);
            for i in start+1..stop {
                gens.insert(next_pairs[i].generators.1);
                gens.insert(next_pairs[i].generators.0);
            }
            debug_assert!(gens.len() > 0);
            let multiplier = hash_table.get_difference(
                lcm, basis.elements[first_generator as usize].monomials[0]);
            self.add_pivot(first_generator, &multiplier, basis, hash_table);

            gens.remove(&first_generator);
            for g in &gens {
                let multiplier = hash_table.get_difference(
                    lcm, basis.elements[*g as usize].monomials[0]);
                self.add_todo(*g, &multiplier, basis, hash_table);
            }
            start = stop;
            gens.clear();
        }
    }

    fn add_pivot(&mut self,
        divisor_idx: BasisLength, multiplier: &[Exponent],
        basis: &Basis, hash_table: &mut HashTable) {

        let mult_mons = hash_table.generate_multiplied_monomials(
            divisor_idx, &multiplier, basis);
        self.pivots.push(
            Row { basis_index : divisor_idx, columns : mult_mons} );
    }

    fn add_todo(&mut self,
        divisor_idx: BasisLength, multiplier: &[Exponent],
        basis: &Basis, hash_table: &mut HashTable) {

        let mult_mons = hash_table.generate_multiplied_monomials(
            divisor_idx, multiplier, basis);
        self.todo.push(
            Row { basis_index : divisor_idx, columns : mult_mons} );
    }

    fn get_reducers(&mut self, basis: &Basis, hash_table: &mut HashTable) {

        // get list of all lms and divmask
        let mut divisor_data_vec: 
            Vec<(DivisorMask,HashTableLength,BasisLength)> = Vec::new();
        for i in 0..basis.elements.len() {
            if basis.elements[i].is_redundant == false {
                divisor_data_vec.push((
                    hash_table.divisor_masks[basis.elements[i].monomials[0] as usize],
                    basis.elements[i].monomials[0],
                    i as BasisLength));
            }
        }
        let mut new_pivot_data: Vec<(BasisLength, ExpVec)> = Vec::new();

        for i in 0..self.todo.len() {
            for j in 0..self.todo[i].columns.len() {
                let c = self.todo[i].columns[j] as usize;
                if hash_table.indices[c] == 0 {
                    hash_table.indices[c] = 1;
                    self.columns.push(c as HashTableLength);
                    match hash_table.find_divisor(c as HashTableLength, &divisor_data_vec, &basis) {
                        Some((divisor_idx, multiplier)) =>
                         { self.add_pivot(divisor_idx, &multiplier, basis, hash_table); //new_pivot_data.push((divisor_idx, multiplier));
                           hash_table.indices[c] = 2 },
                        None => continue,
                    }
                }
            }
        }
        // number of pivots may change, thus use while loop instead of for loop
        // for updated checks on self.pivots.len()
        let mut i = 0;
        while i < self.pivots.len() {
            for j in 0..self.pivots[i].columns.len() {
                let c = self.pivots[i].columns[j] as usize;
                if hash_table.indices[c] == 0 {
                    hash_table.indices[c] = 1;
                    self.columns.push(c as HashTableLength);
                    match hash_table.find_divisor(c as HashTableLength, &divisor_data_vec, &basis) {
                        Some((divisor_idx, multiplier)) =>
                         { self.add_pivot(divisor_idx, &multiplier, basis, hash_table); //new_pivot_data.push((divisor_idx, multiplier));
                           hash_table.indices[c] = 2 },
                        None => continue,
                    }
                }
            }
            i += 1;
        }
        // for np in new_pivot_data {
        //     self.add_pivot(np.0, &np.1, basis, hash_table);
        // }
        // let mut curr_nr_pivs = 0;
        // while curr_nr_pivs < self.pivots.len() {
        //     let mut new_pivot_data: Vec<(BasisLength, ExpVec)> = Vec::new();
        //     for pivots in &self.pivots[curr_nr_pivs..self.pivots.len()] {
        //         for c in &pivots.columns {
        //             if hash_table.indices[*c as usize] == 0 {
        //                 hash_table.indices[*c as usize] = 1;
        //                 self.columns.push(*c);
        //                 match hash_table.find_divisor(*c, &divisor_data_vec, &basis) {
        //                     Some((divisor_idx, multiplier)) =>
        //                      { new_pivot_data.push((divisor_idx, multiplier));
        //                        hash_table.indices[*c as usize] = 2 },
        //                     None => continue,
        //                 }
        //             }
        //         }
        //     }
        //     curr_nr_pivs = self.pivots.len();
        //     for np in new_pivot_data {
        //         self.add_pivot(np.0, &np.1, basis, hash_table);
        //     }
        // }
    }

    fn convert_hashes_to_columns(&mut self, hash_table: &mut HashTable) {

        self.nr_known_pivots = self.pivots.len();
        // set colum index for corresponding monomial hash in hash table
        self.columns.sort_by(|a,b| hash_table.cmp_monomials_by_index_then_drl(*b, *a));
        for i in 0..self.columns.len() {
            hash_table.indices[self.columns[i] as usize] = i as HashTableLength;
        }
        // map hashes to columns in matrix
        for i in 0..self.todo.len() {
            for j in 0..self.todo[i].columns.len() {
                self.todo[i].columns[j] =
                    hash_table.indices[self.todo[i].columns[j] as usize];
            }
        }
        for i in 0..self.pivots.len() {
            for j in 0..self.pivots[i].columns.len() {
                self.pivots[i].columns[j] =
                    hash_table.indices[self.pivots[i].columns[j] as usize];
            }
        }
        // reset indices
        for i in 0..self.columns.len() {
            hash_table.indices[self.columns[i] as usize] = 0;
        }
    }

    fn get_density(&mut self) {
        let mat_size = (self.todo.len()+self.pivots.len()) * self.columns.len();
        let mut nr_nonzero_elements: usize = self.pivots.iter().map(|x| x.columns.len()).sum::<usize>();
        nr_nonzero_elements += self.todo.iter().map(|x| x.columns.len()).sum::<usize>();
        self.density = nr_nonzero_elements as f64 / mat_size as f64 * 100.0;
    }

    fn link_pivots_to_columns(&mut self) {
        self.pivot_lookup = vec!(usize::MAX; self.columns.len());

        for i in 0..self.pivots.len() {
            self.pivot_lookup[self.pivots[i].columns[0] as usize] = i;
        }
        self.get_density();
        print!(" {:7} x {:<7} {:8.2}%", self.todo.len()+self.pivots.len(),
            self.columns.len(), self.density);
        stdout().flush().unwrap();
    }

    // pub fn final_basis_reduction(&mut self, basis: &mut Basis,
    //     hash_table: &mut HashTable) {
    //
    //     // preprocessing
    //     self.get_non_redundant_basis_elements(basis, hash_table);
    //     self.get_reducers(basis, hash_table);
    //     self.convert_hashes_to_columns(hash_table);
    //     self.pivots.sort_by(|a,b| b.columns[0].cmp(&a.columns[0]));
    //     self.link_pivots_to_columns();
    //
    //     // reduction
    //     self.reduce
    // }

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
        let reducer = &self.pivots[self.pivot_lookup[col_idx]];
        let reducer_coefficients =
            &basis.elements[reducer.basis_index as usize].coefficients;
        let reducer_columns = &reducer.columns;
        debug_assert!(
            reducer_columns.len() == reducer_coefficients.len());

        let multiplier = dense_row[col_idx];

        // update dense row applying multiplied reducer
        reducer_columns
            .iter().zip(reducer_coefficients).for_each(|(a,b)|
                multiply_add_with_check(
                    &mut dense_row[*a as usize], multiplier, *b as DenseRowCoefficient, characteristic_2));
    }

    fn update_interreduced_pivot(&mut self, dense_row: DenseRow, col_idx: usize, basis: &mut Basis) {

        let (cols, cfs) = generate_sparse_row_from_dense_row(
            dense_row, col_idx, basis.characteristic as DenseRowCoefficient);

        let pivot_idx = self.pivot_lookup[col_idx];
        self.pivots[pivot_idx].columns = cols;
        basis.elements[self.pivots[pivot_idx].basis_index as usize].coefficients = cfs;
    }

    fn add_new_pivot(&mut self, dense_row: DenseRow, col_idx: usize, basis: &mut Basis) {

        let (cols, cfs) = generate_sparse_row_from_dense_row(
            dense_row, col_idx, basis.characteristic as DenseRowCoefficient);

        self.pivots.push(
            Row {
                basis_index: basis.elements.len() as BasisLength,
                columns: cols,});
        self.pivot_lookup[col_idx] = self.pivots.len()-1;
        basis.elements.push(
            Element {
                coefficients: cfs,
                monomials: Vec::new(),
                is_redundant: false,});
    }

    fn reduce_row(&mut self, idx: usize, basis: &mut Basis) {

        let row = &self.todo[idx];
        let mut dense_row: DenseRow = vec!(0; self.columns.len());

        let cfs = &basis.elements[row.basis_index as usize].coefficients;
        debug_assert!(cfs.len() == row.columns.len());

        let characteristic = basis.characteristic as DenseRowCoefficient;

        let start_column = row.columns[0] as usize;
        let last_column  = self.columns.len();

        for (i,c) in row.columns.iter().enumerate() {
            dense_row[*c as usize] = cfs[i] as DenseRowCoefficient;
        }

        let mut new_pivot_index = 0;
        for i in start_column..last_column {
            if dense_row[i] != 0 {
                dense_row[i] %= characteristic;
                if dense_row[i] != 0 {
                    if self.pivot_lookup[i] != usize::MAX {
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
        if row.columns.len() > 1 {
            let mut dense_row: DenseRow = vec!(0; self.columns.len());

            let cfs = &basis.elements[row.basis_index as usize].coefficients;
            debug_assert!(cfs.len() == row.columns.len());

            let characteristic = basis.characteristic as DenseRowCoefficient;

            let pivot_index  = row.columns[0] as usize;
            let start_column = row.columns[0+1] as usize;
            let last_column  = self.columns.len();

            for (i,c) in row.columns.iter().enumerate() {
                dense_row[*c as usize] = cfs[i] as DenseRowCoefficient;
            }

            for i in start_column..last_column {
                if dense_row[i] != 0 {
                    dense_row[i] %= characteristic;
                    if dense_row[i] != 0 {
                        if self.pivot_lookup[i] != usize::MAX {
                            self.apply_reducer(&mut dense_row, i, basis);
                        }
                    }
                }
            }
            self.update_interreduced_pivot(dense_row, pivot_index, basis);
        }
    }

    pub fn reduce(&mut self, basis: &mut Basis) {

        // set previous basis length before adding new elements / pivots
        basis.previous_length = basis.elements.len() as BasisLength;

        // println!("reducers for this matrix");
        // for p in &self.pivots {
        //     println!("red {}", p.basis_index);
        // }

        // find new pivots, reduce todo rows correspondingly
        (0..self.todo.len()).for_each(|i| {
            self.reduce_row(i, basis);
        });

        let nr_known_pivots = self.nr_known_pivots;
        // sort newly found pivots by decreasing column index
        // to prepare interreduction process
        self.pivots[nr_known_pivots..].sort_by(|a,b| a.columns[0].cmp(&b.columns[0]));
        for (i,r) in self.pivots[nr_known_pivots..].iter().enumerate() {
            self.pivot_lookup[r.columns[0] as usize] = i + nr_known_pivots;
        }

        // interreduce newly found pivots
        for i in nr_known_pivots..self.pivots.len() {
            self.interreduce_row(i, basis);
        }
    }

    pub fn postprocessing(&mut self, basis: &mut Basis, hash_table: &HashTable) {

        print!(" {:9} new {:9} zero",
            basis.elements.len() as i64 -basis.previous_length as i64,
            self.todo.len() as i64 - basis.elements.len() as i64
                + basis.previous_length as i64);
        stdout().flush().unwrap();
        // change column indices to monomial hash table positions
        self.pivots[self.nr_known_pivots..].iter_mut().for_each(|a|
            a.columns.iter_mut().for_each(|b| *b = self.columns[*b as usize]));

        // copy monomial data for new elements to basis
        self.pivots[self.nr_known_pivots..].iter().for_each(|a|
            basis.elements[a.basis_index as usize].monomials = a.columns.clone());

        basis.elements[(basis.previous_length as usize)..]
            .sort_by(|a,b| hash_table.cmp_monomials_by_drl(b.monomials[0], a.monomials[0]));
    }
}

fn multiply_add_with_check(
    coeff: &mut DenseRowCoefficient,
    multiplier: DenseRowCoefficient,
    reducer: DenseRowCoefficient,
    characteristic_squared: DenseRowCoefficient) {

    *coeff -= multiplier * reducer;
    *coeff += (*coeff >> 63) & characteristic_squared;
}

fn generate_sparse_row_from_dense_row(
    dense_row: DenseRow, col_idx: usize, characteristic: DenseRowCoefficient)
    -> (MonomVec, CoeffVec) {

    let mut cols: MonomVec = Vec::new();
    let mut cfs: CoeffVec = Vec::new();

    let lc = dense_row[col_idx] as Coefficient;

    if lc != 1 {
        let inv = modular_inverse(lc, characteristic as Characteristic) as DenseRowCoefficient;

        for (i, c) in dense_row[col_idx..].iter().enumerate() {
            if *c != 0 {
                cfs.push(((inv * *c) % characteristic) as Coefficient);
                cols.push((i+col_idx) as HashTableLength);
            }
        }
    } else {
        for (i, c) in dense_row[col_idx..].iter().enumerate() {
            if *c != 0 {
                cfs.push(*c as Coefficient);
                cols.push((i+col_idx) as HashTableLength);
            }
        }
    }
    return (cols, cfs);
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_density() {
        let fc : Characteristic = 65521;
        let cfs : Vec<CoeffVec> = vec![vec![-2,65523], vec![1, -3],
        vec![1, -1], vec![1, 1]];
        let exps : Vec<Vec<ExpVec>> = vec![vec![vec![0,3,1], vec![1,1,0]],
        vec![vec![0,2,0], vec![1,1,0]], vec![vec![0,0,2], vec![1,0,0]],
        vec![vec![0,0,1], vec![0,0,0]]];
        let mut hash_table = HashTable::new(&exps);
        let basis = Basis::new::<i32>(&mut hash_table, fc, cfs, exps);

        let mult: ExpVec = vec![0,0,2];
        let mut matrix = Matrix::new();

        matrix.add_todo(3, &mult, &basis, &mut hash_table);
        matrix.get_reducers(&basis, &mut hash_table);
        matrix.convert_hashes_to_columns(&mut hash_table);
        assert_eq!(matrix.columns.len(), 8);
        matrix.pivots.sort_by(|a,b| a.columns[0].cmp(&b.columns[0]));
        matrix.link_pivots_to_columns();
        matrix.get_density();
        assert_eq!(matrix.density, 25.0);
    }
    #[test]
    fn test_multiply_add_with_check() {
        let mut coeff: DenseRowCoefficient = 100;
        let multiplier: DenseRowCoefficient = 2;
        let reducer: DenseRowCoefficient = 31;
        let characteristic_squared: DenseRowCoefficient = 101*101;

        // without negative overflow
        multiply_add_with_check(&mut coeff, multiplier, reducer, characteristic_squared);
        assert_eq!(coeff, 38);
        // with negative overflow
        multiply_add_with_check(&mut coeff, multiplier, reducer, characteristic_squared);
        assert_eq!(coeff, 10177);
    }

    #[test]
    fn test_generate_sparse_row_from_dense_row() {
        let dense_row: DenseRow = [0,2,3,0,0,23].to_vec();
        let characteristic: DenseRowCoefficient = 101;
        let (cols, cfs) = generate_sparse_row_from_dense_row(
            dense_row, 1, characteristic);
        assert_eq!(cols, [1,2,5]);
        assert_eq!(cfs, [1,52,62]);
    }
    #[test]
    fn test_interreduce_row() {
        let fc : Characteristic = 65521;
        let cfs : Vec<CoeffVec> = vec![vec![-2,65523], vec![1, -3],
        vec![1, -1], vec![1, 1]];
        let exps : Vec<Vec<ExpVec>> = vec![vec![vec![0,3,1], vec![1,1,0]],
        vec![vec![0,2,0], vec![1,1,0]], vec![vec![0,0,2], vec![1,0,0]],
        vec![vec![0,0,1], vec![0,0,0]]];
        let mut hash_table = HashTable::new(&exps);
        let mut basis = Basis::new::<i32>(&mut hash_table, fc, cfs, exps);

        let mult: ExpVec = vec![0,0,2];
        let mut matrix = Matrix::new();

        matrix.add_todo(3, &mult, &basis, &mut hash_table);
        matrix.get_reducers(&basis, &mut hash_table);
        matrix.convert_hashes_to_columns(&mut hash_table);
        matrix.pivots.sort_by(|a,b| a.columns[0].cmp(&b.columns[0]));
        matrix.link_pivots_to_columns();
        matrix.reduce_row(0, &mut basis);
        matrix.interreduce_row(0, &mut basis);
        assert_eq!(matrix.pivots.len(), 7);
        assert_eq!(matrix.pivots[0].columns, [0,7]);
        assert_eq!(matrix.pivots[0].basis_index, 0);
        assert_eq!(
            basis.elements[matrix.pivots[0].basis_index as usize].coefficients,
            [1,21840]);
    }
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

        let mult: ExpVec = vec![0,0,2];
        let mut matrix = Matrix::new();

        matrix.add_todo(3, &mult, &basis, &mut hash_table);
        matrix.get_reducers(&basis, &mut hash_table);
        matrix.convert_hashes_to_columns(&mut hash_table);
        matrix.pivots.sort_by(|a,b| a.columns[0].cmp(&b.columns[0]));
        matrix.link_pivots_to_columns();
        matrix.reduce_row(0, &mut basis);
        assert_eq!(matrix.pivots.len(), 7);
        assert_eq!(matrix.pivots[6].columns, [6,7]);
        assert_eq!(matrix.pivots[6].basis_index, 4);
        assert_eq!(basis.elements[matrix.pivots[6].basis_index as usize].coefficients, [1,43681]);
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
        let mut matrix = Matrix::new();

        matrix.add_todo(0, &mult, &basis, &mut hash_table);
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
        let mut matrix = Matrix::new();

        matrix.add_todo(0, &mult, &basis, &mut hash_table);
        matrix.get_reducers(&basis, &mut hash_table);

        assert_eq!(matrix.todo.len(), 1);
        assert_eq!(matrix.todo[0].basis_index, 0);
        assert_eq!(hash_table.exponents[matrix.todo[0].columns[0] as usize], [1,1,1]);
        assert_eq!(hash_table.exponents[matrix.todo[0].columns[1] as usize], [1,1,0]);
        assert_eq!(matrix.pivots.len(), 2);
        assert_eq!(matrix.pivots[0].basis_index, 0);
        assert_eq!(hash_table.exponents[matrix.pivots[0].columns[0] as usize], [1,1,1]);
        assert_eq!(hash_table.exponents[matrix.pivots[0].columns[1] as usize], [1,1,0]);
        assert_eq!(matrix.pivots[1].basis_index, 2);
        assert_eq!(hash_table.exponents[matrix.pivots[1].columns[0] as usize], [1,1,0]);
        assert_eq!(hash_table.exponents[matrix.pivots[1].columns[1] as usize], [0,2,0]);
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
        let mut matrix = Matrix::new();

        matrix.add_todo(0, &mult, &basis, &mut hash_table);
        assert_eq!(matrix.todo[0].basis_index, 0);
        assert_eq!(matrix.todo[0].columns, [1,8]);
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
        let mut matrix = Matrix::new();

        matrix.add_pivot(0, &mult, &basis, &mut hash_table);
        assert_eq!(matrix.pivots[0].basis_index, 0);
        assert_eq!(matrix.pivots[0].columns, [1,8]);
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

        assert_eq!(matrix.columns.len(), 1);
        assert_eq!(matrix.todo.len(), 1);
        assert_eq!(matrix.todo[0].basis_index, 1);
        assert_eq!(matrix.todo[0].columns, [4,5]);
        assert_eq!(matrix.pivots.len(), 1);
        assert_eq!(matrix.pivots[0].basis_index, 0);
        assert_eq!(matrix.pivots[0].columns, [4,6]);
    }
}
