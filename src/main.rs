mod primitives;
mod matrix;
mod io;
mod hash_table;
mod pairs;
mod arithmetic;
mod basis;

use crate::io::{
    read_file,
    Config,
};

use crate::basis::{
    Basis,
};

use crate::hash_table::{
    HashTable,
};

use crate::pairs::{
    PairSet,
};

use crate::matrix::{
    Matrix,
};

fn main() {
    let config = Config::new();
    let (variables, characteristic, coefficients, exponents) = read_file(config);
    let mut hash_table = HashTable::new(&exponents);
    let mut basis = Basis::new::<i32>(&mut hash_table, characteristic, coefficients, exponents);
    let mut pairs = PairSet::new();
    loop {
        if basis.previous_length < basis.elements.len() {
            pairs.update(&mut basis, &mut hash_table);
            basis.update_data(&hash_table);
        }
        if basis.is_constant() || pairs.is_empty() {
            break;
        }
        let mut matrix = Matrix::new();
        matrix.preprocessing(&basis, &mut pairs, &mut hash_table);
        matrix.reduce(&mut basis);
        matrix.postprocessing(&mut basis, &hash_table);
        for i in basis.previous_length..basis.elements.len() {
            println!("lm[{}] = {:?} (#mons {}) -- lc {}",
            i, hash_table.monomials[basis.elements[i].monomials[0]].exponents,
            basis.elements[i].monomials.len(),
            basis.elements[i].coefficients[0]);
        }
    }
    println!("((( final basis )))");

        for (i,e) in basis.elements.into_iter().filter(|x| x.is_redundant == false).enumerate() {
            println!("lm[{}] = {:?} (#mons {} -- redundant? {})",
            i, hash_table.monomials[e.monomials[0]].exponents,
            e.monomials.len(), e.is_redundant);
        }

    println!("done");


    // for c in coefficients {
    // println!("cfs {:?}", c);
    // }
    // let exp: ExpVec = vec![1,1,1];
    // println!("and char is {}", characteristic);
    // println!("size {}", std::mem::size_of::<usize>());
    // let mut hash_table: HashTable = HashTable::new(&exponents);
    // let mut basis: Basis = Basis:new(
    // for i in (0..exponents.len()).step_by(variables.len()) {
    //     println!("{:?}", &exponents[i..i+variables.len()]);
    //     map.insert(&exponents[i..i+variables.len()]);
    // }
    // for val in map.values() {
    //     println!("{val}");
    // let basis = basis::generate_initial_basis(
    //     characteristic,
    //     lengths,
    //     coefficients,
    //     exponents);
    // let basis = basis::Basis::<u32> {
    //         previous_length: 0,
    //         is_constant: false,
    //         maximum_total_degree: 0,
    //         is_redundant: vec![false; lengths.len()],
    //         monomials: monomials,
    //         coefficients: coefficients,
    //     }
}
