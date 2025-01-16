mod types;
mod io;
mod hash_table;
mod meta_data;
mod arithmetic;
mod basis;

use crate::io::{
    read_file,
    Config,
};
use crate::hash_table::{
    HashTable,
};

use crate::types::*;

fn main() {
    let config = Config::new();
    let (variables, characteristic, coefficients, exponents) = read_file(config);
    for c in coefficients {
    println!("cfs {:?}", c);
    }
    let exp: ExpVec = vec![1,1,1];
    println!("and char is {}", characteristic);
    println!("size {}", std::mem::size_of::<usize>());
    let mut hash_table: HashTable = HashTable::new(&exponents);
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
