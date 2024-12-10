mod io;
mod arithmetic;
use crate::io::file_handling::{
    read_file,
    Config
};

use std::collections::HashMap;

fn main() {
    let config = Config::new();
    let (variables, characteristic, lengths, coefficients, exponents) = read_file(config);
    println!("cfs {:?}", coefficients);
    println!("len {:?}", lengths);
    println!("and char is {}", characteristic);
    let mut map: HashMap<Vec<i32>, i32> = HashMap::new();
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
