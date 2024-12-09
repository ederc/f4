use crate::io::file_handling::{
    read_file,
    Config
};

mod io;

fn main() {
    let config = Config::new();
    let (variables, characteristic, lengths, coefficients, exponents) = read_file(config);
    println!("cfs {:?}", coefficients);
    println!("len {:?}", lengths);
    println!("and char is {}", characteristic);
}
