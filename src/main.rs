mod primitives;
mod matrix;
mod io;
mod hash_table;
mod pairs;
mod arithmetic;
mod basis;

use std::time::{
    Duration,
    Instant,
};

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
    let overall_time = Instant::now();
    let mut pre_overall_time = Duration::new(0,0);
    let mut up_overall_time  = Duration::new(0,0);
    let mut la_overall_time  = Duration::new(0,0);
    let config = Config::new();
    let (variables, characteristic, coefficients, exponents) = read_file(config);
    let mut hash_table = HashTable::new(&exponents);
    let mut basis = Basis::new::<i32>(&mut hash_table, characteristic, coefficients, exponents);
    let mut pairs = PairSet::new();
    println!("deg     sel   pairs        matrix        density                 data               round time");
    println!("----------------------------------------------------------------------------------------------");
    loop {
        let rd_time = Instant::now();
        if basis.previous_length < basis.elements.len() {
            let up_time = Instant::now();
            pairs.update(&mut basis, &mut hash_table);
            basis.update_data(&hash_table);
            up_overall_time += up_time.elapsed();
        }
        if basis.is_constant() || pairs.is_empty() {
            break;
        }
        let mut matrix = Matrix::new();
        let pre_time = Instant::now();
        matrix.preprocessing(&basis, &mut pairs, &mut hash_table);
        pre_overall_time += pre_time.elapsed();
        let la_time = Instant::now();
        matrix.reduce(&mut basis);
        la_overall_time += la_time.elapsed();
        matrix.postprocessing(&mut basis, &hash_table);
        println!("{:13.3} sec ", rd_time.elapsed().as_secs_f64());
    }

    println!("----------------------------------------------------------------------------------------------");
    println!("length of basis: {}", basis.elements
        .into_iter()
        .filter(|x| x.is_redundant == false)
        .collect::<Vec<_>>()
        .len());
    println!("----------------------------------------------------------------------------------------------");
    println!("overall time:        {:13.3} sec", overall_time.elapsed().as_secs_f64());
    println!("----------------------------------------------------------------------------------------------");
    println!("update time:         {:13.3} sec {:6.1}%", up_overall_time.as_secs_f64(),
        up_overall_time.as_secs_f64() / overall_time.elapsed().as_secs_f64() * 100.0);
    println!("preprocessing time:  {:13.3} sec {:6.1}%", pre_overall_time.as_secs_f64(),
        pre_overall_time.as_secs_f64() / overall_time.elapsed().as_secs_f64() * 100.0);
    println!("linear algebra time: {:13.3} sec {:6.1}%", la_overall_time.as_secs_f64(),
        la_overall_time.as_secs_f64() / overall_time.elapsed().as_secs_f64() * 100.0);
    println!("----------------------------------------------------------------------------------------------");

    // for (i,e) in basis.elements.into_iter().filter(|x| x.is_redundant == false).enumerate() {
    //     println!("lm[{}] = {:?} (#mons {} -- redundant? {})",
    //     i, hash_table.monomials[e.monomials[0]].exponents,
    //     e.monomials.len(), e.is_redundant);
    // }
}
