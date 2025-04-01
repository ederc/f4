mod primitives;
mod matrix;
mod io;
mod hash_table;
mod pairs;
mod update;
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

use crate::update::{
    LcmVec,
    update_lcms,
};

use crate::matrix::{
    Matrix,
};

fn main() {
    let overall_time = Instant::now();
    let mut pre_overall_time = Duration::new(0,0);
    let mut post_overall_time = Duration::new(0,0);
    let mut up_overall_time  = Duration::new(0,0);
    let mut la_overall_time  = Duration::new(0,0);
    let config = Config::new();
    let (variables, characteristic, coefficients, exponents) = read_file(config);
    let mut hash_table = HashTable::new(&exponents);
    let mut basis = Basis::new::<i32>(&mut hash_table, characteristic, coefficients, exponents);
    let mut pairs = PairSet::new();
    let mut lcm_vec = LcmVec::new();

    println!("deg     sel   pairs        matrix        density              new  data             round time");
    println!("----------------------------------------------------------------------------------------------");
    loop {
        let rd_time = Instant::now();
        if (basis.previous_length as usize) < basis.elements.len() {
            let up_time = Instant::now();
            update_lcms(&mut lcm_vec, &mut basis, &mut hash_table);
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
        let post_time = Instant::now();
        matrix.postprocessing(&mut basis, &hash_table);
        post_overall_time += post_time.elapsed();
        println!("{:13.3} sec ", rd_time.elapsed().as_secs_f64());
    }

    // // final reduction
    // let rd_time = Instant::now();
    // let final_reduction_time = Instant::now();
    // let mut matrix = Matrix::new();
    // matrix.final_basis_reduction(&mut basis, &mut hash_table);
    //
    // // printing info
    // println!("{:13.3} sec ", rd_time.elapsed().as_secs_f64());
    println!("----------------------------------------------------------------------------------------------");
    println!("length of basis: {}", basis.elements
        .into_iter()
        .filter(|x| x.is_redundant == false)
        .collect::<Vec<_>>()
        .len());
    println!("----------------------------------------------------------------------------------------------");
    println!("hash table capacity: {}", hash_table.length.ilog2());
    println!("hash table insert data: {} -> {} ({:3.3}%) -> {} ({:3.3}%)", hash_table.nr_in,
        hash_table.nr_in_ex,
        (hash_table.nr_in_ex as f64) / (hash_table.nr_in as f64) * 100.0,
        hash_table.nr_in_new,
        (hash_table.nr_in_new as f64) / (hash_table.nr_in as f64) * 100.0);
    println!("----------------------------------------------------------------------------------------------");
    println!("overall time:        {:13.3} sec", overall_time.elapsed().as_secs_f64());
    println!("----------------------------------------------------------------------------------------------");
    println!("update time:         {:13.3} sec {:6.1}%", up_overall_time.as_secs_f64(),
        up_overall_time.as_secs_f64() / overall_time.elapsed().as_secs_f64() * 100.0);
    println!("preprocessing time:  {:13.3} sec {:6.1}%", pre_overall_time.as_secs_f64(),
        pre_overall_time.as_secs_f64() / overall_time.elapsed().as_secs_f64() * 100.0);
    println!("linear algebra time: {:13.3} sec {:6.1}%", la_overall_time.as_secs_f64(),
        la_overall_time.as_secs_f64() / overall_time.elapsed().as_secs_f64() * 100.0);
    println!("postprocessingt time: {:12.3} sec {:6.1}%", post_overall_time.as_secs_f64(),
        post_overall_time.as_secs_f64() / overall_time.elapsed().as_secs_f64() * 100.0);
    println!("----------------------------------------------------------------------------------------------");

    // for (i,e) in basis.elements.into_iter().filter(|x| x.is_redundant == false).enumerate() {
    //     println!("lm[{}] = {:?} (#mons {} -- redundant? {})",
    //     i, hash_table.monomials[e.monomials[0]].exponents,
    //     e.monomials.len(), e.is_redundant);
    // }
}
