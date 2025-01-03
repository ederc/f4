use crate::types::*;

pub struct Basis<T> {
    previous_length: BasisLength,
    is_constant: bool,
    maximum_total_degree: Degree,
    is_redundant: Vec<bool>,
    monomials: Vec<Vec<Vec<Exponent>>>,
    coefficients: Vec<Vec<T>>,
}
