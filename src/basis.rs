use crate::types::*;

pub struct Basis {
    previous_length: BasisLength,
    is_constant: bool,
    maximum_total_degree: Degree,
    is_redundant: Vec<bool>,
    monomials: Vec<MonomVec>,
    coefficients: Vec<CoeffVec>,
}

// impl<'a> Basis<'a> {
//     pub fn new<T>(coefficients:
//         let len = coefficients.len();
//         let mut basis = Basis {
//             previous_length      = 0,
//             is_constant          = false,
//             maximum_total_degree = 0,

