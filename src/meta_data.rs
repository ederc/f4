use crate::types::*;

pub struct MetaData {
    pub characteristic: Characteristic,
    pub nr_pairs_reduced: usize,
    pub nr_redundant_elements: usize,
    pub nr_input_generators: usize,
}
