use crate::primitives::*;

pub struct MetaData {
    pub characteristic: Characteristic,
    pub nr_pairs_reduced: usize,
    pub nr_redundant_elements: usize,
    pub nr_input_generators: usize,
}

impl MetaData {
    pub fn new(
        characteristic: Characteristic,
        nr_input_generators: usize) -> MetaData {
        
        return MetaData {
                characteristic : characteristic,
                nr_pairs_reduced : 0,
                nr_redundant_elements : 0,
                nr_input_generators : nr_input_generators,
                };
    }
}
