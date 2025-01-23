// Exponent stuff
pub type Exponent = u32;
pub type ExpVec = Vec<Exponent>;

// Hash table stuff
pub type HashValue = usize;
pub type HashTableLength = usize;
pub type DivisorMask = usize;
pub type MonomVec = Vec<HashTableLength>;

// Field stuff
pub type Characteristic = i32;
pub type Coefficient = Characteristic;
pub type CoeffVec = Vec<Coefficient>;

// Basis stuff
pub type Degree = u32;
pub type BasisLength = usize;

// Pair stuff
pub enum Criterion {
    Keep,
    Chain,
    Product,
}

pub struct Pair {
    pub lcm: HashTableLength,
    pub generators: (BasisLength, BasisLength),
    pub criterion: Criterion,
}

pub type PairVec = Vec<Pair>;

// Meta data stuff
pub struct MetaData {
    pub characteristic: Characteristic,
    pub nr_pairs_reduced: usize,
    pub nr_redundant_elements: usize,
    pub nr_input_generators: usize,
}
