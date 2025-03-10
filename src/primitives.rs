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

// Matrix stuff
pub type DenseRowCoefficient = i128;
pub type DenseRow = Vec<DenseRowCoefficient>;
