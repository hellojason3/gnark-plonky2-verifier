use ff::PrimeField;

#[derive(PrimeField)]
#[PrimeFieldModulus = "21888242871839275222246405745257275088548364400416034343698204186575808495617"]
#[PrimeFieldGenerator = "7"]
#[PrimeFieldReprEndianness = "little"]
pub struct FrBn254([u64; 4]);

impl FrBn254 {
    /// The Montgomery form (R = 2^256), least significant limb first.
    pub(crate) fn montgomery_limbs(&self) -> [u64; 4] {
        self.0
    }

    pub(crate) fn from_montgomery_limbs(limbs: [u64; 4]) -> Self {
        Self(limbs)
    }
}
