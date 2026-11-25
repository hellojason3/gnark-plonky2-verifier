use ark_ff::fields::{Field, Fp256, MontBackend, MontConfig};
use ark_ff::BigInt;
use lazy_static::lazy_static;

use crate::fr_bn254::FrBn254;
use crate::poseidon_bn254_constants::C_CONSTANTS;
use crate::poseidon_bn254_constants::M_MATRIX;
use crate::poseidon_bn254_constants::P_MATRIX;
use crate::poseidon_bn254_constants::S_CONSTANTS;

pub const RATE: usize = 3;
pub const WIDTH: usize = 4;
pub const FULL_ROUNDS: usize = 8;
pub const PARTIAL_ROUNDS: usize = 56;
pub const GOLDILOCKS_ELEMENTS: usize = 3;

pub type PoseidonState = [FrBn254; WIDTH];

#[derive(MontConfig)]
#[modulus = "21888242871839275222246405745257275088548364400416034343698204186575808495617"]
#[generator = "5"]
struct ArithmeticConfig;

/// The same field as `FrBn254` with faster multiplication and squaring. Both
/// hold the Montgomery form (R = 2^256) in four little-endian limbs.
type Element = Fp256<MontBackend<ArithmeticConfig, 4>>;

fn to_element(value: &FrBn254) -> Element {
    Element::new_unchecked(BigInt(value.montgomery_limbs()))
}

fn from_element(value: &Element) -> FrBn254 {
    FrBn254::from_montgomery_limbs(value.0 .0)
}

struct Constants {
    round: Vec<Element>,
    sparse: Vec<Element>,
    mds: [[Element; WIDTH]; WIDTH],
    pre_sparse: [[Element; WIDTH]; WIDTH],
}

fn to_matrix(rows: &[Vec<FrBn254>]) -> [[Element; WIDTH]; WIDTH] {
    core::array::from_fn(|j| core::array::from_fn(|i| to_element(&rows[j][i])))
}

lazy_static! {
    static ref CONSTANTS: Constants = Constants {
        round: C_CONSTANTS.iter().map(to_element).collect(),
        sparse: S_CONSTANTS.iter().map(to_element).collect(),
        mds: to_matrix(&M_MATRIX),
        pre_sparse: to_matrix(&P_MATRIX),
    };
}

pub fn permution(state: &mut PoseidonState) {
    let constants: &Constants = &CONSTANTS;
    let mut elements: [Element; WIDTH] = core::array::from_fn(|i| to_element(&state[i]));
    let half = FULL_ROUNDS / 2;
    let partial_start = (half + 1) * WIDTH;
    let tail_start = partial_start + PARTIAL_ROUNDS;

    add_constants(&mut elements, &constants.round[..WIDTH]);
    for i in 0..half - 1 {
        full_round(&mut elements, &constants.round[(i + 1) * WIDTH..(i + 2) * WIDTH], &constants.mds);
    }
    full_round(&mut elements, &constants.round[half * WIDTH..partial_start], &constants.pre_sparse);

    let partial_constants = &constants.round[partial_start..tail_start];
    for (constant, sparse) in partial_constants.iter().zip(constants.sparse.chunks_exact(WIDTH * 2 - 1)) {
        let first = exp5(&elements[0]) + constant;
        let mut new_first = sparse[0] * first;
        for j in 1..WIDTH {
            new_first += sparse[j] * elements[j];
        }
        for k in 1..WIDTH {
            elements[k] += first * sparse[WIDTH + k - 1];
        }
        elements[0] = new_first;
    }

    for i in 0..half - 1 {
        full_round(&mut elements, &constants.round[tail_start + i * WIDTH..tail_start + (i + 1) * WIDTH], &constants.mds);
    }
    full_round(&mut elements, &[], &constants.mds);

    for (out, element) in state.iter_mut().zip(&elements) {
        *out = from_element(element);
    }
}

fn add_constants(state: &mut [Element; WIDTH], round_constants: &[Element]) {
    for (element, constant) in state.iter_mut().zip(round_constants) {
        *element += constant;
    }
}

fn exp5(x: &Element) -> Element {
    x.square().square() * x
}

fn full_round(state: &mut [Element; WIDTH], round_constants: &[Element], matrix: &[[Element; WIDTH]; WIDTH]) {
    for element in state.iter_mut() {
        *element = exp5(element);
    }
    add_constants(state, round_constants);
    let input = *state;
    for (i, out) in state.iter_mut().enumerate() {
        let mut sum = matrix[0][i] * input[0];
        for j in 1..WIDTH {
            sum += matrix[j][i] * input[j];
        }
        *out = sum;
    }
}

#[cfg(test)]
mod tests {
    use std::ops::{AddAssign, MulAssign};

    use ff::Field as _;

    use super::*;

    /// The permutation as first written, on the ff-derived field.
    fn permution_reference(state: &mut PoseidonState) {
        fn ark(state: &mut PoseidonState, it: usize) {
            for i in 0..WIDTH {
                state[i].add_assign(&C_CONSTANTS[it + i]);
            }
        }
        fn exp5(mut x: FrBn254) -> FrBn254 {
            let aux = x;
            x = x.square();
            x = x.square();
            x.mul_assign(&aux);
            x
        }
        fn mix(state: &mut PoseidonState, constant_matrix: &[Vec<FrBn254>]) {
            let mut result: PoseidonState = [FrBn254::ZERO; WIDTH];
            for (i, result_element) in result.iter_mut().enumerate() {
                for j in 0..WIDTH {
                    let mut mul = constant_matrix[j][i];
                    mul.mul_assign(&state[j]);
                    result_element.add_assign(&mul);
                }
            }
            *state = result;
        }
        fn full_rounds(state: &mut PoseidonState, first: bool) {
            for i in 0..FULL_ROUNDS / 2 - 1 {
                for element in state.iter_mut() {
                    *element = exp5(*element);
                }
                if first {
                    ark(state, (i + 1) * WIDTH);
                } else {
                    ark(state, (FULL_ROUNDS / 2 + 1) * WIDTH + PARTIAL_ROUNDS + i * WIDTH);
                }
                mix(state, &M_MATRIX);
            }
            for element in state.iter_mut() {
                *element = exp5(*element);
            }
            if first {
                ark(state, (FULL_ROUNDS / 2) * WIDTH);
                mix(state, &P_MATRIX);
            } else {
                mix(state, &M_MATRIX);
            }
        }
        ark(state, 0);
        full_rounds(state, true);
        for i in 0..PARTIAL_ROUNDS {
            state[0] = exp5(state[0]);
            state[0].add_assign(&C_CONSTANTS[(FULL_ROUNDS / 2 + 1) * WIDTH + i]);
            let mut new_state0 = FrBn254::ZERO;
            for j in 0..WIDTH {
                let mut mul = S_CONSTANTS[(WIDTH * 2 - 1) * i + j];
                mul.mul_assign(&state[j]);
                new_state0.add_assign(&mul);
            }
            for k in 1..WIDTH {
                let mut mul = state[0];
                mul.mul_assign(&S_CONSTANTS[(WIDTH * 2 - 1) * i + WIDTH + k - 1]);
                state[k].add_assign(&mul);
            }
            state[0] = new_state0;
        }
        full_rounds(state, false);
    }

    #[test]
    fn permutation_equals_reference_on_chained_states() {
        let mut state = [FrBn254::from(1u64), FrBn254::from(2u64), FrBn254::from(3u64), FrBn254::from(4u64)];
        for round in 0..20_000u64 {
            let mut expected = state;
            permution_reference(&mut expected);
            permution(&mut state);
            assert!(state == expected, "permutation differs at round {round}");
            state[(round % 4) as usize] += FrBn254::from(round).square();
        }
    }

    #[test]
    fn permutation_equals_reference_on_edge_states() {
        let top = -FrBn254::ONE;
        for mut state in [[FrBn254::ZERO; WIDTH], [top; WIDTH], [FrBn254::ONE, top, FrBn254::ZERO, top]] {
            let mut expected = state;
            permution_reference(&mut expected);
            permution(&mut state);
            assert!(state == expected);
        }
    }
}
