use num_bigint::{BigUint, ToBigUint};
use std::{collections::BTreeMap, convert::TryInto};

const BETA: i32 = 2;
const GAMMA1: i32 = 5;

// This computation counts the complete distribution of all possible
// signatures in vanilla Dilithium and tweaked Dilithium.  We iterate
// over all possible values of y1, y2, y1p, c, and cp, and compute z
// using the two different methods.  We check whether (z, c) is a valid
// signature, and if so we count *all* the possible hash functions that
// would generate that signature.
//
// View the hash function as a partially programmed random oracle.
// We program the following entries in the internal RO table, i.e.:
//
//   (y1, y2)  |-> c
//   (y1p, y2) |-> cp
//
// Let #Y be the order of the set of one y-coefficient.
// Let #B be the amount of possible c values.
//
// Now we count the complete amount of possible random oracles that
// can exist under these contraints.  There are three cases:
//
//   - Case 1 [ y1 != y1p ]: In this case, we have queried the RO
//     on two different inputs.  That means that there are two inputs
//     in the RO list.  The amount of possible inputs that are left
//     is (#Y)^2 minus two.  For each of these inputs we can choose #B
//     outputs.  So the amount of functions that we count is
//     (#B)^(#Y - 2).
//
//   - Case 2 [ y1 == y1p && c == cp ]: A collision in y1 occurs.
//     In this case, we have queried the RO twice, but with the
//     same inputs.  That means that there is only one entry in the RO
//     list.  The amount of possible inputs that can still be added is
//     (#Y)^2 minus one.
//
//   - Case 3: [ y1 == y1p && c == cp ]: A collision in y1 occurs,
//     but c != cp.  The RO is a function, so this cannot happen. We
//     do not count this case.

const ORD_B: i32 = 2 * BETA - 1;
const ORD_Y: i32 = 2 * GAMMA1;

type CounterMap = BTreeMap<(i32, i32, i32), num_bigint::BigUint>;

fn z_is_in_bounds(z: i32) -> bool {
    z.abs() < GAMMA1 - BETA
}

fn dilithium_ztrick(y1: i32, y2: i32, y1p: i32) -> CounterMap {
    assert!(-GAMMA1 < y1 && y1 <= GAMMA1);
    assert!(-GAMMA1 < y2 && y2 <= GAMMA1);
    assert!(-GAMMA1 < y1p && y1p <= GAMMA1);

    let mut zc = BTreeMap::new();
    let mut count_signature = |z1, z2, c, v| {
        let counter = zc.entry((z1, z2, c)).or_default();
        *counter += v;
    };

    // Iterate over h.  H is a tuple of two mappings from Y to B, because
    // at most we will be doing two iterations.
    // We already know the values of the inputs, they are (y1, y2) and
    // (y1p, y2), so we only need to iterate over two challenge values.
    for (c, cp) in itertools::iproduct!(-BETA + 1..BETA, -BETA + 1..BETA) {
        let mut abort = false;
        let mut z1 = y1 + c;
        let mut z2 = y2 + c;
        if !z_is_in_bounds(z1) {
            // z1 is not in bounds; resample y1.
            z1 = y1p + cp;
            z2 = y2 + cp;
            if !z_is_in_bounds(z1) || !z_is_in_bounds(z2) {
                // After two iterations, still no valid signature.
                abort = true;
            }
        } else if !z_is_in_bounds(z2) {
            // We already saw z1, so both y1 and y2 will be resampled.
            abort = true;
        }

        if abort {
            // This signature will not be output.
            continue;
        }

        if (y1, y2) != (y1p, y2) {
            // Case 1
            let programmed = (ORD_Y * ORD_Y - 2).try_into().unwrap();
            let count = ORD_B.to_biguint().unwrap().pow(programmed);
            count_signature(z1, z2, cp, count);
            continue;
        } else if (y1, y2) == (y1p, y2) && c == cp {
            // Case 2
            let programmed = ((ORD_Y * ORD_Y) - 1).try_into().unwrap();
            let count = ORD_B.to_biguint().unwrap().pow(programmed);
            count_signature(z1, z2, cp, count);
            continue;
        } else if (y1, y2) == (y1p, y2) && c != cp {
            // Case 3
            continue;
        } else {
            unreachable!();
        }
    }

    zc
}

fn dilithium_vanilla(y1: i32, y2: i32, y1p: i32) -> CounterMap {
    assert!(-GAMMA1 < y1 && y1 <= GAMMA1);
    assert!(-GAMMA1 < y2 && y2 <= GAMMA1);
    assert!(-GAMMA1 < y1p && y1p <= GAMMA1);

    let mut zc = BTreeMap::new();
    let mut count_signature = |z1, z2, c, v: BigUint| {
        let counter = zc.entry((z1, z2, c)).or_default();
        *counter += v;
    };

    // Iterate over h.  H is a tuple of two mappings from Y to B, because
    // at most we will be doing two iterations.
    // We already know the values of the inputs, they are (y1, y2) and
    // (y1p, y2), so we only need to iterate over two challenge values.
    for (c, cp) in itertools::iproduct!(-BETA + 1..BETA, -BETA + 1..BETA) {
        let mut abort = false;
        let z1 = y1 + c;
        let z2 = y2 + c;
        let _z1p = y1p + cp;
        let _z2p = y2 + cp;
        if !z_is_in_bounds(z1) || !z_is_in_bounds(z2) {
            abort = true;
        }

        if abort {
            // This signature will not be output.
            continue;
        }
        
        if (y1, y2) != (y1p, y2) {
            // Case 1.
            let programmed = (ORD_Y * ORD_Y) as u32 - 2;
            let count = ORD_B.to_biguint().unwrap().pow(programmed);
            count_signature(z1, z2, cp, count);
        } else if (y1, y2) == (y1p, y2) && c == cp {
            // Case 2.
            let programmed = (ORD_Y * ORD_Y) as u32 - 1;
            let count = ORD_B.to_biguint().unwrap().pow(programmed);
            count_signature(z1, z2, cp, count);
            continue;
        } else if (y1, y2) == (y1p, y2) && c != cp {
            // Case 3.
            continue;
        } else {
            unreachable!();
        }
    }

    zc
}

fn is_uniform<I, K, V>(iter: &mut I) -> bool
where
    I: Iterator<Item = (K, V)>,
    K: std::fmt::Debug,
    V: std::fmt::Debug + PartialEq,
{
    let (_, baseline) = iter.next().expect("empty map");
    while let Some((k, v)) = iter.next() {
        if v != baseline {
            eprintln!(
                "value at {:?} ({:?}) is not same as baseline ({:?})",
                k, v, baseline
            );
            return false;
        }
    }
    true
}

fn merge_results(map: &mut CounterMap, other: &CounterMap) {
    for (k, v) in other.iter() {
        let entry = map.entry(*k).or_default();
        *entry += v;
    }
}

fn main() {
    let mut vanilla_results = BTreeMap::<_, BigUint>::new();
    for (y1, y2, y1p) in itertools::iproduct!(
        (-GAMMA1 + 1..=GAMMA1),
        (-GAMMA1 + 1..=GAMMA1),
        (-GAMMA1 + 1..=GAMMA1)
    ) {
        merge_results(&mut vanilla_results, &dilithium_vanilla(y1, y2, y1p));
    }

    let mut ztrick_results = BTreeMap::new();
    for (y1, y2, y1p) in itertools::iproduct!(
        (-GAMMA1 + 1..=GAMMA1),
        (-GAMMA1 + 1..=GAMMA1),
        (-GAMMA1 + 1..=GAMMA1)
    ) {
        merge_results(&mut ztrick_results, &dilithium_ztrick(y1, y2, y1p));
    }

    dbg!(GAMMA1, BETA);
    eprintln!("vanilla:");
    dbg!(is_uniform(&mut vanilla_results.iter()));

    eprintln!("ztrick:");
    dbg!(is_uniform(&mut ztrick_results.iter()));

    // Print the amount of possible signatures for this tuple.
    eprintln!("(z1, z2, c): [vanilla], [ztrick]");
    for k in vanilla_results.keys() {
        eprintln!(
            "{:?}: {}, {}",
            k,
            vanilla_results.get(k).map(|x| x).unwrap(),
            ztrick_results.get(k).map(|x| x).unwrap()
        );
    }
}
