use std::collections::BTreeMap;

const BETA: i32 = 10;
const GAMMA1: i32 = 100;

fn z_is_in_bounds(z: i32) -> bool {
    z.abs() < GAMMA1 - BETA
}

fn dilithium_ztrick(y1: i32, y2: i32, y1p: i32) -> BTreeMap<(i32, i32, i32), u8> {
    assert!(-GAMMA1 < y1 && y1 <= GAMMA1);
    assert!(-GAMMA1 < y2 && y2 <= GAMMA1);
    assert!(-GAMMA1 < y1p && y1p <= GAMMA1);

    let mut zc = BTreeMap::new();
    let mut count_signature = |z1, z2, c| {
        let counter = zc.entry((z1, z2, c)).or_default();
        *counter += 1;
    };

    // Iterate over h.  H is a tuple of two mappings from Y to B, because
    // at most we will be doing two iterations.
    // We already know the values of the inputs, they are (y1, y2) and
    // (y1p, y2), so we only need to iterate over two challenge values.
    for (c, cp) in itertools::iproduct!(-BETA..=BETA, -BETA..=BETA) {
        let mut z1 = y1 + c;
        let mut z2 = y2 + c;
        if !z_is_in_bounds(z1) {
            // Sample y1 again.
            if (y1, y2) == (y1p, y2) {
                // Hash function has to return the same challenge.
                z1 = y1p + c;
            } else {
                // Hash function just maps to a new random value.
                z1 = y1p + cp;
                z2 = y2 + cp;
            }
        } else if !z_is_in_bounds(z2) {
            // Abort completely.
            continue;
        } else {
            // Good singature in one iteration.  Output.
            count_signature(z1, z2, c);
            continue;
        }

        // z-check from the second iteration.
        if !z_is_in_bounds(z1) {
            // Abort completely.
            continue;
        }
        if !z_is_in_bounds(z2) {
            // Abort completely.
            continue;
        }

        // Good signature in two iterations.  Ouput.
        count_signature(z1, z2, cp);
    }

    return zc;
}

fn is_uniform(map: &BTreeMap<(i32, i32, i32), u8>) -> bool {
    let baseline = map.values().next().expect("empty map");

    let z1_range = -GAMMA1 + BETA + 1..GAMMA1 - BETA;
    let z2_range = -GAMMA1 + BETA + 1..GAMMA1 - BETA;
    let c_range = -BETA..=BETA;

    for key in itertools::iproduct!(z1_range, z2_range, c_range) {
        if map.get(&key) != Some(baseline) {
            eprintln!(
                "value at {:?} ({:?}) is not same as baseline ({:?})",
                key,
                map.get(&key),
                baseline
            );
            return false;
        }
    }
    true
}

fn main() {
    let mut results = BTreeMap::new();
    for (y1, y2, y1p) in itertools::iproduct!(
        (-GAMMA1 + 1..=GAMMA1),
        (-GAMMA1 + 1..=GAMMA1),
        (-GAMMA1 + 1..=GAMMA1)
    ) {
        results.extend(dilithium_ztrick(y1, y2, y1p));
    }

    dbg!(results.len());
    dbg!(&results);
    dbg!(is_uniform(&results));
}
