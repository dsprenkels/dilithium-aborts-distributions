use std::collections::BTreeMap;

const BETA: i32 = 2;
const GAMMA1: i32 = 4;

type CounterMap = BTreeMap<(i32, i32, i32), i32>;

fn z_is_in_bounds(z: i32) -> bool {
    z.abs() < GAMMA1 - BETA
}

fn dilithium_ztrick(y1: i32, y2: i32, y1p: i32) -> CounterMap {
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
    for (c, cp) in itertools::iproduct!(-BETA + 1..BETA, -BETA + 1..BETA) {
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

            // z-check from the second iteration.
            if z_is_in_bounds(z1) && z_is_in_bounds(z2) {
                // Good signature in two iterations.  Ouput.
                count_signature(z1, z2, cp);
            }
        } else if !z_is_in_bounds(z2) {
            // Abort completely.
            continue;
        } else {
            // Good singature in one iteration.  Output.
            count_signature(z1, z2, c);
            continue;
        }
    }

    return zc;
}

fn dilithium_vanilla(y1: i32, y2: i32, y1p: i32) -> CounterMap {
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
    for (c, cp) in itertools::iproduct!(-BETA + 1..BETA, -BETA + 1..BETA) {
        let z1 = y1 + c;
        let z2 = y2 + c;
        let _z1p = y1p + cp;
        let _z2p = y2 + cp;
        if z_is_in_bounds(z1) && z_is_in_bounds(z2) {
            // Good singature in one iteration.  Output.
            count_signature(z1, z2, c);
            continue;
        }
    }

    return zc;
}

fn is_uniform(map: &CounterMap) -> bool {
    let mut iter = map.iter();
    let (_, baseline) = iter.next().expect("empty map");
    for (k, v) in iter {
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
    let mut vanilla_results = BTreeMap::new();
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
    dbg!(vanilla_results.len());
    // dbg!(&vanilla_results);
    dbg!(is_uniform(&vanilla_results));

    eprintln!("ztrick:");
    dbg!(ztrick_results.len());
    // dbg!(&ztrick_results);
    dbg!(is_uniform(&ztrick_results));

    eprintln!("(z1, z2, c): [vanilla], [ztrick]");
    for k in vanilla_results.keys() {
        eprintln!(
            "{:?}: {}, {}",
            k,
            vanilla_results.get(k).map(|x| *x).unwrap_or_default(),
            ztrick_results.get(k).map(|x| *x).unwrap_or_default()
        );
    }
}
