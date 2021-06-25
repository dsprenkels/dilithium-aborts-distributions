use std::collections::BTreeMap;

const BETA: i64 = 2;
const GAMMA1: i64 = 7;
const ORD_B: i64 = 2 * BETA - 1;
const ORD_Y: i64 = 2 * GAMMA1;

type CounterMap = BTreeMap<(i64, i64, i64), i64>;

fn z_is_in_bounds(z: i64) -> bool {
    z.abs() < GAMMA1 - BETA
}

fn dilithium_ztrick(y1: i64, y2: i64, y1p: i64) -> CounterMap {
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
            z1 = y1p + cp;
            z2 = y2 + cp;
            if !z_is_in_bounds(z1) || !z_is_in_bounds(z2) {
                abort = true;
            }
        } else if !z_is_in_bounds(z2) {
            abort = true;
        }

        if abort {
            continue;
        } else if (y1, y2) != (y1p, y2) {
            count_signature(z1, z2, cp, ORD_B.pow(ORD_Y as u32 - 2));
        } else if (y1, y2) == (y1p, y2) && c == cp {
            count_signature(z1, z2, cp, ORD_B.pow(ORD_Y as u32 - 1));
            continue;
        } else if (y1, y2) == (y1p, y2) && c != cp {
            continue;
        } else {
            unreachable!();
        }
    }

    zc
}

fn dilithium_vanilla(y1: i64, y2: i64, y1p: i64) -> CounterMap {
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
        let z1 = y1 + c;
        let z2 = y2 + c;
        let _z1p = y1p + cp;
        let _z2p = y2 + cp;
        if !z_is_in_bounds(z1) || !z_is_in_bounds(z2) {
            abort = true;
        }

        if abort {
            continue;
        } else if (y1, y2) != (y1p, y2) {
            count_signature(z1, z2, cp, ORD_B.pow(ORD_Y as u32 - 2));
        } else if (y1, y2) == (y1p, y2) && c == cp {
            count_signature(z1, z2, cp, ORD_B.pow(ORD_Y as u32 - 1));
            continue;
        } else if (y1, y2) == (y1p, y2) && c != cp {
            continue;
        } else {
            unreachable!();
        }
    }

    zc
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
    let mut vanilla_results = BTreeMap::<_, i64>::new();
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

    dbg!(&ztrick_results);

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
