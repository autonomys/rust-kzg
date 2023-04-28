use criterion::{criterion_group, criterion_main, Criterion};
use kzg_bench::benches::kzg::{bench_commit_to_poly, bench_compute_proof_single};

use zkcrypto::fftsettings::ZkFFTSettings;
use zkcrypto::kzg_proofs::{generate_trusted_setup, KZGSettings};
use zkcrypto::kzg_types::{ZkG1Projective, ZkG2Projective};
use zkcrypto::poly::ZPoly;
use zkcrypto::zkfr::blsScalar;

fn bench_commit_to_poly_(c: &mut Criterion) {
    bench_commit_to_poly::<
        blsScalar,
        ZkG1Projective,
        ZkG2Projective,
        ZPoly,
        ZkFFTSettings,
        KZGSettings,
    >(c, &generate_trusted_setup)
}

fn bench_compute_proof_single_(c: &mut Criterion) {
    bench_compute_proof_single::<
        blsScalar,
        ZkG1Projective,
        ZkG2Projective,
        ZPoly,
        ZkFFTSettings,
        KZGSettings,
    >(c, &generate_trusted_setup)
}

criterion_group! {
    name = benches;
    config = Criterion::default().sample_size(10);
    targets = bench_commit_to_poly_, bench_compute_proof_single_
}

criterion_main!(benches);
