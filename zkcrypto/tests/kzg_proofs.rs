#[cfg(test)]
mod tests {
    use kzg_bench::tests::kzg_proofs::*;
    use zkcrypto::fftsettings::ZkFFTSettings;
    use zkcrypto::kzg_proofs::{generate_trusted_setup, KZGSettings};
    use zkcrypto::kzg_types::{ZkG1Projective, ZkG2Projective};
    use zkcrypto::poly::ZPoly;
    use zkcrypto::zkfr::blsScalar;

    #[test]
    fn test_proof_single() {
        proof_single::<blsScalar, ZkG1Projective, ZkG2Projective, ZPoly, ZkFFTSettings, KZGSettings>(
            &generate_trusted_setup,
        );
    }

    #[test]
    fn test_commit_to_nil_poly() {
        commit_to_nil_poly::<
            blsScalar,
            ZkG1Projective,
            ZkG2Projective,
            ZPoly,
            ZkFFTSettings,
            KZGSettings,
        >(&generate_trusted_setup);
    }

    #[test]
    #[should_panic(expected = "Poly given is too long")]
    fn test_commit_to_too_long_poly() {
        commit_to_too_long_poly::<
            blsScalar,
            ZkG1Projective,
            ZkG2Projective,
            ZPoly,
            ZkFFTSettings,
            KZGSettings,
        >(&generate_trusted_setup);
    }

    // #[test]
    // fn commit_to_too_long_poly_returns_err_() {
    // commit_to_too_long_poly_returns_err::<blsScalar, ZkG1Projective, ZkG2Projective, ZPoly, ZkFFTSettings, KZGSettings>(&generate_trusted_setup);
    // }

    #[test]
    fn test_proof_multi() {
        proof_multi::<blsScalar, ZkG1Projective, ZkG2Projective, ZPoly, ZkFFTSettings, KZGSettings>(
            &generate_trusted_setup,
        );
    }
}
