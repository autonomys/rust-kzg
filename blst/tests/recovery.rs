 #[path = "./local_tests/local_recovery.rs"]
 pub mod local_recovery;

#[cfg(test)]
mod tests {
    use kzg_bench::tests::recover::*;
    // uncomment to use the local tests
    use crate::local_recovery::{recover_random, recover_simple, recover_g1};

    use blst_rust::types::fft_settings::FsFFTSettings;
    use blst_rust::types::fr::FsFr;
    use blst_rust::types::g1::FsG1;
    use blst_rust::types::poly::{FsPoly, FsPolyG1};

    // Shared tests
    #[test]
    fn recover_simple_() {
        recover_simple::<FsFr, FsFFTSettings, FsPoly, FsPoly>();
    }

    #[test]
    fn recover_random_() {
        recover_random::<FsFr, FsFFTSettings, FsPoly, FsPoly>();
    }

    #[test]
    fn more_than_half_missing_() {
        more_than_half_missing::<FsFr, FsFFTSettings, FsPoly, FsPoly>();
    }

    #[test]
    fn recover_g1_() {
        recover_g1::<FsG1, FsFr, FsFFTSettings, FsPolyG1, FsPolyG1>()
    }
}
