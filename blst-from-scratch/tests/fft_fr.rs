#[cfg(test)]
mod tests {
    use blst_from_scratch::fft_fr::{fft_fr_fast, fft_fr_slow};
    use blst_from_scratch::types::fft_settings::FsFFTSettings;
    use blst_from_scratch::types::fr::FsFr;
    use kzg_bench::tests::fft_fr::{compare_sft_fft, inverse_fft, roundtrip_fft, stride_fft};

    #[test]
    fn compare_sft_fft_() {
        compare_sft_fft::<FsFr, FsFFTSettings>(&fft_fr_slow, &fft_fr_fast);
    }

    #[test]
    fn roundtrip_fft_() {
        roundtrip_fft::<FsFr, FsFFTSettings>();
    }

    #[test]
    fn inverse_fft_() {
        inverse_fft::<FsFr, FsFFTSettings>();
    }

    #[test]
    fn stride_fft_() {
        stride_fft::<FsFr, FsFFTSettings>();
    }
}
