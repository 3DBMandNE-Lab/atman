#![no_main]
//! Fuzz the samples TSV reader. Must never panic on arbitrary bytes.
use libfuzzer_sys::fuzz_target;
use std::io::Write;

fuzz_target!(|data: &[u8]| {
    let mut tmp = match tempfile::NamedTempFile::new() {
        Ok(t) => t,
        Err(_) => return,
    };
    if tmp.write_all(data).is_err() {
        return;
    }
    let _ = atman::io::read_samples(tmp.path());
});
