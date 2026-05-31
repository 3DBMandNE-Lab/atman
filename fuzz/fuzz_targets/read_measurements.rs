#![no_main]
//! Fuzz the canonical measurements TSV reader. The reader must never panic,
//! index out of bounds, or hang on arbitrary bytes — only return `Ok`/`Err`.
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
    let _ = atman::io::read_measurements_long(tmp.path());
});
