//! Python bindings (`lazy_scs._core`), built by maturin with the `python` feature.

use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use pyo3::types::PyBytes;
use std::collections::BTreeMap;

type SolveResult<'py> = (Bound<'py, PyBytes>, u64, u64, usize, BTreeMap<String, u64>);

/// solve(strings: list[bytes], check=False, verify=False)
///     -> (superstring, W, W_ext, rounds, stats)
///
/// W is the certified lower bound W' <= OPT; the superstring has length <= 2W.
/// `verify` re-justifies every count from scratch (adds verify_* entries to stats).
/// The GIL is released while solving.
#[pyfunction]
#[pyo3(signature = (strings, check=false, verify=false))]
fn solve<'py>(py: Python<'py>, strings: Vec<Vec<u8>>, check: bool, verify: bool) -> PyResult<SolveResult<'py>> {
    let (r, ok) = py.allow_threads(|| {
        let r = crate::solve(strings.clone(), check, verify);
        let ok = crate::covers(&r.superstring, &strings);
        (r, ok)
    });
    if !ok {
        return Err(PyRuntimeError::new_err("internal error: the output misses an input string"));
    }
    let stats = r.stats.iter().map(|(k, v)| (k.to_string(), *v)).collect();
    Ok((PyBytes::new(py, &r.superstring), r.w, r.w_ext, r.rounds, stats))
}

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(solve, m)?)?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
