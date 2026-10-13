"""lazy-scs: scalable factor-2 approximation for the shortest common superstring.

    >>> import lazy_scs
    >>> lazy_scs.shortest_superstring(["GATTACA", "TACAGG", "AGGCT"])
    'GATTACAGGCT'
    >>> r = lazy_scs.solve(["GATTACA", "TACAGG", "AGGCT"])
    >>> r.length, r.lower_bound, r.ratio
    (11, 11, 1.0)

The superstring always has length <= 2 * lower_bound <= 2 * OPT, and `ratio` = length / lower_bound
is a per-instance certificate: ratio 1.0 means the output is provably optimal.

Strings may be `str` (encoded as UTF-8) or `bytes`. The result has the same type as the inputs.
The core is written in Rust and releases the GIL while it runs.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Sequence, Union

from ._core import __version__
from ._core import solve as _solve

__all__ = ["Result", "solve", "shortest_superstring", "__version__"]

StrOrBytes = Union[str, bytes]


@dataclass(frozen=True)
class Result:
    superstring: StrOrBytes
    """A superstring of every input (same type as the inputs)."""
    lower_bound: int
    """Certified lower bound W' on the optimal length: lower_bound <= OPT <= len(superstring)."""
    rounds: int
    """Solver rounds (1 + number of times a lazy periodic rule had to be enabled)."""
    stats: dict = field(default_factory=dict)
    """Internal counters (and verify_* entries when solve(..., verify=True))."""

    @property
    def length(self) -> int:
        return len(self.superstring)

    @property
    def ratio(self) -> float:
        """length / lower_bound >= length / OPT: how far from optimal the output can be, at most."""
        return self.length / self.lower_bound if self.lower_bound else 1.0


def _encode(strings: Iterable[StrOrBytes]) -> tuple[list[bytes], bool]:
    items = list(strings)
    if not items:
        return [], True
    kinds = {type(s) for s in items}
    if kinds <= {bytes, bytearray, memoryview}:
        return [bytes(s) for s in items], False
    if kinds <= {str}:
        return [s.encode("utf-8") for s in items], True
    raise TypeError("strings must be all str or all bytes, got " + ", ".join(sorted(k.__name__ for k in kinds)))


def solve(strings: Sequence[StrOrBytes], *, verify: bool = False, check: bool = False) -> Result:
    """Run the factor-2 approximation and return the superstring with its certified lower bound.

    verify: re-justify every count of the lower bound from scratch (independent checker; slower).
    check:  extra internal consistency checks while building (for debugging).
    """
    data, is_text = _encode(strings)
    sup, w, _w_ext, rounds, stats = _solve(data, check=check, verify=verify)
    if verify and (stats.get("verify_unjustified") or stats.get("verify_ext_violations")):
        raise RuntimeError(f"lower-bound verification failed: {stats}")
    if is_text:
        try:
            sup = sup.decode("utf-8")
        except UnicodeDecodeError as e:   # overlaps are found byte-wise; only possible with non-ASCII text
            raise ValueError("the superstring of the UTF-8 encoded inputs is not valid UTF-8; pass bytes "
                             "instead (e.g. [s.encode() for s in strings])") from e
    return Result(superstring=sup, lower_bound=w, rounds=rounds, stats=dict(stats))


def shortest_superstring(strings: Sequence[StrOrBytes]) -> StrOrBytes:
    """A superstring of all `strings` with length at most twice the optimum."""
    return solve(strings).superstring
