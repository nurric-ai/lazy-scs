"""Stage 2b: the 2-approximation driven by compressed counts (lazy rules: extension
rules only so far; raises NeedRule if Section 5.3 needs a periodic rule)."""
from collections import Counter
import scs2
from scs2 import (Walks, Connector, closed_walks, build_layers, euler_superstring,
                  connected_to_eps, reduce_inputs, Result)
from compressed import CompressedCounts


def two_approx_c(strings, texts=()):
    req = reduce_inputs(strings)
    if not req:
        return Result("", 0, {})
    cc = CompressedCounts(req, texts)
    V = cc.sam                                  # membership via the automaton
    up, down = cc.edges()
    walks = closed_walks(up, down)
    groups, layers = build_layers(walks)
    stats = Counter()
    added = Walks(V)
    Connector(groups, layers, added, V, stats).process()
    assert added.cost <= cc.W
    U, D = up + added.up, down + added.down
    assert connected_to_eps(U, D)
    T = euler_superstring(U, D)
    assert len(T) <= 2 * cc.W and all(s in T for s in req)
    stats.update(groups=len(groups), layers=len(layers), states=cc.states)
    return Result(T, cc.W, dict(stats))


def two_approx_c_lazy(strings, max_rounds=1000):
    """Compressed counts + periodic rules enabled per text only when Section 5.3 needs them."""
    texts = set()
    for rnd in range(max_rounds):
        try:
            r = two_approx_c(strings, frozenset(texts))
            r.stats.update(lazy_texts=len(texts), lazy_rounds=rnd + 1)
            return r
        except scs2.NeedRule as e:
            v, p, k = e.args[0]
            key = (p, scs2.least_rotation(v[:p]))
            assert key not in texts, ("rule enabled but still missing", v, p, k)
            texts.add(key)
    raise RuntimeError("no convergence")
