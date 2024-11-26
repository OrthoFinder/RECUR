import multiprocessing as mp
from typing import Tuple

try:
    from importlib.metadata import PackageNotFoundError, version
    __version__ = version(__name__)
except PackageNotFoundError:
    from ._version import __version__ as __version__


def find_balanced_pair(n: int) -> Tuple[int, int]:
    best_pair = (1, n)
    min_difference = abs(n - 1)

    for i in range(1, int(n**0.5) + 1):
        if n % i == 0:
            pair = (i, n // i)
            difference = abs(pair[0] - pair[1])
            if difference < min_difference:
                min_difference = difference
                best_pair = pair

    return best_pair


recur_iqtree_nthreads, iqtree_nthreads = find_balanced_pair(mp.cpu_count())
recur_nthreads = mp.cpu_count()
