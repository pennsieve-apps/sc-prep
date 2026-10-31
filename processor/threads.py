"""
The thread count, fixed once at start-up from THREADS (0 = the cgroup CPU quota, via
safe_inputs.cpu_threads) and set as the BLAS / OpenMP / numba environment variables. Import this
before numpy: the libraries read these variables when they load, so a later change (e.g. a Lambda
event's THREADS) would not take effect and is not attempted; on Lambda the count is the function's
vCPUs.
"""
from __future__ import annotations

import os

from processor.safe_inputs import cpu_threads

VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS", "NUMEXPR_NUM_THREADS")


def configure(value: str | None) -> int:
    """Set the thread variables from THREADS (an invalid value counts as 0; main() rejects it with exit 2)."""
    try:
        n = cpu_threads(int(value or 0))
    except ValueError:
        n = cpu_threads(0)
    for v in VARS:
        os.environ[v] = str(n)
    return n


THREADS = configure(os.environ.get("THREADS"))
