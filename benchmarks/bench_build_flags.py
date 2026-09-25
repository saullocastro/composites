r"""Computational cost of the laminate calculations, to assess the build flags

Times the creation of a laminate with ``laminated_plate`` and the methods and
functions that recompute its stiffness, for laminates with a different number
of plies. It is used to assess the Cython directives and the compiler flags of
``setup.py``, see ``CHANGELOG.md``.

Run it against the build of composites to be assessed, e.g.::

    PYTHONPATH=<path to composites repository> python bench_build_flags.py

"""
import json
import sys
import time

import numpy as np

import composites
from composites import laminated_plate
from composites.core import (laminate_from_LaminationParameters,
                             n_double_laminate, GradABD)

REPEAT = 5
CFRP = (138e9, 9.3e9, 0.3, 4.6e9, 4.6e9, 2.3e9)
STACKS = {
    '8 plies': [0, 45, -45, 90, 90, -45, 45, 0],
    '32 plies': [0, 45, -45, 90]*8,
    '128 plies': [0, 45, -45, 90]*32,
}


def pin_process():
    # NOTE a single core and a high priority reduce the timing noise
    try:
        import psutil
        p = psutil.Process()
        p.cpu_affinity([p.cpu_affinity()[-1]])
        if hasattr(psutil, 'HIGH_PRIORITY_CLASS'):
            p.nice(psutil.HIGH_PRIORITY_CLASS)
        else:
            p.nice(-10)
    except Exception:
        pass


def timeit(func, number):
    best = np.inf
    for _ in range(REPEAT):
        t0 = time.perf_counter()
        for _ in range(number):
            func()
        best = min(best, (time.perf_counter() - t0)/number)
    return best


def main():
    pin_process()
    results = {'composites': composites.__file__, 'cases': {}}
    for name, stack in STACKS.items():
        lam = laminated_plate(stack=stack, plyt=0.125e-3, laminaprop=CFRP)
        mat = lam.plies[0].matlamina
        lp = lam.calc_lamination_parameters()
        h = lam.h
        grad = GradABD()
        funcs = {
            'laminated_plate': (lambda: laminated_plate(
                stack=stack, plyt=0.125e-3, laminaprop=CFRP), 20),
            'calc_constitutive_matrix': (lam.calc_constitutive_matrix, 200),
            'calc_transverse_shear_stiffness': (
                lam.calc_transverse_shear_stiffness, 200),
            'calc_lamination_parameters': (lam.calc_lamination_parameters,
                                           2000),
            'laminate_from_LaminationParameters': (
                lambda: laminate_from_LaminationParameters(h, mat, lp), 2000),
            'calc_LP_grad': (lambda: grad.calc_LP_grad(h, mat, lp), 2000),
            'n_double_laminate': (lambda: n_double_laminate(
                h, 4, np.array([0., 45., -45., 90.]), mat), 2000),
        }
        case = {fname: timeit(f, number)*1e6
                for fname, (f, number) in funcs.items()}
        results['cases'][name] = case
        print(name)
        for fname, t in case.items():
            print('    %-36s %9.2f us' % (fname, t))
        sys.stdout.flush()
    if len(sys.argv) > 1:
        with open(sys.argv[1], 'w') as f:
            json.dump(results, f, indent=2)


if __name__ == '__main__':
    main()
