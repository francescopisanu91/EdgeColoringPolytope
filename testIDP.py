import pulp
import numpy as np
from tqdm import tqdm


def find_integer_decomposition(target, generators, max_coeff=1000, verbose=False):
    target = np.asarray(target, dtype=int)
    generators = [np.asarray(v, dtype=int) for v in generators]

    k = len(generators)
    m = len(target)

    prob = pulp.LpProblem("MonoidMembership", pulp.LpStatusOptimal)

    lambdas = [
        pulp.LpVariable(f"lambda_{i}", lowBound=0, upBound=max_coeff, cat="Integer")
        for i in range(k)
    ]

    # vincoli
    for j in range(m):
        prob += (
            pulp.lpSum(lambdas[i] * generators[i][j] for i in range(k))
            == target[j]
        )
    prob += pulp.lpSum(lambdas) == max_coeff


    # OBIETTIVO: soluzione più piccola
    prob += pulp.lpSum(lambdas)

    status = prob.solve(pulp.PULP_CBC_CMD(msg=verbose))

    if pulp.LpStatus[status] in ("Optimal", "Feasible"):
        coeffs = []
        for l in lambdas:
            val = pulp.value(l)
            if val is None:
                val = 0
            coeffs.append(int(round(val)))
        return True, coeffs

    return False, None

def hasIDP(targets, generators, max_coeff=1000, verbose=False):
    n = len(targets)

    for target in tqdm(
        targets,
        desc="Testing IDP targets",
        total=n,
        unit="target"
    ):
        state, _ = find_integer_decomposition(
            target,
            generators,
            max_coeff=max_coeff,
            verbose=verbose
        )

        if state is False:
            return False

    return True
