from numba import njit
import numpy as np


@njit(cache=True)
def fast_find_optimal_borders(
        candidates_count,
        L,
        min_seg_snps,
        min_seg_bp,
        snps_positions,
        candidate_numbers,
        first_snp_number,
        inf_score,
        total_snps_count,
        b_penalty
):
    score = np.zeros(candidates_count + 1, dtype=np.float64)
    best_boundaries_count = np.zeros(candidates_count + 1, dtype=np.int_)
    has_boundary_cache = np.zeros((candidates_count + 1, candidates_count), dtype=np.bool_)
    unique_positions = None
    if min_seg_snps:
        unique_positions = np.zeros(candidates_count + 1)
        current_index = 0
        for i in range(1, candidates_count + 1):
            if snps_positions[i] != snps_positions[i - 1]:
                current_index += 1
            unique_positions[i] = current_index
    for i in range(candidates_count + 1):
        score[i] = L[0, i]

        kf = -1
        current_optimal = score[i]

        check_optimal = True
        if min_seg_snps or min_seg_bp:
            last_index = candidate_numbers[i] + 1 - first_snp_number if i != candidates_count else -1
            first_index = 0
            if (min_seg_snps and unique_positions[last_index] - unique_positions[first_index] < min_seg_snps) or \
                    (min_seg_bp and snps_positions[last_index] - snps_positions[first_index] < min_seg_bp):
                score[i] -= inf_score
                check_optimal = False
        if check_optimal:
            for k in range(i):
                parameter_penalty = -1 / 2 * (best_boundaries_count[k] + 1) * (np.log(total_snps_count) + 1) * b_penalty

                z_penalty = 0
                if min_seg_snps or min_seg_bp:
                    last_index = candidate_numbers[i] + 1 - first_snp_number if i != candidates_count else -1
                    first_index = candidate_numbers[k] + 1 - first_snp_number
                    if (min_seg_snps and unique_positions[last_index] -
                        unique_positions[first_index] < min_seg_snps) or (
                            min_seg_bp and snps_positions[last_index] - snps_positions[first_index] < min_seg_bp):
                        z_penalty = -inf_score

                likelihood = score[k] + L[k + 1, i]
                candidate = likelihood + parameter_penalty + z_penalty
                if candidate > current_optimal:
                    current_optimal = candidate
                    score[i] = likelihood + z_penalty
                    kf = k
        if kf != -1:
            for j in range(0, kf):
                has_boundary_cache[i, j] = has_boundary_cache[kf, j]
            has_boundary_cache[i, kf] = 1
            for j in range(kf + 1, i):
                has_boundary_cache[i, j] = 0

        best_boundaries_count[i] = has_boundary_cache[i, :].sum()

    boundaries_indexes = [candidate_numbers[j] for j in range(candidates_count) if
                          has_boundary_cache[-1, j]]
    return boundaries_indexes