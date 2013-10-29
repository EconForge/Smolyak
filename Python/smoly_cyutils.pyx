cpdef gen_points(long[:, ::1] indices, dict An):
    from itertools import product
    # cdef list points
    cdef long[::1] el
    cdef int i

    points = []
    for el in indices:
        points.extend(list(product(*[An[i] for i in el])))

        return points
