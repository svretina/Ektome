
def sphere(float x, float y, float z, float b):
    return (x + b) * (x + b) + y * y + z * z

cimport cython
@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing
def loop(float b, float exr2, double[:] x, double[:] y, double[:] z, double[:] xbounds, double[:] ybounds, double[:] zbounds , double[:, :, :] mask):
    cdef Py_ssize_t i, j, k
    for j in range(y.shape[0]):
        if y[j] < ybounds[0] or y[j] > ybounds[1]:
            continue
        for i in range(x.shape[0]):
            if x[i] < xbounds[0] or x[i] > xbounds[1]:
                continue
            for k in range(z.shape[0]):
                if z[k] < zbounds[0] or z[k] > zbounds[1]:
                    continue
                if sphere(x[i], y[j], z[k], b) >= exr2:
                    mask[i,j,k] = 1
    return mask
