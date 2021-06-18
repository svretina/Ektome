#pythran export loop(float64, float64,
#                    float64 [], float64 [] , float64 [],
#                    float64 list, float64 list , float64 list,
#                    (int, int, int))
import numpy as np

def loop(b,exr,x,y,z,xbounds,ybounds,zbounds,mask_shape):
    mask = np.ones(mask_shape)*np.nan
    exr2 = exr * exr
    for j in range(y.shape[0]):
        if y[j] < ybounds[0] or y[j] > ybounds[1]:
            continue
        y2 = y[j] * y[j]
        for i in range(x.shape[0]):
            if x[i] < xbounds[0] or x[i] > xbounds[1]:
                continue
            x2 = (x[i] + b) * (x[i] + b)
            for k in range(z.shape[0]):
                z2 = z[k] * z[k]
                if z[k] < zbounds[0] or z[k] > zbounds[1]:
                    continue
                if x2+y2+z2 >= exr2:
                    mask[i,j,k] = 1
    return mask
