import numpy as np


def cheby_eval(x, n):
    past_val = np.ones_like(x)
    curr_val = x

    if n == 1:
        return past_val

    if n == 2:
        return curr_val

    for i in xrange(3, n + 1):
        temp = 2*x*curr_val - past_val
        past_val = curr_val
        curr_val = temp

    return curr_val


def cheby2n(x, n):
    # computes the chebychev polynomials of the first kind
    dim = x.shape
    results = np.zeros((n+1, ) + dim)
    results[0, ...] = np.ones(dim)
    results[1, ...] = x
    for i in range(2,n+1):
        results[i,...] = 2 * x * results[i-1,...] - results[i-2,...]
    return results


def permute(a):
    """
    Creates all unique combinations of a list a that is passed in.
    Function is based off of a function written by John Lettman:
    TCHS Computer Information Systems.  My thanks to him.
    """

    a.sort() # Sort.

    ## Output the first input sorted.
    yield list(a)

    i = 0
    first = 0
    alen = len(a)

    ## "alen" could also be used for the reference to the last element.

    while(True):
        i = alen - 1

        while(True):
            i -= 1 # i--

            if(a[i] < a[(i + 1)]):
                j = alen - 1

                while(not (a[i] < a[j])): j -= 1 # j--

                a[i], a[j] = a[j], a[i] # swap(a[j], a[i])
                t = a[(i + 1):alen]
                t.reverse()
                a[(i + 1):alen] = t

                # Output current.
                yield list(a)

                break # next.

            if(i == first):
                a.reverse()

                # yield list(a)
                return


def num_grid_points(d, mu):
    if mu == 1:
        return 2*d + 1

    if mu == 2:
        return 1 + 4*d + 4*d*(d-1)/2.

    if mu == 3:
        return 1 + 8*d + 12*d*(d-1)/2. + 8*d*(d-1)*(d-2)/6.


def m_i(i):
    if i < 0:
        raise ValueError('i must be positive')
    elif i == 0:
        return 0
    elif i == 1:
        return 1
    else:
        return 2**(i - 1) + 1
