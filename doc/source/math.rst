Mathematical Background
=======================

This page is still a work in progress.  I will try to provide a very basic description of each of the pieces that are necessary.

Smolyak Grid
------------

One of the standard forms of building an approximation grid is to use a simple tensor-product.  While this seems relatively trivial for low dimensions, but the number of points to evaluate quickly becomes intractable for larger dimensions.  In Bellman (1961) this problem is referred to as the "curse of dimensionality."  Two years after Bellman's paper, Sergey Smolyak introduced a numerical technique where the number of grid points needed to approximate grew polynomially instead of exponentially.  As stated in Judd, Maliar, Maliar, Valero; the idea behind this technique was "that some elements produced by tensor-product rules are more important for representing multidimensional functions than the others."  The tensor-product typically takes as parameters the dimension of the grid and the number of points, :math:`n` to be evaluated at in each dimension which produces a grid with :math:`n^d` points (Note: There are other ways of doing this where the number of points is different in each dimension, but the resulting number of points is similar).  The Smolyak grid takes an "accuracy" parameter :math:`\mu` and the number of dimensions as parameters.  The number of points in the grid is determined by the dimension and :math:`\mu`.  The number of Smolyak grid points at :math:`\mu=1` is :math:`1 + 2d`, at :math:`\mu=2` it is :math:`1+4d+4d(d-1)`, etc...  Notice that the number of grid points grows linearly at :math:`\mu=1`, and quadratically at :math:`\mu=2`.

The standard construction of a Smolyak grid uses nested sets of points.  One typically uses the extrema of the Chebyshev Polynomials, which are known as the Chebyshev-Gauss-Lobatto points.  We will continue our description using these points since the code is implemented using them.  The nested sets require two conditions.  First, that each set :math:`S_i` has :math:`2^{i-1}` points for :math:`i \geq 2` and 1 if :math:`i=1`.  Secondly, that the sets are nested.  The first four nested sets are:

.. math::

    i=1 \text{ } : \text{ }S_1 = \{ 0 \}

    i=2 \text{ } : \text{ }S_2 = \{ -1, 0, 1 \}

    i=3 \text{ } : \text{ }S_3 = \{ -1, -\frac{1}{\sqrt{2}}, 0, \frac{1}{\sqrt{2}}, 1 \}

    i=4 \text{ } : \text{ }S_4 = \{ -1, -\frac{\sqrt{2 + \sqrt{2}}}{2}, -\frac{1}{\sqrt{2}}, -\frac{\sqrt{2 - \sqrt{2}}}{2}, 0, \frac{\sqrt{2 - \sqrt{2}}}{2}, \frac{1}{\sqrt{2}}, \frac{\sqrt{2 + \sqrt{2}}}{2}, 1 \}


One then takes the tensor product of the unidimensional sets and then picks out the products that satisfy :math:`d \leq \sum_{j=1}^d i_j \leq d + \mu` where :math:`i_j` is the index of the unidimensional sets.  For example if the parameters were :math:`d = 2, \mu = 2` then you would build the first four nested sets of points (as shown above) and then take all of the tensor products that satisfied :math:`2 \leq i_1 + i_2 \leq 4` which would give you :math:` S_1 \text{X} S_1, S_1 \text{X} S_2, S_2 \text{X} S_1, S_2 \text{X} S_2`.  Then these would be your points.  We will define the set of grid points to be :math:`\mathcal{H}^{d, \mu} := \bigcup_{d \leq |i^*| \leq d + \mu}  \prod{i=1}^{d+\mu} S_i` where :math:`i^* = \begin{bmatrix} i_1 & \dots i_d \end{bmatrix}`.

You can see how there would be repeated points and hence this method could be improved upon.  This is one of the key results of the Judd, Maliar, Maliar, Valero (2013) paper.  Instead of building nested sets, they build disjoint sets :math:`A_i` such that :math:`A_1 = \{ 0 \}` and :math:`A_i = S_i \backslash S_{i-1}` for all :math:`i \geq 2`.  Then the points are taken from the tensor-products of these sets in the same fashion as described above.

The following is a computationally efficient way of finding these grid points.  The first step is to create the first :math:`\mu + 1` unidimensional disjoint sets :math:`A_i` as described above. Then one should create a vector of possible values i.e. :math:`\begin{bmatrix} 1, 2, \dots, \mu + 1 \end{bmatrix}`.  It is important to note that these possible values only range from 1 to :math:`\mu + 1` since the smallest possible index is 1 (To see this think of given :math:`d` indexes.  The smallest the first :math:`d-1` of them can be is 1 which would sum to :math:`d-1` hence the last index could only be valued up to :math:`\mu + 1`.  We can then check all of the combinations (with replacement) of these to find the sets of indexes that would work.  Once we have these we can permute them to capture all of the indexes that would work.  Then you only take the tensor-products of these sets (Add reference to our code smol_inds and build_grid here).


Smolyak Basis Polynomials
-------------------------
