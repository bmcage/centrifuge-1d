from __future__ import division, print_function

import numpy as np

# copy of scipy 0.13 as this deprecated starting scipy 0.14
# rewrite in terms of 0.14 when we can do as in MonoCubicInterp here.
from scipy.interpolate import KroghInterpolator

def _isscalar(x):
    """Check whether x is if a scalar type, or 0-dim"""
    return np.isscalar(x) or hasattr(x, 'shape') and x.shape == ()

class _Interpolator1D(object):
    """
    Common features in univariate interpolation

    Deal with input data type and interpolation axis rolling.  The
    actual interpolator can assume the y-data is of shape (n, r) where
    `n` is the number of x-points, and `r` the number of variables,
    and use self.dtype as the y-data type.

    Attributes
    ----------
    _y_axis
        Axis along which the interpolation goes in the original array
    _y_extra_shape
        Additional trailing shape of the input arrays, excluding
        the interpolation axis.
    dtype
        Dtype of the y-data arrays. Can be set via set_dtype, which
        forces it to be float or complex.

    Methods
    -------
    __call__
    _prepare_x
    _finish_y
    _reshape_yi
    _set_yi
    _set_dtype
    _evaluate

    """

    __slots__ = ('_y_axis', '_y_extra_shape', 'dtype')

    def __init__(self, xi=None, yi=None, axis=None):
        self._y_axis = axis
        self._y_extra_shape = None
        self.dtype = None
        if yi is not None:
            self._set_yi(yi, xi=xi, axis=axis)

    def __call__(self, x):
        """
        Evaluate the interpolant

        Parameters
        ----------
        x : array-like
            Points to evaluate the interpolant at.

        Returns
        -------
        y : array-like
            Interpolated values. Shape is determined by replacing
            the interpolation axis in the original array with the shape of x.

        """
        x, x_shape = self._prepare_x(x)
        y = self._evaluate(x)
        return self._finish_y(y, x_shape)

    def _evaluate(self, x):
        """
        Actually evaluate the value of the interpolator.
        """
        raise NotImplementedError()

    def _prepare_x(self, x):
        """Reshape input x array to 1-D"""
        x = np.asarray(x)
        x_shape = x.shape
        return x.ravel(), x_shape

    def _finish_y(self, y, x_shape):
        """Reshape interpolated y back to n-d array similar to initial y"""
        y = y.reshape(x_shape + self._y_extra_shape)
        if self._y_axis != 0 and x_shape != ():
            nx = len(x_shape)
            ny = len(self._y_extra_shape)
            s = (list(range(nx, nx + self._y_axis))
                 + list(range(nx)) + list(range(nx+self._y_axis, nx+ny)))
            y = y.transpose(s)
        return y

    def _reshape_yi(self, yi, check=False):
        yi = np.rollaxis(np.asarray(yi), self._y_axis)
        if check and yi.shape[1:] != self._y_extra_shape:
            ok_shape = "%r + (N,) + %r" % (self._y_extra_shape[-self._y_axis:],
                                           self._y_extra_shape[:-self._y_axis])
            raise ValueError("Data must be of shape %s" % ok_shape)
        return yi.reshape((yi.shape[0], -1))

    def _set_yi(self, yi, xi=None, axis=None):
        if axis is None:
            axis = self._y_axis
        if axis is None:
            raise ValueError("no interpolation axis specified")

        yi = np.asarray(yi)

        shape = yi.shape
        if shape == ():
            shape = (1,)
        if xi is not None and shape[axis] != len(xi):
            raise ValueError("x and y arrays must be equal in length along "
                             "interpolation axis.")

        self._y_axis = (axis % yi.ndim)
        self._y_extra_shape = yi.shape[:self._y_axis]+yi.shape[self._y_axis+1:]
        self.dtype = None
        self._set_dtype(yi.dtype)

    def _set_dtype(self, dtype, union=False):
        if np.issubdtype(dtype, np.complexfloating) \
               or np.issubdtype(self.dtype, np.complexfloating):
            self.dtype = np.complex_
        else:
            if not union or self.dtype != np.complex_:
                self.dtype = np.float_


class _Interpolator1DWithDerivatives(_Interpolator1D):
    def derivatives(self, x, der=None):
        """
        Evaluate many derivatives of the polynomial at the point x

        Produce an array of all derivative values at the point x.

        Parameters
        ----------
        x : array-like
            Point or points at which to evaluate the derivatives

        der : None or integer
            How many derivatives to extract; None for all potentially
            nonzero derivatives (that is a number equal to the number
            of points). This number includes the function value as 0th
            derivative.

        Returns
        -------
        d : ndarray
            Array with derivatives; d[j] contains the j-th derivative.
            Shape of d[j] is determined by replacing the interpolation
            axis in the original array with the shape of x.

        Examples
        --------
        >>> KroghInterpolator([0,0,0],[1,2,3]).derivatives(0)
        array([1.0,2.0,3.0])
        >>> KroghInterpolator([0,0,0],[1,2,3]).derivatives([0,0])
        array([[1.0,1.0],
               [2.0,2.0],
               [3.0,3.0]])

        """
        x, x_shape = self._prepare_x(x)
        y = self._evaluate_derivatives(x, der)

        y = y.reshape((y.shape[0],) + x_shape + self._y_extra_shape)
        if self._y_axis != 0 and x_shape != ():
            nx = len(x_shape)
            ny = len(self._y_extra_shape)
            s = ([0] + list(range(nx+1, nx + self._y_axis+1))
                 + list(range(1,nx+1)) +
                 list(range(nx+1+self._y_axis, nx+ny+1)))
            y = y.transpose(s)
        return y

    def derivative(self, x, der=1):
        """
        Evaluate one derivative of the polynomial at the point x

        Parameters
        ----------
        x : array-like
            Point or points at which to evaluate the derivatives

        der : integer, optional
            Which derivative to extract. This number includes the
            function value as 0th derivative.

        Returns
        -------
        d : ndarray
            Derivative interpolated at the x-points.  Shape of d is
            determined by replacing the interpolation axis in the
            original array with the shape of x.

        Notes
        -----
        This is computed by evaluating all derivatives up to the desired
        one (using self.derivatives()) and then discarding the rest.

        """
        x, x_shape = self._prepare_x(x)
        y = self._evaluate_derivatives(x, der+1)
        return self._finish_y(y[der], x_shape)

class PiecewisePolynomial(_Interpolator1DWithDerivatives):
    """Piecewise polynomial curve specified by points and derivatives

    This class represents a curve that is a piecewise polynomial. It
    passes through a list of points and has specified derivatives at
    each point. The degree of the polynomial may vary from segment to
    segment, as may the number of derivatives available. The degree
    should not exceed about thirty.

    Appending points to the end of the curve is efficient.

    Parameters
    ----------
    xi : array-like
        a sorted 1-d array of x-coordinates
    yi : array-like or list of array-likes
        yi[i][j] is the j-th derivative known at xi[i]   (for axis=0)
    orders : list of integers, or integer
        a list of polynomial orders, or a single universal order
    direction : {None, 1, -1}
        indicates whether the xi are increasing or decreasing
        +1 indicates increasing
        -1 indicates decreasing
        None indicates that it should be deduced from the first two xi
    axis : int, optional
        Axis in the yi array corresponding to the x-coordinate values.

    Notes
    -----
    If orders is None, or orders[i] is None, then the degree of the
    polynomial segment is exactly the degree required to match all i
    available derivatives at both endpoints. If orders[i] is not None,
    then some derivatives will be ignored. The code will try to use an
    equal number of derivatives from each end; if the total number of
    derivatives needed is odd, it will prefer the rightmost endpoint. If
    not enough derivatives are available, an exception is raised.

    """

    def __init__(self, xi, yi, orders=None, direction=None, axis=0):
        _Interpolator1DWithDerivatives.__init__(self, axis=axis)

        if axis != 0:
            try:
                yi = np.asarray(yi)
            except ValueError:
                raise ValueError("If yi is a list, then axis must be 0")

            preslice = ((slice(None,None,None),) * (axis % yi.ndim))
            slice0 = preslice + (0,)
            slice1 = preslice + (slice(1, None, None),)
        else:
            slice0 = 0
            slice1 = slice(1, None, None)

        yi0 = np.asarray(yi[slice0])
        self._set_yi(yi0)

        self.xi = [xi[0]]
        self.yi = [self._reshape_yi(yi0)]
        self.n = 1
        self.r = np.prod(self._y_extra_shape, dtype=np.int64)

        self.direction = direction
        self.orders = []
        self.polynomials = []
        self.extend(xi[1:],yi[slice1],orders)

    def _make_polynomial(self,x1,y1,x2,y2,order,direction):
        """Construct the interpolating polynomial object

        Deduces the number of derivatives to match at each end
        from order and the number of derivatives available. If
        possible it uses the same number of derivatives from
        each end; if the number is odd it tries to take the
        extra one from y2. In any case if not enough derivatives
        are available at one end or another it draws enough to
        make up the total from the other end.
        """
        n = order+1
        n1 = min(n//2,len(y1))
        n2 = min(n-n1,len(y2))
        n1 = min(n-n2,len(y1))
        if n1+n2 != n:
            raise ValueError("Point %g has %d derivatives, point %g has %d derivatives, but order %d requested" % (x1, len(y1), x2, len(y2), order))
        if not (n1 <= len(y1) and n2 <= len(y2)):
            raise ValueError("`order` input incompatible with length y1 or y2.")

        xi = np.zeros(n)
        yi = np.zeros((n, self.r), dtype=self.dtype)

        xi[:n1] = x1
        yi[:n1] = y1[:n1].reshape((n1, self.r))
        xi[n1:] = x2
        yi[n1:] = y2[:n2].reshape((n2, self.r))

        return KroghInterpolator(xi,yi,axis=0)

    def append(self, xi, yi, order=None):
        """
        Append a single point with derivatives to the PiecewisePolynomial

        Parameters
        ----------
        xi : float
            Input
        yi : array_like
            `yi` is the list of derivatives known at `xi`
        order : integer or None
            a polynomial order, or instructions to use the highest
            possible order

        """
        yi = self._reshape_yi(yi, check=True)
        self._set_dtype(yi.dtype, union=True)

        if self.direction is None:
            self.direction = np.sign(xi-self.xi[-1])
        elif (xi-self.xi[-1])*self.direction < 0:
            raise ValueError("x coordinates must be in the %d direction: %s" % (self.direction, self.xi))

        self.xi.append(xi)
        self.yi.append(yi)

        if order is None:
            n1 = len(self.yi[-2])
            n2 = len(self.yi[-1])
            n = n1+n2
            order = n-1

        self.orders.append(order)
        self.polynomials.append(self._make_polynomial(
            self.xi[-2], self.yi[-2],
            self.xi[-1], self.yi[-1],
            order, self.direction))
        self.n += 1

    def extend(self, xi, yi, orders=None):
        """
        Extend the PiecewisePolynomial by a list of points

        Parameters
        ----------
        xi : array_like
            A sorted list of x-coordinates.
        yi : list of lists of length N1
            ``yi[i]`` (if ``axis == 0``) is the list of derivatives known
            at ``xi[i]``.
        orders : int or list of ints
            A list of polynomial orders, or a single universal order.
        direction : {None, 1, -1}
            Indicates whether the `xi` are increasing or decreasing.

                +1 indicates increasing

                -1 indicates decreasing

            None indicates that it should be deduced from the first two `xi`.

        """
        if self._y_axis == 0:
            # allow yi to be a ragged list
            for i in range(len(xi)):
                if orders is None or _isscalar(orders):
                    self.append(xi[i],yi[i],orders)
                else:
                    self.append(xi[i],yi[i],orders[i])
        else:
            preslice = (slice(None,None,None),) * self._y_axis
            for i in range(len(xi)):
                if orders is None or _isscalar(orders):
                    self.append(xi[i],yi[preslice + (i,)],orders)
                else:
                    self.append(xi[i],yi[preslice + (i,)],orders[i])

    def _evaluate(self, x):
        if _isscalar(x):
            pos = np.clip(np.searchsorted(self.xi, x) - 1, 0, self.n-2)
            y = self.polynomials[pos](x)
        else:
            m = len(x)
            pos = np.clip(np.searchsorted(self.xi, x) - 1, 0, self.n-2)
            y = np.zeros((m, self.r), dtype=self.dtype)
            if y.size > 0:
                for i in range(self.n-1):
                    c = pos == i
                    y[c] = self.polynomials[i](x[c])
        return y

    def _evaluate_derivatives(self, x, der=None):
        if der is None and self.polynomials:
            der = self.polynomials[0].n
        if _isscalar(x):
            pos = np.clip(np.searchsorted(self.xi, x) - 1, 0, self.n-2)
            y = self.polynomials[pos].derivatives(x,der=der)
        else:
            m = len(x)
            pos = np.clip(np.searchsorted(self.xi, x) - 1, 0, self.n-2)
            y = np.zeros((der,m,self.r), dtype=self.dtype)
            if y.size > 0:
                for i in range(self.n-1):
                    c = pos == i
                    y[:,c] = self.polynomials[i].derivatives(x[c],der=der)
        return y

class PchipInterpolator(PiecewisePolynomial):
    """PCHIP 1-d monotonic cubic interpolation

    x and y are arrays of values used to approximate some function f,
    with ``y = f(x)``.  The interpolant uses monotonic cubic splines
    to find the value of new points.

    Parameters
    ----------
    x : ndarray
        A 1-D array of monotonically increasing real values.  `x` cannot
        include duplicate values (otherwise f is overspecified)
    y : ndarray
        A 1-D array of real values.  `y`'s length along the interpolation
        axis must be equal to the length of `x`. If N-D array, use axis
        parameter to select correct axis.
    axis : int, optional
        Axis in the y array corresponding to the x-coordinate values.

    Notes
    -----
    Assumes x is sorted in monotonic order (e.g. ``x[1] > x[0]``).

    """
    def __init__(self, x, y, axis=0):
        x = np.asarray(x)
        y = np.asarray(y)

        axis = axis % y.ndim

        xp = x.reshape((x.shape[0],) + (1,)*(y.ndim-1))
        yp = np.rollaxis(y, axis)
        self.origyi = yp

        data = np.empty((yp.shape[0], 2) + yp.shape[1:], y.dtype)
        data[:,0] = yp
        data[:,1] = PchipInterpolator._find_derivatives(xp, yp)

        s = list(range(2, y.ndim + 1))
        s.insert(axis, 1)
        s.insert(axis, 0)
        data = data.transpose(s)

        PiecewisePolynomial.__init__(self, x, data, orders=3, direction=None,
                                     axis=axis)

    @staticmethod
    def _edge_case(m0, d1, out):
        m0 = np.atleast_1d(m0)
        d1 = np.atleast_1d(d1)
        mask = (d1 != 0) & (m0 != 0)
        out[mask] = 1.0/(1.0/m0[mask]+1.0/d1[mask])

    @staticmethod
    def _find_derivatives(x, y):
        # Determine the derivatives at the points y_k, d_k, by using
        #  PCHIP algorithm is:
        # We choose the derivatives at the point x_k by
        # Let m_k be the slope of the kth segment (between k and k+1)
        # If m_k=0 or m_{k-1}=0 or sgn(m_k) != sgn(m_{k-1}) then d_k == 0
        # else use weighted harmonic mean:
        #   w_1 = 2h_k + h_{k-1}, w_2 = h_k + 2h_{k-1}
        #   1/d_k = 1/(w_1 + w_2)*(w_1 / m_k + w_2 / m_{k-1})
        #   where h_k is the spacing between x_k and x_{k+1}

        y_shape = y.shape
        if y.ndim == 1:
            # So that _edge_case doesn't end up assigning to scalars
            x = x[:,None]
            y = y[:,None]

        hk = x[1:] - x[:-1]
        mk = (y[1:] - y[:-1]) / hk
        smk = np.sign(mk)
        condition = ((smk[1:] != smk[:-1]) | (mk[1:] == 0) | (mk[:-1] == 0))

        w1 = 2*hk[1:] + hk[:-1]
        w2 = hk[1:] + 2*hk[:-1]
        # values where division by zero occurs will be excluded
        # by 'condition' afterwards
        with np.errstate(divide='ignore'):
            whmean = 1.0/(w1+w2)*(w1/mk[1:] + w2/mk[:-1])
        dk = np.zeros_like(y)
        dk[1:-1][condition] = 0.0
        dk[1:-1][~condition] = 1.0/whmean[~condition]

        # For end-points choose d_0 so that 1/d_0 = 1/m_0 + 1/d_1 unless
        #  one of d_1 or m_0 is 0, then choose d_0 = 0
        PchipInterpolator._edge_case(mk[0],dk[1], dk[0])
        PchipInterpolator._edge_case(mk[-1],dk[-2], dk[-1])

        return dk.reshape(y_shape)

    def root(self, y):
        """
        Evaluate for which `x` we have `f(x) = y`.
        As PchipInterpolator is monotonic, the solution is unique if the
        Interpolator has been constructed with y[i+1] > y[i]

        Parameters
        ----------
        y : array-like
            Point or points at which to evaluate `f^{-1}(y)=x`

        Returns
        -------
        d : ndarray
            root interpolated at the y-points.
        """
        # first determine the correct cubic polynomial
        y, y_shape = self._prepare_x(y)
        if _isscalar(y):
            pos = np.clip(np.searchsorted(self.origyi, y) - 1, 0, self.n-2)
            poly = self.polynomials[pos]
            x = self._poly_inv(poly, y)
        else:
            m = len(y)
            pos = np.clip(np.searchsorted(self.origyi, y) - 1, 0, self.n-2)
            x = np.zeros((m, self.r), dtype=self.dtype)
            if x.size > 0:
                for i in range(self.n-1):
                    c = pos == i
                    if not any(c):
                        continue
                    poly  = self.polynomials[i]
                    ress = self._poly_inv(poly, y[c])
                    x[c] = ress[:, np.newaxis]
        return self._finish_y(x, y_shape)

    @staticmethod
    def _poly_inv(cubic, y):
        """Given a cubic KroghInterpolator polynomial,
           we determine the root f(x) = y, where x must be
           bounded by the edges of the Krogh polynomial

        Parameters
        ----------
        y : array-like
            Point or points at which to evaluate `f^{-1}(y)=x`
        Returns
        -------
        d : ndarray
            root of the cubic polynomial interpolated at the y-points.
        """
        from scipy.interpolate._ppoly import real_roots
        x0 = cubic.xi[0]
        x2 = cubic.xi[2]
        if (cubic.yi[0][0] >= cubic.yi[2][0]):
            raise ValueError("Not a strictly increasing monotone function")
        #convert Krogh c structure to c structure of PPoly
        ourc = np.empty((4,1), cubic.c.dtype)
        ourc[0, 0] = cubic.c[3,0]
        ourc[1, 0] = cubic.c[2,0]
        ourc[2, 0] = cubic.c[1,0]
        ourc[1,0] += cubic.c[3]*(x0-x2)
        ourc = ourc.reshape(4,1,1)
        y = np.asarray(y)
        result = np.empty(y.shape, float)
        for ind, yval in enumerate(y):
            ourc[3, 0, 0] = cubic.c[0,0] - yval
            roots = real_roots(ourc, np.array([x0,x2], float), 0, 0)
            result[ind] = roots[0]
        return result


class MonoCubicInterp(PchipInterpolator):
    """
    extend pchip with fast inversion
    """
    def root(self, y):
        """
        Evaluate for which `x` we have `f(x) = y`.
        As PchipInterpolator is monotonic, the solution is unique if the
        Interpolator has been constructed with y[i+1] > y[i]

        Parameters
        ----------
        y : array-like
            Point or points at which to evaluate `f^{-1}(y)=x`

        Returns
        -------
        d : ndarray
            root interpolated at the y-points.
        """
        # first determine the correct cubic polynomial
        y, y_shape = self._prepare_x(y)
        if _isscalar(y):
            pos = np.clip(np.searchsorted(self.origyi, y) - 1, 0, self.n-2)
            poly = self.polynomials[pos]
            x = self._poly_inv(poly, y)
        else:
            m = len(y)
            pos = np.clip(np.searchsorted(self.origyi, y) - 1, 0, self.n-2)
            x = np.zeros((m, self.r), dtype=self.dtype)
            if x.size > 0:
                for i in range(self.n-1):
                    c = pos == i
                    if not any(c):
                        continue
                    poly = self.polynomials[i]
                    ress = self._poly_inv(poly, y[c])
                    x[c] = ress[:, np.newaxis]
        return self._finish_y(x, y_shape)

    @staticmethod
    def _poly_inv(cubic, y):
        """Given a cubic KroghInterpolator polynomial,
           we determine the root f(x) = y, where x must be
           bounded by the edges of the Krogh polynomial

        Parameters
        ----------
        y : array-like
            Point or points at which to evaluate `f^{-1}(y)=x`
        Returns
        -------
        d : ndarray
            root of the cubic polynomial interpolated at the y-points.
        """
        from scipy.interpolate._ppoly import real_roots
        x0 = cubic.xi[0]
        x2 = cubic.xi[2]
        if (cubic.yi[0][0] >= cubic.yi[2][0]):
            raise ValueError("Not a strictly increasing monotone function")
        #convert Krogh c structure to c structure of PPoly
        ourc = np.empty((4,1), cubic.c.dtype)
        ourc[0, 0] = cubic.c[3,0]
        ourc[1, 0] = cubic.c[2,0]
        ourc[2, 0] = cubic.c[1,0]
        ourc[1,0] += cubic.c[3]*(x0-x2)
        ourc = ourc.reshape(4,1,1)
        y = np.asarray(y)
        result = np.empty(y.shape, float)
        for ind, yval in enumerate(y):
            ourc[3, 0, 0] = cubic.c[0,0] - yval
            roots = real_roots(ourc, np.array([x0,x2], float), 0, 0)
            if len(roots[0]) == 1:
                result[ind] = roots[0]
            elif len(roots[0]) == 0:
                #print ('roots', roots[0], x0, x2, ourc, yval)
                #check to be sure not the extreme points
                #print ('cubxy', cubic.xi, cubic.yi)
                if np.allclose(cubic.yi[0], yval):
                    result[ind] = cubic.xi[0]
                elif np.allclose(cubic.yi[2], yval):
                    result[ind] = cubic.xi[2]
                else:
                    print ('WARNING: no root found in interval')
                    raise Exception ('No Root found')
            else:
                raise Exception ('Multiple roots found')
        return result

class QuadraticBspline(object):
    """
    Implementation of interpolating quadratic B-Spline
    The grid goes from 0 to n=len(input)-1
    The splines go from -1 to n, with B_{-1} to left of interval, and B_n to right
    Intervals goes from 0 to n-1, with only internal intervals
    """
    def __init__(self, xi, yi):
        self.xi = xi
        self.origyi = yi
        #for interpolating origyi, we need at knots
        #   y0 = origy0 and yi = 2*origyi-yi-1
        self.yi = np.empty(len(yi), float)
        self.yi[0] = yi[0]
        for ind in range(1,len(yi)):
            self.yi[ind] = 2*yi[ind]-self.yi[ind-1]
        self.size = len(self.xi) -1;
        self.h = xi[1:] - xi[:-1]
        #we create the spline functions
        self.splines = []

    def quadspline(self, j, x):
        """
        spline j evaluated at scalar x, spline j has peak from j to j+1 point
        """
        if j> self.size or j < -1:
            raise Exception("A Quad requested that does not exist " + str(j))
        elif j == -1:
            if self.xi[0] <= x <= self.xi[1]:
                return (self.xi[1]-x)**2 / (2*self.h[0]**2)
            else:
                return 0.
        elif j == self.size:
            if self.xi[-2] <= x <= self.xi[-1]:
                return (x-self.xi[-2])**2 / (2*self.h[-1]**2)
            else:
                return 0.
        if j==0:
            lbound = self.xi[0]
        else:
            lbound = self.xi[j-1]
        if j == self.size -1:
            rbound = self.xi[-1]
        else:
            rbound = self.xi[j+2]
        if not (lbound <= x <= rbound):
            return 0.
        if x > self.xi[j+1]:
            return (self.xi[j+2]-x)**2 / (2*self.h[j+1]**2)
        elif x > self.xi[j]:
            y = (0.5 + (x-self.xi[j])/self.h[j] \
                    - (x-self.xi[j])**2/self.h[j]**2)
            #if j == self.size-1:
            #    y += (x-self.xi[j])**2 / (2*self.h[j]**2)
            #    print (x,y)
            #elif j == 0:
            #    y += (self.xi[1]-x)**2 / (2*self.h[0]**2)
            return y
        else:
            return (x-self.xi[j-1])**2 / (2*self.h[j-1]**2)

    def quadspline_der1(self, j, x):
        """
        spline j evaluated at scalar x
        """
        if j> self.size or j < -1:
            raise Exception("A Quad requested that does not exist " + str(j))
        elif j == -1:
            if self.xi[0] <= x <= self.xi[1]:
                return -2*(self.xi[1]-x) / (2*self.h[0]**2)
            else:
                return 0.
        elif j == self.size:
            if self.xi[-2] <= x <= self.xi[-1]:
                return 2*(x-self.xi[-2]) / (2*self.h[-1]**2)
            else:
                return 0.
        if j==0:
            lbound = self.xi[0]
        else:
            lbound = self.xi[j-1]
        if j == self.size -1:
            rbound = self.xi[-1]
        else:
            rbound = self.xi[j+2]
        if not (lbound <= x <= rbound):
            return 0.
        if x > self.xi[j+1]:
            return -2*(self.xi[j+2]-x) / (2*self.h[j+1]**2)
        elif x > self.xi[j]:
            return (1/self.h[j]- 2*(x-self.xi[j])/self.h[j]**2)
        else:
            return 2*(x-self.xi[j-1]) / (2*self.h[j-1]**2)

    def quadspline_der2(self, j, x):
        """
        spline j evaluated at scalar x
        """
        if j> self.size or j < -1:
            raise Exception("A Quad requested that does not exist " + str(j))
        elif j == -1:
            if self.xi[0] <= x <= self.xi[1]:
                return 2 / (2*self.h[0]**2)
            else:
                return 0.
        elif j == self.size:
            if self.xi[-2] <= x <= self.xi[-1]:
                return 2 / (2*self.h[-1]**2)
            else:
                return 0.
        if j==0:
            lbound = self.xi[0]
        else:
            lbound = self.xi[j-1]
        if j == self.size -1:
            rbound = self.xi[-1]
        else:
            rbound = self.xi[j+2]
        if not (lbound <= x <= rbound):
            return 0.
        if x > self.xi[j+1]:
            return 2 / (2*self.h[j+1]**2)
        elif x > self.xi[j]:
            return (-2/self.h[j]**2)
        else:
            return 2 / (2*self.h[j-1]**2)

    def __minplusval(self, ind):
        if ind-2 < 0:
            yim = 0
        else:
            yim = self.yi[ind-2]
        if ind >= self.size+1:
            yip = 0
        else:
            yip = self.yi[ind]
        return yim, yip

    def __call__(self, x):
        """
        Evaluate the interpolant

        Parameters
        ----------
        x : array-like
            Points to evaluate the interpolant at.

        Returns
        -------
        y : array-like
            Interpolated values.
        """
        if _isscalar(x):
            ind = np.searchsorted(self.xi, x)
            if (ind == 0 or ind == self.size+1):
                #out of bounds
                y = 0.
            else:
                yim, yip = self.__minplusval(ind)

                y = (yim*self.quadspline(ind-2,x) +
                     self.yi[ind-1]*self.quadspline(ind-1,x) +
                     yip*self.quadspline(ind,x)  )
        else:
           y = np.zeros(len(x), dtype=x.dtype)
           inds = np.searchsorted(self.xi, x)
           for i, ind in enumerate(inds):
                if (ind == 0 or ind == self.size+1):
                    #out of bounds
                    y[i] = 0.
                else:
                   yim, yip = self.__minplusval(ind)
                   y[i] = (yim*self.quadspline(ind-2,x[i]) +
                           self.yi[ind-1]*self.quadspline(ind-1,x[i]) +
                           yip*self.quadspline(ind,x[i])  )
        return y

    def derivative(self, x, der=1):
        """
        Evaluate one derivative of the polynomial at the point x

        Parameters
        ----------
        x : array-like
            Point or points at which to evaluate the derivatives

        der : integer, optional
            Which derivative to extract. This number includes the
            function value as 0th derivative.

        Returns
        -------
        d : ndarray
            Derivative interpolated at the x-points.
        """
        if der == 0:
            return self.__call__(x)
        elif der == 1:
            if _isscalar(x):
                ind = np.searchsorted(self.xi, x)
                if (ind == 0 or ind == self.size+1):
                    #out of bounds
                    y = 0.
                else:
                    yim, yip = self.__minplusval(ind)
                    y = (yim*self.quadspline_der1(ind-2,x) +
                         self.yi[ind  ]*self.quadspline_der1(ind-1 ,x) +
                         yip*self.quadspline_der1(ind,x)  )
            else:
               y = np.zeros(len(x), dtype=x.dtype)
               inds = np.searchsorted(self.xi, x)
               for i, ind in enumerate(inds):
                   if (ind == 0 or ind == self.size+1):
                       #out of bounds
                       y[i] = 0.
                   else:
                       yim, yip = self.__minplusval(ind)
                       y[i] = (yim*self.quadspline_der1(ind-2,x[i]) +
                               self.yi[ind  ]*self.quadspline_der1(ind-1,x[i]) +
                               yip*self.quadspline_der1(ind,x[i])  )
            return y
        elif der == 2:
            if _isscalar(x):
                ind = np.searchsorted(self.xi, x)
                if (ind == 0 or ind == self.size+1):
                    #out of bounds
                    y = 0.
                else:
                    yim, yip = self.__minplusval(ind)
                    y = (yim*self.quadspline_der2(ind-2,x) +
                         self.yi[ind  ]*self.quadspline_der2(ind-1,x) +
                         yip*self.quadspline_der2(ind,x)  )
            else:
               y = np.zeros(len(x), dtype=x.dtype)
               inds = np.searchsorted(self.xi, x)
               for i, ind in enumerate(inds):
                   if (ind == 0 or ind == self.size+1):
                       #out of bounds
                       y[i] = 0.
                   else:
                       yim, yip = self.__minplusval(ind)
                       y[i] = (yim*self.quadspline_der2(ind-2,x[i]) +
                               self.yi[ind  ]*self.quadspline_der2(ind-1,x[i]) +
                               yip*self.quadspline_der2(ind,x[i])  )
            return y
        else:
            if _isscalar(x):
                y = 0
            else:
               y = np.zeros(len(x), dtype=x.dtype)
            return y

class PiecewiseLinear(object):
    """
    Implementation of piecewise linear function with smoothing in the
    neighbourhood of the grid point to stay differentiable.
    """
    def __init__(self, xi, yi, dx_eps=1e-1):
        self.dx_epsilon = dx_eps
        self.xi = xi
        self.yi = yi
        self.dxi = xi[1:] - xi[:-1]
        self.Dyi = np.asarray((yi[1:] - yi[:-1])/self.dxi, dtype=float)

        if not all(self.dxi > 2*dx_eps):
            raise Exception("For piecewise linear interpolation the points "
                            "must be at least at distance 2*dx_eps. Current "
                            "value of dx_eps: {}". format(dx_eps))
        size = np.alen(xi)
        self.size = size

    def _polynomial(self, x, x1, y1, D0, D1):
        """
        Find value of x on polynomial joing two piecewise linear segments,
        where (x1, y1) is the junction of the two linear segments and
        D0 is the derivative of the first and D1 of the second segment.
        """
        t = (x-x1)/self.dx_epsilon
        return y1 + self.dx_epsilon*((1 + t*t)*(D1 - D0) + 2*t*(D1 + D0))/4

    def _polynomial_der1(self, x, x1, y1, D0, D1):
        """
        Find derivative of x on polynomial joing two piecewise linear segments,
        where (x1, y1) is the junction of the two linear segments and
        D0 is the derivative of the first and D1 of the second segment.
        """
        t = (x-x1)/self.dx_epsilon
        return 0.5*((D1 - D0)*t + D1 + D0)

    def _polynomial_der2(self, x, x1, y1, D0, D1):
        """
        Find second derivative of polynomial joing two piecewise linear
        segments, where (x1, y1) is the junction of the two linear segments and
        D0 is the derivative of the first and D1 of the second segment.
        """
        return 0.5*(D1 - D0)/self.dx_epsilon

    def __call__(self, x):
        xi  = self.xi
        yi  = self.yi
        Dyi = self.Dyi
        dx_eps = self.dx_epsilon
        size = self.size
        end = self.size - 1

        y = np.zeros(np.alen(x))
        if np.isscalar(x):
            x = (x, )

        idxs = np.searchsorted(self.xi, x)
        for (i, idx) in enumerate(idxs):
            if idx == 0: # out of bounds
                y[i] = yi[0] + Dyi[0]*(x[i]-xi[0])
            elif idx == size:
                y[i] = yi[-1] + (x[i]-xi[-1])*Dyi[-1]
            elif (idx < size - 1) and (xi[idx] - x[i] < dx_eps):
                # close to the right end of segment - all but last segment
                y[i] = self._polynomial(x[i], xi[idx], yi[idx],
                                        Dyi[idx-1], Dyi[idx])
            elif (idx > 1) and (x[i] - xi[idx-1] < dx_eps):
                # close to the left end of segment - all but first segment
                y[i] = self._polynomial(x[i], xi[idx-1], yi[idx-1],
                                        Dyi[idx-2], Dyi[idx-1])
            else:
                y[i] = yi[idx-1] + Dyi[idx-1]*(x[i]-xi[idx-1])
        return y

    def derivative(self, x, der=1):
        """
        NOTE: Second derivative in the points (x-eps), (x+eps) is not continuous
        function.
        """
        if np.isscalar(x):
            x = (x, )

        if der == 0:
            return self.__call__(x)

        Dy = np.zeros(np.alen(x))
        size = self.size
        xi  = self.xi
        yi  = self.yi
        Dyi = self.Dyi
        dx_eps = self.dx_epsilon

        idxs = np.searchsorted(self.xi, x)

        if der == 1:
            for (i, idx) in enumerate(idxs):
                if idx == 0: # out of bounds
                    Dy[i] = Dyi[0]
                elif idx == size:
                    Dy[i] = Dyi[-1]
                elif (idx < size -1 ) and (xi[idx] - x[i] < dx_eps):
                    # close to the right end of segment - all but last segment
                    Dy[i] = self._polynomial_der1(x[i], xi[idx], yi[idx],
                                                  Dyi[idx-1], Dyi[idx])
                elif (idx > 1) and (x[i] - xi[idx-1] < dx_eps):
                    # close to the left end of segment - all but first segment
                    Dy[i] = self._polynomial_der1(x[i], xi[idx-1], yi[idx-1],
                                                  Dyi[idx-2], Dyi[idx-1])
                else:
                    Dy[i] = Dyi[idx-1]

        elif der == 2:
            for (i, idx) in enumerate(idxs):
                if (idx == 0) or (idx == size): # out of bounds
                    Dy[i] = 0.0
                elif (idx < size - 1) and (xi[idx] - x[i] < dx_eps):
                    # close to the right end of segment - all but last segment
                    Dy[i] = self._polynomial_der2(x[i], xi[idx], yi[idx],
                                                  yi[idx-1], Dyi[idx])
                elif (idx > 1) and (x[i] - xi[idx-1] < dx_eps):
                    # close to the left end of segment - all but first segment
                    x[i]
                    xi[idx-1]
                    yi[idx-1]
                    Dyi[idx-2]
                    Dyi[idx-1]
                    #print(idx, size - 1)

                    Dy[i] = self._polynomial_der2(x[i], xi[idx-1], yi[idx-1],
                                                  Dyi[idx-2], Dyi[idx-1])
                else:
                    Dy[i] = 0.0
        # else all derivatives are 0

        return Dy

if __name__ == "__main__":
    #test of QuadraticBspline
    x = np.linspace(0, 10, 21)
    mcub = MonoCubicInterp
    mcub_lin = mcub(x, x)
    lin = x
    mcub_kwa = mcub(x, np.power(x,2))
    kwa = x**2
    mcub_cub = mcub(x, np.power(x,3))
    cub = x**3
    mcub_limit = mcub(np.array([0,2]), np.array([0,4]))
    h = -np.power(10, x[::-1])/5
    u = 1/np.power(1+ np.power(-0.015 * h, 1.3), 1. - 1./1.3)
    mcub_vGn = mcub(h, u)

    bspl = QuadraticBspline
    bspl_lin = bspl(x, x)
    bspl_kwa = bspl(x, np.power(x,2))
    bspl_cub = bspl(x, np.power(x,3))
    bspl_limit = bspl(np.array([0,2]), np.array([0,4]))
    h = -np.power(10, x[::-1])/5
    u = 1/np.power(1+ np.power(-0.015 * h, 1.3), 1. - 1./1.3)
    bspl_vGn = bspl(h, u)

    import pylab
    pylab.figure(1)
    xaxis = np.linspace(-1000., 1e-1, 100)
    mcuby = mcub_vGn(xaxis)
    bsply = bspl_vGn(xaxis)
    pylab.xlim(xmin=-2001)
    pylab.plot(h, 1/np.power(1+ np.power(-0.015 * h, 1.3), 1. - 1./1.3), 'g-')
    pylab.plot(xaxis, mcuby, 'bo')
    pylab.plot(xaxis, bsply, 'ro')

    pylab.figure(2)
    xaxis = np.linspace(0, 10, 100)
    mcuby = mcub_lin(xaxis)
    bsply = bspl_lin(xaxis)
    pylab.plot(x, lin, 'go')
    pylab.plot(xaxis, mcuby, 'bo')
    pylab.plot(xaxis, bsply, 'ro')
    pylab.plot(xaxis, map(bspl_lin.quadspline, np.ones(len(xaxis),float),xaxis), 'ko')

    # test of PiecewiseLinear
    pylab.figure(3)
    plin = PiecewiseLinear
    plin_vGn = plin(h, u, dx_eps=0.2)
    xaxis = np.linspace(-1000., 1e-1, 10000)
    pylab.xlim(xmin=-2001)
    y=plin_vGn(xaxis)
    pylab.plot(xaxis, y,'-')
    pylab.plot(xaxis.T, plin_vGn.derivative(xaxis, der=1).T, 'o')
    #pylab.plot(xaxis, plin_vGn.derivative(xaxis, der=2), 'o')
    #pylab.plot(xaxis, plin_vGn.derivative(xaxis, der=3))
    pylab.plot(h, 1/np.power(1+ np.power(-0.015 * h, 1.3), 1. - 1./1.3), 'g-')
    pylab.show()
