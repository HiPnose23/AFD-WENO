import numpy as np
import matplotlib.pylab as plt
from numpy.polynomial import polynomial as P
from scipy.special import bernoulli, factorial, gamma


def afd_coefficients(order):
    """
    Compute Alternative Finite Difference (AFD) coefficients as given by
    Merriman (2003) [1_].

    Parameters
    ----------
    order : int
        Order up to which to compute the AFD coefficients.

    Returns
    -------
    a : list
        List containg the AFD coefficients.

    Examples
    --------
    >>> import numpy as np
    >>> import afd_order_new as afd
    >>> afd.afd_coefficients(4)
    array([ 1.        ,  0.        , -0.04166667,  0.        ,  0.00121528])
    >>> np.allclose(afd.afd_coefficients(10),
    ...             [1.0, 0.0, -1.0/24.0, 0.0, 7.0/5760.0, 0.0,
    ...              -31.0/967680.0, 0.0, 127.0/154828800.0, 0.0,
    ...              -73.0/3503554560.0])
    True

    References
    ----------
    .. [1] Merriman, Understanding the Shu-Osher Conservative Finite Difference
           Form, 2003, Journal of Scientific Computing , Vol. 19, No. 1/3,
           p. 309-322.
           https://doi.org/10.1023/a:1025312210724
    """
    # AFD WENO coefficients from (adapted to scipy's bernoulli) Eq. (18) in
    # Merriman (2003) [1_]
    a = np.zeros(order+1)
    B = bernoulli(order+1)
    for n2 in range(0, order+1, 2):
        a[n2] = -(2.**n2 - 2.)*B[n2]/(2.**n2*factorial(n2, exact=True))
    return a


def afd_fd(r, f, A, V, xi, dxi, method=False):
    r"""
    Compute a conservative Alternative Finite Difference (AFD) approximation to
    the flux divergence as

    .. math::

       \frac{1}{V(\xi)} \frac{\mathrm{d}}{\mathrm{d}\xi} A(\xi) f(\xi)
       \approx   \frac{1}{V(\xi) \Delta \xi}
                 \left(  A(\xi + \Delta \xi/2) f(\xi + \Delta \xi/2)
                       - A(\xi - \Delta \xi/2) f(\xi - \Delta \xi/2)
                 \right)
               + \sum_{l=1}^{r}
                 \left(
                   \frac{\partial^{r} g}{\partial \xi^{r}}(\xi + \Delta \xi/2)
                 - \frac{\partial^{r} g}{\partial \xi^{r}}(\xi - \Delta \xi/2)
                 \right)

    where :math:`f(\xi)` is the flux function, :math:`A(\xi)` and :math:`V(\xi)`
    are the area and volume functions, and :math: `g(\xi) = A(\xi) f(\xi)`.

    Parameters
    ----------
    r : int
        Order of the AFD approximation.
    f : callable
        Flux function.
    A : callable
        Area function.
    V : callable
        Volume function.
    xi : float
        Coordinate :math:`\xi`.
    dxi : float
        Cell size :math:`\Delta \xi`.
    method : bool
        Method choice to compute flux divergence (default: False)

    Returns
    -------
    dF : float
        AFD approximation of flux divergence.
    """
    # verify that order is positive
    if r < 0:
        raise ValueError("Works only for r >= 0!")
    # set up stencil for left/right cell interface interpolation
    if r == 0:
        k = 0
    else:
        k = (r + 2)//2
    xi_ = xi + np.arange(-k, +k + 1)*dxi
    # compute AFD coefficients
    afd_coeffs = afd_coefficients(r)
    # left/right cell interface
    xim = xi - 0.5*dxi
    xip = xi + 0.5*dxi
    # interpolate g = A * f over left/right stencil
    if r > 0:
        gm = np.polyfit(xi_[ :-1], A(xi_[ :-1])*f(xi_[ :-1]), 2*k-1)
        gp = np.polyfit(xi_[1:  ], A(xi_[1:  ])*f(xi_[1:  ]), 2*k-1)
    # compute flux divergence
    if method: # method A: add up AFD sum of finite differences
        dF = (A(xip)*f(xip) - A(xim)*f(xim))/(V(xi)*dxi)
        for l in range(1, r+1):
            dF += ( afd_coeffs[l]
                   *(  np.polyval(np.polyder(gp, l), xip)
                     - np.polyval(np.polyder(gm, l), xim)
                    )/(V(xi)*dxi))*dxi**l
    else: # method B: add up AFD flux expansion and finite difference it
        Fm = A(xim)*f(xim)
        Fp = A(xip)*f(xip)
        for l in range(1, r+1):
            Fm += afd_coeffs[l]*np.polyval(np.polyder(gm, l), xim)*dxi**l
            Fp += afd_coeffs[l]*np.polyval(np.polyder(gp, l), xip)*dxi**l
        dF = (Fp - Fm)/(V(xi)*dxi)
    return dF


def EoC(hs, errs):
    """Estimate order of convergence given discretization parameter and
       errors."""
    # exclude errors once they stop decreasing or become zero
    k = np.where((errs[1:] >= errs[:-1]) | (errs[1:] == 0))[0]
    m = k[0] if k.size else errs.size
    if m < 2: # make sure that we have at least two measurements
        m = errs.size
    return float(np.polyfit(np.log(hs[:m]), np.log(errs[:m]), 1)[0])


def test01_cartesian():
    """Cartesian test case 1."""
    import sympy as sp
    # generate flux, area, volume and exact divergence functions using sympy
    xi = sp.symbols('xi')
    A_expr = 0*xi + 1
    V_expr = 0*xi + 1
    f_expr = sp.sin(xi)
    df_expr = sp.diff(A_expr*f_expr, xi)/V_expr
    A = sp.lambdify(xi, A_expr, modules='numpy')
    V = sp.lambdify(xi, V_expr, modules='numpy')
    f = sp.lambdify(xi, f_expr, modules='numpy')
    df_exact = sp.lambdify(xi, df_expr, modules='numpy')
    # evaluate AFD flux divergence and accuracy analysis
    xi = 0.2
    dxi = 0.5**np.arange(0, 15)
    for r in range(0, 10, 2):
        err = np.zeros_like(dxi)
        for i in range(len(dxi)):
            df = afd_fd(r, f, A, V, xi, dxi[i])
            err[i] = abs(df - df_exact(xi))
        eoc = EoC(dxi, err)
        # print(f"{err = }")
        # print(f"{eoc = }")
        plt.loglog(dxi, err, "-+", label=f"AFD {r = } (EoC = {eoc:4.2f})")
    plt.grid()
    plt.xlabel(r"$h$")
    plt.ylabel(r"Error")
    plt.legend()
    plt.title("Cartesian AFD flux divergence")
    plt.show()


def test01_curvilinear(m=0, xi0=0.2, keep_cell=0):
    """Curvilinear test case 1."""
    import sympy as sp
    # generate flux, area, volume and exact divergence functions using sympy
    xi = sp.symbols('xi', latex_name=r'\xi')
    A_expr = xi**m
    V_expr = xi**m
    f_expr = sp.sin(xi)
    df_expr = sp.diff(A_expr*f_expr, xi)/V_expr
    print(df_expr)
    A = sp.lambdify(xi, A_expr, modules='numpy')
    V = sp.lambdify(xi, V_expr, modules='numpy')
    f = sp.lambdify(xi, f_expr, modules='numpy')
    df_exact = sp.lambdify(xi, df_expr, modules='numpy')
    # evaluate AFD flux divergence and accuracy analysis
    xi = xi0
    if m > 0:
      dxi0 = 2.*xi
      if keep_cell > 0:
          dxi0 /= keep_cell
    else:
      dxi0 = 1.
    dxi = 0.5**np.arange(0, 25)
    fig, ax = plt.subplots()
    for r in range(0, 6, 2):
        err = np.zeros_like(dxi)
        for i in range(len(dxi)):
            if m > 0 and keep_cell > 0:
                xi = (keep_cell - 0.5)*dxi[i]
            df = afd_fd(r, f, A, V, xi, dxi[i])
            err[i] = abs(df - df_exact(xi))
        eoc = EoC(dxi, err)
        ax.loglog(dxi, err, "-+", label=f"AFD {r = } (EoC = {eoc:4.2f})")
    ax.grid()
    ax.set_xlabel(r"$h$")
    ax.set_ylabel(r"Error")
    ax.legend()
    return fig, ax


if __name__ == "__main__":
    # Cartesian case
    # test01_cartesian()
    fig, ax = test01_curvilinear(m=0)
    ax.set_title("AFD Cartesian")
    fig.savefig("AFD_Cartesian.pdf")

    # Cylindrical case
    fig, ax = test01_curvilinear(m=1)
    ax.set_title("AFD cylindrical")
    fig.savefig("AFD_cylindrical.pdf")
    fig, ax = test01_curvilinear(m=1, keep_cell=1)
    ax.set_title("AFD cylindrical for first cell")
    fig.savefig("AFD_cylindrical_firstcell.pdf")


    # Spherical case
    fig, ax = test01_curvilinear(m=2)
    ax.set_title("AFD spherical")
    fig.savefig("AFD_spherical.pdf")
    fig, ax = test01_curvilinear(m=2, keep_cell=1)
    ax.set_title("AFD spherical for first cell")
    fig.savefig("AFD_spherical_firstcell.pdf")
