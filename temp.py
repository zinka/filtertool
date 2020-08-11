from sympy import nsolve, exp, Symbol
from sympy.core import symbol
import sympy as sym


def transformation_LP():
    return


def transformation_HP():
    return


def transformation_BP():
    return


def transformation_BS():
    return


def g2MQ_BPF(g, FBW):
    pass


def stepped_impedance(LC, Z0_high, Z0_low, lg_high, lg_low, fc, ladder_type="LC"):
    """
    Function to convert LC values into TX line physical lengths ... see Sec. 5.1.1, page 114
    @param LC: LC or CL values as an array
    @param Z0_high: Z0 of inductor line
    @param Z0_low: Z0 of capacitor line
    @param lg_high: guided wavelength of inductor line @ cut-off frequency
    @param lg_low: guided wavelength of capacitor line @ cut-off frequency
    @param fc: cut-off frequency
    @param ladder_type: "LC" or "CL"
    """
    LC_len_guess = zeros_like(LC)
    if ladder_type == "LC":
        LC_len_guess[1:-1:2] = (lg_high / 2 / pi) * arcsin(
            2 * pi * fc * LC[1:-1:2] / Z0_high
        )
        LC_len_guess[2:-1:2] = (lg_low / 2 / pi) * arcsin(
            2 * pi * fc * LC[2:-1:2] * Z0_low
        )
    if ladder_type == "CL":
        LC_len_guess[1:-1:2] = (lg_low / 2 / pi) * arcsin(
            2 * pi * fc * LC[1:-1:2] * Z0_low
        )
        print(LC_len_guess[1:-1:2])
        LC_len_guess[2:-1:2] = (lg_high / 2 / pi) * arcsin(
            2 * pi * fc * LC[2:-1:2] / Z0_high
        )

    #     lL = Symbol("lL")
    #     lC = Symbol("lC")
    #     f1 = (
    #         Z0_high * sym.sin(2 * pi * lL / lg_high)
    #         + Z0_low * sym.tan(pi * lC / lg_low)
    #         - 2 * pi * fc * LC[1]
    #     )
    #     f2 = (
    #         sym.sin(2 * pi * lC / lg_low) / Z0_low
    #         + 2 * sym.tan(pi * lL / lg_high) / Z0_high
    #         - 2 * pi * fc * LC[2]
    #     )
    #     print([LC_len_guess[1], LC_len_guess[2]])
    #     print(nsolve([f1, f2], [lL, lC], [LC_len_guess[1], LC_len_guess[2]]))

    return LC_len_guess

    # ==========================================================================
    # J. S. Hong, Ch.5, LPF Excercise
    # ==========================================================================

    #     g = Chebyshev_g(LAr=0.1, LAs=40, Ws=6)[1]
    #     LC = g2LC_LPF(g, Z0=50, fc=1e9, ladder_type="LC")
    #     LC_len = stepped_impedance(
    #         LC, Z0_high=93, Z0_low=24, lg_high=118e-3, lg_low=105e-3, fc=1e9
    #     )
    #     print(LC_len)

    g = Butterworth_g(LAs=20, Ws=1.5)
    print("order: ", g[0])
    print("parameters: ", g[1])
    g2S(g=g[1], Wstart=0.0001, Wstop=3, npoints=1001, dB_limit=-70)

    # Pozar's example ... ex 8.6, page 424
    g = Butterworth_g(LAs=20, Ws=1.5)[1]
    LC = g2LC_LPF(g, Z0=50, fc=2.5e9, ladder_type="CL")
    print("LC values:\n", LC)
    LC_len = stepped_impedance(
        LC,
        Z0_high=120,
        Z0_low=20,
        lg_high=118e-3,
        lg_low=105e-3,
        fc=1e9,
        ladder_type="CL",
    )


# def g2KJ_HP(g, LC_values, ZY0=1, ZYL=1, plot):
#
# def g2KJ_BP(g, LC_values, ZY0=1, ZYL=1, plot):
#
# def g2KJ_BS(g, LC_values, ZY0=1, ZYL=1, plot):


# def g2LC_LPF(g, Z0=50, fc=1e9, ladder_type="LC"):  # only for LPFs
#     if ladder_type == "LC":
#         g[1:-1:2] = g[1:-1:2] * Z0 / (2 * pi * fc)
#         g[2:-1:2] = g[2:-1:2] / Z0 / (2 * pi * fc)
#         g[0] = Z0
#         g[-1] = g[-1] * Z0
#         g2LC = g
#     if ladder_type == "CL":
#         g[1:-1:2] = g[1:-1:2] / Z0 / (2 * pi * fc)
#         g[2:-1:2] = g[2:-1:2] * Z0 / (2 * pi * fc)
#         g[0] = Z0
#         g[-1] = g[-1] * Z0
#         g2LC = g
#     return g2LC


def g2LC(g, Z0=50, fc=1e9, ladder_type="LC"):  # only for LPFs
    if ladder_type == "LC":
        g[1:-1:2] = g[1:-1:2] * Z0 / (2 * np.pi * fc)
        g[2:-1:2] = g[2:-1:2] / Z0 / (2 * np.pi * fc)
        g[0] = Z0
        g[-1] = g[-1] * Z0
        g2LC = g
    if ladder_type == "CL":
        g[1:-1:2] = g[1:-1:2] / Z0 / (2 * np.pi * fc)
        g[2:-1:2] = g[2:-1:2] * Z0 / (2 * np.pi * fc)
        g[0] = Z0
        g[-1] = g[-1] * Z0
        g2LC = g
    return g2LC


def stub_LPF():
    # Implement Fig. 5.3 from J. S. Hong
    pass


# ==========================================================================
# J. S. Hong, stub LPF
# ==========================================================================

# ==========================================================================
# J. S. Hong, elliptic LPF
# ==========================================================================

# ==========================================================================
# J. S. Hong, interdigital
# ==========================================================================

# ==========================================================================
# J. S. Hong, combline
# ==========================================================================

# ==========================================================================
# J. S. Hong, pseudo-combline
# ==========================================================================

# ==========================================================================
# J. S. Hong, stub bandpass
# ==========================================================================
