#! /usr/bin/env python

# Author: Srinivasa Rao Zinka (srinivas . zinka [at] gmail . com)
# Copyright (c) 2020 Srinivasa Rao Zinka
# License: New BSD License

from matplotlib import pyplot as plt
from numpy import (
    pi,
    ceil,
    log,
    arange,
    sin,
    arccosh,
    sqrt,
    tanh,
    sinh,
    ones,
    r_,
    arctan,
    arctan2,
    array,
)
from warnings import warn
from pathspec import __updated__
from coupling import cutoff
import skrf as rf

__updated__ = "11-08-2020"


def Butterworth_g(LAs, Ws):
    """
    Function to obtain low pass parameters of Butterworth filters (Sec. 3.2.1, page 40)
    @param LAs: minimum stop-band attenuation in dB at Ws
    @param Ws: stop-band cutoff frequency
    """
    n = ceil(log(10 ** (0.1 * LAs) - 1) / (2 * log(Ws)))
    i = arange(n + 2)
    g = 2 * sin((2 * i - 1) * pi / 2 / n)
    g[0] = 1
    g[-1] = 1
    return n, g


def Chebyshev_g(LAr, LAs, Ws):
    """
    Function to obtain low pass parameters of Chebyshev filters (Sec. 3.2.2, page 40)
    @param LAr: passband ripple in dB
    @param LAs: minimum stop-band attenuation in dB at Ws
    @param Ws: stop-band cutoff frequency
    """
    n = int(
        ceil(
            arccosh(sqrt((10 ** (0.1 * LAs) - 1) / (10 ** (0.1 * LAr) - 1)))
            / arccosh(Ws)
        )
    )
    beta = log(1 / tanh(LAr / 17.37))
    gamma = sinh(beta / 2 / n)
    g = ones(n + 2)
    g[1] = (2 / gamma) * sin(pi / 2 / n)
    g[-1] = 1 / tanh(beta / 4) ** 2 if n % 2 == 0 else 1

    for i in range(2, n + 1):
        num = 4 * sin((2 * i - 1) * pi / 2 / n) * sin((2 * i - 3) * pi / 2 / n)
        den = gamma ** 2 + sin((i - 1) * pi / n) ** 2
        g[i] = num / den / g[i - 1]
    return n, g


def Elliptic_g():
    return


def Gaussian_g():
    return


def g2RLC_LPF(g, Z0=50, fc=1e9, ladder_type="LC"):
    """
    Function to transform g values to LPF LC values
    @param g: low-pass prototype's g values
    @param Z0: characteristic impedance of port1 and 2
    @param fc: cut-off frequency in Hz
    @param ladder_type: "LC" or "CL"
    """
    if ladder_type == "LC":
        g[1:-1:2] = g[1:-1:2] * Z0 / g[0] / (2 * pi * fc)
        g[2:-1:2] = g[2:-1:2] * g[0] / Z0 / (2 * pi * fc)
    if ladder_type == "CL":
        g[2:-1:2] = Z0 / g[0] / (2 * pi * fc) * g[2:-1:2]
        g[1:-1:2] = g[0] / Z0 / (2 * pi * fc) * g[1:-1:2]
    g[0] = g[0] * Z0
    g[-1] = g[-1] * Z0
    return g


def RLC2S_LPF(
    g, Wstart, Wstop, npoints=1001, dB_limit=-100, ladder_type="LC", W=True, plot=False
):
    """
    Function to plot S21 and S11 from g parameters.
    Main purpose of this function is to ensure that the g values are indeed correct.
    @param g: g parameters as a list [g_0, g_1, ..., g_{n+1}]
    @param start: starting frequency (avoid 0)
    @param stop: stop frequency
    @param npoints: number of points for plotting purpose
    """
    if not ladder_type == "LC":
        warn("ladder_type must be LC for the time being")
        return

    freq = rf.Frequency(
        start=Wstart / (2 * pi), stop=Wstop / (2 * pi), unit="Hz", npoints=npoints
    )
    tl_media = rf.DefinedGammaZ0(freq, z0=1, gamma=1j * freq.w / rf.c)

    # GND and port elements
    gnd = rf.Circuit.Ground(freq, name="gnd")
    port1 = rf.Circuit.Port(freq, name="port1", z0=g[0])
    port2 = rf.Circuit.Port(freq, name="port2", z0=g[-1])

    g = g[1:-1]

    # L and C elements for LC ladder
    ckt_elements = []
    ind_flag = True  # L starts first
    for i in range(len(g)):
        ckt_elements.append(
            tl_media.inductor(g[i], name=str(i))
        ) if ind_flag else ckt_elements.append(tl_media.capacitor(g[i], name=str(i)))
        ind_flag = not ind_flag  # toggle between inductors and cpacitors

    # node connections
    cnx = []
    for i in range(int(len(g) / 2) + 1):
        if i == 0:
            cnx.append([(port1, 0), (ckt_elements[i], 0)])  # starting
        elif i == int(len(g) / 2):
            if int(len(g)) % 2 == 0:  # ending ... even order case
                cnx.append(
                    [
                        (ckt_elements[2 * i - 2], 1),
                        (ckt_elements[2 * i - 1], 0),
                        (port2, 0),
                    ]
                )
            else:  # ending ... odd order case
                cnx.append(
                    [
                        (ckt_elements[2 * i - 2], 1),
                        (ckt_elements[2 * i - 1], 0),
                        (ckt_elements[2 * i], 0),
                    ]
                )
                cnx.append(
                    [(ckt_elements[2 * i], 1), (port2, 0), ]
                )

        else:  # in between
            cnx.append(
                [
                    (ckt_elements[2 * i - 2], 1),
                    (ckt_elements[2 * i - 1], 0),
                    (ckt_elements[2 * i], 0),
                ]
            )
    # GND connections ... same for even and odd cases
    gnd_cnx = []
    for i in range(int(len(g) / 2) + 1):
        if i == 0:
            gnd_cnx.append((gnd, 0))
        elif i == int(len(g) / 2):
            gnd_cnx.append((ckt_elements[2 * i - 1], 1))
        else:
            gnd_cnx.append((ckt_elements[2 * i - 1], 1))
    cnx.append(gnd_cnx)

    # create a circuit out of all the elements and connections
    cir = rf.Circuit(cnx)
    ntw = cir.network
    freq = 2 * pi * ntw.frequency.f if W else ntw.frequency.f
    if plot:
        title = (
            "Transfer Functions vs Angular Frequency"
            if W
            else "Transfer Functions vs Frequency"
        )
        xlabel = "w" if W else "f"
        ylabel = "|H(w)|" if W else "|H(f)|"
        rf.plotting.plot_rectangular(
            freq,
            cutoff(rf.mathFunctions.complex2dB(ntw.s[:, 1, 0]), dB_limit),
            x_label=xlabel,
            y_label=ylabel,
            title=title,
            show_legend=True,
            label="S21",
        )
        rf.plotting.plot_rectangular(
            freq,
            cutoff(rf.mathFunctions.complex2dB(ntw.s[:, 0, 0]), dB_limit),
            label="S11",
        )

        if not (port1.z0[0] == port2.z0[0]):
            warn("port2 and port2 impedances are not equal")

        s_ax = plt.gca()
        s_ax.axvspan(0, 1, facecolor="r", alpha=0.2)
        s_ax.axhspan(0, -3, facecolor="g", alpha=0.2)
        #     s_ax.axvspan(0, 1.6, facecolor="r", alpha=0.2)

        plt.show()
    return freq, ntw.s[:, 0, 0], ntw.s[:, 1, 0]


def g2KJ_LPF(g, LC_values, ZY0=1, ZYL=1):
    """
    Lowpass prototype with K or J inverters ... Sec. 3.4.2, page 57
    L or C values are in general obtained from lowpass/highpass/bandpass/bandstop characteristics
    X or B values need to be converted to L or C values (xi = W_0*Li and bi = W_0*C_i) ... Fig.3.20, page 59
    @param g: lowpass prototype's g values
    @param LC_values: L or C values (see Fig. 3.18 a and b)
    @param ZY0: Z0 or Y0
    @param ZYL: ZL or YL
    """
    num1 = r_[ZY0, LC_values]
    num2 = r_[LC_values, ZYL]
    den1 = g[0:-1]
    den2 = r_[g[1:-1], g[-1]]
    KJ = sqrt((num1 * num2) / (den1 * den2))
    return KJ


def endcoupled_BPF_tht_Cg(g, FBW, Y0, W0):
    """
    Function to calculate (initial) electrical lengths and gap coupling values
    @param g: lowpass prototype's g values
    @param FBW: fractional bandwidth
    @param Y0: characteristic admittance of port1 and 2
    @param W0: BPF's center frequency (rad/s)
    """
    J = g2KJ_LPF(g, LC_values=pi / 2 * FBW * ones((g.size - 2,)))  # normalised wrt y0
    B = J / (1 - J * J)
    Cg = B / W0 * Y0
    tht = []
    for i in range(0, len(B) - 1):
        tht.append(pi - 0.5 * (arctan(2 * B[i]) + arctan(2 * B[i + 1])))
    return tht, Cg


def endcoupled_BPF2(tht, lg0, Y0, W0, Cp):
    """
    Function to provide the final corrected electrical lengths
    @param tht: initial electrical lengths of resonators
    @param lg0: guided wavelength
    @param Y0: characteristic admittance of port1 and 2
    @param W0: BPF's center frequency (rad/s)
    @param Cp: Shunt capacitances calculated from EM solver
    """
    final_len = array(tht) * lg0 / 2 / pi
    for i in range(len(tht)):
        tmp1 = W0 * lg0 / Y0 / 2 / pi * Cp[i]
        tmp2 = W0 * lg0 / Y0 / 2 / pi * Cp[i + 1]
        final_len[i] = final_len[i] - tmp1 - tmp2
    return final_len


def parallel_coupled_BPF(g, FBW, Y0):
    """
    Function to calculate Z0e and Z0o for parallel-coupled BPD.
    In addition to calculating widtha dn gaps from Z0 values, one should compensate 
    for microstrip open ends.
    @param g: lowpass prototype's g values
    @param FBW: fractional bandwidth
    @param Y0: characteristic admittance of port1 and 2
    """
    J = g2KJ_LPF(g, LC_values=pi / 2 * FBW * ones((g.size - 2,)))  # normalised wrt y0
    Z0e = (1 / Y0) * (1 + J + J * J)
    Z0o = (1 / Y0) * (1 - J + J * J)
    return Z0e, Z0o


def hairpin_BPF(g, FBW):
    """
    Function to calculate coupling coefficients and external quality factors 
    corresponding to hair-pin BPFs 
    @param g: lowpass prototype's g values
    @param FBW: fractional bandwidth
    """
    den1 = g[1:-1]
    den2 = r_[g[2:-1], g[-1]]
    M = (FBW / sqrt(den1 * den2))[0:-1]  # discarding the last value
    Qi = g[0] * g[1] / FBW
    Qo = g[-2] * g[-1] / FBW
    return Qi, M, Qo

# ==============================================================================


if __name__ == "__main__":

    print("======================== Simulation Starts ========================")

    # ==========================================================================
    # Butterworth and Chebyshev g values
    # ==========================================================================

    #         print(Butterworth_g(LAs=40, Ws=2))  # Sec. 3.2.1, page 40
    #         print(Chebyshev_g(LAr=0.1, LAs=40, Ws=2))  # Sec. 3.2.2, page 42

    # ===========================================================================
    # Plotting H(jw) from g values
    # ===========================================================================

    #         g = Butterworth_g(LAs=50, Ws=2)
    #         print("order: ", g[0])
    #         print("parameters: ", g[1])

    #         g = Chebyshev_g(LAr=0.1, LAs=70, Ws=2)
    #         print("order: ", g[0])
    #         print("parameters: ", g[1])

    #         g2S(g=g[1], Wstart=0.0001, Wstop=3, npoints=1001, dB_limit=-70)

    # ===========================================================================
    # g to RLC values and plot w.r.t f
    # ===========================================================================

    #     g = Chebyshev_g(LAr=0.1, LAs=40, Ws=6)[1]
    #     RLC2S(g, Wstart=0.0001, Wstop=3, npoints=1001, dB_limit=-70, plot=True)
    #
    #     RLC = g2RLC_LPF(g, Z0=50, fc=1e9, ladder_type="LC")
    #     RLC2S(
    #         RLC,
    #         Wstart=0.0001,
    #         Wstop=2 * pi * 3e9,
    #         npoints=1001,
    #         dB_limit=-70,
    #         W=False,
    #         plot=True,
    #     )

    # ==========================================================================
    # J. S. Hong, Ch.4, g to K conversion, (based on example given in Sec. 5.2.3)
    # ==========================================================================

    #     n, g = Chebyshev_g(LAr=0.1, LAs=40, Ws=3)
    #     print(n, g)
    #     KJ = g2KJ(g, LC_values=ones((g.size - 2,)))  # -2 because of g_0 and g_{n+1}
    #     print(KJ)

    # ==========================================================================
    # J. S. Hong, end coupled, Sec. 5.2.1, page 125
    # ==========================================================================

    #     g = Chebyshev_g(LAr=0.1, LAs=20, Ws=3)[1]
    #     print(g)
    #     RLC2S_LPF(g, Wstart=0.0001, Wstop=3, npoints=1001, dB_limit=-70, plot=False)
    #     [tht, Cg] = endcoupled_BPF_tht_Cg(g, 0.028, Y0=1 / 50, W0=2 * pi * 6e9)
    #     print(tht)
    #     print(Cg)
    #     Cp = [0.0049e-12, 0.0457e-12, 0.0457e-12, 0.0049e-12]
    #     endcoupled_BPF2(tht, lg0=18.27e-3, Y0=1 / 50, W0=2 * pi * 6e9, Cp=Cp)

    # ===========================================================================
    # J. S. Hong, parallel-coupled, Sec. 5.2.2, page 129
    # ===========================================================================

    #     g = Chebyshev_g(LAr=0.1, LAs=40, Ws=3)[1]
    #     print(g)
    #     RLC2S_LPF(g, Wstart=0.0001, Wstop=3, npoints=1001, dB_limit=-70, plot=False)
    #     parallel_coupled_BPF(g, FBW=0.15, Y0=1 / 50)

    # ==========================================================================
    # J. S. Hong, hair-pin filter, Sec. 5.2.3, page 131
    # ==========================================================================

    g = Chebyshev_g(LAr=0.1, LAs=40, Ws=3)[1]
    RLC2S_LPF(g, Wstart=0.0001, Wstop=3, npoints=1001, dB_limit=-70, plot=False)
    hairpin_BPF(g, FBW=0.2)

# ==========================================================================
# Notes Section
# ==========================================================================
# TODO obtain filter phase
# TODO stepped, LC ladder, semi lumped LPF
# TODO end and edge coupled filters
# composite filters from Pozar
# Richard transformation and kuroda identies pozar

# refactor everythin and commit

# TODO implement elliptic and Gaussian low-pass prototype filters
# TODO Ws calculations an actual L an C values
# update filter design article to include these traditional filters as well
# TODO what to do with Chebyshev even case where S21 is greater than 0dB?
