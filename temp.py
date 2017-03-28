

def g2LC(g, Z0=50, fc=1e9, ladder_type='LC'): # only for LPFs
    if(ladder_type == 'LC'):
        g[1:-1:2] = g[1:-1:2] * Z0 / (2 * np.pi * fc)
        g[2:-1:2] = g[2:-1:2] / Z0 / (2 * np.pi * fc)
        g[0] = Z0
        g[-1] = g[-1] * Z0
        g2LC = g
    if(ladder_type == 'CL'):
        g[1:-1:2] = g[1:-1:2] / Z0 / (2 * np.pi * fc)
        g[2:-1:2] = g[2:-1:2] * Z0 / (2 * np.pi * fc)
        g[0] = Z0
        g[-1] = g[-1] * Z0
        g2LC = g
    return g2LC

def stepped_impedance(g2LC, Z0_high, Z0_low, lg_high, lg_low, fc, ladder_type='LC', optimize=False):
    LC_ang = np.zeros_like(g2LC)
    if(ladder_type == 'LC'):
        g2LC[1:-1:2] = (lg_high / 2 / np.pi) * np.arcsin(2 *
                                                         np.pi * fc * g2LC[1:-1:2] / Z0_high)
        LC_ang[1:-1:2] = g2LC[1:-1:2] * (2 * np.pi / lg_high)
        g2LC[2:-1:2] = (lg_low / 2 / np.pi) * np.arcsin(2 *
                                                        np.pi * fc * g2LC[2:-1:2] * Z0_low)
        LC_ang[2:-1:2] = g2LC[2:-1:2] * (2 * np.pi / lg_low)
        LC_len = g2LC
    if(ladder_type == 'CL'):
        pass
    if(optimize):
        pass # implement (5.4) from J. S. Hong

    return LC_len, LC_ang

def stub_LPF():
    # Implement Fig. 5.3 from J. S. Hong
    pass
