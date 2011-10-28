# -*- coding: utf-8 -*-


def relhum(t, w, p):
    """
    c NCL function: rh = relhum (t,w,p)

    Parameters
    ----------
    t : array_like
        Temperature [K]
    w : array_like
    mixing ratio [kg/kg]
    p : array_like
        pressure [Pa]




    """
    tp = t
    if tp > 375.16:
        tp = 375.16
    elif tp < 173.16:
        tp = 173.16

    it = tp - 173.16
    t2 = 173.16 + it
    es = (t2 + 1.0 - tp) * table[it + 1] + (tp - t2) * table[it + 2]
    es = es * 1.0e-01

    rh = (w * (p - 0.378 * es) / (0.622 * es)) * 100.0

    if rh < 0.:
        # FIXME: relhum = 0.0001
        pass

    return rh
