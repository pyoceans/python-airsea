#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# relhum.py
#
# purpose:
# author:   Filipe P. A. Fernandes
# e-mail:   ocefpaf@gmail
# web:      http://ocefpaf.tiddlyspot.com/
# created:  22-Apr-2011
# modified: Fri 22 Apr 2011 09:25:57 AM EDT
#
# obs:
#

def relhum(t,w,p):
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

    it = tp-173.16
    t2 = 173.16+it
    es = ( t2 + 1.0 - tp) * table[it+1] + (tp-t2) * table[it+2]
    es = es*1.0e-01

    rh = ( w *(p-0.378*es) / (0.622*es) )* 100.0

    if rh < 0.:
        relhum = 0.0001

    return rh
