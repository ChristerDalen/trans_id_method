import numpy as np
from scipy.optimize import fsolve

# Shared global variables
d = None
xi = None
tau = None

def Chen(Y, R, T, Kc):
    """
    CHEN 
    """
    global d, xi, tau

    cp1, cp2, A, cm1, tm1, tp1 = id_6_param3(Y, R, T)

    cinf = (cp1 * cp2 - cm1**2) / (cp1 + cp2 - 2 * cm1)

    Km = cinf / (Kc * (A - cinf))
    K = cinf / A
    H = (1/3) * (
        (cp1 - cinf) / cinf +
        (cinf - cm1) / (cp1 - cinf) +
        (cp2 - cinf) / (cinf - cm1)
    )

    xi = -np.log(H) / np.sqrt(np.pi**2 + (np.log(H))**2)
    tau = (tm1 - tp1) * np.sqrt(1 - xi**2) / np.pi
    d = 2 * tp1 - tm1

    wu = fsolve(fun, 1.0)[0]

    Gcl = K / np.sqrt((1 - tau**2 * wu**2)**2 + (2 * xi * tau * wu)**2)
    GcGp = Gcl / np.sqrt(1 + 2 * Gcl + Gcl**2)
    GM = 1 / GcGp
    Kcu = Kc * GM

    tm = (1 / wu) * np.sqrt(Kcu**2 * Km**2 - 1)
    dm = (1 / wu) * (np.pi - np.arctan(tm * wu))

    return np.real(Km), np.real(tm), np.real(dm)


def FSL(Y, R, T, kc, fc):
    yp1, yp2, dR, ym1, tm1, tp1, tp2 = id_6_param3(Y, R, T)
    yinf = (yp1 * yp2 - ym1**2) / (yp1 + yp2 - 2 * ym1)

    tr = 1.694
    y0 = 0
    OS = (yp1 - yinf) / (yinf - y0)

    wd = 2 * np.pi / (tp2 - tp1)
    K = (yinf - y0) / dR

    phi = (tp2 - tr) * wd
    xi = np.cos(phi)
    wn = wd / np.sin(phi)

    a = 2 * xi * wn
    b = wn**2

    kp = K / (kc * (1 - K))
    tau = (
        -((1 + kc * kp) * a) +
        np.sqrt(((-(1 + kc * kp) * a))**2 - 4 * b * (1 - (kc**2) * 
(kp**2)))
    ) / (2 * b)

    A = -(1 + kc * kp) / tau
    to = -2 * A / b

    return kp, tau, to


def fun(wu):
    global d, xi, tau
    return -d * wu - np.arctan((2 * xi * tau * wu) / np.sqrt(1 - tau**2 * 
wu**2)) + np.pi


def id_6_param3(Y, R, T):
    dyp = np.max(Y)
    ip = np.argmax(Y)
    tp = T[ip]

    dys = np.max(R)

    Y_tail = Y[ip:]
    iu_rel = np.argmin(Y_tail)
    dyu = np.min(Y_tail)
    tu = T[ip + iu_rel]

    Y_tail2 = Y[ip + iu_rel - 1:]
    ip2_rel = np.argmax(Y_tail2)
    dyp2 = np.max(Y_tail2)
    tp2 = T[ip2_rel + ip + iu_rel - 1]

    return dyp, dyp2, dys, dyu, tu, tp, tp2


def JR(Y, R, T, Kc):
    cp1, cp2, A, cm1, tm1, tp1 = id_6_param3(Y, R, T)
    cinf = (cp1 * cp2 - cm1**2) / (cp1 + cp2 - 2 * cm1)

    v = (cinf - cm1) / (cp1 - cinf)
    xi = -np.log(v) / np.sqrt(np.pi**2 + (np.log(v))**2)

    Km = cinf / (Kc * (A - cinf))
    K = Km * Kc
    tau = (tm1 - tp1) * np.sqrt(1 - xi**2) / np.pi

    gamma1, gamma2, delta = -0.6143, 0.1247, 0.3866
    alpha = 2 * xi * tau * (1 + K) / (delta + gamma1 * K)
    beta = -1 / (delta + gamma1 * K)

    A1 = beta**2 * gamma2 * K + beta * delta
    B1 = 2 * gamma2 * K * alpha * beta + alpha * delta
    C1 = gamma2 * K * alpha**2 - tau**2 * (1 + K)

    tm = (-B1 + np.sqrt(B1**2 - 4 * A1 * C1)) / (2 * A1)
    dm = alpha + beta * tm

    return Km, tm, dm


def Lee(Y, R, T, Kc):
    cp1, cp2, A, cm1, tm1, tp1 = id_6_param3(Y, R, T)
    cinf = (cp1 * cp2 - cm1**2) / (cp1 + cp2 - 2 * cm1)

    v = (cinf - cm1) / (cp1 - cinf)
    xi = -np.log(v) / np.sqrt(np.pi**2 + (np.log(v))**2)
    tau = (tm1 - tp1) * np.sqrt(1 - xi**2) / np.pi

    Km = cinf / (Kc * (A - cinf))
    alpha = xi / tau
    beta = np.sqrt(1 - xi**2) / tau
    v_angle = np.arctan(beta / alpha)

    dm = (v_angle + np.pi / 4) / beta
    for _ in range(10):
        dm = (1 / beta) * (
            v_angle + np.arctan(
                beta * np.exp(-alpha * dm) /
                (Kc * Km * np.sqrt(alpha**2 + beta**2) * np.cos(beta * dm 
- v_angle))
            )
        )

    tm = (1 / alpha) * (1 + Kc * Km * np.exp(alpha * dm) * np.cos(beta * 
dm))

    return Km, tm, dm


def MF(Y, R, T, Kc, Ti, fc):
    global d, xi, tau

    cp1, cp2, A, cm1, tm1, tp1 = id_6_param3(Y, R, T)
    css = (cp1 * cp2 - cm1**2) / (cp1 + cp2 - 2 * cm1)

    K = css / A
    p = -(1 / (2 * np.pi)) * np.log((cp2 - css) / (cp1 - css))
    xi = np.sqrt(p**2 / (1 + p**2))
    tau = (tp2 - tp1) * np.sqrt(1 - xi**2) / (2 * np.pi)

    Sc = (css * np.ones(len(Y)) - Y).sum() * 0.001
    d = Sc / css - 2 * xi * tau

    wc = fsolve(fun, 0.0)[0]
    Kp = Ti / (Kc * Sc) * css

    M = K / np.sqrt((1 - tau**2 * wc**2)**2 + (2 * tau * xi * wc)**2)
    tp = np.sqrt((Kc * Kp)**2 * (1 + Ti**2 * wc**2) - M**2 * Ti**2 * 
wc**2) / (M * wc**2 * Ti)
    dp = (1 / wc) * (np.arctan(wc * Ti) + np.arctan(1 / (tau * wc)))

    return K, tp, dp


def YS_est(Y, R, T, Kc):
    cp1, cp2, A, cm1, tm1, tp1 = id_6_param3(Y, R, T)
    cinf = (cp1 * cp2 - cm1**2) / (cp1 + cp2 - 2 * cm1)

    v = (cinf - cm1) / (cp1 - cinf)
    xi = -np.log(v) / np.sqrt(np.pi**2 + (np.log(v))**2)

    Km = cinf / (Kc * (A - cinf))
    K = Km * Kc

    s1 = xi * np.sqrt(K + 1) + np.sqrt(xi**2 * (K + 1) + K)
    s2 = np.sqrt((1 - xi**2) * (K + 1))

    dm = 2 * (tm1 - tp1) * s2 / (np.pi * s1)
    tm = (tm1 - tp1) * s1 * s2 / np.pi

    return Km, tm, dm


def YS(Y, R, T, Kc):
    cp1, cp2, A, cm1, tm1, tp1 = id_6_param3(Y, R, T)
    cinf = (cp1 * cp2 - cm1**2) / (cp1 + cp2 - 2 * cm1)

    v = (cinf - cm1) / (cp1 - cinf)
    xi = -np.log(v) / np.sqrt(np.pi**2 + (np.log(v))**2)

    Km = cinf / (Kc * (A - cinf))
    K = Km * Kc

    s1 = xi * np.sqrt(K + 1) + np.sqrt(xi**2 * (K + 1) + K - 1)
    s1 = np.real(s1)
    s2 = np.sqrt((1 - xi**2) * (K + 1))

    dm = 2 * (tm1 - tp1) * s2 / (np.pi * s1)
    tm = (tm1 - tp1) * s1 * s2 / np.pi

    return Km, tm, dm

