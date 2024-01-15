import math
import numpy as np

PI = math.pi
S32 = math.sqrt(3) / 2
G = math.sqrt(1 / S32)


def A0(kf):
    return PI * kf**2


def A1(kf):
    a = math.acos((G/2) / kf)  # half-angle of sector
    b = kf * math.sin(a)  # half-len of segment
    c = PI * kf**2 * (a / PI)  # sector area
    d = 0.5 * (2*b) * (G/2)  # triangle area
    e = c - d  # segment area
    return 2 * e


def T0(kf):
    s = -G * S32 + math.sqrt(4 * kf**2 - G**2) / 2
    sx = S32 * s
    sy = s / 2
    a = math.atan2((G/2) + sx, G*S32 + sy)  # half-angle of sector
    b = kf * math.sin(a)  # half-len of segment
    c = PI * kf**2 * (a / PI)  # sector area
    h = kf * math.cos(a)  # sub-triangle height
    d = 0.5 * (2*b) * h  # sub-triangle area
    e = c - d  # segment area
    f = 0.5 * (2*b) * (2*b*S32)  # triangle area
    return f + 3*e


def main():
    for f in np.linspace(0.5, 2, 10):
        kf = G * f
        q = A0(kf) - 3*A1(kf) + 2*T0(kf)
        s = (A0(kf) - A1(kf)) / 2
        t = s - (A1(kf) - T0(kf))
        print(f"{kf/G:.2f} {q} {s} {t}")


if __name__ == '__main__':
    main()
