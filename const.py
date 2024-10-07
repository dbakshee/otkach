import math
import numpy as np

PI = math.pi
S32 = math.sqrt(3) / 2
G = math.sqrt(1 / S32)

class Sector:
    def __init__(self, radius):
        self.radius = radius
        self.angle = None
        self.len = None
        self.height = None

    def withAngle(self, angle):
        self.angle = angle
        self.len = 2 * self.radius * math.sin(angle / 2)
        self.height = 2 * self.radius * math.cos(angle / 2)
        return self

    def withLen(self, len):
        self.len = len
        self.angle = 2 * math.asin(len / 2 / self.radius)
        self.height = math.sqrt(self.radius**2 - (len / 2)**2)
        return self

    def withHeight(self, height):
        self.height = height
        self.angle = 2 * math.acos(height / self.radius)
        self.len = 2 * math.sqrt(self.radius**2 - height**2)
        return self

    def asector(self):
        return PI * self.radius**2 * (self.angle / (2 * PI))

    def atriangle(self):
        return 0.5 * self.len * self.height

    def asegment(self):
        return self.asector() - self.atriangle()


def A0(kf):
    return PI * kf**2


def A1(kf):
    if False:
        a = math.acos((G/2) / kf)  # half-angle of sector
        b = kf * math.sin(a)  # half-len of segment
        c = PI * kf**2 * (a / PI)  # sector area
        d = 0.5 * (2*b) * (G/2)  # triangle area
        e = c - d  # segment area
        return 2 * e
    else:
        a1 = Sector(kf).withHeight(G/2)
        return 2 * a1.asegment()


def T0(kf):
    if False:
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
    else:
        s = -G * S32 + math.sqrt(4 * kf**2 - G**2) / 2
        b = G + 2 * s * S32
        aside = Sector(kf).withLen(b).asegment()
        acenter = Sector(b).withLen(b).atriangle()
        return acenter + 3 * aside

def main():
    first = True
    for f in np.linspace(0.5, 2, 10):
        kf = G * f
        q = A0(kf) - 3*A1(kf) + 2*T0(kf)
        s = (A0(kf) - A1(kf)) / 2
        t = s - (A1(kf) - T0(kf))
        if first:
            first = False
            print(f"{'kf':5s}A0-3A1+2T0 (A0-A1)/2 A0/2-3A1/2+T0")
        print(f"{kf/G:.2f} {q:.4f} {s:.4f} {t:.4f}")


if __name__ == '__main__':
    main()
