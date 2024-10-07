settings.outformat = "pdf";
unitsize(2cm);
real g = 1; // lattice
real kf = g * 1.2; // k fermi
real s32 = sqrt(3) / 2;
real h = g * s32;
real s = -g * s32 + sqrt(4 * kf^2 - g^2) / 2;
real sx = s32 * s;
real sy = s / 2;
real an = atan2(g/2 + sx, h + sy) / pi * 180; // half-angle top arc
real bn = atan2(sy, g + sx) / pi * 180; // half-angle right arc
pen tr = opacity(0.05); // transparent

path twoT0minusA1 =
    arc((g/2, -h), kf, angle1=90+an, angle2=90-an) --
    arc((0,0), kf, angle1=bn, angle2=-bn) --
    arc((g/2,  h), kf, angle1=-90+an, angle2=-90-an) --
    arc((g,0), kf, angle1=180+bn, angle2=180-bn) -- cycle;

// draw

dot((0,0));
dot((g,0));
dot((g/2,h));
dot((g/2,-h));

filldraw(circle((0,0), kf), green+tr, green+white + dashed);
filldraw(circle((g/2,h), kf), green+tr, green+white + dashed);
filldraw(circle((g,0), kf), green+tr, green+white + dashed);
filldraw(circle((g/2,-h), kf), green+tr, green+white + dashed);

fill(twoT0minusA1, green+opacity(0.2));

draw(arc((g/2,-h), r=kf, angle1=90+an, angle2=90-an/2), red, arrow=Arrow(HookHead));
draw(arc((g/2,-h), r=kf, angle1=90-an/2, angle2=90-an), red);

draw(arc((0,0), r=kf, angle1=bn, angle2=-bn/2), red, arrow=Arrow(HookHead));
draw(arc((0,0), r=kf, angle1=-bn/2, angle2=-bn), red);

draw(arc((g/2,h), r=kf, angle1=-90+an, angle2=-90-an), blue+1.5);
draw(arc((g/2,h), r=kf, angle1=-90+an, angle2=-90-an/2), red, arrow=Arrow(HookHead));
draw(arc((g/2,h), r=kf, angle1=-90-an/2, angle2=-90-an), red);

draw(arc((g,0), r=kf, angle1=180+bn, angle2=180-bn/2), red, arrow=Arrow(HookHead));
draw(arc((g,0), r=kf, angle1=180-bn/2, angle2=180-bn), red);

label("$2T_0-A_1$", (g/2,0));
clip(shift((g/2,0)) * scale(1.5) * box((-g,-h),(g,h)));
//clip(shift((g/2,0)) * scale(5) * box((-kf,-h),(kf,h)));
