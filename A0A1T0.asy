settings.outformat = "pdf";
unitsize(1.5cm);
real g = 3; // lattice
real kf = g * 2.2; // k fermi
real s32 = sqrt(3) / 2;
real h = g * s32;
real s = -g * s32 + sqrt(4 * kf^2 - g^2) / 2;
real sx = s32 * s;
real sy = s / 2;
real an = atan2(g/2 + sx, h + sy) / pi * 180; // half-angle top arc
real bn = atan2(sy, g + sx) / pi * 180; // half-angle right arc
pen tr = opacity(0.05); // transparent
real anA0 = 0;
real WA1=15;//8.5;
real WT0=20;
path A0 =
    arc((0, 0), kf, angle1=0, angle2=360) -- cycle;

path A0_A11 =
    arc((WA1, 0), kf, angle1=0, angle2=360) -- cycle;
path A0_A12 =
    arc((WA1+g, 0), kf, angle1=0, angle2=360) -- cycle;
path A0_T01 =
    arc((WT0, h/2), kf, angle1=0, angle2=360) -- cycle;
path A0_T02 =
    arc((WT0+g, h/2), kf, angle1=0, angle2=360) -- cycle;
path A0_T03 =
    arc((WT0+g/2, -h/2), kf, angle1=0, angle2=360) -- cycle;
//path twoT0minusA1 =
//    arc((g/2, -h), kf, angle1=90+an, angle2=90-an) --
//    arc((0,0), kf, angle1=bn, angle2=-bn) --
//    arc((g/2,  h), kf, angle1=-90+an, angle2=-90-an) --
//    arc((g,0), kf, angle1=180+bn, angle2=180-bn) -- cycle;

// draw

dot((0,0));

dot((WA1,0));
dot((WA1+g,0));

dot((WT0,h/2));
dot((WT0+g,h/2));
dot((WT0+g/2,-h/2));
//dot((g,0));
//dot((g/2,h));
//dot((g/2,-h));

filldraw(circle((0,0), kf), red+tr, red+white + dashed);

filldraw(circle((WA1,0), kf), red+tr, red+white + dashed);
filldraw(circle((WA1+g,0), kf), red+tr, red+white + dashed);

fill(A0, red+opacity(0.2));

fill(A0_A11, red+opacity(0.2));
fill(A0_A12, red+opacity(0.2));

fill(A0_T01, red+opacity(0.2));
fill(A0_T02, red+opacity(0.2));
fill(A0_T03, red+opacity(0.2));

draw(arc((0,0), r=kf, angle1=180, angle2=0), red, arrow=Arrow(HookHead));
draw(arc((0,0), r=kf, angle1=0, angle2=-180), red, arrow=Arrow(HookHead));

draw(arc((WA1,0), r=kf, angle1=180, angle2=0), red, arrow=Arrow(HookHead));
draw(arc((WA1,0), r=kf, angle1=0, angle2=-180), red, arrow=Arrow(HookHead));

draw(arc((WA1+g,0), r=kf, angle1=180, angle2=0), red, arrow=Arrow(HookHead));
draw(arc((WA1+g,0), r=kf, angle1=0, angle2=-180), red, arrow=Arrow(HookHead));


//T0
//draw(arc((WT0,h/2), r=kf, angle1=120, angle2=-30), red, arrow=Arrow(HookHead));
draw(arc((WT0,h/2), r=kf, angle1=-30, angle2=-390), red, arrow=Arrow(HookHead));

draw(arc((WT0+g,h/2), r=kf, angle1=210, angle2=-150), red, arrow=Arrow(HookHead));
//draw(arc((WT0+g,h/2), r=kf, angle1=60, angle2=-240), red, arrow=Arrow(HookHead));

draw(arc((WT0+g/2,-h/2), r=kf, angle1=90, angle2=-270), red, arrow=Arrow(HookHead));

label("\large $A_0$", (0,0), NE, 0.5*blue);
label("\large $A_1$", (WA1+g/2,0),  0.5*blue);
label("\large $T_0$", (WT0+g/2,0), N, 0.5*blue);
clip(shift((g/2,0)) * scale(5) * box((-g,-h),(WT0,h)));
