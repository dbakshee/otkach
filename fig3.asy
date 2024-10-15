import graph;
settings.outformat = "pdf";
unitsize(1.3cm);
real kf = 1.4;
real g = 1;

int n = 7;

pair[] z;
for (int i = 0; i < n; ++i) {  z[i] = (i*g, 0);  }

path[] c;
for (int i = 0; i < n; ++i) {  c[i] = circle(z[i], kf);  }

pen pb = 0.8*white;
for (int i = 0; i < n; ++i) {  draw(c[i], pb); dot(z[i], pb); }
{ pair zz=(-1g,0); draw(circle(zz,kf), pb); dot(zz, pb);}
//------------------------------------------------------------------------------
{
    real a = acos(g/2/kf) / pi * 180;
    pen p = 0.6*blue;
    draw(arc(z[1],kf,a,360-a), p+1.3, Arrow(size=10, position=2.5));
    draw(arc(z[2],kf,180+a,180-a), p+1.3, Arrow(size=10, position=1.5));
    label("\Large$A_0{-}A_1$", z[1]+(0,kf), N, p);
}

{
    real a = acos(g/kf) / pi * 180;
    pen p = 0.6*blue;
    draw(arc(z[3],kf,a,360-a), p+1.3, Arrow(size=10, position=2.5));
    draw(arc(z[5],kf,180+a,180-a), p+1.3, Arrow(size=10, position=1.5));
    label("\Large$A_0{-}A_2$", z[3]+(0,kf), N, p);
}

label("\Large\bf b)", (-0.9kf,1.5kf), SE);
limits((-0.9kf,-1.1kf),(4.5g,1.5kf), Crop);
