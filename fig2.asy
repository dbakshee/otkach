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
    pen p = 0.6*blue;
    draw(arc(z[0],kf,0,180), p+1.3, Arrow(size=10, position=0.3));
    draw(arc(z[0],kf,180,360), p+1.3, Arrow(size=10, position=1.3));
    label("\Large$A_0$", z[0]+(0,kf), N, p);
}

{
    real a = acos(g/2/kf) / pi * 180;
    pen p = 0.6*blue;
    draw(arc(z[2],kf,-a,+a), p+1.3, Arrow(size=10, position=1.5));
    draw(arc(z[3],kf,180-a,180+a), p+1.3, Arrow(size=10, position=1.5));
    label("\Large$A_1$", (z[2]+z[3])/2+(0,kf), N, p);
}

{
    real a = acos(g/kf) / pi * 180;
    pen p = 0.6*blue;
    draw(arc(z[3],kf,-a,+a), p+1.3, Arrow(size=10, position=1.5));
    draw(arc(z[5],kf,180-a,180+a), p+1.3, Arrow(size=10, position=1.5));
    label("\Large$A_2$", (z[3]+z[5])/2+(0,kf), N, p);
}

label("\Large\bf a)", (-1.1kf,1.5kf), SE);
limits((-1.1kf,-1.1kf),(4.5g,1.5kf), Crop);
