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

for (int i = 0; i < n; ++i) {  draw(c[i], 0.8white); dot(z[i], 0.8white); }


// draw("$g$",z[1]--z[2], Arrows(HookHead));

// picture a1;
// fill(a1, c[0], red + opacity(0.2));
// clip(a1, c[1]);
// add(a1);
//label("$A_1$", (g/2,0.7kf), red);

// picture a2;
// fill(a2, c[1], green + opacity(0.2));
// clip(a2, c[3]);
// add(a2);
//label("$A_2$", (2g,kf/2), 0.5*green);

//fill(c[4], blue + opacity(0.2));

//xaxis("$k_x$", YEquals(-1.2kf), xmax=5.5g, Arrow);
//yaxis("$k_y$", XEquals(2.5g), ymin=-1.4kf, ymax=1.3kf, Arrow);

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

limits((-1.1kf,-1.1kf),(4.5g,1.5kf), Crop);
