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

for (int i = 0; i < n; ++i) {  draw(c[i]); dot(z[i]); }
draw(circle((-g,0), kf));
dot((-g,0));


draw("$g$",z[1]--z[2], Arrows(HookHead));

picture a1;
fill(a1, c[0], red + opacity(0.2));
clip(a1, c[1]);
add(a1);
label("$A_1$", (g/2,0.7kf), red);

picture a2;
fill(a2, c[1], green + opacity(0.2));
clip(a2, c[3]);
add(a2);
label("$A_2$", (2g,kf/2), 0.5*green);

fill(c[4], blue + opacity(0.2));
label("$A_0$", (4g,0), N, 0.5*blue);

xaxis("$k_x$", YEquals(-1.2kf), xmax=5.5g, Arrow);
yaxis("$k_y$", XEquals(2.5g), ymin=-1.4kf, ymax=1.3kf, Arrow);

limits((-0.5g,-1.5kf),(5.5g,1.3kf), Crop);
