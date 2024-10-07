settings.outformat = "pdf";
unitsize(0.8cm);
real h = 4;
real l = h;
real w = h / 5;
real c = 1;

fill((-c,h)--(l+c,h)--(l+c,h-w)--(-c,h-w)--cycle, 0.9*white);
fill((-c,0)--(l+c,0)--(l+c,w)--(-c,w)--cycle, 0.9*white);
fill((0,h)--(l,h)--(l,0)--(0,0)--cycle, 0.7*white);

draw((-c,h)--(l+c,h));
draw((-c,0)--(l+c,0));
draw((-c,w)--(0,w)--(0,h-w)--(-c,h-w));
draw((l+c,w)--(l,w)--(l,h-w)--(l+c,h-w));

int numv = 10;
real d = l / numv;
for (int i = 1; i <= numv; ++i) {
    real xi = 0 + i*d-d/2;
    draw((xi,h*0.01)--(xi,h*0.99), 0.6*white);
}
draw((l/2-2*d, w/2)--(l/2+2*d, w/2), 0.5*blue+1.5, Arrow);
label("1",(0,w/2),W);
label("2",(l,w/2),E);
label("4",(0,h-w/2),W);
label("3",(l,h-w/2),E);
label("\large $I_{12}$",(l/2,w/2),NW, 0.5*blue);
label("\large $\displaystyle R_{xx}=\frac{V_{34}}{I_{12}}$",(l/2,h/2),N, 0.5*blue);
