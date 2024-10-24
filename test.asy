import graph;
settings.outformat = "pdf";
unitsize(1.8cm);
real kf = 0.8;
real g = 1;
pair g1 = (1g, 0);
pair g2 = (g/2, sqrt(3)/2*g);
int m = 5;
int n = 4;

pair[] z;
for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
        z[i*n + j] = (i-(int)(j/2))*g1 + j*g2; }}

path[] c;
for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
        c[i*n + j] = circle(z[i*n + j], kf); }}

pen pb = 0.8*white;
for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
        draw(c[i*n + j], pb);
        dot(z[i*n + j], pb); }}
//------------------------------------------------------------------------------
{
    picture a0;
    fill(a0, c[0*n+1], red + opacity(0.2));
    add(a0);
    label("$A_0$", z[0*n+1],  0.7*red);
}
//
{
    picture a1;
    fill(a1, c[1*n+1], green + opacity(0.2));
    clip(a1, c[2*n+1]);
    add(a1);
    label("$A_1$", (z[1*n+1]+z[2*n+1])/2, 0.5*green);
}
//
{
    picture t0;
    fill(t0, c[1*n+1], blue + opacity(0.2));
    clip(t0, c[1*n+0]);
    clip(t0, c[2*n+0]);
    add(t0);
    real y = sqrt(g**2 - g**2/4)/3;
    label("$T_0$", (z[1*n+0]+z[2*n+0])/2 + (0, y), 0.5*blue);
}
//
{
    pen pen = orange;
    path u = c[0*n+3];
    path l = c[0*n+2];
    path r = c[1*n+2];
    real tub = intersections(u,l)[0][0]*90;
    real tue = intersections(u,r)[1][0]*90;
    real tlb = intersections(l,r)[0][0]*90;
    real tle = intersections(l,u)[1][0]*90;
    real trb = intersections(r,u)[0][0]*90;
    real tre = intersections(r,l)[0][0]*90;
    fill(
        arc(z[0*n+3], kf, tub, tue)..
        arc(z[1*n+2], kf, trb, tre)..
        arc(z[0*n+2], kf, tlb, tle)..cycle,
        pen + opacity(0.2));
    real y = sqrt(g**2 - g**2/4)/2;
    label("$2A_1-T_0$", (z[0*n+2]+z[1*n+2])/2 + (0, y), 0.5*pen);
}
