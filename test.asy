import graph;
settings.outformat = "pdf";
unitsize(1.8cm);
real kf = 0.8;
real g = 1;
pair g1 = (1g, 0);
pair g2 = (g/2, sqrt(3)/2*g);
int m = 5;
int n = 4;

pair z(int i, int j) {
    return (i-(int)(j/2))*g1 + j*g2;
}

pen pb = 0.8*white;
for (int i = -1; i < m + 1; ++i) {
    for (int j = -1; j < n + 1; ++j) {
        draw(circle(z(i,j), kf), pb);
        dot(z(i,j), pb); }}
//------------------------------------------------------------------------------
{
    pen pen = red;
    pair zc = z(0,1);
    fill(arc(zc, kf, 0, 360)..cycle, pen + opacity(0.2));
    draw(arc(zc, kf, 0, 360), Arrow(position=1.5));
    label("$A_0$", zc,  0.5*pen);
}
//
{
    pen pen = green;
    pair zl = z(1,1);
    pair zr = z(2,1);
    path l = circle(zl, kf);
    path r = circle(zr, kf);
    real tlb = intersections(l,r)[1][0]*90;
    real tle = intersections(l,r)[0][0]*90 + 360;
    real trb = intersections(r,l)[0][0]*90;
    real tre = intersections(r,l)[1][0]*90;
    fill(
        arc(zl, kf, tlb, tle)..
        arc(zr, kf, trb, tre)..
        cycle, pen + opacity(0.2));
    draw(arc(zl, kf, tlb, tle));
    draw(arc(zr, kf, trb, tre), Arrow(position=0.6));
    label("$A_1$", (zl + zr)/2, 0.5*pen);
}
//
{
    pen pen = blue;
    pair zu = z(1,1);
    pair zl = z(1,0);
    pair zr = z(2,0);
    path u = circle(zu, kf);
    path l = circle(zl, kf);
    path r = circle(zr, kf);
    real tub = intersections(u,r)[0][0]*90;
    real tue = intersections(u,l)[1][0]*90;
    real tlb = intersections(l,u)[0][0]*90;
    real tle = intersections(l,r)[0][0]*90;
    real trb = intersections(r,l)[0][0]*90;
    real tre = intersections(r,u)[1][0]*90;
    fill(
        arc(zu, kf, tub, tue)..
        arc(zl, kf, tlb, tle)..
        arc(zr, kf, trb, tre)..
        cycle, pen +  opacity(0.2));
    draw(arc(zu, kf, tub, tue));
    draw(arc(zl, kf, tlb, tle));
    draw(arc(zr, kf, trb, tre), Arrow(position=0.6));
    real y = sqrt(g**2 - g**2/4)/3;
    label("$T_0$", (zl+zr)/2 + (0, y), 0.5*pen);
}
//
{
    pen pen = orange;
    pair zu = z(0,3);
    pair zl = z(0,2);
    pair zr = z(1,2);
    path u = circle(zu, kf);
    path l = circle(zl, kf);
    path r = circle(zr, kf);
    real tub = intersections(u,l)[0][0]*90;
    real tue = intersections(u,r)[1][0]*90;
    real tlb = intersections(l,r)[0][0]*90;
    real tle = intersections(l,u)[1][0]*90;
    real trb = intersections(r,u)[0][0]*90;
    real tre = intersections(r,l)[0][0]*90;
    fill(
        arc(zu, kf, tub, tue)..
        arc(zr, kf, trb, tre)..
        arc(zl, kf, tlb, tle)..cycle,
        pen + opacity(0.2));
    draw(arc(zu, kf, tub, tue));
    draw(arc(zr, kf, trb, tre));
    draw(arc(zl, kf, tlb, tle), Arrow(position=0.5));
    real y = sqrt(g**2 - g**2/4)/2;
    label("$2A_1-T_0$", (zl+zr)/2 + (0, y), 0.5*pen);
}
//
{
    pen pen = cyan;
    pair zu = z(2,3);
    pair zl = z(2,2);
    pair zr = z(3,2);
    path u = circle(zu, kf);
    path l = circle(zl, kf);
    path r = circle(zr, kf);
    real tubl = intersections(u,l)[0][0]*90;
    real tuel = intersections(u,r)[0][0]*90;
    real trbc = intersections(r,u)[1][0]*90;
    real trec = intersections(r,l)[1][0]*90;
    real tlbc = intersections(l,r)[1][0]*90;
    real tlec = intersections(l,u)[0][0]*90;
    real tubr = intersections(u,l)[1][0]*90;
    real tuer = intersections(u,r)[1][0]*90;
    real tlb = intersections(l,r)[0][0]*90;
    real tle = intersections(l,u)[1][0]*90;
    real trb = intersections(r,u)[0][0]*90;
    real tre = intersections(r,l)[0][0]*90;
    // draw(arc(zu, kf, tubl, tuel), Arrow);
    // draw(arc(zr, kf, trbc, trec), Arrow);
    // draw(arc(zl, kf, tlbc, tlec, CCW), Arrow);
    // draw(arc(zu, kf, tubr, tuer), Arrow);
    // draw(arc(zr, kf, trb, tre), Arrow);
    // draw(arc(zl, kf, tlb, tle), Arrow);
    fill(
        arc(zu, kf, tubl, tuel)..
        arc(zr, kf, trbc, trec)..
        arc(zl, kf, tlbc, tlec, CCW)..
        arc(zu, kf, tubr, tuer)..
        arc(zr, kf, trb, tre)..
        arc(zl, kf, tlb, tle)..
        cycle,
        pen + opacity(0.2));
    draw(arc(zu, kf, tubl, tuel));
    draw(arc(zr, kf, trbc, trec));
    draw(arc(zl, kf, tlbc, tlec, CCW));
    draw(arc(zu, kf, tubr, tuer));
    draw(arc(zr, kf, trb, tre));
    draw(arc(zl, kf, tlb, tle), Arrow(position=0.5));
    real y = sqrt(g**2 - g**2/4)/2;
    label("$3A_1-2T_0$", (zl+zr)/2 + (0, y), 0.5*pen);
}
{
    pen pen = yellow;
    pair zl = z(3,1);
    pair zr = z(4,1);
    path l = circle(zl, kf);
    path r = circle(zr, kf);
    real tlb = intersections(l,r)[0][0]*90;
    real tle = intersections(l,r)[1][0]*90;
    real trb = intersections(r,l)[1][0]*90;
    real tre = intersections(r,l)[0][0]*90 + 360;
    fill(
        arc(zl, kf, tlb, tle)..
        arc(zr, kf, trb, tre)..
        cycle,
        pen + opacity(0.2));
    draw(arc(zl, kf, tlb, tle), Arrow(position=1));
    draw(arc(zr, kf, trb, tre), Arrow(position=1));
    label("$2A_0-A_1$", (zl+zr)/2, 0.3*pen);
}
{
    pen pen = magenta;
    pair zl = z(4,2);
    pair zr = z(4,3);
    path l = circle(zl, kf);
    path r = circle(zr, kf);
    real trb = intersections(r,l)[0][0]*90;
    real tre = intersections(r,l)[1][0]*90 + 360;
    real tlb = intersections(l,r)[0][0]*90;
    real tle = intersections(l,r)[1][0]*90;
    fill(
        arc(zl, kf, tlb, tle)..
        arc(zr, kf, trb, tre)..
        cycle,
        pen + opacity(0.4));
    fill(r,
        pen + opacity(0.2));
    draw(arc(zr, kf, 0, 360), Arrow(position=1));
    draw(arc(zr, 0.97*kf, trb, tre), Arrow(position=0.3));
    draw(arc(zl, kf, tlb, tle), Arrow(position=0.5));
    label("$A_0+A_1$", zr, N, 0.5*pen);
}

limits(z(0,1) - (1.1kf,1.1kf),
       z(4,3) + (1.1kf,1.1kf), Crop);
