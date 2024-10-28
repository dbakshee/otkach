import graph;
settings.outformat = "pdf";
unitsize(1.8cm);
real kf = 0.8;
real g = 1;
real gy = sqrt(3)/2*g;
pair g1 = (1g, 0);
pair g2 = (g/2, gy);
int m = 5;
int n = 4;

pair z(int i, int j) {
    return (i-(int)(j/2))*g1 + j*g2;
}

pen pb = 0.8*white;
for (int i = -1; i < m + 2; ++i) {
    for (int j = -1; j < n + 1; ++j) {
        draw(circle(z(i,j), kf), pb);
        dot(z(i,j), pb); }}
//------------------------------------------------------------------------------
{
    string label = "$A_0$";
    pen pen = red;
    pair zc = z(0,1);
    fill(arc(zc, kf, 0, 360)..cycle, pen + opacity(0.2));
    draw(arc(zc, kf, 0, 360), Arrow(position=2.5));
    label(label, zc);
}
//
{
    string label = "$A_1$";
    pen pen = green;
    pair zl = z(2,0);
    pair zr = z(2,1);
    path l = circle(zl, kf);
    path r = circle(zr, kf);
    real tlb = intersections(l,r)[0][0]*90;
    real tle = intersections(l,r)[1][0]*90 + 360;
    real trb = intersections(r,l)[0][0]*90;
    real tre = intersections(r,l)[1][0]*90;
    fill(
        arc(zl, kf, tlb, tle)..
        arc(zr, kf, trb, tre)..
        cycle, pen + opacity(0.2));
    draw(arc(zl, kf, tlb, tle));
    draw(arc(zr, kf, trb, tre), Arrow(position=0.6));
    label(label, (zl + zr)/2);
}
{
    string label = "$A_{1a}$";
    pen pen = green;
    pair zu = z(2,4);
    pair zd = z(2,2);
    real kf1 = g*1.05;
    path u = circle(zu, kf1);
    path d = circle(zd, kf1);
    draw(u, 0.3white+dashed);
    draw(d, 0.3white+dashed);
    real tab = intersections(d,u)[0][0]*90;
    real tae = intersections(d,u)[1][0]*90 + 360;
    real tbb = intersections(u,d)[0][0]*90;
    real tbe = intersections(u,d)[1][0]*90;
    fill(
        arc(zd, kf1, tab, tae)..
        arc(zu, kf1, tbb, tbe)..
        cycle, pen + opacity(0.2));
    draw(arc(zd, kf1, tab, tae), 0.5pen+dashed);
    draw(arc(zu, kf1, tbb, tbe), 0.5pen+dashed, Arrow(position=0.6));
    label(label, (zu + zd)/2);
}
//
{
    string label = "$T_0$";
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
    label(label, (zl+zr)/2 + (0, y));
}
//
{
    string label = "$2A_1{-}T_0$";
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
    label(label, (zl+zr)/2 + (0, y));
}
//
{
    string label = "$3A_1{-}2T_0$";
    pen pen = cyan;
    pair zd = z(3,2);
    pair zl = z(2,3);
    pair zr = z(3,3);
    path d = circle(zd, kf);
    path l = circle(zl, kf);
    path r = circle(zr, kf);
    real[][] dl = intersections(d,l);
    real[][] dr = intersections(d,r);
    real[][] rl = intersections(r,l);
    path res =
        subpath(d, dr[1][0], dl[1][0])..
        subpath(l, dl[1][1], rl[1][1])..
        subpath(r, rl[1][0], dr[0][1])..
        subpath(d, dr[0][0], dl[0][0])..
        subpath(l, dl[0][1], rl[0][1]+4)..
        subpath(r, rl[0][0], dr[1][1])..
        cycle;
    fill(res, pen + opacity(0.2));
    draw(res, Arrow(position=0.6));
    real y = sqrt(g**2 - g**2/4)/2;
    label(label, (zl+zr)/2 - (0, y));
}
{
    string label = "$3A_1{-}2T_0$";
    pen pen = cyan;
    pair zo = z(2,2);
    pair za = z(1,2);
    pair zb = z(1,1);
    pair zc = z(2,1);
    path o = circle(zo, kf);
    path a = circle(za, kf);
    path b = circle(zb, kf);
    path c = circle(zc, kf);
    real tab = intersections(o,a)[0][0]*90;
    real tae = intersections(o,c)[1][0]*90;
    real tbb = intersections(c,o)[0][0]*90;
    real tbe = intersections(c,b)[0][0]*90;
    real tcb = intersections(b,c)[0][0]*90;
    real tce = intersections(b,a)[0][0]*90;
    real tdb = intersections(a,b)[1][0]*90;
    real tde = intersections(a,o)[0][0]*90;
    draw(arc(zo, kf, tab, tae), Arrow(position=0.6));
    draw(arc(zc, kf, tbb, tbe));
    draw(arc(zb, kf, tcb, tce));
    draw(arc(za, kf, tdb, tde, CCW));
    fill(
        arc(zo, kf, tab, tae)..
        arc(zc, kf, tbb, tbe)..
        arc(zb, kf, tcb, tce)..
        arc(za, kf, tdb, tde, CCW)..
        cycle, pen + opacity(0.2));
    real y = sqrt(g**2 - g**2/4)/2;
    label(label, (zb+zc)/2 + (0, y));
}
{
    string label = "$3A_1{-}2T_0$";
    pen pen = cyan;
    pair za = z(0,4);
    pair zb = z(0,3);
    pair zc = z(1,4);
    pair zd = z(1,3);
    path a = circle(za, kf);
    path b = circle(zb, kf);
    path c = circle(zc, kf);
    path d = circle(zd, kf);
    real tab = intersections(b,d)[0][0]*90;
    real tae = intersections(b,a)[1][0]*90;
    real tbb = intersections(a,b)[0][0]*90;
    real tbe = intersections(a,c)[1][0]*90;
    real tcb = intersections(c,a)[1][0]*90;
    real tce = intersections(c,d)[1][0]*90;
    real tdb = intersections(d,c)[0][0]*90;
    real tde = intersections(d,b)[0][0]*90;
    fill(
        arc(zb, kf, tab, tae)..
        arc(za, kf, tbb, tbe)..
        arc(zc, kf, tcb, tce)..
        arc(zd, kf, tdb, tde)..
        cycle, pen + opacity(0.2));
    draw(arc(zb, kf, tab, tae), Arrow(position=1.6));
    draw(arc(za, kf, tbb, tbe));
    draw(arc(zc, kf, tcb, tce), Arrow(position=1.6));
    draw(arc(zd, kf, tdb, tde));
    label(label, (zb+zc)/2);
}
{
    string label = "$2A_0{-}A_1$";
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
    draw(arc(zr, kf, trb, tre), Arrow(position=1.5));
    label(label, (zl+zr)/2);
}
{
    string label = "$A_0{+}A_1$";
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
    label(label, zr, N);
}

xaxis("$k_x$", YEquals(g*sqrt(3)), xmin=-0.5g, xmax=-0.05g, Arrow);
yaxis(XEquals(-0.5g), ymin=g*sqrt(3), ymax=2.5g, Arrow);
label("$k_y$", (-0.5g,2.5g), N);

//xaxis("$k_x$", YEquals(0), xmin=-1g, xmax=g/2, Arrow);
//yaxis("$k_y$", XEquals(0-0.5g), ymin=-gy, ymax=gy, Arrow);

limits(z(0,1) - (1.1g,0.95gy),
       z(4,3) + (1.05kf,1.2kf), Crop);
