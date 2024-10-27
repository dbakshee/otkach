import graph;
settings.outformat = "pdf";
unitsize(1.8cm);
real kf = 0.8;
real g = 1;
real gy = sqrt(3)/2*g;
pair g1 = (1g, 0);
pair g2 = (g/2, gy);

pair z(int i, int j) {
    if (j >= 0) {
        return i*g1 + j*g2 - (int)(j/2)*(g,0);
    } else {
        return i*g1 + j*g2 - (int)((j-1)/2)*(g,0);
    }
}

pen pb = 0.8*white;
for (int i = -1; i < 4; ++i) {
    for (int j = -1; j < 4; ++j) {
        draw(circle(z(i,j), kf), pb);
        dot(z(i,j), pb); }}
//------------------------------------------------------------------------------
{
    string label = "$A_0{-}A_1$";
    pen pen = red;
    pair zl = z(0,0);
    pair zr = z(1,0);
    path l = circle(zl, kf);
    path r = circle(zr, kf);
    real[][] lr = intersections(l,r);
    path pa = subpath(l, lr[0][0], lr[1][0]);
    path pb = subpath(r, lr[0][1], lr[1][1]);
    fill(pa..reverse(pb)..cycle, pen + opacity(0.2));
    draw(pa, Arrow(position=2.5));
    draw(pb, Arrow(position=1.5));
    label(label, zl - g1/4);
}
{
    string label = "$A_1{-}T_0$";
    pen pen = green;
    pair zl = z(1,0);
    pair zr = z(2,0);
    pair zd = z(1,-1);
    path l = circle(zl, kf);
    path r = circle(zr, kf);
    path d = circle(zd, kf);
    real[][] lr = intersections(l,r);
    real[][] ld = intersections(l,d);
    real[][] rd = intersections(r,d);
    path pa = subpath(d, rd[0][1], ld[0][1]);
    path pb = subpath(l, ld[0][0], lr[1][0]);
    path pc = subpath(r, lr[1][1], rd[0][0]);
    fill(pa..pb..pc..cycle, pen + opacity(0.2));
    draw(pa, Arrow(position=0.6));
    draw(pb, Arrow(position=1.6));
    draw(reverse(pc), Arrow(position=0.6));
    label(label, zd + 2g1/6, N);
}
{
    string label = "$2(A_1{-}T_0)$";
    pen pen = green;
    pair zl = z(1,0);
    pair zr = z(2,0);
    pair zu = z(1,1);
    path l = circle(zl, kf);
    path r = circle(zr, kf);
    path u = circle(zu, kf);
    real[][] lr = intersections(l,r);
    real[][] lu = intersections(l,u);
    real[][] ru = intersections(r,u);
    path pa = subpath(l, lr[0][0], lu[1][0]);
    path pb = subpath(u, lu[1][1], ru[1][1]);
    path pc = subpath(r, lr[0][1], ru[1][0]);
    path pd = subpath(l, lu[0][0], lr[0][0]);
    path pe = subpath(u, lu[0][1], ru[0][1]);
    path pf = subpath(r, ru[0][0], lr[0][1]);
    fill(pa..pb..reverse(pc)..cycle, pen + opacity(0.2));
    fill(reverse(pd)..pe..pf..cycle, pen + opacity(0.2));
    draw(pa, Arrow(position=0.7));
    draw(pb, Arrow(position=0.6));
    draw(pc, Arrow(position=0.6));
    draw(pd, Arrow(position=0.6));
    draw(pe, Arrow(position=0.6));
    draw(pf, Arrow(position=1.6));
    label(label, (zl + zr)/2 + g1/4);
}
{
    string label = "$A_0{-}2A_1{+}T_0$";
    pen pen = blue;
    pair zl = z(-1,3);
    pair zr = z(0,3);
    pair zd = z(0,2);
    path l = circle(zl, kf);
    path r = circle(zr, kf);
    path d = circle(zd, kf);
    real[][] lr = intersections(l,r);
    real[][] ld = intersections(l,d);
    real[][] rd = intersections(r,d);
    path pa = subpath(l, ld[0][0], lr[1][0]);
    path pb = subpath(r, lr[1][1], rd[1][0]);
    path pc = subpath(d, ld[0][1], rd[1][1]+4);
    fill(pa..pb..reverse(pc)..cycle, pen + opacity(0.2));
    draw(pa, Arrow(position=1.2));
    draw(pb, Arrow(position=1));
    draw(pc, Arrow(position=2));
    label(label, zd + (0,3g/7));
}
{
    string label = "$A_0{-}3A_1{+}2T_0$";
    pen pen = yellow;
    pair za = z(1,3); //label("a", za);
    pair zb = z(2,3); //label("b", zb);
    pair zc = z(2,2); //label("c", zc);
    pair zd = z(3,2); //label("d", zd);
    path a = circle(za, kf);
    path b = circle(zb, kf);
    path c = circle(zc, kf);
    path d = circle(zd, kf);
    real[][] ac = intersections(a,c);
    real[][] bc = intersections(b,c);
    real[][] dc = intersections(d,c);
    real[][] ab = intersections(a,b);
    real[][] bd = intersections(b,d);
    path pa = subpath(a, ac[0][0], ab[1][0]);
    path pb = subpath(b, ab[1][1], bd[0][0]);
    path pc = subpath(d, bd[0][1], dc[1][0]);
    path pd = subpath(c, ac[0][1], dc[1][1]);
    fill(pa..pb..pc..reverse(pd)..cycle, pen + opacity(0.5));
    draw(pa, Arrow(position=1.2));
    draw(pb);
    draw(pc, Arrow(position=1.5));
    draw(pd, Arrow(position=1.7));
    label(label, zc+(0,3g/7));
}

xaxis("$k_x$", YEquals(-gy), xmin=-1g, xmax=0, Arrow);
yaxis("$k_y$", XEquals(-1g), ymin=-gy, ymax=0, Arrow);

limits(z(-1,0) - (1.2kf,1.2kf),
       z(3,2) + (1.2kf,1.2kf), Crop);
