import graph;
settings.outformat = "pdf";
unitsize(1cm);
real g = 0.2;
real gy = sqrt(3)/2*g;
pair g1 = (1g, 0);
pair g2 = (g/2, gy);

real w = 5;
real wr = 1;
real l = 10;
real lr = 2.5;

pair z(int i, int j) {
    if (j >= 0) {
        return i*g1 + j*g2 - (int)(j/2)*(g,0);
    } else {
        return i*g1 + j*g2 - (int)((j-1)/2)*(g,0);
    }
}
{
    pen pb = 0.8*white + 3;
    for (int i = -12; i <= 12; ++i) {
    for (int j = -14; j <= 14; ++j) {
        if (i == 12 && (j % 2 == 1)) { continue; }
        dot(z(i,j) + (l/2,w/2), pb);
        //fill(circle(z(i,j) + (l/2,w/2), 0.05), pb);
    }}
}
{
    pen pen = black + 2;
    draw((0, w)--(l, w), pen);
    draw((0, w-wr)--(lr, w-wr)--(lr,wr)--(0,wr), pen);
    draw((l, w-wr)--(l-lr, w-wr)--(l-lr,wr)--(l,wr), pen);
    draw((0, 0)--(l, 0), pen);
    path p =
        (0, w)--(l, w)--
        (l, w-wr)--(l-lr, w-wr)--(l-lr,wr)--(l,wr)--
        (l,0)--(0,0)--
        (0, wr)--(lr, wr)--(lr,w-wr)--(0,w-wr)--
        cycle;
    fill(p, pen + opacity(0.1));
    draw((3lr/2,wr/2)--(l-3lr/2,wr/2), black+1, Arrow(size=15));
    label("\Large $J_{12}$", (l/2, wr/2), N);

}
label("\Large\tt 4", (lr/2,w-wr/2));
label("\Large\tt 3", (l-lr/2,w-wr/2));
label("\Large\tt 2", (l-lr/2,wr/2));
label("\Large\tt 1", (lr/2,wr/2));
{
    path b = (lr - wr/2, w - wr/2)--(lr - wr/2, w + wr/2)--
    (l - lr + wr/2, w + wr/2)--(l - lr + wr/2, w - wr/2);
    draw(b, black+1);
    fill(circle(point(b,0.0), 0.1));
    fill(circle(point(b,3.0), 0.1));
    label("\Large $V_{43}$", (l/2,w+wr/2), N);
}
