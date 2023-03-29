# W03: 1201x125 - 1920x800
# W05: 1001x181 - 1600x1158
# W02: 2001x47 - 1600x150
# W02: 2001x80 - 1600x255

vals = "dos rxy21 rxx rxy12 r2t22 r2t11 ryy"

set terminal pngcairo  transparent enhanced font "arial,10" fontscale 1.0 size 1600, 255
unset key
unset colorbox
set view map
do for [val in vals] {
    bin = "w02/w02_".val.".bin"
    png = "w02/w02a_".val.".png"
    set output png
    plot bin matrix binary with image
}
