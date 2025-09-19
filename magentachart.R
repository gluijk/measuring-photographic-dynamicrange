# Create magenta chart for Camera RAW Dynamic range measurement
# www.overfitting.net
# https://www.overfitting.net/2025/07/rango-dinamico-de-un-sensor-de-imagen.html

# library(tiff)
library(png)


################################

# CREATE MAGENTA CHART

# FUll HD chart
DIMX=1920
DIMY=1080

NCOLS=11
NROWS=7

MARGIN=0.7  # allow for outward space
WIDTH=round( min(floor(DIMX/NCOLS), floor(DIMY/NROWS)) * MARGIN)
HEIGHT=WIDTH

DIMXc=WIDTH*NCOLS
DIMYc=HEIGHT*NROWS
chart=array(0, dim=c(DIMYc, DIMXc, 3))

f=1.4  # nonlinear colur scale factor (inverse gamma)
for (f in seq(from=1, to=3, by=0.1)) {
    val=seq(0, 1, length.out=NCOLS*NROWS)^f
    p=1
    for (j in 1:NROWS) {
        for (i in 1:NCOLS) {
            x1=(i-1)*WIDTH+1
            x2= i   *WIDTH
            y1=(j-1)*HEIGHT+1
            y2= j   *HEIGHT
            patch=which(row(chart[,,1])>=y1 & row(chart[,,1])<=y2 &
                        col(chart[,,1])>=x1 & col(chart[,,1])<=x2)
            
            # UniWB for Canon 350D: R=162, G=64 y B=104
            chart[,,1][patch]=val[p]*162/162  # R
            chart[,,2][patch]=val[p]* 64/162  # G
            chart[,,3][patch]=val[p]*104/162  # B
            p=p+1
        }
    }
    
    R=15  # white circles radius
    DIMX=DIMXc + 4*WIDTH + 2*R
    DIMY=DIMYc + 4*HEIGHT + 2*R
    chartfinal=array(0, dim=c(DIMY, DIMX, 3))
    
    OFFSETX=round((DIMX-DIMXc)/2)
    OFFSETY=round((DIMY-DIMYc)/2)
    
    chartfinal[(OFFSETY+1):(OFFSETY+DIMYc), (OFFSETX+1):(OFFSETX+DIMXc),]=chart
    
    # Now we add 4 white circles
    x0=c(R+1+WIDTH,  R+1+WIDTH,     DIMX-R-WIDTH, DIMX-R-WIDTH)
    y0=c(R+1+HEIGHT, DIMY-R-HEIGHT, R+1+HEIGHT,   DIMY-R-HEIGHT)
    for (i in 1:4) {
        indices=which( ((row(chartfinal[,,1])-y0[i])/R)^2 +
                       ((col(chartfinal[,,1])-x0[i])/R)^2 < 1 )
        chartfinal[,,1][indices]=1
        chartfinal[,,2][indices]=1 
        chartfinal[,,3][indices]=1 
    }
    
    # Write chart
    # writeTIFF(chart, "magentachart.tif", bits.per.sample=16)
    writePNG(chartfinal, paste0("magentachart_", f*10, ".png"))
}



