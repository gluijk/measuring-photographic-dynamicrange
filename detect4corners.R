# Detect 4 corners
# www.overfitting.net
# https://www.overfitting.net/2025/07/rango-dinamico-de-un-sensor-de-imagen.html

# library(tiff)
library(png)


################################

# img is a grayscale matrix
name="magentachart_blur_noise.png"
chart=readPNG(name)
DIMX=ncol(chart)
DIMY=nrow(chart)

for (sector in 1:4) {
    if (sector==1) img=chart[1:(DIMY/2), 1:(DIMX/2)]
    if (sector==2) img=chart[1:(DIMY/2), (DIMX/2+1):DIMX]
    if (sector==3) img=chart[(DIMY/2+1):DIMY, (DIMX/2+1):DIMX]
    if (sector==4) img=chart[(DIMY/2+1):DIMY, 1:(DIMX/2)]
    
    # 1. Threshold for top 0.05% brightest pixels
    q <- quantile(img, probs = 0.9995)
    
    # 2. Coordinates of pixels above threshold
    coords <- which(img >= q, arr.ind = TRUE)
    img[coords]=1-img[coords]
    
    # 3. Median coordinates, rounded
    center_y <- round(median(coords[, 1]))
    center_x <- round(median(coords[, 2]))
    
    img[center_y, 1:(DIMX/2)]=1-img[center_y, 1:(DIMX/2)]
    img[1:(DIMY/2), center_x]=1-img[1:(DIMY/2), center_x]
    
    if (sector==1) chart[1:(DIMY/2), 1:(DIMX/2)]=img
    if (sector==2) chart[1:(DIMY/2), (DIMX/2+1):DIMX]=img
    if (sector==3) chart[(DIMY/2+1):DIMY, (DIMX/2+1):DIMX]=img
    if (sector==4) chart[(DIMY/2+1):DIMY, 1:(DIMX/2)]=img
}


# Write chart
writePNG(chart, paste0(name, "_detection.png"))




