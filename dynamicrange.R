# Camera RAW Dynamic range measurement
# www.overfitting.net
# https://www.overfitting.net/2025/07/rango-dinamico-de-un-sensor-de-imagen.html

library(tiff)
library(Cairo)
library(Rcpp)
# library(microbenchmark)


################################

# KEYSTONE CORRECTION (SAME FOR ALL IMAGES)

# Undo distortion function
calculate_keystone = function(xu, yu, xd, yd, DIMX, DIMY) {
    # Calculate the k parameters that define a keystone correction
    # from 4 pairs of (xu,yu) coords + 4 pairs of (xd,yd) coords
    
    # NOTE: we swap the distorted and undistorted trapezoids because
    # we want to model the transformation
    # FROM CORRECTED coords (DST) -> TO UNCORRECTED coords (ORG)
    # DIMX, DIMY not needed for calculations, just to plot the trapezoids
    
    # Solve 8 equations linear system: A * k = b -> k = inv(A) * b
    A=matrix(nrow=8, ncol=8)
    A[1,]=c(xd[1], yd[1], 1, 0,     0,     0, -xd[1]*xu[1], -yd[1]*xu[1])
    A[2,]=c(0,     0,     0, xd[1], yd[1], 1, -xd[1]*yu[1], -yd[1]*yu[1])
    A[3,]=c(xd[2], yd[2], 1, 0,     0,     0, -xd[2]*xu[2], -yd[2]*xu[2])
    A[4,]=c(0,     0,     0, xd[2], yd[2], 1, -xd[2]*yu[2], -yd[2]*yu[2])
    A[5,]=c(xd[3], yd[3], 1, 0,     0,     0, -xd[3]*xu[3], -yd[3]*xu[3])
    A[6,]=c(0,     0,     0, xd[3], yd[3], 1, -xd[3]*yu[3], -yd[3]*yu[3])
    A[7,]=c(xd[4], yd[4], 1, 0,     0,     0, -xd[4]*xu[4], -yd[4]*xu[4])
    A[8,]=c(0,     0,     0, xd[4], yd[4], 1, -xd[4]*yu[4], -yd[4]*yu[4])
    b=as.matrix(c(xu[1], yu[1], xu[2], yu[2], xu[3], yu[3], xu[4], yu[4]))
    k=solve(A, b)  # equivalent to inv(A) * b = solve(A) %*% b
    
    # Check keystone correction
    for (i in 1:4) print(undo_keystone_coords(xd[i], yd[i], k))
    
    # Plot trapezoids
    plot(c(xd, xd[1]), c(yd, yd[1]), type='l', col='red', asp=1,
         xlab='X', ylab='Y', xlim=c(1,DIMX), ylim=c(DIMY,1))
    lines(c(xu, xu[1]), c(yu, yu[1]), type='l', col='blue')
    for (i in 1:4) {
        lines(c(xd[i], xu[i]), c(yd[i], yu[i]), type='l', lty=3, col='darkgray')
    }
    abline(h=c(1,DIMY), v=c(1,DIMX))
    
    return (k)
}

undo_keystone_coords = function(xd, yd, k) {
    # Return the keystone correction of a pair (or list of pairs) of coords
    denom=k[7]*xd + k[8]*yd + 1
    xu=(k[1]*xd + k[2]*yd + k[3]) / denom
    yu=(k[4]*xd + k[5]*yd + k[6]) / denom
    return(c(xu, yu))  # return pair (xu, yu)
}

undo_keystone = function(imgd, k) {
    # Return the keystone correction of an image (matrix)
    DIMX=ncol(imgd)
    DIMY=nrow(imgd)
    imgc=imgd*0
    
    for (x in 1:DIMX) {
        for (y in 1:DIMY) {
            denom=k[7]*x + k[8]*y + 1
            xu=round((k[1]*x + k[2]*y + k[3]) / denom)
            yu=round((k[4]*x + k[5]*y + k[6]) / denom)
            
            if (xu>=1 & xu<=DIMX & yu>=1 & yu<=DIMY)
                imgc[y, x]=imgd[yu, xu]  # nearest neighbour interp
        }
    }
    
    return(imgc)
}

cppFunction('
NumericMatrix undo_keystone_cpp(NumericMatrix imgd, NumericVector k) {
  // Return the keystone correction of an image (matrix)
  int DIMX = imgd.ncol();  // cols = X
  int DIMY = imgd.nrow();  // rows = Y
  NumericMatrix imgc(DIMY, DIMX);  // output image, same size

  for (int x = 0; x < DIMX; x++) {
    for (int y = 0; y < DIMY; y++) {
      double denom = k[6] * x + k[7] * y + 1;
      int xu = round((k[0] * x + k[1] * y + k[2]) / denom);
      int yu = round((k[3] * x + k[4] * y + k[5]) / denom);
      
      if (xu >= 0 && xu < DIMX && yu >= 0 && yu < DIMY) {
        imgc(y, x) = imgd(yu, xu);
      }
    }
  }

  return imgc;
}
')

cppFunction('
NumericMatrix undo_keystone_cpp_old(NumericMatrix imgd, NumericVector k) {
  // Return the keystone correction of an image (matrix)
  int DIMX = imgd.ncol();  // cols = X
  int DIMY = imgd.nrow();  // rows = Y
  NumericMatrix imgc(DIMY, DIMX);  // output image, same size

  for (int x = 0; x < DIMX; x++) {
    for (int y = 0; y < DIMY; y++) {
      double xp = x + 1;  // convert to 1-based indexing
      double yp = y + 1;  // (needed to use k values)

      double denom = k[6] * xp + k[7] * yp + 1;
      double xu = (k[0] * xp + k[1] * yp + k[2]) / denom;
      double yu = (k[3] * xp + k[4] * yp + k[5]) / denom;

      int xup = round(xu) - 1;  // back to 0-based indexing
      int yup = round(yu) - 1;  // (needed to write on imgc)

      if (xup >= 0 && xup < DIMX && yup >= 0 && yup < DIMY) {
        imgc(y, x) = imgd(yup, xup);
      }
    }
  }

  return imgc;
}
')


################################

# SNR CALCULATION OVER THE PATCHES

cppFunction('
List analyze_patches(NumericMatrix imgcrop, int NCOLS, int NROWS, double SAFE) {
  // Optimize SNR calculation over the patches
  int DIMX = imgcrop.ncol();
  int DIMY = imgcrop.nrow();
  
  std::vector<double> Signal;
  std::vector<double> Noise;

  for (int i = 1; i <= NCOLS; i++) {
    for (int j = 1; j <= NROWS; j++) {
      int x1 = round((i - 1) * DIMX / NCOLS + SAFE);
      int x2 = round(i * DIMX / NCOLS - SAFE);
      int y1 = round((j - 1) * DIMY / NROWS + SAFE);
      int y2 = round(j * DIMY / NROWS - SAFE);

      std::vector<double> values;

      for (int y = y1; y <= y2; y++) {
        for (int x = x1; x <= x2; x++) {
          if (y >= 0 && y < DIMY && x >= 0 && x < DIMX)
            values.push_back(imgcrop(y, x));
        }
      }

      if (values.size() == 0) continue;

      // Compute mean and stddev
      double sum = 0.0;
      for (double v : values) sum += v;
      double mean = sum / values.size();

      double var = 0.0;
      for (double v : values) var += (v - mean) * (v - mean);
      double stdev = sqrt(var / values.size());

      // Count saturated pixels
      int sat = 0;
      for (double v : values) if (v > 0.9) sat++;

      if (mean > 0 && 20 * log10(mean / stdev) >= -10 && ((double)sat / values.size()) < 0.01) {
        Signal.push_back(mean);
        Noise.push_back(stdev);

        // Draw black inner rectangle
        for (int y = y1; y <= y2; y++) {
          if (x1 >= 0 && x1 < DIMX) imgcrop(y, x1) = 0;
          if (x2 >= 0 && x2 < DIMX) imgcrop(y, x2) = 0;
        }
        for (int x = x1; x <= x2; x++) {
          if (y1 >= 0 && y1 < DIMY) imgcrop(y1, x) = 0;
          if (y2 >= 0 && y2 < DIMY) imgcrop(y2, x) = 0;
        }

        // Draw white outer rectangle
        for (int y = y1; y <= y2; y++) {
          if (x1-1 >= 0 && x1-1 < DIMX) imgcrop(y, x1-1) = 1;
          if (x2+1 >= 0 && x2+1 < DIMX) imgcrop(y, x2+1) = 1;
        }
        for (int x = x1; x <= x2; x++) {
          if (y1-1 >= 0 && y1-1 < DIMY) imgcrop(y1-1, x) = 1;
          if (y2+1 >= 0 && y2+1 < DIMY) imgcrop(y2+1, x) = 1;
        }
      }
    }
  }

  return List::create(
    Named("Signal") = Signal,
    Named("Noise") = Noise,
    Named("imgcrop") = imgcrop
  );
}
')


################################

# 0. BLACK AND SAT POINTS CALCULATION FROM DARKFRAME AND BLOWN RAW FILES

# It is mandatory to use accurate BLACK (especially) and SAT values
# to properly locate each pixel's level vs saturation

# dcraw -v -D -t 0 -4 -T BLACK.DNG SAT.DNG
imgblack=readTIFF("BLACK.tiff", as.is=TRUE)  # read unmodified integer RAW data
imgsat=readTIFF("SAT.tiff", as.is=TRUE)  # read unmodified integer RAW data

# BLACK level
for (rawchan in c('R', 'G1', 'G2', 'B')) {
    if (rawchan=='R') {
        imgBayerB=imgblack[row(imgblack)%%2 & col(imgblack)%%2]
        imgBayerS=imgsat[row(imgsat)%%2 & col(imgsat)%%2]        
    } else if (rawchan=='G1') {
        imgBayerB=imgblack[row(imgblack)%%2 & !col(imgblack)%%2]
        imgBayerS=imgsat[row(imgsat)%%2 & !col(imgsat)%%2]        
    } else if (rawchan=='G2') {
        imgBayerB=imgblack[!row(imgblack)%%2 & col(imgblack)%%2]
        imgBayerS=imgsat[!row(imgsat)%%2 & col(imgsat)%%2]        
    } else if (rawchan=='B') {
        imgBayerB=imgblack[!row(imgblack)%%2 & !col(imgblack)%%2]
        imgBayerS=imgsat[!row(imgsat)%%2 & !col(imgsat)%%2]        
    }
 
    print(paste0(rawchan, " BLACK: median=",
                 median(imgBayerB), ", mean=", mean(imgBayerB)))
    freq_table=table(as.vector(imgBayerB))
    write.csv2(freq_table, paste0("BLACK_freq_", rawchan, ".csv"))
    
    print(paste0(rawchan, " SAT: median=",
                 median(imgBayerS), ", mean=", mean(imgBayerS)))
    freq_table=table(as.vector(imgBayerS))
    write.csv2(freq_table, paste0("SAT_freq_", rawchan, ".csv"))
}

# The conclusion is that the best:
# BLACK is 255 (median) or 254.85 (floating point mean)
# SAT   is 4095 for all 3 channels
BLACK=mean(imgblack)  # 254.85


################################

# 1. EXTRACT RAW DATA AS 16-bit TIFF FILES READ THEM AND NORMALIZE

# dcraw -v -D -t 0 -4 -T *.DNG  # replace DNG for any camera RAW format

BLACK=256  # Olympus OM-1 black level. Unfortunately we used 256 instead of 255 or 254.85
SAT=4095  # Olympus OM-1 sat level

filepath=getwd()
filenames=list.files(path=filepath, pattern="\\.tiff$",  # pattern="\\.tif{1,2}$",
                      ignore.case=TRUE, full.names=FALSE)
filenamesISO=gsub(".tiff", "", filenames) 
filenamesISO=toupper(gsub("^iso0*", "iso", filenamesISO))

CairoPNG("SNRcurvesBLACK254.85.png", width=1920, height=1080)  # HQ Full HD curves

N=length(filenames)  # number of RAW files to process
for (image in 1:N) {
    NAME=filenames[image]
    cat(paste0('Processing "', NAME, '" ...\n'))
    
    # Read RAW data
    img=readTIFF(NAME, as.is=TRUE)  # read unmodified integer RAW data
    # hist(img, breaks=800, main=paste0('"', NAME, '" RAW histogram'))
    # abline(v=BLACK, col='red')
    
    # Normalize to floating point 0..1 range (negative values are allowed)
    img=img-BLACK
    img=img/(SAT-BLACK)
    # hist(img, breaks=800)
    # abline(v=0, col='red')

    
################################
    
# 2. EXTRACT INDIVIDUAL RAW CHANNEL(S) AND APPLY KEYSTONE CORRECTION
    
    # Keep one RAW channel (criteria checked OK):
    imgBayer=img[row(img)%%2 & col(img)%%2]  # R
    # imgBayer=img[row(img)%%2 & !col(img)%%2]  # G1
    # imgBayer=img[!row(img)%%2 & col(img)%%2]  # G2
    # imgBayer=img[!row(img)%%2 & !col(img)%%2]  # B
    dim(imgBayer)=dim(img)/2
    
    # Save uncorrected but normalized image
    # imgsave=imgBayer
    # imgsave[imgsave<0]=0
    # writeTIFF(imgsave, paste0("uncorrectednormalizedchart_", NAME, ".tif"), bits.per.sample=16)
    # rm(imgsave)

    # Correct keystone distortion
    DIMX=ncol(imgBayer)
    DIMY=nrow(imgBayer)
    if (image==1) {  # calculate the same k for all keystone corrections
        # The 4 points (xu,yu) should be calculated from chart (corner marks)
        # Distorted points (source)
        xu=c(119, 99, 2515, 2473)  # top-left, bottom-left, bottom-right, top-right
        yu=c(170, 1687, 1679, 158)
        
        # Undistorted points (destination):
        
        # top-left
        xtl=(xu[1] + xu[2]) / 2
        ytl=(yu[1] + yu[4]) / 2
        # bottom-right
        xbr=(xu[3] + xu[4]) / 2
        ybr=(yu[2] + yu[3]) / 2
        
        xd=c(xtl, xtl, xbr, xbr)
        yd=c(ytl, ybr, ybr, ytl)
        
        k=calculate_keystone(xu, yu, xd, yd, DIMX, DIMY)
    }
    
    imgc=undo_keystone_cpp(imgBayer, k)

    # benchmark_results=microbenchmark(
    #     Cpp_loops_old=undo_keystone_cpp_old(imgBayer, k),
    #     Cpp_loops=undo_keystone_cpp(imgBayer, k),
    #     times=30
    # )
    # benchmark_results
    # boxplot(benchmark_results)
    
    # Unit: milliseconds
    # expr        min          lq       mean     median         uq        max neval
    # R_loops 24111.0131 24253.91095 24385.2646 24273.2952 24365.5697 25671.9213    15
    # Cpp_loops    98.4014    98.96145   101.4232   100.1011   101.5062   116.5063    15
    
    # Save corrected image
    # imgsave=imgc
    # imgsave[imgsave<0]=0
    # writeTIFF(imgsave, paste0("correctedchart_", NAME, ".tif"), bits.per.sample=16)
    # rm(imgsave)
    
    
################################
    
# 3. READ PATCHES TO FORM 7x11 GRID AND COLLECT (EV,SNR) PAIRS
    
    # Crop patches area
    imgcrop=imgc[round(ytl):round(ybr), round(xtl):round(xbr)]

    # Analyze imgcrop dividing it in NCOLS x NROWS patches leaving a SAFE guard
    # DIMX=ncol(imgcrop)
    # DIMY=nrow(imgcrop)
    # for (i in 1:NCOLS) {
    #     for (j in 1:NROWS) {
    #         x1=round((i-1)*DIMX/NCOLS + SAFE)
    #         x2=round( i   *DIMX/NCOLS - SAFE)
    #         y1=round((j-1)*DIMY/NROWS + SAFE)
    #         y2=round( j   *DIMY/NROWS - SAFE)
    #         patch=which(row(imgcrop)>=y1 & row(imgcrop)<=y2 & 
    #                     col(imgcrop)>=x1 & col(imgcrop)<=x2)
    #         values=imgcrop[patch]
    #         S=mean(values)  # S=mean
    #         N=var(values)^0.5  # N=stdev
    #         
    #         # Ignore patches with negative average values, SNR < -10dB or
    #         # >1% of saturated/nonlinear (>90%) values
    #         if (S>0 & 20*log10(S/N) >= -10 & length(values[values>0.9])/length(values)<0.01) {
    #             Signal=c(Signal,S)
    #             Noise=c(Noise, N)
    #             
    #             # Draw patch rectangle
    #             imgcrop[y1:y2,x1]=0
    #             imgcrop[y1:y2,x2]=0
    #             imgcrop[y1,x1:x2]=0
    #             imgcrop[y2,x1:x2]=0
    #             
    #             imgcrop[y1:y2,(x1-1)]=1
    #             imgcrop[y1:y2,(x2+1)]=1
    #             imgcrop[(y1-1),x1:x2]=1
    #             imgcrop[(y2+1),x1:x2]=1
    #         }
    # 
    #     }
    # }
    calc=analyze_patches(imgcrop, NCOLS=11, NROWS=7, SAFE=50)  # SAFE=80
    Signal=calc$Signal
    Noise=calc$Noise
    SNR=Signal/Noise
    imgcrop=calc$imgcrop
    patches_used=length(Signal)

    # Order from lower to higher SNR values (to plot beautifully)
    neworder=order(SNR)  # SNR will be the independent variable in splines
    Signal=Signal[neworder]
    Noise=Noise[neworder]
    SNR=SNR[neworder]
    
    imgsave=imgcrop
    imgsave[imgsave<0]=0
    writeTIFF(imgsave, paste0("croppedchart_usedpatches_", NAME, ".tif"), bits.per.sample=16)
    rm(imgsave)
    
    # SNR curves in dB: we'll plot blue scatter points and red curves
    if (image==1) {
        plot(log2(Signal), 20*log10(SNR), xlim=c(-14,0), ylim=c(-10,20),
             pch=16, cex=0.5, col='blue',
             main='SNR curves - Olympus OM-1',
             xlab='RAW exposure (EV)', ylab='SNR (dB)')
        abline(h=c(0,12), lty=2)
        abline(v=seq(-14,0,1), lty=2, col='gray')
        axis(side=1, at=-14:0)
    } else {
        points(log2(Signal), 20*log10(SNR), pch=16, cex=0.5, col='blue')        
    }
    
    # # SNR curves in EV
    # plot(log2(Signal), log2(SNR), xlim=c(-14,0), ylim=c(-2,4), 
    #      pch=16, cex=0.5, col='blue',
    #      main=paste0('SNR curves\nOlympus OM-1 at ', NAME),
    #      xlab='RAW exposure (EV)', ylab='SNR (EV)')
    # abline(h=c(0,2), lty=2)
    # abline(v=seq(-14,0,1), lty=2, col='gray')

    
################################

# 4. APPROXIMATION CURVES TO CALCULATE DR VALUES
    
    # Soft cubic splines approximation
    # https://stackoverflow.com/questions/37528590/r-smooth-spline-smoothing-spline-is-not-smooth-but-overfitting-my-data
    # NOTES:
    #   Inverted variables: Signal = f(SNR) -> DR = f(SNR threshold)
    #   Splines already in log transformed domains (softer derivatives)
    spline_fit=smooth.spline(20*log10(SNR), log2(Signal),
                             spar=0.5, nknots=10)  # spar controls smoothness
    valuesx=predict(spline_fit, 20*log10(SNR))$y
    valuesy=20*log10(SNR)
    lines(valuesx, valuesy, col='red', type='l')
    text(max(valuesx), max(valuesy)+0.5, labels=filenamesISO[image], col='red')
    
    # Now calculate the Dynamic Range for a given SNR threshold criteria
    TH_dB=c(12, 0)  # Photographic (12dB) and Engineering (0dB) DR
    DR_EV=c()
    for (threshold in 1:length(TH_dB)) {
        DR_EV=c(DR_EV, -predict(spline_fit, TH_dB[threshold])$y)
        xpos=-DR_EV[threshold]
        ypos=TH_dB[threshold]
        points(xpos, ypos, pch=3, cex=1.5)
        text(xpos, ypos+0.5,
             labels=paste0(round(DR_EV[threshold],1), "EV"))
    }

    # Create dataframe with all calculated DR
    if (image==1) {
        DR_df=data.frame(tiff_file=NAME,
                         DR_EV_12dB=DR_EV[1], DR_EV_0dB=DR_EV[2],
                         patches_used=patches_used)        
    } else {
        new_row=data.frame(tiff_file=NAME,
                           DR_EV_12dB=DR_EV[1], DR_EV_0dB=DR_EV[2],
                           patches_used=patches_used)
        DR_df=rbind(DR_df, new_row)
    }
}

text(-14+0.1, 12+0.5, labels="Photographic DR (SNR>12dB)", adj=0)
text(-14+0.1,  0+0.5, labels="Engineering DR (SNR>0dB)", adj=0)

dev.off()

# Print calculated DR for each ISO
DR_df=DR_df[order(DR_df$tiff_file), ]
print(DR_df)


