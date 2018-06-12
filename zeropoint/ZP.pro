PRO readfiles

  readcol, "tempfiles/files.txt",files,f='a'
  fileno = n_elements(files) - 1

  FOR i=0,fileno DO BEGIN

     readcol, files[i],$
     "IMAGE"+files[i],"XINIT"+files[i],"YINIT"+files[i],"ID"+files[i],$
              "COORDS"+files[i],"LID"+files[i],"sl1"+files[i],$
     "XCENTER"+files[i],"YCENTER"+files[i],"XSHIFT"+files[i],"YSHIFT"+files[i],$
              "XERR"+files[i],"YERR"+files[i],"CIER"+files[i],"CERROR"+files[i],"sl2"+files[i],$
     "MSKY"+files[i],"STDEV"+files[i],"SSKEW"+files[i],"NSKY"+files[i],$
              "NSREJ"+files[i],"SIER"+files[i],"SERROR"+files[i],"sl3"+files[i],$
     "ITIME"+files[i],"XAIRMASS"+files[i],"IFILTER"+files[i],"OTIME"+files[i],"sl4"+files[i],$
     "RAPERT"+files[i],"SUM"+files[i],"AREA"+files[i],"FLUX"+files[i],$
              "MAG"+files[i],"MERR"+files[i],"PIER"+files[i],"PERROR"+files[i],$
     f='a,f,f,i,a,i,a,f,f,f,f,a,a,i,a,a,f,f,f,i,i,i,a,a,f,f,a,a,a,f,f,f,f,f,f,i,a'

     ;print, COORDSB1[1]

  ENDFOR

END


PRO magnitude


print, "     Reading in files..."

  readcol,'b1.dat',$
   IMAGEB1,XINITB1,YINITB1,IDB1,COORDSB1,LIDB1,sl1B1,$
   XCENTERB1,YCENTERB1,XSHIFTB1,YSHIFTB1,XERRB1,YERRB1,CIERB1,CERRORB1,sl2B1,$
   MSKYB1,STDEVB1,SSKEWB1,NSKYB1,NSREJB1,SIERB1,SERRORB1,sl3B1,$
   ITIMEB1,XAIRMASSB1,IFILTERB1,OTIMEB1,sl4B1,$
   RAPERTB1,SUMB1,AREAB1,FLUXB1,MAGB1,MERRB1,PIERB1,PERRORB1,$
  f='a,f,f,i,a,i,a,f,f,f,f,a,a,i,a,a,f,f,f,i,i,i,a,a,f,f,a,a,a,f,f,f,f,f,f,i,a'

  readcol,'b2.dat',$
   IMAGEB2,XINITB2,YINITB2,IDB2,COORDSB2,LIDB2,sl1B2,$
   XCENTERB2,YCENTERB2,XSHIFTB2,YSHIFTB2,XERRB2,YERRB2,CIERB2,CERRORB2,sl2B2,$
   MSKYB2,STDEVB2,SSKEWB2,NSKYB2,NSREJB2,SIERB2,SERRORB2,sl3B2,$
   ITIMEB2,XAIRMASSB2,IFILTERB2,OTIMEB2,sl4B2,$
   RAPERTB2,SUMB2,AREAB2,FLUXB2,MAGB2,MERRB2,PIERB2,PERRORB2,$
  f='a,f,f,i,a,i,a,f,f,f,f,a,a,i,a,a,f,f,f,i,i,i,a,a,f,f,a,a,a,f,f,f,f,f,f,i,a'

  readcol,'r1.dat',$
   IMAGER1,XINITR1,YINITR1,IDR1,COORDSR1,LIDR1,sl1R1,$
   XCENTERR1,YCENTERR1,XSHIFTR1,YSHIFTR1,XERRR1,YERRR1,CIERR1,CERRORR1,sl2R1,$
   MSKYR1,STDEVR1,SSKEWR1,NSKYR1,NSREJR1,SIERR1,SERRORR1,sl3R1,$
   ITIMER1,XAIRMASSR1,IFILTERR1,OTIMER1,sl4R1,$
   RAPERTR1,SUMR1,AREAR1,FLUXR1,MAGR1,MERRR1,PIERR1,PERRORR1,$
  f='a,f,f,i,a,i,a,f,f,f,f,a,a,i,a,a,f,f,f,i,i,i,a,a,f,f,a,a,a,f,f,f,f,f,f,i,a'

  readcol,'r2.dat',$
   IMAGER2,XINITR2,YINITR2,IDR2,COORDSR2,LIDR2,sl1R2,$
   XCENTERR2,YCENTERR2,XSHIFTR2,YSHIFTR2,XERRR2,YERRR2,CIERR2,CERRORR2,sl2R2,$
   MSKYR2,STDEVR2,SSKEWR2,NSKYR2,NSREJR2,SIERR2,SERRORR2,sl3R2,$
   ITIMER2,XAIRMASSR2,IFILTERR2,OTIMER2,sl4R2,$
   RAPERTR2,SUMR2,AREAR2,FLUXR2,MAGR2,MERRR2,PIERR2,PERRORR2,$
  f='a,f,f,i,a,i,a,f,f,f,f,a,a,i,a,a,f,f,f,i,i,i,a,a,f,f,a,a,a,f,f,f,f,f,f,i,a'

  readcol,'v1.dat',$
   IMAGEV1,XINITV1,YINITV1,IDV1,COORDSV1,LIDV1,sl1V1,$
   XCENTERV1,YCENTERV1,XSHIFTV1,YSHIFTV1,XERRV1,YERRV1,CIERV1,CERRORV1,sl2V1,$
   MSKYV1,STDEVV1,SSKEWV1,NSKYV1,NSREJV1,SIERV1,SERRORV1,sl3V1,$
   ITIMEV1,XAIRMASSV1,IFILTERV1,OTIMEV1,sl4V1,$
   RAPERTV1,SUMV1,AREAV1,FLUXV1,MAGV1,MERRV1,PIERV1,PERRORV1,$
  f='a,f,f,i,a,i,a,f,f,f,f,a,a,i,a,a,f,f,f,i,i,i,a,a,f,f,a,a,a,f,f,f,f,f,f,i,a'

  readcol,'v2.dat',$
   IMAGEV2,XINITV2,YINITV2,IDV2,COORDSV2,LIDV2,sl1V2,$
   XCENTERV2,YCENTERV2,XSHIFTV2,YSHIFTV2,XERRV2,YERRV2,CIERV2,CERRORV2,sl2V2,$
   MSKYV2,STDEVV2,SSKEWV2,NSKYV2,NSREJV2,SIERV2,SERRORV2,sl3V2,$
   ITIMEV2,XAIRMASSV2,IFILTERV2,OTIMEV2,sl4V2,$
   RAPERTV2,SUMGV2,AREAV2,FLUXV2,MAGV2,MERRV2,PIERV2,PERRORV2,$
  f='a,f,f,i,a,i,a,f,f,f,f,a,a,i,a,a,f,f,f,i,i,i,a,a,f,f,a,a,a,f,f,f,f,f,f,i,a'

print, "     Opening output file..."

OPENU, lun, 'starmags.txt', /Get_lun

tol = 20
starlimit = n_elements(IMAGEB1) - 1

print, "     Finding matches..."

counter = 0
;B1
  FOR i = 0, starlimit DO BEGIN
;B2
    FOR j = 0, starlimit DO BEGIN
;R1
       FOR k = 0, starlimit DO BEGIN
;R2          
          FOR l = 0, starlimit DO BEGIN
;V1
             FOR m = 0, starlimit DO BEGIN
;V2
                FOR n = 0, starlimit DO BEGIN

      IF (ABS(XCENTERB1[i] - XCENTERB2[j]) LE tol) AND (YCENTERB1[i] - YCENTERB2[j] LE tol) AND $
         (ABS(XCENTERB1[i] - XCENTERR1[k]) LE tol) AND (YCENTERB1[i] - YCENTERR1[k] LE tol) AND $
         (ABS(XCENTERB1[i] - XCENTERR2[l]) LE tol) AND (YCENTERB1[i] - YCENTERR2[l] LE tol) AND $
         (ABS(XCENTERB1[i] - XCENTERV1[m]) LE tol) AND (YCENTERB1[i] - YCENTERV1[m] LE tol) AND $
         (ABS(XCENTERB1[i] - XCENTERV2[n]) LE tol) AND (YCENTERB1[i] - YCENTERV2[n] LE tol) THEN BEGIN

         B = [MAGB1[i],MAGB2[j]]
         R = [MAGR1[k],MAGR2[l]]
         V = [MAGV1[m],MAGV2[n]]
         counter = counter + 1
 
         printf, lun, XCENTERB1[i], YCENTERB1[i], MEAN(B), MEAN(R), MEAN(V),$
                 f='(f9.4,f12.4,f10.4,f10.4,f10.4)'

      ENDIF

            ENDFOR

          ENDFOR

        ENDFOR

      ENDFOR

    ENDFOR
    
 ENDFOR

print, "     Found", counter, " stars", f='(a10,i2,a)' 
print, ""
print, "     NOW, INPUT 'realmags.txt' WITH *CORRESPONDING* REAL MAGNITUDES"
print, "              FOR COMPARISON IN ORDER: V, B-V, V-R"

CLOSE, lun
free_lun, lun

END

PRO zeropoint

  readcol, 'starmags.txt',$
;***put filters here
   X,Y,measureB,measureR,measureV,$
  f='f,f,f,f,f'

  readcol, 'realmags.txt',$
   realV,realB_V,realV_R,$
  f='f,f,f'

;***INPUT MAGNITUDE TRANSFORMATIONS HERE

;realG_R = 1.124*(realB_V) - 0.252
;realG = realV + 0.634*(realB_V) - 0.108
;realR = realG - realG_R

realR = realV - realV_R
realB = realV + realB_V

;
;

;ZP = real - measured + 25

ZP_B = fltarr(2)
ZP_R = fltarr(2)
ZP_V = fltarr(2)

  FOR i = 0, 1 DO BEGIN

    ZP_B[i] = realB[i] - measureB[i] + 25

  ENDFOR

  FOR i = 0, 1 DO BEGIN

    ZP_R[i] = realR[i] - measureR[i] + 25

 ENDFOR

  FOR i = 0, 1 DO BEGIN

    ZP_V[i] = realV[i] - measureV[i] + 25

  ENDFOR

OPENW, lun, 'zeropoints', /Get_lun

  printf, lun, "ZP_B = ", MEAN(ZP_B);, f='(f9.4)'
  printf, lun, "ZP_R = ", MEAN(ZP_R);, f='(f9.4)'
  printf, lun, "ZP_V = ", MEAN(ZP_V);, f='(f9.4)'

CLOSE, lun
free_lun, lun

END
