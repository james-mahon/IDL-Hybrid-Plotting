pro img_cont, a, x, y,dx,dy,nfrm, WINDOW_SCALE = window_scale, $
  ASPECT = aspect, INTERP = interp,POSTSCRIPT=postscript
  ;+
  ; NAME:
  ; IMAGE_CONT
  ; PURPOSE:
  ; Overlay an image and a contour plot.
  ; CATEGORY:
  ; General graphics.
  ; CALLING SEQUENCE:
  ; IMAGE_CONT, A
  ; INPUTS:
  ; A = 2 dimensional array to display.
  ; KEYWORD PARAMETERS:
  ; /WINDOW_SCALE = set to scale the window size to the image size,
  ;   otherwise the image size is scaled to the window size.
  ;   Ignored when outputting to devices with scalable pixels.
  ; /ASPECT = set to retain image's aspect ratio.  Assumes square
  ;   pixels.  If /WINDOW_SCALE is set, the aspect ratio is
  ;   retained.
  ; /INTERP = set to bi-linear interpolate if image is resampled.
  ; OUTPUTS:
  ; No explicit outputs.
  ; COMMON BLOCKS:
  ; none.
  ; SIDE EFFECTS:
  ; The currently selected display is affected.
  ; RESTRICTIONS:
  ; None that are obvious.
  ; PROCEDURE:
  ; If the device has scalable pixels then the image is written over
  ; the plot window.
  ; MODIFICATION HISTORY:
  ; DMS, May, 1988.
  ;       DOLS MAy 2011, add dx and dy in call to label axis cosistently
  ;-

  nclr = !d.n_colors

  if keyword_set(postscript) then begin
    set_plot,'ps
    device,filename='img.eps
    device,/encapsulated
    !p.font=0
    device,/palatino
    device,bits=8
    device,/color
  endif

  ;f_read_coord,'coord.dat',x,y,z,dzg,dzc,nx,ny,nz
  ;ax = y
  ;az = z

  ;a = bytscl(a)
  a = a<254
  maxa=max(a)
  mina=min(a)
  ;help,a,maxa,mina
  ;print,!d.n_colors,nclr
  ;a = (not(bytscl(a)))*nclr/256
  ;print,max(a)
  ;stop
  sz=size(a)
  ;ab = bytscl(a)
  ;a=rebin(a,4*sz(1),4*sz(2))
  ;sz=size(a)


  ;!x.title = 'y (10!u3!n km)'
  ;!y.title = 'z (10!u3!n km)'
  ;!x.title = 'time (s)'
  ;!y.title = 'mode'

  !p.charsize=1.2
  !y.margin = [4,4]
  !x.margin = [10,10]
  ;clrb = findgen(sz(2))*nclr/max(findgen(sz(2)))
  ;;clrb = findgen(sz(2))*255/max(findgen(sz(2)))
  ;for i=sz(1)-15,sz(1)-1 do begin
  ;   a(i,*) = max(clrb)-clrb
  ;;   a(i,*) = clrb
  ;endfor

  on_error,2                      ;Return to caller if an error occurs
  sz = size(a)      ;Size of image
  if sz(0) lt 2 then message, 'Parameter not 2D'

  ;set window used by contour

  ;DOLS label of axis
  ;x1=strtrim(dx,1) & c1=strmid(x1,0,6)
  x1=strtrim(dx,1) & x1=strmid(x1,0,4)
  y1=strtrim(dy,1) & y1=strmid(y1,0,4)

  contour,[[0,0],[1,1]],/nodata, xstyle=1, ystyle = 1,$
    title='HORIZONTAL CUT BUT NOT AT EQUAT PLANE!!!',$
    ;            ytitle='y (60 km)',xtitle='x (60 km)',$
    ytitle='y ('+y1+' km)',xtitle='x (' + x1 + ' km)',$
    ;xrange=[min(x),max(x)],yrange=[min(y),max(y)+1],$  ; dipslay in km or pxl
    xrange=[min(x),max(x)],yrange=[min(y),max(y)],$ ; display in RIo DOLS 1/11/13
    ;            title='t = '+strmid(strtrim(string(0.5*nfrm*100.),2),0,6) + ' (s)',$
    /isotropic
  ;            xrange=[min(x),800],yrange=[min(y),max(y)+1]
  p1 = !P & x1 = !X & y1= !Y

  px = !x.window * !d.x_vsize ;Get size of window in device units
  py = !y.window * !d.y_vsize
  swx = px(1)-px(0)   ;Size in x in device units
  swy = py(1)-py(0)   ;Size in Y
  six = float(sz(1))    ;Image sizes
  siy = float(sz(2))
  aspi = six / siy    ;Image aspect ratio
  aspw = swx / swy    ;Window aspect ratio
  f = aspi / aspw     ;Ratio of aspect ratios

  if (!d.flags and 1) ne 0 then begin ;Scalable pixels?
    if keyword_set(aspect) then begin ;Retain aspect ratio?
      ;Adjust window size
      if f ge 1.0 then swy = swy / f else swx = swx * f
    endif
    print,max(a)
    tv,a,px(0),py(0),xsize = swx, ysize = swy, /device


  endif else begin  ;Not scalable pixels
    if keyword_set(window_scale) then begin ;Scale window to image?
      tv,a,px(0),py(0)  ;Output image
      swx = six   ;Set window size from image
      swy = siy
    endif else begin    ;Scale window
      if keyword_set(aspect) then begin
        if f ge 1.0 then swy = swy / f else swx = swx * f
      endif   ;aspect
      print,max(a)
      tv,poly_2d((a),$  ;Have to resample image
        [[0,0],[six/swx,0]], [[0,siy/swy],[0,0]],$
        keyword_set(interp),swx,swy), $
        px(0),py(0)

    endelse     ;window_scale
  endelse     ;scalable pixels

  axis,xaxis=2,xstyle=1
  axis,yaxis=2,ystyle=1
  ;     axis,yaxis=1,ystyle=1,$
  ;;     yticks=2,$
  ;;     ytickv=[min(y),(min(y)+max(y))/2,max(y)], $
  ;;     ytickname=[string(min(a)),string((min(a)+max(a))/2),string(max(a))], $
  ;;     ytitle='E!dz!n (mV/m)',$
  ;;     yrange = [mina*(2.3e-25/1.6e-19)*1e6,maxa*(2.3e-25/1.6e-19*1e6)]
  ;     ytitle='(u!de!n)!dz!n (km/s)',$
  ;     yrange = [mina,maxa]


  mx = !d.n_colors-1    ;Brightest color
  colors = [mx,mx,mx,0,0,0] ;color vectors
  if !d.name eq 'PS' then colors = mx - colors ;invert line colors for pstscrp
  m = max(a)
  ;contour,[[0,0],[1,1]],/nodata, xstyle=1, ystyle = 1
  ;contour,a,/noerase,/nodata,/yst,$  ;Do the contour
  ;   pos = [px(0),py(0), px(0)+swx,py(0)+swy],/dev,$
  ;          c_color =  colors, charsize=1.5, $
  ;          xtitle='x (km)', ytitle='y (km)', $
  ;          levels=[0.02*m,0.04*m,0.06*m,0.2*m,0.4*m,0.5*m,0.7*m,0.9*m]
  ;;    xrange=[min(ax),max(ax)], xstyle=1
  ;dx = !x.crange(1) - !x.crange(0)
  ;dy = !y.crange(1) - !y.crange(0)
  ;xyouts,0.05*dx+!x.crange(0),0.85*dy+!y.crange(0),tit,/data
  ;!p.multi=[0,1,2]

  !P = p1 & !X = x1 & !Y = y1

  ;contour,[[0,0],[1,1]],/nodata, xstyle=1, ystyle = 1,/noerase

  if keyword_set(postscript) then begin
    device,/close
    set_plot,'x'
    !p.font=-1
  endif

  return
end

;************************************************************************
;************************************************************************
;            MAIN CODE STARTS HERE
;************************************************************************
;************************************************************************




pro Create_movie

close,1

;Plot settings-------------------------------------
device,decompose=0
!p.CHARSIZE=1.2
!p.multi=[0,1,1]
;--------------------------------------------------

;Reading parameters and coordinates----------------
common parameters
print,'First need to run Read_para.pro'
read,'How many processors?',num_proc

Read_coord,nx,ny,nz,x,y,z,dz_grid,dz_cell
;--------------------------------------------------

;Compute center processor--------------------------
dz = dz_cell(0)
nnz = (nz-1) + (num_proc-1)*(nz-2)+1 ;number of grid points on stitched grid
gz = findgen(nnz)*dz
cz = gz(nnz-1)/2.
print,'Center of the vertical grid in km ',cz

middle_proc_bot = num_proc/2+1 ;central processor counted from below
proc_read=1+num_proc-middle_proc_bot;output file number
gz_proc_base=z(nz-1) * (middle_proc_bot-1) + z(1) - (middle_proc_bot-1) * dz

zslc = 0
;--------------------------------------------------

;Prepare window for animation----------------------
xsz = 1600
ysz = 1000
nframe = 100
XINTERANIMATE, SET=[xsz,ysz, nframe], /SHOWLOAD
;--------------------------------------------------

;Set directory-------------------------------------
dir = dir_root + filenom + 'databig/grid/'
;--------------------------------------------------

;Read variables------------------------------------
nt = 10 ;TODO: Temp values, need to update

for nfrm = 1,nframe do begin
  print,'Frame: ',nfrm
  ;string(input,FORMAT='(I2.1)')
  Read_scalar, dir, 'c.np_3d_'+strtrim(string(proc_read,format='(I2.1)'),2),    nfrm,np_xy,np_xz ;bulk particle density
  Read_vector, dir, 'c.up_3d_'+strtrim(string(proc_read,format='(I2.1)'),2),    nfrm,up_xy,up_xz ;bulk particle velocity
  Read_scalar, dir, 'c.temp_p_3d_'+strtrim(string(proc_read,format='(I2.1)'),2),nfrm,tp_xy,tp_xz ;temp
  Read_vector, dir, 'c.b1_3d_'+strtrim(string(proc_read,format='(I2.1)'),2),    nfrm,b1_xy,b1_xz ;b field
;--------------------------------------------------

;TODO: Not sure------------------------------------
  uf=up_xy
  uf2d = reform(uf(*,*,*))
  uf2d2 = sqrt(uf2d(*,*,0)^2 + uf2d(*,*,1)^2 + uf2d(*,*,2)) ;total 3d speed in the xy plane
;--------------------------------------------------

;Redefine variables--------------------------------
  pf2d = reform(np_xy/1e15 * tp_xy * 1.6e-19)
  nf2d = reform(np_xy(*,*))
  b12dxy = reform(b1_xy(*,*,*)) ;vector component of bfield in the xy plane
  ntot2d = np_xy
  Ti2d = tp_xy
  uf2d = up_xy
;--------------------------------------------------

;Axis Reversal-------------------------------------
  xmid = (x(0) + x(nx-1))/2
  ymid = (y(0) + y(ny-1))/2
  ;central point
  
  y1 = reverse(y)
  x1 = reverse(x)
  ;so flow comes from the left and jupiter is at top of image
  
  x = xmid-x1
  y = ymid-y1
  gz = gz(zslc) - gz
  ;center axis on Europa
;--------------------------------------------------

;Disc plotting-------------------------------------
  ;TODO: shift_R, shift_y is forced
  shift_R = 1
  shift_y = 1
  CX_center = max(x) - (shift_R * Rio)
  CY_center = shift_y * max(y)
  print,'************************************************************'
  print,'DISC LOCATION'
  print,'*************'
  print,'SHIFT OF DISC SHIFT_R from outflow boundary(Rio)=',SHIFT_R
  print,'CX_center= max(X)-(shift_R*Rio): max(X)(Rio)=',MAX(X)/Rio
  print,'CENTER OF DISC XY from box center(km/Rio)=',CX_center,CY_center,CX_center/Rio,CY_center/Rio
  print,'minmax (x)(km)=',min(x),max(x),' Rio=',min(x)/Rio ,max(x)/Rio
  print,'minmax (Y)(km)=',min(Y),max(Y),' Rio=',min(Y)/Rio,max(Y)/Rio
  print,'************************************************************'
    theta = findgen(360)
    cx = CX_center/Rio + cos(theta*!dtor)
    cy = CY_center/Rio + sin(theta*!dtor)
;--------------------------------------------------

;Vector Reversal-----------------------------------
  uf2dx = rotate(-1. *b12dxy(*,*,0),2)
  uf2dy = rotate(-1. *b12dxy(*,*,1),2)
  uf2dz = rotate(1. *b12dxy(*,*,2),2)
  
  b2dx = rotate(-1.* b12dxy(*,*,0),2)
  b2dy = rotate(-1.* b12dxy(*,*,1),2)
  b2dz = rotate(1.* b12dxy(*,*,2),2)  
  ntot2d = rotate(ntot2d(*,*),2)
  pf2d = rotate(pf2d(*,*),2)
  ti2d = rotate(ti2d(*,*),2)
  
  uf3dx = fltarr(nx,ny)
  uf3dy = fltarr(nx,ny)
  uf3dz = fltarr(nx,ny)
;--------------------------------------------------

;Create movie with 4 windows-----------------------
  !p.multi=[0,2,2]
  loadct,39,/silent
  im_clip=ntot2d/1e15
  img_cont,bytscl(im_clip),x/Rio,y/Rio,dx,dy,nfrm
  
  loadct,0,/silent
  oplot,cx,cy,psym = 1,color = 255,symsize = 0.5
  
  loadct,39,/silent
  ;TODO: Investigate why color bars aren't working
  ;my_colorbar, range = [min(im_clip),max(im_clip)],$
  ;  position=[!y.window(0),!x.window(1)+0.01, !y.window(1), $
  ;            !x.window(1)+0.03],$
  ;  ncolors = 254,/vertical,/right,$
  ;  title='density (cc)',$
  ;  charsize=1.4,yaxis=1,minor=0
    
  img_cont,bytscl(ti2d),x/RIo,y/RIo,dx,dy,nfrm
  
  loadct,0,/silent
  oplot,cx,cy,psym=1,thick=1,color=255,symsize=0.5
  
  ;my_colorbar, range=[min(ti2d),max(ti2d)],$
  ;  position=[!y.window(0), !x.window(1)+0.01, !y.window(1), $
  ;                 !x.window(1)+0.03],$
  ;  ncolors=254,/vertical,/right,$
  ;  title = 'ion temperature (eV)',$
  ;  charsize=1.4,yaxis=1,minor=0
    
  JUMP_PASS1:
;--------------------------------------------------

;Plotting x component of flow velocity-------------
  loadct,39,/silent
  
  print,'  PLOT VELCOCITY but  check max/min Uf2dx=',max(uf2dx),min(uf2dx)
  
  ;TODO: commented out if block.  Function?
  ;if (answer_v_sat eq 'y') then begin
  ;  print,'PRESCRIBED saturation of velocity at Vsat=',vsat
  ;  print,' but  nfrm  check max/min Uf2dx=',nfrm , max(uf2dx),min(uf2dx)
  ;  if(nfrm eq nframe) then begin   ; for the last plot
  ;     WT= where (uf2dx lt vsat)
  ;      uf2dx(wT)=0.
  ;     endif
  ;  endif
    
  img_cont,bytscl(uf2dx),x/RIo,y/RIo,dx,dy,nfrm
  
  loadct,0,/silent
  
  oplot,cx,cy,psym=1,thick=1,color=255,symsize=0.5
  
  loadct,39,/silent
  ;my_colorbar,position=[!y.window(0), !x.window(1)+0.01, !y.window(1), $
  ;                    !x.window(1)+0.03],ncolors=255,/vertical,/right,$
  ;          minrange=[min(abs(-uf2dx))],$
  ;          maxrange=[max(abs(-uf2dx))],$
  ;          divisions=4,title='flow velocity VX (km/s)',$
  ;          charsize=1.4
;--------------------------------------------------

;Plotting Flowlines--------------------------------
  print,'Plotting flowlines'
  
  Ystart=0.
  ddy = 100.
  lmin = 0
  lmax = 30
  ddy = 24.
  
  if(answer_diversion eq 'y') then begin
    lmin=0  & lmax= 480
    lmin=0  & lmax= 680
    ddy= 1.  ;  for dx= 12km separation of fwl in km
  endif
  
  for l = lmin,lmax do begin
    xx = x(1)
    yy = Ystart+y(ny/2) + (-1*lmax+2*l)* ddy
    xx0 = xx
    yy0 = yy
    
    r = sqrt(xx^2 + yy^2)
    
    while (xx lt max(x)) do begin
      i = floor(xx/dx-x(0)/dx)-1
      ip = ceil(xx/dx-x(0)/dx)-1
      j = floor(yy/dy - y(0)/dy)-1
      jp = ceil(yy/dy - y(0)/dy)-1
      
      uxx = uf2dx(i,j)
      uyy = uf2dy(i,j)
      
      xx1 = xx + uxx * dt
      yy1 = yy + uyy * dt
      r1 = sqrt((xx1-CX_CENTER)^2 + (yy1-CY_CENTER)^2)
      
      IF(answer_diversion eq 'y') then begin
        IF (r1 le Rio) then begin ; to stop fwl that hit Io
          print,'FWL l hit IO ',l
          print,'     #fwl, Yres(km)=',lmax,dy
          print,'    XX0, YY0   =',XX0,YY0
          print,'    XX1 , YY1 ,r1 hit( km)=',XX1,YY1,r1
          print,'    XX1 , YY1 ,r1 hit(rio)=',XX1/Rio,YY1/Rio,r1/Rio
          GOTO,JUMP_OTHER_FWL
        endif
      ENDIF
      
      if(uxx le 1.) then begin ; there I am sure that the flow is very slow at 1km/s
        print,' WARNING, FWL#',l,'  UX IS <1km/s or NEGATIVE  reflection or on disc?'
        print,'              xx,yy, r (Rio)=',xx/Rio,yy/Rio,r/Rio
        print,'              uxx,uyy (km/s)=',uxx,uyy
        GOTO,JUMP_OTHER_FWL
      endif
      
      if((l eq 471) or (l eq 270)) then begin
        plot,[xx,xx1]/RIo,[yy,yy1]/Rio,/data,thick=2.0,linestyle=0,color=255
      endif
      
      xx = xx1
      yy = yy1
      r = r1
      
      endwhile
    JUMP_OTHER_FWL:
    
    endfor
  JUMP_PASS:
  
  mo = m_bkg * 1.6726e-27
  
  b2dx = b2dx * Mo /1.6e-19 / 1e-9
  b2dy = b2dy * Mo /1.6e-19 / 1e-9
  b2dz = b2dz * Mo /1.6e-19 / 1e-9
  
  b12d_tot = b2dz
  ;TODO: what is b12d_tot?
  
  if (answer_sat eq 'y') then begin
    WB = where(b12d_tot ge 401)
    b12d_tot(WB) = 401.
  endif
  
  img_cont,bytscl(b12d_tot),x/RIo,y/RIo,dx,dy,nfrm
  
  loadct,0,/silent
  
  oplot,cx,cy,psym=1,thick=1,color=255,symsize=0.5
  
  loadct,39,/silent
  ;my_colorbar,position=[!y.window(0), !x.window(1)+0.01, !y.window(1), $
  ;                      !x.window(1)+0.03],ncolors=255,/vertical,/right,$
  ;       minrange=[min(abs(b12d_tot))],$
  ;       maxrange=[max(abs(b12d_tot))],divisions=4,$
  ;       title='BZ (nT)',$
  ;       charsize=1.4
         
  JUMP_PASS_COLORBAR:
  JUMP_B:
  
  im = tvrd(true=1)
  tv,im
  
  xinteranimate, frame = nfrm-1, image = im
  
  print,'***********************************************************'
  print,'***********************************************************'
  print,'END OF TIMESTEP nfrm=',nfrm
  print,'     minmax(nf) (cm-3)= ',min(ntot2d)/1e15,max(ntot2d)/1.e15
  print,'     minmax(Ti) (eV)= ',min(Ti2d),max(Ti2d)
  print,'     minmax(ufx) (km/s)= ',min(uf2dx),max(uf2dx)
  print,'     minmax(bz) (nT)= ',min(b12d_tot),max(b12d_tot)
  print,'DISC SHIFTED shift_r(Rio)=?',shift_R
  print,'   center Cx/y_center(Reur)=',CX_center/Rio_km, Cy_center/Rio_km
  print,'***********************************************************'
  print,'  '
  print,'  '
  print,'  '
  print,'  '
endfor

xinteranimate,/keep_pixmaps

print,'WARNING contour plot might not be equatorail plane'
print,'----------------------------------------------------'
print,' total number of processors comm_sz=',num_proc
print,'   equatorial plane  global coord cz (km)=',cz
print,'   middle procfrom bottom & hybrid middle proc(from top) output to be read=',middle_proc_bot,proc_read
print,'   base (nz=2 fortran)  of middle processor (km)=',gz_proc_base

return
end
  