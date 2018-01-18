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