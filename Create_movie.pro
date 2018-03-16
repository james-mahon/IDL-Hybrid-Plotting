pro Create_movie

close,1

;Forced Variables----------------------------------
;  Shift_R, shift_y
;  m_bkg
;--------------------------------------------------

;Plot settings-------------------------------------
device,decompose=0
!p.CHARSIZE=1.2
!p.multi=[0,1,1]
;--------------------------------------------------

;Questions-----------------------------------------
answer_diversion='  '
READ,'Compute diversion from many flowlines (y/n)?', answer_diversion
answer_sat = ' '
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
;
;Set directory-------------------------------------
dir = dir_root + filenom + 'databig/grid/'
;--------------------------------------------------

;Read variables------------------------------------

for nfrm = 1,nframe do begin
  print,'Frame: ',nfrm
  ;string(input,FORMAT='(I2.1)')
  Read_scalar, dir, 'c.np_3d_'+strtrim(string(proc_read,format='(I2.1)'),2),    nfrm,np_xy,np_xz ;bulk particle density
  Read_vector, dir, 'c.up_3d_'+strtrim(string(proc_read,format='(I2.1)'),2),    nfrm,up_xy,up_xz ;bulk particle velocity
  Read_scalar, dir, 'c.temp_p_3d_'+strtrim(string(proc_read,format='(I2.1)'),2),nfrm,tp_xy,tp_xz ;temp
  Read_vector, dir, 'c.b1_3d_'+strtrim(string(proc_read,format='(I2.1)'),2),    nfrm,b1_xy,b1_xz ;b field
;--------------------------------------------------

;Create velocity array-----------------------------
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
  
  ;loadct,0,/silent
  oplot,cx,cy,psym = 1,color = 255,symsize = 0.5
  
  ;loadct,39,/silent
  ;mycolorbar. range = [min(im_clip),max(im_clip)],$
    ;position=[!y.window(0),!x.window(1)+0.01, !y.window(1), $
    ;          !x.window(1)+0.03],$
    ;ncolors = 254,/vertical,/right,$
    ;title='density (cc)',$
    ;charsize=1.4,yaxis=1,minor=0)
   
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
  
  m_bkg = 15.975128
  mo = m_bkg * 1.6726e-27
  
  b2dx = b2dx * Mo /1.6e-19 / 1e-9
  b2dy = b2dy * Mo /1.6e-19 / 1e-9
  b2dz = b2dz * Mo /1.6e-19 / 1e-9
  
  b12d_tot = b2dz
  
  ;if (answer_sat eq 'y') then begin
  ;  WB = where(b12d_tot ge 401)
  ;  b12d_tot(WB) = 401.
  ;endif
  
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
  print,'   center Cx/y_center(Reur)=',CX_center/Rio, Cy_center/Rio
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
  