pro Read_vector,dir,file,numframe,xy,xz

  common parameters

  ;Open file-----------------------------------------
  close 1
  fullpath = dir+file+'.dat'
  openr,1,fullpath,/f77_unformatted
  ;--------------------------------------------------

  ;create array--------------------------------------
  frame = 0l
  xy=fltarr(nx,ny,3,/nozero)
  xz=fltarr(nx,ny,3,/nozero)
  framecount = 1
  ;--------------------------------------------------

  ;Populate array------------------------------------
  while (framecount le numframe) do begin
    frame = 0l
    readu,1,frame
    readu,1,xz,xy
    framecount = framecount + 1
  endwhile
  ;--------------------------------------------------

  print 'image # ', frame

  close,1

  return

end