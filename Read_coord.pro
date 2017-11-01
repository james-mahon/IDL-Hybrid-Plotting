pro Read_coord
;Set file path-------------------------------------
dir_root='/Users/jama3001/Data/' ;Base folder (User specific)
filenom= '2017-Wed-Nov-01/pluto-5/'  ;Date and run number (Dynamic)
dir=dir_root + filenom
full_path=dir +'databig/c.coord.dat' ;Coordinate file (Static)
;--------------------------------------------------

;Read domain dimensions-----------------------------
nx=0l
ny=0l
nz=0l

openr,11,full_path,/f77_unformatted
readu,11,nx
readu,11,ny
readu,11,nz
;--------------------------------------------------

;Construct grid------------------------------------
x=fltarr(nx,/nozero)
y=fltarr(ny,/nozero)
z=fltarr(nz,/nozero)
dz_grid=fltarr(nz,/nozero)
dz_cell=fltarr(nz,/nozero)

readu,11,x
readu,11,y
readu,11,z
readu,11,dz_grid
readu,11,dz_cell
close,11
;--------------------------------------------------
end