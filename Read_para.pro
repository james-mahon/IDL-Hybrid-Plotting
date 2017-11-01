pro Read_para
;Define common parameters--------------------------
COMMON parameters,$
  nx,ny,nz,dx,dy,delz,$
  nt,dt,nout,$
  Ni_max,$
  m_pu,$
  alpha,beta,$
  b0_init,$
  nf_init,$
  vsw,vth
;--------------------------------------------------

;Set file path-------------------------------------
dir_root='/Users/jama3001/Data/' ;Base folder
filenom= '2017-Tue-Oct-24/pluto-1/'  ;Date and run number
dir=dir_root + filenom
full_path=dir +'databig/para.dat' ;Parameter file
;--------------------------------------------------

;Define variable types-----------------------------
nx=0l & ny=0l &  nz=0l
comm_sz=0l & nproc_up=0l  & nproc_down= 0l
dx=0. & dy=0. & delZ=0. & BIG_delz=0.
nt=0l & nout=0l & dt=0.
dtsub_init=0. & ntsub=0l
out_dir='1234567890'
vtop=0. & vbottom=0.
Ni_max=0l
mproton=0.0d & m_pu=0.0d & m_heavy=0.
np_top=0. & np_bottom=0.
b0_top=0. & b0_bottom=0.
vth_top=0. & vth_bottom=0.
alpha=0.0d & beta=0.
RIo=0.
ion_amu=0.
m_pu=0.0d
nf_init=0.
dt_frac=0.
vsw=0.
vth=0.
Ni_tot_frac=0.
dx_frac=0.
nu_init_frac=0.
mrestart=0l
ri0=0.0
part_nout=0l
;--------------------------------------------------

;Read Variables------------------------------------
close,1
openr,1,full_path,/f77_unformatted
      readu,1,nx,ny,nz,dx,dy,delz
      readu,1,nt,dtsub_init,ntsub,dt,nout
      readu,1,out_dir
      readu,1,vtop,vbottom
      readu,1,Ni_max
      readu,1,mproton,m_pu,m_heavy
      readu,1,np_top,np_bottom
      readu,1,b0_top,b0_bottom
      readu,1,vth_top,vth_bottom
      readu,1,alpha,beta
      readu,1,RIo
      readu,1,b0_init
      readu,1,ion_amu
      readu,1,m_pu
      readu,1,nf_init
      readu,1,dt_frac
      readu,1,vsw
      readu,1,vth
      readu,1,Ni_tot_frac
      readu,1,dx_frac
      readu,1,nu_init_frac
      readu,1,mrestart
      readu,1,ri0
      readu,1,part_nout
close,1
;--------------------------------------------------
end