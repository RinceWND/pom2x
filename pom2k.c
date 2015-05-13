C***********************************************************************
C
C     Common blocks for pom2k.f
C
C***********************************************************************
C
C     source_c should agree with source in pom2k.f.
C
      character*40 source_c
      parameter(source_c='pom2k  2006-05-03')
C
C***********************************************************************
C
C     Array sizes:
C
      integer
     $  im             ,imm1           ,imm2           ,jm             ,
     $  jmm1           ,jmm2           ,kb             ,kbm1           ,
     $  kbm2 
C
C***********************************************************************
C
C     Set size of problem here:
C
C     parameter
C -- file2ic (iproblem=3)
C    $ (im=41          ,jm=61          ,kb=16)
C -- seamount (iproblem=1)
C    $ (im=65          ,jm=49          ,kb=21)
C ---- from parameter file generated by runpom2k ----
      include 'grid'
C
C***********************************************************************
C                                     
      parameter
     $ (imm1=im-1      ,imm2=im-2      ,jmm1=jm-1      ,jmm2=jm-2      ,
     $  kbm1=kb-1      ,kbm2=kb-2      )
C
C-----------------------------------------------------------------------
C
C     Scalars:
C
      real
     $  alpha          ,dte            ,dti            ,dti2           ,  
     $  grav           ,hmax           ,kappa          ,pi             ,
     $  ramp           ,rfe            ,rfn            ,rfs            ,
     $  rfw            ,rhoref         ,sbias          ,slmax          ,
     $  small          ,tbias          ,time           ,tprni          ,
     $  umol           ,vmaxl          ,time0

      integer
     $  iint           ,iprint         ,iskp           ,jskp           ,
     $  kl1            ,kl2            ,mode           ,ntp
C
      common/blkcon/ 
     $  alpha          ,dte            ,dti            ,dti2           ,
     $  grav           ,hmax           ,kappa          ,pi             ,
     $  ramp           ,rfe            ,rfn            ,rfs            ,
     $  rfw            ,rhoref         ,sbias          ,slmax          ,
     $  small          ,tbias          ,time           ,tprni          ,
     $  umol           ,vmaxl          ,time0          ,
     $  iint           ,iprint         ,iskp           ,jskp           ,
     $  kl1            ,kl2            ,mode           ,ntp
C
C-----------------------------------------------------------------------
C
C     1-D arrays:
C
      real
     $  dz             ,dzz            ,z              ,zz
C
      common/blk1d/ 
     $  dz(kb)         ,dzz(kb)        ,z(kb)          ,zz(kb) 
C
C-----------------------------------------------------------------------
C
C     2-D arrays:
C
      real
     $  aam2d          ,advua          ,advva          ,adx2d          ,
     $  ady2d          ,art            ,aru            ,arv            ,
     $  cbc            ,cor            ,d              ,drx2d          ,
     $  dry2d          ,dt             ,dum            ,dvm            ,
     $  dx             ,dy             ,east_c         ,east_e         ,
     $  east_u         ,east_v         ,e_atmos        ,egb            ,
     $  egf            ,el             ,elb            ,elf            ,
     $  et             ,etb            ,etf            ,fluxua         ,
     $  fluxva         ,fsm            ,h              ,north_c        ,
     $  north_e        ,north_u        ,north_v        ,psi            ,
     $  rot            ,ssurf          ,swrad          ,vfluxb         ,
     $  tps            ,tsurf          ,ua             ,vfluxf         ,
     $  uab            ,uaf            ,utb            ,utf            ,
     $  va             ,vab            ,vaf            ,
     $  vtb            ,vtf            ,wssurf         ,wtsurf         ,
     $  wubot          ,wusurf         ,wvbot          ,wvsurf         ,
! Arrays for interpolation
     $  swradb         ,tsurfb         ,ssurfb         ,wssurfb        ,
     $  wtsurfb        ,wusurfb        ,wvsurfb        ,
     $  swradf         ,tsurff         ,ssurff         ,wssurff        ,
     $  wtsurff        ,wusurff        ,wvsurff        ,
     $  qqb            ,dqb            ,qqf            ,dqf
C
      common/blk2d/
     $  aam2d(im,jm)   ,advua(im,jm)   ,advva(im,jm)   ,adx2d(im,jm)   ,
     $  ady2d(im,jm)   ,art(im,jm)     ,aru(im,jm)     ,arv(im,jm)     ,
     $  cbc(im,jm)     ,cor(im,jm)     ,d(im,jm)       ,drx2d(im,jm)   ,
     $  dry2d(im,jm)   ,dt(im,jm)      ,dum(im,jm)     ,dvm(im,jm)     ,
     $  dx(im,jm)      ,dy(im,jm)      ,east_c(im,jm)  ,east_e(im,jm)  ,
     $  east_u(im,jm)  ,east_v(im,jm)  ,e_atmos(im,jm) ,egb(im,jm)     ,
     $  egf(im,jm)     ,el(im,jm)      ,elb(im,jm)     ,elf(im,jm)     ,
     $  et(im,jm)      ,etb(im,jm)     ,etf(im,jm)     ,fluxua(im,jm)  ,
     $  fluxva(im,jm)  ,fsm(im,jm)     ,h(im,jm)       ,north_c(im,jm) ,
     $  north_e(im,jm) ,north_u(im,jm) ,north_v(im,jm) ,psi(im,jm)     ,
     $  rot(im,jm)     ,ssurf(im,jm)   ,swrad(im,jm)   ,vfluxb(im,jm)  ,
     $  tps(im,jm)     ,tsurf(im,jm)   ,ua(im,jm)      ,vfluxf(im,jm)  ,
     $  uab(im,jm)     ,uaf(im,jm)     ,utb(im,jm)     ,utf(im,jm)     ,
     $  va(im,jm)      ,vab(im,jm)     ,vaf(im,jm)     ,
     $  vtb(im,jm)     ,vtf(im,jm)     ,wssurf(im,jm)  ,wtsurf(im,jm)  ,
     $  wubot(im,jm)   ,wusurf(im,jm)  ,wvbot(im,jm)   ,wvsurf(im,jm)  ,
! Arrays for interpolation
     $  swradb(im,jm)  ,tsurfb(im,jm)  ,ssurfb(im,jm)  ,wssurfb(im,jm) ,
     $  wtsurfb(im,jm) ,wusurfb(im,jm) ,wvsurfb(im,jm) ,
     $  swradf(im,jm)  ,tsurff(im,jm)  ,ssurff(im,jm)  ,wssurff(im,jm) ,
     $  wtsurff(im,jm) ,wusurff(im,jm) ,wvsurff(im,jm) ,
     $  qqb(im,jm)     ,dqb(im,jm)     ,qqf(im,jm)     ,dqf(im,jm)
C
C-----------------------------------------------------------------------
C
C     3-D arrays:
C
      real 
     $  aam            ,advx           ,advy           ,a              ,
     $  c              ,drhox          ,drhoy          ,dtef           ,
     $  ee             ,gg             ,kh             ,km             ,
     $  kq             ,l              ,q2b            ,q2             ,
     $  q2lb           ,q2l            ,rho            ,rmean          ,
     $  sb             ,sclim          ,s              ,tb             ,
     $  tclim          ,t              ,ub             ,uf             ,
     $  u              ,vb             ,vf             ,v              ,
     $  w              ,zflux
C
      common/blk3d/
     $  aam(im,jm,kb)  ,advx(im,jm,kb) ,advy(im,jm,kb) ,a(im,jm,kb)    ,
     $  c(im,jm,kb)    ,drhox(im,jm,kb),drhoy(im,jm,kb),dtef(im,jm,kb) ,
     $  ee(im,jm,kb)   ,gg(im,jm,kb)   ,kh(im,jm,kb)   ,km(im,jm,kb)   ,
     $  kq(im,jm,kb)   ,l(im,jm,kb)    ,q2b(im,jm,kb)  ,q2(im,jm,kb)   ,
     $  q2lb(im,jm,kb) ,q2l(im,jm,kb)  ,rho(im,jm,kb)  ,rmean(im,jm,kb),
     $  sb(im,jm,kb)   ,sclim(im,jm,kb),s(im,jm,kb)    ,tb(im,jm,kb)   ,
     $  tclim(im,jm,kb),t(im,jm,kb)    ,ub(im,jm,kb)   ,uf(im,jm,kb)   ,
     $  u(im,jm,kb)    ,vb(im,jm,kb)   ,vf(im,jm,kb)   ,v(im,jm,kb)    ,
     $  w(im,jm,kb)    ,zflux(im,jm,kb)
C
C-----------------------------------------------------------------------
C
C     1 and 2-D boundary value arrays:
C
      real
     $  ele            ,eln            ,els            ,elw            ,
     $  sbe            ,sbn            ,sbs            ,sbw            ,
     $  tbe            ,tbn            ,tbs            ,tbw            ,
     $  uabe           ,uabw           ,ube            ,ubw            ,
     $  vabn           ,vabs           ,vbn            ,vbs            ,
     $  uabs           ,uabn           ,vabw           ,vabe           ,
     $  ubs            ,ubn            ,vbw            ,vbe            ,
! Arrays for interpolation
!  prev. month
     $  eleb           ,elnb           ,elsb           ,elwb           ,
     $  sbeb           ,sbnb           ,sbsb           ,sbwb           ,
     $  tbeb           ,tbnb           ,tbsb           ,tbwb           ,
     $  uabeb          ,uabwb          ,ubeb           ,ubwb           ,
     $  vabnb          ,vabsb          ,vbnb           ,vbsb           ,
     $  uabsb          ,uabnb          ,vabwb          ,vabeb          ,
     $  ubsb           ,ubnb           ,vbwb           ,vbeb           ,
!  next. month
     $  elef           ,elnf           ,elsf           ,elwf           ,
     $  sbef           ,sbnf           ,sbsf           ,sbwf           ,
     $  tbef           ,tbnf           ,tbsf           ,tbwf           ,
     $  uabef          ,uabwf          ,ubef           ,ubwf           ,
     $  vabnf          ,vabsf          ,vbnf           ,vbsf           ,
     $  uabsf          ,uabnf          ,vabwf          ,vabef          ,
     $  ubsf           ,ubnf           ,vbwf           ,vbef
C
      common/bdry/
     $  ele(jm)        ,eln(im)        ,els(im)        ,elw(jm)        ,
     $  sbe(jm,kb)     ,sbn(im,kb)     ,sbs(im,kb)     ,sbw(jm,kb)     ,
     $  tbe(jm,kb)     ,tbn(im,kb)     ,tbs(im,kb)     ,tbw(jm,kb)     ,
     $  uabe(jm)       ,uabw(jm)       ,ube(jm,kb)     ,ubw(jm,kb)     ,
     $  vabn(im)       ,vabs(im)       ,vbn(im,kb)     ,vbs(im,kb)     ,
     $  uabs(im)       ,uabn(im)       ,vabw(jm)       ,vabe(jm)       ,
     $  ubs(im,kb)     ,ubn(im,kb)     ,vbw(jm,kb)     ,vbe(jm,kb)     ,
!
     $  eleb(jm)       ,elnb(jm)       ,elsb(im)       ,elwb(jm)       ,
     $  sbeb(jm,kb)    ,sbnb(im,kb)    ,sbsb(im,kb)    ,sbwb(jm,kb)    ,
     $  tbeb(jm,kb)    ,tbnb(im,kb)    ,tbsb(im,kb)    ,tbwb(jm,kb)    ,
     $  uabeb(jm)      ,uabwb(jm)      ,ubeb(jm,kb)    ,ubwb(jm,kb)    ,
     $  vabnb(im)      ,vabsb(im)      ,vbnb(im,kb)    ,vbsb(im,kb)    ,
     $  uabsb(im)      ,uabnb(im)      ,vabwb(jm)      ,vabeb(jm)      ,
     $  ubsb(im,kb)    ,ubnb(im,kb)    ,vbwb(jm,kb)    ,vbeb(jm,kb)    ,
!
     $  elef(jm)       ,elnf(jm)       ,elsf(im)       ,elwf(jm)       ,
     $  sbef(jm,kb)    ,sbnf(im,kb)    ,sbsf(im,kb)    ,sbwf(jm,kb)    ,
     $  tbef(jm,kb)    ,tbnf(im,kb)    ,tbsf(im,kb)    ,tbwf(jm,kb)    ,
     $  uabef(jm)      ,uabwf(jm)      ,ubef(jm,kb)    ,ubwf(jm,kb)    ,
     $  vabnf(im)      ,vabsf(im)      ,vbnf(im,kb)    ,vbsf(im,kb)    ,
     $  uabsf(im)      ,uabnf(im)      ,vabwf(jm)      ,vabef(jm)      ,
     $  ubsf(im,kb)    ,ubnf(im,kb)    ,vbwf(jm,kb)    ,vbef(jm,kb)
C
C-----------------------------------------------------------------------
C
C     Character variables:
C
      character*26
     $  time_start
C
      character*40
     $  source,title
C
      common/blkchar/
     $  time_start     ,source         ,title
C
C-----------------------------------------------------------------------
!
!     NetCDF variables:
!
      integer ncid, varid, ncptime
      integer dim_time, dim_srho, dim_sw, dim_lat, dim_lon
      integer status
      character*128 filename
      integer dbg_lvl, dbg_seq_i, dbg_step
      real    dbg_off
C
C-----------------------------------------------------------------------
C     my path settings variables
      character*224 pth_wrk, pth_grd, pth_flx, pth_bry, pth_out, pth_bkp
     $ , ptf_rst, pth_dbg
      character*32 pfx_dmn, pfx_dbg
      common/pth/ pth_wrk, pth_grd, pth_flx, pth_bry, pth_out, pth_bkp
     $ , ptf_rst, pfx_dmn, pth_dbg, pfx_dbg
C
      integer rf_ts, rf_sts, rf_uv, rf_swrad, rf_wtsur, rf_wsurf, rf_el
      real tsplt_ts, tsplt_uv, tsplt_swrad, tsplt_wsurf, tsplt_wtsur
     $ , tsplt_el
      common/upd_cond/ rf_ts, rf_sts, rf_uv, rf_swrad, rf_wtsur
     $ , rf_wsurf, rf_el, tsplt_ts, tsplt_uv, tsplt_swrad, tsplt_wsurf
     $ , tsplt_wtsur, tsplt_el

      real fac  ! fac - interpolation factor
      integer m ! m   - interpolation month
      common/misc/ m, fac, dbg_lvl, dbg_off, dbg_seq_i, dbg_step
C-----------------------------------------------------------------------
C
C     End of common blocks
C
C-----------------------------------------------------------------------
C