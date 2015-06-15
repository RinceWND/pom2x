      program pom2k
C
C **********************************************************************
C *                                                                    *
C *   The last code change as rcorded in pom2k.change was on           *
C *                                                                    *
C *                     2006-05-03                                     *
C *                                  (adding IC from file)             *
C *                                                                    *
C * FUNCTION    :  This is a version of the three dimensional, time    *
C *                dependent, primitive equation, ocean model          *
C *                developed by Alan Blumberg and George Mellor with   *
C *                subsequent contributions by Leo Oey, Steve Brenner  *
C *                and others. It is now called the Princeton Ocean    *
C *                Model. Two references are:                          *
C *                                                                    *
C *                Blumberg, A.F. and G.L. Mellor; Diagnostic and      *
C *                  prognostic numerical circulation studies of the   *
C *                  South Atlantic Bight, J. Geophys. Res. 88,        *
C *                  4579-4592, 1983.                                  *
C *                                                                    *
C *                Blumberg, A.F. and G.L. Mellor; A description of a  *
C *                  three-dimensional coastal ocean circulation model,*
C *                  Three-Dimensional Coastal Ocean Models, Coastal   *
C *                  and Estuarine Sciences, 4, N.S. Heaps, ed.,       *
C *                  American Geophysical Union, 1-16, 1987.           *
C *                                                                    *
C *                In subroutine profq the model makes use of the      *
C *                turbulence closure sub-model described in:          *
C *                                                                    *
C *                Mellor, G.L. and T. Yamada; Development of a        *
C *                  turbulence closure model for geophysical fluid    *
C *                  problems, Rev. Geophys. Space Phys., 20, No. 4,   *
C *                  851-875, 1982.                                    *
C *            (note recent profq that includes breaking waves)        *
C *                                                                    *
C *                A user's guide is available:                        *
C *                                                                    *
C *                Mellor, G.L.; User's guide for a three-dimensional, *
C *                  primitive equation, numerical ocean model.        *
C *                  Princeton University Report, 1998.                *
C *                                                                    *
C *                In October 2001, the source code underwent          *
C *                revision by John Hunter of the University of        *
C *                Tasmania. Major aspects of the revision were:       *
C *                                                                    *
C *                (1) The revision was based on pom98 updated to      *
C *                    12/9/2001.                                      *
C *                (2) Declaration of all variables.                   *
C *                (3) Rationalisation of the input of all constants.  *
C *                (4) Modifications to the "printer" output.          *
C *                (5) Output to a netCDF file.                        *
C *                (6) Inclusion of surface freshwater flux.           *
C *                (7) Inclusion of atmospheric pressure.              *
C *                (8) Inclusion of an additional problem to check (6) *
C *                    and (7), above.                                 *
C *                (9) Inclusion of option for Smolarkiewicz           *
C *                    advection scheme.                               *
C *                                                                    *
C *                This revised version is functionally almost         *
C *                equivalent to pom98. The output to device 6 from    *
C *                the "seamount" problem should be almost the same,   *
C *                any differences being due to minor format changes   *
C *                and improvements in rounding.                       *
C *                                                                    *
C *                This revision was helped by the following people:   *
C *                Tal Ezer, Peter Holloway, George Mellor, Rich       *
C *                Signell, Ian Webster, Brian Williams and Emma Young.*
C *                                                                    *
C **********************************************************************
C *                                                                    *
C *                                  GENERAL NOTES                     *
C *                                                                    *
C *                1. All units are S.I. (M.K.S.) unless otherwise     *
C *                   stated. NOTE that time is in days from the start *
C *                   of the run.                                      *
C *                                                                    *
C *                2. "b", <nothing> and "f" refers to backward,       *
C *                   central and forward time levels.                 *
C *                                                                    *
C *                3. NetCDF output may be used. In order to omit/use  *
C *                   netCDF, comment/uncomment all statements         *
C *                   carrying the comment "*netCDF*" at the end of    *
C *                   the line (or set netcdf_file='nonetcdf')         *
C *                                                                    *
C *                4. NetCDF is version 3. An attempt has been made to *
C *                   conform to the NetCDF Climate and Forecast (CF)  *
C *                   Metadata Conventions, but this may not yet be    *
C *                   complete (see:                                   *
C *                                                                    *
C *          http://www.cgd.ucar.edu/cms/eaton/cf-metadata/index.html) *
C *                                                                    *
C *                5. In order to use netCDF, the program should be    *
C *                   compiled with the appropriate library. For       *
C *                   example, if using g77, you may need to type:     *
C *                                                                    *
C *                     g77 -o pom2k pom2k.f /usr/lib/libnetcdf.a      *
C *                                                                    *
C *                   You should also have the "include" file of       *
C *                   netCDF subroutines (pom2k.n).                    *
C *                                                                    *
C *                6. In order to use netCDF, you may need to change   *
C *                   the name of the "include" file in the statement: *
C *                                                                    *
C *                     include '/usr/include/netcdf.inc'              *
C *                                                                    *
C *                   in subroutine write_netcdf                       *
C *                                                                    *
C **********************************************************************
C *                                                                    *
C *                                SOFTWARE LICENSING                  *
C *                                                                    *
C *                This program is free software; you can redistribute *
C *                it and/or modify it under the terms of the GNU      *
C *                General Public License as published by the Free     *
C *                Software Foundation, either Version 2 of the        *
C *                license, or (at your option) any later version.     *
C *                                                                    *
C *                This program is distributed in the hope that it     *
C *                will be useful, but without any warranty; without   *
C *                even the implied warranty of merchantability or     *
C *                fitness for a particular purpose. See the GNU       *
C *                General Public License for more details.            *
C *                                                                    *
C *                A copy of the GNU General Public License is         *
C *                available at http://www.gnu.org/copyleft/gpl.html   *
C *                or by writing to the Free Software Foundation, Inc.,*
C *                59 Temple Place - Suite 330, Boston, MA 02111, USA. *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
C     New declarations plus ispi,isp2i:
C
      real aam_init,atot
      real cbcmax,cbcmin,darea
      real days,dte2,dvol
      real eaver
      real horcon
      real ispi,isp2i
      real period,prtd1,prtd2
      real saver,smoth,sw,swtch
      real taver
      real vamax,vtot,tsalt
      real z0b
      real tatm,satm
      integer io(100),jo(100),ko(100)
      integer i,iend,iext,imax,ispadv,isplit,iswtch
      integer j,jmax
      integer k
      integer nadv,nbct,nbcs,npg,nitera,nread
      integer iproblem
      logical lramp
      character*120 netcdf_file

      integer bkp_gap           ! rwnd:
      real slice_s, slice_f     !     :
      type bnd                  !
        logical nth
        logical est
        logical sth
        logical wst
      end type
      type BC
        logical    wnd
        logical    lrd
        logical    srd
        logical    ssf
        logical    vap
        type (bnd) bnd
      end type                  ! rwnd:
      type (BC) bcf
C
C***********************************************************************
C
C     source should agree with source_c in pom2k.c and source_n in pom2k.n.
C
      source='pom2k  2006-05-03'
C
c     if(source.ne.source_c) then
c       write(6,7)
c   7   format(/'Incompatible versions of program and include files ',
c    $          '..... program terminated'/)
c       stop
c     endif
C
C***********************************************************************
C
      small=1.e-9           ! Small value
C
      pi=atan(1.e0)*4.e0    ! PI
C
C***********************************************************************
C
C     Input of filenames and constants:
C
C     NOTE that the array sizes im, jm and kb should be set in pom2k.c
C
C-----------------------------------------------------------------------
C
      title='Run 1                                   ' ! run's title
C
C-----------------------------------------------------------------------
C
      netcdf_file='pom2k.nc'  ! netCDF output file
c     netcdf_file='nonetcdf'  ! disable netCDF output
C
C-----------------------------------------------------------------------
C
C     Problem number:
C
C     iproblem      problem      initialisation
C                    type          subroutine
C
C         1        seamount       seamount
C
C         2        conservation   box
C                  box
C
C         3        IC from file   file2ic
C
      iproblem=1
C
C-----------------------------------------------------------------------
C
C       mode                     description
C
C        2        2-D calculation (bottom stress calculated in advave)
C
C        3        3-D calculation (bottom stress calculated in profu,v)
C
C        4        3-D calculation with t and s held fixed
C
      mode=3
C
C-----------------------------------------------------------------------
C
C     Advection scheme:
C
C      nadv     Advection scheme
C
C        1       Centred scheme, as originally provide in POM
C        2       Smolarkiewicz iterative upstream scheme, based on
C                subroutines provided by Gianmaria Sannino and Vincenzo
C                Artale
C
      nadv=1
C
C-----------------------------------------------------------------------
C
C     Constants for Smolarkiewicz iterative upstream scheme.
C
C     Number of iterations. This should be in the range 1 - 4. 1 is
C     standard upstream differencing; 3 adds 50% CPU time to POM:
C
      nitera=2
C
C     Smoothing parameter. This should preferably be 1, but 0 < sw < 1
C     gives smoother solutions with less overshoot when nitera > 1:
C
      sw=0.5e0
C
C-----------------------------------------------------------------------
C
C     Index to indicate whether run to start from restart file
C     (nread=0: no restart input file; nread=1: restart input file):
C
      nread=0
C
C-----------------------------------------------------------------------
C
C     External (2-D) time step (secs.) according to CFL:
C
      dte=6.e0
C
C-----------------------------------------------------------------------
C
C     <Internal (3-D) time step>/<External (2-D) time step>
C     (dti/dte; dimensionless):
C
      isplit=30
C
C-----------------------------------------------------------------------
C
C     Date and time of start of initial run of model in format (i.e.
C     UDUNITS convention)
C
C       YYYY-MM-DD HH:MM:SS <+/->HH:MM
C
C     where "<+/->HH:MM" is the time zone (positive eastwards from
C     Coordinated Universal Time). NOTE that the climatological time
C     axis (i.e. beginning of year zero, which does not exist in the
C     real-world calendar) has been used here. Insert your own date
C     and time as required:
C
      time_start='2000-01-01 00:00:00 +00:00'
C
C-----------------------------------------------------------------------
C
      days=0.25e0       ! run duration in days
C
C-----------------------------------------------------------------------
C
      prtd1=0.125e0     ! Initial print interval (days)
C
C-----------------------------------------------------------------------
C
      prtd2=1.e0         ! Final print interval (days)
C
C-----------------------------------------------------------------------
C
      swtch=1000.e0      ! Time to switch from prtd1 to prtd2
C
C-----------------------------------------------------------------------
C
      iskp=4             ! Printout skip interval in i
C
C-----------------------------------------------------------------------
C
      jskp=3             ! Printout skip interval in j
C
C-----------------------------------------------------------------------
C
C     Logical for inertial ramp (.true. if inertial ramp to be applied
C     to wind stress and baroclinic forcing, otherwise .false.)
C
      lramp=.false.
C
C-----------------------------------------------------------------------
C
C     Reference density (recommended values: 1025 for seawater,
C     1000 for freswater; S.I. units):
C
      rhoref=1025.e0
C
C-----------------------------------------------------------------------
C
      tbias=0.e0         ! Temperature bias (deg. C)
C
C-----------------------------------------------------------------------
C
      sbias=0.e0         ! Salinity bias
C
C-----------------------------------------------------------------------
C
      grav=9.806e0       ! gravity constant (S.I. units)
C
C-----------------------------------------------------------------------
C
      kappa=0.4e0        ! von Karman's constant
C
C-----------------------------------------------------------------------
C
      z0b=.01e0          ! Bottom roughness (metres)
C
C-----------------------------------------------------------------------
C
      cbcmin=.0025e0     ! Minimum bottom friction coeff.
C
C-----------------------------------------------------------------------
C
      cbcmax=1.e0        ! Maximum bottom friction coeff.
C
C-----------------------------------------------------------------------
C
      horcon=0.2e0       ! Smagorinsky diffusivity coeff.
C
C-----------------------------------------------------------------------
C
C     Inverse horizontal turbulent Prandtl number
C     (ah/am; dimensionless):
C
C     NOTE that tprni=0.e0 yields zero horizontal diffusivity!
C
      tprni=.2e0
C
C-----------------------------------------------------------------------
C
C     Background viscosity used in subroutines profq, proft, profu and
C     profv (S.I. units):
C
      umol=2.e-5
C
C-----------------------------------------------------------------------
C
C     Maximum depth used in radiation boundary condition in subroutine
C     bcond (metres):
C
      hmax=4500.e0
C
C-----------------------------------------------------------------------
C
C     Maximum magnitude of vaf (used in check that essentially tests
C     for CFL violation):
C
      vmaxl=100.e0
C
C-----------------------------------------------------------------------
C
C     Maximum allowable value of:
C
C       <difference of depths>/<sum of depths>
C
C     for two adjacent cells (dimensionless). This is used in subroutine
C     slpmax. If >= 1, then slpmax is not applied:
C
      slmax=2.e0
C
C-----------------------------------------------------------------------
C
C     Integers defining the number of logarithmic layers at the
C     surface and bottom (used by subroutine depth). The number of
C     logarithmic layers are kl1-2 at the surface and kb-kl2-1
C     at the bottom. For no log portions, set kl1=2 and kl2=kb-1:
C
      kl1=6
      kl2=kb-2
C
C-----------------------------------------------------------------------
C
C     Water type, used in subroutine proft.
C
C       ntp    Jerlov water type
C
C        1            i
C        2            ia
C        3            ib
C        4            ii
C        5            iii
C
      ntp=2
C
C-----------------------------------------------------------------------
C
C     Surface temperature boundary condition, used in subroutine proft:
C
C       nbct   prescribed    prescribed   short wave
C              temperature      flux      penetration
C
C        1        no           yes           no
C        2        no           yes           yes
C        3        yes          no            no
C        4        yes          no            yes
C
      nbct=1
C
C-----------------------------------------------------------------------
C
C     Surface salinity boundary condition, used in subroutine proft:
C
C       nbcs   prescribed    prescribed
C               salinity      flux
C
C        1        no           yes
C        3        yes          no
C
C     NOTE that only 1 and 3 are allowed for salinity.
C
      nbcs=1
C
!-----------------------------------------------------------------------
!
!     Pressure gradient scheme:
!
!       nbcs   scheme order
!
!        1        2nd (Mellor and Yamada???)
!        2        4th (McCalpin)
!
      npg = 1
!
C-----------------------------------------------------------------------
C
C     Step interval during which external (2-D) mode advective terms are
C     not updated (dimensionless):
C
      ispadv=5
C
C-----------------------------------------------------------------------
C
C     Constant in temporal filter used to prevent solution splitting
C     (dimensionless):
C
      smoth=0.10e0
C
C-----------------------------------------------------------------------
C
C     Weight used for surface slope term in external (2-D) dynamic
C     equation (a value of alpha = 0.e0 is perfectly acceptable, but the
C     value, alpha=.225e0 permits a longer time step):
C
      alpha=0.225e0
C
C-----------------------------------------------------------------------
C
C     Initial value of aam:
C
      aam_init=500.e0
C
C-----------------------------------------------------------------------
C
C     Backup gap default value                  ! rwnd:
C                                               !     :
      bkp_gap = 60                              !     :
!     Initialise time0 here to avoid overriding
      time0 = 0.e0
!     Interpolation params init
      m   = 1
      fac = 0.
!     Initialise debug stuff
      dbg_step  = 0
      dbg_seq_i = 0
!     Clock vars
      slice_s = 0.
      slice_f = 0.

!  BC flags
      bcf%wnd = .true.
      bcf%lrd = .false.
      bcf%srd = .false.
      bcf%ssf = .true.
      bcf%vap = .false.
C
C     End of input of constants
C***********************************************************************
C
C --- Above are the default parameters, alternatively one can
C --- use parameters from a file created by runscript runpom2k
C
      include 'params'
C
C***********************************************************************
C
C     Calculate some constants:
C
      dti=dte*float(isplit)
      dte2=dte*2
      dti2=dti*2
C
      iend=max0(nint(days*24.e0*3600.e0/dti),2)
      iprint=nint(prtd1*24.e0*3600.e0/dti)
      iswtch=nint(swtch*24.e0*3600.e0/dti)
C
      ispi=1.e0/float(isplit)
      isp2i=1.e0/(2.e0*float(isplit))
C
C-----------------------------------------------------------------------
C
C     Print initial summary:
C
      write(6,'(/,'' source   = '',a40)') source
      write(6,'('' title      = '',a40/)') title
      write(6,'('' iproblem   = '',i10)') iproblem
      write(6,'('' mode       = '',i10)') mode
      write(6,'('' nadv       = '',i10)') nadv
      write(6,'('' nitera     = '',i10)') nitera
      write(6,'('' sw         = '',f10.4)') sw
      write(6,'('' nread      = '',i10)') nread
      write(6,'('' dte        = '',f10.2)') dte
      write(6,'('' dti        = '',f10.1)') dti
      write(6,'('' isplit     = '',i10)') isplit
      write(6,'('' time_start = '',a26)') time_start
      write(6,'('' days       = '',f10.4)') days
      write(6,'('' iend       = '',i10)') iend
      write(6,'('' prtd1      = '',f10.4)') prtd1
      write(6,'('' iprint     = '',i10)') iprint
      write(6,'('' prtd2      = '',f10.4)') prtd2
      write(6,'('' swtch      = '',f10.2)') swtch
      write(6,'('' iswtch     = '',i10)') iswtch
      write(6,'('' iskp, jskp = '',i5'','',i5)') iskp,jskp
      write(6,'('' lramp      = '',l10)') lramp
      write(6,'('' rhoref     = '',f10.3)') rhoref
      write(6,'('' tbias      = '',f10.3)') tbias
      write(6,'('' sbias      = '',f10.3)') sbias
      write(6,'('' grav       = '',f10.4)') grav
      write(6,'('' kappa      = '',f10.4)') kappa
      write(6,'('' z0b        = '',f10.6)') z0b
      write(6,'('' cbcmin     = '',f10.6)') cbcmin
      write(6,'('' cbcmax     = '',f10.6)') cbcmax
      write(6,'('' horcon     = '',f10.3)') horcon
      write(6,'('' tprni      = '',f10.4)') tprni
      write(6,'('' umol       = '',f10.4)') umol
      write(6,'('' hmax       = '',f10.2)') hmax
      write(6,'('' vmaxl      = '',f10.4)') vmaxl
      write(6,'('' slmax      = '',f10.4)') slmax
      write(6,'('' kl1, kl2   = '',i5,'','',i5)') kl1,kl2
      write(6,'('' ntp        = '',i10)') ntp
      write(6,'('' nbct       = '',i10)') nbct
      write(6,'('' nbcs       = '',i10)') nbcs
      write(6,'('' ispadv     = '',i10)') ispadv
      write(6,'('' smoth      = '',f10.4)') smoth
      write(6,'('' alpha      = '',f10.4)') alpha
C
C-----------------------------------------------------------------------
C
C     Initialise boundary arrays:
C
      do i=1,im
        vabn(i)=0.e0
        vabs(i)=0.e0
        uabn(i)=0.e0        ! rwnd:
        uabs(i)=0.e0        !     :
        eln(i)=0.e0
        els(i)=0.e0
        do k=1,kb
          vbn(i,k)=0.e0
          vbs(i,k)=0.e0
          ubn(i,k)=0.e0     ! rwnd:
          ubs(i,k)=0.e0     !     :
          tbn(i,k)=0.e0
          tbs(i,k)=0.e0
          sbn(i,k)=0.e0
          sbs(i,k)=0.e0
        end do
      end do
C
      do j=1,jm
        uabe(j)=0.e0
        uabw(j)=0.e0
        vabw(j)=0.e0        ! rwnd:
        vabe(j)=0.e0        !     :
        ele(j)=0.e0
        elw(j)=0.e0
        do k=1,kb
          ube(j,k)=0.e0
          ubw(j,k)=0.e0
          vbe(j,k)=0.e0     ! rwnd:
          vbw(j,k)=0.e0     !     :
          tbe(j,k)=0.e0
          tbw(j,k)=0.e0
          sbe(j,k)=0.e0
          sbw(j,k)=0.e0
        end do
      end do
C
C-----------------------------------------------------------------------
C
C     Initialise 2-D and 3-D arrays for safety (this may be overwritten
C     later):
C
      do j=1,jm
        do i=1,im
          uab(i,j)=0.e0
          vab(i,j)=0.e0
          elb(i,j)=0.e0
          etb(i,j)=0.e0
          e_atmos(i,j)=0.e0
          vfluxb(i,j)=0.e0
          vfluxf(i,j)=0.e0
          wusurf(i,j)=0.e0
          wvsurf(i,j)=0.e0
          wtsurf(i,j)=0.e0
          wssurf(i,j)=0.e0
          swrad(i,j)=0.e0
          drx2d(i,j)=0.e0
          dry2d(i,j)=0.e0
        end do
      end do
C
      do k=1,kbm1
        do j=1,jm
          do i=1,im
            ub(i,j,k)=0.e0
            vb(i,j,k)=0.e0
            aam(i,j,k)=aam_init         ! rwnd : / killing pom2k_bug as ??? said /
          end do
        end do
      end do
C
C-----------------------------------------------------------------------
C
C     Set up sigma layers:
C
      if(iproblem.lt.3) call depth
C
C-----------------------------------------------------------------------
C
C     Read in grid data, and initial and lateral boundary conditions:
C
      if(iproblem.eq.1) then
        call seamount
      else if(iproblem.eq.2) then
        call box
      else if(iproblem.eq.3) then
        call file2ic
      else if(iproblem.eq.11) then  !rwnd:
        call ncdf2ic                !    :
      else if(iproblem.eq.71) then  !    :
        call cortes                 !    :
      else
        write(6,8)
    8   format(/' Invalid value of iproblem ..... program terminated'/)
        stop
      endif
C
C     Inertial period for temporal filter:
C
      period=(2.e0*pi)/abs(cor(im/2,jm/2))/86400.e0
C
C     Initialise time:
C
!      time0=0.e0       ! Do not initialise time0 here since we already set it in params
      time=0.e0
C
      ncptime=0                     ! rwnd:
C     Initialize read flags:        !     :
      rf_ts    = 0  ! bry T and S   !     :
      rf_sts   = 0  ! SST and SSS   !     :
      rf_uv    = 0  ! bry u and v   !     :
      rf_swrad = 0  !               !     :
      rf_wtsur = 0  !               !     :
      rf_wsurf = 0  !               !     :
      rf_el    = 0  !               !     :
C
C     Initial conditions:
C
C     NOTE that lateral thermodynamic boundary conditions are often set
C     equal to the initial conditions and are held constant thereafter.
C     Users can of course create variable boundary conditions.
C
      do i=1,im
        do j=1,jm
          ua(i,j)=uab(i,j)
          va(i,j)=vab(i,j)
          el(i,j)=elb(i,j)
          et(i,j)=etb(i,j)
          etf(i,j)=et(i,j)
          d(i,j)=h(i,j)+el(i,j)
          dt(i,j)=h(i,j)+et(i,j)
          w(i,j,1)=vfluxf(i,j)
        end do
      end do
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            l(i,j,k)=0.1*dt(i,j)
            q2b(i,j,k)=small
            q2lb(i,j,k)=l(i,j,k)*q2b(i,j,k)
            kh(i,j,k)=l(i,j,k)*sqrt(q2b(i,j,k))
            km(i,j,k)=kh(i,j,k)
            kq(i,j,k)=kh(i,j,k)
            aam(i,j,k)=aam_init
          end do
        end do
      end do
C
      do k=1,kbm1
        do i=1,im
          do j=1,jm
            q2(i,j,k)=q2b(i,j,k)
            q2l(i,j,k)=q2lb(i,j,k)
            t(i,j,k)=tb(i,j,k)
            s(i,j,k)=sb(i,j,k)
            u(i,j,k)=ub(i,j,k)
            v(i,j,k)=vb(i,j,k)
          end do
        end do
      end do
C
      call dens(s,t,rho)
C
      if (npg.eq.1) then
        call baropg
      else
        call baropg_mcc
      end if
C
      do k=1,kbm1
        do j=1,jm
          do i=1,im
            drx2d(i,j)=drx2d(i,j)+drhox(i,j,k)*dz(k)
            dry2d(i,j)=dry2d(i,j)+drhoy(i,j,k)*dz(k)
          end do
        end do
      end do
C
C     Calculate bottom friction coefficient:
C
      do j=1,jm
        do i=1,im
          cbc(i,j)=(kappa/log((1.e0+zz(kbm1))*h(i,j)/z0b))**2
          cbc(i,j)=max(cbcmin,cbc(i,j))
C
C     If the following is invoked, then it is probable that the wrong
C     choice of z0b or vertical spacing has been made:
C
          cbc(i,j)=min(cbcmax,cbc(i,j))
        end do
      end do
C
C     Calculate external (2-D) CFL time step:
C
      do j=1,jm
        do i=1,im
          tps(i,j)=0.5e0/sqrt(1.e0/dx(i,j)**2+1.e0/dy(i,j)**2)
     $               /sqrt(grav*(h(i,j)+small))*fsm(i,j)
        end do
      end do
C
C-----------------------------------------------------------------------
C
C     The following data are needed for a seamless restart. if nread=1,
C     data had been created by a previous run (see write(71) at end of
C     this program). nread=0 denotes a first time run.
C
      if(nread.eq.1) then
        open(70, file=trim(pth_wrk)//trim(ptf_rst), action='read'
     $          ,form='unformatted')
        read(70) time0,
     $           wubot,wvbot,aam2d,ua,uab,va,vab,el,elb,et,etb,egb,
     $           utb,vtb,u,ub,w,v,vb,t,tb,s,sb,rho,
     $           adx2d,ady2d,advua,advva,
     $           km,kh,kq,l,q2,q2b,aam,q2l,q2lb
        close(70)
      end if
C
      do j=1,jm
        do i=1,im
          d(i,j)=h(i,j)+el(i,j)
          dt(i,j)=h(i,j)+et(i,j)
        end do
      end do
C
      time=time0
C
C-----------------------------------------------------------------------
C
C     Initial conditions:
C
!     Update boundary conditions:
      call upd_mnth
      if (nbct.eq.2.or.nbct.eq.4) call flux(3)                 ! rwnd: swrad
      if (nbct.eq.1.or.nbct.eq.2.or.nbct.eq.4) call flux(4)    ! rwnd: wtsurf
      if (nbct.eq.3.or.nbct.eq.4) call bry(45)                 ! 43 - SST and SSS; 45 - 1st layer T and S
      if (bcf%wnd) call flux(52)     ! rwnd: wusurf
!<      call bry(2)       ! rwnd: el
!<      call bry(3)       !     : u & ua
      call bry(4)       !     : ts
C
      call ic2ncdf                  ! rwnd:
C
      ncptime = ncptime+1           ! rwnd:
      call ncflush(ncptime)         !     :
      !call ncTgtFlush(49)           !     :
      filename = trim(pth_wrk)//trim(pth_out)//
     $             trim(title)//"_tgt.csv"
      open(49, file=filename)
      write(49,'(f7.3,";",f10.5,";",f10.5)')time,u(70,155,1),v(70,155,1)
C
C-----------------------------------------------------------------------
C
C     Initialise netCDF output and output initial set of data:
C
        if(netcdf_file.ne.'nonetcdf') then
c     call write_netcdf(netcdf_file,1)                        ! *netCDF*
c     call write_netcdf(netcdf_file,2)                        ! *netCDF*
        endif
C
C-----------------------------------------------------------------------
!
      if (dbg_lvl.eq.10) then
        call debug_write_xy(aam2d, "aam2d", 0)
        call debug_write_xy(advua, "advua", 0)
        call debug_write_xy(advva, "advva", 0)
        call debug_write_xy(adx2d, "adx2d", 0)
        call debug_write_xy(ady2d, "ady2d", 0)
        call debug_write_xy(d,     "d"    , 0)
        call debug_write_xy(dt,    "dt"   , 0)
        call debug_write_xy(drx2d, "drx2d", 0)
        call debug_write_xy(dry2d, "dry2d", 0)
        call debug_write_xy(egf,   "egf",   0)
        call debug_write_xy(el,    "el",    0)
        call debug_write_xy(elb,   "elb",   0)
        call debug_write_xy(elf,   "elf",   0)
        call debug_write_xy(et,    "et",    0)
        call debug_write_xy(etb,   "etb",   0)
        call debug_write_xy(etf,   "etf",   0)
        call debug_write_xy(egf,   "egf",   0)
        call debug_write_xy(fluxua,"fluxua",0)
        call debug_write_xy(fluxva,"fluxva",0)
      end if
!
!      write(*,*) time, " [", iprint, "]"    ! rwnd: TODO: remove from `production` or always output to file
      call cpu_time(slice_s)
C
      do 9000 iint=1,iend      !  Begin internal (3-D) mode
C
!        write(*,*) "** ",iint," / ", iend, " **"
C
        time=dti*float(iint)/86400.e0+time0
!
        call upd_mnth   ! rwnd: Update month for BC reading
C
        if(lramp) then
          ramp=time/period
          if(ramp.gt.1.e0) ramp=1.e0
        else
          ramp=1.e0
        endif
C
!
        if (time.ge.dbg_off) dbg_step = dbg_step + 1
        dbg_seq_i = 0
!
C       write(6,2) mode,iint,time
C   2   format(' mode,iint,time =',2i5,f9.2)
C
C-----------------------------------------------------------------------
C
C     Set time dependent, surface and lateral boundary conditions.
C     The latter will be used in subroutine bcond. Users may
C     wish to create a subroutine to supply wusurf, wvsurf, wtsurf,
C     wssurf, swrad and vflux.
        if (nbct.eq.2.or.nbct.eq.4) call flux(3)            ! rwnd: swrad
        if (nbct.eq.1.or.nbct.eq.2.or.nbct.eq.4) call flux(4)    ! rwnd: wtsurf
        if (nbct.eq.3.or.nbct.eq.4) call bry(45)  ! 43 - SST and SSS
        if (bcf%wnd) call flux(52) !   51 - 2000 4xDaily; 52 - Climatological Daily
!<        call bry(2)     ! elevation BC
!<        call bry(3)     ! u,v BC
        call bry(4)     ! t,s BC
!
        if ((dbg_lvl.eq.10).and.(time.ge.dbg_off)) then
          call debug_write_xy(swrad,  "swrad",  dbg_step)
          call debug_write_xy(wtsurf, "wtsurf", dbg_step)
          call debug_write_xy(tsurf,  "tsurf",  dbg_step)
          call debug_write_xy(ssurf,  "ssurf",  dbg_step)
          call debug_write_xy(wusurf, "wusurf", dbg_step)
          call debug_write_xy(wvsurf, "wvsurf", dbg_step)
          call debug_write_xy(elb,    "elb",    dbg_step)
          call debug_write_x(uabs,   "uabs",   dbg_step)
          call debug_write_x(vabs,   "vabs",   dbg_step)
          call debug_write_xz(ubs,    "ubs",    dbg_step)
          call debug_write_xz(vbs,    "vbs",    dbg_step)
          call debug_write_xz(tbs,    "tbs",    dbg_step)
          call debug_write_xz(sbs,    "sbs",    dbg_step)
        end if
!
C
C     Introduce simple wind stress. Value is negative for westerly or
C     southerly winds. The following wind stress has been tapered
C     along the boundary to suppress numerically induced oscilations
C     near the boundary (Jamart and Ozer, J.G.R., 91, 10621-10631).
C     To make a healthy surface Ekman layer, it would be well to set
C     kl1=9.
C
      if (iproblem.ne.11) then
        do j=2,jmm1
          do i=2,imm1
c
      if((iproblem.ne.3).and.(iproblem.ne.11)
     $ .and.(iproblem.ne.71)) then     ! constant wind read in file2ic ! RWND[+]/ additional condition for iproblem={11,71} /
c
c           wusurf(i,j)=ramp*(1.e-4*cos(pi*(j-1)/jmm1))
            wusurf(i,j)=1.00*(1.e-4*cos(pi*(j-1)/jmm1))
     $                    *.25e0*(dvm(i,j+1)+dvm(i-1,j+1)
     $                          +dvm(i-1,j)+dvm(i,j))
C --- no wind ----
c           wusurf(i,j)=0.e0
            wvsurf(i,j)=0.e0
       endif
            e_atmos(i,j)=0.e0
            vfluxf(i,j)=0.e0
C
C     Set w(i,j,1)=vflux(i,j).ne.0 if one wishes non-zero flow across
C     the sea surface. See calculation of elf(i,j) below and subroutines
C     vertvl, advt1 (or advt2). If w(1,j,1)=0, and, additionally, there
C     is no net flow across lateral boundaries, the basin volume will be
C     constant; if also vflux(i,j).ne.0, then, for example, the average
C     salinity will change and, unrealistically, so will total salt.
C
            w(i,j,1)=vfluxf(i,j)
C
C     Set wtsurf to the sensible heat, the latent heat (which involves
C     only the evaporative component of vflux) and the long wave
C     radiation:
C
            wtsurf(i,j)=0.e0
C
C     Set swrad to the short wave radiation:
C
            swrad(i,j)=0.e0
C
C     To account for change in temperature of flow crossing the sea
C     surface (generally quite small compared to latent heat effect)
C
            tatm=t(i,j,1)+tbias    ! an approximation
            wtsurf(i,j)=wtsurf(i,j)+vfluxf(i,j)*(tatm-t(i,j,1)-tbias)
C
C     Set the salinity of water vapor/precipitation which enters/leaves
C     the atmosphere (or e.g., an ice cover)
C
            satm=0.e0
            wssurf(i,j)=            vfluxf(i,j)*(satm-s(i,j,1)-sbias)
C
          end do
        end do
      end if
C
C-----------------------------------------------------------------------
C
C     Set lateral viscosity:
C
C     If mode=2 then initial values of aam2d are used. If one wishes
C     to use Smagorinsky lateral viscosity and diffusion for an
C     external (2-D) mode calculation, then appropiate code can be
C     adapted from that below and installed just before the end of the
C     "if(mode.eq.2)" loop in subroutine advave.
C
C     Calculate Smagorinsky lateral viscosity:
C
C       ( hor visc = horcon*dx*dy*sqrt((du/dx)**2+(dv/dy)**2
C                                     +.5*(du/dy+dv/dx)**2) )
C
        if(mode.ne.2) then
          call advct(a,c,ee)
          if (npg.eq.1) then
            call baropg
          else
            call baropg_mcc
          end if
C
!
          if ((dbg_lvl.eq.10).and.(time.ge.dbg_off)) then
            call debug_write_xy(advx,  "advx",  dbg_step)   ! rwnd:
            call debug_write_xy(advy,  "advy",  dbg_step)   !     :
            call debug_write_xy(drhox, "drhox", dbg_step)   !     :
            call debug_write_xy(drhoy, "drhoy", dbg_step)   !     :

            call debug_write_xy(uaf,   "uaf",   dbg_step)   ! rwnd:
          end if
!
          do k=1,kbm1
            do j=2,jmm1
              do i=2,imm1
                aam(i,j,k)=horcon*dx(i,j)*dy(i,j)
     $                      *sqrt( ((u(i+1,j,k)-u(i,j,k))/dx(i,j))**2
     $                            +((v(i,j+1,k)-v(i,j,k))/dy(i,j))**2
     $                      +.5e0*(.25e0*(u(i,j+1,k)+u(i+1,j+1,k)
     $                                   -u(i,j-1,k)-u(i+1,j-1,k))
     $                      /dy(i,j)
     $                      +.25e0*(v(i+1,j,k)+v(i+1,j+1,k)
     $                             -v(i-1,j,k)-v(i-1,j+1,k))
     $                      /dx(i,j)) **2)
              end do
            end do
          end do
!
          if ((dbg_lvl.eq.10).and.(time.ge.dbg_off)) then
            call debug_write_3d(aam,   "aam",   dbg_step)   !     :
          end if
!
C
C     Form vertical averages of 3-D fields for use in external (2-D)
C     mode:
C
          do j=1,jm
            do i=1,im
              adx2d(i,j)=0.e0
              ady2d(i,j)=0.e0
              drx2d(i,j)=0.e0
              dry2d(i,j)=0.e0
              aam2d(i,j)=0.e0
            end do
          end do
C
          do k=1,kbm1
            do j=1,jm
              do i=1,im
                adx2d(i,j)=adx2d(i,j)+advx(i,j,k)*dz(k)
                ady2d(i,j)=ady2d(i,j)+advy(i,j,k)*dz(k)
                drx2d(i,j)=drx2d(i,j)+drhox(i,j,k)*dz(k)
                dry2d(i,j)=dry2d(i,j)+drhoy(i,j,k)*dz(k)
                aam2d(i,j)=aam2d(i,j)+aam(i,j,k)*dz(k)
              end do
            end do
          end do
C
          call advave(tps)
C
          do j=1,jm
            do i=1,im
              adx2d(i,j)=adx2d(i,j)-advua(i,j)
              ady2d(i,j)=ady2d(i,j)-advva(i,j)
            end do
          end do
C
        endif
C
        do j=1,jm
          do i=1,im
            egf(i,j)=el(i,j)*ispi
          end do
        end do
C
        do j=1,jm
          do i=2,im
            utf(i,j)=ua(i,j)*(d(i,j)+d(i-1,j))*isp2i
          end do
        end do
        do j=2,jm
          do i=1,im
            vtf(i,j)=va(i,j)*(d(i,j)+d(i,j-1))*isp2i
          end do
        end do
!
        if ((dbg_lvl.eq.10).and.(time.ge.dbg_off)) then
          call debug_write_xy(advua, "advua", dbg_step)   ! rwnd:
          call debug_write_xy(advva, "advva", dbg_step)   !     :
          call debug_write_xy(adx2d, "adx2d", dbg_step)   !     :
          call debug_write_xy(ady2d, "ady2d", dbg_step)   !     :
          call debug_write_xy(egf,   "egf",   dbg_step)   !     :
          call debug_write_xy(utf,   "utf",   dbg_step)   !     :
          call debug_write_xy(vtf,   "vtf",   dbg_step)   !     :
        end if
!
C
C-----------------------------------------------------------------------
C
        do 8000 iext=1,isplit    ! Begin external (2-D) mode
C
C         write(6,3) iext,time
C   3     format(' iext,time =',i5,f9.2)
C
          do j=2,jm
            do i=2,im
              fluxua(i,j)=.25e0*(d(i,j)+d(i-1,j))
     $                     *(dy(i,j)+dy(i-1,j))*ua(i,j)
              fluxva(i,j)=.25e0*(d(i,j)+d(i,j-1))
     $                     *(dx(i,j)+dx(i,j-1))*va(i,j)
            end do
          end do
!
        if ((dbg_lvl.eq.10).and.(time.ge.dbg_off)) then
          call debug_write_xy(fluxua,"fluxua",dbg_step)   ! rwnd:
          call debug_write_xy(fluxva,"fluxva",dbg_step)   !     :
          call debug_write_xy(art,   "art",   dbg_step)   !     :
          call debug_write_xy(elb,   "elb",   dbg_step)   !     :
          call debug_write_xy(vfluxf,"vfluxf",dbg_step)   !     :
        end if
!
C
C     NOTE addition of surface freshwater flux, w(i,j,1)=vflux, compared
C     with pom98.f. See also modifications to subroutine vertvl.
C
          do j=2,jmm1
            do i=2,imm1
              elf(i,j)=elb(i,j)
     $                  +dte2*(-(fluxua(i+1,j)-fluxua(i,j)
     $                          +fluxva(i,j+1)-fluxva(i,j))/art(i,j)
     $                          -vfluxf(i,j))
            end do
          end do
C
!
          if ((dbg_lvl.ge.3).and.(time.ge.dbg_off)) then
            call debug_write_txt_xy(elf,   "elf",   dbg_step)   ! rwnd:
            call debug_write_txt_xy(fluxua,"fluxua",dbg_step)   ! rwnd:
            call debug_write_txt_xy(fluxva,"fluxva",dbg_step)   !     :
            call debug_write_txt_xy(va,    "va",    dbg_step)   !     :
            call debug_write_txt_xy(elb,   "elb",   dbg_step)   !     :
          end if
!
          call bcond(1)
!
          if ((dbg_lvl.eq.10).and.(time.ge.dbg_off)) then
            call debug_write_xy(elf,"elf",dbg_step)   ! rwnd:
          end if
!
          if(mod(iext,ispadv).eq.0) call advave(tps)
!
          if ((dbg_lvl.eq.10).and.(time.ge.dbg_off)) then
            call debug_write_xy(adx2d,"adx2d",dbg_step)   ! rwnd:
            call debug_write_xy(advua,"advua",dbg_step)   ! rwnd:
            call debug_write_xy(va,"va",dbg_step)   ! rwnd:
            call debug_write_xy(el,"el",dbg_step)   ! rwnd:
            call debug_write_xy(elb,"elb",dbg_step)   ! rwnd:
            call debug_write_xy(e_atmos,"e_atmos",dbg_step)   ! rwnd:
            call debug_write_xy(drx2d,"drx2d",dbg_step)   ! rwnd:
            call debug_write_xy(wusurf,"wusurf",dbg_step)   ! rwnd:
            call debug_write_xy(wubot,"wubot",dbg_step)   ! rwnd:
            call debug_write_xy(uab,"uab",dbg_step)   ! rwnd:
            call debug_write_xy(uaf,"uaf",dbg_step)   ! rwnd:
          end if
!
C
          do j=2,jmm1
            do i=2,im
              uaf(i,j)=adx2d(i,j)+advua(i,j)
     $                  -aru(i,j)*.25e0
     $                    *(cor(i,j)*d(i,j)*(va(i,j+1)+va(i,j))
     $                     +cor(i-1,j)*d(i-1,j)*(va(i-1,j+1)+va(i-1,j)))
     $                  +.25e0*grav*(dy(i,j)+dy(i-1,j))
     $                    *(d(i,j)+d(i-1,j))
     $                    *((1.e0-2.e0*alpha)
     $                       *(el(i,j)-el(i-1,j))
     $                      +alpha*(elb(i,j)-elb(i-1,j)
     $                             +elf(i,j)-elf(i-1,j))
     $                      +e_atmos(i,j)-e_atmos(i-1,j))
     $                  +drx2d(i,j)+aru(i,j)*(wusurf(i,j)-wubot(i,j))
            end do
          end do
C
          do j=2,jmm1
            do i=2,im
              uaf(i,j)=((h(i,j)+elb(i,j)+h(i-1,j)+elb(i-1,j))
     $                    *aru(i,j)*uab(i,j)
     $                  -4.e0*dte*uaf(i,j))
     $                 /((h(i,j)+elf(i,j)+h(i-1,j)+elf(i-1,j))
     $                     *aru(i,j))
            end do
          end do
C
          do j=2,jm
            do i=2,imm1
              vaf(i,j)=ady2d(i,j)+advva(i,j)
     $                  +arv(i,j)*.25e0
     $                    *(cor(i,j)*d(i,j)*(ua(i+1,j)+ua(i,j))
     $                     +cor(i,j-1)*d(i,j-1)*(ua(i+1,j-1)+ua(i,j-1)))
     $                  +.25e0*grav*(dx(i,j)+dx(i,j-1))
     $                    *(d(i,j)+d(i,j-1))
     $                    *((1.e0-2.e0*alpha)*(el(i,j)-el(i,j-1))
     $                      +alpha*(elb(i,j)-elb(i,j-1)
     $                             +elf(i,j)-elf(i,j-1))
     $                      +e_atmos(i,j)-e_atmos(i,j-1))
     $                  +dry2d(i,j)+arv(i,j)*(wvsurf(i,j)-wvbot(i,j))
            end do
          end do
C
          do j=2,jm
            do i=2,imm1
              vaf(i,j)=((h(i,j)+elb(i,j)+h(i,j-1)+elb(i,j-1))
     $                    *vab(i,j)*arv(i,j)
     $                  -4.e0*dte*vaf(i,j))
     $                 /((h(i,j)+elf(i,j)+h(i,j-1)+elf(i,j-1))
     $                     *arv(i,j))
            end do
          end do
C
!
          if ((dbg_lvl.eq.10).and.(time.ge.dbg_off)) then
            call debug_write_xy(uaf,"uaf",dbg_step)   ! rwnd:
            call debug_write_xy(vaf,"vaf",dbg_step)   ! rwnd:
          end if
!
          call bcond(2)
!
          if ((dbg_lvl.eq.10).and.(time.ge.dbg_off)) then
            call debug_write_xy(uaf,"uaf",dbg_step)   ! rwnd:
            call debug_write_xy(vaf,"vaf",dbg_step)   ! rwnd:
          end if
!
C
          if(iext.eq.(isplit-2))then
            do j=1,jm
              do i=1,im
                etf(i,j)=.25e0*smoth*elf(i,j)
              end do
            end do
C
          else if(iext.eq.(isplit-1)) then
C
            do j=1,jm
              do i=1,im
                etf(i,j)=etf(i,j)+.5e0*(1.-.5e0*smoth)*elf(i,j)
              end do
            end do
C
          else if(iext.eq.isplit) then
C
            do j=1,jm
              do i=1,im
                etf(i,j)=(etf(i,j)+.5e0*elf(i,j))*fsm(i,j)
              end do
            end do
C
!
          if ((dbg_lvl.eq.10).and.(time.ge.dbg_off)) then
            call debug_write_xy(etf,"etf",dbg_step)   ! rwnd:
          end if
!
          endif
C
C     Stop if velocity condition violated (generally due to CFL
C     criterion not being satisfied):
C
          vamax=0.e0
C
          do j=1,jm
            do i=1,im
              if(abs(vaf(i,j)).ge.vamax) then
                vamax=abs(vaf(i,j))
	        imax=i
	        jmax=j
              endif
            end do
          end do
C
          if(vamax.le.vmaxl) then
C
C     Apply filter to remove time split and reset time sequence:
C
            do j=1,jm
              do i=1,im
                ua(i,j)=ua(i,j)
     $                   +.5e0*smoth*(uab(i,j)-2.e0*ua(i,j)+uaf(i,j))
                va(i,j)=va(i,j)
     $                   +.5e0*smoth*(vab(i,j)-2.e0*va(i,j)+vaf(i,j))
                el(i,j)=el(i,j)
     $                   +.5e0*smoth*(elb(i,j)-2.e0*el(i,j)+elf(i,j))
                elb(i,j)=el(i,j)
                el(i,j)=elf(i,j)
                d(i,j)=h(i,j)+el(i,j)
                uab(i,j)=ua(i,j)
                ua(i,j)=uaf(i,j)
                vab(i,j)=va(i,j)
                va(i,j)=vaf(i,j)
              end do
            end do
!
            if ((dbg_lvl.eq.10).and.(time.ge.dbg_off)) then
              call debug_write_xy(ua, "ua", dbg_step)   ! rwnd:
              call debug_write_xy(va, "va", dbg_step)   ! rwnd:
              call debug_write_xy(uab,"uab",dbg_step)   ! rwnd:
              call debug_write_xy(vab,"vab",dbg_step)   ! rwnd:
              call debug_write_xy(el, "el", dbg_step)   ! rwnd:
              call debug_write_xy(elb,"elb",dbg_step)   ! rwnd:
              call debug_write_xy(d,  "d",  dbg_step)   ! rwnd:
            end if
!
C
            if(iext.ne.isplit) then
              do j=1,jm
                do i=1,im
                  egf(i,j)=egf(i,j)+el(i,j)*ispi
                end do
              end do
              do j=1,jm
                do i=2,im
                  utf(i,j)=utf(i,j)+ua(i,j)*(d(i,j)+d(i-1,j))*isp2i
                end do
              end do
              do j=2,jm
                do i=1,im
                  vtf(i,j)=vtf(i,j)+va(i,j)*(d(i,j)+d(i,j-1))*isp2i
                end do
              end do
            endif
!
            if ((dbg_lvl.ge.3).and.(time.ge.dbg_off)) then
              call debug_write_txt_xy(ua,"ua",dbg_step)   ! rwnd:
              call debug_write_txt_xy(va,"va",dbg_step)   ! rwnd:
            end if
!
C
          endif
C
 8000 continue        ! End of external (2-D) mode
C
C-----------------------------------------------------------------------
C
        if(vamax.le.vmaxl) then
C
C     Continue with internal (3-D) mode calculation:
C
          if((iint.ne.1.or.time0.ne.0.e0).and.mode.ne.2) then
C
C     Adjust u(z) and v(z) such that depth average of (u,v) = (ua,va):
C
            do j=1,jm
              do i=1,im
                tps(i,j)=0.e0
              end do
            end do
C
            do k=1,kbm1
              do j=1,jm
                do i=1,im
                  tps(i,j)=tps(i,j)+u(i,j,k)*dz(k)
                end do
              end do
            end do
C
            do k=1,kbm1
              do j=1,jm
                do i=2,im
                  u(i,j,k)=(u(i,j,k)-tps(i,j))+
     $                     (utb(i,j)+utf(i,j))/(dt(i,j)+dt(i-1,j))
                end do
              end do
            end do
C
            do j=1,jm
              do i=1,im
                tps(i,j)=0.e0
              end do
            end do
C
            do k=1,kbm1
              do j=1,jm
                do i=1,im
                  tps(i,j)=tps(i,j)+v(i,j,k)*dz(k)
                end do
              end do
            end do
C
            do k=1,kbm1
              do j=2,jm
                do i=1,im
                  v(i,j,k)=(v(i,j,k)-tps(i,j))+
     $                     (vtb(i,j)+vtf(i,j))/(dt(i,j)+dt(i,j-1))
                end do
              end do
            end do
!
            if ((dbg_lvl.eq.10).and.(time.ge.dbg_off)) then
              call debug_write_3d(u,"u",dbg_step)   ! rwnd:
              call debug_write_3d(v,"v",dbg_step)   ! rwnd:
            end if
!
C
C     vertvl calculates w from u, v, dt (h+et), etf and etb:
C
            call vertvl(a,c)
            call bcond(5)
!
            if ((dbg_lvl.eq.10).and.(time.ge.dbg_off)) then
              call debug_write_3d(w,    "w",    dbg_step)   ! rwnd:
            end if
!
C
C
            do k=1,kb
              do j=1,jm
                do i=1,im
                  uf(i,j,k)=0.e0
                  vf(i,j,k)=0.e0
                end do
              end do
            end do
C
C     Calculate q2f and q2lf using uf, vf, a and c as temporary
C     variables:
C
            call advq(q2b,q2,uf,a,c)
            call advq(q2lb,q2l,vf,a,c)
            call profq(a,c,tps,dtef)
            call bcond(6)
!
            if ((dbg_lvl.eq.10).and.(time.ge.dbg_off)) then
              call debug_write_3d(w,    "w",    dbg_step)   ! rwnd:
              call debug_write_3d(km,   "km",   dbg_step)   ! rwnd:
              call debug_write_3d(kh,   "kh",   dbg_step)   ! rwnd:
              call debug_write_3d(uf,   "q2",   dbg_step)   ! rwnd:
              call debug_write_3d(vf,   "q2l",  dbg_step)   ! rwnd:
            end if
!
C
            do k=1,kb
              do j=1,jm
                do i=1,im
                  q2(i,j,k)=q2(i,j,k)
     $                       +.5e0*smoth*(uf(i,j,k)+q2b(i,j,k)
     $                                    -2.e0*q2(i,j,k))
                  q2l(i,j,k)=q2l(i,j,k)
     $                       +.5e0*smoth*(vf(i,j,k)+q2lb(i,j,k)
     $                                    -2.e0*q2l(i,j,k))
                  q2b(i,j,k)=q2(i,j,k)
                  q2(i,j,k)=uf(i,j,k)
                  q2lb(i,j,k)=q2l(i,j,k)
                  q2l(i,j,k)=vf(i,j,k)
                end do
              end do
            end do
!
            if ((dbg_lvl.eq.10).and.(time.ge.dbg_off)) then
              call debug_write_3d(q2b, "q2b", dbg_step)   ! rwnd:
              call debug_write_3d(q2,  "q2" , dbg_step)   ! rwnd:
              call debug_write_3d(q2lb,"q2lb",dbg_step)   ! rwnd:
              call debug_write_3d(q2l, "q2l", dbg_step)   ! rwnd:
            end if
!
C
C     Calculate tf and sf using uf, vf, a and c as temporary variables:
C
            if(mode.ne.4) then
C
              if(nadv.eq.1) then
C
                call advt1(tb,t,tclim,uf,a,c)
                call advt1(sb,s,sclim,vf,a,c)
C
              else if(nadv.eq.2) then
C
                call advt2(tb,t,tclim,uf,a,c,nitera,sw)
                call advt2(sb,s,sclim,vf,a,c,nitera,sw)
C
              else
C
                write(6,9)
    9           format(/'Invalid value for nadv ..... ',
     $                 'program terminated'/)
                stop
C
              endif
!
              if ((dbg_lvl.eq.10).and.(time.ge.dbg_off)) then
                call debug_write_3d(uf, "t", dbg_step)   ! rwnd:
                call debug_write_3d(vf, "s", dbg_step)   ! rwnd:
              end if
!
C
              call proft(uf,wtsurf,tsurf,nbct,tps)
              call proft(vf,wssurf,ssurf,nbcs,tps)
              call bcond(4)
!
              if ((dbg_lvl.eq.10).and.(time.ge.dbg_off)) then
                call debug_write_3d(uf, "t", dbg_step)   ! rwnd:flag0
                call debug_write_3d(vf, "s", dbg_step)   ! rwnd:
              end if
!
C
              do k=1,kb
                do j=1,jm
                  do i=1,im
                    t(i,j,k)=t(i,j,k)
     $                        +.5e0*smoth*(uf(i,j,k)+tb(i,j,k)
     $                                     -2.e0*t(i,j,k))
                    s(i,j,k)=s(i,j,k)
     $                        +.5e0*smoth*(vf(i,j,k)+sb(i,j,k)
     $                                     -2.e0*s(i,j,k))
                    tb(i,j,k)=t(i,j,k)
                    t(i,j,k)=uf(i,j,k)
                    sb(i,j,k)=s(i,j,k)
                    s(i,j,k)=vf(i,j,k)
                  end do
                end do
              end do
!
              if ((dbg_lvl.eq.10).and.(time.ge.dbg_off)) then
                call debug_write_3d(t, "t", dbg_step)   ! rwnd:
                call debug_write_3d(s, "s", dbg_step)   ! rwnd:
                call debug_write_3d(tb,"tb",dbg_step)   ! rwnd:
                call debug_write_3d(sb,"sb",dbg_step)   ! rwnd:
              end if
!
C
              call dens(s,t,rho)
!
              if ((dbg_lvl.eq.10).and.(time.ge.dbg_off)) then
                call debug_write_3d(rho, "rho", dbg_step)   ! rwnd:
              end if
!
C
            endif
C
C     Calculate uf and vf:
C
            call advu
            call advv
            call profu
            call profv
            call bcond(3)
!
              if ((dbg_lvl.eq.10).and.(time.ge.dbg_off)) then
                call debug_write_3d(uf, "uf", dbg_step)   ! rwnd:
                call debug_write_3d(vf, "vf", dbg_step)   ! rwnd:
              end if
!
C
            do j=1,jm
              do i=1,im
                tps(i,j)=0.e0
              end do
            end do
C
            do k=1,kbm1
              do j=1,jm
                do i=1,im
                  tps(i,j)=tps(i,j)
     $                      +(uf(i,j,k)+ub(i,j,k)-2.e0*u(i,j,k))*dz(k)
                end do
              end do
            end do
C
            do k=1,kbm1
              do j=1,jm
                do i=1,im
                  u(i,j,k)=u(i,j,k)
     $                      +.5e0*smoth*(uf(i,j,k)+ub(i,j,k)
     $                                   -2.e0*u(i,j,k)-tps(i,j))
                end do
              end do
            end do
C
            do j=1,jm
              do i=1,im
                tps(i,j)=0.e0
              end do
            end do
C
            do k=1,kbm1
              do j=1,jm
                do i=1,im
                  tps(i,j)=tps(i,j)
     $                      +(vf(i,j,k)+vb(i,j,k)-2.e0*v(i,j,k))*dz(k)
                end do
              end do
            end do
C
            do k=1,kbm1
              do j=1,jm
                do i=1,im
                  v(i,j,k)=v(i,j,k)
     $                      +.5e0*smoth*(vf(i,j,k)+vb(i,j,k)
     $                                   -2.e0*v(i,j,k)-tps(i,j))
                end do
              end do
            end do
C
            do k=1,kb
              do j=1,jm
                do i=1,im
                  ub(i,j,k)=u(i,j,k)
                  u(i,j,k)=uf(i,j,k)
                  vb(i,j,k)=v(i,j,k)
                  v(i,j,k)=vf(i,j,k)
                end do
              end do
            end do
!
            if ((dbg_lvl.eq.10).and.(time.ge.dbg_off)) then
              call debug_write_3d(u, "u", dbg_step)   ! rwnd:
              call debug_write_3d(v, "v", dbg_step)   ! rwnd:
              call debug_write_3d(ub,"ub",dbg_step)   ! rwnd:
              call debug_write_3d(vb,"vb",dbg_step)   ! rwnd:
            end if
!
C
          endif
C
          do j=1,jm
            do i=1,im
              egb(i,j)=egf(i,j)
              etb(i,j)=et(i,j)
              et(i,j)=etf(i,j)
              dt(i,j)=h(i,j)+et(i,j)
              utb(i,j)=utf(i,j)
              vtb(i,j)=vtf(i,j)
              vfluxb(i,j)=vfluxf(i,j)
            end do
          end do
!
          if ((dbg_lvl.eq.10).and.(time.ge.dbg_off)) then
            call debug_write_xy(egb,   "egb",   dbg_step)   ! rwnd:flag0
            call debug_write_xy(etb,   "etb",   dbg_step)   ! rwnd:
            call debug_write_xy(et,    "et",    dbg_step)   ! rwnd:
            call debug_write_xy(dt,    "dt",    dbg_step)   ! rwnd:
            call debug_write_xy(utb,   "utb",   dbg_step)   ! rwnd:
            call debug_write_xy(vtb,   "vtb",   dbg_step)   ! rwnd:
            call debug_write_xy(vfluxb,"vfluxb",dbg_step)   ! rwnd:
          end if
!
C
        endif
C
C-----------------------------------------------------------------------
C
C     Beginning of print section:
C
        if(iint.ge.iswtch) iprint=nint(prtd2*24.e0*3600.e0/dti)
C
        if(mod(iint,3).eq.0.or.vamax.gt.vmaxl) then
!          call ncTgtFlush(49)
          write(49, *) time, ";", u(70,155,1), ";", v(70,155,1)
        end if
        if(mod(iint,iprint).eq.0.or.vamax.gt.vmaxl) then
C
          call cpu_time(slice_f)
          write(6,4) time,iint,iext,iprint,(slice_f-slice_s)/60.
    4     format(/
     $    '**************************************************',
     $    '**************************************************',
     $    '*************************'//
     $    ' time =',f9.4,', iint =',i8,', iext =',i8,', iprint =',i8,//
     $    ' exectime =',f9.4,' minutes')
C
C     Select print statements in printall as desired:
C
          ncptime = ncptime+1           ! rwnd:
          call ncflush(ncptime)         !     :
C
          vtot=0.e0
          atot=0.e0
          taver=0.e0
          saver=0.e0
          eaver=0.e0
          do k=1,kbm1
            do j=1,jm
              do i=1,im
                darea=dx(i,j)*dy(i,j)*fsm(i,j)
                dvol=darea*dt(i,j)*dz(k)
                vtot=vtot+dvol
                taver=taver+tb(i,j,k)*dvol
                saver=saver+sb(i,j,k)*dvol
              end do
            end do
          end do
C
          do j=1,jm
            do i=1,im
              darea=dx(i,j)*dy(i,j)*fsm(i,j)
              atot=atot+darea
              eaver=eaver+et(i,j)*darea
            end do
          end do
C
          taver=taver/vtot
          saver=saver/vtot
          eaver=eaver/atot
          tsalt=(saver+sbias)*vtot
C
          write(6,5) vtot,atot,eaver,taver,saver,tsalt
    5     format('vtot = ',e16.7,'   atot = ',e16.7,
     $           '  eaver =',e16.7/'taver =',e16.7,
     $           '   saver =',e16.7,'  tsalt =',e16.7)
C
C     Write netCDF output:
C
            if(netcdf_file.ne.'nonetcdf') then
c         call write_netcdf(netcdf_file,2)                    ! *netCDF*
            endif
C
          if(vamax.gt.vmaxl) then
C
            write(6,4) time,iint,iext,iprint
C
            ncptime = ncptime+1         ! rwnd:
            call ncflush(ncptime)       !     :
C
            write(6,6) vamax,imax,jmax
    6       format(///////////////////
     $             '************************************************'/
     $             '************ abnormal job end ******************'/
     $             '************* user terminated ******************'/
     $             '************************************************'/
     $             ' vamax =',e12.3,'   imax,jmax =',2i5)
C
C     Close netCDF file:
C
              if(netcdf_file.ne.'nonetcdf') then
c           call write_netcdf(netcdf_file,3)                  ! *netCDF*
              endif
C
            stop
C
          endif

          call cpu_time(slice_s)
C
        endif
C
C     End of print section
C
C----------------------------------------------------------
C
C     Back up if needed
C
        if (mod(time,float(bkp_gap)).eq.0) then
          if (iint.ne.iend) then
            write(*,*) "[*] Backing up..."
            open(71,file=trim(pth_wrk)//trim(pth_bkp)//
     $              trim(title)//'.restart.bkp',
     $              form='unformatted',position='rewind')
            write(71) time,
     $      wubot,wvbot,aam2d,ua,uab,va,vab,el,elb,et,etb,egb,
     $      utb,vtb,u,ub,w,v,vb,t,tb,s,sb,rho,adx2d,ady2d,advua,advva,
     $      km,kh,kq,l,q2,q2b,aam,q2l,q2lb
            close(71)
          end if
        end if
C-----------------------------------------------------------------------
C
 9000 continue       !  End of internal (3-D) mode
C
C-----------------------------------------------------------------------
C
!      call ncTgtFlush(49)
      write(49, *) time, ";", u(70,155,1), ";", v(70,155,1)
      close(49)
!
      write(6,4) time,iint,iext,iprint
C
C     Save this data for a seamless restart:
C
      write(*,*) "[ ] Writing restart file..."                      ! rwnd:
      open(71,file=trim(pth_wrk)//trim(pth_out)                     !     :
     $             //trim(title)//'.restart',                       !     :
     $             form='unformatted',position='rewind')            !     :
      write(71) time,                                               !     :
     $  wubot,wvbot,aam2d,ua,uab,va,vab,el,elb,et,etb,egb,          !     :
     $  utb,vtb,u,ub,w,v,vb,t,tb,s,sb,rho,adx2d,ady2d,advua,advva,  !     :
     $  km,kh,kq,l,q2,q2b,aam,q2l,q2lb                              !     :
      close(71)                                                     !     :
!
      write(*,*) "[O] Simulations completed!"                       ! rwnd:
C
C     Close netCDF file:
C
        if(netcdf_file.ne.'nonetcdf') then
c     call write_netcdf(netcdf_file,3)                        ! *netCDF*
        endif
C
      stop
C
      end
C
C     End of main program
C
C-----------------------------------------------------------------------
C
      subroutine advave(curv2d)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Calculates horizontal advection and diffusion.      *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real curv2d(im,jm)
      integer i,j
C
C     u-advection and diffusion:
C
C     Advective fluxes:
C
      do j=1,jm
        do i=1,im
          advua(i,j)=0.e0
        end do
      end do
C
      do j=2,jm
        do i=2,imm1
          fluxua(i,j)=.125e0*((d(i+1,j)+d(i,j))*ua(i+1,j)
     $                       +(d(i,j)+d(i-1,j))*ua(i,j))
     $                      *(ua(i+1,j)+ua(i,j))
        end do
      end do
C
      do j=2,jm
        do i=2,im
          fluxva(i,j)=.125e0*((d(i,j)+d(i,j-1))*va(i,j)
     $                       +(d(i-1,j)+d(i-1,j-1))*va(i-1,j))
     $                      *(ua(i,j)+ua(i,j-1))
        end do
      end do
C
C     Add viscous fluxes:
C
      do j=2,jm
        do i=2,imm1
          fluxua(i,j)=fluxua(i,j)
     $                 -d(i,j)*2.e0*aam2d(i,j)*(uab(i+1,j)-uab(i,j))
     $                   /dx(i,j)
        end do
      end do
C
      do j=2,jm
        do i=2,im
          tps(i,j)=.25e0*(d(i,j)+d(i-1,j)+d(i,j-1)+d(i-1,j-1))
     $              *(aam2d(i,j)+aam2d(i,j-1)
     $                +aam2d(i-1,j)+aam2d(i-1,j-1))
     $              *((uab(i,j)-uab(i,j-1))
     $                 /(dy(i,j)+dy(i-1,j)+dy(i,j-1)+dy(i-1,j-1))
     $               +(vab(i,j)-vab(i-1,j))
     $                 /(dx(i,j)+dx(i-1,j)+dx(i,j-1)+dx(i-1,j-1)))
          fluxua(i,j)=fluxua(i,j)*dy(i,j)
          fluxva(i,j)=(fluxva(i,j)-tps(i,j))*.25e0
     $                 *(dx(i,j)+dx(i-1,j)+dx(i,j-1)+dx(i-1,j-1))
        end do
      end do
C
      do j=2,jmm1
        do i=2,imm1
          advua(i,j)=fluxua(i,j)-fluxua(i-1,j)
     $                +fluxva(i,j+1)-fluxva(i,j)
        end do
      end do
C
C     u-advection and diffusion:
C
      do j=1,jm
        do i=1,im
          advva(i,j)=0.e0
        end do
      end do
C
C     Advective fluxes:
C
      do j=2,jm
        do i=2,im
          fluxua(i,j)=.125e0*((d(i,j)+d(i-1,j))*ua(i,j)
     $                       +(d(i,j-1)+d(i-1,j-1))*ua(i,j-1))
     $                      *(va(i-1,j)+va(i,j))
        end do
      end do
C
      do j=2,jmm1
        do i=2,im
          fluxva(i,j)=.125e0*((d(i,j+1)+d(i,j))*va(i,j+1)
     $                       +(d(i,j)+d(i,j-1))*va(i,j))
     $                      *(va(i,j+1)+va(i,j))
        end do
      end do
C
C     Add viscous fluxes:
C
      do j=2,jmm1
        do i=2,im
          fluxva(i,j)=fluxva(i,j)
     $                 -d(i,j)*2.e0*aam2d(i,j)*(vab(i,j+1)-vab(i,j))
     $                   /dy(i,j)
        end do
      end do
C
      do j=2,jm
        do i=2,im
          fluxva(i,j)=fluxva(i,j)*dx(i,j)
          fluxua(i,j)=(fluxua(i,j)-tps(i,j))*.25e0
     $                 *(dy(i,j)+dy(i-1,j)+dy(i,j-1)+dy(i-1,j-1))
        end do
      end do
C
      do j=2,jmm1
        do i=2,imm1
          advva(i,j)=fluxua(i+1,j)-fluxua(i,j)
     $                +fluxva(i,j)-fluxva(i,j-1)
        end do
      end do
C
      if(mode.eq.2) then
C
        do j=2,jmm1
          do i=2,imm1
            wubot(i,j)=-0.5e0*(cbc(i,j)+cbc(i-1,j))
     $                  *sqrt(uab(i,j)**2
     $                        +(.25e0*(vab(i,j)+vab(i,j+1)
     $                                 +vab(i-1,j)+vab(i-1,j+1)))**2)
     $                  *uab(i,j)
          end do
        end do
C
        do j=2,jmm1
          do i=2,imm1
            wvbot(i,j)=-0.5e0*(cbc(i,j)+cbc(i,j-1))
     $                  *sqrt(vab(i,j)**2
     $                        +(.25e0*(uab(i,j)+uab(i+1,j)
     $                                +uab(i,j-1)+uab(i+1,j-1)))**2)
     $                  *vab(i,j)
          end do
        end do
C
        do j=2,jmm1
          do i=2,imm1
            curv2d(i,j)=.25e0
     $                   *((va(i,j+1)+va(i,j))*(dy(i+1,j)-dy(i-1,j))
     $                    -(ua(i+1,j)+ua(i,j))*(dx(i,j+1)-dx(i,j-1)))
     $                   /(dx(i,j)*dy(i,j))
          end do
        end do
C
        do j=2,jmm1
          do i=3,imm1
            advua(i,j)=advua(i,j)-aru(i,j)*.25e0
     $                  *(curv2d(i,j)*d(i,j)
     $                    *(va(i,j+1)+va(i,j))
     $                    +curv2d(i-1,j)*d(i-1,j)
     $                    *(va(i-1,j+1)+va(i-1,j)))
          end do
        end do
C
        do j=3,jmm1
          do i=2,imm1
            advva(i,j)=advva(i,j)+arv(i,j)*.25e0
     $                  *(curv2d(i,j)*d(i,j)
     $                    *(ua(i+1,j)+ua(i,j))
     $                    +curv2d(i,j-1)*d(i,j-1)
     $                    *(ua(i+1,j-1)+ua(i,j-1)))
          end do
        end do
C
      endif
C
      return
C
      end
C
      subroutine advct(xflux,yflux,curv)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Calculates the horizontal portions of momentum      *
C *                advection well in advance of their use in advu and  *
C *                advv so that their vertical integrals (created in   *
C *                the main program) may be used in the external (2-D) *
C *                mode calculation.                                   *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real xflux(im,jm,kb),yflux(im,jm,kb)
      real curv(im,jm,kb)
      real dtaam
      integer i,j,k
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            curv(i,j,k)=0.e0
            advx(i,j,k)=0.e0
            xflux(i,j,k)=0.e0
            yflux(i,j,k)=0.e0
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            curv(i,j,k)=.25e0*((v(i,j+1,k)+v(i,j,k))
     $                         *(dy(i+1,j)-dy(i-1,j))
     $                         -(u(i+1,j,k)+u(i,j,k))
     $                         *(dx(i,j+1)-dx(i,j-1)))
     $                       /(dx(i,j)*dy(i,j))
          end do
        end do
      end do
C
C     Calculate x-component of velocity advection:
C
C     Calculate horizontal advective fluxes:
C
      do k=1,kbm1
        do j=1,jm
          do i=2,imm1
            xflux(i,j,k)=.125e0*((dt(i+1,j)+dt(i,j))*u(i+1,j,k)
     $                           +(dt(i,j)+dt(i-1,j))*u(i,j,k))
     $                         *(u(i+1,j,k)+u(i,j,k))
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=2,jm
          do i=2,im
            yflux(i,j,k)=.125e0*((dt(i,j)+dt(i,j-1))*v(i,j,k)
     $                           +(dt(i-1,j)+dt(i-1,j-1))*v(i-1,j,k))
     $                         *(u(i,j,k)+u(i,j-1,k))
          end do
        end do
      end do
C
C    Add horizontal diffusive fluxes:
C
      do k=1,kbm1
        do j=2,jm
          do i=2,imm1
            xflux(i,j,k)=xflux(i,j,k)
     $                    -dt(i,j)*aam(i,j,k)*2.e0
     $                    *(ub(i+1,j,k)-ub(i,j,k))/dx(i,j)
            dtaam=.25e0*(dt(i,j)+dt(i-1,j)+dt(i,j-1)+dt(i-1,j-1))
     $             *(aam(i,j,k)+aam(i-1,j,k)
     $               +aam(i,j-1,k)+aam(i-1,j-1,k))
            yflux(i,j,k)=yflux(i,j,k)
     $                    -dtaam*((ub(i,j,k)-ub(i,j-1,k))
     $                            /(dy(i,j)+dy(i-1,j)
     $                              +dy(i,j-1)+dy(i-1,j-1))
     $                            +(vb(i,j,k)-vb(i-1,j,k))
     $                            /(dx(i,j)+dx(i-1,j)
     $                              +dx(i,j-1)+dx(i-1,j-1)))
C
            xflux(i,j,k)=dy(i,j)*xflux(i,j,k)
            yflux(i,j,k)=.25e0*(dx(i,j)+dx(i-1,j)
     $                          +dx(i,j-1)+dx(i-1,j-1))*yflux(i,j,k)
          end do
        end do
      end do
C
C     Do horizontal advection:
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            advx(i,j,k)=xflux(i,j,k)-xflux(i-1,j,k)
     $                   +yflux(i,j+1,k)-yflux(i,j,k)
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=2,jmm1
          do i=3,imm1
            advx(i,j,k)=advx(i,j,k)
     $                   -aru(i,j)*.25e0
     $                     *(curv(i,j,k)*dt(i,j)
     $                        *(v(i,j+1,k)+v(i,j,k))
     $                       +curv(i-1,j,k)*dt(i-1,j)
     $                        *(v(i-1,j+1,k)+v(i-1,j,k)))
          end do
        end do
      end do
C
C-----------------------------------------------------------------------
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            advy(i,j,k)=0.e0
            xflux(i,j,k)=0.e0
            yflux(i,j,k)=0.e0
          end do
        end do
      end do
C
C     Calculate y-component of velocity advection:
C
C     Calculate horizontal advective fluxes:
C
      do k=1,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=.125e0*((dt(i,j)+dt(i-1,j))*u(i,j,k)
     $                           +(dt(i,j-1)+dt(i-1,j-1))*u(i,j-1,k))
     $                         *(v(i,j,k)+v(i-1,j,k))
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=2,jmm1
          do i=1,im
            yflux(i,j,k)=.125e0*((dt(i,j+1)+dt(i,j))*v(i,j+1,k)
     $                           +(dt(i,j)+dt(i,j-1))*v(i,j,k))
     $                         *(v(i,j+1,k)+v(i,j,k))
          end do
        end do
      end do
C
C    Add horizontal diffusive fluxes:
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,im
            dtaam=.25e0*(dt(i,j)+dt(i-1,j)+dt(i,j-1)+dt(i-1,j-1))
     $             *(aam(i,j,k)+aam(i-1,j,k)
     $               +aam(i,j-1,k)+aam(i-1,j-1,k))
            xflux(i,j,k)=xflux(i,j,k)
     $                    -dtaam*((ub(i,j,k)-ub(i,j-1,k))
     $                            /(dy(i,j)+dy(i-1,j)
     $                              +dy(i,j-1)+dy(i-1,j-1))
     $                            +(vb(i,j,k)-vb(i-1,j,k))
     $                            /(dx(i,j)+dx(i-1,j)
     $                              +dx(i,j-1)+dx(i-1,j-1)))
            yflux(i,j,k)=yflux(i,j,k)
     $                    -dt(i,j)*aam(i,j,k)*2.e0
     $                    *(vb(i,j+1,k)-vb(i,j,k))/dy(i,j)
C
            xflux(i,j,k)=.25e0*(dy(i,j)+dy(i-1,j)
     $                          +dy(i,j-1)+dy(i-1,j-1))*xflux(i,j,k)
            yflux(i,j,k)=dx(i,j)*yflux(i,j,k)
          end do
        end do
      end do
C
C     Do horizontal advection:
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            advy(i,j,k)=xflux(i+1,j,k)-xflux(i,j,k)
     $                   +yflux(i,j,k)-yflux(i,j-1,k)
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=3,jmm1
          do i=2,imm1
            advy(i,j,k)=advy(i,j,k)
     $                   +arv(i,j)*.25e0
     $                     *(curv(i,j,k)*dt(i,j)
     $                        *(u(i+1,j,k)+u(i,j,k))
     $                       +curv(i,j-1,k)*dt(i,j-1)
     $                        *(u(i+1,j-1,k)+u(i,j-1,k)))
          end do
        end do
      end do
C
      return
C
      end
C
      subroutine advq(qb,q,qf,xflux,yflux)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Calculates horizontal advection and diffusion, and  *
C *                vertical advection for turbulent quantities.        *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real qb(im,jm,kb),q(im,jm,kb),qf(im,jm,kb)
      real xflux(im,jm,kb),yflux(im,jm,kb)
      integer i,j,k
C
C     Do horizontal advection:
C
      do k=2,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=.125e0*(q(i,j,k)+q(i-1,j,k))
     $                    *(dt(i,j)+dt(i-1,j))*(u(i,j,k)+u(i,j,k-1))
            yflux(i,j,k)=.125e0*(q(i,j,k)+q(i,j-1,k))
     $                    *(dt(i,j)+dt(i,j-1))*(v(i,j,k)+v(i,j,k-1))
          end do
        end do
      end do
C
C     Do horizontal diffusion:
C
      do k=2,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=xflux(i,j,k)
     $                    -.25e0*(aam(i,j,k)+aam(i-1,j,k)
     $                            +aam(i,j,k-1)+aam(i-1,j,k-1))
     $                          *(h(i,j)+h(i-1,j))
     $                          *(qb(i,j,k)-qb(i-1,j,k))*dum(i,j)
     $                          /(dx(i,j)+dx(i-1,j))
            yflux(i,j,k)=yflux(i,j,k)
     $                    -.25e0*(aam(i,j,k)+aam(i,j-1,k)
     $                            +aam(i,j,k-1)+aam(i,j-1,k-1))
     $                          *(h(i,j)+h(i,j-1))
     $                          *(qb(i,j,k)-qb(i,j-1,k))*dvm(i,j)
     $                          /(dy(i,j)+dy(i,j-1))
            xflux(i,j,k)=.5e0*(dy(i,j)+dy(i-1,j))*xflux(i,j,k)
            yflux(i,j,k)=.5e0*(dx(i,j)+dx(i,j-1))*yflux(i,j,k)
          end do
        end do
      end do
C
C     Do vertical advection, add flux terms, then step forward in time:
C
      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            qf(i,j,k)=(w(i,j,k-1)*q(i,j,k-1)-w(i,j,k+1)*q(i,j,k+1))
     $                 *art(i,j)/(dz(k)+dz(k-1))
     $                 +xflux(i+1,j,k)-xflux(i,j,k)
     $                 +yflux(i,j+1,k)-yflux(i,j,k)
            qf(i,j,k)=((h(i,j)+etb(i,j))*art(i,j)
     $                 *qb(i,j,k)-dti2*qf(i,j,k))
     $                /((h(i,j)+etf(i,j))*art(i,j))
          end do
        end do
      end do
C
      return
C
      end
C
      subroutine advt1(fb,f,fclim,ff,xflux,yflux)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Integrates conservative scalar equations.           *
C *                                                                    *
C *                This is centred scheme, as originally provide in    *
C *                POM (previously called advt).                       *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real fb(im,jm,kb),f(im,jm,kb),fclim(im,jm,kb),ff(im,jm,kb)
      real xflux(im,jm,kb),yflux(im,jm,kb)
      integer i,j,k
C
      do j=1,jm
        do i=1,im
           f(i,j,kb)=f(i,j,kbm1)
           fb(i,j,kb)=fb(i,j,kbm1)
        end do
      end do
C
C     Do advective fluxes:
C
      do k=1,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=.25e0*((dt(i,j)+dt(i-1,j))
     $                          *(f(i,j,k)+f(i-1,j,k))*u(i,j,k))
            yflux(i,j,k)=.25e0*((dt(i,j)+dt(i,j-1))
     $                          *(f(i,j,k)+f(i,j-1,k))*v(i,j,k))
          end do
        end do
      end do
C
C     Add diffusive fluxes:
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            fb(i,j,k)=fb(i,j,k)-fclim(i,j,k)
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=xflux(i,j,k)
     $                    -.5e0*(aam(i,j,k)+aam(i-1,j,k))
     $                         *(h(i,j)+h(i-1,j))*tprni
     $                         *(fb(i,j,k)-fb(i-1,j,k))*dum(i,j)
     $                         /(dx(i,j)+dx(i-1,j))
            yflux(i,j,k)=yflux(i,j,k)
     $                    -.5e0*(aam(i,j,k)+aam(i,j-1,k))
     $                         *(h(i,j)+h(i,j-1))*tprni
     $                         *(fb(i,j,k)-fb(i,j-1,k))*dvm(i,j)
     $                         /(dy(i,j)+dy(i,j-1))
            xflux(i,j,k)=.5e0*(dy(i,j)+dy(i-1,j))*xflux(i,j,k)
            yflux(i,j,k)=.5e0*(dx(i,j)+dx(i,j-1))*yflux(i,j,k)
          end do
        end do
      end do
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            fb(i,j,k)=fb(i,j,k)+fclim(i,j,k)
          end do
        end do
      end do
C
C     Do vertical advection:
C
      do j=2,jmm1
        do i=2,imm1
          zflux(i,j,1)=f(i,j,1)*w(i,j,1)*art(i,j)
          zflux(i,j,kb)=0.e0
        end do
      end do
C
      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            zflux(i,j,k)=.5e0*(f(i,j,k-1)+f(i,j,k))*w(i,j,k)*art(i,j)
          end do
        end do
      end do
C
C     Add net horizontal fluxes and then step forward in time:
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            ff(i,j,k)=xflux(i+1,j,k)-xflux(i,j,k)
     $                 +yflux(i,j+1,k)-yflux(i,j,k)
     $                 +(zflux(i,j,k)-zflux(i,j,k+1))/dz(k)
C
            ff(i,j,k)=(fb(i,j,k)*(h(i,j)+etb(i,j))*art(i,j)
     $                 -dti2*ff(i,j,k))
     $                 /((h(i,j)+etf(i,j))*art(i,j))
          end do
        end do
      end do
C
      return
C
      end
C
      subroutine advt2(fb,f,fclim,ff,xflux,yflux,nitera,sw)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Integrates conservative scalar equations.           *
C *                                                                    *
C *                This is a first-order upstream scheme, which        *
C *                reduces implicit diffusion using the Smolarkiewicz  *
C *                iterative upstream scheme with an antidiffusive     *
C *                velocity.                                           *
C *                                                                    *
C *                It is based on the subroutines of Gianmaria Sannino *
C *                (Inter-university Computing Consortium, Rome, Italy)*
C *                and Vincenzo Artale (Italian National Agency for    *
C *                New Technology and Environment, Rome, Italy),       *
C *                downloaded from the POM FTP site on 1 Nov. 2001.    *
C *                The calculations have been simplified by removing   *
C *                the shock switch option. It should be noted that    *
C *                this implementation does not include cross-terms    *
C *                which are in the original formulation.              *
C *                                                                    *
C *                fb,f,fclim,ff . as used in subroutine advt1         *
C *                xflux,yflux ... working arrays used to save memory  *
C *                nitera ........ number of iterations. This should   *
C *                                be in the range 1 - 4. 1 is         *
C *                                standard upstream differencing;     *
C *                                3 adds 50% CPU time to POM.         *
C *                sw ............ smoothing parameter. This should    *
C *                                preferably be 1, but 0 < sw < 1     *
C *                                gives smoother solutions with less  *
C *                                overshoot when nitera > 1.          *
C *                                                                    *
C *                Reference:                                          *
C *                                                                    *
C *                Smolarkiewicz, P.K.; A fully multidimensional       *
C *                  positive definite advection transport algorithm   *
C *                  with small implicit diffusion, Journal of         *
C *                  Computational Physics, 54, 325-362, 1984.         *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real fb(im,jm,kb),f(im,jm,kb),fclim(im,jm,kb),ff(im,jm,kb)
      real xflux(im,jm,kb),yflux(im,jm,kb)
      real sw
      integer nitera
      real fbmem(im,jm,kb),eta(im,jm)
      real xmassflux(im,jm,kb),ymassflux(im,jm,kb),zwflux(im,jm,kb)
      integer i,j,k,itera
C
C     Calculate horizontal mass fluxes:
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            xmassflux(i,j,k)=0.e0
            ymassflux(i,j,k)=0.e0
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,im
            xmassflux(i,j,k)=0.25e0*(dy(i-1,j)+dy(i,j))
     $                             *(dt(i-1,j)+dt(i,j))*u(i,j,k)
          end do
        end do
C
        do j=2,jm
          do i=2,imm1
            ymassflux(i,j,k)=0.25e0*(dx(i,j-1)+dx(i,j))
     $                             *(dt(i,j-1)+dt(i,j))*v(i,j,k)
          end do
        end do
      end do
C
      do j=1,jm
        do i=1,im
          fb(i,j,kb)=fb(i,j,kbm1)
          eta(i,j)=etb(i,j)
        end do
      end do
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            zwflux(i,j,k)=w(i,j,k)
            fbmem(i,j,k)=fb(i,j,k)
          end do
        end do
      end do
C
C     Start Smolarkiewicz scheme:
C
      do itera=1,nitera
C
C     Upwind advection scheme:
C
        do k=1,kbm1
          do j=2,jm
            do i=2,im
              xflux(i,j,k)=0.5e0
     $                      *((xmassflux(i,j,k)+abs(xmassflux(i,j,k)))
     $                        *fbmem(i-1,j,k)+
     $                        (xmassflux(i,j,k)-abs(xmassflux(i,j,k)))
     $                        *fbmem(i,j,k))
C
              yflux(i,j,k)=0.5e0
     $                      *((ymassflux(i,j,k)+abs(ymassflux(i,j,k)))
     $                        *fbmem(i,j-1,k)+
     $                        (ymassflux(i,j,k)-abs(ymassflux(i,j,k)))
     $                        *fbmem(i,j,k))
            end do
          end do
        end do
C
        do j=2,jmm1
          do i=2,imm1
            zflux(i,j,1)=0.e0
            if(itera.eq.1) zflux(i,j,1)=w(i,j,1)*f(i,j,1)*art(i,j)
            zflux(i,j,kb)=0.e0
          end do
        end do
C
        do k=2,kbm1
          do j=2,jmm1
            do i=2,imm1
              zflux(i,j,k)=0.5e0
     $                      *((zwflux(i,j,k)+abs(zwflux(i,j,k)))
     $                       *fbmem(i,j,k)+
     $                        (zwflux(i,j,k)-abs(zwflux(i,j,k)))
     $                       *fbmem(i,j,k-1))
              zflux(i,j,k)=zflux(i,j,k)*art(i,j)
            end do
          end do
        end do
C
C     Add net advective fluxes and step forward in time:
C
        do k=1,kbm1
          do j=2,jmm1
            do i=2,imm1
              ff(i,j,k)=xflux(i+1,j,k)-xflux(i,j,k)
     $                 +yflux(i,j+1,k)-yflux(i,j,k)
     $                 +(zflux(i,j,k)-zflux(i,j,k+1))/dz(k)
              ff(i,j,k)=(fbmem(i,j,k)*(h(i,j)+eta(i,j))*art(i,j)
     $                   -dti2*ff(i,j,k))/((h(i,j)+etf(i,j))*art(i,j))
            end do
          end do
        end do
C
C     Calculate antidiffusion velocity:
C
        call smol_adif(xmassflux,ymassflux,zwflux,ff,sw)
C
        do j=1,jm
          do i=1,im
            eta(i,j)=etf(i,j)
            do k=1,kb
              fbmem(i,j,k)=ff(i,j,k)
            end do
          end do
        end do
C
C     End of Smolarkiewicz scheme
C
      end do
C
C     Add horizontal diffusive fluxes:
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            fb(i,j,k)=fb(i,j,k)-fclim(i,j,k)
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=2,jm
          do i=2,im
            xmassflux(i,j,k)=0.5e0*(aam(i,j,k)+aam(i-1,j,k))
            ymassflux(i,j,k)=0.5e0*(aam(i,j,k)+aam(i,j-1,k))
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=2,jm
          do i=2,im
           xflux(i,j,k)=-xmassflux(i,j,k)*(h(i,j)+h(i-1,j))*tprni
     $                   *(fb(i,j,k)-fb(i-1,j,k))*dum(i,j)
     $                   *(dy(i,j)+dy(i-1,j))*0.5e0/(dx(i,j)+dx(i-1,j))
           yflux(i,j,k)=-ymassflux(i,j,k)*(h(i,j)+h(i,j-1))*tprni
     $                   *(fb(i,j,k)-fb(i,j-1,k))*dvm(i,j)
     $                   *(dx(i,j)+dx(i,j-1))*0.5e0/(dy(i,j)+dy(i,j-1))
          end do
        end do
      end do
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            fb(i,j,k)=fb(i,j,k)+fclim(i,j,k)
          end do
        end do
      end do
C
C     Add net horizontal fluxes and step forward in time:
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            ff(i,j,k)=ff(i,j,k)-dti2*(xflux(i+1,j,k)-xflux(i,j,k)
     $                               +yflux(i,j+1,k)-yflux(i,j,k))
     $                           /((h(i,j)+etf(i,j))*art(i,j))
          end do
        end do
      end do
C
      return
C
      end
C
      subroutine advu
C **********************************************************************
C *                                                                    *
C * ROUTINE NAME:  advu                                                *
C *                                                                    *
C * FUNCTION    :  Does horizontal and vertical advection of           *
C *                u-momentum, and includes coriolis, surface slope    *
C *                and baroclinic terms.                               *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      integer i,j,k
C
C     Do vertical advection:
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            uf(i,j,k)=0.e0
          end do
        end do
      end do
C
      do k=2,kbm1
        do j=1,jm
          do i=2,im
            uf(i,j,k)=.25e0*(w(i,j,k)+w(i-1,j,k))
     $                     *(u(i,j,k)+u(i,j,k-1))
          end do
        end do
      end do
C
C     Combine horizontal and vertical advection with coriolis, surface
C     slope and baroclinic terms:
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            uf(i,j,k)=advx(i,j,k)
     $                 +(uf(i,j,k)-uf(i,j,k+1))*aru(i,j)/dz(k)
     $                 -aru(i,j)*.25e0
     $                   *(cor(i,j)*dt(i,j)
     $                      *(v(i,j+1,k)+v(i,j,k))
     $                     +cor(i-1,j)*dt(i-1,j)
     $                       *(v(i-1,j+1,k)+v(i-1,j,k)))
     $                 +grav*.125e0*(dt(i,j)+dt(i-1,j))
     $                   *(egf(i,j)-egf(i-1,j)+egb(i,j)-egb(i-1,j)
     $                     +(e_atmos(i,j)-e_atmos(i-1,j))*2.e0)
     $                   *(dy(i,j)+dy(i-1,j))
     $                 +drhox(i,j,k)
          end do
        end do
      end do
C
C     Step forward in time:
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            uf(i,j,k)=((h(i,j)+etb(i,j)+h(i-1,j)+etb(i-1,j))
     $                 *aru(i,j)*ub(i,j,k)
     $                 -2.e0*dti2*uf(i,j,k))
     $                /((h(i,j)+etf(i,j)+h(i-1,j)+etf(i-1,j))
     $                  *aru(i,j))
          end do
        end do
      end do
C
      return
C
      end
C
      subroutine advv
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Does horizontal and vertical advection of           *
C *                v-momentum, and includes coriolis, surface slope    *
C *                and baroclinic terms.                               *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      integer i,j,k
C
C     Do vertical advection:
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            vf(i,j,k)=0.e0
          end do
        end do
      end do
C
      do k=2,kbm1
        do j=2,jm
          do i=1,im
            vf(i,j,k)=.25e0*(w(i,j,k)+w(i,j-1,k))
     $                     *(v(i,j,k)+v(i,j,k-1))
          end do
        end do
      end do
C
C     Combine horizontal and vertical advection with coriolis, surface
C     slope and baroclinic terms:
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            vf(i,j,k)=advy(i,j,k)
     $                 +(vf(i,j,k)-vf(i,j,k+1))*arv(i,j)/dz(k)
     $                 +arv(i,j)*.25e0
     $                   *(cor(i,j)*dt(i,j)
     $                      *(u(i+1,j,k)+u(i,j,k))
     $                     +cor(i,j-1)*dt(i,j-1)
     $                       *(u(i+1,j-1,k)+u(i,j-1,k)))
     $                 +grav*.125e0*(dt(i,j)+dt(i,j-1))
     $                   *(egf(i,j)-egf(i,j-1)+egb(i,j)-egb(i,j-1)
     $                     +(e_atmos(i,j)-e_atmos(i,j-1))*2.e0)
     $                   *(dx(i,j)+dx(i,j-1))
     $                 +drhoy(i,j,k)
          end do
        end do
      end do
C
C     Step forward in time:
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            vf(i,j,k)=((h(i,j)+etb(i,j)+h(i,j-1)+etb(i,j-1))
     $                 *arv(i,j)*vb(i,j,k)
     $                 -2.e0*dti2*vf(i,j,k))
     $                /((h(i,j)+etf(i,j)+h(i,j-1)+etf(i,j-1))
     $                  *arv(i,j))
          end do
        end do
      end do
C
      return
C
      end
C
      subroutine areas_masks
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Calculates areas and masks.                         *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      integer i,j
C
C     Calculate areas of "t" and "s" cells:
C
      do j=1,jm
        do i=1,im
          art(i,j)=dx(i,j)*dy(i,j)
        end do
      end do
C
C     Calculate areas of "u" and "v" cells:
C
      do j=2,jm
        do i=2,im
          aru(i,j)=.25e0*(dx(i,j)+dx(i-1,j))*(dy(i,j)+dy(i-1,j))
          arv(i,j)=.25e0*(dx(i,j)+dx(i,j-1))*(dy(i,j)+dy(i,j-1))
        end do
      end do
C
      do j=1,jm
        aru(1,j)=aru(2,j)
        arv(1,j)=arv(2,j)
      end do
C
      do i=1,im
        aru(i,1)=aru(i,2)
        arv(i,1)=arv(i,2)
      end do
C
C     Initialise and set up free surface mask:
C
!      do j=1,jm
!        do i=1,im
!          fsm(i,j)=0.e0
!          dum(i,j)=0.e0
!          dvm(i,j)=0.e0
!          if(h(i,j).gt.1.e0) fsm(i,j)=1.e0
!        end do
!      end do
!C
!C     Set up velocity masks:
!C
!      do j=2,jm
!        do i=2,im
!          dum(i,j)=fsm(i,j)*fsm(i-1,j)
!          dvm(i,j)=fsm(i,j)*fsm(i,j-1)
!        end do
!      end do
!!      dum(2:im,1) = 1.     ! rwnd:
!!      dvm(1,2:jm) = 1.     !     :
!
!    I guess I was right about i=j=1 borders. The code below is from pom08-wad.
      do j=1,jm
        do i=1,im
          fsm(i,j) = 1.
          dum(i,j) = 1.
          dvm(i,j) = 1.
          if (h(i,j)<=1.) then
            h(i,j)   = 1.
            fsm(i,j) = 0.
            dum(i,j) = 0.
            dvm(i,j) = 0.
          end if
        end do
      end do
!    Close southwest corner.
!      fsm(1,1) = 0.
!
      do j=1,jmm1
        do i=1,im
          if (fsm(i,j)==0..and.fsm(i,j+1)/=0.) dvm(i,j+1) = 0.
        end do
      end do
      do j=1,jm
        do i=1,im-1
          if (fsm(i,j)==0..and.fsm(i+1,j)/=0.) dum(i+1,j) = 0.
        end do
      end do
!
C
      return
C
      end
C
      subroutine baropg
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Calculates  baroclinic pressure gradient.           *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      integer i,j,k
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            rho(i,j,k)=rho(i,j,k)-rmean(i,j,k)
          end do
        end do
      end do
C
C     Calculate x-component of baroclinic pressure gradient:
C
      do j=2,jmm1
        do i=2,imm1
          drhox(i,j,1)=.5e0*grav*(-zz(1))*(dt(i,j)+dt(i-1,j))
     $                  *(rho(i,j,1)-rho(i-1,j,1))
        end do
      end do
C
      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhox(i,j,k)=drhox(i,j,k-1)
     $                    +grav*.25e0*(zz(k-1)-zz(k))
     $                      *(dt(i,j)+dt(i-1,j))
     $                      *(rho(i,j,k)-rho(i-1,j,k)
     $                        +rho(i,j,k-1)-rho(i-1,j,k-1))
     $                    +grav*.25e0*(zz(k-1)+zz(k))
     $                      *(dt(i,j)-dt(i-1,j))
     $                      *(rho(i,j,k)+rho(i-1,j,k)
     $                        -rho(i,j,k-1)-rho(i-1,j,k-1))
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhox(i,j,k)=.25e0*(dt(i,j)+dt(i-1,j))
     $                        *drhox(i,j,k)*dum(i,j)
     $                        *(dy(i,j)+dy(i-1,j))
          end do
        end do
      end do
C
C     Calculate y-component of baroclinic pressure gradient:
C
      do j=2,jmm1
        do i=2,imm1
          drhoy(i,j,1)=.5e0*grav*(-zz(1))*(dt(i,j)+dt(i,j-1))
     $                  *(rho(i,j,1)-rho(i,j-1,1))
        end do
      end do
C
      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhoy(i,j,k)=drhoy(i,j,k-1)
     $                    +grav*.25e0*(zz(k-1)-zz(k))
     $                      *(dt(i,j)+dt(i,j-1))
     $                      *(rho(i,j,k)-rho(i,j-1,k)
     $                        +rho(i,j,k-1)-rho(i,j-1,k-1))
     $                    +grav*.25e0*(zz(k-1)+zz(k))
     $                      *(dt(i,j)-dt(i,j-1))
     $                      *(rho(i,j,k)+rho(i,j-1,k)
     $                        -rho(i,j,k-1)-rho(i,j-1,k-1))
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhoy(i,j,k)=.25e0*(dt(i,j)+dt(i,j-1))
     $                        *drhoy(i,j,k)*dvm(i,j)
     $                        *(dx(i,j)+dx(i,j-1))
          end do
        end do
      end do
C
      do k=1,kb
        do j=2,jmm1
          do i=2,imm1
            drhox(i,j,k)=ramp*drhox(i,j,k)
            drhoy(i,j,k)=ramp*drhoy(i,j,k)
          end do
        end do
      end do
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            rho(i,j,k)=rho(i,j,k)+rmean(i,j,k)
          end do
        end do
      end do
C
      return
C
      end
C
      subroutine baropg_mcc
!-----------------------------------------------------------------------
!
!  FUNCTION: Calculates  baroclinic pressure gradient (from sbPOM)
!            4th order correction terms, following McCalpin
!
!-----------------------------------------------------------------------
      implicit none
      include 'pom2k.c'
      integer i,j,k
      real d4(im,jm),ddx(im,jm),drho(im,jm,kb),rhou(im,jm,kb)

      do k=1,kb
        do j=1,jm
          do i=1,im
            rho(i,j,k)=rho(i,j,k)-rmean(i,j,k)
          end do
        end do
      end do

! compute terms correct to 4th order
      do i=1,im
        do j=1,jm
          ddx(i,j)=0.
          d4(i,j)=0.
        end do
      end do
      do k=1,kb
        do j=1,jm
          do i=1,im
            rhou(i,j,k)=0.
            drho(i,j,k)=0.
          end do
        end do
      end do

! compute DRHO, RHOU, DDX and D4
      do j=1,jm
        do i=2,im
          do k=1,kbm1
            drho(i,j,k)=(rho(i,j,k)-rho(i-1,j,k))*dum(i,j)
            rhou(i,j,k)=0.5*(rho(i,j,k)+rho(i-1,j,k))*dum(i,j)
          end do
          ddx(i,j)=(d(i,j)-d(i-1,j))*dum(i,j)
          d4(i,j)=.5*(d(i,j)+d(i-1,j))*dum(i,j)
        end do
      end do

        do j=1,jm
          do i=2,imm1
            do k=1,kbm1
              drho(i,j,k)=drho(i,j,k) - (1./24.)*
     $                    (dum(i+1,j)*(rho(i+1,j,k)-rho(i,j,k))-
     $                    2*(rho(i,j,k)-rho(i-1,j,k))+
     $                    dum(i-1,j)*(rho(i-1,j,k)-rho(i-2,j,k)))
              rhou(i,j,k)=rhou(i,j,k) + (1./16.)*
     $                    (dum(i+1,j)*(rho(i,j,k)-rho(i+1,j,k))+
     $                    dum(i-1,j)*(rho(i-1,j,k)-rho(i-2,j,k)))
            end do
            ddx(i,j)=ddx(i,j)-(1./24.)*
     $               (dum(i+1,j)*(d(i+1,j)-d(i,j))-
     $               2*(d(i,j)-d(i-1,j))+
     $               dum(i-1,j)*(d(i-1,j)-d(i-2,j)))
            d4(i,j)=d4(i,j)+(1./16.)*
     $              (dum(i+1,j)*(d(i,j)-d(i+1,j))+
     $              dum(i-1,j)*(d(i-1,j)-d(i-2,j)))
          end do
        end do

! calculate x-component of baroclinic pressure gradient
      do j=2,jmm1
        do i=2,imm1
          drhox(i,j,1)=grav*(-zz(1))*d4(i,j)*drho(i,j,1)
        end do
      end do

      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhox(i,j,k)=drhox(i,j,k-1)
     $                   +grav*0.5e0*dzz(k-1)*d4(i,j)
     $                   *(drho(i,j,k-1)+drho(i,j,k))
     $                   +grav*0.5e0*(zz(k-1)+zz(k))*ddx(i,j)
     $                   *(rhou(i,j,k)-rhou(i,j,k-1))
          end do
        end do
      end do

      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhox(i,j,k)=.25e0*(dt(i,j)+dt(i-1,j))
     $                        *drhox(i,j,k)*dum(i,j)
     $                        *(dy(i,j)+dy(i-1,j))
          end do
        end do
      end do

! compute terms correct to 4th order
      do i=1,im
        do j=1,jm
          ddx(i,j)=0.
          d4(i,j)=0.
        end do
      end do
      do k=1,kb
        do j=1,jm
          do i=1,im
            rhou(i,j,k)=0.
            drho(i,j,k)=0.
          end do
        end do
      end do

! compute DRHO, RHOU, DDX and D4
      do j=2,jm
        do i=1,im
          do k=1,kbm1
            drho(i,j,k)=(rho(i,j,k)-rho(i,j-1,k))*dvm(i,j)
            rhou(i,j,k)=.5*(rho(i,j,k)+rho(i,j-1,k))*dvm(i,j)
          end do
          ddx(i,j)=(d(i,j)-d(i,j-1))*dvm(i,j)
          d4(i,j)=.5*(d(i,j)+d(i,j-1))*dvm(i,j)
        end do
      end do

        do j=2,jmm1
          do i=1,im
            do k=1,kbm1
              drho(i,j,k)=drho(i,j,k)-(1./24.)*
     $                    (dvm(i,j+1)*(rho(i,j+1,k)-rho(i,j,k))-
     $                    2*(rho(i,j,k)-rho(i,j-1,k))+
     $                    dvm(i,j-1)*(rho(i,j-1,k)-rho(i,j-2,k)))
              rhou(i,j,k)=rhou(i,j,k)+(1./16.)*
     $                    (dvm(i,j+1)*(rho(i,j,k)-rho(i,j+1,k))+
     $                    dvm(i,j-1)*(rho(i,j-1,k)-rho(i,j-2,k)))
            end do
            ddx(i,j)=ddx(i,j)-(1./24)*
     $               (dvm(i,j+1)*(d(i,j+1)-d(i,j))-
     $               2*(d(i,j)-d(i,j-1))+
     $               dvm(i,j-1)*(d(i,j-1)-d(i,j-2)))
            d4(i,j)=d4(i,j)+(1./16.)*
     $              (dvm(i,j+1)*(d(i,j)-d(i,j+1))+
     $              dvm(i,j-1)*(d(i,j-1)-d(i,j-2)))
          end do
        end do

! calculate y-component of baroclinic pressure gradient
      do j=2,jmm1
        do i=2,imm1
          drhoy(i,j,1)=grav*(-zz(1))*d4(i,j)*drho(i,j,1)
        end do
      end do

      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhoy(i,j,k)=drhoy(i,j,k-1)
     $                   +grav*0.5e0*dzz(k-1)*d4(i,j)
     $                   *(drho(i,j,k-1)+drho(i,j,k))
     $                   +grav*0.5e0*(zz(k-1)+zz(k))*ddx(i,j)
     $                   *(rhou(i,j,k)-rhou(i,j,k-1))
          end do
        end do
      end do

      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhoy(i,j,k)=.25e0*(dt(i,j)+dt(i,j-1))
     $                        *drhoy(i,j,k)*dvm(i,j)
     $                        *(dx(i,j)+dx(i,j-1))
          end do
        end do
      end do

      do k=1,kb
        do j=2,jmm1
          do i=2,imm1
            drhox(i,j,k)=ramp*drhox(i,j,k)
            drhoy(i,j,k)=ramp*drhoy(i,j,k)
          end do
        end do
      end do

      do k=1,kb
        do j=1,jm
          do i=1,im
            rho(i,j,k)=rho(i,j,k)+rmean(i,j,k)
          end do
        end do
      end do

      return
      end

!
      subroutine bcond_orig(idx) ! orig
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Applies open boundary conditions.                   *
C *                                                                    *
C *                Closed boundary conditions are automatically        *
C *                enabled through specification of the masks, dum,    *
C *                dvm and fsm, in which case the open boundary        *
C *                conditions, included below, will be overwritten.    *
C *                                                                    *
C *                                The C-Grid                          *
C *                                **********                          *
C *                                                                    *
C *                The diagram is for the case where u and v are the   *
C *                primary boundary conditions together with t and     *
C *                s (co-located with el)                              *
C *                                                                    *
C *                All interpolations are centered in space except     *
C *                those at lateral open boundary where an upstream    *
C *                Horizontal locations of e(el), t and s (etc.) are   *
C *                coincident.                                         *
C *                                                                    *
C *                People not acquainted with sigma coordinates have   *
C *                often asked what kind of boundary condition is      *
C *                applied along closed horizontal boundaries.         *
C *                Although the issue is not as important as it might  *
C *                be  for z-level grids, a direct answer is "half-    *
C *                slip" which, of course, is between free slip and    *
C *                non-slip.                                           *
C
C
C East and West end points for the C-grid in POM.
C
C                      west
C
C           v(1,j+1)=0           v(2,j+1) . . . .
C
C     ----<---<----<-----
C     |                 |
C  u(1,j)   el(1,j)   u(2,j)=BC  el(2,j)   u(3,j) . . . .
C             |                   |
C             -----<----<----<-----
C
C           v(1,j)=0              v(2,j) . . . .
C
C                                                    east
C
C                              . . . .  v(im-1,j+1)           v(im,j+1)=0
C
C
C                 . . .  .  u(im-1,j)   el(im-1,j)  u(im,j)=BC  el(im,j)
C                                            |                   |
C                                            ----->----->---->----
C
C                              . . . .   v(im-1,j)             v(im,j)=0
C
C  Notes:
C    1. The suffixes, f  or af, have been deleted.
C    2. All variables NOT designated as boundary condition (=BC) or set to
C zero or obtained from an interior point are calculated points.
C    3. u(1,j) is never used but is obtained from the interior point for
C cosmetic output. Its counterpart, u(im+1,j), does not exist.
C    4. v=0 at i=1 and i=im are used as open inflow BC's unless specified
C otherwise.
C    5. The south and north extremal points are obtained from the above by
C permuting u to v, v to u, i to j and j to i.


C **********************************************************************

      implicit none
C
      include 'pom2k.c'
C
      integer idx
      real ga,u1,wm
      integer i,j,k
C
      if(idx.eq.1) then
C
C-----------------------------------------------------------------------
C
C     External (2-D) boundary conditions:
C
C     In this example, the governing boundary conditions are a radiation
C     condition on uaf in the east and in the west, and vaf in the north
C     and south. The tangential velocities are set to zero on both
C     boundaries. These are only one set of possibilities and may not
C     represent a choice which yields the most physically realistic
C     result.
C
C     Elevation (in this application, elevation is not a primary
C     boundary condition):
!     Уровень на границе (в данной реализации не является "подаваемым" условием);
C
!     Уровень внутренних граничных точек для красоты копируется во внешние.
        do j=1,jm
          elf(1,j)=elf(2,j)     ! (ю_внутр -> ю_внешн)
          elf(im,j)=elf(imm1,j) ! (с_внутр -> с_внешн)
        end do
C
        do i=1,im
          elf(i,1)=elf(i,2)     ! (з_внутр -> з_внешн)
          elf(i,jm)=elf(i,jmm1) ! (в_внутр -> в_внешн)
        end do
C
!     Наложение маски свободной поверхности
        do j=1,jm
          do i=1,im
            elf(i,j)=elf(i,j)*fsm(i,j)
          end do
        end do
C
        return  ! выход из процедуры
C
      else if(idx.eq.2) then
C
C     External (2-D) velocity:
!     Баротропные скорости на границе:
C
        do j=2,jmm1
C
C     East:
C
!     Условие излучения на восточной границе для баротропного режима
!     Нормальная компонента интеграла скорости = интеграл скорости на восточной границе
!     (если мы учитываем уровень на границе, то) +фазовая скорость(частота?)*(разница_в_уровне_между_внутренней_граничной_точкой_и_граничным_условием)
          uaf(im,j)=uabe(j)
     $               +rfe*sqrt(grav/h(imm1,j))
     $                         *(el(imm1,j)-ele(j))
          uaf(im,j)=ramp*uaf(im,j) ! Применяем инерционный фильтр, если установлен.
          vaf(im,j)=0.e0           ! Тангенциальная компонента баротропной скорости отсутствует (неизвестна).
C
C     West:
C
!     То же самое на западной границе
          uaf(2,j)=uabw(j)
     $              -rfw*sqrt(grav/h(2,j))
     $                        *(el(2,j)-elw(j))
          uaf(2,j)=ramp*uaf(2,j)
          uaf(1,j)=uaf(2,j)
          vaf(1,j)=0.e0
C
        end do
C
        do i=2,imm1
C
C     North:
C
!     Северная
          vaf(i,jm)=vabn(i)
     $               +rfn*sqrt(grav/h(i,jmm1))
     $                         *(el(i,jmm1)-eln(i))
          vaf(i,jm)=ramp*vaf(i,jm)
          uaf(i,jm)=0.e0
C
C     South:
C
!     Южная
          vaf(i,2)=vabs(i)
     $              -rfs*sqrt(grav/h(i,2))
     $                        *(el(i,2)-els(i))
          vaf(i,2)=ramp*vaf(i,2)
          vaf(i,1)=vaf(i,2)
          uaf(i,1)=0.e0 !uabs(i)    ! Изн. `0.e0`. Тангенциальная компонента интеграла скорости известна.
C
        end do
C
!     Наложение маски свободной поверхности
        do j=1,jm
          do i=1,im
            uaf(i,j)=uaf(i,j)*dum(i,j)
            vaf(i,j)=vaf(i,j)*dvm(i,j)
          end do
        end do
C
        return
C
      else if(idx.eq.3) then
C
C-----------------------------------------------------------------------
C
C     Internal (3-D) boundary conditions:
C
C     Velocity (radiation conditions; smoothing is used in the direction
C     tangential to the boundaries):
C
!     Бароклинные скорости на границе (выполняется сглаживание тангенциальной компоненты):
!
        do k=1,kbm1
          do j=2,jmm1
C
C     East:
C
!     Восточная
            ga=sqrt(h(im,j)/hmax)   ! Некоторый коэффициент (тем меньше, чем меньше глубина в точке)
!     Условие излучения на восточной границе для бароклинного режима
!     Нормальная компонента скорости = "взвешенное" среднее скоростей трёх "соседних" внутренних точек у восточной границы * коэф.
!                                     + "взвешенное" среднее скоростей трёх "соседних" внешних точек на восточной границе * (1-коэф.)
            uf(im,j,k)=ga*(.25e0*u(imm1,j-1,k)+.5e0*u(imm1,j,k)
     $                     +.25e0*u(imm1,j+1,k))
     $                  +(1.e0-ga)*(.25e0*u(im,j-1,k)+.5e0*u(im,j,k)
     $                    +.25e0*u(im,j+1,k))
            vf(im,j,k)=0.e0
C
C     West:
C
!     Западная
            ga=sqrt(h(1,j)/hmax)
            uf(2,j,k)=ga*(.25e0*u(3,j-1,k)+.5e0*u(3,j,k)
     $                    +.25e0*u(3,j+1,k))
     $                 +(1.e0-ga)*(.25e0*u(2,j-1,k)+.5e0*u(2,j,k)
     $                   +.25e0*u(2,j+1,k))
            uf(1,j,k)=uf(2,j,k)
            vf(1,j,k)=0.e0
          end do
        end do
C
        do k=1,kbm1
          do i=2,imm1
C
C     North:
C
            ga=sqrt(h(i,jm)/hmax)
            vf(i,jm,k)=ga*(.25e0*v(i-1,jmm1,k)+.5e0*v(i,jmm1,k)
     $                     +.25e0*v(i+1,jmm1,k))
     $                  +(1.e0-ga)*(.25e0*v(i-1,jm,k)+.5e0*v(i,jm,k)
     $                    +.25e0*v(i+1,jm,k))
            uf(i,jm,k)=0.e0
C
C     South:
C
            ga=sqrt(h(i,1)/hmax)
            vf(i,2,k)=ga*(.25e0*v(i-1,3,k)+.5e0*v(i,3,k)
     $                    +.25e0*v(i+1,3,k))
     $                  +(1.e0-ga)*(.25e0*v(i-1,2,k)+.5e0*v(i,2,k)
     $                    +.25e0*v(i+1,2,k))
!     $                 +(1.e0-ga)*(.25e0*vbs(i-1,k)+.5e0*vbs(i,k)   ! rwnd:
!     $                   +.25e0*vbs(i+1,k))                         ! rwnd:
            vf(i,1,k)=vf(i,2,k)
            uf(i,1,k)=0.e0  !ubs(i,k)   ! Изн. `0.e0`. Тангенциальная компонента скорости на границе известна.
          end do
        end do
C
!     Наложение маски свободной поверхности
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*dum(i,j)
              vf(i,j,k)=vf(i,j,k)*dvm(i,j)
            end do
          end do
        end do
C
        return
C
      else if(idx.eq.4) then
C
C     Temperature and salinity boundary conditions (using uf and vf,
C     respectively):
C
!     ГУ температуры и солёности:
        do k=1,kbm1
          do j=1,jm
C
C     East:
C
!     Расчёт условного параметра (2*u*(dt/dx)?)
            u1=2.e0*u(im,j,k)*dti/(dx(im,j)+dx(imm1,j))
!     Если скорость отрицательная, то вток - берём температуру и солёность из граничных условий, ...
            if(u1.le.0.e0) then
              uf(im,j,k)=t(im,j,k)-u1*(tbe(j,k)-t(im,j,k))
              vf(im,j,k)=s(im,j,k)-u1*(sbe(j,k)-s(im,j,k))
!     ... иначе отток - берём температуру и солёность из внутренней граничной точки.
            else
              uf(im,j,k)=t(im,j,k)-u1*(t(im,j,k)-t(imm1,j,k))
              vf(im,j,k)=s(im,j,k)-u1*(s(im,j,k)-s(imm1,j,k))
              ! И если это не поверхность и не "дно", то учитываем вертикальную скорость.
              if(k.ne.1.and.k.ne.kbm1) then
                wm=.5e0*(w(imm1,j,k)+w(imm1,j,k+1))*dti
     $              /((zz(k-1)-zz(k+1))*dt(imm1,j))
                uf(im,j,k)=uf(im,j,k)-wm*(t(imm1,j,k-1)-t(imm1,j,k+1))
                vf(im,j,k)=vf(im,j,k)-wm*(s(imm1,j,k-1)-s(imm1,j,k+1))
              endif
            endif
C
C     West:
C
!     Западная
            u1=2.e0*u(2,j,k)*dti/(dx(1,j)+dx(2,j))
            if(u1.ge.0.e0) then
              uf(1,j,k)=t(1,j,k)-u1*(t(1,j,k)-tbw(j,k))
              vf(1,j,k)=s(1,j,k)-u1*(s(1,j,k)-sbw(j,k))
            else
              uf(1,j,k)=t(1,j,k)-u1*(t(2,j,k)-t(1,j,k))
              vf(1,j,k)=s(1,j,k)-u1*(s(2,j,k)-s(1,j,k))
              if(k.ne.1.and.k.ne.kbm1) then
                wm=.5e0*(w(2,j,k)+w(2,j,k+1))*dti
     $              /((zz(k-1)-zz(k+1))*dt(2,j))
                uf(1,j,k)=uf(1,j,k)-wm*(t(2,j,k-1)-t(2,j,k+1))
                vf(1,j,k)=vf(1,j,k)-wm*(s(2,j,k-1)-s(2,j,k+1))
              endif
            endif
          end do
        end do
C
        do k=1,kbm1
          do i=1,im
C
C     North:
C
!     Северная
            u1=2.e0*v(i,jm,k)*dti/(dy(i,jm)+dy(i,jmm1))
            if(u1.le.0.e0) then
              uf(i,jm,k)=t(i,jm,k)-u1*(tbn(i,k)-t(i,jm,k))
              vf(i,jm,k)=s(i,jm,k)-u1*(sbn(i,k)-s(i,jm,k))
            else
              uf(i,jm,k)=t(i,jm,k)-u1*(t(i,jm,k)-t(i,jmm1,k))
              vf(i,jm,k)=s(i,jm,k)-u1*(s(i,jm,k)-s(i,jmm1,k))
              if(k.ne.1.and.k.ne.kbm1) then
                wm=.5e0*(w(i,jmm1,k)+w(i,jmm1,k+1))*dti
     $              /((zz(k-1)-zz(k+1))*dt(i,jmm1))
                uf(i,jm,k)=uf(i,jm,k)-wm*(t(i,jmm1,k-1)-t(i,jmm1,k+1))
                vf(i,jm,k)=vf(i,jm,k)-wm*(s(i,jmm1,k-1)-s(i,jmm1,k+1))
              endif
            endif
C
C     South:
C
!     Южная
            u1=2.e0*v(i,2,k)*dti/(dy(i,1)+dy(i,2))
            if(u1.ge.0.e0) then
              uf(i,1,k)=t(i,1,k)-u1*(t(i,1,k)-tbs(i,k))
              vf(i,1,k)=s(i,1,k)-u1*(s(i,1,k)-sbs(i,k))
            else
              uf(i,1,k)=t(i,1,k)-u1*(t(i,2,k)-t(i,1,k))
              vf(i,1,k)=s(i,1,k)-u1*(s(i,2,k)-s(i,1,k))
              if(k.ne.1.and.k.ne.kbm1) then
                wm=.5e0*(w(i,2,k)+w(i,2,k+1))*dti
     $              /((zz(k-1)-zz(k+1))*dt(i,2))
                uf(i,1,k)=uf(i,1,k)-wm*(t(i,2,k-1)-t(i,2,k+1))
                vf(i,1,k)=vf(i,1,k)-wm*(s(i,2,k-1)-s(i,2,k+1))
              endif
            endif
          end do
        end do
C
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)
            end do
          end do
        end do
C
        return
C
      else if(idx.eq.5) then
C
C     Vertical velocity boundary conditions:
C
!     Условие вертикальной скорости на границе:
!     Отсутствует. Просто накладывается маска свободной поверхности.
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              w(i,j,k)=w(i,j,k)*fsm(i,j)
            end do
          end do
        end do
C
        return
C
      else if(idx.eq.6) then
C
C     q2 and q2l boundary conditions:
C
!     ГУ турбулентной кинетической энергии и турб. масштаба:
        do k=1,kb
          do j=1,jm
C
C     East:
C
            u1=2.e0*u(im,j,k)*dti/(dx(im,j)+dx(imm1,j))
!     Если вток, то используем в качестве граничного, значение близкое к нулю.
            if(u1.le.0.e0) then
              uf(im,j,k)=q2(im,j,k)-u1*(small-q2(im,j,k))
              vf(im,j,k)=q2l(im,j,k)-u1*(small-q2l(im,j,k))
            else
!     Если отток - для расчёта используются внутренние точки.
              uf(im,j,k)=q2(im,j,k)-u1*(q2(im,j,k)-q2(imm1,j,k))
              vf(im,j,k)=q2l(im,j,k)-u1*(q2l(im,j,k)-q2l(imm1,j,k))
            endif
C
C     West:
C
            u1=2.e0*u(2,j,k)*dti/(dx(1,j)+dx(2,j))
            if(u1.ge.0.e0) then
              uf(1,j,k)=q2(1,j,k)-u1*(q2(1,j,k)-small)
              vf(1,j,k)=q2l(1,j,k)-u1*(q2l(1,j,k)-small)
            else
              uf(1,j,k)=q2(1,j,k)-u1*(q2(2,j,k)-q2(1,j,k))
              vf(1,j,k)=q2l(1,j,k)-u1*(q2l(2,j,k)-q2l(1,j,k))
            endif
          end do
        end do
C
        do k=1,kb
          do i=1,im
C
C     North:
C
            u1=2.e0*v(i,jm,k)*dti/(dy(i,jm)+dy(i,jmm1))
            if(u1.le.0.e0) then
              uf(i,jm,k)=q2(i,jm,k)-u1*(small-q2(i,jm,k))
              vf(i,jm,k)=q2l(i,jm,k)-u1*(small-q2l(i,jm,k))
            else
              uf(i,jm,k)=q2(i,jm,k)-u1*(q2(i,jm,k)-q2(i,jmm1,k))
              vf(i,jm,k)=q2l(i,jm,k)-u1*(q2l(i,jm,k)-q2l(i,jmm1,k))
            endif
C
C     South:
C
            u1=2.e0*v(i,2,k)*dti/(dy(i,1)+dy(i,2))
            if(u1.ge.0.e0) then
              uf(i,1,k)=q2(i,1,k)-u1*(q2(i,1,k)-small)
              vf(i,1,k)=q2l(i,1,k)-u1*(q2l(i,1,k)-small)
            else
              uf(i,1,k)=q2(i,1,k)-u1*(q2(i,2,k)-q2(i,1,k))
              vf(i,1,k)=q2l(i,1,k)-u1*(q2l(i,2,k)-q2l(i,1,k))
            endif
          end do
        end do
C
!     Наложение маски свободной поверхности с предотвращением получения нулевых значений.
        do k=1,kb
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)+1.e-10
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)+1.e-10
            end do
          end do
        end do
C
        return
C
      endif
C
      end
C
      subroutine bcond_inflow(idx)
C **********************************************************************
C
C **********************************************************************

      implicit none
C
      include 'pom2k.c'
C
      integer idx
      real ga,u1,wm
      integer i,j,k
C
      if(idx.eq.1) then
C
C-----------------------------------------------------------------------
C
C     External (2-D) boundary conditions:
C
C     Elevation (in this application, elevation is not a primary
C     boundary condition):
!     Уровень на границе (в данной реализации не является "подаваемым" условием);
!
!     Inflow condition for elevation.
        do j=1,jm
          elf(1,j)=elf(2,j)     ! (ю_внутр -> ю_внешн)
          elf(im,j)=elf(imm1,j) ! (с_внутр -> с_внешн)
        end do
!
        do i=1,im
          elf(i,1)=elf(i,2)     ! (з_внутр -> з_внешн)
          elf(i,jm)=elf(i,jmm1) ! (в_внутр -> в_внешн)
        end do
C
!     Наложение маски свободной поверхности
        do j=1,jm
          do i=1,im
            elf(i,j)=elf(i,j)*fsm(i,j)
          end do
        end do
C
        return  ! выход из процедуры
C
      else if(idx.eq.2) then
C
C     External (2-D) velocity:
!     Баротропные скорости на границе:
C
        do j=2,jmm1
C
C     East:
C
!     Inflow condition for barotropic velocities
          uaf(im,j)=2*(uabe(j)-sqrt(grav/h(imm1,j))*ele(j))
     $              /(h(im,j)+elf(im,j)+h(imm1,j)+elf(imm1,j))
          uaf(im,j)=ramp*uaf(im,j) ! Применяем инерционный фильтр, если установлен. Nedded...? or not?
          vaf(im,j)=0.e0           ! (set) Тангенциальная компонента баротропной скорости отсутствует (неизвестна).
C
C     West:
C
!     То же самое на западной границе
          uaf(2,j)=2*(uabw(j)-sqrt(grav/h(2,j))*elw(j))
     $             /(h(1,j)+elf(1,j)+h(2,j)+elf(2,j))
          uaf(2,j)=ramp*uaf(2,j)
          uaf(1,j)=uaf(2,j)
          vaf(1,j)=0.e0             ! (set)
C
        end do
C
        do i=2,imm1
C
C     North:
C
!     Северная
          vaf(i,jm)=2*(vabn(i)-sqrt(grav/h(i,jmm1))*eln(i))
     $              /(h(i,jm)+elf(i,jm)+h(i,jmm1)+elf(i,jmm1))
          vaf(i,jm)=ramp*vaf(i,jm)
          uaf(i,jm)=0.e0            ! (set)
C
C     South:
C
!     Южная
          vaf(i,2)=2*(vabs(i)-sqrt(grav/h(i,2))*eln(i))
     $             /(h(i,1)+elf(i,1)+h(i,2)+elf(i,2))
          vaf(i,2)=ramp*vaf(i,2)
          vaf(i,1)=vaf(i,2)
          uaf(i,1)=uabs(i)    ! Изн. `0.e0`. Тангенциальная компонента интеграла скорости известна.
C
        end do
C
!     Наложение маски свободной поверхности
        do j=1,jm
          do i=1,im
            uaf(i,j)=uaf(i,j)*dum(i,j)
            vaf(i,j)=vaf(i,j)*dvm(i,j)
          end do
        end do
C
        return
C
      else if(idx.eq.3) then
C
C-----------------------------------------------------------------------
C
C     Internal (3-D) boundary conditions:
!     Бароклинные скорости на границе:
!
        do k=1,kbm1
          do j=2,jmm1
C
C     East:
C
!     Восточная
            uf(im,j,k)=0.e0             ! (BC: unknown)
            vf(im,j,k)=0.e0             ! (set)
C
C     West:
C
!     Западная
            uf(2,j,k)=0.e0              ! (BC: unknown)
            uf(1,j,k)=uf(2,j,k)
            vf(1,j,k)=0.e0              ! (set)
          end do
        end do
C
        do k=1,kbm1
          do i=2,imm1
C
C     North:
C
            vf(i,jm,k)=0.e0             ! (BC: unknown)
            uf(i,jm,k)=0.e0             ! (set)
C
C     South:
C
            vf(i,2,k)=vbs(i,k)          ! (BC)
            vf(i,1,k)=vf(i,2,k)
            uf(i,1,k)=ubs(i,k)   ! Изн. `0.e0`. Тангенциальная компонента скорости на границе известна.
          end do
        end do
C
!     Наложение маски свободной поверхности
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*dum(i,j)
              vf(i,j,k)=vf(i,j,k)*dvm(i,j)
            end do
          end do
        end do
C
        return
C
      else if(idx.eq.4) then
C
C     Temperature and salinity boundary conditions (using uf and vf,
C     respectively):
C
!     ГУ температуры и солёности:
        do k=1,kbm1
          do j=1,jm
C
C     East:
C
!     Расчёт условного параметра (2*u*(dt/dx)?)
            u1=2.e0*u(im,j,k)*dti/(dx(im,j)+dx(imm1,j))
!     Если скорость отрицательная, то вток - берём температуру и солёность из граничных условий, ...
            if(u1.le.0.e0) then
              uf(im,j,k)=t(im,j,k)-u1*(tbe(j,k)-t(im,j,k))
              vf(im,j,k)=s(im,j,k)-u1*(sbe(j,k)-s(im,j,k))
!     ... иначе отток - берём температуру и солёность из внутренней граничной точки.
            else
              uf(im,j,k)=t(im,j,k)-u1*(t(im,j,k)-t(imm1,j,k))
              vf(im,j,k)=s(im,j,k)-u1*(s(im,j,k)-s(imm1,j,k))
              ! И если это не поверхность и не "дно", то учитываем вертикальную скорость.
              if(k.ne.1.and.k.ne.kbm1) then
                wm=.5e0*(w(imm1,j,k)+w(imm1,j,k+1))*dti
     $              /((zz(k-1)-zz(k+1))*dt(imm1,j))
                uf(im,j,k)=uf(im,j,k)-wm*(t(imm1,j,k-1)-t(imm1,j,k+1))
                vf(im,j,k)=vf(im,j,k)-wm*(s(imm1,j,k-1)-s(imm1,j,k+1))
              endif
            endif
C
C     West:
C
!     Западная
            u1=2.e0*u(2,j,k)*dti/(dx(1,j)+dx(2,j))
            if(u1.ge.0.e0) then
              uf(1,j,k)=t(1,j,k)-u1*(t(1,j,k)-tbw(j,k))
              vf(1,j,k)=s(1,j,k)-u1*(s(1,j,k)-sbw(j,k))
            else
              uf(1,j,k)=t(1,j,k)-u1*(t(2,j,k)-t(1,j,k))
              vf(1,j,k)=s(1,j,k)-u1*(s(2,j,k)-s(1,j,k))
              if(k.ne.1.and.k.ne.kbm1) then
                wm=.5e0*(w(2,j,k)+w(2,j,k+1))*dti
     $              /((zz(k-1)-zz(k+1))*dt(2,j))
                uf(1,j,k)=uf(1,j,k)-wm*(t(2,j,k-1)-t(2,j,k+1))
                vf(1,j,k)=vf(1,j,k)-wm*(s(2,j,k-1)-s(2,j,k+1))
              endif
            endif
          end do
        end do
C
        do k=1,kbm1
          do i=1,im
C
C     North:
C
!     Северная
            u1=2.e0*v(i,jm,k)*dti/(dy(i,jm)+dy(i,jmm1))
            if(u1.le.0.e0) then
              uf(i,jm,k)=t(i,jm,k)-u1*(tbn(i,k)-t(i,jm,k))
              vf(i,jm,k)=s(i,jm,k)-u1*(sbn(i,k)-s(i,jm,k))
            else
              uf(i,jm,k)=t(i,jm,k)-u1*(t(i,jm,k)-t(i,jmm1,k))
              vf(i,jm,k)=s(i,jm,k)-u1*(s(i,jm,k)-s(i,jmm1,k))
              if(k.ne.1.and.k.ne.kbm1) then
                wm=.5e0*(w(i,jmm1,k)+w(i,jmm1,k+1))*dti
     $              /((zz(k-1)-zz(k+1))*dt(i,jmm1))
                uf(i,jm,k)=uf(i,jm,k)-wm*(t(i,jmm1,k-1)-t(i,jmm1,k+1))
                vf(i,jm,k)=vf(i,jm,k)-wm*(s(i,jmm1,k-1)-s(i,jmm1,k+1))
              endif
            endif
C
C     South:
C
!     Южная
            u1=2.e0*v(i,2,k)*dti/(dy(i,1)+dy(i,2))
            if(u1.ge.0.e0) then
              uf(i,1,k)=t(i,1,k)-u1*(t(i,1,k)-tbs(i,k))
              vf(i,1,k)=s(i,1,k)-u1*(s(i,1,k)-sbs(i,k))
            else
              uf(i,1,k)=t(i,1,k)-u1*(t(i,2,k)-t(i,1,k))
              vf(i,1,k)=s(i,1,k)-u1*(s(i,2,k)-s(i,1,k))
              if(k.ne.1.and.k.ne.kbm1) then
                wm=.5e0*(w(i,2,k)+w(i,2,k+1))*dti
     $              /((zz(k-1)-zz(k+1))*dt(i,2))
                uf(i,1,k)=uf(i,1,k)-wm*(t(i,2,k-1)-t(i,2,k+1))
                vf(i,1,k)=vf(i,1,k)-wm*(s(i,2,k-1)-s(i,2,k+1))
              endif
            endif
          end do
        end do
C
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)
            end do
          end do
        end do
C
        return
C
      else if(idx.eq.5) then
C
C     Vertical velocity boundary conditions:
C
!     Условие вертикальной скорости на границе:
!     Отсутствует. Просто накладывается маска свободной поверхности.
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              w(i,j,k)=w(i,j,k)*fsm(i,j)
            end do
          end do
        end do
C
        return
C
      else if(idx.eq.6) then
C
C     q2 and q2l boundary conditions:
C
!     ГУ турбулентной кинетической энергии и турб. масштаба:
        do k=1,kb
          do j=1,jm
C
C     East:
C
            u1=2.e0*u(im,j,k)*dti/(dx(im,j)+dx(imm1,j))
!     Если вток, то используем в качестве граничного, значение близкое к нулю.
            if(u1.le.0.e0) then
              uf(im,j,k)=q2(im,j,k)-u1*(small-q2(im,j,k))
              vf(im,j,k)=q2l(im,j,k)-u1*(small-q2l(im,j,k))
            else
!     Если отток - для расчёта используются внутренние точки.
              uf(im,j,k)=q2(im,j,k)-u1*(q2(im,j,k)-q2(imm1,j,k))
              vf(im,j,k)=q2l(im,j,k)-u1*(q2l(im,j,k)-q2l(imm1,j,k))
            endif
C
C     West:
C
            u1=2.e0*u(2,j,k)*dti/(dx(1,j)+dx(2,j))
            if(u1.ge.0.e0) then
              uf(1,j,k)=q2(1,j,k)-u1*(q2(1,j,k)-small)
              vf(1,j,k)=q2l(1,j,k)-u1*(q2l(1,j,k)-small)
            else
              uf(1,j,k)=q2(1,j,k)-u1*(q2(2,j,k)-q2(1,j,k))
              vf(1,j,k)=q2l(1,j,k)-u1*(q2l(2,j,k)-q2l(1,j,k))
            endif
          end do
        end do
C
        do k=1,kb
          do i=1,im
C
C     North:
C
            u1=2.e0*v(i,jm,k)*dti/(dy(i,jm)+dy(i,jmm1))
            if(u1.le.0.e0) then
              uf(i,jm,k)=q2(i,jm,k)-u1*(small-q2(i,jm,k))
              vf(i,jm,k)=q2l(i,jm,k)-u1*(small-q2l(i,jm,k))
            else
              uf(i,jm,k)=q2(i,jm,k)-u1*(q2(i,jm,k)-q2(i,jmm1,k))
              vf(i,jm,k)=q2l(i,jm,k)-u1*(q2l(i,jm,k)-q2l(i,jmm1,k))
            endif
C
C     South:
C
            u1=2.e0*v(i,2,k)*dti/(dy(i,1)+dy(i,2))
            if(u1.ge.0.e0) then
              uf(i,1,k)=q2(i,1,k)-u1*(q2(i,1,k)-small)
              vf(i,1,k)=q2l(i,1,k)-u1*(q2l(i,1,k)-small)
            else
              uf(i,1,k)=q2(i,1,k)-u1*(q2(i,2,k)-q2(i,1,k))
              vf(i,1,k)=q2l(i,1,k)-u1*(q2l(i,2,k)-q2l(i,1,k))
            endif
          end do
        end do
C
!     Наложение маски свободной поверхности с предотвращением получения нулевых значений.
        do k=1,kb
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)+1.e-10
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)+1.e-10
            end do
          end do
        end do
C
        return
C
      endif
C
      end
C
      subroutine bcond_sbc(idx)
C **********************************************************************
C   Stupid BCs
C **********************************************************************
!
      implicit none
C
      include 'pom2k.c'
C
      integer idx
      real ga,u1,wm
      integer i,j,k
C
      if(idx.eq.1) then
C
C-----------------------------------------------------------------------
C
C     External (2-D) boundary conditions:
C
        do j=1,jm
          elf( 1,j) = elf(  2 ,j)
          elf(im,j) = elf(imm1,j)
        end do
        do i=1,im
          elf(i, 1) = elf(i,  2 )
          elf(i,jm) = elf(i,jmm1)
        end do
!     Наложение маски свободной поверхности
        do j=1,jm
          do i=1,im
            elf(i,j)=elf(i,j)*fsm(i,j)
          end do
        end do
!
        return  ! выход из процедуры
!
      else if(idx.eq.2) then
!
!     External (2-D) velocity:
!     Баротропные скорости на границе:
!
        do j=1,jm
!     East:
          uaf(im,j) = 2.*uabe(j)
     $               /(h(im,j)+elf(im,j)+h(imm1,j)+elf(imm1,j))
          vaf(im,j) = va(im,j)+ramp*(vabe(j)-va(im,j))
!     West:
          uaf(2,j) = 2.*uabw(j)
     $               /(h(1,j)+elf(1,j)+h(2,j)+elf(2,j))
          vaf(1,j) = va(2,j)+ramp*(vabw(j)-va(2,j))
          uaf(1,j) = uaf(2,j)
        end do
!
        do i=2,imm1
!     North:
          vaf(i,jm) = 2.*vabn(i)
     $               /(h(i,jm)+elf(i,jm)+h(i,jmm1)+elf(i,jmm1))
          uaf(i,jm) = ua(i,jm)+ramp*(uabn(i)-ua(i,jm))
!     South:
          vaf(i,2) = 2.*vabs(i)
     $               /(h(i,1)+elf(i,1)+h(i,2)+elf(i,2))
          uaf(i,1) = ua(i,2)+ramp*(uabs(i)-ua(i,2))
          vaf(i,1) = vaf(i,2)
        end do
C
!     Наложение маски свободной поверхности
        do j=1,jm
          do i=1,im
            uaf(i,j)=uaf(i,j)*dum(i,j)
            vaf(i,j)=vaf(i,j)*dvm(i,j)
          end do
        end do
C
        return
C
      else if(idx.eq.3) then
C
C-----------------------------------------------------------------------
C
C     Internal (3-D) boundary conditions:
!     Бароклинные скорости на границе:
!
        do k=1,kbm1
          do j=1,jm
!     East:
            uf(im,j,k) = u(im,j,k)+ramp*(ube(j,k)-u(im,j,k))
            vf(im,j,k) = v(im,j,k)+ramp*(vbe(j,k)-v(im,j,k))
!     West:
            uf(2,j,k) = u(2,j,k)+ramp*(ubw(j,k)-u(2,j,k))
            vf(1,j,k) = v(2,j,k)+ramp*(vbw(j,k)-v(2,j,k))
            uf(1,j,k) = uf(2,j,k)
          end do
        end do
!
        do k=1,kbm1
          do i=1,im
!     North:
            vf(i,jm,k) = v(i,jm,k)+ramp*(vbn(i,k)-v(i,jm,k))
            uf(i,jm,k) = u(i,jm,k)+ramp*(ubn(i,k)-u(i,jm,k))
!     South:
            vf(i,2,k) = v(i,2,k)+ramp*(vbs(i,k)-v(i,2,k))
            uf(i,1,k) = u(i,2,k)+ramp*(ubs(i,k)-u(i,2,k))
            vf(i,1,k) = vf(i,2,k)
          end do
        end do
!     Наложение маски свободной поверхности
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*dum(i,j)
              vf(i,j,k)=vf(i,j,k)*dvm(i,j)
            end do
          end do
        end do
!
        return
!
      else if(idx.eq.4) then
!
!     Temperature and salinity boundary conditions (using uf and vf,
!     respectively):
!
!     ГУ температуры и солёности:
        do k=1,kbm1
          do j=1,jm
!
!     East:
!
            uf(im,j,k) = t(im,j,k)-dti/(dx(im,j)+dx(imm1,j))
     $                  *(
     $                    (u(im,j,k)+abs(u(im,j,k)))
     $                    *(t(im,j,k)-t(imm1,j,k))
     $                   +(u(im,j,k)-abs(u(im,j,k)))
     $                    *(tbe(j,k)-t(im,j,k))
     $                   )
            vf(im,j,k) = s(im,j,k)-dti/(dx(im,j)+dx(imm1,j))
     $                  *(
     $                    (u(im,j,k)+abs(u(im,j,k)))
     $                    *(s(im,j,k)-s(imm1,j,k))
     $                   +(u(im,j,k)-abs(u(im,j,k)))
     $                    *(sbe(j,k)-s(im,j,k))
     $                   )
!     West:
            uf(1,j,k) = t(1,j,k)-dti/(dx(1,j)+dx(2,j))
     $                  *(
     $                    (u(2,j,k)+abs(u(2,j,k)))
     $                    *(t(1,j,k)-tbw(j,k))
     $                   +(u(2,j,k)-abs(u(2,j,k)))
     $                    *(t(2,j,k)-t(1,j,k))
     $                   )
            vf(1,j,k) = s(1,j,k)-dti/(dx(1,j)+dx(2,j))
     $                  *(
     $                    (u(2,j,k)+abs(u(2,j,k)))
     $                    *(s(1,j,k)-sbw(j,k))
     $                   +(u(2,j,k)-abs(u(2,j,k)))
     $                    *(s(2,j,k)-s(1,j,k))
     $                   )
          end do
        end do
C
        do k=1,kbm1
          do i=1,im
!     North:
            uf(i,jm,k) = t(i,jm,k)-dti/(dy(i,jm)+dy(i,jmm1))
     $                  *(
     $                    (v(i,jm,k)+abs(v(i,jm,k)))
     $                    *(t(i,jm,k)-t(i,jmm1,k))
     $                   +(v(i,jm,k)-abs(v(i,jm,k)))
     $                    *(tbn(i,k)-t(i,jm,k))
     $                   )
            vf(i,jm,k) = s(i,jm,k)-dti/(dy(i,jm)+dy(i,jmm1))
     $                  *(
     $                    (v(i,jm,k)+abs(v(i,jm,k)))
     $                    *(s(i,jm,k)-s(i,jmm1,k))
     $                   +(v(i,jm,k)-abs(v(i,jm,k)))
     $                    *(sbn(i,k)-s(i,jm,k))
     $                   )
!     South:
            uf(i,1,k) = t(i,1,k)-dti/(dy(i,1)+dy(i,2))
     $                  *(
     $                    (v(i,2,k)+abs(v(i,2,k)))
     $                    *(t(i,1,k)-tbs(i,k))
     $                   +(v(i,2,k)-abs(v(i,2,k)))
     $                    *(t(i,2,k)-t(i,1,k))
     $                   )
            vf(i,1,k) = s(i,1,k)-dti/(dy(i,1)+dy(i,2))
     $                  *(
     $                    (v(i,2,k)+abs(v(i,2,k)))
     $                    *(s(i,1,k)-sbs(i,k))
     $                   +(v(i,2,k)-abs(v(i,2,k)))
     $                    *(s(i,2,k)-s(i,1,k))
     $                   )
          end do
        end do
!
        return
C
      else if(idx.eq.5) then
C
C     Vertical velocity boundary conditions:
C
!     Условие вертикальной скорости на границе:
!     Отсутствует. Просто накладывается маска свободной поверхности.
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              w(i,j,k)=w(i,j,k)*fsm(i,j)
            end do
          end do
        end do
C
        return
C
      else if(idx.eq.6) then
C
C     q2 and q2l boundary conditions:
C
!     ГУ турбулентной кинетической энергии и турб. масштаба:
        do k=1,kb
          do j=1,jm
C
C     East:
C
            uf(1,j,k) = 1.e-10
            vf(1,j,k) = 1.e-10
C
C     West:
C
            uf(im,j,k) = 1.e-10
            vf(im,j,k) = 1.e-10
!
          end do
        end do
C
        do k=1,kb
          do i=1,im
C
C     North:
C
            uf(i,1,k) = 1.e-10
            vf(i,1,k) = 1.e-10
C
C     South:
C
            uf(i,jm,k) = 1.e-10
            vf(i,jm,k) = 1.e-10
!
          end do
        end do
C
!     Наложение маски свободной поверхности с предотвращением получения нулевых значений.
!       Left from original pom2k
        do k=1,kb
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)+1.e-10
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)+1.e-10
            end do
          end do
        end do
C
        return
C
      endif
C
      end
!
      subroutine bcond_clamped(idx)
C **********************************************************************
C   My clamped BCs
C **********************************************************************

      implicit none
C
      include 'pom2k.c'
C
      integer idx
      real ga,u1,wm
      integer i,j,k
C
      if(idx.eq.1) then
C
C-----------------------------------------------------------------------
C
C     External (2-D) boundary conditions:
C
        do j=2,jmm1
!     West:
          elf(2,j)    = el(2,j)+ramp*(elw(j)-el(2,j))
          elf(1,j)    = elf(2,j)
!     East:
          elf(imm1,j) = el(imm1,j)+ramp*(ele(j)-el(imm1,j))
          elf( im ,j) = elf(imm1,j)
        end do
!
        do i=2,imm1
!     South:
          elf(i,2)    = el(i,2)+ramp*(els(i)-el(i,2))
          elf(i,1)    = elf(i,2)
!     North:
          elf(i,jmm1) = el(i,jmm1)+ramp*(eln(i)-el(i,jmm1))
          elf(i, jm ) = elf(i,jmm1)
        end do
!
        elf(1, 1) = (elf(1,  2 )+el(2, 1))*.5   ! Why elf+el? And not elf+elf?
        elf(1,jm) = (elf(1,jmm1)+el(2,jm))*.5
!       Two next lines were not present in RaSBCs.
        elf(im, 1) = (elf(im,  2 )+el(imm1, 1))*.5
        elf(im,jm) = (elf(im,jmm1)+el(imm1,jm))*.5
!
!        do j=1,jm
!          elf(1,j)=elf(2,j)     ! (з_внутр -> з_внешн)
!          elf(im,j)=elf(imm1,j) ! (в_внутр -> в_внешн)
!        end do
!
!        do i=1,im
!          elf(i,1)=elf(i,2)     ! (ю_внутр -> ю_внешн)
!          elf(i,jm)=elf(i,jmm1) ! (с_внутр -> с_внешн)
!        end do
!
C
!     Наложение маски свободной поверхности
        do j=1,jm
          do i=1,im
            elf(i,j)=elf(i,j)*fsm(i,j)
          end do
        end do
C
        return  ! выход из процедуры
C
      else if(idx.eq.2) then
C
C     External (2-D) velocity:
!     Баротропные скорости на границе:
!
        do j=2,jmm1
!
!     East:
!
          uaf(im,j) = ua(im,j)+ramp*(uabe(j)-ua(im,j))
          vaf(im,j) = ramp*vabe(j)
!     West:
!
          uaf(2,j) = ramp*uabw(j)
          vaf(1,j) = ramp*vabw(j)
          uaf(1,j) = uaf(2,j)
!
        end do
!
        do i=2,imm1
!
!     North:
!
          vaf(i,jm) = ramp*vabn(i)
          uaf(i,jm) = ramp*uabn(i)
!
!     South:
!
          vaf(i,2) = ramp*vabs(i)
          uaf(i,1) = ramp*uabs(i)
          vaf(i,1) = vaf(i,2)
!
        end do
C
!     Наложение маски свободной поверхности
        do j=1,jm
          do i=1,im
            uaf(i,j)=uaf(i,j)*dum(i,j)
            vaf(i,j)=vaf(i,j)*dvm(i,j)
          end do
        end do
C
        return
C
      else if(idx.eq.3) then
C
C-----------------------------------------------------------------------
C
C     Internal (3-D) boundary conditions:
!     Бароклинные скорости на границе:
!
        do k=1,kbm1
          do j=2,jmm1
!
!     East:
!
            uf(im,j,k) = ramp*ube(j,k)
            vf(im,j,k) = ramp*vbe(j,k)
!
!     West:
!
            uf(2,j,k) = ramp*ubw(j,k)
            vf(1,j,k) = ramp*vbw(j,k)
            uf(1,j,k) = uf(2,j,k)
!
          end do
        end do
!
        do k=1,kbm1
          do i=2,imm1
!
!     North:
!
            vf(i,jm,k) = ramp*vbn(i,k)
            uf(i,jm,k) = ramp*ubn(i,k)
!
!     South:
!
            vf(i,2,k) = ramp*vbs(i,k)
            uf(i,1,k) = ramp*ubs(i,k)
            vf(i,1,k) = vf(i,2,k)
!
          end do
        end do
C
!     Наложение маски свободной поверхности
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*dum(i,j)
              vf(i,j,k)=vf(i,j,k)*dvm(i,j)
            end do
          end do
        end do
C
        return
C
      else if(idx.eq.4) then
!
!     Temperature and salinity boundary conditions (using uf and vf,
!     respectively):
!
!     ГУ температуры и солёности:
        do k=1,kbm1
          do j=2,jmm1
!
!     East:
!
            uf(im,j,k) = ramp*tbe(j,k)
            vf(im,j,k) = ramp*sbs(j,k)
!
!     West:
!
            uf(1,j,k) = ramp*tbw(j,k)
            vf(1,j,k) = ramp*sbw(j,k)
!
          end do
        end do
C
        do k=1,kbm1
          do i=2,imm1
!
!     North:
!
            uf(i,jm,k) = ramp*tbn(i,k)
            vf(i,jm,k) = ramp*sbn(i,k)
!
!     South:
!
            uf(i,1,k) = ramp*tbs(i,k)
            vf(i,1,k) = ramp*sbs(i,k)
!
          end do
        end do
C
        do k=1,kbm1
!     Average corner cells
          uf(1, 1,k) = (uf(2, 1,k)+uf(1,   2,k))*.5
          vf(1, 1,k) = (vf(2, 1,k)+vf(1,   2,k))*.5
          uf(1,jm,k) = (uf(2,jm,k)+uf(1,jmm1,k))*.5
          vf(1,jm,k) = (vf(2,jm,k)+vf(1,jmm1,k))*.5
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)
            end do
          end do
        end do
C
        return
C
      else if(idx.eq.5) then
C
C     Vertical velocity boundary conditions:
C
!     Условие вертикальной скорости на границе:
!     Отсутствует. Просто накладывается маска свободной поверхности.
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              w(i,j,k)=w(i,j,k)*fsm(i,j)
            end do
          end do
        end do
C
        return
C
      else if(idx.eq.6) then
C
C     q2 and q2l boundary conditions:
C
!     ГУ турбулентной кинетической энергии и турб. масштаба:
        do k=1,kb
          do j=1,jm
C
C     East:
C
            uf(1,j,k) = 1.e-10
            vf(1,j,k) = 1.e-10
C
C     West:
C
            uf(im,j,k) = 1.e-10
            vf(im,j,k) = 1.e-10
!
          end do
        end do
C
        do k=1,kb
          do i=1,im
C
C     North:
C
            uf(i,1,k) = 1.e-10
            vf(i,1,k) = 1.e-10
C
C     South:
C
            uf(i,jm,k) = 1.e-10
            vf(i,jm,k) = 1.e-10
!
          end do
        end do
C
!     Наложение маски свободной поверхности с предотвращением получения нулевых значений.
!       Left from original pom2k
        do k=1,kb
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)+1.e-10
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)+1.e-10
            end do
          end do
        end do
C
        return
C
      endif
C
      end
C
      subroutine bcond(idx) ! bcond_ras(idx)
C **********************************************************************
C   Rochford and Shulman
C **********************************************************************

      implicit none
C
      include 'pom2k.c'
C
      integer idx
      real ga,u1,wm
      integer i,j,k
C
      if(idx.eq.1) then
C
C-----------------------------------------------------------------------
C
C     External (2-D) boundary conditions:
C
        do j=2,jmm1
!     West:
          ga = sqrt( (h(1,j)+h(2,j))*.5*grav )
     $        *dte/( (dx(1,j)+dx(2,j))*.5 )         ! Change 1->2 and 2->3 maybe...?
          elf(1,j) =     ga*(.25*el(2,j-1)+.5*el(2,j)+.25*el(2,j+1))
     $              +(1-ga)*(.25*el(1,j-1)+.5*el(1,j)+.25*el(1,j+1))
!     East:
          ga = sqrt( (h(im,j)+h(imm1,j))*.5*grav )
     $        *dte/( (dx(im,j)+dx(imm1,j))*.5 )
          elf(im,j) =     ga*(.25*el(imm1,j-1)+.5*el(imm1,j)
     $                        +.25*el(imm1,j+1))
     $               +(1-ga)*(.25*el( im ,j-1)+.5*el( im ,j)
     $                        +.25*el( im ,j+1))
        end do
!
        do i=2,imm1
!     South:
          ga = sqrt( (h(i,1)+h(i,2))*.5*grav )
     $        *dte/( (dy(i,1)+dy(i,2))*.5 )         ! Change 1->2 and 2->3 maybe...?
          elf(i,1) =     ga*(.25*el(i-1,2)+.5*el(i,2)+.25*el(i+1,2))
     $              +(1-ga)*(.25*el(i-1,1)+.5*el(i,1)+.25*el(i+1,1))
!     North:
          ga = sqrt( (h(i,jm)+h(i,jmm1))*.5*grav )
     $        *dte/( (dy(i,jm)+dy(i,jmm1))*.5 )
          elf(i,jm) =     ga*(.25*el(i-1,jmm1)+.5*el(i,jmm1)
     $                        +.25*el(i+1,jmm1))
     $               +(1-ga)*(.25*el(i-1,jm)+.5*el(i,jm)
     $                        +.25*el(i+1,jm))
        end do
!
        elf(1, 1) = (elf(1,  2 )+el(2, 1))*.5   ! Why elf+el? And not elf+elf?
        elf(1,jm) = (elf(1,jmm1)+el(2,jm))*.5
!       Two next lines were not present in RaSBCs.
        elf(im, 1) = (elf(im,  2 )+el(imm1, 1))*.5
        elf(im,jm) = (elf(im,jmm1)+el(imm1,jm))*.5
!
!        do j=1,jm
!          elf(1,j)=elf(2,j)     ! (з_внутр -> з_внешн)
!          elf(im,j)=elf(imm1,j) ! (в_внутр -> в_внешн)
!        end do
!
!        do i=1,im
!          elf(i,1)=elf(i,2)     ! (ю_внутр -> ю_внешн)
!          elf(i,jm)=elf(i,jmm1) ! (с_внутр -> с_внешн)
!        end do
!
C
!     Наложение маски свободной поверхности
        do j=1,jm
          do i=1,im
            elf(i,j)=elf(i,j)*fsm(i,j)
          end do
        end do
C
        return  ! выход из процедуры
C
      else if(idx.eq.2) then
C
C     External (2-D) velocity:
!     Баротропные скорости на границе:
!
        do j=2,jmm1
!
!     East:
!
          !u1 = ( h(im,j)+h(imm1,j) )*( dy(im,j)+dy(imm1,j) )*.25
          ga = sqrt( grav/h(imm1,j) )!*u1
          uaf(im,j)=( uabe(j)*ramp-ga*( ele(j)*ramp-el(imm1,j) ) )!/u1
!
          wm = ( ua(imm1,j-1)+ua(imm1,j) )*.5
          !u1 = ( h(im,j-1)+h(im,j) )*( dx(im,j-1)+dx(im,j) )*.25
          u1 = ramp*vabe(j)!/u1
          ga = dte/(dx(im,j-1)+dx(im,j)+dx(imm1,j-1)+dx(imm1,j))*2.
          vaf(im,j) = va(im,j)-ga*
     $               ( (wm+abs(wm))*(va(im,j)-va(imm1,j))
     $                +(wm-abs(wm))*(u1      -va( im ,j)) )
!          uaf(im,j) = 0.
!          vaf(im,j) = 0.
!
!          uaf(im,j) = uaf(imm1,j)
!
!     West:
!
          !u1 = ( h(1,j)+h(2,j) )*( dy(1,j)+dy(2,j) )*.25
          ga = sqrt( grav/h(2,j) )!*u1
          uaf(2,j) = ( uabw(j)*ramp-ga*( el(2,j)-elw(j)*ramp ) )!/u1
!          uaf(2,j) = 0.
!
          wm = ( ua(2,j-1)+ua(2,j) )*.5
          !u1 = ( h(1,j-1)+h(1,j) )*( dx(1,j-1)+dx(1,j) )*.25
          u1 = ramp*vabw(j)!/u1
          ga = dte/(dx(1,j-1)+dx(1,j)+dx(2,j-1)+dx(2,j))*2.
          vaf(1,j) = va(1,j)-ga*
     $               ( (wm+abs(wm))*(va(1,j)-u1     )
     $                +(wm-abs(wm))*(va(2,j)-va(1,j)) )
!
          uaf(1,j) = uaf(2,j)
!
        end do
!
        do i=2,imm1
!
!     North:
!
          !u1 = ( h(i,jm)+h(i,jmm1) )*( dx(i,jm)+dx(i,jmm1) )*.25
          ga = sqrt( grav/h(i,jmm1) )!*u1
          vaf(i,jm)=( vabn(i)*ramp-ga*( eln(i)*ramp-el(i,jmm1) ) )!/u1
!          vaf(i,jm) = 0.
!
          wm = ( va(i-1,jmm1)+va(i,jmm1) )*.5
          !u1 = ( h(i-1,jm)+h(i,jm) )*( dy(i-1,jm)+dy(i,jm) )*.25
          u1 = ramp*uabn(i)!/u1
          ga = dte/(dy(i-1,jm)+dy(i,jm)+dy(i-1,jmm1)+dy(i,jmm1))*2.
          uaf(i,jm) = ua(i,jm)-ga*
     $               ( (wm+abs(wm))*(ua(i,jm)-ua(i,jmm1))
     $                +(wm-abs(wm))*(u1      -ua(i, jm )) )
!
!          vaf(i,jm) = vaf(i,jmm1)
!
!     South:
!
          !u1 = ( h(i,1)+h(i,2) )*( dx(i,1)+dx(i,2) )*.25
          ga = sqrt( grav/h(i,2) )!*u1
          vaf(i,2) = ( vabs(i)*ramp-ga*( el(i,2)-els(i)*ramp ) )!/u1
!          vaf(i,2) = 0.
!
          wm = ( va(i-1,2)+va(i,2) )*.5
          !u1 = ( h(i-1,1)+h(i,1) )*( dy(i-1,1)+dy(i,1) )*.25
          u1 = ramp*uabs(i)!/u1
          ga = dte/(dy(i-1,1)+dy(i,1)+dy(i-1,2)+dy(i,2))*2.
          uaf(i,1) = ua(i,1)-ga*
     $               ( (wm+abs(wm))*(ua(i,1)-u1     )
     $                +(wm-abs(wm))*(ua(i,2)-ua(i,1)) )
!
          vaf(i,1) = vaf(i,2)
!
        end do
C
!     Наложение маски свободной поверхности
        do j=1,jm
          do i=1,im
            uaf(i,j)=uaf(i,j)*dum(i,j)
            vaf(i,j)=vaf(i,j)*dvm(i,j)
          end do
        end do
C
        return
C
      else if(idx.eq.3) then
C
C-----------------------------------------------------------------------
C
C     Internal (3-D) boundary conditions:
!     Бароклинные скорости на границе:
!
        do k=1,kbm1
          do j=2,jmm1
!
!     East:
!
            ga = sqrt(h(imm1,j)*grav*5.e-3)*dti/dx(imm1,j)
            uf(im,j,k) =     ga*(.25*u(imm1,j-1,k)+.5*u(imm1,j,k)
     $                          +.25*u(imm1,j+1,k))
     $                  +(1-ga)*(.25*u( im ,j-1,k)+.5*u( im ,j,k)
     $                          +.25*u( im ,j+1,k))
!
            u1 = 2.*dti/(dx(im,j-1)+dx(im,j)+dx(imm1,j-1)+dx(imm1,j))
            wm = ( u(im,j,k)+u(im,j-1,k) )*.5
            vf(im,j,k) = v(im,j,k) - u1*
     $                 ( (wm+abs(wm))*(u(im,j,k)-u(imm1,j,k))
     $                  +(wm-abs(wm))*(0.e0     -u( im ,j,k)) )
!            uf(im,j,k) = 0.
!            vf(im,j,k) = 0.
C
C     West:
C
            ga = sqrt(h(2,j)*grav*5.e-3)*dti/dx(2,j)
            uf(2,j,k) =     ga*(.25*u(3,j-1,k)+.5*u(3,j,k)
     $                         +.25*u(3,j+1,k))
     $                 +(1-ga)*(.25*u(2,j-1,k)+.5*u(2,j,k)
     $                         +.25*u(2,j+1,k))
!            uf(2,j,k) = 0.
!
            u1 = 2.*dti/(dx(1,j-1)+dx(1,j)+dx(2,j-1)+dx(2,j))
            wm = ( u(2,j,k)+u(2,j-1,k) )*.5
            vf(1,j,k) = v(1,j,k) - u1*
     $                 ( (wm+abs(wm))*(u(1,j,k)-0.e0    )
     $                  +(wm-abs(wm))*(u(2,j,k)-u(1,j,k)) )
!
            uf(1,j,k) = uf(2,j,k)
!
          end do
        end do
!     Fill corners for numerical safety
        uf(1, 1,:) = (uf(1,   2,:)+u(2, 1,:))*.5
        uf(1,jm,:) = (uf(1,jmm1,:)+u(2,jm,:))*.5
C
        do k=1,kbm1
          do i=2,imm1
C
C     North:
C
            ga = sqrt(h(i,jmm1)*grav*5.e-3)*dti/dy(i,jmm1)
            vf(i,jm,k) =     ga*(.25*v(i-1,jmm1,k)+.5*v(i,jmm1,k)
     $                           +.25*v(i+1,jmm1,k))
     $                  +(1-ga)*(.25*v(i-1, jm ,k)+.5*v(i, jm ,k)
     $                           +.25*v(i+1, jm ,k))
!            vf(i,jm,k) = 0.
!
            u1 = 2.*dti/(dy(i-1,jm)+dy(i,jm)+dy(i-1,jmm1)+dy(i,jmm1))
            wm = ( v(i,jm,k)+v(i-1,jm,k) )*.5
            uf(i,jm,k) = u(i,jm,k) - u1*
     $                 ( (wm+abs(wm))*(u(i,jm,k)-u(i,jmm1,k))
     $                  +(wm-abs(wm))*(0.e0     -u(i, jm ,k)) )
C
C     South:
C
            ga = sqrt(h(i,2)*grav*5.e-3)*dti/dy(i,2)
            vf(i,2,k) =     ga*(.25*v(i-1,3,k)+.5*v(i,3,k)
     $                          +.25*v(i+1,3,k))
     $                 +(1-ga)*(.25*v(i-1,2,k)+.5*v(i,2,k)
     $                          +.25*v(i+1,2,k))
!            vf(i,2,k) =     ga*(.25*v(i-1,2,k)+.5*v(i,2,k)
!     $                          +.25*v(i+1,2,k))
!     $                 +(1-ga)*(.25*vbs(i-1,k)+.5*vbs(i,k)
!     $                          +.25*vbs(i+1,k))
!            vf(i,2,k) = 0.
!
            u1 = 2.*dti/(dy(i-1,1)+dy(i,1)+dy(i-1,2)+dy(i,2))
            wm = ( v(i,2,k)+v(i-1,2,k) )*.5
            uf(i,1,k) = u(i,1,k) - u1*
     $                 ( (wm+abs(wm))*(u(i,1,k)-ubs(i,k))     ! 0.e0 -> ubs(i,k)
     $                  +(wm-abs(wm))*(u(i,2,k)-u(i,1,k)) )
!
            vf(i,1,k) = vf(i,2,k)
!
          end do
        end do
!     Fill corners for numerical safety
        vf( 1,1,:) = (vf(   2,1,:)+v(2, 1,:))*.5
        vf(im,1,:) = (vf(imm1,1,:)+v(im,2,:))*.5
C
!     Наложение маски свободной поверхности
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*dum(i,j)
              vf(i,j,k)=vf(i,j,k)*dvm(i,j)
            end do
          end do
        end do
C
        return
C
      else if(idx.eq.4) then
C
C     Temperature and salinity boundary conditions (using uf and vf,
C     respectively):
C
!     ГУ температуры и солёности:
        do k=1,kbm1
          do j=2,jmm1
C
C     East:
C
            ga = dti/(dx(im,j)+dx(imm1,j))
            uf(im,j,k) = t(im,j,k) - ga *
     $                 ( ( u(imm1,j,k)-abs(u(imm1,j,k)) )
     $                  *(t(im,j,k)-tbe(j,k))
     $                  +( u(imm1,j,k)+abs(u(imm1,j,k)) )
     $                  *(t(imm1,j,k)-t(im,j,k)) )
            vf(im,j,k) = s(im,j,k) - ga *
     $                 ( ( u(imm1,j,k)-abs(v(imm1,j,k)) )
     $                  *(s(im,j,k)-sbe(j,k))
     $                  +( u(imm1,j,k)+abs(v(imm1,j,k)) )
     $                  *(s(imm1,j,k)-s(im,j,k)) )
C
C     West:
C
            ga = dti/(dx(1,j)+dx(2,j))
            uf(1,j,k) = t(1,j,k) - ga *
     $                 ( ( u(2,j,k)+abs(u(2,j,k)) )
     $                  *(t(1,j,k)-tbw(j,k))
     $                  +( u(2,j,k)-abs(u(2,j,k)) )
     $                  *(t(2,j,k)-t(1,j,k)) )
            vf(1,j,k) = s(1,j,k) - ga *
     $                 ( ( u(2,j,k)+abs(v(2,j,k)) )
     $                  *(s(1,j,k)-sbw(j,k))
     $                  +( u(2,j,k)-abs(v(2,j,k)) )
     $                  *(s(2,j,k)-s(1,j,k)) )
          end do
        end do
C
        do k=1,kbm1
          do i=2,imm1
C
C     North:
C
            ga = dti/(dy(i,jm)+dy(i,jmm1))
            uf(i,jm,k) = t(i,jm,k) - ga *
     $                 ( ( v(i,jmm1,k)-abs(v(i,jmm1,k)) )
     $                  *(tbn(i,k)-t(i,jm,k))   !!!!! Revise!!!!
     $                  +( v(i,jmm1,k)+abs(v(i,jmm1,k)) )
     $                  *(t(i,jm,k)-t(i,jmm1,k)) ) !!
            vf(i,jm,k) = s(i,jm,k) - ga *
     $                 ( ( v(i,jmm1,k)-abs(v(i,jmm1,k)) )
     $                  *(sbn(i,k)-s(i,jm,k))   !!!!! R!
     $                  +( v(i,jmm1,k)+abs(v(i,jmm1,k)) )
     $                  *(s(i,jm,k)-s(i,jmm1,k)) ) !!
C
C     South:
C
            ga = dti/(dy(i,1)+dy(i,2))
            uf(i,1,k) = t(i,1,k) - ga *
     $                 ( ( v(i,2,k)+abs(v(i,2,k)) )
     $                  *(t(i,1,k)-tbs(i,k))
     $                  +( v(i,2,k)-abs(v(i,2,k)) )
     $                  *(t(i,2,k)-t(i,1,k)) )
            vf(i,1,k) = s(i,1,k) - ga *
     $                 ( ( v(i,2,k)+abs(v(i,2,k)) )
     $                  *(s(i,1,k)-sbs(i,k))
     $                  +( v(i,2,k)-abs(v(i,2,k)) )
     $                  *(s(i,2,k)-s(i,1,k)) )
          end do
        end do
C
        do k=1,kbm1
!     Average corner cells
          uf(1, 1,k) = (uf(2, 1,k)+uf(1,   2,k))*.5
          vf(1, 1,k) = (vf(2, 1,k)+vf(1,   2,k))*.5
          uf(1,jm,k) = (uf(2,jm,k)+uf(1,jmm1,k))*.5
          vf(1,jm,k) = (vf(2,jm,k)+vf(1,jmm1,k))*.5
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)
            end do
          end do
        end do
C
        return
C
      else if(idx.eq.5) then
C
C     Vertical velocity boundary conditions:
C
!     Условие вертикальной скорости на границе:
        do k=1,kbm1
          do j=1,jm
            w( 1,j,k) = w(   2,j,k)*fsm( 1,j)
            w(im,j,k) = w(imm1,j,k)*fsm(im,j)
          end do
          do i=2,im
            w(i, 1,k) = w(i,   2,k)*fsm(i,   1)
            w(i,jm,k) = w(i,jmm1,k)*fsm(i,jmm1)
          end do
        end do
C
        return
C
      else if(idx.eq.6) then
C
C     q2 and q2l boundary conditions:
C
!     ГУ турбулентной кинетической энергии и турб. масштаба:
        do k=1,kb
          do j=1,jm
C
C     East:
C
            uf(1,j,k) = 1.e-10
            vf(1,j,k) = 1.e-10
C
C     West:
C
            uf(im,j,k) = 1.e-10
            vf(im,j,k) = 1.e-10
!
          end do
        end do
C
        do k=1,kb
          do i=1,im
C
C     North:
C
            uf(i,1,k) = 1.e-10
            vf(i,1,k) = 1.e-10
C
C     South:
C
            uf(i,jm,k) = 1.e-10
            vf(i,jm,k) = 1.e-10
!
          end do
        end do
C
!     Наложение маски свободной поверхности с предотвращением получения нулевых значений.
!       Left from original pom2k
        do k=1,kb
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)+1.e-10
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)+1.e-10
            end do
          end do
        end do
C
        return
C
      endif
C
      end
C
      subroutine bcondorl(idx)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  This is an optional subroutine replacing  bcond and *
C *                using Orlanski's scheme (J. Comp. Phys. 21, 251-269,*
C *                1976), specialized for the seamount problem. To     *
C *                make it work for the seamount problem, I (G.M.)     *
C *                have had to add an extra condition on an "if"       *
C *                statement in the t and s open boundary conditions,  *
C *                which involves the sign of the normal velocity.     *
C *                Thus:                                               *
C *                                                                    *
C *            if(cl.eq.0.e0.and.ubw(j,k).ge.0.e0) uf(1,j,k)=tbw(j,k), *
C *                                                                    *
C *                plus 3 others of the same kind.                     *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      integer idx
      real cl,denom
      integer i,j,k
C
      if(idx.eq.1) then
C
C-----------------------------------------------------------------------
C
C     External (2-D) boundary conditions:
C
C     In this example the governing boundary conditions are a radiation
C     condition on uaf(im,j) in the east and an inflow uaf(2,j) in the
C     west. The tangential velocities are set to zero on both
C     boundaries. These are only one set of possibilities and may not
C     represent a choice which yields the most physically realistic
C     result.
C
C     Elevation (in this application, elevation is not a primary
C     boundary condition):
C
        do  j=1,jm
          elf(1,j)=elf(2,j)
          elf(im,j)=elf(imm1,j)
        end do
C
        do j=1,jm
          do i=1,im
            elf(i,j)=elf(i,j)*fsm(i,j)
          end do
        end do
C
        return
C
      else if(idx.eq.2) then
C
C     External (2-D) velocity:
C
        do j=2,jmm1
C
C     West:
C
          uaf(2,j)=ramp*uabw(j)-sqrt(grav/h(2,j))*(el(2,j)-elw(j))
          uaf(1,j)=uaf(2,j)
          vaf(1,j)=0.e0
C
C     East:
C
          uaf(im,j)=ramp*uabe(j)
     $               +sqrt(grav/h(imm1,j))*(el(imm1,j)-ele(j))
          vaf(im,j)=0.e0
C
        end do
C
        do j=1,jm
          do i=1,im
            uaf(i,j)=uaf(i,j)*dum(i,j)
            vaf(i,j)=vaf(i,j)*dvm(i,j)
          end do
        end do
C
        return
C
      else if(idx.eq.3) then
C
C-----------------------------------------------------------------------
C
C     Internal (3-D) boundary conditions:
C
C     Eastern and western radiation boundary conditions according to
C     Orlanski's explicit scheme:
C
        do k=1,kbm1
          do j=2,jmm1
C
C     West:
C
            denom=(uf(3,j,k)+ub(3,j,k)-2.e0*u(4,j,k))
            if(denom.eq.0.e0)denom=0.01e0
            cl=(ub(3,j,k)-uf(3,j,k))/denom
            if(cl.gt.1.e0) cl=1.e0
            if(cl.lt.0.e0) cl=0.e0
            uf(2,j,k)=(ub(2,j,k)*(1.e0-cl)+2.e0*cl*u(3,j,k))
     $                 /(1.e0+cl)
            uf(1,j,k)=uf(2,j,k)
            vf(1,j,k)=0.e0
C
C     East:
C
            denom=(uf(im-1,j,k)+ub(im-1,j,k)-2.e0*u(im-2,j,k))
            if(denom.eq.0.e0)denom=0.01e0
            cl=(ub(im-1,j,k)-uf(im-1,j,k))/denom
            if(cl.gt.1.e0) cl=1.e0
            if(cl.lt.0.e0) cl=0.e0
            uf(im,j,k)=(ub(im,j,k)*(1.e0-cl)+2.e0*cl*u(im-1,j,k))
     $                  /(1.e0+cl)
            vf(im,j,k)=0.e0
          end do
        end do
C
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*dum(i,j)
              vf(i,j,k)=vf(i,j,k)*dvm(i,j)
            end do
          end do
        end do
C
        return
C
      else if(idx.eq.4) then
C
C     Temperature and salinity boundary conditions (using uf and vf,
C     respectively):
C
        do k=1,kbm1
          do j=1,jm
C
C     West:
C
            ubw(j,k)=ub(2,j,k)
            denom=(uf(2,j,k)+tb(2,j,k)-2.e0*t(3,j,k))
            if(denom.eq.0.e0) denom=0.01e0
            cl=(tb(2,j,k)-uf(2,j,k))/denom
            if(cl.gt.1.e0) cl=1.e0
            if(cl.lt.0.e0) cl=0.e0
            uf(1,j,k)=(tb(1,j,k)*(1.e0-cl)+2.e0*cl*t(2,j,k))/(1.e0+cl)
            if(cl.eq.0.e0.and.ubw(j,k).ge.0.e0) uf(1,j,k)=tbw(j,k)
C
            denom=(vf(2,j,k)+sb(2,j,k)-2.e0*s(3,j,k))
            if(denom.eq.0.e0) denom=0.01e0
            cl=(sb(2,j,k)-vf(2,j,k))/denom
            if(cl.gt.1.e0) cl=1.e0
            if(cl.lt.0.e0) cl=0.e0
            vf(1,j,k)=(sb(1,j,k)*(1.e0-cl)+2.e0*cl*s(2,j,k))/(1.e0+cl)
            if(cl.eq.0.e0.and.ubw(j,k).ge.0.e0) vf(1,j,k)=sbw(j,k)
C
C     East:
C
            ube(j,k)=ub(im,j,k)
            denom=(uf(im-1,j,k)+tb(im-1,j,k)-2.e0*t(im-2,j,k))
            if(denom.eq.0.e0) denom=0.01e0
            cl=(tb(im-1,j,k)-uf(im-1,j,k))/denom
            if(cl.gt.1.e0) cl=1.e0
            if(cl.lt.0.e0) cl=0.e0
            uf(im,j,k)=(tb(im,j,k)*(1.e0-cl)+2.e0*cl*t(im-1,j,k))
     $                  /(1.e0+cl)
            if(cl.eq.0.e0.and.ube(j,k).le.0.e0) uf(im,j,k)=tbe(j,k)
C
            denom=(vf(im-1,j,k)+sb(im-1,j,k)-2.e0*s(im-2,j,k))
            if(denom.eq.0.e0) denom=0.01e0
            cl=(sb(im-1,j,k)-vf(im-1,j,k))/denom
            if(cl.gt.1.e0) cl=1.e0
            if(cl.lt.0.e0) cl=0.e0
            vf(im,j,k)=(sb(im,j,k)*(1.e0-cl)+2.e0*cl*s(im-1,j,k))
     $                  /(1.e0+cl)
            if(cl.eq.0.e0.and.ube(j,k).le.0.e0) vf(im,j,k)=sbe(j,k)
C
          end do
        end do
C
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)
            end do
          end do
        end do
C
        return
C
      else if(idx.eq.5) then
C
C     Vertical velocity boundary conditions:
C
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              w(i,j,k)=w(i,j,k)*fsm(i,j)
            end do
          end do
        end do
C
        return
C
      else if(idx.eq.6) then
C
C     q2 and q2l boundary conditions:
C
        do k=1,kb
C
          do j=1,jm
            uf(im,j,k)=1.e-10
            vf(im,j,k)=1.e-10
            uf(1,j,k)=1.e-10
            vf(1,j,k)=1.e-10
          end do
C
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)
            end do
          end do
        end do
C
        return
C
      endif
C
      end
C
      subroutine box
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Sets up conservation box problem.                   *
C *                                                                    *
C *                This basin uses the same grid as the seamount       *
C *                problem, but it has a flat bottom, is surrounded by *
C *                walls and is initialised with uniform salinity and  *
C *                temperature. It is forced by a surface input of     *
C *                water of the same temperature and salinity as the   *
C *                water in the basin. Therefore, the temperature and  *
C *                salinity in the basin should not change, and the    *
C *                free surface should fall at a rate vflux. It is also*
C *                forced by a steady atmospheric pressure field which *
C *                depresses the southwestern half of the model by 1 m *
C *                and elevates the northeastern half of the model by  *
C *                1 m.                                                *
C *                                                                    *
C *                Since this problem defines its own fixed e_atmos,   *
C *                tatm, satm and e_atmos, comment out corresponding   *
C *                declarations after the do 9000 statement in main    *
C *                program.                                            *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real depth,delx,tatm,satm
      integer i,j,k
C
C     Water depth:
C
      depth=4500.e0
C
C     Grid size:
C
      delx=8000.e0
C
C     Set up grid dimensions, areas of free surface cells, and
C     Coriolis parameter:
C
      do j=1,jm
        do i=1,im
C
C     For constant grid size:
C
C         dx(i,j)=delx
C         dy(i,j)=delx
C
C     For variable grid size:
C
          dx(i,j)=delx-delx*sin(pi*float(i)/float(im))/2.e0
          dy(i,j)=delx-delx*sin(pi*float(j)/float(jm))/2.e0
C
          cor(i,j)=1.e-4
C
        end do
      end do
C
C     Calculate horizontal coordinates of grid points and rotation
C     angle.
C
C     NOTE that this is introduced solely for the benefit of any post-
C     processing software, and in order to conform with the requirements
C     of the NetCDF Climate and Forecast (CF) Metadata Conventions.
C
C     There are four horizontal coordinate systems, denoted by the
C     subscripts u, v, e and c ("u" is a u-point, "v" is a v-point,
C     "e" is an elevation point and "c" is a cell corner), as shown
C     below. In addition, "east_*" is an easting and "north_*" is a
C     northing. Hence the coordinates of the "u" points are given by
C     (east_u,north_u).
C
C     Also, if the centre point of the cell shown below is at
C     (east_e(i,j),north_e(i,j)), then (east_u(i,j),north_u(i,j)) are
C     the coordinates of the western of the two "u" points,
C     (east_v(i,j),north_v(i,j)) are the coordinates of the southern of
C     the two "v" points, and (east_c(i,j),north_c(i,j)) are the
C     coordinates of the southwestern corner point of the cell. The
C     southwest corner of the entire grid is at
C     (east_c(1,1),north_c(1,1)).
C
C                      |              |
C                    --c------v-------c--
C                      |              |
C                      |              |
C                      |              |
C                      |              |
C                      u      e       u
C                      |              |
C                      |              |
C                      |              |
C                      |              |
C                    --c------v-------c--
C                      |              |
C
C
C     NOTE that the following calculation of east_c and north_c only
C     works properly for a rectangular grid with east and north aligned
C     with i and j, respectively:
C
      do j=1,jm
        east_c(1,j)=0.e0
        do i=2,im
          east_c(i,j)=east_c(i-1,j)+dx(i-1,j)
        end do
      end do
C
      do i=1,im
        north_c(i,1)=0.e0
        do j=2,jm
          north_c(i,j)=north_c(i,j-1)+dy(i,j-1)
        end do
      end do
C
C     The following works properly for any grid:
C
C     Elevation points:
C
      do j=1,jm-1
        do i=1,im-1
          east_e(i,j)=(east_c(i,j)+east_c(i+1,j)
     $                  +east_c(i,j+1)+east_c(i+1,j+1))/4.e0
          north_e(i,j)=(north_c(i,j)+north_c(i+1,j)
     $                   +north_c(i,j+1)+north_c(i+1,j+1))/4.e0
        end do
      end do
C
C     Extrapolate ends:
C
      do i=1,im-1
        east_e(i,jm)
     $    =((east_c(i,jm)+east_c(i+1,jm))*3.e0
     $       -east_c(i,jm-1)-east_c(i+1,jm-1))/4.e0
        north_e(i,jm)
     $    =((north_c(i,jm)+north_c(i+1,jm))*3.e0
     $       -north_c(i,jm-1)-north_c(i+1,jm-1))/4.e0
      end do
C
      do j=1,jm-1
        east_e(im,j)
     $    =((east_c(im,j)+east_c(im,j+1))*3.e0
     $       -east_c(im-1,j)-east_c(im-1,j+1))/4.e0
        north_e(im,j)
     $    =((north_c(im,j)+north_c(im,j+1))*3.e0
     $       -north_c(im-1,j)-north_c(im-1,j+1))/4.e0
      end do
C
      east_e(im,jm)=east_e(im-1,jm)+east_e(im,jm-1)
     $               -(east_e(im-2,jm)+east_e(im,jm-2))/2.e0
      north_e(im,jm)=north_e(im-1,jm)+north_e(im,jm-1)
     $               -(north_e(im-2,jm)+north_e(im,jm-2))/2.e0
C
C     u-points:
C
      do j=1,jm-1
        do i=1,im
          east_u(i,j)=(east_c(i,j)+east_c(i,j+1))/2.e0
          north_u(i,j)=(north_c(i,j)+north_c(i,j+1))/2.e0
        end do
      end do
C
C     Extrapolate ends:
C
      do i=1,im
        east_u(i,jm)=(east_c(i,jm)*3.e0-east_c(i,jm-1))/2.e0
        north_u(i,jm)=(north_c(i,jm)*3.e0-north_c(i,jm-1))/2.e0
      end do
C
C     v-points:
C
      do j=1,jm
        do i=1,im-1
          east_v(i,j)=(east_c(i,j)+east_c(i+1,j))/2.e0
          north_v(i,j)=(north_c(i,j)+north_c(i+1,j))/2.e0
        end do
      end do
C
C     Extrapolate ends:
C
      do j=1,jm
        east_v(im,j)=(east_c(im,j)*3.e0-east_c(im-1,j))/2.e0
        north_v(im,j)=(north_c(im,j)*3.e0-north_c(im-1,j))/2.e0
      end do
C
C     rot is the angle (radians, anticlockwise) of the i-axis relative
C     to east, averaged to a cell centre:
C
C     (NOTE that the following calculation of rot only works properly
C     for this particular rectangular grid)
C
      do j=1,jm
        do i=1,im
          rot(i,j)=0.e0
        end do
      end do
C
C     Define depth:
C
      do i=1,im
        do j=1,jm
          h(i,j)=depth
        end do
      end do
C
C     Close the north and south boundaries:
C
      do i=1,im
        h(i,1)=1.e0
        h(i,jm)=1.e0
      end do
C
C     Close the east and west boundaries:
C
      do j=1,jm
        h(1,j)=1.e0
        h(im,j)=1.e0
      end do
C
C     Calculate areas and masks:
C
      call areas_masks
C
C     Adjust bottom topography so that cell to cell variations
C     in h do not exceed parameter slmax:
C
      if(slmax.lt.1.e0) call slpmax
C
C     Set tbias and sbias here for test (tbias and sbias would
C     normally only be set in the main program):
C
      tbias=10.e0
      sbias=20.e0
      write(6,1) tbias,sbias
    1 format(/' tbias and sbias changed in subroutine box to:'/
     $         2f10.3//)
C
C     Set initial conditions:
C
      do k=1,kbm1
        do j=1,jm
          do i=1,im
            tb(i,j,k)=20.e0-tbias
            sb(i,j,k)=35.e0-sbias
            tclim(i,j,k)=tb(i,j,k)
            sclim(i,j,k)=sb(i,j,k)
          end do
        end do
      end do
C
C     Initialise uab and vab as necessary
C     (NOTE that these have already been initialised to zero in the
C     main program):
C
      do j=1,jm
        do i=1,im
C     No conditions necessary for this problem
        end do
      end do
C
C     Set surface boundary conditions, e_atmos, vflux, wusurf,
C     wvsurf, wtsurf, wssurf and swrad, as necessary
C     (NOTE:
C      1. These have all been initialised to zero in the main program.
C      2. The temperature and salinity of inflowing water must be
C         defined relative to tbias and sbias.):
C
      do j=1,jm
        do i=1,im
          if(i+j-57.le.0) then
            e_atmos(i,j)=1.e0
          else
            e_atmos(i,j)=-1.e0
          endif
C
C     Ensure atmospheric pressure cannot make water depth go negative:
C
          e_atmos(i,j)=min(e_atmos(i,j),h(i,j))
C
          vfluxf(i,j)=-0.0001e0
C
C     See main program, just after "Begin numerical integration", for
C     an explanation of these terms:
C
          tatm=20.e0
          satm=35.e0
C
        end do
      end do
C
C     Initialise elb, etb, dt and aam2d:
C
      do j=1,jm
        do i=1,im
          elb(i,j)=-e_atmos(i,j)
          etb(i,j)=-e_atmos(i,j)
          dt(i,j)=h(i,j)-e_atmos(i,j)
          aam2d(i,j)=aam(i,j,1)
        end do
      end do
C
      call dens(sb,tb,rho)
C
C     Generated horizontally averaged density field (in this
C     application, the initial condition for density is a function
C     of z (the vertical cartesian coordinate) -- when this is not
C     so, make sure that rmean has been area averaged BEFORE transfer
C     to sigma coordinates):
C
      do k=1,kbm1
        do j=1,jm
          do i=1,im
            rmean(i,j,k)=rho(i,j,k)
          end do
        end do
      end do
C
C     Set lateral boundary conditions, for use in subroutine bcond
C     (in this problem, all lateral boundaries are closed through
C     the specification of the masks fsm, dum and dvm):
C
      rfe=1.e0
      rfw=1.e0
      rfn=1.e0
      rfs=1.e0
C
C     Set thermodynamic boundary conditions (for the seamount
C     problem, and other possible applications, lateral thermodynamic
C     boundary conditions are set equal to the initial conditions and
C     are held constant thereafter - users may, of course, create
C     variable boundary conditions):
C
      do k=1,kbm1
C
        do j=1,jm
          tbe(j,k)=tb(im,j,k)
          tbw(j,k)=tb(1,j,k)
          sbe(j,k)=sb(im,j,k)
          sbw(j,k)=sb(1,j,k)
        end do
C
        do i=1,im
          tbn(i,k)=tb(i,jm,k)
          tbs(i,k)=tb(i,1,k)
          sbn(i,k)=sb(i,jm,k)
          sbs(i,k)=sb(i,1,k)
        end do
C
      end do
C
      return
C
      end
C
      subroutine cortes                 ! rwnd :
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Sets up coriolis test.                              *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real depth,delx,tatm,satm
      real rad,re,dlat,dlon
      integer i,j,k
      rad=0.01745329
      re=6371.E3
C
C     Sigma layers:
C
      do k=1,kb
        z(k) = -float(k-1)/float(kb-1)
      end do
      do k=1,kbm1
        dz(k) = z(k)-z(k+1)
        zz(k) = (z(k)+z(k+1))/2
      end do
      dz(kb) = 0.
      zz(kb) = zz(kbm1)-dz(kbm1)
      do k=1,kbm1
        dzz(k) = zz(k)-zz(k+1)
      end do
      dzz(kb) = 0.
C
C     Water depth:
C
      depth=4500.e0
C
C     Grid size:
C
      delx=8000.e0
C
C     Set up grid dimensions, areas of free surface cells, and
C     Coriolis parameter:
C
      do j=1,jm
        do i=1,im

          north_e(i,j) = float(j)*180/float(jm+1)-90.
          east_e(i,j) = float(i)*360/float(im+1)-180.
          cor(i,j)=2.*7.29E-5*sin(north_e(i,j)*rad)
C
        end do
      end do
C
        do j=1,jm
          do i=2,im-1
            dx(i,j)=0.5*rad*re*sqrt(((east_e(i+1,j)-east_e(i-1,j))
     1 *cos(north_e(i,j)*rad))**2+(north_e(i+1,j)-north_e(i-1,j))**2)
          end do
            dx(1,j)=dx(2,j)
            dx(im,j)=dx(im-1,j)
        end do
C
        do i=1,im
          do j=2,jm-1
            dy(i,j)=0.5*rad*re*sqrt(((east_e(i,j+1)-east_e(i,j-1))
     1 *cos(north_e(i,j)*rad))**2+(north_e(i,j+1)-north_e(i,j-1))**2)
          end do
            dy(i,1)=dy(i,2)
            dy(i,jm)=dy(i,jm-1)
        end do
C
C --- the following grids are needed only for netcdf plotting
C
C     Corner of cell points:
C
      do j=2,jm
        do i=2,im
          east_c(i,j)=(east_e(i,j)+east_e(i-1,j)
     $                  +east_e(i,j-1)+east_e(i-1,j-1))/4.e0
          north_c(i,j)=(north_e(i,j)+north_e(i-1,j)
     $                   +north_e(i,j-1)+north_e(i-1,j-1))/4.e0
        end do
      end do
C
C
C     Extrapolate ends (approx.):
C
      do i=2,im
        east_c(i,1)=2.*east_c(i,2)-east_c(i,3)
        north_c(i,1)=2.*north_c(i,2)-north_c(i,3)
      end do
        east_c(1,1)=2.*east_c(2,1)-east_c(3,1)
C
      do j=2,jm
        east_c(1,j)=2.*east_c(2,j)-east_c(3,j)
        north_c(1,j)=2.*north_c(2,j)-north_c(3,j)
      end do
        north_c(1,1)=2.*north_c(1,2)-north_c(1,3)
C
C     u-points:
C
      do j=1,jm-1
        do i=1,im
          east_u(i,j)=(east_c(i,j)+east_c(i,j+1))/2.e0
          north_u(i,j)=(north_c(i,j)+north_c(i,j+1))/2.e0
        end do
      end do
C
C     Extrapolate ends:
C
      do i=1,im
        east_u(i,jm)=(east_c(i,jm)*3.e0-east_c(i,jm-1))/2.e0
        north_u(i,jm)=(north_c(i,jm)*3.e0-north_c(i,jm-1))/2.e0
      end do
C
C     v-points:
C
      do j=1,jm
        do i=1,im-1
          east_v(i,j)=(east_c(i,j)+east_c(i+1,j))/2.e0
          north_v(i,j)=(north_c(i,j)+north_c(i+1,j))/2.e0
        end do
      end do
C
C     Extrapolate ends:
C
      do j=1,jm
        east_v(im,j)=(east_c(im,j)*3.e0-east_c(im-1,j))/2.e0
        north_v(im,j)=(north_c(im,j)*3.e0-north_c(im-1,j))/2.e0
      end do
C
C     rot is the angle (radians, anticlockwise) of the i-axis relative
C     to east, averaged to a cell centre: (only needed for CDF plotting)
C
      do j=1,jm
        do i=1,im-1
          rot(i,j)=0.
          dlat=north_e(i+1,j)-north_e(i,j)
          dlon= east_e(i+1,j)- east_e(i,j)
           if(dlon.ne.0.) rot(i,j)=atan(dlat/dlon)
        end do
       rot(im,j)=rot(im-1,j)
      end do
C
C     Define depth:
C
      do i=1,im
        do j=1,jm
          h(i,j)=depth
        end do
      end do
C
C     Close the north and south boundaries:
C
      do i=1,im
        h(i,1)=1.e0
        h(i,jm)=1.e0
      end do
C
C     Close the east and west boundaries:
C
      do j=1,jm
        h(1,j)=1.e0
        h(im,j)=1.e0
      end do
C
C     Calculate areas and masks:
C
      call areas_masks
C
C     Adjust bottom topography so that cell to cell variations
C     in h do not exceed parameter slmax:
C
      if(slmax.lt.1.e0) call slpmax
C
C     Set tbias and sbias here for test (tbias and sbias would
C     normally only be set in the main program):
C
      tbias=10.e0
      sbias=20.e0
      write(6,1) tbias,sbias
    1 format(/' tbias and sbias changed in subroutine cortes to:'/
     $         2f10.3//)
C
C     Set initial conditions:
C
      do k=1,kbm1
        do j=1,jm
          do i=1,im
            tb(i,j,k)=20.e0-tbias
            sb(i,j,k)=35.e0-sbias
            tclim(i,j,k)=tb(i,j,k)
            sclim(i,j,k)=sb(i,j,k)
          end do
        end do
      end do
C
C     Initialise uab and vab as necessary
C     (NOTE that these have already been initialised to zero in the
C     main program):
C
      do j=1,jm
        do i=1,im
C     No conditions necessary for this problem
        end do
      end do
C
C     Set surface boundary conditions, e_atmos, vflux, wusurf,
C     wvsurf, wtsurf, wssurf and swrad, as necessary
C     (NOTE:
C      1. These have all been initialised to zero in the main program.
C      2. The temperature and salinity of inflowing water must be
C         defined relative to tbias and sbias.):
C
C
C     See main program, just after "Begin numerical integration", for
C     an explanation of these terms:
C
          tatm=20.e0
          satm=35.e0
C
C     Initialise elb, etb, dt and aam2d:
C
      do j=1,jm
        do i=1,im
          elb(i,j)=-.5*(j-1)/jm+.5
          etb(i,j)=elb(i,j)
          dt(i,j)=h(i,j)+etb(i,j)
          aam2d(i,j)=aam(i,j,1)
        end do
      end do
C
      call dens(sb,tb,rho)
C
C     Generated horizontally averaged density field (in this
C     application, the initial condition for density is a function
C     of z (the vertical cartesian coordinate) -- when this is not
C     so, make sure that rmean has been area averaged BEFORE transfer
C     to sigma coordinates):
C
      do k=1,kbm1
        do j=1,jm
          do i=1,im
            rmean(i,j,k)=rho(i,j,k)
          end do
        end do
      end do
C
C     Set lateral boundary conditions, for use in subroutine bcond
C     (in this problem, all lateral boundaries are closed through
C     the specification of the masks fsm, dum and dvm):
C
      rfe=1.e0
      rfw=1.e0
      rfn=1.e0
      rfs=1.e0
C
C     Set thermodynamic boundary conditions (for the seamount
C     problem, and other possible applications, lateral thermodynamic
C     boundary conditions are set equal to the initial conditions and
C     are held constant thereafter - users may, of course, create
C     variable boundary conditions):
C
      do k=1,kbm1
C
        do j=1,jm
          tbe(j,k)=tb(im,j,k)
          tbw(j,k)=tb(1,j,k)
          sbe(j,k)=sb(im,j,k)
          sbw(j,k)=sb(1,j,k)
        end do
C
        do i=1,im
          tbn(i,k)=tb(i,jm,k)
          tbs(i,k)=tb(i,1,k)
          sbn(i,k)=sb(i,jm,k)
          sbs(i,k)=sb(i,1,k)
        end do
C
      end do
C
      return
C
      end
C
      subroutine dens(si,ti,rhoo)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Calculates (density-1000.)/rhoref.                  *
C *                                                                    *
C *                (see: Mellor, G.L., 1991, J. Atmos. Oceanic Tech.,  *
C *                609-611.)                                           *
C *                                                                    *
C *                ti is potential temperature                         *
C *                                                                    *
C *                If using 32 bit precision, it is recommended that   *
C *                cr,p,rhor,sr,tr,tr2,tr3 and tr4 be made double      *
C *                precision, and the "e"s in the constants be changed *
C *                to "d"s.                                            *
C *                                                                    *
C * NOTE: if pressure is not used in dens, buoyancy term (boygr)       *
C *       in profq must be changed (see note in profq)                 *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real si(im,jm,kb),ti(im,jm,kb),rhoo(im,jm,kb)
      real cr,p,rhor,sr,tr,tr2,tr3,tr4
      integer i,j,k
C
      do k=1,kbm1
        do j=1,jm
          do i=1,im
C
            tr=ti(i,j,k)+tbias
            sr=si(i,j,k)+sbias
            tr2=tr*tr
            tr3=tr2*tr
            tr4=tr3*tr
C
C     Approximate pressure in units of bars:
C
            p=grav*rhoref*(-zz(k)* h(i,j))*1.e-5
C
            rhor=-0.157406e0+6.793952e-2*tr
     $            -9.095290e-3*tr2+1.001685e-4*tr3
     $            -1.120083e-6*tr4+6.536332e-9*tr4*tr
C
            rhor=rhor+(0.824493e0-4.0899e-3*tr
     $            +7.6438e-5*tr2-8.2467e-7*tr3
     $            +5.3875e-9*tr4)*sr
     $            +(-5.72466e-3+1.0227e-4*tr
     $            -1.6546e-6*tr2)*abs(sr)**1.5
     $            +4.8314e-4*sr*sr
C
            cr=1449.1e0+.0821e0*p+4.55e0*tr-.045e0*tr2
     $          +1.34e0*(sr-35.e0)
            rhor=rhor+1.e5*p/(cr*cr)*(1.e0-2.e0*p/(cr*cr))
C
            rhoo(i,j,k)=rhor/rhoref*fsm(i,j)
C
          end do
        end do
      end do
C
      return
C
      end
C
      subroutine depth
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Establishes the vertical sigma grid with log        *
C *                distributions at the top and bottom and a linear    *
C *                distribution in between. The number of layers of    *
C *                reduced thickness are kl1-2 at the surface and      *
C *                kb-kl2-1 at the bottom. kl1 and kl2 are defined in  *
C *                the main program. For no log portions, set kl1=2    *
C *                and kl2=kb-1.                                       *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real delz
      integer kdz(12)
      integer k
C
      data kdz/1,1,2,4,8,16,32,64,128,256,512,1024/
C
      z(1)=0.e0
C
      do k=2,kl1
        z(k)=z(k-1)+kdz(k-1)
      end do
C
      delz=z(kl1)-z(kl1-1)
C
      do k=kl1+1,kl2
        z(k)=z(k-1)+delz
      end do
C
      do k=kl2+1,kb
        dz(k)=float(kdz(kb-k+1))*delz/float(kdz(kb-kl2))
        z(k)=z(k-1)+dz(k)
      end do
C
      do k=1,kb
        z(k)=-z(k)/z(kb)
      end do
C
      do k=1,kb-1
        zz(k)=0.5e0*(z(k)+z(k+1))
      end do
C
      zz(kb)=2.e0*zz(kb-1)-zz(kb-2)
C
      do k=1,kb-1
        dz(k)=z(k)-z(k+1)
        dzz(k)=zz(k)-zz(k+1)
      end do
C
      dz(kb)=0.e0
      dzz(kb)=0.e0
C
      write(6,1)
    1 format(/2x,'k',7x,'z',9x,'zz',9x,'dz',9x,'dzz',/)
C
      do k=1,kb
        write(6,2) k,z(k),zz(k),dz(k),dzz(k)
    2   format((' ',i5,4f10.3))
      end do
C
      write(6,3)
    3 format(//)
C
      return
C
      end
C
      subroutine file2ic
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Sets up my own problem.                             *
C *                                                                    *
C * This example read IC from IC.dat file, generated by GRID.f in      *
C * GRID-DATA directory. Only minimal number of fields are read,       *
C * while others are calculated here.                                  *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real rad,re,dlat,dlon,cff
      integer i,j,k
      character*5 field
      rad=0.01745329
      re=6371.E3
C
      write(6,'(/,'' Read grid and initial conditions '',/)')
C
C--- 1D ---
      read(40,'(a5)') field
      write(6,'(a5)') field
       read(40,'(8E12.5)') z
      read(40,'(a5)') field
      write(6,'(a5)') field
       read(40,'(8E12.5)') zz
      read(40,'(a5)') field
      write(6,'(a5)') field
       read(40,'(8E12.5)') dz
      read(40,'(a5)') field
      write(6,'(a5)') field
       read(40,'(8E12.5)') dzz
C--- 2D ---
      read(40,'(a5)') field
      write(6,'(a5)') field
       read(40,'(8E12.5)') east_e
      read(40,'(a5)') field
      write(6,'(a5)') field
       read(40,'(8E12.5)') north_e
      read(40,'(a5)') field
      write(6,'(a5)') field
       read(40,'(8E12.5)') h
C--- 3D ---
      read(40,'(a5)') field
      write(6,'(a5)') field
       read(40,'(8E12.5)') t
      read(40,'(a5)') field
      write(6,'(a5)') field
       read(40,'(8E12.5)') s
      read(40,'(a5)') field
      write(6,'(a5)') field
       read(40,'(8E12.5)') rmean
C--- Constant wind stress read here
C (for time dep. read in loop 9000 & interpolate in time)
      read(40,'(a5)') field
      write(6,'(a5)') field
       read(40,'(8E12.5)') wusurf
      read(40,'(a5)') field
      write(6,'(a5)') field
       read(40,'(8E12.5)') wvsurf
C
C --- print vertical grid distribution
C
      write(6,2)
    2 format(/2x,'k',7x,'z',9x,'zz',9x,'dz',9x,'dzz',/)
      write(6,'(''  '',/)')
      do k=1,kb
        write(6,3) k,z(k),zz(k),dz(k),dzz(k)
    3   format((' ',i5,4f10.3))
      end do
      write(6,'(''  '',//)')
C
C --- calc. surface & lateral BC from climatology
C
        do j=1,jm
          do i=1,im
             tsurf(i,j)=t(i,j,1)
             ssurf(i,j)=s(i,j,1)
            do k=1,kb
              tclim(i,j,k)=t(i,j,k)
              sclim(i,j,k)=s(i,j,k)
            end do
          end do
        end do
C
C                    --- EAST & WEST BCs ---
        do j=1,jm
              ele(j)=0.
              elw(j)=0.
C --- other vel. BCs (fixed in time) can be specified here
              uabe(j)=0.
              uabw(j)=0.
            do k=1,kb
              ubw(j,k)=0.
              ube(j,k)=0.
              tbw(j,k)=tclim(1,j,k)
              sbw(j,k)=sclim(1,j,k)
              tbe(j,k)=tclim(im,j,k)
              sbe(j,k)=sclim(im,j,k)
            end do
        end do
C                    --- NORTH & SOUTH BCs ---
        do i=1,im
              els(i)=0.
              eln(i)=0.
              vabs(i)=0.
              vabn(i)=0.
            do k=1,kb
              vbs(i,k)=0.
              vbn(i,k)=0.
              tbs(i,k)=tclim(i,1,k)
              sbs(i,k)=sclim(i,1,k)
              tbn(i,k)=tclim(i,jm,k)
              sbn(i,k)=sclim(i,jm,k)
            end do
        end do
C
C     Set initial conditions:
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            tb(i,j,k)=t(i,j,k)
            sb(i,j,k)=s(i,j,k)
            ub(i,j,k)=0.
            vb(i,j,k)=0.
          end do
        end do
      end do
C
      call dens(sb,tb,rho)
C
C --- calc. Curiolis Parameter
C
        do j=1,jm
          do i=1,im
            cor(i,j)=2.*7.29E-5*sin(north_e(i,j)*rad)
            aam2d(i,j)=aam(i,j,1)
            elb(i,j)=0.
            etb(i,j)=0.
            dt(i,j)=h(i,j)
          end do
        end do
C
        do j=1,jm
          do i=2,im-1
            dx(i,j)=0.5*rad*re*sqrt(((east_e(i+1,j)-east_e(i-1,j))
     1 *cos(north_e(i,j)*rad))**2+(north_e(i+1,j)-north_e(i-1,j))**2)
          end do
            dx(1,j)=dx(2,j)
            dx(im,j)=dx(im-1,j)
        end do
C
        do i=1,im
          do j=2,jm-1
            dy(i,j)=0.5*rad*re*sqrt(((east_e(i,j+1)-east_e(i,j-1))
     1 *cos(north_e(i,j)*rad))**2+(north_e(i,j+1)-north_e(i,j-1))**2)
          end do
            dy(i,1)=dy(i,2)
            dy(i,jm)=dy(i,jm-1)
        end do
C
C     Calculate areas and masks:
C
      call areas_masks
C
        do j=1,jm
          do i=1,im
             tsurf(i,j)=t(i,j,1)
             ssurf(i,j)=s(i,j,1)
            do k=1,kb
              tclim(i,j,k)=t(i,j,k)
              sclim(i,j,k)=s(i,j,k)
            end do
          end do
        end do
C
C                    --- EAST & WEST BCs ---
        do j=1,jm
            ele(j)=0.
            elw(j)=0.
C --- other vel. BCs (fixed in time) can be specified here
            do k=1,kb
              ubw(j,k)=0                ! RWND
              ube(j,k)=0                !
              uabe(j) =0                ! RWND
              uabw(j) =0                !  Calculate vertical averages
              tbw(j,k)=tclim(1,j,k)
              sbw(j,k)=sclim(1,j,k)
              tbe(j,k)=tclim(im,j,k)
              sbe(j,k)=sclim(im,j,k)
            end do
        end do
C                    --- NORTH & SOUTH BCs ---
        do i=1,im
              els(i)=0. !elb(i,1)   ! RWND
              eln(i)=0. !elb(i,jm)  !
            do k=1,kb
              vbs(i,k)=0                ! RWND
              vbn(i,k)=0                !
              vabs(j) =0                ! RWND
              vabn(j) =0                !  Calculate vertical means
              tbs(i,k)=tclim(i,1,k)
              sbs(i,k)=sclim(i,1,k)
              tbn(i,k)=tclim(i,jm,k)
              sbn(i,k)=sclim(i,jm,k)
            end do
        end do
C
C     Set initial conditions:
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            tb(i,j,k)=t(i,j,k)
            sb(i,j,k)=s(i,j,k)
            ub(i,j,k)=u(i,j,k)
            vb(i,j,k)=v(i,j,k)
          end do
        end do
      end do
C
C --- the following grids are needed only for netcdf plotting
C
C     Corner of cell points:
C
      do j=2,jm
        do i=2,im
          east_c(i,j)=(east_e(i,j)+east_e(i-1,j)
     $                  +east_e(i,j-1)+east_e(i-1,j-1))/4.e0
          north_c(i,j)=(north_e(i,j)+north_e(i-1,j)
     $                   +north_e(i,j-1)+north_e(i-1,j-1))/4.e0
        end do
      end do
C
C
C     Extrapolate ends (approx.):
C
      do i=2,im
        east_c(i,1)=2.*east_c(i,2)-east_c(i,3)
        north_c(i,1)=2.*north_c(i,2)-north_c(i,3)
      end do
        east_c(1,1)=2.*east_c(2,1)-east_c(3,1)
C
      do j=2,jm
        east_c(1,j)=2.*east_c(2,j)-east_c(3,j)
        north_c(1,j)=2.*north_c(2,j)-north_c(3,j)
      end do
        north_c(1,1)=2.*north_c(1,2)-north_c(1,3)
C
C     u-points:
C
      do j=1,jm-1
        do i=1,im
          east_u(i,j)=(east_c(i,j)+east_c(i,j+1))/2.e0
          north_u(i,j)=(north_c(i,j)+north_c(i,j+1))/2.e0
        end do
      end do
C
C     Extrapolate ends:
C
      do i=1,im
        east_u(i,jm)=(east_c(i,jm)*3.e0-east_c(i,jm-1))/2.e0
        north_u(i,jm)=(north_c(i,jm)*3.e0-north_c(i,jm-1))/2.e0
      end do
C
C     v-points:
C
      do j=1,jm
        do i=1,im-1
          east_v(i,j)=(east_c(i,j)+east_c(i+1,j))/2.e0
          north_v(i,j)=(north_c(i,j)+north_c(i+1,j))/2.e0
        end do
      end do
C
C     Extrapolate ends:
C
      do j=1,jm
        east_v(im,j)=(east_c(im,j)*3.e0-east_c(im-1,j))/2.e0
        north_v(im,j)=(north_c(im,j)*3.e0-north_c(im-1,j))/2.e0
      end do
C
C     rot is the angle (radians, anticlockwise) of the i-axis relative
C     to east, averaged to a cell centre: (only needed for CDF plotting)
C
      do j=1,jm
        do i=1,im-1
          rot(i,j)=0.
          dlat=north_e(i+1,j)-north_e(i,j)
          dlon= east_e(i+1,j)- east_e(i,j)
           if(dlon.ne.0.) rot(i,j)=atan(dlat/dlon)
        end do
       rot(im,j)=rot(im-1,j)
      end do
C
C     Set lateral boundary conditions, for use in subroutine bcond
C     set all=0 for closed BCs.
C     Values=0 for vel BC only, =1 is combination of vel+elev.
      rfe=0.e0
      rfw=0.e0
      rfn=0.e0
      rfs=0.e0
C
      return
      end
C
      subroutine profq(sm,sh,dh,cc)
C **********************************************************************
C *                                        Updated: Sep. 24, 2003      *
C * FUNCTION    :  Solves for q2 (twice the turbulent kinetic energy), *
C *                q2l (q2 x turbulent length scale), km (vertical     *
C *                kinematic viscosity) and kh (vertical kinematic     *
C *                diffusivity), using a simplified version of the     *
C *                level 2 1/2 model of Mellor and Yamada (1982).      *
C * In this version, the Craig-Banner sub-model whereby breaking wave  *
C * tke is injected into the surface is included. However, we use an   *
C * analytical solution to the near surface tke equation to solve for  *
C * q2 at the surface giving the same result as C-B diffusion. The new *
C * scheme is simpler and more robust than the latter scheme.          *
C *                                                                    *
C * References                                                         *
C *   Craig, P. D. and M. L. Banner, Modeling wave-enhanced turbulence *
C *     in the ocean surface layer. J. Phys. Oceanogr., 24, 2546-2559, *
C *     1994.                                                          *
C *   Ezer, T., On the seasonal mixed-layer simulated by a basin-scale *
C *     ocean model and the Mellor-Yamada turbulence scheme,           *
C *     J. Geophys. Res., 105(C7), 16,843-16,855, 2000.                *
C *   Mellor, G.L. and T. Yamada, Development of a turbulence          *
C *     closure model for geophysical fluid fluid problems,            *
C *     Rev. Geophys. Space Phys., 20, 851-875, 1982.                  *
C *   Mellor, G. L., One-dimensional, ocean surface layer modeling,    *
C *     a problem and a solution. J. Phys. Oceanogr., 31(3), 790-809,  *
C *     2001.                                                          *
C *   Mellor, G.L. and A. Blumberg, Wave breaking and ocean surface    *
C *     thermal response, J. Phys. Oceanogr., 2003.                    *
C *   Stacey, M. W., Simulations of the wind-forced near-surface       *
C *     circulation in Knight Inlet: a parameterization of the         *
C *     roughness length. J. Phys. Oceanogr., 29, 1363-1367, 1999.     *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real sm(im,jm,kb),sh(im,jm,kb),cc(im,jm,kb)
      real gh(im,jm,kb),boygr(im,jm,kb),dh(im,jm),stf(im,jm,kb)
      real prod(im,jm,kb),kn(im,jm,kb)
      real a1,a2,b1,b2,c1
      real coef1,coef2,coef3,coef4,coef5
      real const1,e1,e2,ghc
      real p,sef,sp,tp
      real l0(im,jm)
      real cbcnst,surfl,shiw
      real utau2, df0,df1,df2
C
      integer i,j,k,ki
C
      equivalence (prod,kn)
C
      data a1,b1,a2,b2,c1/0.92e0,16.6e0,0.74e0,10.1e0,0.08e0/
      data e1/1.8e0/,e2/1.33e0/
      data sef/1.e0/
      data cbcnst/100./surfl/2.e5/shiw/0.0/
C
      do j=1,jm
        do i=1,im
          dh(i,j)=h(i,j)+etf(i,j)
        end do
      end do
C
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            a(i,j,k)=-dti2*(kq(i,j,k+1)+kq(i,j,k)+2.e0*umol)*.5e0
     $                /(dzz(k-1)*dz(k)*dh(i,j)*dh(i,j))
            c(i,j,k)=-dti2*(kq(i,j,k-1)+kq(i,j,k)+2.e0*umol)*.5e0
     $                /(dzz(k-1)*dz(k-1)*dh(i,j)*dh(i,j))
          end do
        end do
      end do
C
C-----------------------------------------------------------------------
C
C     The following section solves the equation:
C
C       dti2*(kq*q2')' - q2*(2.*dti2*dtef+1.) = -q2b
C
C     Surface and bottom boundary conditions:
C
      const1=(16.6e0**(2.e0/3.e0))*sef
C
C initialize fields that are not calculated on all boundaries
C but are later used there
      do i=1,im
        ee(i,jm,1)=0.
        gg(i,jm,1)=0.
        l0(i,jm)=0.
      end do
      do j=1,jm
        ee(im,j,1)=0.
        gg(im,j,1)=0.
        l0(im,j)=0.
      end do
      do i=1,im
      do j=1,jm
       do k=2,kbm1
        prod(i,j,k)=0.
       end do
      end do
      end do
C
      do j=1,jmm1
        do i=1,imm1
          utau2=sqrt((.5e0*(wusurf(i,j)+wusurf(i+1,j)))**2
     $                  +(.5e0*(wvsurf(i,j)+wvsurf(i,j+1)))**2)
C Wave breaking energy- a variant of Craig & Banner (1994)
C see Mellor and Blumberg, 2003.
          ee(i,j,1)=0.e0
          gg(i,j,1)=(15.8*cbcnst)**(2./3.)*utau2
C Surface length scale following Stacey (1999).
          l0(i,j)=surfl*utau2/grav
C
          uf(i,j,kb)=sqrt((.5e0*(wubot(i,j)+wubot(i+1,j)))**2
     $                   +(.5e0*(wvbot(i,j)+wvbot(i,j+1)))**2)*const1
        end do
      end do
C
C    Calculate speed of sound squared:
C
      do k=1,kbm1
        do j=1,jm
          do i=1,im
            tp=t(i,j,k)+tbias
            sp=s(i,j,k)+sbias
C
C     Calculate pressure in units of decibars:
C
            p=grav*rhoref*(-zz(k)* h(i,j))*1.e-4
            cc(i,j,k)=1449.1e0+.00821e0*p+4.55e0*tp -.045e0*tp**2
     $                 +1.34e0*(sp-35.0e0)
            cc(i,j,k)=cc(i,j,k)
     $                 /sqrt((1.e0-.01642e0*p/cc(i,j,k))
     $                   *(1.e0-0.40e0*p/cc(i,j,k)**2))
          end do
        end do
      end do
C
C     Calculate buoyancy gradient:
C
C
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            q2b(i,j,k)=abs(q2b(i,j,k))
            q2lb(i,j,k)=abs(q2lb(i,j,k))
            boygr(i,j,k)=grav*(rho(i,j,k-1)-rho(i,j,k))
     $                    /(dzz(k-1)* h(i,j))
C *** NOTE: comment out next line if dens does not include pressure
     $      +(grav**2)*2.e0/(cc(i,j,k-1)**2+cc(i,j,k)**2)
          end do
        end do
      end do
C
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            l(i,j,k)=abs(q2lb(i,j,k)/q2b(i,j,k))
            if(z(k).gt.-0.5) l(i,j,k)=max(l(i,j,k),kappa*l0(i,j))
            gh(i,j,k)=(l(i,j,k)**2)*boygr(i,j,k)/q2b(i,j,k)
            gh(i,j,k)=min(gh(i,j,k),.028e0)
          end do
        end do
      end do
C
      do j=1,jm
        do i=1,im
          l(i,j,1)=kappa*l0(i,j)
          l(i,j,kb)=0.e0
          gh(i,j,1)=0.e0
          gh(i,j,kb)=0.e0
        end do
      end do
C
C    Calculate production of turbulent kinetic energy:
C
      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            prod(i,j,k)=km(i,j,k)*.25e0*sef
     $                   *((u(i,j,k)-u(i,j,k-1)
     $                      +u(i+1,j,k)-u(i+1,j,k-1))**2
     $                     +(v(i,j,k)-v(i,j,k-1)
     $                      +v(i,j+1,k)-v(i,j+1,k-1))**2)
     $                   /(dzz(k-1)*dh(i,j))**2
C   Add shear due to internal wave field
     $             -shiw*km(i,j,k)*boygr(i,j,k)
            prod(i,j,k)=prod(i,j,k)+kh(i,j,k)*boygr(i,j,k)
          end do
        end do
      end do
C
C  NOTE: Richardson # dep. dissipation correction (Mellor, 2001; Ezer, 2000),
C  depends on ghc the critical number (empirical -6 to -2) to increase mixing.
      ghc=-6.0e0
      do k=1,kb
        do j=1,jm
          do i=1,im
            stf(i,j,k)=1.e0
C It is unclear yet if diss. corr. is needed when surf. waves are included.
c           if(gh(i,j,k).lt.0.e0)
c    $        stf(i,j,k)=1.0e0-0.9e0*(gh(i,j,k)/ghc)**1.5e0
c           if(gh(i,j,k).lt.ghc) stf(i,j,k)=0.1e0
            dtef(i,j,k)=sqrt(abs(q2b(i,j,k)))*stf(i,j,k)
     $                   /(b1*l(i,j,k)+small)
          end do
        end do
      end do
C
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            gg(i,j,k)=1.e0/(a(i,j,k)+c(i,j,k)*(1.e0-ee(i,j,k-1))
     $                      -(2.e0*dti2*dtef(i,j,k)+1.e0))
            ee(i,j,k)=a(i,j,k)*gg(i,j,k)
            gg(i,j,k)=(-2.e0*dti2*prod(i,j,k)+c(i,j,k)*gg(i,j,k-1)
     $                 -uf(i,j,k))*gg(i,j,k)
          end do
        end do
      end do
C
      do k=1,kbm1
        ki=kb-k
        do j=1,jm
          do i=1,im
            uf(i,j,ki)=ee(i,j,ki)*uf(i,j,ki+1)+gg(i,j,ki)
          end do
        end do
      end do
C
C-----------------------------------------------------------------------
C
C     The following section solves the equation:
C
C       dti2(kq*q2l')' - q2l*(dti2*dtef+1.) = -q2lb
C
      do j=1,jm
        do i=1,im
          ee(i,j,2)=0.e0
          gg(i,j,2)=0.e0
          vf(i,j,kb)=0.e0
        end do
      end do
C
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            dtef(i,j,k)=dtef(i,j,k)
     $                   *(1.e0+e2*((1.e0/abs(z(k)-z(1))
     $                               +1.e0/abs(z(k)-z(kb)))
     $                                *l(i,j,k)/(dh(i,j)*kappa))**2)
            gg(i,j,k)=1.e0/(a(i,j,k)+c(i,j,k)*(1.e0-ee(i,j,k-1))
     $                      -(dti2*dtef(i,j,k)+1.e0))
            ee(i,j,k)=a(i,j,k)*gg(i,j,k)
            gg(i,j,k)=(dti2*(-prod(i,j,k)*l(i,j,k)*e1)
     $                 +c(i,j,k)*gg(i,j,k-1)-vf(i,j,k))*gg(i,j,k)
          end do
        end do
      end do
C
      do k=1,kb-2
        ki=kb-k
        do j=1,jm
          do i=1,im
            vf(i,j,ki)=ee(i,j,ki)*vf(i,j,ki+1)+gg(i,j,ki)
          end do
        end do
      end do
C
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            if(uf(i,j,k).le.small.or.vf(i,j,k).le.small) then
              uf(i,j,k)=small
              vf(i,j,k)=0.1*dt(i,j)*small
            endif
          end do
        end do
      end do
C
C-----------------------------------------------------------------------
C
C     The following section solves for km and kh:
C
      coef4=18.e0*a1*a1+9.e0*a1*a2
      coef5=9.e0*a1*a2
C
C     Note that sm and sh limit to infinity when gh approaches 0.0288:
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            coef1=a2*(1.e0-6.e0*a1/b1*stf(i,j,k))
            coef2=3.e0*a2*b2/stf(i,j,k)+18.e0*a1*a2
            coef3=a1*(1.e0-3.e0*c1-6.e0*a1/b1*stf(i,j,k))
            sh(i,j,k)=coef1/(1.e0-coef2*gh(i,j,k))
            sm(i,j,k)=coef3+sh(i,j,k)*coef4*gh(i,j,k)
            sm(i,j,k)=sm(i,j,k)/(1.e0-coef5*gh(i,j,k))
          end do
        end do
      end do
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            kn(i,j,k)=l(i,j,k)*sqrt(abs(q2(i,j,k)))
            kq(i,j,k)=(kn(i,j,k)*.41e0*sh(i,j,k)+kq(i,j,k))*.5e0
            km(i,j,k)=(kn(i,j,k)*sm(i,j,k)+km(i,j,k))*.5e0
            kh(i,j,k)=(kn(i,j,k)*sh(i,j,k)+kh(i,j,k))*.5e0
          end do
        end do
      end do
C cosmetics: make boundr. values as interior
C (even if not used, printout otherwise may show strange values)
      do k=1,kb
        do i=1,im
           km(i,jm,k)=km(i,jmm1,k)*fsm(i,jm)
           kh(i,jm,k)=kh(i,jmm1,k)*fsm(i,jm)
           km(i,1,k)=km(i,2,k)*fsm(i,1)
           kh(i,1,k)=kh(i,2,k)*fsm(i,1)
        end do
        do j=1,jm
           km(im,j,k)=km(imm1,j,k)*fsm(im,j)
           kh(im,j,k)=kh(imm1,j,k)*fsm(im,j)
           km(1,j,k)=km(2,j,k)*fsm(1,j)
           kh(1,j,k)=kh(2,j,k)*fsm(1,j)
        end do
      end do
C
      return
C
      end
C
c ---------------------------------------------------------------------
C
      subroutine proft(f,wfsurf,fsurf,nbc,dh)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Solves for vertical diffusion of temperature and    *
C *                salinity using method described by Richmeyer and    *
C *                Morton.                                             *
C *                                                                    *
C *                Irradiance parameters are from Paulson and Simpson. *
C *                                                                    *
C *                See:                                                *
C *                                                                    *
C *                Richtmeyer R.D., and K.W. Morton, 1967. Difference  *
C *                  Methods for Initial-Value Problems, 2nd edition,  *
C *                  Interscience, New York, 198-201.                  *
C *                                                                    *
C *                Paulson, C. A., and J. Simpson, 1977: Irradiance    *
C *                  measurements in the upper ocean, J. Phys.         *
C *                  Oceanogr., 7, 952-956.                            *
C *                                                                    *
C *                NOTES:                                              *
C *                                                                    *
C *                (1) wfsurf and swrad are negative values when water *
C *                    column is warming or salt is being added.       *
C *                                                                    *
C *                (2) nbc may only be 1 and 3 for salinity.           *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real f(im,jm,kb),wfsurf(im,jm)
      real fsurf(im,jm),dh(im,jm)
      integer nbc
      real rad(im,jm,kb),r(5),ad1(5),ad2(5)
      integer i,j,k,ki
C
C-----------------------------------------------------------------------
C
C     Irradiance parameters after Paulson and Simpson:
C
C       ntp               1      2       3       4       5
C   Jerlov type           i      ia      ib      ii     iii
C
      data r   /       .58e0,  .62e0,  .67e0,  .77e0,  .78e0 /
      data ad1 /       .35e0,  .60e0,  1.0e0,  1.5e0,  1.4e0 /
      data ad2 /       23.e0,  20.e0,  17.e0,  14.e0,  7.9e0 /
C
C-----------------------------------------------------------------------
C
C     Surface boundary condition:
C
C       nbc   prescribed    prescribed   short wave
C             temperature      flux      penetration
C             or salinity               (temperature
C                                           only)
C
C        1        no           yes           no
C        2        no           yes           yes
C        3        yes          no            no
C        4        yes          no            yes
C
C     NOTE that only 1 and 3 are allowed for salinity.
C
C-----------------------------------------------------------------------
C
C     The following section solves the equation:
C
C       dti2*(kh*f')'-f=-fb
C
      do j=1,jm
        do i=1,im
          dh(i,j)=h(i,j)+etf(i,j)
        end do
      end do
C
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            a(i,j,k-1)=-dti2*(kh(i,j,k)+umol)
     $                  /(dz(k-1)*dzz(k-1)*dh(i,j)*dh(i,j))
            c(i,j,k)=-dti2*(kh(i,j,k)+umol)
     $                  /(dz(k)*dzz(k-1)*dh(i,j)*dh(i,j))
          end do
        end do
      end do
C
C     Calculate penetrative radiation. At the bottom any unattenuated
C     radiation is deposited in the bottom layer:
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            rad(i,j,k)=0.e0
          end do
        end do
      end do
C
      if(nbc.eq.2.or.nbc.eq.4) then
C
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              rad(i,j,k)=swrad(i,j)
     $                    *(r(ntp)*exp(z(k)*dh(i,j)/ad1(ntp))
     $                      +(1.e0-r(ntp))*exp(z(k)*dh(i,j)/ad2(ntp)))
            end do
          end do
        end do
C
      endif
C
      if(nbc.eq.1) then
C
        do j=1,jm
          do i=1,im
            ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.e0)
            gg(i,j,1)=-dti2*wfsurf(i,j)/(-dz(1)*dh(i,j))-f(i,j,1)
            gg(i,j,1)=gg(i,j,1)/(a(i,j,1)-1.e0)
          end do
        end do
C
      else if(nbc.eq.2) then
C
        do j=1,jm
          do i=1,im
            ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.e0)
            gg(i,j,1)=dti2*(wfsurf(i,j)+rad(i,j,1)-rad(i,j,2))
     $                 /(dz(1)*dh(i,j))
     $                   -f(i,j,1)
            gg(i,j,1)=gg(i,j,1)/(a(i,j,1)-1.e0)
          end do
        end do
C
      else if(nbc.eq.3.or.nbc.eq.4) then
C
        do j=1,jm
          do i=1,im
            ee(i,j,1)=0.e0
            gg(i,j,1)=fsurf(i,j)
          end do
        end do
C
      endif
C
      do k=2,kbm2
        do j=1,jm
          do i=1,im
            gg(i,j,k)=1.e0/(a(i,j,k)+c(i,j,k)*(1.e0-ee(i,j,k-1))-1.e0)
            ee(i,j,k)=a(i,j,k)*gg(i,j,k)
            gg(i,j,k)=(c(i,j,k)*gg(i,j,k-1)-f(i,j,k)
     $                 +dti2*(rad(i,j,k)-rad(i,j,k+1))
     $                   /(dh(i,j)*dz(k)))
     $                 *gg(i,j,k)
          end do
        end do
      end do
C
C     Bottom adiabatic boundary condition:
C
      do j=1,jm
        do i=1,im
          f(i,j,kbm1)=(c(i,j,kbm1)*gg(i,j,kbm2)-f(i,j,kbm1)
     $                 +dti2*(rad(i,j,kbm1)-rad(i,j,kb))
     $                   /(dh(i,j)*dz(kbm1)))
     $                 /(c(i,j,kbm1)*(1.e0-ee(i,j,kbm2))-1.e0)
        end do
      end do
C
      do k=2,kbm1
        ki=kb-k
        do j=1,jm
          do i=1,im
          f(i,j,ki)=(ee(i,j,ki)*f(i,j,ki+1)+gg(i,j,ki))
          end do
        end do
      end do
C
      return
C
      end
C
      subroutine profu
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Solves for vertical diffusion of x-momentum using   *
C *                method described by Richmeyer and Morton.           *
C *                                                                    *
C *                See:                                                *
C *                                                                    *
C *                Richtmeyer R.D., and K.W. Morton, 1967. Difference  *
C *                  Methods for Initial-Value Problems, 2nd edition,  *
C *                  Interscience, New York, 198-201.                  *
C *                                                                    *
C *                NOTE that wusurf has the opposite sign to the wind  *
C *                speed.                                              *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
      real dh(im,jm)
      integer i,j,k,ki
C
C     The following section solves the equation:
C
C       dti2*(km*u')'-u=-ub
C
      do j=1,jm
        do i=1,im
          dh(i,j)=1.e0
        end do
      end do
C
      do j=2,jm
        do i=2,im
          dh(i,j)=(h(i,j)+etf(i,j)+h(i-1,j)+etf(i-1,j))*.5e0
        end do
      end do
C
      do k=1,kb
        do j=2,jm
          do i=2,im
            c(i,j,k)=(km(i,j,k)+km(i-1,j,k))*.5e0
          end do
        end do
      end do
C
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            a(i,j,k-1)=-dti2*(c(i,j,k)+umol)
     $                  /(dz(k-1)*dzz(k-1)*dh(i,j)*dh(i,j))
            c(i,j,k)=-dti2*(c(i,j,k)+umol)
     $                /(dz(k)*dzz(k-1)*dh(i,j)*dh(i,j))
          end do
        end do
      end do
C
      do j=1,jm
        do i=1,im
          ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.e0)
          gg(i,j,1)=(-dti2*wusurf(i,j)/(-dz(1)*dh(i,j))
     $               -uf(i,j,1))
     $               /(a(i,j,1)-1.e0)
        end do
      end do
C
      do k=2,kbm2
        do j=1,jm
          do i=1,im
            gg(i,j,k)=1.e0/(a(i,j,k)+c(i,j,k)*(1.e0-ee(i,j,k-1))-1.e0)
            ee(i,j,k)=a(i,j,k)*gg(i,j,k)
            gg(i,j,k)=(c(i,j,k)*gg(i,j,k-1)-uf(i,j,k))*gg(i,j,k)
          end do
        end do
      end do
C
      do j=2,jmm1
        do i=2,imm1
          tps(i,j)=0.5e0*(cbc(i,j)+cbc(i-1,j))
     $              *sqrt(ub(i,j,kbm1)**2
     $                +(.25e0*(vb(i,j,kbm1)+vb(i,j+1,kbm1)
     $                         +vb(i-1,j,kbm1)+vb(i-1,j+1,kbm1)))**2)
          uf(i,j,kbm1)=(c(i,j,kbm1)*gg(i,j,kbm2)-uf(i,j,kbm1))
     $                  /(tps(i,j)*dti2/(-dz(kbm1)*dh(i,j))-1.e0
     $                    -(ee(i,j,kbm2)-1.e0)*c(i,j,kbm1))
          uf(i,j,kbm1)=uf(i,j,kbm1)*dum(i,j)
        end do
      end do
C

      do k=2,kbm1
        ki=kb-k
        do j=2,jmm1
          do i=2,imm1
            uf(i,j,ki)=(ee(i,j,ki)*uf(i,j,ki+1)+gg(i,j,ki))*dum(i,j)
          end do
        end do
      end do
C
      do j=2,jmm1
        do i=2,imm1
          wubot(i,j)=-tps(i,j)*uf(i,j,kbm1)
        end do
      end do
C
      return
C
      end
C
      subroutine profv
C **********************************************************************
C                                                                      *
C * FUNCTION    :  Solves for vertical diffusion of y-momentum using   *
C *                method described by Richmeyer and Morton.           *
C *                                                                    *
C *                See:                                                *
C *                                                                    *
C *                Richtmeyer R.D., and K.W. Morton, 1967. Difference  *
C *                  Methods for Initial-Value Problems, 2nd edition,  *
C *                  Interscience, New York, 198-201.                  *
C *                                                                    *
C *                NOTE that wvsurf has the opposite sign to the wind  *
C *                speed.                                              *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
      real dh(im,jm)
      integer i,j,k,ki
C
C     The following section solves the equation:
C
C       dti2*(km*u')'-u=-ub
C
      do j=1,jm
        do i=1,im
          dh(i,j)=1.e0
        end do
      end do
C
      do j=2,jm
        do i=2,im
          dh(i,j)=.5e0*(h(i,j)+etf(i,j)+h(i,j-1)+etf(i,j-1))
        end do
      end do
C
      do k=1,kb
        do j=2,jm
          do i=2,im
            c(i,j,k)=(km(i,j,k)+km(i,j-1,k))*.5e0
          end do
        end do
      end do
C
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            a(i,j,k-1)=-dti2*(c(i,j,k)+umol)
     $                  /(dz(k-1)*dzz(k-1)*dh(i,j)*dh(i,j))
            c(i,j,k)=-dti2*(c(i,j,k)+umol)
     $                /(dz(k)*dzz(k-1)*dh(i,j)*dh(i,j))
          end do
        end do
      end do
C
      do j=1,jm
        do i=1,im
          ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.e0)
          gg(i,j,1)=(-dti2*wvsurf(i,j)/(-dz(1)*dh(i,j))-vf(i,j,1))
     $               /(a(i,j,1)-1.e0)
        end do
      end do
C
      do k=2,kbm2
        do j=1,jm
          do i=1,im
            gg(i,j,k)=1.e0/(a(i,j,k)+c(i,j,k)*(1.e0-ee(i,j,k-1))-1.e0)
            ee(i,j,k)=a(i,j,k)*gg(i,j,k)
            gg(i,j,k)=(c(i,j,k)*gg(i,j,k-1)-vf(i,j,k))*gg(i,j,k)
          end do
        end do
      end do
C
      do j=2,jmm1
        do i=2,imm1
          tps(i,j)=0.5e0*(cbc(i,j)+cbc(i,j-1))
     $              *sqrt((.25e0*(ub(i,j,kbm1)+ub(i+1,j,kbm1)
     $                            +ub(i,j-1,kbm1)+ub(i+1,j-1,kbm1)))**2
     $                    +vb(i,j,kbm1)**2)
          vf(i,j,kbm1)=(c(i,j,kbm1)*gg(i,j,kbm2)-vf(i,j,kbm1))
     $                  /(tps(i,j)*dti2/(-dz(kbm1)*dh(i,j))-1.e0
     $                    -(ee(i,j,kbm2)-1.e0)*c(i,j,kbm1))
          vf(i,j,kbm1)=vf(i,j,kbm1)*dvm(i,j)
        end do
      end do
C
      do k=2,kbm1
        ki=kb-k
        do j=2,jmm1
          do i=2,imm1
            vf(i,j,ki)=(ee(i,j,ki)*vf(i,j,ki+1)+gg(i,j,ki))*dvm(i,j)
          end do
        end do
      end do
C
      do j=2,jmm1
        do i=2,imm1
          wvbot(i,j)=-tps(i,j)*vf(i,j,kbm1)
        end do
      end do
C
      return
C
      end
C
      subroutine seamount
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Sets up for seamount problem.                       *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real delh,delx,elejmid,elwjmid,ra,vel
      integer i,j,k
C
C     Set delh > 1.0 for an island or delh < 1.0 for a seamount:
C
      delh=0.9e0
C
C     Grid size:
C
      delx=8000.e0
C
C     Radius island or seamount:
C
      ra=25000.e0
C
C     Current velocity:
C
      vel=0.2e0
C
C     Set up grid dimensions, areas of free surface cells, and
C     Coriolis parameter:
C
      do j=1,jm
        do i=1,im
C
C     For constant grid size:
C
C         dx(i,j)=delx
C         dy(i,j)=delx
C
C     For variable grid size:
C
          dx(i,j)=delx-delx*sin(pi*float(i)/float(im))/2.e0
          dy(i,j)=delx-delx*sin(pi*float(j)/float(jm))/2.e0
C
          cor(i,j)=1.e-4
C
        end do
      end do
C
C     Calculate horizontal coordinates of grid points and rotation
C     angle.
C
C     NOTE that this is introduced solely for the benefit of any post-
C     processing software, and in order to conform with the requirements
C     of the NetCDF Climate and Forecast (CF) Metadata Conventions.
C
C     There are four horizontal coordinate systems, denoted by the
C     subscripts u, v, e and c ("u" is a u-point, "v" is a v-point,
C     "e" is an elevation point and "c" is a cell corner), as shown
C     below. In addition, "east_*" is an easting and "north_*" is a
C     northing. Hence the coordinates of the "u" points are given by
C     (east_u,north_u).
C
C     Also, if the centre point of the cell shown below is at
C     (east_e(i,j),north_e(i,j)), then (east_u(i,j),north_u(i,j)) are
C     the coordinates of the western of the two "u" points,
C     (east_v(i,j),north_v(i,j)) are the coordinates of the southern of
C     the two "v" points, and (east_c(i,j),north_c(i,j)) are the
C     coordinates of the southwestern corner point of the cell. The
C     southwest corner of the entire grid is at
C     (east_c(1,1),north_c(1,1)).
C
C                      |              |
C                    --c------v-------c--
C                      |              |
C                      |              |
C                      |              |
C                      |              |
C                      u      e       u
C                      |              |
C                      |              |
C                      |              |
C                      |              |
C                    --c------v-------c--
C                      |              |
C
C
C     NOTE that the following calculation of east_c and north_c only
C     works properly for a rectangular grid with east and north aligned
C     with i and j, respectively:
C
      do j=1,jm
        east_c(1,j)=0.e0
        do i=2,im
          east_c(i,j)=east_c(i-1,j)+dx(i-1,j)
        end do
      end do
C
      do i=1,im
        north_c(i,1)=0.e0
        do j=2,jm
          north_c(i,j)=north_c(i,j-1)+dy(i,j-1)
        end do
      end do
C
C     The following works properly for any grid:
C
C     Elevation points:
C
      do j=1,jm-1
        do i=1,im-1
          east_e(i,j)=(east_c(i,j)+east_c(i+1,j)
     $                  +east_c(i,j+1)+east_c(i+1,j+1))/4.e0
          north_e(i,j)=(north_c(i,j)+north_c(i+1,j)
     $                   +north_c(i,j+1)+north_c(i+1,j+1))/4.e0
        end do
      end do
C
C     Extrapolate ends:
C
      do i=1,im-1
        east_e(i,jm)
     $    =((east_c(i,jm)+east_c(i+1,jm))*3.e0
     $       -east_c(i,jm-1)-east_c(i+1,jm-1))/4.e0
        north_e(i,jm)
     $    =((north_c(i,jm)+north_c(i+1,jm))*3.e0
     $       -north_c(i,jm-1)-north_c(i+1,jm-1))/4.e0
      end do
C
      do j=1,jm-1
        east_e(im,j)
     $    =((east_c(im,j)+east_c(im,j+1))*3.e0
     $       -east_c(im-1,j)-east_c(im-1,j+1))/4.e0
        north_e(im,j)
     $    =((north_c(im,j)+north_c(im,j+1))*3.e0
     $       -north_c(im-1,j)-north_c(im-1,j+1))/4.e0
      end do
C
      east_e(im,jm)=east_e(im-1,jm)+east_e(im,jm-1)
     $               -(east_e(im-2,jm)+east_e(im,jm-2))/2.e0
      north_e(im,jm)=north_e(im-1,jm)+north_e(im,jm-1)
     $               -(north_e(im-2,jm)+north_e(im,jm-2))/2.e0
C
C     u-points:
C
      do j=1,jm-1
        do i=1,im
          east_u(i,j)=(east_c(i,j)+east_c(i,j+1))/2.e0
          north_u(i,j)=(north_c(i,j)+north_c(i,j+1))/2.e0
        end do
      end do
C
C     Extrapolate ends:
C
      do i=1,im
        east_u(i,jm)=(east_c(i,jm)*3.e0-east_c(i,jm-1))/2.e0
        north_u(i,jm)=(north_c(i,jm)*3.e0-north_c(i,jm-1))/2.e0
      end do
C
C     v-points:
C
      do j=1,jm
        do i=1,im-1
          east_v(i,j)=(east_c(i,j)+east_c(i+1,j))/2.e0
          north_v(i,j)=(north_c(i,j)+north_c(i+1,j))/2.e0
        end do
      end do
C
C     Extrapolate ends:
C
      do j=1,jm
        east_v(im,j)=(east_c(im,j)*3.e0-east_c(im-1,j))/2.e0
        north_v(im,j)=(north_c(im,j)*3.e0-north_c(im-1,j))/2.e0
      end do
C
C     rot is the angle (radians, anticlockwise) of the i-axis relative
C     to east, averaged to a cell centre:
C
C     (NOTE that the following calculation of rot only works properly
C     for this particular rectangular grid)
C
      do j=1,jm
        do i=1,im
          rot(i,j)=0.e0
        end do
      end do
C
C     Define depth:
C
      do i=1,im
        do j=1,jm
C
          h(i,j)=4500.e0*(1.e0-delh
     $                          *exp(-((east_c(i,j)
     $                                   -east_c((im+1)/2,j))**2
     $                                +(north_c(i,j)
     $                                   -north_c(i,(jm+1)/2))**2)
     $                                /ra**2))
          if(h(i,j).lt.1.e0) h(i,j)=1.e0
C
        end do
      end do
C
C     Close the north and south boundaries to form a channel:
C
      do i=1,im
        h(i,1)=1.e0
        h(i,jm)=1.e0
      end do
C
C     Calculate areas and masks:
C
      call areas_masks
C
C     Adjust bottom topography so that cell to cell variations
C     in h do not exceed parameter slmax:
C
      if(slmax.lt.1.e0) call slpmax
C
C     Set initial conditions:
C
      do k=1,kbm1
        do j=1,jm
          do i=1,im
            tb(i,j,k)=5.e0+15.e0*exp(zz(k)*h(i,j)/1000.e0)-tbias
            sb(i,j,k)=35.e0-sbias
            tclim(i,j,k)=tb(i,j,k)
            sclim(i,j,k)=sb(i,j,k)
            ub(i,j,k)=vel*dum(i,j)
          end do
        end do
      end do
C
C     Initialise uab and vab as necessary
C     (NOTE that these have already been initialised to zero in the
C     main program):
C
      do j=1,jm
        do i=1,im
          uab(i,j)=vel*dum(i,j)
        end do
      end do
C
C     Set surface boundary conditions, e_atmos, vflux, wusurf,
C     wvsurf, wtsurf, wssurf and swrad, as necessary
C     (NOTE:
C      1. These have all been initialised to zero in the main program.
C      2. The temperature and salinity of inflowing water must be
C         defined relative to tbias and sbias.):
C
      do j=1,jm
        do i=1,im
C     No conditions necessary for this problem
        end do
      end do
C
C     Initialise elb, etb, dt and aam2d:
C
      do j=1,jm
        do i=1,im
          elb(i,j)=-e_atmos(i,j)
          etb(i,j)=-e_atmos(i,j)
          dt(i,j)=h(i,j)-e_atmos(i,j)
          aam2d(i,j)=aam(i,j,1)
        end do
      end do
C
      call dens(sb,tb,rho)
C
C     Generated horizontally averaged density field (in this
C     application, the initial condition for density is a function
C     of z (the vertical cartesian coordinate) -- when this is not
C     so, make sure that rmean has been area averaged BEFORE transfer
C     to sigma coordinates):
C
      do k=1,kbm1
        do j=1,jm
          do i=1,im
            rmean(i,j,k)=rho(i,j,k)
          end do
        end do
      end do
C
C     Set lateral boundary conditions, for use in subroutine bcond
C     (in the seamount problem, the east and west boundaries are open,
C     while the south and north boundaries are closed through the
C     specification of the masks fsm, dum and dvm):
C
      rfe=1.e0
      rfw=1.e0
      rfn=1.e0
      rfs=1.e0
C
      do j=2,jmm1
        uabw(j)=uab(2,j)
        uabe(j)=uab(imm1,j)
C
C     Set geostrophically conditioned elevations at the boundaries:
C
        ele(j)=ele(j-1)-cor(imm1,j)*uab(imm1,j)/grav*dy(imm1,j-1)
        elw(j)=elw(j-1)-cor(2,j)*uab(2,j)/grav*dy(2,j-1)
      end do
C
C     Adjust boundary elevations so that they are zero in the middle
C     of the channel:
C
      elejmid=ele(jmm1/2)
      elwjmid=elw(jmm1/2)
      do j=2,jmm1
        ele(j)=(ele(j)-elejmid)*fsm(im,j)
        elw(j)=(elw(j)-elwjmid)*fsm(2,j)
      end do
C
C     Set thermodynamic boundary conditions (for the seamount
C     problem, and other possible applications, lateral thermodynamic
C     boundary conditions are set equal to the initial conditions and
C     are held constant thereafter - users may, of course, create
C     variable boundary conditions):
C
      do k=1,kbm1
C
        do j=1,jm
          tbe(j,k)=tb(im,j,k)
          tbw(j,k)=tb(1,j,k)
          sbe(j,k)=sb(im,j,k)
          sbw(j,k)=sb(1,j,k)
        end do
C
        do i=1,im
          tbn(i,k)=tb(i,jm,k)
          tbs(i,k)=tb(i,1,k)
          sbn(i,k)=sb(i,jm,k)
          sbs(i,k)=sb(i,1,k)
        end do
C
      end do
C
      return
C
      end
C
      subroutine slpmax
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Limits the maximum of:                              *
C *                                                                    *
C *                  <difference of depths>/<sum of depths>            *
C *                                                                    *
C *                for two adjacent cells. The maximum possible value  *
C *                is unity.                                           *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real mean,del
      integer i,j,loop
C
      do loop=1,10
C
C     Sweep right:
C
        do j=2,jm-1
C
          do i=2,im-1
            if(fsm(i,j).ne.0.e0.and.fsm(i+1,j).ne.0.e0) then
              if(abs(h(i+1,j)-h(i,j))/(h(i,j)+h(i+1,j)).ge.slmax) then
                mean=(h(i+1,j)+h(i,j))/2.e0
                del=sign(slmax,h(i+1,j)-h(i,j))
                h(i+1,j)=mean*(1.e0+del)
                h(i,j)=mean*(1.e0-del)
              endif
            endif
          end do
C
C    Sweep left:
C
          do i=im-1,2,-1
            if(fsm(i,j).ne.0.e0.and.fsm(i+1,j).ne.0.e0) then
              if(abs(h(i+1,j)-h(i,j))/(h(i,j)+h(i+1,j)).ge.slmax) then
                mean=(h(i+1,j)+h(i,j))/2.e0
                del=sign(slmax,h(i+1,j)-h(i,j))
                h(i+1,j)=mean*(1.e0+del)
                h(i,j)=mean*(1.e0-del)
              endif
            endif
          end do
C
        end do
C
C   Sweep up:
C
        do i=2,im-1
C
          do j=2,jm-1
            if(fsm(i,j).ne.0.e0.and.fsm(i,j+1).ne.0.e0) then
              if(abs(h(i,j+1)-h(i,j))/(h(i,j)+h(i,j+1)).ge.slmax) then
                mean=(h(i,j+1)+h(i,j))/2.e0
                del=sign(slmax,h(i,j+1)-h(i,j))
                h(i,j+1)=mean*(1.e0+del)
                h(i,j)=mean*(1.e0-del)
              endif
            endif
          end do
C
C   Sweep down:
C
          do j=jm-1,2,-1
            if(fsm(i,j).ne.0.e0.and.fsm(i,j+1).ne.0.e0) then
              if(abs(h(i,j+1)-h(i,j))/(h(i,j)+h(i,j+1)).ge.slmax) then
                mean=(h(i,j+1)+h(i,j))/2.e0
                del=sign(slmax,h(i,j+1)-h(i,j))
                h(i,j+1)=mean*(1.e0+del)
                h(i,j)=mean*(1.e0-del)
              endif
            endif
          end do
C
        end do
C
      end do
C
      return
C
      end
C
      subroutine smol_adif(xmassflux,ymassflux,zwflux,ff,sw)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Calculates the antidiffusive velocity used to       *
C *                reduce the numerical diffusion associated with the  *
C *                upstream differencing scheme.                       *
C *                                                                    *
C *                This is based on a subroutine of Gianmaria Sannino  *
C *                (Inter-university Computing Consortium, Rome, Italy)*
C *                and Vincenzo Artale (Italian National Agency for    *
C *                New Technology and Environment, Rome, Italy),       *
C *                downloaded from the POM FTP site on 1 Nov. 2001.    *
C *                The calculations have been simplified by removing   *
C *                the shock switch option.                            *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real ff(im,jm,kb)
      real xmassflux(im,jm,kb),ymassflux(im,jm,kb),zwflux(im,jm,kb)
      real sw
      real mol,abs_1,abs_2
      real value_min,epsilon
      real udx,u2dt,vdy,v2dt,wdz,w2dt
      integer i,j,k
C
      parameter (value_min=1.e-9,epsilon=1.0e-14)
C
C     Apply temperature and salinity mask:
C
      do k=1,kb
        do i=1,im
          do j=1,jm
            ff(i,j,k)=ff(i,j,k)*fsm(i,j)
          end do
        end do
      end do
C
C     Recalculate mass fluxes with antidiffusion velocity:
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,im
            if(ff(i,j,k).lt.value_min.or.
     $         ff(i-1,j,k).lt.value_min) then
              xmassflux(i,j,k)=0.e0
            else
              udx=abs(xmassflux(i,j,k))
              u2dt=dti2*xmassflux(i,j,k)*xmassflux(i,j,k)*2.e0
     $              /(aru(i,j)*(dt(i-1,j)+dt(i,j)))
              mol=(ff(i,j,k)-ff(i-1,j,k))
     $             /(ff(i-1,j,k)+ff(i,j,k)+epsilon)
              xmassflux(i,j,k)=(udx-u2dt)*mol*sw
              abs_1=abs(udx)
              abs_2=abs(u2dt)
              if(abs_1.lt.abs_2) xmassflux(i,j,k)=0.e0
            end if
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=2,jm
          do i=2,imm1
            if(ff(i,j,k).lt.value_min.or.
     $         ff(i,j-1,k).lt.value_min) then
              ymassflux(i,j,k)=0.e0
            else
             vdy=abs(ymassflux(i,j,k))
             v2dt=dti2*ymassflux(i,j,k)*ymassflux(i,j,k)*2.e0
     $             /(arv(i,j)*(dt(i,j-1)+dt(i,j)))
             mol=(ff(i,j,k)-ff(i,j-1,k))
     $            /(ff(i,j-1,k)+ff(i,j,k)+epsilon)
             ymassflux(i,j,k)=(vdy-v2dt)*mol*sw
             abs_1=abs(vdy)
             abs_2=abs(v2dt)
             if(abs_1.lt.abs_2) ymassflux(i,j,k)=0.e0
            end if
          end do
        end do
      end do
C
      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            if(ff(i,j,k).lt.value_min.or.
     $         ff(i,j,k-1).lt.value_min) then
              zwflux(i,j,k)=0.e0
            else
              wdz=abs(zwflux(i,j,k))
              w2dt=dti2*zwflux(i,j,k)*zwflux(i,j,k)/(dzz(k-1)*dt(i,j))
              mol=(ff(i,j,k-1)-ff(i,j,k))
     $             /(ff(i,j,k)+ff(i,j,k-1)+epsilon)
              zwflux(i,j,k)=(wdz-w2dt)*mol*sw
              abs_1=abs(wdz)
              abs_2=abs(w2dt)
              if(abs_1.lt.abs_2)zwflux(i,j,k)=0.e0
            end if
          end do
        end do
      end do
C
      return
C
      end
C
      subroutine vertvl(xflux,yflux)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Calculates vertical velocity.                       *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real xflux(im,jm,kb),yflux(im,jm,kb)
      integer i,j,k
C
C     Reestablish boundary conditions:
C
      do k=1,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=.25e0*(dy(i,j)+dy(i-1,j))
     $                    *(dt(i,j)+dt(i-1,j))*u(i,j,k)
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=2,jm
          do i=2,im
            yflux(i,j,k)=.25e0*(dx(i,j)+dx(i,j-1))
     $                    *(dt(i,j)+dt(i,j-1))*v(i,j,k)
          end do
        end do
      end do
C
C     NOTE that, if one wishes to include freshwater flux, the
C     surface velocity should be set to vflux(i,j). See also
C     change made to 2-D volume conservation equation which
C     calculates elf.
C
        do j=2,jmm1
          do i=2,imm1
            w(i,j,1)=0.5*(vfluxb(i,j)+vfluxf(i,j))
          end do
        end do
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            w(i,j,k+1)=w(i,j,k)
     $                +dz(k)*((xflux(i+1,j,k)-xflux(i,j,k)
     $                        +yflux(i,j+1,k)-yflux(i,j,k))
     $                        /(dx(i,j)*dy(i,j))
     $                        +(etf(i,j)-etb(i,j))/dti2)
          end do
        end do
      end do
C
      return
C
      end
C
      subroutine ncdf2ic
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Sets up my own problem.                             *
C *                                                                    *
C * This example reads IC from NetCDF files, generated by <deleted>    *
C * Only minimal number of fields are read,                            *
C * while others are calculated here.  TODO: Read as much as we can    *
C *                                                                    *
C **********************************************************************
C
      use netcdf
      implicit none
C
      include 'pom2k.c'
C
      real datr(im, jm, kbm1, 2)
      real datu(imm1, jm, kbm1, 2)
      real datv(im, jmm1, kbm1, 2)
!
      real rad,re,dlat,dlon,cff
      integer i,j,k,mm
      real lom(12)   ! length of month
      data lom /31,28.25,31,30,31,30,31,31,30,31,30,31/
      rad=0.01745329
      re=6371.E3
C
      write(6,'(/,'' Read grid and initial conditions '',/)')
C
C--- 1D ---
      filename = trim(pth_wrk)//trim(pth_grd)//
     $           trim(pfx_dmn)//"roms_lvl.nc"   ! `roms_lvl` is actually `roms_bry` but the new one which has incorrect velocities but has new sigma-levels.
      write(*,*) "\\",trim(filename)
      call check( nf90_open(filename, NF90_NOWRITE, ncid) )
      call check( nf90_inq_varid(ncid, "Cs_w", varid) )
      call check( nf90_get_var(ncid, varid, z(kb:1:-1)) )
      write(*, *) "[O] z retrieved"
      call check( nf90_inq_varid(ncid, "Cs_r", varid) )
      call check( nf90_get_var(ncid, varid, zz(kbm1:1:-1)) )
      zz(kb) = zz(kbm1)+zz(kbm1)-zz(kbm1-1)
      write(*, *) "[O] zz retrieved"
      do k=1,kbm1
        dz(k)  =  z(k)- z(k+1)
        dzz(k) = zz(k)-zz(k+1)
      end do
      dz(kb)  = 0.
      dzz(kb) = 0.
      write(*, *) "[O] z derivatives (dz,dzz) calculated"
      call check( nf90_close(ncid) )

      filename = trim(pth_wrk)//trim(pth_grd)//
     $           trim(pfx_dmn)//"roms_ini.nc"
      write(*,*) "\\",trim(filename)
      call check( nf90_open(filename, NF90_NOWRITE, ncid) )
C--- 3D ---
      call check( nf90_inq_varid(ncid, "temp", varid) )
      call check( nf90_get_var(ncid, varid, t(:,:,kbm1:1:-1)) )
      t(:,:,kb) = t(:,:,kbm1)
      write(*, *) "[O] potential temperature retrieved"
      call check( nf90_inq_varid(ncid, "salt", varid) )
      call check( nf90_get_var(ncid, varid, s(:,:,kbm1:1:-1)) )
      s(:,:,kb) = s(:,:,kbm1)
      write(*, *) "[O] salinity retrieved"
      call check( nf90_close(ncid) )
      
      filename = trim(pth_wrk)//trim(pth_grd)//
     $           trim(pfx_dmn)//"roms_grd.nc"
      write(*,*) "\\",trim(filename)
      call check( nf90_open(filename, NF90_NOWRITE, ncid) )
C--- 2D ---
      call check( nf90_inq_varid(ncid, "lon_rho", varid) )      ! Read east_u, east_v as well, maybe?
      call check( nf90_get_var(ncid, varid, east_e) )
      write(*, *) "[O] east_e retrieved"
      call check( nf90_inq_varid(ncid, "lat_rho", varid) )
      call check( nf90_get_var(ncid, varid, north_e) )
      write(*, *) "[O] north_e retrieved"
      call check( nf90_inq_varid(ncid, "mask_rho", varid) )
      call check( nf90_get_var(ncid, varid, fsm) )
      write(*, *) "[O] free surface mask retrieved"
      call check( nf90_inq_varid(ncid, "h", varid) )
      call check( nf90_get_var(ncid, varid, h) )
      do i=1,im
        do j=1,jm
            if (fsm(i,j)==0) h(i,j) = 1.
        end do
      end do
      write(*, *) "[O] bottom depth retrieved"
      call check( nf90_close(ncid) )

!      filename = trim(pth_wrk)//trim(pth_flx)//
!     $           trim(pfx_dmn)//"roms_frc.nc"
!      write(*,*) "\\",trim(filename)
!      call check( nf90_open(filename, NF90_NOWRITE, ncid) )
C--- (Constant) Wind stress
!      call check( nf90_inq_varid(ncid, "sustr", varid) )
!      call check( nf90_get_var(ncid, varid, wusurf(2:im,:)) )
!      wusurf(1,:) = wusurf(2,:)*.5
!      wusurf(:,:) = wusurf(:,:)/rhoref
!      write(*, *) "[O] Zonal component of momentum flux retrieved"
!      call check( nf90_inq_varid(ncid, "svstr", varid) )
!      call check( nf90_get_var(ncid, varid, wvsurf(:,2:jm)) )
!      wvsurf(:,1) = wvsurf(:,2)*.5
!      wvsurf(:,:) = wvsurf(:,:)/rhoref
!      write(*, *) "[O] Meridional component of momentum flux retrieved"

C      Simulate from zero elevation to avoid artificial waves during spin-up

!      call check( nf90_close(ncid) )

!    Read from roms_clm and interpolate in time depending on time_offset
      filename = trim(pth_wrk)//trim(pth_bry)//
     $           trim(pfx_dmn)//"roms_clm.nc"
      write(*,*) "\\",trim(filename)
      call check( nf90_open(filename, NF90_NOWRITE, ncid) )
!    Get month for time0
      call upd_mnth

      call check( nf90_inq_varid(ncid, "temp", varid) )

      if (m.ne.1) then
        call check( nf90_get_var(ncid, varid, datr(:,:,kbm1:1:-1,:),
     $                          (/1,1,1,m-1/), (/im,jm,kbm1,2/)) )
      else
        call check( nf90_get_var(ncid, varid, datr(:,:,kbm1:1:-1,1),
     $                          (/1,1,1,12/),(/im,jm,kbm1,1/)) )
        call check( nf90_get_var(ncid, varid, datr(:,:,kbm1:1:-1,2),
     $                          (/1,1,1,1/), (/im,jm,kbm1,1/)) )
      end if

!    FSM is not defined yet! Be careful not to apply it!
      t(:,:,1:kbm1) = datr(:,:,1:kbm1,1)+fac*(datr(:,:,1:kbm1,2)
     $               -datr(:,:,1:kbm1,1))

      call check( nf90_inq_varid(ncid, "salt", varid) )

      if (m.ne.1) then
        call check( nf90_get_var(ncid, varid, datr,
     $                          (/1,1,1,m-1/), (/im,jm,kbm1,2/)) )
      else
        call check( nf90_get_var(ncid, varid, datr(:,:,:,1),
     $                          (/1,1,1,12/),(/im,jm,kbm1,1/)) )
        call check( nf90_get_var(ncid, varid, datr(:,:,:,2),
     $                          (/1,1,1,1/), (/im,jm,kbm1,1/)) )
      end if

      do k=1,kbm1
        s(:,:,k) = datr(:,:,kb-k,1)+fac*(datr(:,:,kb-k,2)
     $            -datr(:,:,kb-k,1))
      end do
!
!   Read elevation
!
      elb = 0.
      el  = 0.
!      call check( nf90_inq_varid(ncid, "zeta", varid) )
!
!      if (m.ne.1) then
!        call check( nf90_get_var(ncid, varid, datr,
!     $                          (/1,1,1,m-1/), (/im,jm,1,2/)) )
!      else
!        call check( nf90_get_var(ncid, varid, datr(:,:,1,1),
!     $                          (/1,1,1,12/),(/im,jm,1,1/)) )
!        call check( nf90_get_var(ncid, varid, datr(:,:,1,2),
!     $                          (/1,1,1,1/), (/im,jm,1,1/)) )
!      end if
!
!      elb(:,:) = datr(:,:,1,1)+fac*(datr(:,:,1,2)
!     $            -datr(:,:,1,1))
!      el = elb      ! rwnd: Is this correct?
!
!    Calculate annual means
!
      call check( nf90_inq_varid(ncid, "temp", varid) )
      do mm=1,12
        call check( nf90_get_var(ncid, varid, datr(:,:,kbm1:1:-1,1),
     $                          (/1,1,1,m/),(/im,jm,kbm1,1/)) )
        do i=1,im
          do j=1,jm
            do k=1,kb
              tclim(i,j,k) = tclim(i,j,k)+datr(i,j,k,1)*lom(mm)
            end do
          end do
        end do
      end do
      tclim = tclim/365.25
      call check( nf90_inq_varid(ncid, "salt", varid) )
      do mm=1,12
        call check( nf90_get_var(ncid, varid, datr(:,:,kbm1:1:-1,1),
     $                          (/1,1,1,m/),(/im,jm,kbm1,1/)) )
        do i=1,im
          do j=1,jm
            do k=1,kb
              sclim(i,j,k) = sclim(i,j,k)+datr(i,j,k,1)*lom(mm)
            end do
          end do
        end do
      end do
      sclim = sclim/365.25
!      call check( nf90_close(ncid) )
!
!    Read annual means
!
!      filename = trim(pth_wrk)//trim(pth_grd)//
!     $           trim(pfx_dmn)//"pom_ann.nc"
!      write(*,*) "\\",trim(filename)
!      call check( nf90_open(filename, NF90_NOWRITE, ncid) )

!      call check( nf90_inq_varid(ncid, "t_ann", varid) )
!      call check( nf90_get_var(ncid, varid, tclim) )
!      call check( nf90_inq_varid(ncid, "s_ann", varid) )
!      call check( nf90_get_var(ncid, varid, sclim) )
!      call check( nf90_inq_varid(ncid, "r_ann", varid) )
!      call check( nf90_get_var(ncid, varid, rmean) )
!    Test rmean calculations
!      call dens(sclim,tclim,rho)
!      rho = rmean-rho
!      cff = 0.
!      dlon = 0.
!      do i=1,im
!        do j=1,jm
!          do k=1,kb
!            cff = cff+rho(i,j,k)
!            if (dlon<abs(rho(i,j,k))) dlon=rho(i,j,k)
!          end do
!        end do
!      end do
!      write(*,*) "Max difference in densities:  ", dlon
!      cff = cff/im/jm/kb
!      write(*,*) "Mean difference in densities: ", cff
      call check( nf90_close(ncid) )
      tclim = t
      sclim = s
      call dens(sclim,tclim,rmean)
!
      write(*, *) "[+] Finished reading IC."
C
C --- print vertical grid distribution
C
      write(6,2)
    2 format(/2x,'k',7x,'z',9x,'zz',9x,'dz',9x,'dzz',/)
      write(6,'(''  '',/)')
      do k=1,kb
        write(6,3) k,z(k),zz(k),dz(k),dzz(k)
    3   format((' ',i5,4f10.3))
      end do
      write(6,'(''  '',//)')
C
C --- calc. Curiolis Parameter
C
        do j=1,jm
          do i=1,im
            cor(i,j)=2.*7.29E-5*sin(north_e(i,j)*rad)
            aam2d(i,j)=aam(i,j,1)   ! RWND: aam is initialized with aam_init already
            !elb(i,j)=0.            ! RWND (Read few lines above or initialized to zero)
            etb(i,j)=elb(i,j)       ! RWND //PV: 0.
            dt(i,j)=h(i,j)+etb(i,j) ! RWND //PV: h(i,j)
          end do
        end do
C
        do j=1,jm
          do i=2,im-1
            dx(i,j)=0.5*rad*re*sqrt(((east_e(i+1,j)-east_e(i-1,j))
     1 *cos(north_e(i,j)*rad))**2+(north_e(i+1,j)-north_e(i-1,j))**2)
          end do
            dx(1,j)=dx(2,j)
            dx(im,j)=dx(im-1,j)
        end do
C
        do i=1,im
          do j=2,jm-1
            dy(i,j)=0.5*rad*re*sqrt(((east_e(i,j+1)-east_e(i,j-1))
     1 *cos(north_e(i,j)*rad))**2+(north_e(i,j+1)-north_e(i,j-1))**2)
          end do
            dy(i,1)=dy(i,2)
            dy(i,jm)=dy(i,jm-1)
        end do
C
C     Calculate areas and masks:
C
      call areas_masks
!
!      Comment the code below to avoid forced closed boundary.
!      do j=1,jm
!        fsm( 1, j) = 0.0
!        fsm(im, j) = 0.0
!      end do
!      do i=1,im
!        fsm( i, 1) = 0.0
!        fsm( i,jm) = 0.0
!      end do
!!      Recalculate masks for u and v
!!      TODO: this can be done simplier by touching only near-boundary cells
!      do j=2,jm
!        do i=2,im
!          dum(i,j)=fsm(i,j)*fsm(i-1,j)
!          dvm(i,j)=fsm(i,j)*fsm(i,j-1)
!        end do
!      end do
!!!!!        call slpmax
C
C --- calc. surface & lateral BC from climatology
C
!        do j=1,jm
!          do i=1,im
!             tsurf(i,j)=t(i,j,1)
!             ssurf(i,j)=s(i,j,1)
!            do k=1,kb
!              tclim(i,j,k)=t(i,j,k)
!              sclim(i,j,k)=s(i,j,k)
!            end do
!          end do
!        end do
C
C                    --- EAST & WEST BCs ---
        do j=1,jm
            ele(j)=0.
            elw(j)=0.
C --- other vel. BCs (fixed in time) can be specified here
            do k=1,kb
              ubw(j,k)=0.                ! RWND
              ube(j,k)=0.                !
              uabe(j) =0.                ! RWND
              uabw(j) =0.                !  Calculate vertical averages
              tbw(j,k)=tclim(1,j,k)
              sbw(j,k)=sclim(1,j,k)
              tbe(j,k)=tclim(im,j,k)
              sbe(j,k)=sclim(im,j,k)
            end do
        end do
C                    --- NORTH & SOUTH BCs ---
        do i=1,im
              els(i)=0. !elb(i,1)   ! RWND
              eln(i)=0. !elb(i,jm)  !
            do k=1,kb
              vbs(i,k)=0.                ! RWND
              vbn(i,k)=0.                !
              vabs(j) =0.                ! RWND
              vabn(j) =0.                !  Calculate vertical means
              tbs(i,k)=tclim(i,1,k)
              sbs(i,k)=sclim(i,1,k)
              tbn(i,k)=tclim(i,jm,k)
              sbn(i,k)=sclim(i,jm,k)
            end do
        end do
C
C     Set initial conditions:
C       and apply free-surface mask ! rwnd:
      do k=1,kb
        t(:,:,k) = t(:,:,k)*fsm(:,:)
        s(:,:,k) = s(:,:,k)*fsm(:,:)
        do j=1,jm
          do i=1,im
            tb(i,j,k)=t(i,j,k)
            sb(i,j,k)=s(i,j,k)
            ub(i,j,k)=u(i,j,k)
            vb(i,j,k)=v(i,j,k)
          end do
        end do
      end do
C
      call dens(sb,tb,rho)
!      rmean = rho   ! remove the line to avoid rmean overriding
C
C
C --- the following grids are needed only for netcdf plotting
C
C     Corner of cell points:
C
      do j=2,jm
        do i=2,im
          east_c(i,j)=(east_e(i,j)+east_e(i-1,j)
     $                  +east_e(i,j-1)+east_e(i-1,j-1))/4.e0
          north_c(i,j)=(north_e(i,j)+north_e(i-1,j)
     $                   +north_e(i,j-1)+north_e(i-1,j-1))/4.e0
        end do
      end do
C
C
C     Extrapolate ends (approx.):
C
      do i=2,im
        east_c(i,1)=2.*east_c(i,2)-east_c(i,3)
        north_c(i,1)=2.*north_c(i,2)-north_c(i,3)
      end do
        east_c(1,1)=2.*east_c(2,1)-east_c(3,1)
C
      do j=2,jm
        east_c(1,j)=2.*east_c(2,j)-east_c(3,j)
        north_c(1,j)=2.*north_c(2,j)-north_c(3,j)
      end do
        north_c(1,1)=2.*north_c(1,2)-north_c(1,3)
C
C     u-points:
C
      do j=1,jm-1
        do i=1,im
          east_u(i,j)=(east_c(i,j)+east_c(i,j+1))/2.e0
          north_u(i,j)=(north_c(i,j)+north_c(i,j+1))/2.e0
        end do
      end do
C
C     Extrapolate ends:
C
      do i=1,im
        east_u(i,jm)=(east_c(i,jm)*3.e0-east_c(i,jm-1))/2.e0
        north_u(i,jm)=(north_c(i,jm)*3.e0-north_c(i,jm-1))/2.e0
      end do
C
C     v-points:
C
      do j=1,jm
        do i=1,im-1
          east_v(i,j)=(east_c(i,j)+east_c(i+1,j))/2.e0
          north_v(i,j)=(north_c(i,j)+north_c(i+1,j))/2.e0
        end do
      end do
C
C     Extrapolate ends:
C
      do j=1,jm
        east_v(im,j)=(east_c(im,j)*3.e0-east_c(im-1,j))/2.e0
        north_v(im,j)=(north_c(im,j)*3.e0-north_c(im-1,j))/2.e0
      end do
C
C     rot is the angle (radians, anticlockwise) of the i-axis relative
C     to east, averaged to a cell centre: (only needed for CDF plotting)
C
      do j=1,jm
        do i=1,im-1
          rot(i,j)=0.
          dlat=north_e(i+1,j)-north_e(i,j)
          dlon= east_e(i+1,j)- east_e(i,j)
           if(dlon.ne.0.) rot(i,j)=atan(dlat/dlon)
        end do
       rot(im,j)=rot(im-1,j)
      end do
C
C     Set lateral boundary conditions, for use in subroutine bcond
C     set all=0 for closed BCs.
C     Values=0 for vel BC only, =1 is combination of vel+elev.
      rfe=0.e0
      rfw=0.e0
      rfn=0.e0
      rfs=0.e0  ! Meaningless with RaS boundary conditions.
C
      return

      contains
      subroutine check(status)
        integer, intent ( in) :: status
        if(status /= nf90_noerr) then
          stop "Stopped"
        end if
      end subroutine check
      end
C
      subroutine ncdf2ic_orig
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Sets up my own problem.                             *
C *                                                                    *
C * This example reads IC from NetCDF files, generated by <deleted>    *
C * Only minimal number of fields are read,                            *
C * while others are calculated here.  TODO: Read as much as we can    *
C *                                                                    *
C **********************************************************************
C
      use netcdf
      implicit none
C
      include 'pom2k.c'
C
      real dat(imm1, kbm1)
      real datr(im, jm, kbm1, 2)
      real datu(imm1, jm, kbm1, 2)
      real datv(im, jmm1, kbm1, 2)
C
      real rad,re,dlat,dlon,cff
      integer i,j,k
      character*5 field
      rad=0.01745329
      re=6371.E3
C
      write(6,'(/,'' Read grid and initial conditions '',/)')
C
C--- 1D ---
      filename = trim(pth_wrk)//trim(pth_grd)//
     $           trim(pfx_dmn)//"pom_grd.nc"
      write(*,*) filename
      call check( nf90_open(filename, NF90_NOWRITE, ncid) )
      write(*, *) ncid
      call check( nf90_inq_varid(ncid, "z", varid) )
      call check( nf90_get_var(ncid, varid, z) )
      write(*, *) "[O] z retrieved"
      call check( nf90_inq_varid(ncid, "zz", varid) )
      call check( nf90_get_var(ncid, varid, zz) )
      write(*, *) "[O] zz retrieved"
      call check( nf90_inq_varid(ncid, "dz", varid) )
      call check( nf90_get_var(ncid, varid, dz) )
      write(*, *) "[O] dz retrieved"
      call check( nf90_inq_varid(ncid, "dzz", varid) )
      call check( nf90_get_var(ncid, varid, dzz) )
      write(*, *) "[O] dzz retrieved"
C--- 2D ---
      call check( nf90_inq_varid(ncid, "lon_rho", varid) )
      call check( nf90_get_var(ncid, varid, east_e) )
      write(*, *) "[O] east_e retrieved"
      call check( nf90_inq_varid(ncid, "lat_rho", varid) )
      call check( nf90_get_var(ncid, varid, north_e) )
      write(*, *) "[O] north_e retrieved"
      call check( nf90_inq_varid(ncid, "h", varid) )
      call check( nf90_get_var(ncid, varid, h) )
      write(*, *) "[O] bottom depth retrieved"
C--- 3D ---
      call check( nf90_inq_varid(ncid, "T", varid) )
      call check( nf90_get_var(ncid, varid, t) )
      write(*, *) "[O] potential temperature retrieved"
      call check( nf90_inq_varid(ncid, "S", varid) )
      call check( nf90_get_var(ncid, varid, s) )
      write(*, *) "[O] salinity retrieved"
      call check( nf90_inq_varid(ncid, "rmean", varid) )
      call check( nf90_get_var(ncid, varid, rmean) )
      write(*, *) "[O] mean density retrieved"
C--- (Constant) Wind stress
      call check( nf90_inq_varid(ncid, "WUsur", varid) )
      call check( nf90_get_var(ncid, varid, wusurf) )
      write(*, *) "[O] Zonal component of momentum flux retrieved"
      call check( nf90_inq_varid(ncid, "WVsur", varid) )
      call check( nf90_get_var(ncid, varid, wvsurf) )
      write(*, *) "[O] Meridional component of momentum flux retrieved"
      call check( nf90_inq_varid(ncid, "swrad", varid) )
      call check( nf90_get_var(ncid, varid, swrad) )
      write(*, *) "[O] Short wave radiation retrieved"
!      call check( nf90_inq_varid(ncid, "wtsur", varid) )
!      call check( nf90_get_var(ncid, varid, wtsurf) )
!      write(*, *) "[O] Heat flux retrieved"
!      call check( nf90_inq_varid(ncid, "U", varid) )
!      call check( nf90_get_var(ncid, varid, u) )
!      write(*, *) "[O] Longitudal component of current velocity ",
!     $"retrieved"
!      do i=1,im
!        do j=1,jm
!          do k=1,kb
!            uab(i,j) = uab(i,j)+u(i,j,k)*dz(k)
!          end do
!        end do
!      end do
!      write(*, *) "[O] Integrals of longitudal component of current",
!     $"velocity computed"
!      call check( nf90_inq_varid(ncid, "V", varid) )
!      call check( nf90_get_var(ncid, varid, v) )
!      write(*, *) "[O] Latitudal component of current velocity ",
!     $"retrieved"
!      do i=1,im
!        do j=1,jm
!          do k=1,kb
!            vab(i,j) = vab(i,j)+v(i,j,k)*dz(k)
!          end do
!        end do
!      end do
!      write(*, *) "[O] Integrals of latitudal component of current ",
!     $"velocity computed"
!      call check( nf90_inq_varid(ncid, "SSH", varid) )
!      call check( nf90_get_var(ncid, varid, elb) )
!      el = elb  ! RWND
!      write(*, *) "[O] Sea surface height retrieved"

C      Simulate from zero elevation to avoid artificial waves during spin-up

      call check( nf90_close(ncid) )

!    Read from roms_clm and interpolate in time depending on time_offset
      filename = trim(pth_wrk)//trim(pth_bry)//
     $           trim(pfx_dmn)//"roms_clm.nc"
      call check( nf90_open(filename, NF90_NOWRITE, ncid) )
!    Get month for time0
      call upd_mnth

      call check( nf90_inq_varid(ncid, "temp", varid) )

      if (m.ne.1) then
        call check( nf90_get_var(ncid, varid, datr(:,:,kbm1:1:-1,:),
     $                          (/1,1,1,m-1/), (/im,jm,kbm1,2/)) )
      else
        call check( nf90_get_var(ncid, varid, datr(:,:,kbm1:1:-1,1),
     $                          (/1,1,1,12/),(/im,jm,kbm1,1/)) )
        call check( nf90_get_var(ncid, varid, datr(:,:,kbm1:1:-1,2),
     $                          (/1,1,1,1/), (/im,jm,kbm1,1/)) )
      end if

!    FSM is not defined yet! Be careful not to apply it!
        t(:,:,1:kbm1) = datr(:,:,1:kbm1,1)+fac*(datr(:,:,1:kbm1,2)
     $             -datr(:,:,1:kbm1,1))

      call check( nf90_inq_varid(ncid, "salt", varid) )

      if (m.ne.1) then
        call check( nf90_get_var(ncid, varid, datr,
     $                          (/1,1,1,m-1/), (/im,jm,kbm1,2/)) )
      else
        call check( nf90_get_var(ncid, varid, datr(:,:,:,1),
     $                          (/1,1,1,12/),(/im,jm,kbm1,1/)) )
        call check( nf90_get_var(ncid, varid, datr(:,:,:,2),
     $                          (/1,1,1,1/), (/im,jm,kbm1,1/)) )
      end if

      do k=1,kbm1
        s(:,:,k) = datr(:,:,kb-k,1)+fac*(datr(:,:,kb-k,2)
     $            -datr(:,:,kb-k,1))
      end do

      call check( nf90_close(ncid) )
!
!    Read annual means
!
      filename = trim(pth_wrk)//trim(pth_bry)//
     $           trim(pfx_dmn)//"pom_ann.nc"
      call check( nf90_open(filename, NF90_NOWRITE, ncid) )
      call check( nf90_inq_varid(ncid, "t_ann", varid) )
      call check( nf90_get_var(ncid, varid, tclim) )
      call check( nf90_inq_varid(ncid, "s_ann", varid) )
      call check( nf90_get_var(ncid, varid, sclim) )
      call check( nf90_inq_varid(ncid, "r_ann", varid) )
      call check( nf90_get_var(ncid, varid, rmean) )
      call check( nf90_close(ncid) )
!
      write(*, *) "[+] Finished reading IC."
C
C --- print vertical grid distribution
C
      write(6,2)
    2 format(/2x,'k',7x,'z',9x,'zz',9x,'dz',9x,'dzz',/)
      write(6,'(''  '',/)')
      do k=1,kb
        write(6,3) k,z(k),zz(k),dz(k),dzz(k)
    3   format((' ',i5,4f10.3))
      end do
      write(6,'(''  '',//)')
C
C --- calc. Curiolis Parameter
C
        do j=1,jm
          do i=1,im
            cor(i,j)=2.*7.29E-5*sin(north_e(i,j)*rad)
            aam2d(i,j)=aam(i,j,1)   ! RWND: aam is initialized with aam_init already
            !elb(i,j)=0.            ! RWND (Read few lines above or initialized to zero)
            etb(i,j)=elb(i,j)       ! RWND //PV: 0.
            dt(i,j)=h(i,j)+etb(i,j) ! RWND //PV: h(i,j)
          end do
        end do
C
        do j=1,jm
          do i=2,im-1
            dx(i,j)=0.5*rad*re*sqrt(((east_e(i+1,j)-east_e(i-1,j))
     1 *cos(north_e(i,j)*rad))**2+(north_e(i+1,j)-north_e(i-1,j))**2)
          end do
            dx(1,j)=dx(2,j)
            dx(im,j)=dx(im-1,j)
        end do
C
        do i=1,im
          do j=2,jm-1
            dy(i,j)=0.5*rad*re*sqrt(((east_e(i,j+1)-east_e(i,j-1))
     1 *cos(north_e(i,j)*rad))**2+(north_e(i,j+1)-north_e(i,j-1))**2)
          end do
            dy(i,1)=dy(i,2)
            dy(i,jm)=dy(i,jm-1)
        end do
C
C     Calculate areas and masks:
C
      call areas_masks
!
!      Comment the code below to avoid forced closed boundary.
!      do j=1,jm
!        fsm( 1, j) = 0.0
!        fsm(im, j) = 0.0
!      end do
!      do i=1,im
!        fsm( i, 1) = 0.0
!        fsm( i,jm) = 0.0
!      end do
!!      Recalculate masks for u and v
!!      TODO: this can be done simplier by touching only near-boundary cells
!      do j=2,jm
!        do i=2,im
!          dum(i,j)=fsm(i,j)*fsm(i-1,j)
!          dvm(i,j)=fsm(i,j)*fsm(i,j-1)
!        end do
!      end do
!!!!!        call slpmax
C
C --- calc. surface & lateral BC from climatology
C
        do j=1,jm
          do i=1,im
             tsurf(i,j)=t(i,j,1)
             ssurf(i,j)=s(i,j,1)
!            do k=1,kb
!              tclim(i,j,k)=t(i,j,k)
!              sclim(i,j,k)=s(i,j,k)
!            end do
          end do
        end do
C
C                    --- EAST & WEST BCs ---
        do j=1,jm
            ele(j)=0.
            elw(j)=0.
C --- other vel. BCs (fixed in time) can be specified here
            do k=1,kb
              ubw(j,k)=0                ! RWND
              ube(j,k)=0                !
              uabe(j) =0                ! RWND
              uabw(j) =0                !  Calculate vertical averages
              tbw(j,k)=tclim(1,j,k)
              sbw(j,k)=sclim(1,j,k)
              tbe(j,k)=tclim(im,j,k)
              sbe(j,k)=sclim(im,j,k)
            end do
        end do
C                    --- NORTH & SOUTH BCs ---
        do i=1,im
              els(i)=0. !elb(i,1)   ! RWND
              eln(i)=0. !elb(i,jm)  !
            do k=1,kb
              vbs(i,k)=0                ! RWND
              vbn(i,k)=0                !
              vabs(j) =0                ! RWND
              vabn(j) =0                !  Calculate vertical means
              tbs(i,k)=tclim(i,1,k)
              sbs(i,k)=sclim(i,1,k)
              tbn(i,k)=tclim(i,jm,k)
              sbn(i,k)=sclim(i,jm,k)
            end do
        end do
C
C     Set initial conditions:
C       and apply free-surface mask ! rwnd:
      do k=1,kb
        t(:,:,k) = t(:,:,k)*fsm(:,:)
        s(:,:,k) = s(:,:,k)*fsm(:,:)
        do j=1,jm
          do i=1,im
            tb(i,j,k)=t(i,j,k)
            sb(i,j,k)=s(i,j,k)
            ub(i,j,k)=u(i,j,k)
            vb(i,j,k)=v(i,j,k)
          end do
        end do
      end do
C
      call dens(sb,tb,rho)
!      rmean = rho   ! remove the line to avoid rmean overriding
C
C
C --- the following grids are needed only for netcdf plotting
C
C     Corner of cell points:
C
      do j=2,jm
        do i=2,im
          east_c(i,j)=(east_e(i,j)+east_e(i-1,j)
     $                  +east_e(i,j-1)+east_e(i-1,j-1))/4.e0
          north_c(i,j)=(north_e(i,j)+north_e(i-1,j)
     $                   +north_e(i,j-1)+north_e(i-1,j-1))/4.e0
        end do
      end do
C
C
C     Extrapolate ends (approx.):
C
      do i=2,im
        east_c(i,1)=2.*east_c(i,2)-east_c(i,3)
        north_c(i,1)=2.*north_c(i,2)-north_c(i,3)
      end do
        east_c(1,1)=2.*east_c(2,1)-east_c(3,1)
C
      do j=2,jm
        east_c(1,j)=2.*east_c(2,j)-east_c(3,j)
        north_c(1,j)=2.*north_c(2,j)-north_c(3,j)
      end do
        north_c(1,1)=2.*north_c(1,2)-north_c(1,3)
C
C     u-points:
C
      do j=1,jm-1
        do i=1,im
          east_u(i,j)=(east_c(i,j)+east_c(i,j+1))/2.e0
          north_u(i,j)=(north_c(i,j)+north_c(i,j+1))/2.e0
        end do
      end do
C
C     Extrapolate ends:
C
      do i=1,im
        east_u(i,jm)=(east_c(i,jm)*3.e0-east_c(i,jm-1))/2.e0
        north_u(i,jm)=(north_c(i,jm)*3.e0-north_c(i,jm-1))/2.e0
      end do
C
C     v-points:
C
      do j=1,jm
        do i=1,im-1
          east_v(i,j)=(east_c(i,j)+east_c(i+1,j))/2.e0
          north_v(i,j)=(north_c(i,j)+north_c(i+1,j))/2.e0
        end do
      end do
C
C     Extrapolate ends:
C
      do j=1,jm
        east_v(im,j)=(east_c(im,j)*3.e0-east_c(im-1,j))/2.e0
        north_v(im,j)=(north_c(im,j)*3.e0-north_c(im-1,j))/2.e0
      end do
C
C     rot is the angle (radians, anticlockwise) of the i-axis relative
C     to east, averaged to a cell centre: (only needed for CDF plotting)
C
      do j=1,jm
        do i=1,im-1
          rot(i,j)=0.
          dlat=north_e(i+1,j)-north_e(i,j)
          dlon= east_e(i+1,j)- east_e(i,j)
           if(dlon.ne.0.) rot(i,j)=atan(dlat/dlon)
        end do
       rot(im,j)=rot(im-1,j)
      end do
C
C     Set lateral boundary conditions, for use in subroutine bcond
C     set all=0 for closed BCs.
C     Values=0 for vel BC only, =1 is combination of vel+elev.
      rfe=0.e0
      rfw=0.e0
      rfn=1.e0
      rfs=1.e0  ! Meaningless with RaS boundary conditions.
C
      return

      contains
      subroutine check(status)
        integer, intent ( in) :: status
        if(status /= nf90_noerr) then
          stop "Stopped"
        end if
      end subroutine check
      end
C
      subroutine ic2ncdf
C **********************************************************************
C *                                                                    *
C *                         POM2K SOURCE CODE                          *
C *                                                                    *
C * ROUTINE NAME:  ic2ncdf                                             *
C *                                                                    *
C * FUNCTION    :  Creates a netCDF file for further output            *
C *                                                                    *
C **********************************************************************
C
        use netcdf
        implicit none
C
        include 'pom2k.c'
C
        integer i, j, k, dim_strim, dim_auxuv
        real bot_depth(im,jm,kb)
C
        filename = trim(pth_wrk)//trim(pth_out)//
     $             trim(title)//"_data.nc"
        call check( nf90_create(filename, NF90_CLASSIC_MODEL, ncid) )
        call check( nf90_def_dim(ncid, "time", NF90_UNLIMITED
     $          , dim_time) )
        call check( nf90_def_dim(ncid, "s_rho", kb, dim_srho) )
        call check( nf90_def_dim(ncid, "s_w", kb, dim_sw) )
        call check( nf90_def_dim(ncid, "s_trim", 3, dim_strim) )
        call check( nf90_def_dim(ncid, "aux_uv", 2, dim_auxuv) )
        call check( nf90_def_dim(ncid, "latitude", jm, dim_lat) )
        call check( nf90_def_dim(ncid, "longitude", im, dim_lon) )
        call check( nf90_def_var(ncid, "Time", NF90_DOUBLE,
     $          (/ dim_time /), varid) )
        call check( nf90_def_var(ncid, "Level", NF90_DOUBLE,
     $          (/ dim_sw /), varid) )
        if (mode.ge.3) then
          call check( nf90_def_var(ncid, "Depth", NF90_DOUBLE,
     $          (/ dim_lon, dim_lat, dim_sw /), varid) )
        end if
!        call check( nf90_def_var(ncid, "FSM", NF90_DOUBLE,
!     $          (/ dim_lon, dim_lat /), varid) )
!        call check( nf90_def_var(ncid, "CFL", NF90_DOUBLE,
!     $          (/ dim_lon, dim_lat /), varid) )
        call check( nf90_def_var(ncid, "Tavg", NF90_DOUBLE,
     $          (/ dim_time /), varid) )
!        call check( nf90_def_var(ncid, "Vtot", NF90_DOUBLE,
!     $          (/ dim_time /), varid) )
        call check( nf90_def_var(ncid, "Eavg", NF90_DOUBLE,
     $          (/ dim_time /), varid) )
        call check( nf90_def_var(ncid, "Qavg", NF90_DOUBLE,
     $          (/ dim_time /), varid) )
        call check( nf90_def_var(ncid, "Mtot", NF90_DOUBLE,
     $          (/ dim_time /), varid) )
!        call check( nf90_def_var(ncid, "tgt-point_1", NF90_DOUBLE,
!     $          (/ dim_auxuv, dim_time /), varid) ) ! u,v
        call check( nf90_def_var(ncid, "Latitude", NF90_DOUBLE,
     $          (/ dim_lon, dim_lat /), varid) )
        call check( nf90_def_var(ncid, "Longitude", NF90_DOUBLE,
     $          (/ dim_lon, dim_lat /), varid) )
        if (mode.ge.3) then
!          call check( nf90_def_var(ncid, "U", NF90_DOUBLE,
!     $          (/ dim_lon, dim_lat, dim_srho, dim_time /), varid) )
!          call check( nf90_def_var(ncid, "V", NF90_DOUBLE,
!     $          (/ dim_lon, dim_lat, dim_srho, dim_time /), varid) )

!          call check( nf90_def_var(ncid, "W", NF90_DOUBLE,
!     $          (/ dim_lon, dim_lat, dim_srho, dim_time /), varid) )
! TODO: Is it normal to leave it with dim_srho and not change to dim_sw?
!          call check( nf90_def_var(ncid, "T", NF90_DOUBLE,
!     $          (/ dim_lon, dim_lat, dim_srho, dim_time /), varid) )
!          call check( nf90_def_var(ncid, "S", NF90_DOUBLE,
!     $          (/ dim_lon, dim_lat, dim_srho, dim_time /), varid) )
!          call check( nf90_def_var(ncid, "RHO", NF90_DOUBLE,
!     $          (/ dim_lon, dim_lat, dim_srho, dim_time /), varid) )
!          call check( nf90_def_var(ncid, "AAM",
!     $ NF90_DOUBLE, (/ dim_lon, dim_lat, dim_srho, dim_time /), varid) )
!          call check( nf90_def_var(ncid, "KM",
!     $ NF90_DOUBLE, (/ dim_lon, dim_lat, dim_srho, dim_time /), varid) )
!          call check( nf90_def_var(ncid, "Q2",
!     $ NF90_DOUBLE, (/ dim_lon, dim_lat, dim_srho, dim_time /), varid) )
        end if
!
        if (mode.ge.3) then
        call check( nf90_def_var(ncid, "SSU", NF90_DOUBLE,
     $          (/ dim_lon, dim_lat, dim_strim, dim_time /), varid) )
        call check( nf90_def_var(ncid, "SSV", NF90_DOUBLE,
     $          (/ dim_lon, dim_lat, dim_strim, dim_time /), varid) )
        call check( nf90_def_var(ncid, "SST", NF90_DOUBLE,
     $          (/ dim_lon, dim_lat, dim_strim, dim_time /), varid) )
        call check( nf90_def_var(ncid, "SSS", NF90_DOUBLE,
     $          (/ dim_lon, dim_lat, dim_strim, dim_time /), varid) )
        call check( nf90_def_var(ncid, "SSR", NF90_DOUBLE,
     $          (/ dim_lon, dim_lat, dim_strim, dim_time /), varid) )
        end if
        call check( nf90_def_var(ncid, "EL",  NF90_DOUBLE,
     $          (/ dim_lon, dim_lat, dim_time /), varid) )
        call check( nf90_def_var(ncid, "UA",  NF90_DOUBLE,
     $          (/ dim_lon, dim_lat, dim_time /), varid) )
        call check( nf90_def_var(ncid, "VA",  NF90_DOUBLE,
     $          (/ dim_lon, dim_lat, dim_time /), varid) )
!        call check( nf90_def_var(ncid, "T_S",
!     $ NF90_DOUBLE, (/ dim_lon, dim_sw, dim_time /), varid) )
!        call check( nf90_def_var(ncid, "S_S",
!     $ NF90_DOUBLE, (/ dim_lon, dim_sw, dim_time /), varid) )
!        call check( nf90_def_var(ncid, "U_S",
!     $ NF90_DOUBLE, (/ dim_lon, dim_sw, dim_time /), varid) )
!        call check( nf90_def_var(ncid, "V_S",
!     $ NF90_DOUBLE, (/ dim_lon, dim_sw, dim_time /), varid) )
!        call check( nf90_def_var(ncid, "ELS",
!     $ NF90_DOUBLE, (/ dim_lon, dim_time /), varid) )
        call check( nf90_def_var(ncid, "PSI_nw",
     $ NF90_DOUBLE, (/ dim_lon, dim_lat, dim_time /), varid) )
        call check( nf90_def_var(ncid, "PSI_ew",
     $ NF90_DOUBLE, (/ dim_lon, dim_lat, dim_time /), varid) )
        call check( nf90_enddef(ncid) )
        if (mode.ge.3) then
        ! Calculate Depth array
          do k=1,kb
            do j=1,jm
              do i=1,im
                bot_depth(i,j,k) = z(k)*h(i,j)
              end do
            end do
          end do
          call check( nf90_inq_varid(ncid, "Depth", varid) )
          call check( nf90_put_var(ncid, varid, bot_depth) )
        end if
        call check( nf90_inq_varid(ncid, "Latitude", varid) )
        call check( nf90_put_var(ncid, varid, north_e) )
        call check( nf90_inq_varid(ncid, "Longitude", varid) )
        call check( nf90_put_var(ncid, varid, east_e) )
        call check( nf90_inq_varid(ncid, "Level", varid) )
        call check( nf90_put_var(ncid, varid, z) )
!        call check( nf90_inq_varid(ncid, "FSM", varid) )
!        call check( nf90_put_var(ncid, varid, fsm) )
!        call check( nf90_inq_varid(ncid, "CFL", varid) )
!        call check( nf90_put_var(ncid, varid, tps) )
        write(*, *) "NetCDF output has been initialized (",ncid,")."
        call check( nf90_close(ncid) )
C
        return
C
        contains
          subroutine check(status)
            integer, intent ( in) :: status
            if(status /= nf90_noerr) then
              write(*,*) "Error status: ", status
              stop "Stopped"
            end if
          end subroutine check
C
      end
C
      subroutine ncTgtFlush(idx)
!
        implicit none
!
        include 'pom2k.c'
!
        integer :: idx
!
        write(*,*) "...print..."
        write(idx, *) time, ";", u(70,155,1), ";", v(70,155,1)
!
      end
!
      subroutine ncflush(ptime)
C **********************************************************************
C *                                                                    *
C *                         POM2K SOURCE CODE                          *
C *                                                                    *
C * ROUTINE NAME:  ncflush                                             *
C *                                                                    *
C * FUNCTION    :  Writes a set of outputs to netCDF file              *
C *                                                                    *
C **********************************************************************
C
        use netcdf
        implicit none
C
        include 'pom2k.c'
C
        integer :: ptime
        logical :: NOK
        integer :: count
        integer :: i, j, k
        real :: vtot, tavg, atot, eavg, qavg, qtot, mtot
        real :: darea, dvol
        integer nlyrs
        parameter (nlyrs = 3)
        integer :: lyrs(nlyrs)
        data lyrs /1,2,15/
C
        count = 0
        NOK = .true.
        filename = trim(pth_wrk)//trim(pth_out)//
     $             trim(title)//'_data.nc'
        call check( nf90_open(filename, NF90_WRITE, ncid) )
        call check( nf90_inq_varid(ncid, "Time", varid) )
        do while (NOK)
          count = count+1
          status = nf90_put_var(ncid, varid, time, (/ptime/))
          if ((status.eq.nf90_noerr).or.(count.gt.3)) NOK = .false.
        end do
C
          vtot=0.e0
          atot=0.e0
          qtot=0
          mtot=0
          tavg=0.e0
          qavg=0
C
          do j=1,jm
            do i=1,im
              darea=dx(i,j)*dy(i,j)*fsm(i,j)
              do k=1,kbm1
                dvol=darea*dt(i,j)*dz(k)
                mtot=mtot+(rho(i,j,k)*rhoref+1000)*dvol
                vtot=vtot+dvol
                tavg=tavg+tb(i,j,k)*dvol
                qavg=qavg+q2(i,j,k)*dvol
              end do
              atot=atot+darea
              eavg=eavg+et(i,j)*darea
            end do
          end do
C
          eavg=eavg/atot
          tavg=tavg/vtot
          qtot=qavg
          qavg=qtot/vtot
C
        call check( nf90_inq_varid(ncid, "Tavg", varid) )
        call check( nf90_put_var(ncid, varid, tavg, (/ptime/)) )
!        call check( nf90_inq_varid(ncid, "Vtot", varid) )
!        call check( nf90_put_var(ncid, varid, vtot, (/ptime/)) )
        call check( nf90_inq_varid(ncid, "Eavg", varid) )
        call check( nf90_put_var(ncid, varid, eavg, (/ptime/)) )
!        call check( nf90_inq_varid(ncid, "Etot", varid) )
!        call check( nf90_put_var(ncid, varid, eavg+qavg, (/ptime/)) )
        call check( nf90_inq_varid(ncid, "Qavg", varid) )
        call check( nf90_put_var(ncid, varid, qavg, (/ptime/)) )
        call check( nf90_inq_varid(ncid, "Mtot", varid) )
        call check( nf90_put_var(ncid, varid, mtot, (/ptime/)) )
!        call check( nf90_inq_varid(ncid, "tgt-point_1", varid) )
!        call check( nf90_put_var(ncid, varid,u(70,155,1),(/1,ptime/)))
!        call check( nf90_put_var(ncid, varid,v(70,155,1),(/2,ptime/)))
        if (mode.ge.3) then
!          call check( nf90_inq_varid(ncid, "U", varid) )
!          call check( nf90_put_var(ncid, varid, u, (/1,1,1,ptime/)
!     $     ,(/im,jm,kb,1/)) )
!          call check( nf90_inq_varid(ncid, "V", varid) )
!          call check( nf90_put_var(ncid, varid, v, (/1,1,1,ptime/)
!     $     ,(/im,jm,kb,1/)) )
!          call check( nf90_inq_varid(ncid, "W", varid) )
!          call check( nf90_put_var(ncid, varid, w, (/1,1,1,ptime/)
!     $     ,(/im,jm,kb,1/)) )
!          call check( nf90_inq_varid(ncid, "T", varid) )
!          call check( nf90_put_var(ncid, varid, t, (/1,1,1,ptime/)
!     $     ,(/im,jm,kb,1/)) )
!          call check( nf90_inq_varid(ncid, "S", varid) )
!          call check( nf90_put_var(ncid, varid, s, (/1,1,1,ptime/)
!     $     ,(/im,jm,kb,1/)) )
!          call check( nf90_inq_varid(ncid, "RHO", varid) )
!          call check( nf90_put_var(ncid, varid, rho, (/1,1,1,ptime/)
!     $     ,(/im,jm,kb,1/)) )
!          call check( nf90_inq_varid(ncid, "AAM", varid) )
!          call check( nf90_put_var(ncid, varid, aam, (/1,1,1,ptime/)
!     $     ,(/im,jm,kb,1/)) )
!          call check( nf90_inq_varid(ncid, "KM", varid) )
!          call check( nf90_put_var(ncid, varid, km, (/1,1,1,ptime/)
!     $     ,(/im,jm,kb,1/)) )
!          call check( nf90_inq_varid(ncid, "Q2", varid) )
!          call check( nf90_put_var(ncid, varid, q2, (/1,1,1,ptime/)
!     $     ,(/im,jm,kb,1/)) )
        end if
!
        if (mode.ge.3) then
        call check( nf90_inq_varid(ncid, "SSU", varid) )
        call check( nf90_put_var(ncid, varid,   u(:,:,lyrs)
     $   ,(/1,1,1,ptime/),(/im,jm,nlyrs,1/)) )
        call check( nf90_inq_varid(ncid, "SSV", varid) )
        call check( nf90_put_var(ncid, varid,   v(:,:,lyrs)
     $   ,(/1,1,1,ptime/),(/im,jm,nlyrs,1/)) )
        call check( nf90_inq_varid(ncid, "SST", varid) )
        call check( nf90_put_var(ncid, varid,   t(:,:,lyrs)
     $   ,(/1,1,1,ptime/),(/im,jm,nlyrs,1/)) )
        call check( nf90_inq_varid(ncid, "SSS", varid) )
        call check( nf90_put_var(ncid, varid,   s(:,:,lyrs)
     $   ,(/1,1,1,ptime/),(/im,jm,nlyrs,1/)) )
        call check( nf90_inq_varid(ncid, "SSR", varid) )
        call check( nf90_put_var(ncid, varid, rho(:,:,lyrs)
     $   ,(/1,1,1,ptime/),(/im,jm,nlyrs,1/)) )
        end if
        call check( nf90_inq_varid(ncid, "EL", varid) )
        call check( nf90_put_var(ncid, varid, elb, (/1,1,ptime/)
     $   ,(/im,jm,1/)) )
        call check( nf90_inq_varid(ncid, "UA", varid) )
        call check( nf90_put_var(ncid, varid, uab, (/1,1,ptime/)
     $   ,(/im,jm,1/)) )
        call check( nf90_inq_varid(ncid, "VA", varid) )
        call check( nf90_put_var(ncid, varid, vab, (/1,1,ptime/)
     $   ,(/im,jm,1/)) )
!        call check( nf90_inq_varid(ncid, "T_S", varid) )
!        call check( nf90_put_var(ncid, varid, tbs, (/1,1,ptime/)
!     $   ,(/im,kb,1/)) )
!        call check( nf90_inq_varid(ncid, "S_S", varid) )
!        call check( nf90_put_var(ncid, varid, sbs, (/1,1,ptime/)
!     $   ,(/im,kb,1/)) )
!        call check( nf90_inq_varid(ncid, "U_S", varid) )
!        call check( nf90_put_var(ncid, varid, ubs, (/1,1,ptime/)
!     $   ,(/im,kb,1/)) )
!        call check( nf90_inq_varid(ncid, "V_S", varid) )
!        call check( nf90_put_var(ncid, varid, vbs, (/1,1,ptime/)
!     $   ,(/im,kb,1/)) )
!        call check( nf90_inq_varid(ncid, "ELS", varid) )
!        call check( nf90_put_var(ncid, varid, els, (/1,ptime/)
!     $   ,(/im,1/)) )
!        call check( nf90_inq_varid(ncid, "aux1", varid) )
!        call check( nf90_put_var(ncid, varid, sb, (/1,1,1,ptime/)
!     $   ,(/im,jm,kb,1/)) )
!        call check( nf90_inq_varid(ncid, "aux2", varid) )
!        call check( nf90_put_var(ncid, varid, wusurf, (/1,1,ptime/)
!     $   ,(/im,jm,1/)) )
        call check( nf90_close(ncid) )
Cc
        call findpsi2nc(ptime)
C
        return
C
        contains
          subroutine check(status)
            integer, intent ( in) :: status
!            if (DBG) write(*,*) status
            if(status /= nf90_noerr) then
              stop "Stopped"
            end if
          end subroutine check
C
      end
C
      subroutine findpsi2nc(tind)
C **********************************************************************
C *                                                                    *
C * ROUTINE NAME:  findpsi2nc                                          *
C *                                                                    *
C * FUNCTION    :  Calculates the stream function, first assuming      *
C *                zero on the southern boundary and then, using the   *
C *                values on the western boundary, the stream function *
C *                is calculated again. If the elevation field is near *
C *                steady state, the two calculations should agree;    *
C *                otherwise not.                                      *
C *                                                                    *
C **********************************************************************
C
      use netcdf
      implicit none
C
      include 'pom2k.c'
C
      integer i,j, tind
C
      filename = trim(pth_wrk)//trim(pth_out)//trim(title)//"_data.nc"

      call check( nf90_open(filename, NF90_WRITE, ncid) )

      do j=1,jm
        do i=1,im
          psi(i,j)=0.e0
        end do
      end do
C
C     Sweep northward:
C
      do j=2,jmm1
        do i=2,im
          psi(i,j+1)=psi(i,j)
     $                +.25e0*uab(i,j)*(d(i,j)+d(i-1,j))
     $                  *(dy(i,j)+dy(i-1,j))
        end do
      end do
C
      call check( nf90_inq_varid(ncid, "PSI_nw", varid) )
      call check( nf90_put_var(ncid, varid, psi, (/1,1,tind/)) )
C
C    Sweep eastward:
C
      do j=2,jm
        do i=2,imm1
          psi(i+1,j)=psi(i,j)
     $                -.25e0*vab(i,j)*(d(i,j)+d(i,j-1))
     $                  *(dx(i,j)+dx(i,j-1))
        end do
      end do
C
      call check( nf90_inq_varid(ncid, "PSI_ew", varid) )
      call check( nf90_put_var(ncid, varid, psi, (/1,1,tind/)) )
      call check( nf90_close(ncid) )
C
      return
C
      contains
          subroutine check(status)
            integer, intent ( in) :: status
!            if (DBG) write(*,*) status
            if(status /= nf90_noerr) then
              stop "Stopped"
            end if
          end subroutine check
C
      end
C
      subroutine flux(idx)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Reads (if necessary) forcing dataset                *
C *                and interpolates in time.                           *
C *                                                                    *
C **********************************************************************
C
      use netcdf
      implicit none
C
      include 'pom2k.c'
!
      integer, intent (in) :: idx
      integer :: i,j
      real :: qq, dq, sst

      real :: hf_fac = -4.1876e6 ! Heat flux covertion factor
C
      select case (idx)

        case (4) ! Long wave radiation and other fluxes

          if (m.ne.rf_wtsur) then
C
            rf_wtsur = m
C
            filename = trim(pth_wrk)//trim(pth_flx)//
     $                 trim(pfx_dmn)//"roms_frc.nc"
            call check( nf90_open(filename, NF90_NOWRITE, ncid) )
            call check( nf90_inq_varid(ncid, "shflux", varid) )
            if (m.ne.1) then
              call check(nf90_get_var(ncid, varid, qqb,(/1,1,m-1/)))
              call check( nf90_get_var(ncid, varid, qqf, (/1,1,m/)))
            else
              call check( nf90_get_var(ncid, varid, qqb,(/1,1,12/)))
              call check( nf90_get_var(ncid, varid, qqf, (/1,1,1/)))
            end if
            call check( nf90_inq_varid(ncid, "dQdSST", varid) )
            if (m.ne.1) then
              call check(nf90_get_var(ncid, varid, dqb,(/1,1,m-1/)))
              call check( nf90_get_var(ncid, varid, dqf, (/1,1,m/)))
            else
              call check(nf90_get_var(ncid, varid, dqb, (/1,1,12/)))
              call check( nf90_get_var(ncid, varid, dqf, (/1,1,1/)))
            end if
            ! We need to read tsurf again, since we cannot be sure whether it has been read or not.
            call check( nf90_inq_varid(ncid, "SST", varid) )
            if (m.ne.1) then
              call check(nf90_get_var(ncid, varid, tsurfb, (/1,1,m-1/)))
              call check( nf90_get_var(ncid, varid, tsurff, (/1,1,m/)))
            else
              call check( nf90_get_var(ncid, varid, tsurfb, (/1,1,12/)))
              call check( nf90_get_var(ncid, varid, tsurff, (/1,1,1/)))
            end if

            call check( nf90_close(ncid) )
C
            write(*,*) "Read wtsurf:  ", m
C
          end if

          do i=1,im
            do j=1,jm
              ! Interpolate
              qq  = qqb(i,j)   +fac*(qqf(i,j)   -qqb(i,j)   )
              dq  = dqb(i,j)   +fac*(dqf(i,j)   -dqb(i,j)   )
              sst = tsurfb(i,j)+fac*(tsurff(i,j)-tsurfb(i,j))
              ! Convert wtsurf from W/m^2 to K m/s. TODO: convert the qq,dq,sst while reading, maybe?
              wtsurf(i,j) = (qq+dq*(t(i,j,1)-sst))/hf_fac-swrad(i,j)  ! Be sure to read swrad with flux(3) before calling flux(4).
!              write(*,*) wtsurf(i,j)
            end do
          end do
!
          return
C
        case (3) ! Short wave radiation flux
C
          if (m.ne.rf_swrad) then
C
            rf_swrad = m
C
            filename = trim(pth_wrk)//trim(pth_flx)//
     $                 trim(pfx_dmn)//"roms_frc.nc"
            call check( nf90_open(filename, NF90_NOWRITE, ncid) )
            call check( nf90_inq_varid(ncid, "swrad", varid) )
            if (m.ne.1) then
              call check(nf90_get_var(ncid, varid, swradb, (/1,1,m-1/)))
              call check( nf90_get_var(ncid, varid, swradf, (/1,1,m/)) )
            else
              call check( nf90_get_var(ncid, varid, swradb, (/1,1,12/)))
              call check( nf90_get_var(ncid, varid, swradf, (/1,1,1/)) )
            end if
            call check( nf90_close(ncid) )

            ! Convert swrad from W/m^2 to K m/s.
            swradb = swradb/hf_fac
C
            write(*,*) "Read swrad:   ", m

          end if
!       Interpolate
          swrad = swradb+fac*(swradf-swradb)

          return
C
        case (52) ! momentum flux
C
          if (m.ne.rf_wsurf) then
C
            rf_wsurf = m
C
            filename = trim(pth_wrk)//trim(pth_flx)//
     $                 trim(pfx_dmn)//"roms_frc.nc"
            call check( nf90_open(filename, NF90_NOWRITE, ncid) )
            call check( nf90_inq_varid(ncid, "sustr", varid) )
            if (m.ne.1) then
              call check( nf90_get_var(ncid, varid, wusurfb(2:im,:)
     $                   , (/1,1,m-1/), (/imm1,jm,1/)) )
              call check( nf90_get_var(ncid, varid, wusurff(2:im,:)
     $                   , (/1,1,m/),   (/imm1,jm,1/)) )
            else
              call check( nf90_get_var(ncid, varid, wusurfb(2:im,:)
     $                   , (/1,1,12/), (/imm1,jm,1/)) )
              call check( nf90_get_var(ncid, varid, wusurff(2:im,:)
     $                   , (/1,1,1/),  (/imm1,jm,1/)) )
            end if
            call check( nf90_inq_varid(ncid, "svstr", varid) )
            if (m.ne.1) then
              call check( nf90_get_var(ncid, varid, wvsurfb(:,2:jm)
     $                   , (/1,1,m-1/), (/im,jmm1,1/)) )
              call check( nf90_get_var(ncid, varid, wvsurff(:,2:jm)
     $                   , (/1,1,m/),   (/im,jmm1,1/)) )
            else
              call check( nf90_get_var(ncid, varid, wvsurfb(:,2:jm)
     $                   , (/1,1,12/), (/im,jmm1,1/)) )
              call check( nf90_get_var(ncid, varid, wvsurff(:,2:jm)
     $                   , (/1,1,1/),  (/im,jmm1,1/)) )
            end if
            call check( nf90_close(ncid) )

            ! Fill in boundary values (and taper them, maybe?)
            wusurfb(:,1) = wusurfb(:,2)!*.5
            wusurff(:,1) = wusurff(:,2)!*.5
            wvsurfb(1,:) = wvsurfb(2,:)!*.5
            wvsurff(1,:) = wvsurff(2,:)!*.5

            ! Convert s?str (in N/m^2) to w?surf (in m^2/s^2).
            wusurfb = -wusurfb/rhoref
            wvsurfb = -wvsurfb/rhoref
            wusurff = -wusurff/rhoref
            wvsurff = -wvsurff/rhoref

            write(*,*) "Read w*surf:  ", m

          end if

          wusurf = wusurfb+fac*(wusurff-wusurfb)
          wvsurf = wvsurfb+fac*(wvsurff-wvsurfb)

          return

      end select
C
      return
C
        contains
          subroutine check(status)
            integer, intent ( in) :: status
!            if (DBG) write(*,*) status
            if(status /= nf90_noerr) then
              stop "Stopped"
            end if
          end subroutine check
C
      end
!
      subroutine flux_obs(idx)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Reads (if necessary) forcing dataset                *
C *          !TODO and interpolates in time.                           *
C *                                                                    *
C **********************************************************************
C
      use netcdf
      implicit none
C
      integer, intent (in) :: idx
      integer :: i,j
      integer :: tind
      character*2 :: mm
      include 'pom2k.c'
      real :: qq, dq, sst
      dimension :: qq(im, jm), dq(im, jm), sst(im, jm)
C
      select case (idx)

        case (4) ! Long wave radiation and other fluxes
C
C          tsplt_flx = 1
C
          tind = int(time/tsplt_wtsur)+1
          do
            if (tind.gt.365) then
              tind = tind-365
            else
              exit
            end if
          end do
C
          if (tind.ne.rf_wtsur) then
C
            rf_wtsur = tind
C
            ! The code below assumes a 365 days long year
            mm = '01'
            if (tind.lt.32) then
              mm  = '01'
            else
              if (tind.lt.60) then
                mm  = '02'
                tind = tind-31
              else
                if (tind.lt.91) then
                  mm  = '03'
                  tind = tind-59
                else
                  if (tind.lt.121) then
                    mm  = '04'
                    tind = tind-90
                  else
                    if (tind.lt.152) then
                      mm  = '05'
                      tind = tind-120
                    else
                      if (tind.lt.182) then
                        mm  = '06'
                        tind = tind-151
                      else
                        if (tind.lt.213) then
                          mm  = '07'
                          tind = tind-181
                        else
                          if (tind.lt.244) then
                            mm  = '08'
                            tind = tind-212
                          else
                            if (tind.lt.274) then
                              mm  = '09'
                              tind = tind-243
                           else
                              if (tind.lt.305) then
                                mm  = '10'
                                tind = tind-273
                              else
                                if (tind.lt.335) then
                                  mm  = '11'
                                  tind = tind-304
                                else
                                  if (tind.lt.366) then
                                    mm  = '12'
                                    tind = tind-334
                                  end if
                                end if
                              end if
                            end if
                          end if
                        end if
                      end if
                    end if
                  end if
                end if
              end if
            end if

            filename = trim(pth_wrk)//trim(pth_flx)//
     $                 trim(pfx_dmn)//"pom_flx_"//mm//".nc"
            call check( nf90_open(filename, NF90_NOWRITE, ncid) )
            call check( nf90_inq_varid(ncid, "hf", varid) )
            call check( nf90_get_var(ncid, varid, qq, (/1,1,tind/)))
            call check( nf90_inq_varid(ncid, "dq", varid) )
            call check( nf90_get_var(ncid, varid, dq, (/1,1,tind/)))
            call check( nf90_inq_varid(ncid, "sst", varid) )
            call check( nf90_get_var(ncid, varid, sst, (/1,1,tind/)))
            do i=1,im
              do j=1,jm
                wtsurf(i,j) = -(qq(i,j)+dq(i,j)
     $                        *(t(i,j,1)-sst(i,j)))
     $                        /4.1876e6-swrad(i,j)
C                write(*,*) wtsurf(i,j)
              end do
            end do
            call check( nf90_close(ncid) )
C
            write(*,*) "Read wtsurf:  ", tind, "@", mm
C
          end if
C
        case (3) ! Short wave radiation flux
C
          tind = int(time/tsplt_swrad)+1
          do
            if (tind.gt.365) then
              tind = tind-365
            else
              exit
            end if
          end do
C
          if (tind.ne.rf_swrad) then
C
            rf_swrad = tind
C
            ! The code below assumes a 365 days long year
            mm = '01'
            if (tind.lt.32) then
              mm  = '01'
            else
              if (tind.lt.60) then
                mm  = '02'
                tind = tind-31
              else
                if (tind.lt.91) then
                  mm  = '03'
                  tind = tind-59
                else
                  if (tind.lt.121) then
                    mm  = '04'
                    tind = tind-90
                  else
                    if (tind.lt.152) then
                      mm  = '05'
                      tind = tind-120
                    else
                      if (tind.lt.182) then
                        mm  = '06'
                        tind = tind-151
                      else
                        if (tind.lt.213) then
                          mm  = '07'
                          tind = tind-181
                        else
                          if (tind.lt.244) then
                            mm  = '08'
                            tind = tind-212
                          else
                            if (tind.lt.274) then
                              mm  = '09'
                              tind = tind-243
                           else
                              if (tind.lt.305) then
                                mm  = '10'
                                tind = tind-273
                              else
                                if (tind.lt.335) then
                                  mm  = '11'
                                  tind = tind-304
                                else
                                  if (tind.lt.366) then
                                    mm  = '12'
                                    tind = tind-334
                                  end if
                                end if
                              end if
                            end if
                          end if
                        end if
                      end if
                    end if
                  end if
                end if
              end if
            end if

            filename = trim(pth_wrk)//trim(pth_flx)//
     $                 trim(pfx_dmn)//"pom_flx_"//mm//".nc"
            call check( nf90_open(filename, NF90_NOWRITE, ncid) )
            call check( nf90_inq_varid(ncid, "swrad", varid) )
            call check( nf90_get_var(ncid, varid, swrad, (/1,1,tind/)) )
            do i=1,im
              do j=1,jm
                swrad(i,j) = -swrad(i,j)/4.1876e6
              end do
            end do
            call check( nf90_close(ncid) )
C
            write(*,*) "Read swrad:   ", tind, "@", mm
C
          end if
C
        case (51) ! Momentum flux
C
          tind = int(time/tsplt_wsurf)+1
C
          do
            if (tind.gt.365) then
              tind = tind-365
            else
              exit
            end if
          end do
C
          if (tind.ne.rf_wsurf) then
C
            rf_wsurf = tind
C
            if (tind.gt.1464) then
              tind = 1464
              write(*,*) "w*surf limit of 1464 records",
     $                   "has been reached *spam*"
            end if

            filename = "/home/pkharlamov/POM/pom_in/"//
     $                 "PACO-15_pom_wusur.nc"
            call check( nf90_open(filename, NF90_NOWRITE, ncid) )
            call check( nf90_inq_varid(ncid, "uwnd", varid) )
            call check( nf90_get_var(ncid, varid, wusurf, (/1,1,tind/)))
            call check( nf90_close(ncid) )
C
            filename = "/home/pkharlamov/POM/pom_in/"//
     $                 "PACO-15_pom_wvsur.nc"
            call check( nf90_open(filename, NF90_NOWRITE, ncid) )
            call check( nf90_inq_varid(ncid, "uwnd", varid) )
            call check( nf90_get_var(ncid, varid, wvsurf, (/1,1,tind/)))
            call check( nf90_close(ncid) )

            write(*,*) "Read w*surf:  ", tind

          end if

        case (52) ! Daily clim. momentum flux
C
          tsplt_wsurf = 1
          tind = int(time/tsplt_wsurf)+1
C
          do
            if (tind.gt.365) then
              tind = tind-365
            else
              exit
            end if
          end do
C
          if (tind.ne.rf_wsurf) then
C
            rf_wsurf = tind
C
            mm = '01'
            if (tind.lt.32) then
              mm  = '01'
            else
              if (tind.lt.60) then
                mm  = '02'
                tind = tind-31
              else
                if (tind.lt.91) then
                  mm  = '03'
                  tind = tind-59
                else
                  if (tind.lt.121) then
                    mm  = '04'
                    tind = tind-90
                  else
                    if (tind.lt.152) then
                      mm  = '05'
                      tind = tind-120
                    else
                      if (tind.lt.182) then
                        mm  = '06'
                        tind = tind-151
                      else
                        if (tind.lt.213) then
                          mm  = '07'
                          tind = tind-181
                        else
                          if (tind.lt.244) then
                            mm  = '08'
                            tind = tind-212
                          else
                            if (tind.lt.274) then
                              mm  = '09'
                              tind = tind-243
                           else
                              if (tind.lt.305) then
                                mm  = '10'
                                tind = tind-273
                              else
                                if (tind.lt.335) then
                                  mm  = '11'
                                  tind = tind-304
                                else
                                  if (tind.lt.366) then
                                    mm  = '12'
                                    tind = tind-334
                                  end if
                                end if
                              end if
                            end if
                          end if
                        end if
                      end if
                    end if
                  end if
                end if
              end if
            end if
C
            filename = trim(pth_wrk)//trim(pth_flx)//
     $                 trim(pfx_dmn)//"pom_flx_"//mm//".nc"
            call check( nf90_open(filename, NF90_NOWRITE, ncid) )
            call check( nf90_inq_varid(ncid, "wusurf", varid) )
            call check( nf90_get_var(ncid, varid, wusurf, (/1,1,tind/)))
            call check( nf90_inq_varid(ncid, "wvsurf", varid) )
            call check( nf90_get_var(ncid, varid, wvsurf, (/1,1,tind/)))
            call check( nf90_close(ncid) )

            write(*,*) "Read w*surf:  ", tind, "@", mm

          end if

      end select
C
      return
C
        contains
          subroutine check(status)
            integer, intent ( in) :: status
!            if (DBG) write(*,*) status
            if(status /= nf90_noerr) then
              stop "Stopped"
            end if
          end subroutine check
C
      end
C
      subroutine bry_obs(idx)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Reads (if necessary) boundary conditions (t,s)      *
C *          !TODO and interpolates in time.                           *
C *                                                                    *
C **********************************************************************
C
      use netcdf
      implicit none
C
      integer, intent ( in) :: idx
      integer :: i,j,k
      integer :: tind
      character*2 :: mm

      include 'pom2k.c'

      real :: datXZm(im,kbm1), datXmZm(imm1,kbm1)
      real :: v3f(im,jm,kbm1)

C
      select case (idx)
C
          case (2) ! elevation
C
          tind = int(time/tsplt_el)+1
          do
            if (tind.gt.365) then
              tind = tind-365
            else
              exit
            end if
          end do
C
          if (tind.ne.rf_el) then
C
            rf_el = tind
C
            mm = '01'
            if (tind.lt.32) then
              mm  = '01'
            else
              if (tind.lt.60) then
                mm  = '02'
                tind = tind-31
              else
                if (tind.lt.91) then
                  mm  = '03'
                  tind = tind-59
                else
                  if (tind.lt.121) then
                    mm  = '04'
                    tind = tind-90
                  else
                    if (tind.lt.152) then
                      mm  = '05'
                      tind = tind-120
                    else
                      if (tind.lt.182) then
                        mm  = '06'
                        tind = tind-151
                      else
                        if (tind.lt.213) then
                          mm  = '07'
                          tind = tind-181
                        else
                          if (tind.lt.244) then
                            mm  = '08'
                            tind = tind-212
                          else
                            if (tind.lt.274) then
                              mm  = '09'
                              tind = tind-243
                            else
                              if (tind.lt.305) then
                                mm  = '10'
                                tind = tind-273
                              else
                                if (tind.lt.335) then
                                  mm  = '11'
                                  tind = tind-304
                                else
                                  if (tind.lt.366) then
                                    mm  = '12'
                                    tind = tind-334
                                  end if
                                end if
                              end if
                            end if
                          end if
                        end if
                      end if
                    end if
                  end if
                end if
              end if
            end if

            filename = trim(pth_wrk)//trim(pth_bry)//
     $                 trim(pfx_dmn)//"pom_bry_"//mm//".nc"
            call check( nf90_open(filename, NF90_NOWRITE, ncid) )
C
!            call check( nf90_inq_varid(ncid, "Z_N", varid) )
!            call check( nf90_get_var(ncid, varid, eln,(/1,tind/)))
!            call check( nf90_inq_varid(ncid, "Z_E", varid) )
!            call check( nf90_get_var(ncid, varid, ele,(/1,tind/)) )
            call check( nf90_inq_varid(ncid, "Z_S", varid) )
            call check( nf90_get_var(ncid, varid, els,(/1,tind/)) )
!            call check( nf90_inq_varid(ncid, "Z_W", varid) )
!            call check( nf90_get_var(ncid, varid, elw,(/1,tind/)) )
            call check( nf90_close(ncid) )
C
!            do i=1,im
!              elf( i, jmm1) = eln(i)
!              elf( i, 2   ) = els(i)
!            end do
!            do j=1,jm
!              elf(imm1, j ) = elw(i)
!              elf(2   , j ) = ele(i)
!            end do
C
            write(*,*) "Read el BCs:  ", tind, "@", mm
C
          end if
C
        case (3) ! u and v
C
          tind = int(time/tsplt_uv)+1
          do
            if (tind.gt.365) then
              tind = tind-365
            else
              exit
            end if
          end do
C
          if (tind.ne.rf_uv) then
C
            rf_uv = tind
C
            mm = '01'
            if (tind.lt.32) then
              mm  = '01'
            else
              if (tind.lt.60) then
                mm  = '02'
                tind = tind-31
              else
                if (tind.lt.91) then
                  mm  = '03'
                  tind = tind-59
                else
                  if (tind.lt.121) then
                    mm  = '04'
                    tind = tind-90
                  else
                    if (tind.lt.152) then
                      mm  = '05'
                      tind = tind-120
                    else
                      if (tind.lt.182) then
                        mm  = '06'
                        tind = tind-151
                      else
                        if (tind.lt.213) then
                          mm  = '07'
                          tind = tind-181
                        else
                          if (tind.lt.244) then
                            mm  = '08'
                            tind = tind-212
                          else
                            if (tind.lt.274) then
                              mm  = '09'
                              tind = tind-243
                            else
                              if (tind.lt.305) then
                                mm  = '10'
                                tind = tind-273
                              else
                                if (tind.lt.335) then
                                  mm  = '11'
                                  tind = tind-304
                                else
                                  if (tind.lt.366) then
                                    mm  = '12'
                                    tind = tind-334
                                  end if
                                end if
                              end if
                            end if
                          end if
                        end if
                      end if
                    end if
                  end if
                end if
              end if
            end if

            filename = trim(pth_wrk)//trim(pth_bry)//
     $                 trim(pfx_dmn)//"pom_bry_"//mm//".nc"
            call check( nf90_open(filename, NF90_NOWRITE, ncid) )
C
            call check( nf90_inq_varid(ncid, "U_S", varid) )
            call check( nf90_get_var(ncid, varid, datXmZm,(/1,1,tind/),
     $                  (/imm1,kbm1,1/)))
            call check( nf90_inq_varid(ncid, "V_S", varid) )
            call check( nf90_get_var(ncid, varid, vbs,(/1,1,tind/),
     $                  (/im,kbm1,1/)))
!
!     Shift u-data to the right by 1 cell
!     
            do i=1,imm1
              do k=1,kbm1
                ubs(i+1,k) = datXmZm(i,k)
              end do
              ubs(i+1,kb) = 0.
            end do
!
!     Blank the unfilled values for safety
!
            do k=1,kb
              ubs(1,k) = ubs(2,k)
            end do
!
!       Calculate vertical mean BCs for 2D mode
!
            do i=1,im
              uabs(i) = 0.
              vabs(i) = 0.
              do k=1,kb
                uabs(i) = uabs(i)+ubs(i,k)*dz(k)
                vabs(i) = vabs(i)+vbs(i,k)*dz(k)
              end do
            end do
            call check( nf90_close(ncid) )

            write(*,*) "Read vel BCs: ", tind, "@", mm
C
          end if
C
        case (4) ! temperature and salinity
C
          tind = int(time/tsplt_ts)+1
          do
            if (tind.gt.365) then
              tind = tind-365
            else
              exit
            end if
          end do
C
          if (tind.ne.rf_ts) then
C
            rf_ts = tind
C
            mm = '01'
            if (tind.lt.32) then
              mm  = '01'
            else
              if (tind.lt.60) then
                mm  = '02'
                tind = tind-31
              else
                if (tind.lt.91) then
                  mm  = '03'
                  tind = tind-59
                else
                  if (tind.lt.121) then
                    mm  = '04'
                    tind = tind-90
                  else
                    if (tind.lt.152) then
                      mm  = '05'
                      tind = tind-120
                    else
                      if (tind.lt.182) then
                        mm  = '06'
                        tind = tind-151
                      else
                        if (tind.lt.213) then
                          mm  = '07'
                          tind = tind-181
                        else
                          if (tind.lt.244) then
                            mm  = '08'
                            tind = tind-212
                          else
                            if (tind.lt.274) then
                              mm  = '09'
                              tind = tind-243
                            else
                              if (tind.lt.305) then
                                mm  = '10'
                                tind = tind-273
                              else
                                if (tind.lt.335) then
                                  mm  = '11'
                                  tind = tind-304
                                else
                                  if (tind.lt.366) then
                                    mm  = '12'
                                    tind = tind-334
                                  end if
                                end if
                              end if
                            end if
                          end if
                        end if
                      end if
                    end if
                  end if
                end if
              end if
            end if
C
!           filename = trim(pth_wrk)//trim(pth_bry)//
!     $                 trim(pfx_dmn)//"pom_clm_"//mm//".nc"
!            call check( nf90_open(filename, NF90_NOWRITE, ncid) )
!            call check( nf90_inq_varid(ncid, "t", varid) )
!            call check( nf90_get_var(ncid, varid, v3f,(/1,1,1,tind/)) )
!            do i=1,im
!              do k=1,kbm1
!                tbs(i,k) = v3f(i, 1,k)
!                tbn(i,k) = v3f(i,jm,k)
!              end do
!              tbs(i,kb) = tbs(i,kbm1)
!              tbn(i,kb) = tbn(i,kbm1)
!            end do
!            do j=1,jm
!              do k=1,kbm1
!                tbw(j,k) = v3f( 1,j,k)
!                tbe(j,k) = v3f(im,j,k)
!              end do
!              tbw(j,kb) = tbw(j,kbm1)
!              tbe(j,kb) = tbe(j,kbm1)
!            end do
!            call check( nf90_close(ncid) )
C
C           Test to understand where bry t and s are from.
C           Override climatological with boundary
            filename = trim(pth_wrk)//trim(pth_bry)//
     $                 trim(pfx_dmn)//"pom_bry_"//mm//".nc"
            call check( nf90_open(filename, NF90_NOWRITE, ncid) )
            call check( nf90_inq_varid(ncid, "T_S", varid) )
            call check( nf90_get_var(ncid, varid, tbs,(/1,1,tind/),
     $                  (/im,kbm1,1/)))
!            call check( nf90_inq_varid(ncid, "T_N", varid) )
!            call check( nf90_get_var(ncid, varid, tmpXZm,(/1,1,tind/)) )
!            do i=1,im
!              do k=1,kbm1
!                tbn(i,k) = tmpXZm(i,k)
!              end do
!              tbn(i,kb) = tmpXZm(i,kbm1)
!            end do
!            call check( nf90_inq_varid(ncid, "T_E", varid) )
!            call check( nf90_get_var(ncid, varid, tbe,
!     $                              (/1,1,time_index/) ) )
!            call check( nf90_inq_varid(ncid, "T_W", varid) )
!            call check( nf90_get_var(ncid, varid, tmpYZm,(/1,1,tind/)) )
!            do j=1,im
!              do k=1,kbm1
!                tbw(j,k) = tmpYZm(j,k)
!              end do
!              tbw(j,kb) = tmpYZm(j,kbm1)
!            end do
            call check( nf90_inq_varid(ncid, "S_S", varid) )
            call check( nf90_get_var(ncid, varid, sbs,(/1,1,tind/),
     $                  (/im,kbm1,1/)))
!            call check( nf90_inq_varid(ncid, "S_N", varid) )
!            call check( nf90_get_var(ncid, varid, tmpXZm,(/1,1,tind/)) )
!            do i=1,im
!              do k=1,kbm1
!                sbn(i,k) = tmpXZm(i,k)
!              end do
!              sbn(i,kb) = tmpXZm(i,kbm1)
!            end do
!            call check( nf90_inq_varid(ncid, "S_E", varid) )
!            call check( nf90_get_var(ncid, varid, sbe,
!     $                              (/1,1,time_index/) ) )
!            call check( nf90_inq_varid(ncid, "S_W", varid) )
!            call check( nf90_get_var(ncid, varid, tmpYZm,(/1,1,tind/)) )
!            do j=1,im
!              do k=1,kbm1
!                sbw(j,k) = tmpYZm(j,k)
!              end do
!              sbw(j,kb) = tmpYZm(j,kbm1)
!            end do
            call check( nf90_close(ncid) )

            write(*,*) "Read TS BCs:  ", tind, "@", mm
C
          end if
C
          case (43)

            tind = int(time/tsplt_ts)+1
            do
              if (tind.gt.365) then
                tind = tind-365
              else
                exit
              end if
            end do
C
            if (tind.ne.rf_sts) then
C
              rf_sts = tind
C
              mm = '01'
              if (tind.lt.32) then
                mm  = '01'
              else
                if (tind.lt.60) then
                  mm  = '02'
                  tind = tind-31
                else
                  if (tind.lt.91) then
                    mm  = '03'
                    tind = tind-59
                  else
                    if (tind.lt.121) then
                      mm  = '04'
                      tind = tind-90
                    else
                      if (tind.lt.152) then
                        mm  = '05'
                        tind = tind-120
                      else
                        if (tind.lt.182) then
                          mm  = '06'
                          tind = tind-151
                        else
                          if (tind.lt.213) then
                            mm  = '07'
                            tind = tind-181
                          else
                            if (tind.lt.244) then
                              mm  = '08'
                              tind = tind-212
                            else
                              if (tind.lt.274) then
                                mm  = '09'
                                tind = tind-243
                              else
                                if (tind.lt.305) then
                                  mm  = '10'
                                  tind = tind-273
                                else
                                  if (tind.lt.335) then
                                    mm  = '11'
                                    tind = tind-304
                                  else
                                    if (tind.lt.366) then
                                      mm  = '12'
                                      tind = tind-334
                                    end if
                                  end if
                                end if
                              end if
                            end if
                          end if
                        end if
                      end if
                    end if
                  end if
                end if
              end if

              filename = trim(pth_wrk)//trim(pth_bry)//
     $                 trim(pfx_dmn)//"pom_flx_"//mm//".nc"
              call check( nf90_open(filename, NF90_NOWRITE, ncid) )
              call check( nf90_inq_varid(ncid, "sst", varid) )
              call check( nf90_get_var(ncid, varid, tsurf,(/1,1,tind/)))
              call check( nf90_inq_varid(ncid, "sss", varid) )
              call check( nf90_get_var(ncid, varid, ssurf,(/1,1,tind/)))
              call check( nf90_close(ncid) )
              write(*,*) "Read ST/S BCs:", tind, "@", mm

            end if
C
      end select
C
      return
C
        contains
          subroutine check(status)
            integer, intent ( in) :: status
!            if (DBG) write(*,*) status
            if(status /= nf90_noerr) then
              stop "Stopped"
            end if
          end subroutine check
C
      end
!
      subroutine bry(idx)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Reads (if necessary) boundary conditions (t,s)      *
C *                and interpolates in time.                           *
C *                                                                    *
C **********************************************************************
C
      use netcdf
      implicit none
C
      integer, intent ( in) :: idx
      integer :: i,j,k

      include 'pom2k.c'

C
      select case (idx)
C
        case (2) ! elevation
!
          if (m.ne.rf_el) then
!     If we move to the next month...
            rf_el = m

            filename = trim(pth_wrk)//trim(pth_bry)//
     $                 trim(pfx_dmn)//"roms_bry.nc"
            call check( nf90_open(filename, NF90_NOWRITE, ncid) )
            call check( nf90_inq_varid(ncid, "zeta_south", varid) )
!     ...read neccessary fields
!     South:
            if (m.ne.1) then
              call check( nf90_get_var(ncid, varid, elsb,
     $                          (/1,m-1/), (/im,1/)) )
              call check( nf90_get_var(ncid, varid, elsf,
     $                          (/1,m/), (/im,1/)) )
            else
              call check( nf90_get_var(ncid, varid, elsb,
     $                          (/1,12/),(/im,1/)) )
              call check( nf90_get_var(ncid, varid, elsf,
     $                          (/1,1/), (/im,1/)) )
            end if

            call check( nf90_close(ncid) )

            write(*,*) "Read el BCs:  ", m

          end if
!     Perform interpolation
!     South:
          ! FSM is already defined! Feel free to apply it! (But DON'T!!! Really! Don't! It just boundary conditions.)
          els(:) = (elsb(:)+fac*(elsf(:)-elsb(:)))*fsm(:,1)

          return
C
        case (3) ! u and v
C
          if (m.ne.rf_uv) then
!        If we move to the next month...
            rf_uv = m

            filename = trim(pth_wrk)//trim(pth_bry)//
     $                 trim(pfx_dmn)//"roms_bry.nc"
            call check( nf90_open(filename, NF90_NOWRITE, ncid) )
!        ...read neccessary fields
!     South:
            call check( nf90_inq_varid(ncid, "u_south", varid) )
            if (m.ne.1) then
              call check( nf90_get_var(ncid, varid,
     $         ubsb(2:im,kbm1:1:-1),(/1,1,m-1/),(/imm1,kbm1,1/)) )
              call check( nf90_get_var(ncid, varid,
     $         ubsf(2:im,kbm1:1:-1),(/1,1,m/),  (/imm1,kbm1,1/)) )
            else
              call check( nf90_get_var(ncid, varid,
     $         ubsb(2:im,kbm1:1:-1),(/1,1,12/),(/imm1,kbm1,1/)) )
              call check( nf90_get_var(ncid, varid,
     $         ubsf(2:im,kbm1:1:-1),(/1,1,1/), (/imm1,kbm1,1/)) )
            end if
            call check( nf90_inq_varid(ncid, "v_south", varid) )
            if (m.ne.1) then
              call check( nf90_get_var(ncid, varid,
     $         vbsb(:,kbm1:1:-1),(/1,1,m-1/), (/im,kbm1,1/)) )
              call check( nf90_get_var(ncid, varid,
     $         vbsf(:,kbm1:1:-1),(/1,1,m/),   (/im,kbm1,1/)) )
            else
              call check( nf90_get_var(ncid, varid,
     $         vbsb(:,kbm1:1:-1),(/1,1,12/),(/im,kbm1,1/)) )
              call check( nf90_get_var(ncid, varid,
     $         vbsf(:,kbm1:1:-1),(/1,1,1/), (/im,kbm1,1/)) )
            end if

            call check( nf90_close(ncid) )

            write(*,*) "Read vel BCs: ", m

          end if
!     Perform interpolation
!     South:
!          do k=1,kbm1   ! FSM is already defined! Feel free to apply it!
            ubs(:,1:kbm1) = ubsb(:,1:kbm1)
     $                     +fac*(ubsf(:,1:kbm1)-ubsb(:,1:kbm1))
            vbs(:,1:kbm1) = vbsb(:,1:kbm1)
     $                     +fac*(vbsf(:,1:kbm1)-vbsb(:,1:kbm1))
!          do k=1,kbm1              ! This doesn't work, since it completely wipes out u
!            ubs(:,k) = ubs(:,k)*dum(:,1)
!            vbs(:,k) = vbs(:,k)*dvm(:,2) ! index 1, maybe?
!          end do
!          end do

!     Calculate vertical mean BCs for 2D mode
!
          do i=1,im
            uabs(i) = 0.
            vabs(i) = 0.
            do k=1,kb
              uabs(i) = uabs(i)+ubs(i,k)*dz(k)
              vabs(i) = vabs(i)+vbs(i,k)*dz(k)
            end do
          end do

          return
C
        case (4) ! temperature and salinity
C
          if (m.ne.rf_ts) then
!        If we move to the next month...
            rf_ts = m

            filename = trim(pth_wrk)//trim(pth_bry)//
     $                 trim(pfx_dmn)//"roms_clm.nc"
            call check( nf90_open(filename, NF90_NOWRITE, ncid) )
            call check( nf90_inq_varid(ncid, "temp", varid) )
!    Read climatological BCs
            if (m.ne.1) then
!            south
              call check( nf90_get_var(ncid, varid,
     $         tbsb(:,kbm1:1:-1),(/1,1,1,m-1/), (/im,1,kbm1,1/)) )
              call check( nf90_get_var(ncid, varid,
     $         tbsf(:,kbm1:1:-1),(/1,1,1,m/),   (/im,1,kbm1,1/)) )
!            north
              call check( nf90_get_var(ncid, varid,
     $         tbnb(:,kbm1:1:-1),(/1,jm,1,m-1/), (/im,1,kbm1,1/)) )
              call check( nf90_get_var(ncid, varid,
     $         tbnf(:,kbm1:1:-1),(/1,jm,1,m/),   (/im,1,kbm1,1/)) )
!            west
              call check( nf90_get_var(ncid, varid,
     $         tbwb(:,kbm1:1:-1),(/1,1,1,m-1/), (/1,jm,kbm1,1/)) )
              call check( nf90_get_var(ncid, varid,
     $         tbwf(:,kbm1:1:-1),(/1,1,1,m/),   (/1,jm,kbm1,1/)) )
!            east
              call check( nf90_get_var(ncid, varid,
     $         tbeb(:,kbm1:1:-1),(/im,1,1,m-1/), (/1,jm,kbm1,1/)) )
              call check( nf90_get_var(ncid, varid,
     $         tbef(:,kbm1:1:-1),(/im,1,1,m/),   (/1,jm,kbm1,1/)) )
            else
!            south
              call check( nf90_get_var(ncid, varid,
     $         tbsb(:,kbm1:1:-1),(/1,1,1,12/),(/im,1,kbm1,1/)) )
              call check( nf90_get_var(ncid, varid,
     $         tbsf(:,kbm1:1:-1),(/1,1,1,1/), (/im,1,kbm1,1/)) )
!            north
              call check( nf90_get_var(ncid, varid,
     $         tbnb(:,kbm1:1:-1),(/1,jm,1,12/),(/im,1,kbm1,1/)) )
              call check( nf90_get_var(ncid, varid,
     $         tbnf(:,kbm1:1:-1),(/1,jm,1,1/), (/im,1,kbm1,1/)) )
!            east
              call check( nf90_get_var(ncid, varid,
     $         tbeb(:,kbm1:1:-1),(/1,1,1,12/),(/1,jm,kbm1,1/)) )
              call check( nf90_get_var(ncid, varid,
     $         tbef(:,kbm1:1:-1),(/1,1,1,1/), (/1,jm,kbm1,1/)) )
!            west
              call check( nf90_get_var(ncid, varid,
     $         tbeb(:,kbm1:1:-1),(/im,1,1,12/),(/1,jm,kbm1,1/)) )
              call check( nf90_get_var(ncid, varid,
     $         tbef(:,kbm1:1:-1),(/im,1,1,1/), (/1,jm,kbm1,1/)) )
            end if
!
            call check( nf90_inq_varid(ncid, "salt", varid) )
!    Read climatological salinity BCs
            if (m.ne.1) then
!            south
              call check( nf90_get_var(ncid, varid,
     $         sbsb(:,kbm1:1:-1),(/1,1,1,m-1/), (/im,1,kbm1,1/)) )
              call check( nf90_get_var(ncid, varid,
     $         sbsf(:,kbm1:1:-1),(/1,1,1,m/),   (/im,1,kbm1,1/)) )
!            north
              call check( nf90_get_var(ncid, varid,
     $         sbnb(:,kbm1:1:-1),(/1,jm,1,m-1/), (/im,1,kbm1,1/)) )
              call check( nf90_get_var(ncid, varid,
     $         sbnf(:,kbm1:1:-1),(/1,jm,1,m/),   (/im,1,kbm1,1/)) )
!            west
              call check( nf90_get_var(ncid, varid,
     $         sbwb(:,kbm1:1:-1),(/1,1,1,m-1/), (/1,jm,kbm1,1/)) )
              call check( nf90_get_var(ncid, varid,
     $         sbwf(:,kbm1:1:-1),(/1,1,1,m/),   (/1,jm,kbm1,1/)) )
!            east
              call check( nf90_get_var(ncid, varid,
     $         sbeb(:,kbm1:1:-1),(/im,1,1,m-1/), (/1,jm,kbm1,1/)) )
              call check( nf90_get_var(ncid, varid,
     $         sbef(:,kbm1:1:-1),(/im,1,1,m/),   (/1,jm,kbm1,1/)) )
            else
!            south
              call check( nf90_get_var(ncid, varid,
     $         sbsb(:,kbm1:1:-1),(/1,1,1,12/),(/im,1,kbm1,1/)) )
              call check( nf90_get_var(ncid, varid,
     $         sbsf(:,kbm1:1:-1),(/1,1,1,1/), (/im,1,kbm1,1/)) )
!            north
              call check( nf90_get_var(ncid, varid,
     $         sbnb(:,kbm1:1:-1),(/1,jm,1,12/),(/im,1,kbm1,1/)) )
              call check( nf90_get_var(ncid, varid,
     $         sbnf(:,kbm1:1:-1),(/1,jm,1,1/), (/im,1,kbm1,1/)) )
!            east
              call check( nf90_get_var(ncid, varid,
     $         sbeb(:,kbm1:1:-1),(/1,1,1,12/),(/1,jm,kbm1,1/)) )
              call check( nf90_get_var(ncid, varid,
     $         sbef(:,kbm1:1:-1),(/1,1,1,1/), (/1,jm,kbm1,1/)) )
!            west
              call check( nf90_get_var(ncid, varid,
     $         sbeb(:,kbm1:1:-1),(/im,1,1,12/),(/1,jm,kbm1,1/)) )
              call check( nf90_get_var(ncid, varid,
     $         sbef(:,kbm1:1:-1),(/im,1,1,1/), (/1,jm,kbm1,1/)) )
            end if
!
            call check( nf90_close(ncid) )
!
            if (.false.) then   ! Disable it for now
            filename = trim(pth_wrk)//trim(pth_bry)//
     $                 trim(pfx_dmn)//"roms_bry.nc"
            call check( nf90_open(filename, NF90_NOWRITE, ncid) )
!        ...read neccessary fields
!     South:
            call check( nf90_inq_varid(ncid, "temp_south", varid) )
            if (m.ne.1) then
              call check( nf90_get_var(ncid, varid,
     $         tbsb(:,kbm1:1:-1),(/1,1,m-1/), (/im,kbm1,1/)) )
              call check( nf90_get_var(ncid, varid,
     $         tbsf(:,kbm1:1:-1),(/1,1,m/),   (/im,kbm1,1/)) )
            else
              call check( nf90_get_var(ncid, varid,
     $         tbsb(:,kbm1:1:-1),(/1,1,12/),(/im,kbm1,1/)) )
              call check( nf90_get_var(ncid, varid,
     $         tbsf(:,kbm1:1:-1),(/1,1,1/), (/im,kbm1,1/)) )
            end if
            call check( nf90_inq_varid(ncid, "salt_south", varid) )
            if (m.ne.1) then
              call check( nf90_get_var(ncid, varid,
     $         sbsb(:,kbm1:1:-1),(/1,1,m-1/), (/im,kbm1,1/)) )
              call check( nf90_get_var(ncid, varid,
     $         sbsf(:,kbm1:1:-1),(/1,1,m/),   (/im,kbm1,1/)) )
            else
              call check( nf90_get_var(ncid, varid,
     $         sbsb(:,kbm1:1:-1),(/1,1,12/),(/im,kbm1,1/)) )
              call check( nf90_get_var(ncid, varid,
     $         sbsf(:,kbm1:1:-1),(/1,1,1/), (/im,kbm1,1/)) )
            end if

            call check( nf90_close(ncid) )
            end if

            write(*,*) "Read TS BCs:  ", m

          end if
!     Perform interpolation
!     South:
            tbs(:,1:kbm1) = tbsb(:,1:kbm1)
     $                     +fac*(tbsf(:,1:kbm1)-tbsb(:,1:kbm1))
            sbs(:,1:kbm1) = sbsb(:,1:kbm1)
     $                     +fac*(sbsf(:,1:kbm1)-sbsb(:,1:kbm1))
!     North:
            tbn(:,1:kbm1) = tbnb(:,1:kbm1)
     $                     +fac*(tbnf(:,1:kbm1)-tbnb(:,1:kbm1))
            sbn(:,1:kbm1) = sbnb(:,1:kbm1)
     $                     +fac*(sbnf(:,1:kbm1)-sbnb(:,1:kbm1))
!     East:
            tbe(:,1:kbm1) = tbeb(:,1:kbm1)
     $                     +fac*(tbef(:,1:kbm1)-tbeb(:,1:kbm1))
            sbe(:,1:kbm1) = sbeb(:,1:kbm1)
     $                     +fac*(sbef(:,1:kbm1)-sbeb(:,1:kbm1))
!     West:
            tbw(:,1:kbm1) = tbwb(:,1:kbm1)
     $                     +fac*(tbwf(:,1:kbm1)-tbwb(:,1:kbm1))
            sbw(:,1:kbm1) = sbwb(:,1:kbm1)
     $                     +fac*(sbwf(:,1:kbm1)-sbwb(:,1:kbm1))
          return
C
        case (43)

          if (m.ne.rf_sts) then
!        If we move to the next month...
            rf_sts = m

            filename = trim(pth_wrk)//trim(pth_bry)//
     $                 trim(pfx_dmn)//"roms_frc.nc"
            call check( nf90_open(filename, NF90_NOWRITE, ncid) )
!        ...read neccessary fields
!     South:
            call check( nf90_inq_varid(ncid, "SST", varid) )
            if (m.ne.1) then
              call check( nf90_get_var(ncid, varid,
     $         tsurfb,(/1,1,m-1/), (/im,jm,1/)) )
              call check( nf90_get_var(ncid, varid,
     $         tsurff,(/1,1,m  /), (/im,jm,1/)) )
            else
              call check( nf90_get_var(ncid, varid,
     $         tsurfb,(/1,1,12/), (/im,jm,1/)) )
              call check( nf90_get_var(ncid, varid,
     $         tsurff,(/1,1, 1/), (/im,jm,1/)) )
            end if
            call check( nf90_inq_varid(ncid, "SSS", varid) )
            if (m.ne.1) then
              call check( nf90_get_var(ncid, varid,
     $         ssurfb,(/1,1,m-1/), (/im,jm,1/)) )
              call check( nf90_get_var(ncid, varid,
     $         ssurff,(/1,1,m  /), (/im,jm,1/)) )
            else
              call check( nf90_get_var(ncid, varid,
     $         ssurfb,(/1,1,12/), (/im,jm,1/)) )
              call check( nf90_get_var(ncid, varid,
     $         ssurff,(/1,1, 1/), (/im,jm,1/)) )
            end if

            call check( nf90_close(ncid) )

            write(*,*) "Read ST/S BCs:", m

          end if
!     Perform interpolation
!     South:
!          do k=1,kbm1   ! FSM is already defined! Feel free to apply it!
            tsurf = (tsurfb+fac*(tsurff-tsurfb))*fsm
            ssurf = (ssurfb+fac*(ssurff-ssurfb))*fsm
!          end do
          return
!
        case (45)

          if (m.ne.rf_sts) then
!        If we move to the next month...
            rf_sts = m

            filename = trim(pth_wrk)//trim(pth_bry)//
     $                 trim(pfx_dmn)//"roms_clm.nc"
            call check(nf90_open(filename, NF90_NOWRITE, ncid))
!        ...read neccessary fields
!     South:
            call check( nf90_inq_varid(ncid, "temp", varid) )
            if (m.ne.1) then
              call check( nf90_get_var(ncid, varid, tsurfb,
     $         (/1,1,kbm1,m-1/), (/im,jm,1,1/)) )
              call check( nf90_get_var(ncid, varid, tsurff,
     $         (/1,1,kbm1,m  /), (/im,jm,1,1/)) )
            else
              call check( nf90_get_var(ncid, varid,
     $         tsurfb,(/1,1,kbm1,12/), (/im,jm,1,1/)) )
              call check( nf90_get_var(ncid, varid,
     $         tsurff,(/1,1,kbm1, 1/), (/im,jm,1,1/)) )
            end if
            call check( nf90_inq_varid(ncid, "salt", varid) )
            if (m.ne.1) then
            call check( nf90_get_var(ncid, varid,
     $         ssurfb,(/1,1,kbm1,m-1/), (/im,jm,1,1/)) )
            call check( nf90_get_var(ncid, varid,
     $         ssurff,(/1,1,kbm1,m  /), (/im,jm,1,1/)) )
            else
              call check( nf90_get_var(ncid, varid,
     $         ssurfb,(/1,1,kbm1,12/), (/im,jm,1,1/)) )
              call check( nf90_get_var(ncid, varid,
     $         ssurff,(/1,1,kbm1, 1/), (/im,jm,1,1/)) )
            end if

            call check( nf90_close(ncid) )

            write(*,*) "Read 1l TS BCs:", m

          end if
!     Perform interpolation
!     South:
!          do k=1,kbm1   ! FSM is already defined! Feel free to apply it!
          tsurf = (tsurfb+fac*(tsurff-tsurfb))*fsm
          ssurf = (ssurfb+fac*(ssurff-ssurfb))*fsm
!          end do
          return
C
      end select
C
      return
C
        contains
          subroutine check(status)
            integer, intent ( in) :: status
!            if (DBG) write(*,*) status
            if(status /= nf90_noerr) then
              stop "Stopped"
            end if
          end subroutine check
C
      end
c     include 'pom2k.n'                                       ! *netCDF*
C
C     End of source code
C
C-----------------------------------------------------------------------
C
      subroutine debug_write_txt_xy(var, caption, step)

        implicit none
        include 'pom2k.c'

        real, intent (in)             :: var(im,jm)
        character(len=*), intent (in) :: caption
        integer, intent (in)          :: step
        character*4                   :: path, prefix
        character*64                  :: fmt

        integer                       :: lft,rht,top,bot,hst,vst

        lft = 1
        rht = im
        top = jm
        bot = 1
        hst = 1
        vst =-1
        write(fmt,*) '('
        write(fmt,'(i3)') (rht-lft+1)
        write(fmt,*) '(F10.5,";"))'
        write(path,'(I4.4)') step
        call system('mkdir -p '//trim(pth_wrk)
     $                         //trim(pth_dbg)//trim(path))
        write(prefix, '(I3.3,".")') dbg_seq_i
        dbg_seq_i = dbg_seq_i + 1
        open(42, file=trim(pth_wrk)//trim(pth_dbg)
     $                //path//"/"//trim(prefix)
     $                //trim(pfx_dbg)//"_"//caption//".csv")
        write(42, fmt) var(lft:rht:hst,top:bot:vst)
        close(42)

      end
!!!!
      subroutine debug_write_xy(var, caption, step)

        implicit none
        include 'pom2k.c'

        real, intent (in)             :: var(im,jm)
        character(len=*), intent (in) :: caption
        integer, intent (in)          :: step
        character*4                   :: path, prefix

        write(path,'(I4.4)') step
        call system('mkdir -p '//trim(pth_wrk)
     $                         //trim(pth_dbg)//trim(path))
        write(prefix, '(I3.3,".")') dbg_seq_i
        dbg_seq_i = dbg_seq_i + 1
        open(42, file=trim(pth_wrk)//trim(pth_dbg)
     $                //path//"/"//trim(prefix)
     $                //trim(pfx_dbg)//"_"//caption
     $          ,form='unformatted')
        write(42) var
        close(42)

      end

      subroutine debug_write_xz(var, caption, step)

        implicit none
        include 'pom2k.c'

        real, intent (in)             :: var(im,kb)
        character(len=*), intent (in) :: caption
        integer, intent (in)          :: step
        character*4                   :: path, prefix

        write(path,'(I4.4)') step
        call system('mkdir -p '//trim(pth_wrk)
     $                         //trim(pth_dbg)//trim(path))
        write(prefix, '(I3.3,".")') dbg_seq_i
        dbg_seq_i = dbg_seq_i + 1
        open(42, file=trim(pth_wrk)//trim(pth_dbg)
     $                //path//"/"//trim(prefix)
     $                //trim(pfx_dbg)//"_"//caption
     $          ,form='unformatted')
        write(42) var
        close(42)

      end

      subroutine debug_write_x(var, caption, step)

        implicit none
        include 'pom2k.c'

        real, intent (in)             :: var(im)
        character(len=*), intent (in) :: caption
        integer, intent (in)          :: step
        character*4                   :: path, prefix

        write(path,'(I4.4)') step
        call system('mkdir -p '//trim(pth_wrk)
     $                         //trim(pth_dbg)//trim(path))
        write(prefix, '(I3.3,".")') dbg_seq_i
        dbg_seq_i = dbg_seq_i + 1
        open(42, file=trim(pth_wrk)//trim(pth_dbg)
     $                //path//"/"//trim(prefix)
     $                //trim(pfx_dbg)//"_"//caption
     $          ,form='unformatted')
        write(42) var
        close(42)

      end

      subroutine debug_write_3d(var, caption, step)

        implicit none
        include 'pom2k.c'

        real, intent (in)             :: var(im,jm,kb)
        character(len=*), intent (in) :: caption
        integer, intent (in)          :: step
        character*4                   :: path, prefix

        write(path,'(I4.4)') step
        call system('mkdir -p '//trim(pth_wrk)
     $                         //trim(pth_dbg)//trim(path))
        write(prefix, '(I3.3,".")') dbg_seq_i
        dbg_seq_i = dbg_seq_i + 1
        open(42, file=trim(pth_wrk)//trim(pth_dbg)
     $                //path//"/"//trim(prefix)
     $                //trim(pfx_dbg)//"_"//caption
     $          ,form='unformatted',position='rewind')
        write(42) var
        close(42)

      end

      subroutine upd_mnth

        include 'pom2k.c'

        real :: tind, b, e

        tind = mod(time,365.)
C
        if (tind.lt.15) then
          m = 1
          b = -15.5
          e = 15.
        else
          if (tind.lt.44) then
            m = 2
            b = 15.
            e = 44.
          else
            if (tind.lt.73.5) then
              m = 3
              b = 44.
              e = 73.5
            else
              if (tind.lt.104) then
                m = 4
                b = 73.5
                e = 104.
              else
                if (tind.lt.134.5) then
                  m = 5
                  b = 104.
                  e = 134.5
                else
                  if (tind.lt.165) then
                    m = 6
                    b = 134.5
                    e = 165.
                  else
                    if (tind.lt.195.5) then
                      m = 7
                      b = 165.
                      e = 195.5
                    else
                      if (tind.lt.226.5) then
                        m = 8
                        b = 195.5
                        e = 226.5
                      else
                        if (tind.lt.258) then
                          m = 9
                          b = 226.5
                          e = 258.
                        else
                          if (tind.lt.288.5) then
                            m = 10
                            b = 258.
                            e = 288.5
                          else
                            if (tind.lt.319) then
                              m = 11
                              b = 288.5
                              e = 319.
                            else
                              if (tind.lt.349.5) then
                                m = 12
                                b = 319.
                                e = 349.5
                              else
                                m = 1
                                b = 349.5
                                e = 380.5
                              end if
                            end if
                          end if
                        end if
                      end if
                    end if
                  end if
                end if
              end if
            end if
          end if
        end if
!
        fac = (tind-b)/(e-b)

        return

      end