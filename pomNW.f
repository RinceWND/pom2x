      program pom08
!                                                                      !
!----------------------------------------------------------------------!
!     For recent changes search for !lyo:_20080415:                    !
!                                                                      !
!     Main changes:                                                    !
!                                                                      !
!     (1). George Mellor found bugs in subr.profq (also in pom2k); use !
!          the new one in this version.                                !
!     (2). nsmolar=1 case (for wet-and-dry runs with nwad=1) is        !
!          incomplete, and should NOT be used                          !
!                                                                      !
!                       ... l.oey (Apr/15/2008)                        !
!                                                                      !
!----------------------------------------------------------------------!
!     This POM code has wetting and drying (WAD) capability            !
!                                                                      !
!     Set nwad=1 for WAD,set =0 to turn off WAD                        !
!     Set BOTH nwad=0 & nsmolar=0 to recover the original non-WAD POM  !
!                                                                      !
!     For details see:                                                 !
!                                                                      !
!     O2005 = Oey [2005; OceanModelling,  9, 133-150]                  !
!     O2006 = Oey [2006; OceanModelling, 13, 176-195]                  !
!     O2007 = Oey, Ezer et al [2007; Ocean Dynamics, 57, 205-221]      !
!                                                                      !
!     or, http://www.aos.princeton.edu/WWWPUBLIC/PROFS/                !
!     also, .../WWWPUBLIC/PROFS/pom98_wad_release_descriptions.html    !
!                                                                      !
!     For changes that I made, search for the following keyword:       !
!                                                                      !
      !lyo:!wad:  or simply "wad:"                                     !
      !tne:!wad:  (changes taken from T.Ezer's pom98_wad seamount case)!
!                                                                      !
!     Three (3) files are needed to run the model test cases:          !
!                                                                      !
!     pom2k.f (code), pomNW.c (common blocks etc) & pom2k.n (netcdf)   !
!                                                                      !
!     see the runscript runpom08* provided with these files.           !
!                                                                      !
!     "DOUBLE  PRECISION (i.e. REAL*8)" option can be used to          !
!     compile, with or without WAD; e.g. w/Intel ifort, use -r8:       !
!                                                                      !
!     ifort pom2k.f -o a.out -r8                                       !
!                                                                      !
!     the -r8 option is the same as specifying -double precision_size 64 or        !
!     -auto_double.  However, the code works with or without -r8       !
!                                                                      !
!     A link to the netCDF library also needs to be specified. For this!
!     and other details, see the runscript:                            !
!                                                                      !
!     runpom08*                                                        !
!                                                                      !
!     NOTE: for nsmolar=1, the PARAMETER (IM=???,JM=???) in subroutines!
!     wadadvt2d & wadsmoladif must match that specified in pomNW.c     !
!     This is done using "include 'grid'"                              !
!                                                                      !
!     Send bugs and/or improvements to: lyo@princeton.edu              !
!                                                                      !
!                       ... l.oey (May/02/2007)                        !
!                                 (Jul/24,30/2007)                     !
!                                 (Aug/12,18/2007)                     !
!                                                                      !
!     Support from the Minerals Management Service (MMS) for making    !
!     the development of WAD in POM possible is acknowledged.          !
!                                                                      !
!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
clyo:Notes:                                                            !
c     This version incorporates WAD into pom2k.f, and also corrects    !
c     the bug found by Charles Tang, Glen Carter and Alain Caya (and   !
c     corrected in the 2006-05-03 version of pom2k.f).                 !
c     Search for "clyo:"for all the changes I have made.               !
c                                                                      !
c                       ... lyo (Jun/14/2006)                          !
c                                                                      !
C **********************************************************************
C *                                                                    *
C *   The last code change as rcorded in pomNW.change was on           *
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
      use date_utility
      implicit none
C
!lyo:!wad:Add WAD variables in pomNW.c
      include 'pomNW.c'

      integer, external :: create_output
C
C     New declarations plus ispi,isp2i:
C
      double precision aam_init,atot
      double precision cbcmax,cbcmin,darea
      double precision days,dte2,dvol
      double precision eaver
      double precision horcon
      double precision ispi,isp2i
      double precision period,prtd1,prtd2,fsplt
      double precision saver,smoth,sw,swtch
      double precision taver,time0
      double precision vamax,vtot,tsalt
      double precision z0b
      double precision tatm,satm
      double precision wadsmoth  !lyo:!wad:
      double precision cflmin    !lyo:_20080415:
      integer fprint
      integer io(100),jo(100),ko(100)
      integer i,iend,iext,imax,ispadv,isplit,iswtch
      integer j,jmax
      integer k
      integer nadv,nbct,nbcs,nitera,nread
      integer iproblem
      integer nsmolar  !lyo:!wad:
      logical lramp
      character*120 netcdf_file

      double precision    bkp_gap
      integer ncid   !rwnd: Ncdf file id.
      integer nccnt  !    : Ncdf file counter.
      character*256 filename
      character*26  timestamp
!     target point coordinates
      integer tgt_lon, tgt_lat, tgt_sig
      double precision    slice_b, slice_e  
      
      integer :: date_start(3)
      
!     Formatting parameters
      character*4 :: sBOLD, sRESET
      character*5 :: sRED
      parameter( sBOLD = char(27)//"[1m",
     $           sRED  = char(27)//"[31m",
     $           sRESET= char(27)//"[0m" )
C
C***********************************************************************
C
C     source should agree with source_c in pomNW.c and source_n in
C     pom2k.n.
C
      source='pom08  2008-04-18'
C
      if(source.ne.source_c) then
        write(6,7)
    7   format(/'Incompatible versions of program and include files ',
     $          '..... program terminated; in pom2k.f'/)
        stop
      endif
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
C     NOTE that the array sizes im, jm and kb should be set in
C     pomNW.c or the 'grid' file created by runpom08
C
C-----------------------------------------------------------------------
C
      title='Run 1                                   ' ! run's title
C
C-----------------------------------------------------------------------
C
      netcdf_file='pom2k.nc'  ! netCDF output file
c     netcdf_file='nonetcdf'      ! disable netCDF output
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
C         4n       WAD prob#n, n=1,2 or 3     !lyo:!wad:
C
      iproblem=41
C
!                                                                      !
!----------------------------------------------------------------------!
!                                                                      !
!     WAD Control Parameters:                                          !
!                                                                      !
!lyo:!wad:set nwad=1 for WAD run, =0 to recover non-WAD run with min   !
!     depths set at say 10m (defined in routine called depending on    !
!                            value of iproblem)                        !
!                                                                      !
      nwad=0    !rwnd: PV=1
!                                                                      !
!lyo:!wad:set nsmolar=1 if water depth is to be solved by Smolarkiewicz!
!     See O2005 or O2006.                                              !
!     NOTE: nwad=1 w/nsmolar=1 may require finer grid AND even a larger!
!     hc, below, say hc=0.1 instead of default hc=0.05                 !
!                                                                      !
      nsmolar=0
!                                                                      !
!lyo:!wad:Define hhi0, later hhi will be set to it:                    !
!                                                                      !
!     hhi = HIghest water [above MSL] ever expected, this should       !
!           be the observed_highest + say 1 meter to make sure;        !
!           here, MSL = Mean Sea Level
!                                                                      !
!     Choose hhi0=20 (below) assuming that the highest tide or surge   !
!     is not expected to rise above 20m wrt MSL                        !
!                                                                      !
      hhi0=20.e0
!                                                                      !
!lyo:!wad:Define hc in meters:                                         !
!                                                                      !
!     hc  = Thinnest fluid layer (i.e. water depth) below which cell   !
!           is considered dry; default is 0.05m but can be changed in  !
!           "params" below; used only if nwad=1                        !
!                                                                      !
      hc=0.05e0
!                                                                      !
!lyo:!wad:Set wadsmoth.gt.0 to (cosmetically) smooth out isolated,     !
!     thin-fluid wet cells. Smoothing is done if the wet cell is thin, !
!     thickness <= hc*(1+wadsmoth); recommended non-zero wadsmoth is   !
!     ~0.05 but it should not be >~0.1                                 !
      wadsmoth=0.e0
!                                                                      !
!----------------------------------------------------------------------!
!                                                                      !
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
      nitera=3  !=2   !lyo:!wad:
C
C     Smoothing parameter. This should preferably be 1, but 0 < sw < 1
C     gives smoother solutions with less overshoot when nitera > 1:
C
      sw=1.e0  !=0.5e0   !lyo:!wad:
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
!                                                                      !
!----------------------------------------------------------------------!
!                                                                      !
!lyo:!wad:For 3-d runs w/rivers (strong stratification)                !
!     strong currents can (and physically they should) develop         !
!     in thin film of water over nearly dry cells, so dte and          !
!     isplit (below) may need to be smaller.  For dx~dy~1km,           !
!     I have used values as small as DTE=3.0 and ISPLIT=5.             !
!                                                                      !
!     Note#1: Since alpha=0 when nwad=1 (see below), the corresponding !
!             DTE needs to be smaller than the DTE used for nwad=0.    !
!     Note#2: ISPLIT should also be smaller than its value for no WAD; !
!             see Appendix A of Oey [2006; Ocean Modell. 13, 176-195]. !
!                                                                      !
!     For example, for the seamount problem that Tal Ezer has tested   !
!     (im,jm=65,49),  he found that DTE=60 & ISPLIT=30 worked for      !
!     nwad=0 but blew up for nwad=1. However, by decreasing ISPLIT     !
!     =5->10, the model with nwad=1 works for DTE=20->60, except that  !
!     isolated cells w/thin fluid layers (~5cm) appear, surrounded by  !
!     cells which are mostly dry; this probably is to be expected for  !
!     the relatively large dx&dy >2~4km used.  Using smaller DTE=15    !
!     & ISPLIT=5, or DTE=6 & ISPLIT=10 (as used below) can remove the  !
!     isolated thin-fluid cells.                                       !
!----------------------------------------------------------------------!
C
      dte=6.e0
C
C-----------------------------------------------------------------------
C
C     <Internal (3-D) time step>/<External (2-D) time step>
C     (dti/dte; dimensionless):
C
      isplit=10 !tne:!wad:dte=6.0 & isplit=30 work w/nwad=0 & hmax=200m
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
      time_end  ='not-set'
C
C-----------------------------------------------------------------------
C
      days=0.025d0       ! run duration in days
C
!-----------------------------------------------------------------------
!
      fsplt=1.           ! Output file splitting (days)
!
C-----------------------------------------------------------------------
C
      prtd1=0.0125d0     ! Initial print interval (days)
C
C-----------------------------------------------------------------------
C
      prtd2=1.d0         ! Final print interval (days)
C
C-----------------------------------------------------------------------
C
      swtch=1000.d0      ! Time to switch from prtd1 to prtd2
C
C-----------------------------------------------------------------------
C
      iskp=1             ! Printout skip interval in i
C
C-----------------------------------------------------------------------
C
      jskp=1             ! Printout skip interval in j
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
      zsh=.01e0          ! Bottom Log-layer shift (metres) !lyo:!wad:
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
!lyo:!wad:I would use horcon=0.1 (un-related to WAD)
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
C     equation (a value of alph0 = 0.e0 is perfectly acceptable, but the
C     value, alph0=.225e0 permits a longer time step):
C
      alph0=0.225e0 !lyo:!wad:use alph0 instead of alpha - defined later
C
C-----------------------------------------------------------------------
C
C     Initial value of aam:
C
      aam_init=500.e0
!                                                                      !
!----------------------------------------------------------------------!
!lyo:!wad:Tidal amplitude for iproblem=41 w/ or w/o WAD:               !
!                                                                      !
!lyo:_20080415:
      tidamp=0.e0
      if (iproblem .eq. 41) tidamp=9.e0
!----------------------------------------------------------------------!
!                                                                      !
!     initialise time0 (May be overriden with params)
      time0 = 0.
!     initialise tgt coordinates
      tgt_lon = 1
      tgt_lat = 1
      tgt_sig = 1
!     Clock init
      slice_b = 0.
      slice_e = 0.
!     Climatological cycle
      m0        = 1
      mi        = 1
      clm_cycle = 12

      nccnt = 0
!
!     Initialise bry read flags
!
      rf_clm   = 0  ! note that this one is climate and not bry.
      rf_rmn   = 0  ! this one is rmean
      rf_el    = 0
      rf_uv    = 0
      rf_ts    = 0
      rf_sts   = 0
      rf_swrad = 0
      rf_wtsur = 0
      rf_wsurf = 0
!
!     Initialise BC flags
!
      BC%ipl = .true.
      BC%wnd = .true.
      BC%lrd = .true.   ! TODO: implement lrad
      BC%srd = .true.
      BC%ssf = .false.
      BC%vap = .false.
      BC%clm = .false.
      BC%bnd%nth = .false.
      BC%bnd%est = .false.
      BC%bnd%sth = .false.
      BC%bnd%wst = .false.

      IC%el = .false.
      IC%u  = .false.
C
C     End of input of constants
C***********************************************************************
C
C --- Above are the default parameters, alternatively one can
C --- use parameters from a file created by runscript runpom2k_pow_wad !clyo:wad:
C
      include 'params'
C
clyo:wad:beg:
c     Overwrite some input constants: see "params" above in runpom08
clyo:wad:end:
c
C***********************************************************************
C
C     Calculate some constants:
C
      dti=dte*float(isplit)
      dte2=dte*2
      dti2=dti*2
!                                                                      !
!----------------------------------------------------------------------!
!lyo:!wad:                                                             !
!                                                                      !
!     Dry-cell T/S relaxation using leapfrog-trapezoidal:              !
      dtrat=dti/(0.5*86400.)  !1/2 day relaxation time-scale
      cwetrlx1=(1.-dtrat)/(1.+dtrat); cwetrlx2=2.*dtrat/(1.+dtrat)
!                                                                      !
      alpha=alph0*float(1-nwad) !lyo:!wad: must=0 for WAD
!                                                                      !
      hhi=hhi0*float(nwad)
!                                                                      !
      nsmolar=nsmolar*nwad
!                                                                      !
!     Check nsmolar:                  !lyo:_20080415:                  !
!                                                                      !
      if( (nsmolar.ne.0).and.(mode.ne.2) )then
      write(6,'('' Stopped; incorrect inputs of '')')
      write(6,'('' nsmolar   = '',i10)') nsmolar
      write(6,'('' mode      = '',i10)') mode
      write(6,'('' For nsmolar ne 0, the present version works  '')')
      write(6,'('' only if mode = 2. So either run in barotropic '')')
      write(6,'('' mode = 2, or set nsmolar=0 and resubmit '')')
      stop
      endif
!                                                                      !
!     Check iproblem & nwad compatibility:                             !
!                                                                      !
      if ( iproblem<11 .or. iproblem>19 ) then
        if( (iproblem.ne.41).and.(nwad.eq.1) )then !lyo:_20080415:
        write(6,'('' iproblem   = '',i10)') iproblem
        write(6,'('' nwad       = '',i10)') nwad
        write(6,'('' iproblem must = 41 for nwad = 1 '')') !lyo:_20080415:
        write(6,'('' Stopped: iproblem & nwad are not compatible '')')
        stop
        endif
      end if
!                                                                      !
!----------------------------------------------------------------------!
!                                                                      !
!
      if (time_end/='not-set') then         ! rwnd: If end date is set, override days variable.
        days=Days_in_Between(time_start, time_end)
      end if
      iend=max0(nint(days*24.e0*3600.e0/dti),2)
      iprint=nint(prtd1*24.e0*3600.e0/dti)
      iswtch=nint(swtch*24.e0*3600.e0/dti)
      fprint=nint(fsplt*24.e0*3600.e0/dti)  ! rwnd:
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
      write(6,'('' nwad       = '',i10)') nwad        !lyo:!wad:
      write(6,'('' nsmolar    = '',i10)') nsmolar     !lyo:!wad:
      write(6,'('' hhi        = '',f10.3)') hhi       !lyo:!wad:
      write(6,'('' hc         = '',f10.3)') hc        !lyo:!wad:
      write(6,'('' wadsmoth   = '',f10.3)') wadsmoth  !lyo:!wad:
      write(6,'('' mode       = '',i10)') mode
      write(6,'('' nadv       = '',i10)') nadv
      write(6,'('' nitera     = '',i10)') nitera
      write(6,'('' sw         = '',f10.4)') sw
      write(6,'('' nread      = '',i10)') nread
      write(6,'('' dte        = '',f10.2)') dte
      write(6,'('' dti        = '',f10.1)') dti
      write(6,'('' isplit     = '',i10)') isplit
      write(6,'('' time_start = '',a26)') time_start
      write(6,'('' time_end   = '',a26)') time_end
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
      write(6,'('' zsh        = '',f10.6)') zsh         !lyo:!wad:
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
      write(6,'('' tidamp     = '',f10.4)') tidamp
C
C-----------------------------------------------------------------------
C
C     Initialise boundary arrays:
C
      do i=1,im
        vabn(i)=0.e0
        vabs(i)=0.e0
        eln(i)=0.e0     !lyo:!wad:note:keep=0.e0 inst. of -hhi
        els(i)=0.e0     ! .. redefined later in wadseamount etc.
        do k=1,kb
          vbn(i,k)=0.e0
          vbs(i,k)=0.e0
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
        ele(j)=0.e0     !lyo:!wad:note:keep=0.e0 inst. of -hhi
        elw(j)=0.e0     ! .. redefined later in wadseamount etc.
        do k=1,kb
          ube(j,k)=0.e0
          ubw(j,k)=0.e0
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
          elb(i,j)=0.e0     !lyo:!wad:note:keep=0.e0 inst. of -hhi
          etb(i,j)=0.e0     ! .. redefined later in wadseamount etc.
          e_atmos(i,j)=0.e0 !lyo:!wad:note:=0.e0 is correct
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
      else if(iproblem.eq.11) then
        call ncdf2ic
      else if(iproblem.eq.12) then  !rwnd:unstratified case
        call ncdf2ic_unstrat
      else if(iproblem.eq.13) then  !rwnd:box case to test curvilineal grid
        call ncdf2ic_box
      else if(iproblem.eq.14) then  !rwnd:pom grid generated ics
        call ncdf2ic_pom
      else if(iproblem.eq.41) then  !lyo:!wad:
        call wadseamount
      else
        write(6,8)
    8   format(/' Invalid value of iproblem ..... program terminated'/)
        stop
      endif
!
!lyo:_20080415:
!     Make sure that wetmask = fsm for non-WAD runs:
      if (nwad.eq.0) then
         wetmask(:,:)=fsm(:,:)
         endif
!
C
C     Inertial period for temporal filter:
C
!lyo:_20080415:
      if (cor(im/2,jm/2).ne.0.0) then
         period=(2.e0*pi)/abs(cor(im/2,jm/2))/86400.e0
      else
         period=1.0  !set to 1day
      endif
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
          el(i,j)=elb(i,j)  !lyo:!wad:note:already defined
          et(i,j)=etb(i,j)  !          in wadseamount w/hhi etc
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
      if (horcon.le.0.0) aam(:,:,:)=-horcon  !lyo:_20080415:
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
      call baropg
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
!lyo:!wad:cbc_change:begins:-------------------------------------------!
!                                                                      !
!----------------------------------------------------------------------!
!     Calculate bottom friction coefficient:                           !
!     ... see O2006, Appendix A; also Mellor,JPO,2002:Osc.Blayer...    !
!----------------------------------------------------------------------!
!                                                                      !
      do j=1,jm
        do i=1,im
          cbc(i,j)=cbcmin
!         cbc(i,j)=(kappa/log((1.e0+zz(kbm1))*h(i,j)/z0b))**2
          cbc(i,j)=(kappa/log((zsh+(1.e0+zz(kbm1))*d(i,j))/z0b))**2
!lyo:!wad:lyo's originals:---------------------------------------------!
!vers.1:  cbc(i,j)=(kappa/log((zsh+(zz(kbm1)-z(kb))*d(i,j))/z0b))**2   !
!vers.2:  cbc(i,j)=max(cbcmin,                                         !
!    $        (kappa/log((zsh+(zz(kbm1)-z(kb))*d(i,j))/z0b))**2)       !
!lyo:!wad:lyo's originals:---------------------------------------------!
!
          cbc(i,j)=max(cbcmin,cbc(i,j))
C
C     If the following is invoked, then it is probable that the wrong
C     choice of z0b or vertical spacing has been made:
C
          cbc(i,j)=min(cbcmax,cbc(i,j))
        end do
      end do
!lyo:!wad:cbc_change:ends:---------------------------------------------!
C
C     Calculate external (2-D) CFL time step:
C
!lyo:_20080415:
      cflmin=1.e10
      do j=1,jm
        do i=1,im
          tps(i,j)=0.5e0/sqrt(1.e0/dx(i,j)**2+1.e0/dy(i,j)**2)
     $               /sqrt(grav*(h(i,j)+small))*fsm(i,j)
          if (fsm(i,j).ne.0.0) cflmin=min(cflmin,tps(i,j)) !lyo:_20080415:
        end do
      end do
      
      write(*,*) "Minimal tps: ",minval(tps, tps>0.)
      write(*,*) "Difference between minimal tps and chosen dte:",
     $           minval(tps)-dte
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
     $          ,wetmask,wmarsh   !lyo:!wad:
        close(70)
        
        if (time_end/='not-set') then   ! rwnd: In case of seamless reastart with specific end date, recalculate days left.
          days = Days_in_between(time_start, time_end)-time0
          if (days<=0) then
            write(*,*) sBOLD,sRED,"[X]",sRESET,
     $        " Simulation had already exceeded specified end date."
            write(*,*) "    Please, specify time_end variable ",
     $                 "past the ",
     $        DateTime_Build(Days_since_to_DateTime(time_start, time0))
            write(*,*) sBOLD,sRED,
     $                 "[[ SIMULATION TERMINATED ]]",sRESET
!            write(*,*) char(27),"[37;44m","White on blue background...",sRESET
            stop
          end if
        end if
        
      end if
!
      do j=1,jm
        do i=1,im
          d(i,j)=h(i,j)+el(i,j)
          dt(i,j)=h(i,j)+et(i,j)
        end do
      end do
!                                                                      !
!----------------------------------------------------------------------!
!lyo:!wad:Update cbc:                                                  !
!                                                                      !
      cbc(:,:)=cbcmin      !lyo:_20080415:
      if (mode.ne.2) then  !lyo:_20080415:
      do j=1,jm
        do i=1,im
!         cbc(i,j)=cbcmin  !lyo:_20080415:
          cbc(i,j)=(kappa/log((zsh+(1.e0+zz(kbm1))*d(i,j))/z0b))**2
          cbc(i,j)=max(cbcmin,cbc(i,j))
          cbc(i,j)=min(cbcmax,cbc(i,j))
        end do
      end do
      endif
!                                                                      !
!----------------------------------------------------------------------!
!
!   TODO: Updating month here leads to incorrect m0 value if restarting a file with climate warping
      call upd_mnth(time0, BC%ipl)
!   Store month offset for further warping.
      m0 = mi
!   Prevent writing zero time when it is not meant to be zero.
      time = time0
!
!-----------------------------------------------------------------------
!
      nccnt = int(iint/fprint)
      ncid = create_output(nccnt)  ! rwnd:
      call ncflush(ncid
     $            ,int(modulo(float(iint)/float(iprint)
     $                       ,float(fprint)/float(iprint)))+1)
!
      write(*,*) "Creating target point output..."
      filename = trim(pth_wrk)//trim(pth_out)//
     $             trim(title)//"_tgt.csv"
      open(49, file=filename)
      call time2date(time, time_start, timestamp)
      write(49,*) timestamp, ";",
     $ u(tgt_lon,tgt_lat,tgt_sig), ";", v(tgt_lon,tgt_lat,tgt_sig)
C
C-----------------------------------------------------------------------
C
C     Initial conditions:
C
C     Select print statements in printall as desired:
C
!      call printall
C
C-----------------------------------------------------------------------
C
C     Initialise netCDF output and output initial set of data:
C
!        if(netcdf_file.ne.'nonetcdf') then
!      call write_netcdf(netcdf_file,1)                        ! *netCDF*
!      call write_netcdf(netcdf_file,2)                        ! *netCDF*
!        endif
C
C-----------------------------------------------------------------------
!     Set timer
      call cpu_time(slice_b)
C
      do 9000 iint=1,iend      !  Begin internal (3-D) mode
C
        time=dti*float(iint)/86400.e0+time0
C
        if(lramp) then
          ramp=time/period
          if(ramp.gt.1.e0) ramp=1.e0
        else
          ramp=1.e0
        endif
C
C       write(6,2) mode,iint,time
C   2   format(' mode,iint,time =',2i5,f9.2)
C
C-----------------------------------------------------------------------
C
C     Set time dependent, surface and lateral boundary conditions.
C     The latter will be used in subroutine bcond. Users may
C     wish to create a subroutine to supply wusurf, wvsurf, wtsurf,
C     wssurf, swrad and vflux.
!
!     Update month
      call upd_mnth(time, BC%ipl)
      call clm_warp

      if (iproblem>=11 .and. iproblem<=19) then
!
      else
C     Introduce simple wind stress. Value is negative for westerly or
C     southerly winds. The following wind stress has been tapered
C     along the boundary to suppress numerically induced oscilations
C     near the boundary (Jamart and Ozer, J.G.R., 91, 10621-10631).
C     To make a healthy surface Ekman layer, it would be well to set
C     kl1=9.
C
        do j=2,jmm1
          do i=2,imm1
c
      if(iproblem.ne.3) then     ! constant wind read in file2ic
c
c           wusurf(i,j)=ramp*(1.e-4*cos(pi*(j-1)/jmm1))
!           wusurf(i,j)=1.00*(1.e-4*cos(pi*(j-1)/jmm1))
!    $                    *.25e0*(dvm(i,j+1)+dvm(i-1,j+1)
!    $                          +dvm(i-1,j)+dvm(i,j))
C --- no wind ----
            wusurf(i,j)=0.e0  !lyo:_20080415:
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
C     salinity will change and, undouble precisionistically, so will total salt.
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
      
      if (BC%clm) then
        call bry(0)   ! Update Rmean
        call bry(1)   ! Update TSclim
      end if
      if (BC%wnd) call flux(5)  ! Update wind stress
      if (BC%ele) call bry(2)   ! Update boundary elevation
      if (BC%clm) call bry(4)   ! Update boundary TS
      if (BC%vel) call bry(3)   ! Update current velocity BCs.
C
clyo:
!     call powdriver(iprint,nread,z0b,cbcmin,iend/iprint,fsm)
c
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
          call baropg
C
          if (horcon.gt.0.0) then !lyo:_20080415:
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
          endif  !if (horcon.gt.0.0) then !lyo:_20080415:
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
clyo:beg:changes by Tang et al.
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
clyo:end:
C
C-----------------------------------------------------------------------
C
clyomoving:decide to skip wriv & wmarsh for now; wriv in particular
c     needs to make consistent with vflux


        do 8000 iext=1,isplit    ! Begin external (2-D) mode
C
C         write(6,3) iext,time
C   3     format(' iext,time =',i5,f9.2)
C
!                                                                      !
!----------------------------------------------------------------------!
!lyo:!wad:Water depth can be solved by Smolarkiewicz to prevent oscil. !
!                                                                      !
!     nsmolar=1 case (for wet-and-dry runs with nwad=1) is             !
!     incomplete, and should NOT be used                               !
!                                                                      !
      if (nsmolar.eq.1 .and. mode.eq.2) then  !lyo:_20080415:
      dd0(:,:)=h(:,:)+ el(:,:)*fsm(:,:)
      ddb(:,:)=h(:,:)+elb(:,:)*fsm(:,:)
      tps(:,:)=0.0           !aam2d(:,:)
      call wadadvt2d(ddb,dd0,ddf,nitera,sw,dx,dy,ua,va,art,aru,arv,
     1    fsm,dum,dvm,tps,dte2)
      elf(:,:)=ddf(:,:)-h(:,:)
      else
!
          do j=2,jm
            do i=2,im
              fluxua(i,j)=.25e0*(d(i,j)+d(i-1,j))
     $                     *(dy(i,j)+dy(i-1,j))*ua(i,j)
              fluxva(i,j)=.25e0*(d(i,j)+d(i,j-1))
     $                     *(dx(i,j)+dx(i,j-1))*va(i,j)
            end do
          end do
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
!                                                                      !
      endif
!                                                                      !
!----------------------------------------------------------------------!
C
          call bcond(1)
!                                                                      !
!----------------------------------------------------------------------!
!lyo:!wad:Check wet/dry on ELF:                                        !
!                                                                      !
      if (nwad.eq.1) then
      do j=1,jm; do i=1,im
      DWET(I,J)=H(I,J)+ELF(I,J)*fsm(i,j)
      enddo; enddo
      do j=1,jm; do i=1,im
      wetmask(i,j)=fsm(i,j)
      if (dwet(i,j).le.hc) wetmask(i,j)=0.0
      enddo; enddo
      if (wadsmoth.gt.0.0) then
!     Smooth isolated wet spots:                                       !
      do j=1,jm; do i=1,im
      if( (wetmask(i,j) .gt.
     $    (wetmask(i-1,j)+wetmask(i+1,j)+wetmask(i,j-1)+wetmask(i,j+1)))
     $    .and. (dwet(i,j).le.hc*(1.e0+wadsmoth)) ) then
!    $    .and. (dwet(i,j).le.hc*1.01) ) then
         wetmask(i,j)=0.0
!        elf(i,j)=elb(i,j); dwet(i,j)=h(i,j)+elf(i,j)*fsm(i,j)
         endif
      enddo; enddo
      endif
      endif
!                                                                      !
!----------------------------------------------------------------------!
!                                                                      !
          if(mod(iext,ispadv).eq.0) call advave(tps)
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
          call bcond(2)
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
          endif
!                                                                      !
!----------------------------------------------------------------------!
!lyo:!wad:Check wet/dry on UAF & VAF:                                  !
!                                                                      !
      if (nwad.eq.1) then
      do j=1,jm; do i=2,im
      if (0.5*(dwet(i,j)+dwet(i-1,j)).le.hc) then
      uaf(i,j)=0.0
      else                      !Set flux=0 if coming from dry cell:
      if (wetmask(i-1,j).eq.0.0 .and. uaf(i,j).gt.0.0) uaf(i,j)=0.0
      if (wetmask(i,j)  .eq.0.0 .and. uaf(i,j).lt.0.0) uaf(i,j)=0.0
      endif
      enddo; enddo
!
      do j=2,jm; do i=1,im
      if (0.5*(dwet(i,j)+dwet(i,j-1)).le.hc) then
      vaf(i,j)=0.0
      else                      !Set flux=0 if coming from dry cell:
      if (wetmask(i,j-1).eq.0.0 .and. vaf(i,j).gt.0.0) vaf(i,j)=0.0
      if (wetmask(i,j)  .eq.0.0 .and. vaf(i,j).lt.0.0) vaf(i,j)=0.0
      endif
      enddo; enddo
      endif
!                                                                      !
!----------------------------------------------------------------------!
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
!                                                                      !
!----------------------------------------------------------------------!
!lyo:!wad:Update cbc:                                                  !
!                                                                      !
      if (mode.ne.2) then  !lyo:_20080415:
      do j=1,jm
        do i=1,im
          cbc(i,j)=cbcmin
          cbc(i,j)=(kappa/log((zsh+(1.e0+zz(kbm1))*d(i,j))/z0b))**2
          cbc(i,j)=max(cbcmin,cbc(i,j))
          cbc(i,j)=min(cbcmax,cbc(i,j))
        end do
      end do
      endif
!                                                                      !
!----------------------------------------------------------------------!
!
            if(iext.ne.isplit) then
clyo:beg:changes by Tang et al.
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
clyo:end:
            endif
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
C
C     vertvl calculates w from u, v, dt (h+et), etf and etb:
C
            call vertvl(a,c)
            call bcond(5)
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
C
C     Calculate tf and sf using uf, vf, a and c as temporary variables:
C
            if(mode.ne.4) then
!                                                                      !
!----------------------------------------------------------------------!
!lyo:!wad:The followings prepare Oey's [1996] and Oey & Chen's [1992]  !
!     implementations of river and temperature surface conditions      !
!     i.e. using wriv -- but not implemented yet; they are not         !
!     directly for WAD, but ssurfwet is required later for wad         !
!                                                                      !
            do j=1,jm
               do i=1,im
                 ssurfwet(i,j)=SSURF(i,j)*(1.-WRIV(i,j)/( WRIV(i,j)+
     $                                                     1.E-28)  )
               enddo
            enddo
!                                                                      !
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
C
              call proft(uf,wtsurf,tsurf,nbct,tps)
              call proft(vf,wssurf,ssurf,nbcs,tps)
              call bcond(4)
!                                                                      !
!----------------------------------------------------------------------!
!lyo:!wad:Prepare dry-cell T/S "initial conditions" for next time-step !
!                                                                      !
           if (nwad.eq.1) then
              do k=1,kb-1; do j=1,jm; do i=1,im
                 if (wetmask(i,j).eq.0.0) then
                    uf(i,j,k)=(cwetrlx1*tb(i,j,k)+cwetrlx2*tsurf(i,j)
     $                        *fsm(i,j))
                    vf(i,j,k)=(cwetrlx1*sb(i,j,k)+cwetrlx2*ssurfwet(i,j)
     $                        *fsm(i,j))
                 endif
                 enddo; enddo; enddo
           endif
!----------------------------------------------------------------------!
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
C
              call dens(s,t,rho)
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
!                                                                      !
!----------------------------------------------------------------------!
!lyo:!wad:Check wet/dry on UF & VF:                                    !
!                                                                      !
           if (nwad.eq.1) then
              do k=1,kb-1
                 do j=1,jm; do i=2,im
                    if (wetmask(i,j)*wetmask(i-1,j).eq.0.0)uf(i,j,k)=0.0
                    enddo; enddo
                 do j=2,jm; do i=1,im
                    if (wetmask(i,j)*wetmask(i,j-1).eq.0.0)vf(i,j,k)=0.0
                    enddo; enddo
              enddo
           endif
!----------------------------------------------------------------------!
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
C
        endif
C
C-----------------------------------------------------------------------
C
C     Beginning of print section:
C
        if(iint.ge.iswtch) iprint=nint(prtd2*24.e0*3600.e0/dti)
C
!     Target point output
        if(mod(iint,3).eq.0.or.vamax.gt.vmaxl) then    ! Print it every three internal steps
!          call ncTgtFlush(49)
          call time2date(time, time_start, timestamp)
          write(49,*) timestamp, ";",
     $ u(tgt_lon,tgt_lat,tgt_sig), ";", v(tgt_lon,tgt_lat,tgt_sig)
        end if
!
        if(mod(iint,iprint).eq.0.or.vamax.gt.vmaxl) then
C
!     Stop timer
          call cpu_time(slice_e)
          write(6,4) time,iint,iext,iprint,(slice_e-slice_b)/60.
    4     format(/
     $    '**************************************************',
     $    '**************************************************',
     $    '*************************'//
     $    ' time =',f9.4,', iint =',i8,', iext =',i8,', iprint =',i8,//
     $    ' Time elapsed since prev. output: ',f9.4,' min')
C
C     Select print statements in printall as desired:
C
!          call printall
C
!          call wadout  !lyo:!wad:
          
          if (int(iint/fprint)/=nccnt) then
              nccnt = int(iint/fprint)
              write(*,*) "Output file change (", nccnt, ")"
              call ncclose(ncid)
              write(*,*) "Previous file closed."
              ncid = create_output(nccnt)
          end if
          call ncflush(ncid
     $                ,int(modulo(float(iint)/float(iprint)
     $                           ,float(fprint)/float(iprint)))+1)
!
          vtot=0.e0
          atot=0.e0
          taver=0.e0
          saver=0.e0
          eaver=0.e0
          do k=1,kbm1
            do j=1,jm
              do i=1,im
                darea=dx(i,j)*dy(i,j)*wetmask(i,j) !lyo:!wad:
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
              darea=dx(i,j)*dy(i,j)*wetmask(i,j) !lyo:!wad:
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
!            if(netcdf_file.ne.'nonetcdf') then
!          call write_netcdf(netcdf_file,2)                    ! *netCDF*
!            endif
C
          if(vamax.gt.vmaxl) then
C
            write(6,4) time,iint,iext,iprint
C
!            call printall
            call ncflush(ncid
     $                  ,int(modulo(float(iint)/float(iprint)
     $                            ,float(fprint)/float(iprint)))+2) ! Increment by one in order to not to overwrite the last output.
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
!              if(netcdf_file.ne.'nonetcdf') then
!            call write_netcdf(netcdf_file,3)                  ! *netCDF*
!              endif
C
            stop
C
          endif
!     Reset timer
          call cpu_time(slice_b)
!
        endif
!
        if (mod(time,bkp_gap).eq.0) then
          if (iint.ne.iend) then
            write(*,*) "[*] Backing up..."
            open(71,file=trim(pth_wrk)//trim(pth_bkp)
     $                 //trim(title)//'.restart.bkp',
     $              form='unformatted',position='rewind')
            write(71) time,
     $        wubot,wvbot,aam2d,ua,uab,va,vab,el,elb,et,etb,egb,
     $        utb,vtb,u,ub,w,v,vb,t,tb,s,sb,rho,adx2d,ady2d,advua,advva,
     $        km,kh,kq,l,q2,q2b,aam,q2l,q2lb,
     $        wetmask,wmarsh    ! WAD
            close(71)
          end if
        end if
C
C     End of print section
C
C-----------------------------------------------------------------------
C
 9000 continue       !  End of internal (3-D) mode
C
C-----------------------------------------------------------------------
C
!   Finish target point output
!      call ncTgtFlush(49)
      call ncclose(ncid)
      call time2date(time, time_start, timestamp)
      write(49,*) timestamp, ";",
     $ t(tgt_lon,tgt_lat,tgt_sig), ";", s(tgt_lon,tgt_lat,tgt_sig)
      close(49)
!
      write(6,4) time,iint,iext,iprint
!
!            call printall !lyo:_20080415:final printing
!
C
C     Set levels for output:
C
!      ko(1)=1
!      ko(2)=2
!      ko(3)=kb/2
!      ko(4)=kb-1
!      ko(5)=kb
!C
!      call prxyz('Vertical velocity, w                    ',
!     $           time,w       ,im,iskp,jm,jskp,kb,ko,5,-1.e0)
C
C     call prxyz('Turbulent kinetic energy x 2, q2        ',
C    $           time,q2      ,im,iskp,jm,jskp,kb,ko,5,-1.e0)
C
C     Save this data for a seamless restart:
C
clyo:
!     call powsave
c
      write(*,*) "[ ] Writing restart file..."
      open(71,file=trim(pth_wrk)//trim(pth_bkp)
     $           //trim(title)//'.restart.bkp',
     $        form='unformatted',position='rewind')
      write(71) time,
     $  wubot,wvbot,aam2d,ua,uab,va,vab,el,elb,et,etb,egb,
     $  utb,vtb,u,ub,w,v,vb,t,tb,s,sb,rho,adx2d,ady2d,advua,advva,
     $  km,kh,kq,l,q2,q2b,aam,q2l,q2lb
     $ ,wetmask,wmarsh   !lyo:!wad:
      close(71)
C
C     Close netCDF file:
C
!        if(netcdf_file.ne.'nonetcdf') then
!      call write_netcdf(netcdf_file,3)                        ! *netCDF*
!        endif
C
!lyo:!wad:print final successful completion:
!     Execution time total:
      call cpu_time(slice_e)
      write(6,10) time
   10 format(/2x,'JOB SUCCESSFULLY COMPLT.; time = ',1P1e13.5,' days')
!
      call elapsed_time_print(slice_e)
!
      stop
C
      end
C
C     End of main program
C
C-----------------------------------------------------------------------
C
      subroutine elapsed_time_print(raw)

        double precision, intent(in) :: raw
        double precision             :: sec
        integer          :: tmp, min, h, d

        sec = mod(raw, 60.)
        tmp = int(raw-sec)
        min = mod(tmp/60., 60.)
        tmp = (tmp-min*60)/60
        h   = mod(tmp/60., 24.)
        tmp = (tmp-h*60)/60
        d   = mod(tmp/24., 24.)

        write(*,"('Time elapsed: ',i4,'d ',i2,'h ',i2,'min ',f9.4,'s')")
     $            d, h, min, sec
        write(*,*) "Debug raw: ", raw

      end subroutine
!
      subroutine advave(curv2d)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Calculates horizontal advection and diffusion.      *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pomNW.c'
C
      double precision curv2d(im,jm)
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
      include 'pomNW.c'
C
      double precision xflux(im,jm,kb),yflux(im,jm,kb)
      double precision curv(im,jm,kb)
      double precision dtaam
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
      include 'pomNW.c'
C
      double precision qb(im,jm,kb),q(im,jm,kb),qf(im,jm,kb)
      double precision xflux(im,jm,kb),yflux(im,jm,kb)
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
      include 'pomNW.c'
C
      double precision
     $             fb(im,jm,kb),f(im,jm,kb),fclim(im,jm,kb),ff(im,jm,kb)
      double precision xflux(im,jm,kb),yflux(im,jm,kb)
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
      include 'pomNW.c'
C
      double precision
     $             fb(im,jm,kb),f(im,jm,kb),fclim(im,jm,kb),ff(im,jm,kb)
      double precision xflux(im,jm,kb),yflux(im,jm,kb)
      double precision sw
      integer nitera
      double precision fbmem(im,jm,kb),eta(im,jm)
      double precision
     $          xmassflux(im,jm,kb),ymassflux(im,jm,kb),zwflux(im,jm,kb)
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
      include 'pomNW.c'
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
      include 'pomNW.c'
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
      include 'pomNW.c'
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
C     Initialise and set up free surface mask if no WAD: !lyo:!wad:
C
      if (nwad.eq.0) then !lyo:!wad:
!
!     do j=1,jm
!       do i=1,im
!         fsm(i,j)=0.e0
!         dum(i,j)=0.e0
!         dvm(i,j)=0.e0
!         if(h(i,j).gt.1.e0) fsm(i,j)=1.e0
!       end do
!     end do
C
C     Set up velocity masks:
C
!     do j=2,jm
!       do i=2,im
!         dum(i,j)=fsm(i,j)*fsm(i-1,j)
!         dvm(i,j)=fsm(i,j)*fsm(i,j-1)
!       end do
!     end do
!
!lyo:_20080415:the above original pom2k version does NOT define dum & dvm
!     along i=j=1; use the folliwngs (from subroutine wadh).
      do j=1,jm; do i=1,im
      FSM(I,J)=1.; DUM(I,J)=1.; DVM(I,J)=1.
      IF(h(I,J).LE.1.0)THEN
      h(I,J)=1.0; FSM(I,J)=0.; DUM(I,J)=0.; DVM(I,J)=0.
      ENDIF
      enddo; enddo
      DO J=1,JM-1; DO I=1,IM
      IF(FSM(I,J).EQ.0..AND.FSM(I,J+1).NE.0.)DVM(I,J+1)=0.
      enddo; enddo
      DO J=1,JM;   DO I=1,IM-1
      IF(FSM(I,J).EQ.0..AND.FSM(I+1,J).NE.0.)DUM(I+1,J)=0.
      enddo; enddo
!
      endif   !if (nwad.eq.0) then       !lyo:!wad:
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
      include 'pomNW.c'
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
      subroutine bcond(idx)
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
C *                The diagram below and some of the text was provided *
C *                by D.-S. Ko. It is for the case where u and v are   *
C *                the primary boundary conditions together with t and *
C *                s (co-located with e). e = el itself is rather      *
C *                unimportant and is substituted from an adjacent     *
C *                interior point.                                     *
C *                                                                    *
C *                The dotted line (.......) bounds the interior       *
C *                (non-boundary) grid points. In general, only those  *
C *                variables in the interior are computed and          *
C *                variables at open boundary have to be specified.    *
C *                All interpolations are centered in space except     *
C *                those at lateral open boundary where an upstream    *
C *                scheme is usually used. "BU" indictates a line of   *
C *                points where the boundary conditions should be      *
C *                specified.                                          *
C *                                                                    *
C *                Horizontal locations of e(el), t and s (etc.) are   *
C *                coincident. Unused points are indicated by "NU" and *
C *                "*". However, for attractive output, adjacent       *
C *                interior values may be filled in at these points.   *
C *                                                                    *
C *                People not acquainted with sigma coordinates have   *
C *                often asked what kind of boundary condition is      *
C *                applied along closed horizontal boundaries.         *
C *                Although the issue is not as important as it might  *
C *                be  for z-level grids, a direct answer is "half-    *
C *                slip" which, of course, is between free slip and    *
C *                non-slip.                                           *
C *                                                                    **********
C-------------------------------------------------------------------------------*
C     |                               N O R T H                                 *
C     |                                                                         *
C     |    1     2     3           i-1   i    i+1         im-2  im-1   im       *
C-----+-------------------------------------------------------------------------*
C     |   NU BC BC                                                    BC BC     *
C     |    v  v  v                                                     v  v     *
C     |                                                                         *
C     |BC> u* e  u  e  u  e  .  .  u  e  u  e  u  e  .  .  u  e  u  e  u  e  <BC*
C     |    |     |     |           |     |     |           |     |     |        *
C  jm |BC> +--v--+--v--+--v--   .  +--v--+--v--+--v--   .  +--v--+--v--+--v- <BC*
C     |    |     | ....|...........|.....|.....|...........|.....|.... |        *
C     |    u* e  u :e  u  e  .  .  u  e  u  e  u  e  .  .  u  e  u  e: u  e     *
C     |    |     | :   |           |     |     |           |     |   : |        *
C jm-1|    +--v--+--v--+--v--   .  +--v--+--v--+--v--   .  +--v--+--v--+--v-    *
C     |    |     | :   |           |     |     |           |     |   : |        *
C     |    u* e  u :e  u  e  .  .  u  e  u  e  u  e  .  .  u  e  u  e: u  e     *
C     |    |     | :   |           |     |     |           |     |   : |        *
C jm-2|    +--v--+--v--+--v--   .  +--v--+--v--+--v--   .  +--v--+--v--+--v-    *
C     |            :                                                 :          *
C W   |       .  . :.  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .: .  .    E*
C E   |            :                    INTERIOR                     :         A*
C S   |       .    :.     .     .     .     .     .     .     .     .:    .    S*
C T   |            :                                                 :         T*
C     |    |     | :   |           |     |     |           |     |   : |        *
C     |    u* e  u :e  u  e  .  .  u  e  u  e  u  e  .  .  u  e  u  e: u  e     *
C     |    |     | :   |           |     |     |           |     |   : |        *
C   3 |    +--v--+--v--+--v--   .  +--v--+--v--+--v--   .  +--v--+--v--+--v-    *
C     |    |     | :   |           |     |     |           |     |   : |        *
C     |    u* e  u :e  u  e  .  .  u  e  u  e  u  e  .  .  u  e  u  e: u  e     *
C     |    |     | ....|...........|.....|.....|...........|.....|.... |        *
C   2 |BC> +--v--+--v--+--v--   .  +--v--+--v--+--v--   .  +--v--+--v--+--v- <BC*
C     |    |     |     |           |     |     |           |     |     |        *
C     |BC> u* e  u  e  u  e  .  .  u  e  u  e  u  e  .  .  u  e  u  e  u  e  <BC*
C     |    |     |     |           |     |     |           |     |     |      --*
C   1 |NU> +--v*-+--v*-+--v*-   .  +--v*-+--v*-+--v*-   .  +--v*-+--v*-+--v* <NU*
C     |                                                                         *
C     |    ^  ^  ^                                                     ^  ^     *
C     |   NU BC BC                                                    BC BC     *
C-----+-------------------------------------------------------------------------*
C     |    1     2     3           i-1   i     i+1         im-2  im-1  im       *
C     |                                                                         *
C     |                                S O U T H                                *
C-------------------------------------------------------------------------------*
C                                                                               *
C *******************************************************************************
C
      implicit none
C
      include 'pomNW.c'
C
      integer idx
      double precision ga,u1,wm,etide  !lyo:!wad:add etide
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
C     represent a choice which yields the most physically double precisionistic
C     result.
C
C     Elevation (in this application, elevation is not a primary
C     boundary condition):
C
        do j=1,jm
          elf(1,j)=elf(2,j)
          elf(im,j)=elf(imm1,j)
        end do
C
        do i=1,im
          elf(i,1)=elf(i,2)
          elf(i,jm)=elf(i,jmm1)
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
!tne:!wad:
C Very large (9m) tide-like BC (imposed by tidal inflow vel. not elev.)
C      note that ELW already include hhi (see SEAMOUNT)
!lyo:!wad:Specify time-dependent tidal tidal amplitude w/period=1day,
!         & amplitude=9m; note that tide is to be added to the geostro-
!         phically balanced ELW below.  So ETIDE is referenced wrt MSL:
!         [but perhaps it will be more direct to specify elevation??]
!        etide=-9.e0*sin(2.e0*pi*time/1.e0)
         etide=-tidamp*sin(2.e0*pi*time/1.e0) !minus so ebb near time=0
!
        do j=2,jmm1
C
C     East:
C
!tne:!wad:- for wad replace H with D
          uaf(im,j)=uabe(j)
     $               +rfe*sqrt(grav/d(imm1,j))
     $                         *(el(imm1,j)-ele(j))
          uaf(im,j)=ramp*uaf(im,j)
          vaf(im,j)=0.e0
C
C     West:
C
!tne:!wad:- for wad replace H with D add tide
          uaf(2,j)=uabw(j)
     $              -rfw*sqrt(grav/d(2,j))
     $                        *(el(2,j)-elw(j)-etide)
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
!tne:!wad:- for wad replace H with D
          vaf(i,jm)=vabn(i)
     $               +rfn*sqrt(grav/d(i,jmm1))
     $                         *(el(i,jmm1)-eln(i))
          vaf(i,jm)=ramp*vaf(i,jm)
          uaf(i,jm)=0.e0
C
C     South:
C
!tne:!wad:- for wad replace H with D
          vaf(i,2)=vabs(i)
     $              -rfs*sqrt(grav/d(i,2))
     $                        *(el(i,2)-els(i))
          vaf(i,2)=ramp*vaf(i,2)
          vaf(i,1)=vaf(i,2)
          uaf(i,1)=0.e0
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
C     Velocity (radiation conditions; smoothing is used in the direction
C     tangential to the boundaries):
C
        do k=1,kbm1
          do j=2,jmm1
C
C     East:
C
            ga=sqrt(d(im,j)/hmax)  !lyo:!wad:replace h w/d
            uf(im,j,k)=ga*(.25e0*u(imm1,j-1,k)+.5e0*u(imm1,j,k)
     $                     +.25e0*u(imm1,j+1,k))
     $                  +(1.e0-ga)*(.25e0*u(im,j-1,k)+.5e0*u(im,j,k)
     $                    +.25e0*u(im,j+1,k))
            vf(im,j,k)=0.e0
C
C     West:
C
            ga=sqrt(d(1,j)/hmax)  !lyo:!wad:replace h w/d
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
            ga=sqrt(d(i,jm)/hmax)  !lyo:!wad:replace h w/d
            vf(i,jm,k)=ga*(.25e0*v(i-1,jmm1,k)+.5e0*v(i,jmm1,k)
     $                     +.25e0*v(i+1,jmm1,k))
     $                  +(1.e0-ga)*(.25e0*v(i-1,jm,k)+.5e0*v(i,jm,k)
     $                    +.25e0*v(i+1,jm,k))
            uf(i,jm,k)=0.e0
C
C     South:
C
            ga=sqrt(d(i,1)/hmax)  !lyo:!wad:replace h w/d
            vf(i,2,k)=ga*(.25e0*v(i-1,3,k)+.5e0*v(i,3,k)
     $                    +.25e0*v(i+1,3,k))
     $                 +(1.e0-ga)*(.25e0*v(i-1,2,k)+.5e0*v(i,2,k)
     $                   +.25e0*v(i+1,2,k))
            vf(i,1,k)=vf(i,2,k)
            uf(i,1,k)=0.e0
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
C     East:
C
            u1=2.e0*u(im,j,k)*dti/(dx(im,j)+dx(imm1,j))
            if(u1.le.0.e0) then
              uf(im,j,k)=t(im,j,k)-u1*(tbe(j,k)-t(im,j,k))
              vf(im,j,k)=s(im,j,k)-u1*(sbe(j,k)-s(im,j,k))
            else
              uf(im,j,k)=t(im,j,k)-u1*(t(im,j,k)-t(imm1,j,k))
              vf(im,j,k)=s(im,j,k)-u1*(s(im,j,k)-s(imm1,j,k))
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
          do j=1,jm
C
C     East:
C
            u1=2.e0*u(im,j,k)*dti/(dx(im,j)+dx(imm1,j))
            if(u1.le.0.e0) then
              uf(im,j,k)=q2(im,j,k)-u1*(small-q2(im,j,k))
              vf(im,j,k)=q2l(im,j,k)-u1*(small-q2l(im,j,k))
            else
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
      include 'pomNW.c'
C
      integer idx
      double precision cl,denom
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
C     represent a choice which yields the most physically double precisionistic
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
      include 'pomNW.c'
C
      double precision depth,delx,tatm,satm
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
!
!lyo:!wad:!pom2k_bug:tsurf and ssurf were never defined, but should be:
             tsurf(i,j)=tb(i,j,1)
             ssurf(i,j)=sb(i,j,1)
!
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
!                                                                      !
!----------------------------------------------------------------------!
!lyo:!wad: Set up pdens before 1st call dens; used also in profq:      !
      do k=1,kbm1; do j=1,jm; do i=1,im
         pdens(i,j,k)=grav*rhoref*(-zz(k)*max(h(i,j)-hhi,0.e0))*1.e-5
         enddo; enddo; enddo
!                                                                      !
      call dens(sb,tb,rho)
!                                                                      !
!----------------------------------------------------------------------!
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
      include 'pomNW.c'
C
      double precision si(im,jm,kb),ti(im,jm,kb),rhoo(im,jm,kb)
      double precision cr,p,rhor,sr,tr,tr2,tr3,tr4  !,hij !lyo:!wad:define hij
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
!
!lyo:!wad:Keep "p" defined as in original pom w/o hhi (i.e. "h" defined
!         wrt MSL), but set p=0 for cells where input "h" indicates
!         "dry" (i.e. where h-hhi < 0); i.e. use "pdens":
!
            p=pdens(i,j,k)
!           p=grav*rhoref*(-zz(k)* h(i,j))*1.e-5
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
      include 'pomNW.c'
C
      double precision delz
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
      subroutine findpsi
C **********************************************************************
C *                                                                    *
C * ROUTINE NAME:  findpsi                                             *
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
      implicit none
C
      include 'pomNW.c'
C
      integer i,j
C
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
      call prxy('Streamfunction, psi from u              ',
     $          time,psi,im,iskp,jm,jskp,0.d0)
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
      call prxy('Streamfunction, psi from v              ',
     $          time,psi,im,iskp,jm,jskp,0.d0)
C
      return
C
      end
C
!----------------------------------------------------------------------!
!lyo:_20080415:This version (of file2ic) is basically still the old one!
!     from pom2k - I have not converted it to setup for double precisionistic      !
!     WAD runs - but see subr. wadseamount for an example of how-to.   !
!----------------------------------------------------------------------!
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
      include 'pomNW.c'
C
      double precision rad,re,dlat,dlon,cff
      integer i,j,k,m
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
       read(40,'(8E12.5)') h  !lyo:!wad:note:read +- topography
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
!
!lyo:_20080415:Move "calc. Coriolis etc" through "call areas_masks"
!     from below (after call dens) to here; otherwise, fsm will not
!     be defined but is used in subroutine dens. (This is a pom2k_bug).
C
C --- calc. Coriolis Parameter
C
        do j=1,jm
          do i=1,im
            cor(i,j)=2.*7.29E-5*sin(north_e(i,j)*rad)
!           aam2d(i,j)=aam(i,j,1) !lyo:_20080415:pom2k_bug:aam not defined
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
!                                                                      !
!----------------------------------------------------------------------!
!lyo:!wad: Set up pdens before 1st call dens; used also in profq:      !
      do k=1,kbm1; do j=1,jm; do i=1,im
         pdens(i,j,k)=grav*rhoref*(-zz(k)*max(h(i,j)-hhi,0.e0))*1.e-5
         enddo; enddo; enddo
!                                                                      !
      call dens(sb,tb,rho)
!                                                                      !
!----------------------------------------------------------------------!
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
      subroutine printall
C **********************************************************************
C *                                                                    *
C *                         POM2K SOURCE CODE                          *
C *                                                                    *
C * ROUTINE NAME:  printall                                            *
C *                                                                    *
C * FUNCTION    :  Prints a set of outputs to device 6                 *
C *                                                                    *
C *                Edit as approriate.                                 *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      integer io(100),jo(100),ko(100)
C
      include 'pomNW.c'
C
C     2-D horizontal fields:
C
          call prxy('Depth-averaged u, uab                   ',
     $              time,uab,im,iskp,jm,jskp,0.d0)
C
          call prxy('Depth-averaged v, vab                   ',
     $              time,vab,im,iskp,jm,jskp,0.d0)
C
          call prxy('Surface elevation, elb                  ',
     $              time,elb,im,iskp,jm,jskp,0.d0)
C
C     Calculate and print streamfunction:
C
          call findpsi
C
          if(mode.ne.2) then
C
C     2-D horizontal sections of 3-D fields:
C
C     Set levels for output:
C
            ko(1)=1
            ko(2)=kb/2
            ko(3)=kb-1
C
            call prxyz('x-velocity, u                           ',
     $                 time,u    ,im,iskp,jm,jskp,kb,ko,3,0.d0 )
C
            call prxyz('y-velocity, v                           ',
     $                 time,v    ,im,iskp,jm,jskp,kb,ko,3,0.d0 )
C
            ko(1)=2
            call prxyz('z-velocity, w                           ',
     $                 time,w    ,im,iskp,jm,jskp,kb,ko,3,0.d0 )
            ko(1)=1
C
            call prxyz('Potential temperature, t                ',
     $                 time,t    ,im,iskp,jm,jskp,kb,ko,3,1.d-2)
C
            call prxyz('Salinity, s                              ',
     $                 time,s    ,im,iskp,jm,jskp,kb,ko,3,1.d-2)
C
            call prxyz('(density-1000)/rhoref, rho              ',
     $                 time,rho  ,im,iskp,jm,jskp,kb,ko,3,1.d-5)
C
c           call prxyz('Turbulent kinetic energy x 2, q2        ',
c    $                 time,q2   ,im,iskp,jm,jskp,kb,ko,3,0.d0 )
C
c           call prxyz('Turbulent length scale, l               ',
c    $                 time,l    ,im,iskp,jm,jskp,kb,ko,3,0.d0 )
C
            call prxyz('Horizontal kinematic viscosity, aam     ',
     $                 time,aam  ,im,iskp,jm,jskp,kb,ko,3,0.d0 )
C
            call prxyz('Vertical kinematic viscosity, km        ',
     $                 time,km   ,im,iskp,jm,jskp,kb,ko,3,0.d0 )
C
c           call prxyz('Vertical kinematic diffusivity, kh      ',
c    $                 time,kh   ,im,iskp,jm,jskp,kb,ko,3,0.d0 )
C
C     Vertical sections of 3-D fields, normal to j-axis:
C
C     Set sections for output:
C
            jo(1)=1
            jo(2)=jm/2
            jo(3)=jm-1
C
            call prxz('x-velocity, u                           ',
     $                time,u    ,im,iskp,jm,kb,jo,3,0.d0 ,dt,zz)
C
            call prxz('y-velocity, v                           ',
     $                time,v    ,im,iskp,jm,kb,jo,3,0.d0 ,dt,zz)
C
            call prxz('z-velocity, w                           ',
     $                time,w    ,im,iskp,jm,kb,jo,3,0.d0 ,dt,z )
C
            call prxz('Potential temperature, t                ',
     $                time,t    ,im,iskp,jm,kb,jo,3,1.d-2,dt,zz)
C
            call prxz('Salinity, s                             ',
     $                time,s    ,im,iskp,jm,kb,jo,3,1.d-2,dt,zz)
C
            call prxz('(density-1000)/rhoref, rho              ',
     $                time,rho  ,im,iskp,jm,kb,jo,3,1.d-5,dt,zz)
C
c           call prxz('Turbulent kinetic energy x 2, q2        ',
c    $                time,q2   ,im,iskp,jm,kb,jo,3,0.d0 ,dt,z )
C
c           call prxz('Turbulent length scale, l               ',
c    $                time,l    ,im,iskp,jm,kb,jo,3,0.d0 ,dt,z )
C
c           call prxz('Horizontal kinematic viscosity, aam     ',
c    $                time,aam  ,im,iskp,jm,kb,jo,3,0.d0 ,dt,zz)
C
c           call prxz('Vertical kinematic viscosity, km        ',
c    $                time,km   ,im,iskp,jm,kb,jo,3,0.d0 ,dt,z )
C
c           call prxz('Vertical kinematic diffusivity, kh      ',
c    $                time,kh   ,im,iskp,jm,kb,jo,3,0.d0 ,dt,z )
C
C     Vertical sections of 3-D fields, normal to i-axis:
C
C     Set sections for output:
C
            io(1)=1
            io(2)=im/2
            io(3)=im-1
C
            call pryz('x-velocity, u                           ',
     $                time,u    ,im,jm,jskp,kb,io,3,0.d0 ,dt,zz)
C
            call pryz('y-velocity, v                           ',
     $                time,v    ,im,jm,jskp,kb,io,3,0.d0 ,dt,zz)
C
            call pryz('z-velocity, w                           ',
     $                time,w    ,im,jm,jskp,kb,io,3,0.d0 ,dt,zz)
C
            call pryz('Potential temperature, t                ',
     $                time,t    ,im,jm,jskp,kb,io,3,1.d-2,dt,zz)
C
c           call pryz('Salinity x rho / rhoref, s              ',
c    $                time,s    ,im,jm,jskp,kb,io,3,1.d-2,dt,zz)
C
c           call pryz('(density-1000)/rhoref, rho              ',
c    $                time,rho  ,im,jm,jskp,kb,io,3,1.d-5,dt,zz)
C
c           call pryz('Turbulent kinetic energy x 2, q2        ',
c    $                time,q2   ,im,jm,jskp,kb,io,3,0.d0 ,dt,z )
C
c           call pryz('Turbulent length scale, l               ',
c    $                time,l    ,im,jm,jskp,kb,io,3,0.d0 ,dt,z )
C
c           call pryz('Horizontal kinematic viscosity, aam     ',
c    $                time,aam  ,im,jm,jskp,kb,io,3,0.d0 ,dt,zz)
C
c           call pryz('Vertical kinematic viscosity, km        ',
c    $                time,km   ,im,jm,jskp,kb,io,3,0.d0 ,dt,z )
C
c           call pryz('Vertical kinematic diffusivity, kh      ',
c    $                time,kh   ,im,jm,jskp,kb,io,3,0.d0 ,dt,z )
C
          endif
C
      return
C
      end
C
      subroutine profq(sm,sh,dh,cc) !lyo:_20080415:see notes at beg.
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
      include 'pomNW.c'
C
      double precision sm(im,jm,kb),sh(im,jm,kb),cc(im,jm,kb)
      double precision gh(im,jm,kb),boygr(im,jm,kb),dh(im,jm)
     $                ,stf(im,jm,kb)
      double precision prod(im,jm,kb),kn(im,jm,kb)
      double precision a1,a2,b1,b2,c1
      double precision coef1,coef2,coef3,coef4,coef5
      double precision const1,e1,e2,ghc
      double precision p,sef,sp,tp
      double precision l0(im,jm)
      double precision cbcnst,surfl,shiw
      double precision utau2, df0,df1,df2
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
!
!lyo:!wad:Keep "p" defined as in original pom w/o hhi (i.e. "h" defined
!         wrt MSL), but set p=0 for cells where input "h" indicates
!         "dry" (i.e. where h-hhi < 0); i.e. use "pdens":
!
            p=pdens(i,j,k)*10.e0
!           p=grav*rhoref*(-zz(k)* h(i,j))*1.e-4
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
!lyo:!wad:Note in MAIN, call profq is before call proft/dens. So T&S
!         (hence rho) are defined on sigma-grid using (h+etf),
!         so use "h+et" here:
!    $                    /(dzz(k-1)* h(i,j))
     $                    /(dzz(k-1)* (h(i,j)+et(i,j)))
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
          vf(i,j,1)=0.
          vf(i,j,kb)=0.
          ee(i,j,2)=0.e0
          gg(i,j,2)=-kappa*z(2)*dh(i,j)*q2(i,j,2)
          vf(i,j,kb-1)=kappa*(1+z(kbm1))*dh(i,j)*q2(i,j,kbm1)
        end do
      end do
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            dtef(i,j,k)=dtef(i,j,k)
     $                   *(1.e0+e2*((1.e0/abs(z(k)-z(1))
     $                               +1.e0/abs(z(k)-z(kb)))
     $                                *l(i,j,k)/(dh(i,j)*kappa))**2)
          end do
        end do
      end do
      do k=3,kbm1
        do j=1,jm
          do i=1,im
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
C The following is to counter the problem of the ratio of two
C small numbers (l = q2l/q2). Two options are included below.
      do k=2,kbm1
        do j=1,jm
          do i=1,im
c           if(uf(i,j,k).le.small.or.vf(i,j,k).le.small) then
c             uf(i,j,k)=small
c             vf(i,j,k)=0.1*dt(i,j)*small
c           endif
          uf(i,j,k)=abs(uf(i,j,k))
          vf(i,j,k)=abs(vf(i,j,k))
          end do
        end do
      end do
C
C-----------------------------------------------------------------------
c
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
C  There are 2 options for kq which, unlike km and kh, was
C  was not derived by Mellor and Yamada but was purely
C  empirical based on neutral boundary layer data.
C  The choice is whether or not it should be subject to
C  the stability factor, sh. Generally, there is not a great
C  difference in output.
C
      do k=1,kb
        do j=1,jm
          do i=1,im
!                                                                      !
!----------------------------------------------------------------------!
!lyo:!wad:See O2006, Appendix A:                                       !
!                                                                      !
!           kn(i,j,k)=l(i,j,k)*sqrt(abs(q2(i,j,k)))
            kn(i,j,k)=(kappa*zsh+l(i,j,k))*sqrt(abs(q2(i,j,k)))
!                                                                      !
!----------------------------------------------------------------------!
!
!           kq(i,j,k)=(kn(i,j,k)*.41e0*sh(i,j,k)+kq(i,j,k))*.5e0
            kq(i,j,k)=(kn(i,j,k)*.20+kq(i,j,k))*.5e0
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
      include 'pomNW.c'
C
      double precision f(im,jm,kb),wfsurf(im,jm)
      double precision fsurf(im,jm),dh(im,jm)
      integer nbc
      double precision rad(im,jm,kb),r(5),ad1(5),ad2(5)
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
      include 'pomNW.c'
      double precision dh(im,jm)
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
      include 'pomNW.c'
      double precision dh(im,jm)
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
      subroutine prxy(label,time,a,im,iskp,jm,jskp,scala)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Writes a horizontal 2-D field.                      *
C *                                                                    *
C *                label ....... label for output                      *
C *                time ........ time (days)                           *
C *                a(im,jm,kb).. array to be printed                   *
C *                iskp ........ skipping interval for i               *
C *                jskp ........ skipping interval for j               *
C *                scala ....... < 0 for floating point numbers output *
C *                              0 for integer output, divisor for a   *
C *                                based on magnitudes of |a| values   *
C *                              > 0 for integer output, divisor for a *
C *                                given by scala                      *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      integer im,jm
      double precision a(im,jm)
      double precision time,scala
      integer iskp,jskp
      character label*(*)
      double precision amx,scale
      integer i,ib,ie,j,jwr,cols
C
      if(scala.ge.0.e0) then
        cols=24
      else
        cols=12
      endif
C
      if (scala.lt.0.e0) scale = 1.e0
      if (scala.eq.0.e0) then
        amx=1.e-12
        do j=1,jm,jskp
          do i=1,im,iskp
            amx=max(abs(a(i,j)),amx)
          end do
        end do
          scale=10.e0**(int(log10(amx)+100.e0)-103)
        endif
      if(scala.gt.0.e0) scale=scala
C
      write(6,1) label
    1 format(1x,a40/)
      write(6,2) time,scale
    2 format(' Time = ',f9.4,' days    multiply all values by ',1pe8.2)
C
      do ib=1,im,cols*iskp
C
        ie=ib+(cols-1)*iskp
        if(ie.gt.im) ie=im
C
        if(scala.ge.0.e0) then
          write(6,3) (i,i=ib,ie,iskp)
    3     format(/,2x,24i5,/)
        else
          write(6,4) (i,i=ib,ie,iskp)
    4     format(/,12i10,/)
        endif
C
        do j=1,jm,jskp
          jwr=jm+1-j
          if(scala.ge.0.e0) then
            write(6,5) jwr,(nint(a(i,jwr)/scale),i=ib,ie,iskp)
    5       format(1x,i3,24i5)
          else
            write(6,6) jwr,(a(i,jwr),i=ib,ie,iskp)
    6       format(1x,i2,12(e10.2))
          endif
        end do
C
        write(6,7)
    7   format(//)
C
      end do
C
      return
C
      end
C
      subroutine prxyz(label,time,a,im,iskp,jm,jskp,kb,ko,nko,scala)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Writes horizontal layers of a 3-D field with        *
C *                integers or floating point numbers.                 *
C *                                                                    *
C *                label ....... label for output                      *
C *                time ........ time (days)                           *
C *                a(im,jm,kb).. array to be printed                   *
C *                iskp ........ skipping interval for i               *
C *                jskp ........ skipping interval for j               *
C *                ko .......... 1-D array of k-indices for output     *
C *                nko ......... number of elements in ko              *
C *                scala ....... < 0 for floating point numbers output *
C *                              0 for integer output, divisor for a   *
C *                                based on magnitudes of |a| values   *
C *                              > 0 for integer output, divisor for a *
C *                                given by scala                      *
C *                                                                    *
C *                (NOTE that this combines functions of old prxyz and *
C *                 eprxyz)                                            *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      integer im,jm,kb
      double precision a(im,jm,kb)
      double precision time,scala
      integer ko(*)
      integer iskp,jskp,nko
      character label*(*)
      double precision amx,scale
      integer i,ib,ie,j,jwr,k,iko,cols
C
      if(scala.ge.0.e0) then
        cols=24
      else
        cols=12
      endif
C
      if (scala.lt.0.e0) scale = 1.e0
      if (scala.eq.0.e0) then
        amx=1.e-12
        do iko=1,nko
          k=ko(iko)
          do j=1,jm,jskp
            do i=1,im,iskp
              amx=max(abs(a(i,j,k)),amx)
            end do
          end do
        end do
          scale=10.e0**(int(log10(amx)+100.e0)-103)
        endif
      if(scala.gt.0.e0) scale=scala
C
      write(6,1) label
    1 format(1x,a40/)
      write(6,2) time,scale
    2 format(' Time = ',f9.4,' days    multiply all values by ',1pe8.2)
C
      do iko=1,nko
C
        k=ko(iko)
C
        write(6,3) k
    3   format(3x,/' Layer k = ',i2)
C
        do ib=1,im,cols*iskp
C
          ie=ib+(cols-1)*iskp
          if(ie.gt.im) ie=im
C
          if(scala.ge.0.e0) then
            write(6,4) (i,i=ib,ie,iskp)
    4       format(/,2x,24i5,/)
          else
            write(6,5) (i,i=ib,ie,iskp)
    5       format(/,12i10,/)
          endif
C
          do j=1,jm,jskp
            jwr=jm+1-j
            if(scala.ge.0.e0) then
              write(6,6) jwr,(nint(a(i,jwr,k)/scale),i=ib,ie,iskp)
    6         format(1x,i3,24i5)
            else
              write(6,7) jwr,(a(i,jwr,k),i=ib,ie,iskp)
    7         format(1x,i2,12(e10.2))
            endif
          end do
C
          write(6,8)
    8     format(//)
C
        end do
C
      end do
C
      return
C
      end
C
      subroutine prxz(label,time,a,im,iskp,jm,kb,jo,njo,scala,dt,zz)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Writes vertical section of a 3-D field, in the      *
C *                x- or i-direction .                                 *
C *                                                                    *
C *                label ....... label for output                      *
C *                time ........ time (days)                           *
C *                a(im,jm,kb).. array to be printed                   *
C *                iskp ........ skipping interval for i               *
C *                jo .......... 1-D array of j-indices for output     *
C *                njo ......... number of elements in jo              *
C *                scala ....... < 0 for floating point numbers output *
C *                              0 for integer output, divisor for a   *
C *                                based on magnitudes of |a| values   *
C *                              > 0 for integer output, divisor for a *
C *                                given by scala                      *
C *                dt(im,jm) ... total depth                           *
C *                zz(kb) ...... sigma coordinate                      *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      integer im,jm,kb
      double precision a(im,jm,kb),dt(im,jm),zz(kb)
      double precision time,scala
      integer jo(*)
      integer iskp,njo
      character label*(*)
      double precision amx,scale
      integer i,ib,ie,j,k,ijo,cols
C
      if(scala.ge.0.e0) then
        cols=24
      else
        cols=12
      endif
C
      if (scala.lt.0.e0) scale = 1.e0
      if (scala.eq.0.e0) then
        amx=1.e-12
        do  k=1,kb
          do ijo=1,njo
            j=jo(ijo)
            do i=1,im,iskp
              amx=max(abs(a(i,j,k)),amx)
            end do
          end do
        end do
          scale=10.e0**(int(log10(amx)+100.e0)-103)
        endif
      if(scala.gt.0.e0) scale=scala
C
      write(6,1) label
    1 format(1x,a40/)
      write(6,2) time,scale
    2 format(' Time = ',f9.4,' days    multiply all values by ',1pe8.2)
C
      do ijo=1,njo
C
        j=jo(ijo)
C
        write(6,3) j
    3   format(3x,/' Section j =',i3)
C
        do ib=1,im,cols*iskp
C
          ie=ib+(cols-1)*iskp
          if(ie.gt.im) ie=im
C
          if(scala.ge.0.e0) then
            write(6,4) (i,i=ib,ie,iskp)
    4       format(/,'    i =  ',24i5,/)
          else
            write(6,5) (i,i=ib,ie,iskp)
    5       format(/,'    i =  ',12i10,/)
          endif
C
          write(6,6) (nint(dt(i,j)),i=ib,ie,iskp)
    6     format(8x,'d =',24i5.0,/,'     z or zz')
C
          do k=1,kb
            if(scala.ge.0.e0) then
              write(6,7) k,zz(k),(nint(a(i,j,k)/scale),i=ib,ie,iskp)
    7         format(1x,i2,2x,f6.3,24i5)
            else
              write(6,8) k,zz(k),(a(i,j,k),i=ib,ie,iskp)
    8         format(1x,i2,2x,f6.3,12(e10.2))
            endif
          end do
C
          write(6,9)
    9     format(//)
C
        end do
C
      end do
C
      return
C
      end
C
      subroutine pryz(label,time,a,im,jm,jskp,kb,io,nio,scala,dt,zz)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Writes vertical section of a 3-D field, in the      *
C *                y- or j-direction.                                  *
C *                                                                    *
C *                label ....... label for output                      *
C *                time ........ time (days)                           *
C *                a(im,jm,kb).. array to be printed                   *
C *                jskp ........ skipping interval for j               *
C *                io .......... 1-D array of i-indices for output     *
C *                nio ......... number of elements in io              *
C *                scala ....... < 0 for floating point numbers output *
C *                              0 for integer output, divisor for a   *
C *                                based on magnitudes of |a| values   *
C *                              > 0 for integer output, divisor for a *
C *                                given by scala                      *
C *                dt(im,jm) ... total depth                           *
C *                zz(kb) ...... sigma coordinate                      *
C *                                                                    *
C **********************************************************************
C
      implicit none
      integer im,jm,kb
      double precision a(im,jm,kb),dt(im,jm),zz(kb)
      double precision time,scala
      integer io(*)
      integer jskp,nio
      character label*(*)
      double precision amx,scale
      integer i,j,jb,je,k,iio,cols
C
      if(scala.ge.0.e0) then
        cols=24
      else
        cols=12
      endif
C
      if (scala.lt.0.e0) scale = 1.e0
      if (scala.eq.0.e0) then
        amx=1.e-12
        do  k=1,kb
          do j=1,jm,jskp
            do iio=1,nio
              i=io(iio)
              amx=max(abs(a(i,j,k)),amx)
            end do
          end do
        end do
          scale=10.e0**(int(log10(amx)+100.e0)-103)
        endif
      if(scala.gt.0.e0) scale=scala
C
      write(6,1) label
    1 format(1x,a40/)
      write(6,2) time,scale
    2 format(' Time = ',f9.4,' days    multiply all values by ',1pe8.2)
C
      do iio=1,nio
C
        i=io(iio)
C
        write(6,3) i
    3   format(3x,/' Section i =',i3)
C
        do jb=1,jm,cols*jskp
C
          je=jb+(cols-1)*jskp
          if(je.gt.jm) je=jm
C
          if(scala.ge.0.e0) then
            write(6,4) (j,j=jb,je,jskp)
    4       format(/,'    j =  ',24i5,/)
          else
            write(6,5) (j,j=jb,je,jskp)
    5       format(/,'    j =  ',12i10,/)
          endif
C
          write(6,6) (nint(dt(i,j)),j=jb,je,jskp)
    6     format(8x,'d =',24i5.0,/,'     z or zz')
C
          do k=1,kb
            if(scala.ge.0.e0) then
              write(6,7) k,zz(k),(nint(a(i,j,k)/scale),j=jb,je,jskp)
    7         format(1x,i2,2x,f6.3,24i5)
            else
              write(6,8) k,zz(k),(a(i,j,k),j=jb,je,jskp)
    8         format(1x,i2,2x,f6.3,12(e10.2))
            endif
          end do
C
          write(6,9)
    9     format(//)
C
        end do
C
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
      include 'pomNW.c'
C
      double precision delh,delx,elejmid,elwjmid,ra,vel
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
      hmax=4500.e0   !tne:!wad: deep open ocean (orig val - not wad-related)
      do i=1,im
        do j=1,jm
C
          h(i,j)=hmax*(1.e0-delh
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
!
!lyo:!wad:!pom2k_bug:tsurf and ssurf were never defined, but should be:
             tsurf(i,j)=tb(i,j,1)
             ssurf(i,j)=sb(i,j,1)
!
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
!                                                                      !
!----------------------------------------------------------------------!
!lyo:!wad: Set up pdens before 1st call dens; used also in profq:      !
      do k=1,kbm1; do j=1,jm; do i=1,im
         pdens(i,j,k)=grav*rhoref*(-zz(k)*max(h(i,j)-hhi,0.e0))*1.e-5
         enddo; enddo; enddo
!                                                                      !
      call dens(sb,tb,rho)
!                                                                      !
!----------------------------------------------------------------------!
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
      include 'pomNW.c'
C
      double precision mean,del
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
      include 'pomNW.c'
C
      double precision ff(im,jm,kb)
      double precision
     $          xmassflux(im,jm,kb),ymassflux(im,jm,kb),zwflux(im,jm,kb)
      double precision sw
      double precision mol,abs_1,abs_2
      double precision value_min,epsilon
      double precision udx,u2dt,vdy,v2dt,wdz,w2dt
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
      include 'pomNW.c'
C
      double precision xflux(im,jm,kb),yflux(im,jm,kb)
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
c
!lyo:!wad:                                                             !
!                                                                      !
!----------------------------------------------------------------------!
!                                                                      !
!     Here are WAD-related subroutines in alphabatical order;          !
!     all names begin with "wad"                                       !
!                                                                      !
!----------------------------------------------------------------------!
!                                                                      !
      subroutine wadadvt2d(fb,f,ff,nitera,sw,dx,dy,u,v,art,aru,arv,
     1    fsm,dum,dvm,aam,dti2)
!                                                                      !
!     2-d version of Smolarkiewicz from /home/lyo/pom/pom2k/           !
!     pom2k.f's subroutine advt2:                                      !
!                                                                      !
!     This version gets rid of common block, i.e. the routine is now   !
!     self-contained, except that im & jm need to be specified below   !
!     if "include 'grid'" is NOT used; also, 2-d calculation, i.e,     !
!     solve for  D:                                                    !
!                                                                      !
!     d(D)/dt + d(UD)/dx + d(VD)/dy = Diffusion_of_(D)                 !
!                                                                      !
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
      double precision sw,dti2
      integer nitera,im,jm,kb  !lyo:!wad:kb not required but defined
                               !         in "grid" below
c
!     PARAMETER (IM=65,JM=49)
!     PARAMETER (IM=131,JM=99)
      include 'grid'
c
      double precision fb(im,jm),f(im,jm),ff(im,jm)
      double precision dx(im,jm),dy(im,jm),u(im,jm),v(im,jm),aam(im,jm)
      double precision fsm(im,jm),dum(im,jm),dvm(im,jm)
      double precision art(im,jm),aru(im,jm),arv(im,jm)

c     integer ix,jx
c     PARAMETER (IX=62,JX=18)
      double precision xflux(im,jm),yflux(im,jm)
      double precision fbmem(im,jm)
      double precision xmassflux(im,jm),ymassflux(im,jm)
      integer i,j,k,itera,imm1,jmm1
C
C     Calculate horizontal mass fluxes:
C
      imm1=im-1; jmm1=jm-1
c
        do j=1,jm
          do i=1,im
            xmassflux(i,j)=0.e0
            ymassflux(i,j)=0.e0
          end do
        end do
C
        do j=2,jmm1
          do i=2,im
            xmassflux(i,j)=0.5e0*(dy(i-1,j)+dy(i,j))*u(i,j)
          end do
        end do
C
        do j=2,jm
          do i=2,imm1
            ymassflux(i,j)=0.5e0*(dx(i,j-1)+dx(i,j))*v(i,j)
          end do
        end do
C
        do j=1,jm
          do i=1,im
            fbmem(i,j)=fb(i,j)
          end do
        end do
C
C     Start Smolarkiewicz scheme:
C
      do itera=1,nitera
C
C     Upwind advection scheme:
C
          do j=2,jm
            do i=2,im
              xflux(i,j)=0.5e0
     $                      *((xmassflux(i,j)+abs(xmassflux(i,j)))
     $                        *fbmem(i-1,j)+
     $                        (xmassflux(i,j)-abs(xmassflux(i,j)))
     $                        *fbmem(i,j))
C
              yflux(i,j)=0.5e0
     $                      *((ymassflux(i,j)+abs(ymassflux(i,j)))
     $                        *fbmem(i,j-1)+
     $                        (ymassflux(i,j)-abs(ymassflux(i,j)))
     $                        *fbmem(i,j))
            end do
          end do
C
C     Add net advective fluxes and step forward in time:
C
          do j=2,jmm1
            do i=2,imm1
              ff(i,j)=xflux(i+1,j)-xflux(i,j)
     $                 +yflux(i,j+1)-yflux(i,j)
              ff(i,j)=(fbmem(i,j)*art(i,j)
     $                   -dti2*ff(i,j))/(art(i,j))
            end do
          end do
C
C     Calculate antidiffusion velocity:
C
      call wadsmoladif(xmassflux,ymassflux,ff,sw,fsm,aru,arv,dti2)
c     call wadsmoladif(xmassflux,ymassflux,ff,sw,fsm,aru,arv,dti2,im,jm)
C
        do j=1,jm
          do i=1,im
              fbmem(i,j)=ff(i,j)
          end do

        end do
C
C     End of Smolarkiewicz scheme
C
      end do
C
C     Add horizontal diffusive fluxes:
C
        do j=2,jm
          do i=2,im
            xmassflux(i,j)=0.5e0*(aam(i,j)+aam(i-1,j))
            ymassflux(i,j)=0.5e0*(aam(i,j)+aam(i,j-1))
          end do
        end do
C
        do j=2,jm
          do i=2,im
           xflux(i,j)=-xmassflux(i,j)
     $                   *(fb(i,j)-fb(i-1,j))*dum(i,j)
     $                   *(dy(i,j)+dy(i-1,j))/(dx(i,j)+dx(i-1,j))
           yflux(i,j)=-ymassflux(i,j)
     $                   *(fb(i,j)-fb(i,j-1))*dvm(i,j)
     $                   *(dx(i,j)+dx(i,j-1))/(dy(i,j)+dy(i,j-1))
          end do
        end do
C
C     Add net horizontal fluxes and step forward in time:
C
        do j=2,jmm1
          do i=2,imm1
            ff(i,j)=ff(i,j)-dti2*(xflux(i+1,j)-xflux(i,j)
     $                               +yflux(i,j+1)-yflux(i,j))
     $                           /(art(i,j))
          end do
        end do
C
      return
C
      end
C
!                                                                      !
!----------------------------------------------------------------------!
!                                                                      !
      subroutine wadout
!                                                                      !
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  WAD Outputs (temporary subroutine to be merged with *
! *                the rest of the code later)                         *
! *                                                                    *
C *                Edit as approriate.                                 *
C *                                                                    *
! **********************************************************************
!                                                                      !
      implicit none
C
      integer i,j,k
C
      include 'pomNW.c'
C
C     2-D horizontal fields:
C
          call prxy('wetmask; =0 are land or dry, =1 are wet ',
     $              time,wetmask,im,iskp,jm,jskp,1.d0)
C
          do j=1,jm; do i=1,im
             tps(i,j)=fsm(i,j)-wetmask(i,j)
             enddo; enddo
          call prxy('fsm-wetmask; WAD: =1 are dry cells      ',
     $              time,tps,im,iskp,jm,jskp,1.d0)
C
          do j=1,jm; do i=1,im
             tps(i,j)=(elb(i,j)+hhi)*wetmask(i,j)
             enddo; enddo
          call prxy('elb+hhi (m; w.r.t MSL)                  ',
     $              time,tps,im,iskp,jm,jskp,0.d0)
c
          do j=1,jm; do i=1,im
             tps(i,j)=d(i,j)*fsm(i,j)
             enddo; enddo
          call prxy('Total Water Depth (m)                   ',
     $              time,tps,im,iskp,jm,jskp,0.d0)
c
ctne: -------------- save for matlab plot
c
       write(88,'(2i5,4f8.3)')im-2,jm-2,time*24.,99.9,99.9,99.9
      do j=2,jm-1; do i=2,im-1
       write(88,'(2i5,4f8.3)')i,j,tps(i,j),uab(i,j),vab(i,j),d(i,j)
      enddo; enddo
c
      return
      end
!                                                                      !
!----------------------------------------------------------------------!
!                                                                      !
      subroutine wadh
!                                                                      !
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  Sets up WAD definitions of topography [see Fig.1 of *
! *                O2006]                                              *
! *                                                                    *
! **********************************************************************
!                                                                      !
!     Input (through common block) is "h" defined as follows:          !
!                                                                      !
!     "h" is wrt MSL (Mean Sea Level) at which level h=0, and h<0 for  !
!     depth below the MSL, and h>0 for "land" above the MSL (i.e. for  !
!     marshes, wetlands, hills & mountains etc.)                       !
!                                                                      !
!     Also, nwad, hhi, slmax                                           !
!                                                                      !
!     Outputs:                                                         !
!                                                                      !
!     h --> (i.e. becomes) h_old_pom + hhi (which can be zero)         !
!       i.e. "h" is wrt to some high land-datum (Absolute Land Boundary!
!       or ALB), "hhi" meters above the MSL. All cells other than the  !
!       ALB cells are either always wet or can potentially become wet. !
!                                                                      !
!     Also, fsm, dum, dvm, wetmask, cell areas etc, and hmax           !
!       see more detailed definitions inside the subroutine            !
!                                                                      !
!----------------------------------------------------------------------!
      implicit none
!                                                                      !
      include 'pomNW.c'
!                                                                      !
      integer npos,nneg
      double precision hwatmin
      integer i,j
!                                                                      !
!----------------------------------------------------------------------!
!     Do a rudimentary check on "h:"                                   !
!                                                                      !
      npos=0; nneg=0
      do j=1,jm; do i=1,im
      if (h(i,j).ge.0.0) then
      npos=+1
      else
      nneg=-1
      endif
      if (npos*nneg .lt. 0) go to 101
      enddo; enddo
!                                                                      !
      write(6,'('' Stopped in subr. wadh; incorrect defn of h:'')')
      write(6,'('' h is one-sign only; see comments in subr.wadh'',/)')
      write(6,'('' npos,nneg ='',2i5,/)') npos,nneg
      call prxy('Undisturbed water depth, h',time,h,im,1,jm,1,0.d0)
      stop
!                                                                      !
 101  continue
!                                                                      !
      hkeep(:,:)=h(:,:) !Keep original for plots etc
!                                                                      !
!     HMIN & HC (all +ve) are defined as:                              !
!                                                                      !
!     h = hmin ---> absolute land area [never flooded], FSM=0.0;       !
!     h >= hc  ---> water depth wrt High water level, FSM=1.0;         !
!                                                                      !
!     Note that hc is defined in MAIN --                               !
!     hmin=0.01; hc=0.05; hco=hc       !hco is used to avoid round-
      hmin=0.01; hco=hc                !hco is used to avoid round-
      hc=hc*1.0001                     !off in initial elevation
!     hc=hc*1.01                       !tne:!wad:- using a larger
!                 multiplying factor, "1.01" instead of "1.0001",      !
!                 can help eliminate singular wet spots due to roundoff!
!                                                                      !
      write(6,'('' Outputs from subr. wadh ---------------------'',/)')
      write(6,'('' hmin,hc,hco ='',1p3e13.4,/)') hmin,hc,hco
!                                                                      !
!     Reverse sign of "h", and shift reference if hhi.ne.0, see below..!
!     (Note that hhi is defined in MAIN)                               !
      h(:,:)=-h(:,:)+hhi
!                                                                      !
      if (nwad.eq.1) then
!                                                                      !
!     If nwad=1, then hhi.ne.0.0, and line "h(:,:)=-h(:,:)+hhi" above  !
!     already shifts the reference level from MSL to ALB, the Absolute !
!     Land Boundary, according to O2006, Fig.1:                        !
!                                                                      !
      do j=1,jm; do i=1,im
      FSM(I,J)=1.; DUM(I,J)=1.; DVM(I,J)=1.
      IF(h(I,J).LT.0.0)THEN !Define absolute land areas (never flood):
                            ! temporarily set h=-1.0 (some -ve number)
      h(I,J)=-1.0; FSM(I,J)=0.; DUM(I,J)=0.; DVM(I,J)=0.
      ENDIF
      enddo; enddo
      DO J=1,JM-1; DO I=1,IM
      IF(FSM(I,J).EQ.0..AND.FSM(I,J+1).NE.0.)DVM(I,J+1)=0.
      enddo; enddo
      DO J=1,JM;   DO I=1,IM-1
      IF(FSM(I,J).EQ.0..AND.FSM(I+1,J).NE.0.)DUM(I+1,J)=0.
      enddo; enddo
!                                                                      !
!     The above yields the following definitions for h:                !
!                                                                      !
!     h=-1, absolute land areas [never flooded], FSM=0.0, always       !
!     h>=0, wet but potentially dry areas,  FSM also=1.0, always       !
!                                                                      !
!     Note that slpmax touches fsm=1 points only - wet or potentially  !
!     WAD cells; it changes "h" so a wet cell might become dry; so it  !
!     must be called here before wetmask is defined:                   !
!                                                                      !
!     Adjust bottom topography so that cell to cell variations         !
!     in h do not exceed parameter slmax:                              !
!                                                                      !
      if(slmax.lt.1.e0) call slpmax
!                                                                      !
!     We now define a wetmask, assuming that at t=0, the free surface  !
!     is at the MSL (i.e. elevation=0 wrt MSL).  If this run does not  !
!     begin with time=0 (or if other special free-surface distribution !
!     is desired), then wetmask needs to be redefined or read in e.g.  !
!     from a RESTART file, as done later in MAIN if NREAD.ne.0         !
!                                                                      !
      wetmask(:,:)=fsm(:,:)
      do j=1,jm; do i=1,im
      if (h(i,j).lt.hhi) wetmask(i,j)=0.0
      enddo; enddo
!                                                                      !
!     Finalize definitions of "h":                                     !
!                                                                      !
!     h=hmin, absolute land areas [never flooded], FSM=0.0, always     !
!     h>=hc,  wet but potentially dry areas,  FSM also=1.0, always     !
!                                   ...  but wetmask can be 0 or 1     !
!                                                                      !
      do j=1,jm; do i=1,im
      if (h(i,j).lt.0.0) then
      h(i,j)=hmin
      else
      h(i,j)=h(i,j)+hc
      endif
      enddo; enddo
!                                                                      !
      call areas_masks !Calc. cells' areas etc but do not alter fsm...
                       !if nwad=1
!                                                                      !
      else    !nwad.eq.0
!                                                                      !
      DO J=1,JM; DO I=1,IM
      IF (h(I,J).LE.0.0) h(I,J)=hmin
      ENDDO; ENDDO
!                                                                      !
!     Set minimum water depth to hwatmin:                              !
!     hwatmin=min water depth, must be >1 to match subr.areas_masks'   !
!     definition of land, and large enough to avoid grid cells from    !
!     becoming dry                                                     !
!                                                                      !
      hwatmin=10.0
      DO J=1,JM; DO I=1,IM
      IF (h(I,J).GT.hmin.and.h(I,J).LT.hwatmin) h(I,J)=hwatmin
      ENDDO; ENDDO
!                                                                      !
      call areas_masks      !Calc. cells' areas etc & fsm if nwad=0
!                                                                      !
!     Adjust bottom topography so that cell to cell variations         !
!     in h do not exceed parameter slmax:                              !
!                                                                      !
      if(slmax.lt.1.e0) call slpmax
!                                                                      !
      wetmask(:,:)=fsm(:,:) !wetmask is same as fsm if nwad=0
!                                                                      !
      endif   !if (nwad.eq.1) then...else...
!                                                                      !
!     Print to check:                                                  !
!                                                                      !
      CALL PRXY(' Input h ',TIME, hkeep(1,1),IM,1,JM,1,0.d0)
      CALL PRXY(' h after wadh ',TIME, h(1,1),IM,1,JM,1,0.d0)
      CALL PRXY(' FSM ',TIME, FSM(1,1),IM,1,JM,1,1.d0)
      CALL PRXY(' WETMASK ',TIME, WETMASK(1,1),IM,1,JM,1,1.d0)
      tps(:,:)=FSM(:,:)-WETMASK(:,:)
      CALL PRXY('FSM-WETMASK',TIME,tps(1,1),IM,1,JM,1,1.d0)
!                                                                      !
!     Double-check masks for nwad=0:                                   !
!                                                                      !
      if (nwad.eq.0) then
      do j=1,jm; do i=1,im
      if (fsm(i,j).ne.wetmask(i,j)) then
      write(6,'('' Stopped, nwad ='',i4)') nwad
      write(6,'('' fsm.ne.wetmask @ (i,j) ='',2i4,/)') i,j
      stop
      endif
      enddo; enddo
      endif
!                                                                      !
      return
      end
!                                                                      !
!----------------------------------------------------------------------!
!                                                                      !
      subroutine wadseamount
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Sets up for seamount problem w/wad.                 *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pomNW.c'
C
      double precision delh,delx,elejmid,elwjmid,ra,vel,hland  !lyo:!wad:define hland
      integer nvar                                 !lyo:!wad:define nvar
      integer i,j,k
C
C     Set delh > 1.0 for an island or delh < 1.0 for a seamount:
C
      delh=1.15e0   !tne:!wad: island (for wad) =0.9e0 for orig.seamount
C
C     Grid size:
C
      delx=8000.e0
c
!lyo:!wad:Define variable or uniform grid option:
      nvar=1  !=1 for variable grid; =0 otherwise
!lyo:!wad:Special test case for smaller delx=4km:
      if ((nvar.eq.0).or.(im.eq.131.and.jm.eq.99)) delx=delx*0.5
C
C     Radius island or seamount:
C
      ra=50000.e0   !tne:!wad:
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
!lyo:!wad:variable or uniform grid--
          dx(i,j)=delx-float(nvar)*delx*sin(pi*float(i)/float(im))/2.e0
          dy(i,j)=delx-float(nvar)*delx*sin(pi*float(j)/float(jm))/2.e0
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
!tne:!wad:!lyo:!wad:
!     Note that unlike original seamount, h<0 is now water below MSL, and
!     h>0 is above the MSL and is either land (if h>hland, see below)
!     or is potential wet-and-dry region
!
      hmax=50.e0; hland=5.e0 !note all pnts are potential WAD if we set
                             !hland >= -hmax*(1.e0-delh) (=7.5m);
                             !i.e. "if" below is NOT satisfied
!
      do i=1,im
        do j=1,jm
C
          h(i,j)=-hmax*(1.e0-delh
     $                          *exp(-((east_c(i,j)
     $                                   -east_c((im+1)/2,j))**2
     $                                +(north_c(i,j)
     $                                   -north_c(i,(jm+1)/2))**2)
     $                                /ra**2))
!lyo:!wad:Make region near top of seamount absolute land region:
          if(h(i,j).gt.hland) h(i,j)=hhi+1.e0
C
        end do
      end do
C
C     Close the north and south boundaries to form a channel:
C
      do i=1,im
        h(i,1)=hhi+1.e0  !tne:!wad:
        h(i,jm)=hhi+1.e0 !tne:!wad:
      end do
C
C     Calculate areas and masks:
C
!lyo:!wad:Define "h" consistent with WAD and calculate wetmask, cell areas
!         and fsm etc.
!     call areas_masks
      call wadh
C
C     Adjust bottom topography so that cell to cell variations
C     in h do not exceed parameter slmax:
C
!     if(slmax.lt.1.e0) call slpmax  !lyo:wad:now done in wadh above
C
C     Set initial conditions:
C
      do k=1,kbm1
        do j=1,jm
          do i=1,im
!lyo:wad:
!           tb(i,j,k)=5.e0+15.e0*exp(zz(k)*h(i,j)/1000.e0)-tbias
             if (h(i,j).gt.hhi) then !lyo:wad:stratified water below MSL:
                tb(i,j,k)=5.e0+15.e0*exp(zz(k)*(h(i,j)-hhi)/1000.e0)
             else                    !lyo:wad:well-mixed water above MSL:
                tb(i,j,k)=5.e0+15.e0
             endif
            tb(i,j,k)=tb(i,j,k)-tbias
!
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
!
!lyo:!wad:!pom2k_bug:tsurf and ssurf were never defined, but should be:
             tsurf(i,j)=tb(i,j,1)
             ssurf(i,j)=sb(i,j,1)
!
        end do
      end do
C
C     Initialise elb, etb, dt and aam2d:
C
      do j=1,jm
        do i=1,im
          elb(i,j)=(-e_atmos(i,j)-hhi)*fsm(i,j)      !lyo:!wad:
!lyo:!wad:note:The following "if" is not satisfied if nwad=0 since     !
!     then fsm=wetmask at all cells                                    !
!     if (fsm(i,j).ne.0.0.and.wetmask(i,j).eq.0.0) ! dry but not land
      if (fsm(i,j).ne.wetmask(i,j))                ! dry but not land
     &elb(i,j)=-h(i,j)+hco !note slightly smaller hco<hc to ensure that
                           !initially, cell is wet with elb+h=hco < hc
          etb(i,j)=elb(i,j)               !lyo:!wad:
          dt(i,j)=h(i,j)+elb(i,j)         !lyo:!wad:
          aam2d(i,j)=aam(i,j,1)
        end do
      end do
!                                                                      !
!----------------------------------------------------------------------!
!lyo:!wad:Make sure that initial wetmask remains compatible with       !
!     initial elb.  Also define marsh evaporation, +ve upward or       !
!     evaporate.  Also initialize wriv for rivers.                     !
!lyo:!wad:note:wriv is added here but the implementation following     !
!     Oey [1996; JPO, 26, 145-175; but see p.154 in particular] is not !
!     complete in this code.  E-mail me if want to implement it.       !
!                                                                      !
      if (nwad.eq.1) then
      do j=1,jm; do i=1,im
      wetmask(i,j)=fsm(i,j)
      if ((h(i,j)+elb(i,j)).le.hc) wetmask(i,j)=0.0
      enddo; enddo
      endif
!                                                                      !
      do j=1,jm; do i=1,im
      wmarsh(i,j)=0.0; wriv(i,j)=0.0
      if (fsm(i,j).ne.wetmask(i,j)) then           ! dry but not land
!     2 examples of wmarsh:
!--   wmarsh(i,j)=fsm(i,j)*10./art(i,j)
!--   wmarsh(i,j)=fsm(i,j)*(4.e-7)  !Gill's [1982] value
      endif
      enddo; enddo
!                                                                      !
!----------------------------------------------------------------------!
!lyo:!wad: Set up pdens before 1st call dens; used also in profq:      !
      do k=1,kbm1; do j=1,jm; do i=1,im
         pdens(i,j,k)=grav*rhoref*(-zz(k)*max(h(i,j)-hhi,0.e0))*1.e-5
         enddo; enddo; enddo
!                                                                      !
      call dens(sb,tb,rho)
!                                                                      !
!----------------------------------------------------------------------!
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
!lyo:!wad:Note keep (temporarily) all el* defined wrt MSL - ele, elw,
!         eln & els were initialized to be =0 in MAIN; then adjust
!         them later (below) with hhi:
!
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
!
!tne:!wad:!lyo:wad:
      do i=1,im
         eln(i)=(0.-hhi)*fsm(i,jm)
         els(i)=(0.-hhi)*fsm(i,2)
      enddo
      do j=1,jm
         ele(j)=(ele(j)-hhi)*fsm(im,j)
         elw(j)=(elw(j)-hhi)*fsm(2,j)
      enddo
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
!                                                                      !
!----------------------------------------------------------------------!
!                                                                      !
      subroutine wadsmoladif(xmassflux,ymassflux,ff,sw,fsm,aru,arv,dti2)
!                                                                      !
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
      integer im,jm,kb  !lyo:!wad:kb not required but defined
                        !         in "grid" below
c
!     PARAMETER (IM=65,JM=49)
!     PARAMETER (IM=131,JM=99)
      include 'grid'
c
      double precision ff(im,jm)
      double precision xmassflux(im,jm),ymassflux(im,jm)
      double precision fsm(im,jm),aru(im,jm),arv(im,jm)
      double precision sw,dti2
      double precision mol,abs_1,abs_2
      double precision value_min,epsilon
      double precision udx,u2dt,vdy,v2dt,wdz,w2dt
      integer i,j,k,imm1,jmm1
C
      parameter (value_min=1.e-9,epsilon=1.0e-14)
C
c
      imm1=im-1; jmm1=jm-1
c
C     Apply temperature and salinity mask:
C
        do i=1,im
          do j=1,jm
            ff(i,j)=ff(i,j)*fsm(i,j)
          end do
        end do
C
C     Recalculate mass fluxes with antidiffusion velocity:
C
        do j=2,jmm1
          do i=2,im
            if(ff(i,j).lt.value_min.or.
     $         ff(i-1,j).lt.value_min) then
              xmassflux(i,j)=0.e0
            else
              udx=abs(xmassflux(i,j))
              u2dt=dti2*xmassflux(i,j)*xmassflux(i,j)
     $              /(aru(i,j))
              mol=(ff(i,j)-ff(i-1,j))
     $             /(ff(i-1,j)+ff(i,j)+epsilon)
              xmassflux(i,j)=(udx-u2dt)*mol*sw

              abs_1=abs(udx)
              abs_2=abs(u2dt)
              if(abs_1.lt.abs_2) xmassflux(i,j)=0.e0
            end if
          end do
        end do
C
        do j=2,jm
          do i=2,imm1
            if(ff(i,j).lt.value_min.or.
     $         ff(i,j-1).lt.value_min) then
              ymassflux(i,j)=0.e0
            else
             vdy=abs(ymassflux(i,j))
             v2dt=dti2*ymassflux(i,j)*ymassflux(i,j)
     $             /(arv(i,j))
             mol=(ff(i,j)-ff(i,j-1))
     $            /(ff(i,j-1)+ff(i,j)+epsilon)
             ymassflux(i,j)=(vdy-v2dt)*mol*sw
             abs_1=abs(vdy)
             abs_2=abs(v2dt)
             if(abs_1.lt.abs_2) ymassflux(i,j)=0.e0
            end if
          end do
        end do
C
      return
C
      end
!                                                                      !
!lyo:!wad:END of WAD-related subroutines.                              !
!                                                                      !
!----------------------------------------------------------------------!
!                                                                      !
!
!rwnd: My procedures
!
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
      include 'pomNW.c'
C
      double precision datr(im, jm, kbm1, 2)
      double precision datu(imm1, jm, kbm1, 2)
      double precision datv(im, jmm1, kbm1, 2)
      character (len = 256) :: filename
!
      double precision rad,re,dlat,dlon,cff
      integer i,j,k,mm,ncid,varid
      double precision lom(12)   ! length of month
      data lom /31,28.25,31,30,31,30,31,31,30,31,30,31/
      rad=0.01745329
      re=6371.E3
!
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
!      h = -h       ! If using nwad=0 this line as well as `call wadh` isn't needed.
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
!      call upd_mnth(time, BC%ipl)
!      call clm_warp
!
!    Get temperature IC
!
      call check( nf90_inq_varid(ncid, "temp", varid) )

      if (BC%ipl) then
!    Read and interpolate between two months.
        if (mi.ne.1) then
          call check( nf90_get_var(ncid, varid, datr(:,:,kbm1:1:-1,:),
     $                            (/1,1,1,mi-1/), (/im,jm,kbm1,2/)) )
        else
          call check( nf90_get_var(ncid, varid, datr(:,:,kbm1:1:-1,1),
     $                            (/1,1,1,12/),(/im,jm,kbm1,1/)) )
          call check( nf90_get_var(ncid, varid, datr(:,:,kbm1:1:-1,2),
     $                            (/1,1,1,1/), (/im,jm,kbm1,1/)) )
        end if
        t(:,:,1:kbm1) = datr(:,:,1:kbm1,1)+fac*(datr(:,:,1:kbm1,2)
     $                 -datr(:,:,1:kbm1,1))

      else
!    Read a single month.
        call check( nf90_get_var(ncid, varid, t(:,:,kbm1:1:-1),
     $                          (/1,1,1,mi/), (/im,jm,kbm1,1/)) )
      end if
!
!    Get salinity IC
!
      call check( nf90_inq_varid(ncid, "salt", varid) )

      if (BC%ipl) then
!    Read and interpolate between two months.
        if (mi.ne.1) then
          call check( nf90_get_var(ncid, varid, datr,
     $                            (/1,1,1,mi-1/), (/im,jm,kbm1,2/)) )
        else
          call check( nf90_get_var(ncid, varid, datr(:,:,:,1),
     $                            (/1,1,1,12/),(/im,jm,kbm1,1/)) )
          call check( nf90_get_var(ncid, varid, datr(:,:,:,2),
     $                            (/1,1,1,1/), (/im,jm,kbm1,1/)) )
        end if

        do k=1,kbm1
          s(:,:,k) = datr(:,:,kb-k,1)+fac*(datr(:,:,kb-k,2)
     $              -datr(:,:,kb-k,1))
        end do

      else
!    Read a single month.
        call check( nf90_get_var(ncid, varid, s(:,:,kbm1:1:-1),
     $                          (/1,1,1,mi/), (/im,jm,kbm1,1/)) )
      end if
!
!   Read elevation
!
      elb = 0.
      el  = 0.

      if (IC%el) then

        call check( nf90_inq_varid(ncid, "zeta", varid) )

        if (BC%ipl) then

          if (mi.ne.1) then
            call check( nf90_get_var(ncid, varid, datr,
     $                              (/1,1,1,mi-1/), (/im,jm,1,2/)) )
          else
            call check( nf90_get_var(ncid, varid, datr(:,:,1,1),
     $                              (/1,1,1,12/),(/im,jm,1,1/)) )
            call check( nf90_get_var(ncid, varid, datr(:,:,1,2),
     $                              (/1,1,1,1/), (/im,jm,1,1/)) )
          end if

          elb(:,:) = datr(:,:,1,1)+fac*(datr(:,:,1,2)-datr(:,:,1,1))
        else
          call check( nf90_get_var(ncid, varid, elb,
     $                            (/1,1,mi/),(/im,jm,1/)) )
        end if
        el = elb      ! rwnd: Is this correct?

      end if

      call check( nf90_close(ncid) )
!----------------------------------------------------------------------!
!lyo:!wad: Set up pdens before 1st call dens; used also in profq:      !
      do k=1,kbm1; do j=1,jm; do i=1,im
         pdens(i,j,k)=grav*rhoref*(-zz(k)*max(h(i,j)-hhi,0.e0))*1.e-5
      enddo; enddo; enddo
!
! rmean must be read in bry
! Otherwise you MUST read or calculate rmean yourself. But for that you MUST use t(z) and s(z) that are functions of z-coordiante, not sigma.
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
            etb(i,j)=elb(i,j)       ! RWND //PV: 0.
            dt(i,j)=h(i,j)+elb(i,j) ! RWND //PV: h(i,j)
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
      !call wadh
      call areas_masks
!
! Call these before copy TS to boundaries.
!
      call bry(0)
      call bry(1)
      call bry(11)
      call bry(12)
      if (BC%wnd) call flux(5)
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
!!      TODO: this can be done simpler by touching only near-boundary cells
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
!     Get tsurf and ssurf
!      if (nbct==3 .and. nbcs==3) then   ! <-- nbc* are undefined here. Defined only in main program scope.
!        call bry(45)  ! Get fsurf anyway since it must be initialised.
!      else
!        write(*,*) '[!] WARNING! Cases except from nbc* = 3',
!     $               ' read ICOADS sst and sss.'
!        call bry(43)
!      end if
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
      end subroutine ncdf2ic
!
!
!rwnd: UNSTRATIFIED TEST CASE
!
      subroutine ncdf2ic_unstrat
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Sets up unstratified case.                          *
C *                                                                    *
C * This example reads IC from NetCDF files, generated by ROMS_tools.  *
C * Only minimal number of fields are read,                            *
C * while others are calculated here.  TODO: Read as much as we can    *
C *                                                                    *
C **********************************************************************
C
      use netcdf
      implicit none
C
      include 'pomNW.c'
C
      double precision datr(im, jm, kbm1, 2)
      double precision datu(imm1, jm, kbm1, 2)
      double precision datv(im, jmm1, kbm1, 2)
      character (len = 256) :: filename
!
      double precision rad,re,dlat,dlon,cff
      integer i,j,k,mm,ncid,varid
      double precision lom(12)   ! length of month
      data lom /31,28.25,31,30,31,30,31,31,30,31,30,31/
      rad=0.01745329
      re=6371.E3
!
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
!    Set temperature IC
!
      t = 15.
!
!    Set salinity IC
!
      s = 35.
!
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
!    Make flat bottom (debug reason)
!      h = hmax
      fsm(:,1) = 0
      fsm(1,:) = 0
!      h = -h       ! If using nwad=0 this line as well as `call wadh` isn't needed.
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
!      call upd_mnth(time, BC%ipl)
!      call clm_warp
!
!   Read elevation
!
      elb = 0.
      el  = 0.
!
      call check( nf90_close(ncid) )
!     Override tclim and sclim with IC. Or better comment out annual mean calculations.
      tclim = t
      sclim = s

      tsurf(:,:) = tclim(:,:,1)
      ssurf(:,:) = sclim(:,:,1)

!----------------------------------------------------------------------!
!lyo:!wad: Set up pdens before 1st call dens; used also in profq:      !
      do k=1,kbm1; do j=1,jm; do i=1,im
         pdens(i,j,k)=grav*rhoref*(-zz(k)*max(h(i,j)-hhi,0.e0))*1.e-5
      enddo; enddo; enddo
!
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
            dt(i,j)=h(i,j)+elb(i,j) ! RWND //PV: h(i,j)
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
      !call wadh
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
        !t(:,:,k) = t(:,:,k)*fsm(:,:)    ! whut? why?
        !s(:,:,k) = s(:,:,k)*fsm(:,:)    !
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
      if (BC%wnd) call flux(5)
!     Get tsurf and ssurf
!      if (nbct==3 .and. nbcs==3) then
!        call bry(45)  ! Get fsurf anyway since it must be initialised.
!      else
!        write(*,*) '[!] WARNING! Cases except from nbc* = 3',
!     $               ' read ICOADS sst and sss.'
!        call bry(43)
!      end if
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
!
      subroutine ncdf2ic_pom
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
      use date_utility
      implicit none
C
      include 'pomNW.c'
C
      character (len = 256) :: filename
      integer :: YYYY, MM, DD, hh, ii, ss
!
      double precision rad,re,dlat,dlon,dlnt,cff
      integer i,j,k,ncid,varid
      double precision lom(12)   ! length of month
      data lom /31,28.25,31,30,31,30,31,31,30,31,30,31/
      rad=0.01745329
      re=6371.E3
!
      write(6,'(/,'' Read grid and initial conditions '',/)')
C
C--- 1D ---
      filename = trim(pth_wrk)//trim(pth_grd)//
     $           trim(pfx_dmn)//"pom_grd.nc"
      write(*,*) "\\",trim(filename)
      call check( nf90_open(filename, NF90_NOWRITE, ncid) )
      call check( nf90_inq_varid(ncid, "z", varid) )
      call check( nf90_get_var(ncid, varid, z) )
      write(*, *) "[O] z retrieved"
      call check( nf90_inq_varid(ncid, "dz", varid) )
      call check( nf90_get_var(ncid, varid, dz) )
      write(*, *) "[O] dz retrieved"
      call check( nf90_inq_varid(ncid, "zz", varid) )
      call check( nf90_get_var(ncid, varid, zz) )
      write(*, *) "[O] zz retrieved"
      call check( nf90_inq_varid(ncid, "dzz", varid) )
      call check( nf90_get_var(ncid, varid, dzz) )
      write(*, *) "[O] dzz retrieved"
      call check( nf90_inq_varid(ncid, "h", varid) )
      call check( nf90_get_var(ncid, varid, h) )
      write(*, *) "[O] depth retrieved"
      call check( nf90_inq_varid(ncid, "fsm", varid) )
      call check( nf90_get_var(ncid, varid, fsm) )
      write(*, *) "[O] free surface mask retrieved"
      call check( nf90_inq_varid(ncid, "longitude", varid) )
      call check( nf90_get_var(ncid, varid, east_e) )
      write(*, *) "[O] longitude retrieved"
      call check( nf90_inq_varid(ncid, "latitude", varid) )
      call check( nf90_get_var(ncid, varid, north_e) )
      write(*, *) "[O] latitude retrieved"
      call check( nf90_close(ncid) )
      
!     Override hmax parameter
      hmax = maxval(h)

      filename = trim(pth_wrk)//trim(pth_grd)//
     $           trim(pfx_dmn)//"pom_clm.nc"
      write(*,*) "\\",trim(filename)      
      call check( nf90_open(filename, NF90_NOWRITE, ncid) )
C--- 3D ---
      call check( nf90_inq_varid(ncid, "Tclim", varid) )
      write(*,*) mi
      call check( nf90_get_var(ncid, varid, t, (/1,1,1,mi/),
     &                                         (/im,jm,kb,1/)) )
      write(*, *) "[O] potential temperature retrieved"
      call check( nf90_inq_varid(ncid, "Sclim", varid) )
      call check( nf90_get_var(ncid, varid, s, (/1,1,1,mi/),
     &                                         (/im,jm,kb,1/)) )
      write(*, *) "[O] salinity retrieved"
      call check( nf90_inq_varid(ncid, "Rmean", varid) )
      call check( nf90_get_var(ncid, varid, rmean, (/1,1,1,mi/),
     &                                         (/im,jm,kb,1/)) )
      write(*, *) "[O] rmean retrieved"
      call check( nf90_inq_varid(ncid, "el", varid) )
      call check( nf90_get_var(ncid, varid, elb, (/1,1,mi/),
     &                                         (/im,jm,1/)) )
      write(*, *) "[O] elevation retrieved"
      call check( nf90_close(ncid) )
      
      filename = trim(pth_wrk)//trim(pth_grd)//
     $           trim(pfx_dmn)//"pom_frc.nc"
      write(*,*) "\\",trim(filename)
      call check( nf90_open(filename, NF90_NOWRITE, ncid) )
      call check( nf90_inq_varid(ncid, "wU", varid) )
      call check( nf90_get_var(ncid, varid, wusurf, (/1,1,mi/),
     &                                            (/im,jm,1/)) )
      write(*, *) "[O] wusurf retrieved"
      call check( nf90_inq_varid(ncid, "wV", varid) )
      call check( nf90_get_var(ncid, varid, wvsurf, (/1,1,mi/),
     &                                            (/im,jm,1/)) )
      write(*, *) "[O] wvsurf retrieved"
      call check( nf90_close(ncid) )
!
! Make sure that cells with fsm=0. have depth h=1. for areas_masks subroutines to perform correctly.
! ...or just comment out `call areas_masks` and calculate areas and masks manually.
      do i=1,im
        do j=1,jm
          if (fsm(i,j)==0.) h(i,j)=1.
        end do
      end do

!   Simulate from zero elevation to avoid artificial waves during spin-up
!
!   Read elevation
!
!      elb = 0.
!      el  = 0.

      tclim = t
      sclim = s

!----------------------------------------------------------------------!
!lyo:!wad: Set up pdens before 1st call dens; used also in profq:      !
      do k=1,kbm1; do j=1,jm; do i=1,im
        pdens(i,j,k)=grav*rhoref*(-zz(k)*max(h(i,j)-hhi,0.e0))*1.e-5
      enddo; enddo; enddo
!
      write(*, *) "[+] Initial conditions have been read."
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
            dt(i,j)=h(i,j)+elb(i,j) ! RWND //PV: h(i,j)
          end do
        end do
C
        do j=1,jm
          do i=2,im-1
            dx(i,j)=0.5*rad*re*sqrt(((east_e(i+1,j)-east_e(i-1,j))
     $ *cos(north_e(i,j)*rad))**2+(north_e(i+1,j)-north_e(i-1,j))**2)
          end do
            dx(1,j)=dx(2,j)
            dx(im,j)=dx(im-1,j)
        end do
C
        do i=1,im
          do j=2,jm-1
            dy(i,j)=0.5*rad*re*sqrt(((east_e(i,j+1)-east_e(i,j-1))
     $ *cos(north_e(i,j)*rad))**2+(north_e(i,j+1)-north_e(i,j-1))**2)
          end do
            dy(i,1)=dy(i,2)
            dy(i,jm)=dy(i,jm-1)
        end do
C
C     Calculate areas and masks:
C
      !call wadh
      call areas_masks
!
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
      call bry(0)   ! always read rmean
      call bry(1)   ! always read TSclim
      call bry(11)  ! always get vertical...
      call bry(12)  ! ...and horizontal boundaries for TS
      if (BC%wnd) call flux(5)
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
!
      do i=1,im
        do j=1,jm
! --- calc. tilting angle of curvilinear grid
          if (j==jm) then
            dlon = (east_e(i,j)-east_e(i,j-1))*cos(north_e(i,j)*rad)
            dlat = north_e(i,j)-north_e(i,j-1)
          else
            dlon = (east_e(i,j+1)-east_e(i,j))*cos(north_e(i,j)*rad)
            dlat = north_e(i,j+1)-north_e(i,j)
          endif
          dlnt = (dlon**2+dlat**2)**.5
          rot(i,j) = asin(dlon/dlnt)
        end do
      end do
!
!     Set lateral boundary conditions, for use in subroutine bcond
!     set all=0 for closed BCs.
!     Values=0 for vel BC only, =1 is combination of vel+elev.
      rfe=1.e0
      rfw=1.e0
      rfn=1.e0
      rfs=1.e0
!
      return

      contains
      subroutine check(status)
        integer, intent ( in) :: status
        if(status /= nf90_noerr) then
          write(*,*) status
          stop "Stopped"
        end if
      end subroutine check
      end subroutine ncdf2ic_pom
      
      subroutine ncdf2ic_box
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Sets up my own boxproblem.                          *
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
      include 'pomNW.c'
C
      double precision datr(im, jm, kbm1, 2)
      double precision datu(imm1, jm, kbm1, 2)
      double precision datv(im, jmm1, kbm1, 2)
      character (len = 256) :: filename
!
      double precision rad,re,dlat,dlon,cff
      integer i,j,k,mm,ncid,varid
      double precision lom(12)   ! length of month
      data lom /31,28.25,31,30,31,30,31,31,30,31,30,31/
      rad=0.01745329
      re=6371.E3
!
      write(6,'(/,'' Read grid and initial conditions '',/)')
C
C--- 1D ---
      filename = trim(pth_wrk)//trim(pth_grd)//
     $           trim(pfx_dmn)//"pom_grd.nc"
      write(*,*) "\\",trim(filename)
      call check( nf90_open(filename, NF90_NOWRITE, ncid) )
      call check( nf90_inq_varid(ncid, "z", varid) )
      call check( nf90_get_var(ncid, varid, z) )
      write(*, *) "[O] z retrieved"
      call check( nf90_inq_varid(ncid, "dz", varid) )
      call check( nf90_get_var(ncid, varid, dz) )
      write(*, *) "[O] dz retrieved"
      call check( nf90_inq_varid(ncid, "zz", varid) )
      call check( nf90_get_var(ncid, varid, zz) )
      write(*, *) "[O] zz retrieved"
      call check( nf90_inq_varid(ncid, "dzz", varid) )
      call check( nf90_get_var(ncid, varid, dzz) )
      write(*, *) "[O] dzz retrieved"
      write(*, *) "[!] depth ignored"
      write(*, *) "[!] free surface mask ignored"
      call check( nf90_inq_varid(ncid, "longitude", varid) )
      call check( nf90_get_var(ncid, varid, east_e) )
      write(*, *) "[O] longitude retrieved"
      call check( nf90_inq_varid(ncid, "latitude", varid) )
      call check( nf90_get_var(ncid, varid, north_e) )
      write(*, *) "[O] latitude retrieved"
      call check( nf90_close(ncid) )
!
      fsm = 1.
      fsm( 1,:) = 0.
      fsm(im,:) = 0.
      fsm(:, 1) = 0.
      fsm(:,jm) = 0.
      
      hmax = 1000.
      
      h = hmax
      
      t = 10.
      s = 34.
!
!   Read elevation
!
      elb = 0.
      el  = 0.
      
      elb(:,:) = north_e(:,:)*.01
      el = elb

      tclim = t
      sclim = s

!----------------------------------------------------------------------!
!lyo:!wad: Set up pdens before 1st call dens; used also in profq:      !
      do k=1,kbm1; do j=1,jm; do i=1,im
         pdens(i,j,k)=grav*rhoref*(-zz(k)*max(h(i,j)-hhi,0.e0))*1.e-5
      enddo; enddo; enddo
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
            dt(i,j)=h(i,j)+elb(i,j) ! RWND //PV: h(i,j)
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
C --- calc. surface & lateral BC from climatology
C
        tsurf = t(:,:,1)
        ssurf = s(:,:,1)
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
      rmean = rho   
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
      end subroutine ncdf2ic_box
!
      function create_output(nprint) result(ncid)
C **********************************************************************
C *                                                                    *
C *                         POM2K SOURCE CODE                          *
C *                                                                    *
C * ROUTINE NAME:  create_output                                       *
C *                                                                    *
C * FUNCTION    :  Creates a netCDF file for further output            *
C *                                                                    *
C **********************************************************************
C
        use netcdf
        implicit none
C
        include 'pomNW.c'
C
        integer i, j, k, ncid, varid, fprint, nprint
        double precision bot_depth(im,jm,kb)
        character(len=256) filename
        integer dim_srho, dim_sw, dim_strim, dim_auxuv
        integer dim_lat, dim_lon, dim_time
        integer ktrim
C
        write(filename, '(3a,''.'',i4.4,''.nc'')') trim(pth_wrk),
     $             trim(pth_out),trim(title),nprint

        ktrim = kbm1

        call check( nf90_create(filename, NF90_CLASSIC_MODEL, ncid) )
        call check( nf90_def_dim(ncid, "time", NF90_UNLIMITED
     $          , dim_time) )
        call check( nf90_def_dim(ncid, "s_rho", kb, dim_srho) )
        call check( nf90_def_dim(ncid, "s_w", kb, dim_sw) )
        call check( nf90_def_dim(ncid, "s_trim", ktrim, dim_strim) )
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
        call check( nf90_def_var(ncid, "rot", NF90_DOUBLE,
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

          call check( nf90_def_var(ncid, "U", NF90_DOUBLE,
     $          (/ dim_lon, dim_lat, dim_strim, dim_time /), varid) )
          !call check( nf90_put_att(ncid, varid, "_FillValue", 0.) )

          call check( nf90_def_var(ncid, "V", NF90_DOUBLE,
     $          (/ dim_lon, dim_lat, dim_strim, dim_time /), varid) )
          !call check( nf90_put_att(ncid, varid, "_FillValue", 0.) )

          call check( nf90_def_var(ncid, "T", NF90_DOUBLE,
     $          (/ dim_lon, dim_lat, dim_strim, dim_time /), varid) )
          !call check( nf90_put_att(ncid, varid, "_FillValue", 0.) )

          call check( nf90_def_var(ncid, "S", NF90_DOUBLE,
     $          (/ dim_lon, dim_lat, dim_strim, dim_time /), varid) )
          !call check( nf90_put_att(ncid, varid, "_FillValue", 0.) )

          call check( nf90_def_var(ncid, "R", NF90_DOUBLE,
     $          (/ dim_lon, dim_lat, dim_strim, dim_time /), varid) )
          !call check( nf90_put_att(ncid, varid, "_FillValue", 0.) )

!          call check( nf90_def_var(ncid, "Rmean", NF90_DOUBLE,
!     $          (/ dim_lon, dim_lat, dim_strim, dim_time /), varid) )
          !call check( nf90_put_att(ncid, varid, "_FillValue", 0.) )

!          call check( nf90_def_var(ncid, "wusurf", NF90_DOUBLE,
!     $          (/ dim_lon, dim_lat, dim_time /), varid) )

!          call check( nf90_def_var(ncid, "wvsurf", NF90_DOUBLE,
!     $          (/ dim_lon, dim_lat, dim_time /), varid) )

        end if

        call check( nf90_def_var(ncid, "EL",  NF90_DOUBLE,
     $          (/ dim_lon, dim_lat, dim_time /), varid) )
        !call check( nf90_put_att(ncid, varid, "_FillValue", 0.) )

        call check( nf90_def_var(ncid, "UA",  NF90_DOUBLE,
     $          (/ dim_lon, dim_lat, dim_time /), varid) )
        !call check( nf90_put_att(ncid, varid, "_FillValue", 0.) )

        call check( nf90_def_var(ncid, "VA",  NF90_DOUBLE,
     $          (/ dim_lon, dim_lat, dim_time /), varid) )
        !call check( nf90_put_att(ncid, varid, "_FillValue", 0.) )

!        call check( nf90_def_var(ncid, "T_S",
!     $ NF90_DOUBLE, (/ dim_lon, dim_sw, dim_time /), varid) )
!        call check( nf90_def_var(ncid, "S_S",
!     $ NF90_DOUBLE, (/ dim_lon, dim_sw, dim_time /), varid) )
!        call check( nf90_def_var(ncid, "T_W",
!     $ NF90_DOUBLE, (/ dim_lat, dim_sw, dim_time /), varid) )
!        call check( nf90_def_var(ncid, "S_W",
!     $ NF90_DOUBLE, (/ dim_lat, dim_sw, dim_time /), varid) )
!        call check( nf90_def_var(ncid, "ELS",
!     $ NF90_DOUBLE, (/ dim_lon, dim_time /), varid) )
!        call check( nf90_def_var(ncid, "ELW",
!     $ NF90_DOUBLE, (/ dim_lat, dim_time /), varid) )

        call check( nf90_def_var(ncid, "PSI_nw",
     $ NF90_DOUBLE, (/ dim_lon, dim_lat, dim_time /), varid) )
        call check( nf90_def_var(ncid, "PSI_ew",
     $ NF90_DOUBLE, (/ dim_lon, dim_lat, dim_time /), varid) )

!    Debug vars
!        call check( nf90_def_var(ncid, "wusurf",
!     $ NF90_DOUBLE, (/ dim_lon, dim_lat, dim_time /), varid) )
!        call check( nf90_def_var(ncid, "wvsurf",
!     $ NF90_DOUBLE, (/ dim_lon, dim_lat, dim_time /), varid) )

        call check( nf90_enddef(ncid) )
!
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
        call check( nf90_inq_varid(ncid, "rot", varid) )
        call check( nf90_put_var(ncid, varid, rot) )
!        call check( nf90_inq_varid(ncid, "FSM", varid) )
!        call check( nf90_put_var(ncid, varid, fsm) )
!        call check( nf90_inq_varid(ncid, "CFL", varid) )
!        call check( nf90_put_var(ncid, varid, tps) )
        write(*, *) "NetCDF output has been initialized (",ncid,")."
        write(*, *) "Filepath: ",trim(filename)
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
!
      subroutine findpsi2nc(ncid, tind)
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
      include 'pomNW.c'
C
      integer i,j, tind, ncid, varid
      character(len=256) filename
C
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
      subroutine ncflush(ncid, ri)
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
        include 'pomNW.c'
C
        logical :: NOK
        integer :: count
        integer :: i, j, k, ncid, varid, status
        double precision :: vtot, tavg, atot, eavg, qavg, qtot, mtot
        double precision :: darea,dvol,mtot2(im),mtot3,vtot2(im),vtot3
        double precision :: volAcc(im,kbm1),masAcc(im,kbm1),mtotAcc(im)
        integer nlyrs, fi, ri ! fi - file index, ri - record index
        parameter (nlyrs = kbm1)
        integer :: lyrs(nlyrs)
        !data lyrs /1,2,15/
        character(len=256) filename
!
        write(*,*) "[-] Writing record #", ri
        count = 0
        NOK = .true.
        do i=1,nlyrs
          lyrs(i) = i
        end do
     
        call check( nf90_inq_varid(ncid, "Time", varid) )
!        do while (NOK)
!          count = count+1
        call check( nf90_put_var(ncid, varid, time, (/ri/)) )
!          if ((status.eq.nf90_noerr).or.(count.gt.3)) NOK = .false.
!        end do
C
          vtot=0.d0
          atot=0.d0
          qtot=0.d0
          mtot=0.d0
          tavg=0.d0
          qavg=0.d0
          eavg=0.d0
C
!          volAcc = 0.
!          masAcc = 0.
          do j=1,jm
            do i=1,im
              darea=dx(i,j)*dy(i,j)*wetmask(i,j)
              do k=1,kbm1
                dvol=darea*dt(i,j)*dz(k)
!                volAcc(i,k) = volAcc(i,k)+dvol
                mtot=mtot+(rho(i,j,k)*rhoref+1000.)*dvol
!                masAcc(i,k) = masAcc(i,k)+(rho(i,j,k)*rhoref+1000.)*dvol
                vtot=vtot+dvol
                tavg=tavg+tb(i,j,k)*dvol
                qavg=qavg+q2(i,j,k)*dvol
              end do
              atot=atot+darea
              eavg=eavg+et(i,j)*darea
            end do
          end do
          
!          mtot2 = 0.
!          vtot2 = 0.
!          do k=1,kbm1
!            do i=1,im
!              vtot2(i) = vtot2(i)+volAcc(i,k)
!              mtot2(i) = mtot2(i)+masAcc(i,k)
!            end do
!          end do
!          
!          mtot3 = 0.
!          vtot3 = 0.
!          do i=1,im
!            mtot3 = mtot3+mtot2(i)
!            vtot3 = vtot3+vtot2(i)
!          end do
          
!          write(*,*) "Total vol  (old): ",vtot
!          write(*,*) "Total vol  (new): ",vtot3
!          vtot3 = vtot3-vtot
!          write(*,*) "Total volume discreapncy between ",
!     $               "two sum methods is: ",vtot3
!          write(*,*) "Total mass (old): ",mtot
!          write(*,*) "Total mass (new): ",mtot3
!          mtot3 = mtot3-mtot
!          write(*,*) "Total mass discreapncy between ",
!     $               "two sum methods is: ",mtot3
C
          eavg=eavg/atot
          tavg=tavg/vtot
          qtot=qavg
          qavg=qtot/vtot
C
        call check( nf90_inq_varid(ncid, "Tavg", varid) )
        call check( nf90_put_var(ncid, varid, tavg, (/ri/)) )
!        call check( nf90_inq_varid(ncid, "Vtot", varid) )
!        call check( nf90_put_var(ncid, varid, vtot, (/ptime/)) )
        call check( nf90_inq_varid(ncid, "Eavg", varid) )
        call check( nf90_put_var(ncid, varid, eavg, (/ri/)) )
!        call check( nf90_inq_varid(ncid, "Etot", varid) )
!        call check( nf90_put_var(ncid, varid, eavg+qavg, (/ptime/)) )
        call check( nf90_inq_varid(ncid, "Qavg", varid) )
        call check( nf90_put_var(ncid, varid, qavg, (/ri/)) )
        call check( nf90_inq_varid(ncid, "Mtot", varid) )
        call check( nf90_put_var(ncid, varid, mtot, (/ri/)) )
!
        if (mode.ge.3) then
          call check( nf90_inq_varid(ncid, "U", varid) )
          call check( nf90_put_var(ncid, varid,   u(:,:,lyrs)
     $     ,(/1,1,1,ri/),(/im,jm,nlyrs,1/)) )
          call check( nf90_inq_varid(ncid, "V", varid) )
          call check( nf90_put_var(ncid, varid,   v(:,:,lyrs)
     $     ,(/1,1,1,ri/),(/im,jm,nlyrs,1/)) )
          call check( nf90_inq_varid(ncid, "T", varid) )
          call check( nf90_put_var(ncid, varid,   t(:,:,lyrs)
     $     ,(/1,1,1,ri/),(/im,jm,nlyrs,1/)) )
          call check( nf90_inq_varid(ncid, "S", varid) )
          call check( nf90_put_var(ncid, varid,   s(:,:,lyrs)
     $     ,(/1,1,1,ri/),(/im,jm,nlyrs,1/)) )
          call check( nf90_inq_varid(ncid, "R", varid) )
          call check( nf90_put_var(ncid, varid, rho(:,:,lyrs)
     $     ,(/1,1,1,ri/),(/im,jm,nlyrs,1/)) )
        end if
        
        call check( nf90_inq_varid(ncid, "EL", varid) )
        call check( nf90_put_var(ncid, varid, elb, (/1,1,ri/)
     $   ,(/im,jm,1/)) )
        call check( nf90_inq_varid(ncid, "UA", varid) )
        call check( nf90_put_var(ncid, varid, uab, (/1,1,ri/)
     $   ,(/im,jm,1/)) )
        call check( nf90_inq_varid(ncid, "VA", varid) )
        call check( nf90_put_var(ncid, varid, vab, (/1,1,ri/)
     $   ,(/im,jm,1/)) )
Cc
        call findpsi2nc(ncid, ri)
        call check( nf90_sync(ncid) )
C
        return
C
        contains
          subroutine check(status)
            integer, intent ( in) :: status
        if(status /= nf90_noerr) then
          select case (status)
              case (nf90_ebadid)
                write(*,*) "Not a netCDF id"
              case (nf90_eexist)
                write(*,*) "NetCDF file exists && NC_NOCLOBBER"
              case (nf90_einval)
                write(*,*) "Invalid Argument"
              case (nf90_eperm)
                write(*,*) "Write to read-only"
              case (nf90_enotindefine)
                write(*,*) "Operation not allowed in data mode"
              case (nf90_eindefine)
                write(*,*) "Operation not allowed in define mode"
              case (nf90_einvalcoords)
                write(*,*) "Index exeeds dimension bound"
              case (nf90_emaxdims)
                write(*,*) "NC_MAX_DIMS exceeded"
              case (nf90_enameinuse)
                write(*,*) "String match to name in use"
              case (nf90_enotatt)
                write(*,*) "Attribute not found"
              case (nf90_emaxatts)
                write(*,*) "NC_MAX_ATTRS exceeded"
              case (nf90_ebadtype)
                write(*,*) "Not a netCDF data type"
              case (nf90_ebaddim)
                write(*,*) "Invalid dimension or name"
              case (nf90_eunlimpos)
                write(*,*) "NC_UNLIMITED in the wrong index"
              case (nf90_emaxvars)
                write(*,*) "NC_MAX_VARS exceeded"
              case (nf90_enotvar)
                write(*,*) "Variable not found"
              case (nf90_eglobal)
                write(*,*) "Action prohibited on NC_GLOBAL varid"
              case (nf90_enotnc)
                write(*,*) "Not a netcdf file"
              case (nf90_ests)
                write(*,*) "In Fortran, string too short"
              case (nf90_emaxname)
                write(*,*) "NC_MAX_NAME exceeded"
              case (nf90_eunlimit)
                write(*,*) "NC_UNLIMITED size already in use"
              case (nf90_enorecvars)
                write(*,*) "nc_rec op when there are no record vars"
              case (nf90_echar)
                write(*,*) "Attempt to convert between text & numbers"
              case (nf90_eedge)
                write(*,*) "Edge+start exceeds dimension bound"
              case (nf90_estride)
                write(*,*) "Illegal stride"
              case (nf90_ebadname)
                write(*,*) "Attribute or variable name contains "
     $                    ,"illegal characters"
              case (nf90_erange)
                write(*,*) "Math result not representable"
              case (nf90_enomem)
                write(*,*) "Memory allocation (malloc) failure"
              case (nf90_evarsize)
                write(*,*) "One or more variable sizes violate "
     $                    ,"format constraints"
              case (nf90_edimsize)
                write(*,*) "Invalid dimension size"
              case (nf90_etrunc)
                write(*,*) "File likely truncated or possibly corrupted"
              case default
                write(*,*) "Undefined netCDF error (",status,")"
          end select
          stop "Stopped"
        end if
          end subroutine check
C
      end
!
      subroutine flux(idx)
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  Reads (if necessary) forcing dataset                *
! *                and interpolates in time.                           *
! *                                                                    *
! **********************************************************************
!
      use netcdf
      implicit none
!
      include 'pomNW.c'
!
      integer, intent (in) :: idx
      integer :: i,j, ncid, varid
      double precision :: qq, dq, sst
      character(len=256) filename

      double precision :: hf_fac = -4.1876e6 ! Heat flux convertion factor
!
      select case (idx)

        case (4) ! Long wave radiation and other fluxes

          if (mi.ne.rf_wtsur) then
!
            rf_wtsur = mi
!
            write(*,*) "[@] TODO: Implement long wave radiation",
     $                 " BC reading."
!
          end if

          return
!
        case (3) ! Short wave radiation flux
!
          if (mi.ne.rf_swrad) then
!
            rf_swrad = mi
!
            write(*,*) "[@] TODO: Implement short wave radiation",
     $                 " BC reading."
!
          end if

          return
!
        case (5) ! momentum flux
!
          if (mi.ne.rf_wsurf) then
!
            rf_wsurf = mi
!
            filename = trim(pth_wrk)//trim(pth_flx)//
     $                 trim(pfx_dmn)//"pom_frc.nc"
            call check( nf90_open(filename, NF90_NOWRITE, ncid) )
!
            call check( nf90_inq_varid(ncid, "wU", varid) )

            call check( nf90_get_var(ncid, varid, wusurf
     $                   ,(/1,1,mi/), (/im,jm,1/)) )

            if (BC%ipl) then
              if (mi.ne.1) then
                call check( nf90_get_var(ncid, varid, wusurfb
     $                   ,(/1,1,mi-1/), (/im,jm,1/)) )
              else
                call check( nf90_get_var(ncid, varid, wusurfb
     $                   ,(/1,1,12/),   (/im,jm,1/)) )
              end if
            end if

            call check( nf90_inq_varid(ncid, "wV", varid) )

            call check( nf90_get_var(ncid, varid, wvsurf
     $                   , (/1,1,mi/),   (/im,jmm1,1/)) )

            if (BC%ipl) then
              if (mi.ne.1) then
                call check( nf90_get_var(ncid, varid, wvsurfb
     $                   , (/1,1,mi-1/), (/im,jm,1/)) )
              else
                call check( nf90_get_var(ncid, varid, wvsurfb
     $                   , (/1,1,12/),   (/im,jm,1/)) )
              end if
            end if

            call check( nf90_close(ncid) )

            ! Fill in boundaries
            wusurf(1,:) = wusurf(2,:)
            wvsurf(:,1) = wvsurf(:,2)
            ! Taper wind stress along the boundary
            do i = 2, imm1
              do j = 2, jmm1
                wusurf(i,j) = .25*wusurf(i,j)
     $                       *(fsm(i,j+1)+fsm(i+1,j)
     $                        +fsm(i,j-1)+fsm(i-1,j)) ! TODO: or use fsm? Fsm may be preferrable for i=1, j=1
                wvsurf(i,j) = .25*wvsurf(i,j)
     $                       *(fsm(i,j+1)+fsm(i+1,j)
     $                        +fsm(i,j-1)+fsm(i-1,j)) !
              end do
            end do
            do i = 2, imm1
              wusurf(i, 1) = wusurf(i, 1)/3
     $                      *(fsm(i,   2)+fsm(i+1, 1)
     $                       +fsm(i-1, 1))
              wusurf(i,jm) = wusurf(i,jm)/3
     $                      *(fsm(i,jmm1)+fsm(i+1,jm)
     $                       +fsm(i-1,jm))
            end do
            do j = 2, jmm1
              wvsurf( 1,j) = wvsurf( 1,j)/3
     $                      *(fsm( 1,j+1)+fsm(   2,j)
     $                       +fsm( 1,j-1))
              wvsurf(im,j) = wvsurf(im,j)/3
     $                      *(fsm(im,j+1)+fsm(imm1,j)
     $                       +fsm(im,j-1))
            end do
            wusurf( 1, 1) = .5*wusurf( 1, 1)
     $                     *(fsm( 1,   2)+fsm(   2,   1))
            wusurf(im, 1) = .5*wusurf(im, 1)
     $                     *(fsm(im,   2)+fsm(imm1,   1))
            wusurf( 1,jm) = .5*wusurf( 1,jm)
     $                     *(fsm( 2,  jm)+fsm(   2,jmm1))
            wusurf(im,jm) = .5*wusurf(im,jm)
     $                     *(fsm(im,jmm1)+fsm(imm1,  jm))
            wvsurf( 1, 1) = .5*wvsurf( 1, 1)
     $                     *(fsm( 1,   2)+fsm(   2,   1))
            wvsurf(im, 1) = .5*wvsurf(im, 1)
     $                     *(fsm(im,   2)+fsm(imm1,   1))
            wvsurf( 1,jm) = .5*wvsurf( 1,jm)
     $                     *(fsm( 2,  jm)+fsm(   2,jmm1))
            wvsurf(im,jm) = .5*wvsurf(im,jm)
     $                     *(fsm(im,jmm1)+fsm(imm1,  jm))

!
            if (BC%ipl) then
                
              wusurff = wusurf
              wvsurff = wvsurf

              ! Taper for *b fields too.

              ! Fill in boundaries
              wusurfb(1,:) = wusurfb(2,:)
              wvsurfb(:,1) = wvsurfb(:,2)
              ! Taper wind stress along the boundary
              do i = 2, imm1
                do j = 2, jmm1
                  wusurfb(i,j) = .25*wusurfb(i,j)
     $                          *(fsm(i,j+1)+fsm(i+1,j)
     $                           +fsm(i,j-1)+fsm(i-1,j)) ! TODO: or use fsm?
                  wvsurfb(i,j) = .25*wvsurfb(i,j)
     $                          *(fsm(i,j+1)+fsm(i+1,j)
     $                           +fsm(i,j-1)+fsm(i-1,j)) !
                end do
              end do
              do i = 2, imm1
                wusurfb(i, 1) = wusurfb(i, 1)/3
     $                         *(fsm(i,   2)+fsm(i+1, 1)
     $                          +fsm(i-1, 1))
                wusurfb(i,jm) = wusurfb(i,jm)/3
     $                         *(fsm(i,jmm1)+fsm(i+1,jm)
     $                          +fsm(i-1,jm))
              end do
              do j = 2, jmm1
                wvsurfb( 1,j) = wvsurfb( 1,j)/3
     $                         *(fsm( 1,j+1)+fsm(   2,j)
     $                          +fsm( 1,j-1))
                wvsurfb(im,j) = wvsurfb(im,j)/3
     $                         *(fsm(im,j+1)+fsm(imm1,j)
     $                          +fsm(im,j-1))
              end do
              wusurfb( 1, 1) = .5*wusurfb( 1, 1)
     $                        *(fsm( 1,   2)+fsm(   2,   1))
              wusurfb(im, 1) = .5*wusurfb(im, 1)
     $                        *(fsm(im,   2)+fsm(imm1,   1))
              wusurfb( 1,jm) = .5*wusurfb( 1,jm)
     $                        *(fsm( 2,  jm)+fsm(   2,jmm1))
              wusurfb(im,jm) = .5*wusurfb(im,jm)
     $                        *(fsm(im,jmm1)+fsm(imm1,  jm))
              wvsurfb( 1, 1) = .5*wvsurfb( 1, 1)
     $                        *(fsm( 1,   2)+fsm(   2,   1))
              wvsurfb(im, 1) = .5*wvsurfb(im, 1)
     $                        *(fsm(im,   2)+fsm(imm1,   1))
              wvsurfb( 1,jm) = .5*wvsurfb( 1,jm)
     $                        *(fsm( 2,  jm)+fsm(   2,jmm1))
              wvsurfb(im,jm) = .5*wvsurfb(im,jm)
     $                        *(fsm(im,jmm1)+fsm(imm1,  jm))

            end if

            write(*,*) "[-] Read w*surf:  ", mi

          end if

          if (BC%ipl) then
            wusurf = wusurfb+fac*(wusurff-wusurfb)
            wvsurf = wvsurfb+fac*(wvsurff-wvsurfb)
          end if

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
      subroutine flux_roms(idx)
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
      include 'pomNW.c'
!
      integer, intent (in) :: idx
      integer :: i,j, ncid, varid
      double precision :: qq, dq, sst
      character(len=256) filename

      double precision :: hf_fac = -4.1876e6 ! Heat flux convertion factor
C
      select case (idx)

        case (4) ! Long wave radiation and other fluxes

          if (mi.ne.rf_wtsur) then
C
            rf_wtsur = mi
C
            filename = trim(pth_wrk)//trim(pth_flx)//
     $                 trim(pfx_dmn)//"roms_frc.nc"
            call check( nf90_open(filename, NF90_NOWRITE, ncid) )
            call check( nf90_inq_varid(ncid, "shflux", varid) )
            if (mi.ne.1) then
              call check(nf90_get_var(ncid, varid, qqb,(/1,1,mi-1/)))
              call check( nf90_get_var(ncid, varid, qqf, (/1,1,mi/)))
            else
              call check( nf90_get_var(ncid, varid, qqb,(/1,1,12/)))
              call check( nf90_get_var(ncid, varid, qqf, (/1,1,1/)))
            end if
            call check( nf90_inq_varid(ncid, "dQdSST", varid) )
            if (mi.ne.1) then
              call check(nf90_get_var(ncid, varid, dqb,(/1,1,mi-1/)))
              call check( nf90_get_var(ncid, varid, dqf, (/1,1,mi/)))
            else
              call check(nf90_get_var(ncid, varid, dqb, (/1,1,12/)))
              call check( nf90_get_var(ncid, varid, dqf, (/1,1,1/)))
            end if
            ! We need to read tsurf again, since we cannot be sure whether it has been read or not.
            call check( nf90_inq_varid(ncid, "SST", varid) )
            if (mi.ne.1) then
              call check(nf90_get_var(ncid, varid, tsurfb,(/1,1,mi-1/)))
              call check( nf90_get_var(ncid, varid, tsurff, (/1,1,mi/)))
            else
              call check( nf90_get_var(ncid, varid, tsurfb, (/1,1,12/)))
              call check( nf90_get_var(ncid, varid, tsurff, (/1,1,1/)))
            end if

            call check( nf90_close(ncid) )
C
            write(*,*) "Read wtsurf:  ", mi
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
          if (mi.ne.rf_swrad) then
C
            rf_swrad = mi
C
            filename = trim(pth_wrk)//trim(pth_flx)//
     $                 trim(pfx_dmn)//"roms_frc.nc"
            call check( nf90_open(filename, NF90_NOWRITE, ncid) )
            call check( nf90_inq_varid(ncid, "swrad", varid) )
            if (mi.ne.1) then
              call check(nf90_get_var(ncid, varid, swradb,(/1,1,mi-1/)))
              call check( nf90_get_var(ncid, varid, swradf,(/1,1,mi/)) )
            else
              call check( nf90_get_var(ncid, varid, swradb, (/1,1,12/)))
              call check( nf90_get_var(ncid, varid, swradf, (/1,1,1/)) )
            end if
            call check( nf90_close(ncid) )

            ! Convert swrad from W/m^2 to K m/s.
            swradb = swradb/hf_fac
C
            write(*,*) "Read swrad:   ", mi

          end if
!       Interpolate
          swrad = swradb+fac*(swradf-swradb)

          return
C
        case (5) ! momentum flux
C
          if (mi.ne.rf_wsurf) then
C
            rf_wsurf = mi
C
            filename = trim(pth_wrk)//trim(pth_flx)//
     $                 trim(pfx_dmn)//"roms_frc.nc"
            call check( nf90_open(filename, NF90_NOWRITE, ncid) )
!
            call check( nf90_inq_varid(ncid, "sustr", varid) )

            call check( nf90_get_var(ncid, varid, wusurf(2:im,:)
     $                   , (/1,1,mi/),   (/imm1,jm,1/)) )

            if (BC%ipl) then
              if (mi.ne.1) then
                call check( nf90_get_var(ncid, varid, wusurfb(2:im,:)
     $                   , (/1,1,mi-1/), (/imm1,jm,1/)) )
              else
                call check( nf90_get_var(ncid, varid, wusurfb(2:im,:)
     $                   , (/1,1,12/), (/imm1,jm,1/)) )
              end if
            end if

            call check( nf90_inq_varid(ncid, "svstr", varid) )

            call check( nf90_get_var(ncid, varid, wvsurf(:,2:jm)
     $                   , (/1,1,mi/),   (/im,jmm1,1/)) )

            if (BC%ipl) then
              if (mi.ne.1) then
                call check( nf90_get_var(ncid, varid, wvsurfb(:,2:jm)
     $                   , (/1,1,mi-1/), (/im,jmm1,1/)) )
              else
                call check( nf90_get_var(ncid, varid, wvsurfb(:,2:jm)
     $                   , (/1,1,12/), (/im,jmm1,1/)) )
              end if
            end if

            call check( nf90_close(ncid) )

            ! Fill in boundaries
            wusurf(1,:) = wusurf(2,:)
            wvsurf(:,1) = wvsurf(:,2)
            ! Taper wind stress along the boundary
            do i = 2, imm1
              do j = 2, jmm1
                wusurf(i,j) = .25*wusurf(i,j)
     $                       *(fsm(i,j+1)+fsm(i+1,j)
     $                        +fsm(i,j-1)+fsm(i-1,j)) ! TODO: or use fsm? Fsm may be preferrable for i=1, j=1
                wvsurf(i,j) = .25*wvsurf(i,j)
     $                       *(fsm(i,j+1)+fsm(i+1,j)
     $                        +fsm(i,j-1)+fsm(i-1,j)) !
              end do
            end do
            do i = 2, imm1
              wusurf(i, 1) = wusurf(i, 1)/3
     $                      *(fsm(i,   2)+fsm(i+1, 1)
     $                       +fsm(i-1, 1))
              wusurf(i,jm) = wusurf(i,jm)/3
     $                      *(fsm(i,jmm1)+fsm(i+1,jm)
     $                       +fsm(i-1,jm))
            end do
            do j = 2, jmm1
              wvsurf( 1,j) = wvsurf( 1,j)/3
     $                      *(fsm( 1,j+1)+fsm(   2,j)
     $                       +fsm( 1,j-1))
              wvsurf(im,j) = wvsurf(im,j)/3
     $                      *(fsm(im,j+1)+fsm(imm1,j)
     $                       +fsm(im,j-1))
            end do
            wusurf( 1, 1) = .5*wusurf( 1, 1)
     $                     *(fsm( 1,   2)+fsm(   2,   1))
            wusurf(im, 1) = .5*wusurf(im, 1)
     $                     *(fsm(im,   2)+fsm(imm1,   1))
            wusurf( 1,jm) = .5*wusurf( 1,jm)
     $                     *(fsm( 2,  jm)+fsm(   2,jmm1))
            wusurf(im,jm) = .5*wusurf(im,jm)
     $                     *(fsm(im,jmm1)+fsm(imm1,  jm))
            wvsurf( 1, 1) = .5*wvsurf( 1, 1)
     $                     *(fsm( 1,   2)+fsm(   2,   1))
            wvsurf(im, 1) = .5*wvsurf(im, 1)
     $                     *(fsm(im,   2)+fsm(imm1,   1))
            wvsurf( 1,jm) = .5*wvsurf( 1,jm)
     $                     *(fsm( 2,  jm)+fsm(   2,jmm1))
            wvsurf(im,jm) = .5*wvsurf(im,jm)
     $                     *(fsm(im,jmm1)+fsm(imm1,  jm))

!
            if (BC%ipl) then
                
              wusurff = wusurf
              wvsurff = wvsurf

              ! Taper for *b fields too.

              ! Fill in boundaries
              wusurfb(1,:) = wusurfb(2,:)
              wvsurfb(:,1) = wvsurfb(:,2)
              ! Taper wind stress along the boundary
              do i = 2, imm1
                do j = 2, jmm1
                  wusurfb(i,j) = .25*wusurfb(i,j)
     $                          *(fsm(i,j+1)+fsm(i+1,j)
     $                           +fsm(i,j-1)+fsm(i-1,j)) ! TODO: or use fsm?
                  wvsurfb(i,j) = .25*wvsurfb(i,j)
     $                          *(fsm(i,j+1)+fsm(i+1,j)
     $                           +fsm(i,j-1)+fsm(i-1,j)) !
                end do
              end do
              do i = 2, imm1
                wusurfb(i, 1) = wusurfb(i, 1)/3
     $                         *(fsm(i,   2)+fsm(i+1, 1)
     $                          +fsm(i-1, 1))
                wusurfb(i,jm) = wusurfb(i,jm)/3
     $                         *(fsm(i,jmm1)+fsm(i+1,jm)
     $                          +fsm(i-1,jm))
              end do
              do j = 2, jmm1
                wvsurfb( 1,j) = wvsurfb( 1,j)/3
     $                        *(fsm( 1,j+1)+fsm(   2,j)
     $                         +fsm( 1,j-1))
                wvsurfb(im,j) = wvsurfb(im,j)/3
     $                        *(fsm(im,j+1)+fsm(imm1,j)
     $                         +fsm(im,j-1))
              end do
              wusurfb( 1, 1) = .5*wusurfb( 1, 1)
     $                        *(fsm( 1,   2)+fsm(   2,   1))
              wusurfb(im, 1) = .5*wusurfb(im, 1)
     $                        *(fsm(im,   2)+fsm(imm1,   1))
              wusurfb( 1,jm) = .5*wusurfb( 1,jm)
     $                        *(fsm( 2,  jm)+fsm(   2,jmm1))
              wusurfb(im,jm) = .5*wusurfb(im,jm)
     $                        *(fsm(im,jmm1)+fsm(imm1,  jm))
              wvsurfb( 1, 1) = .5*wvsurfb( 1, 1)
     $                        *(fsm( 1,   2)+fsm(   2,   1))
              wvsurfb(im, 1) = .5*wvsurfb(im, 1)
     $                        *(fsm(im,   2)+fsm(imm1,   1))
              wvsurfb( 1,jm) = .5*wvsurfb( 1,jm)
     $                        *(fsm( 2,  jm)+fsm(   2,jmm1))
              wvsurfb(im,jm) = .5*wvsurfb(im,jm)
     $                        *(fsm(im,jmm1)+fsm(imm1,  jm))

              ! Convert s?str (in N/m^2) to w?surf (in m^2/s^2).
              wusurfb = -wusurfb/rhoref
              wvsurfb = -wvsurfb/rhoref
              wusurff = -wusurff/rhoref
              wvsurff = -wvsurff/rhoref
            else
              wusurf  = -wusurf/rhoref
              wvsurf  = -wvsurf/rhoref
            end if

            write(*,*) "Read w*surf:  ", mi

          end if

          if (BC%ipl) then
            wusurf = wusurfb+fac*(wusurff-wusurfb)
            wvsurf = wvsurfb+fac*(wvsurff-wvsurfb)
          end if

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
      subroutine bry_roms(idx)
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
      integer :: i,j,k, ncid,varid

      character(len=256) filename

      include 'pomNW.c'

C
      select case (idx)
!
          case (0)            ! Read rmean from pom grid provided file.

          if (mi.ne.rf_rmn) then

            rf_rmn = mi

            filename = trim(pth_wrk)//trim(pth_grd)
     $                 //trim(pfx_dmn)//"pom_clm.nc"
            write(*,*) "Reading BC: ", trim(filename)
            call check( nf90_open(filename, NF90_NOWRITE, ncid) )

            call check( nf90_inq_varid(ncid, "Rmean", varid) )
            call check( nf90_get_var(ncid, varid,
     $           rmean,(/1,1,1,mi/),(/im,jm,kbm1,1/)) )

            if (BC%ipl) then
              if (mi==1) then
                call check( nf90_get_var(ncid, varid,
     $                      rmeanb,(/1,1,1,12/),(/im,jm,kbm1,1/)) )
              else
                call check( nf90_get_var(ncid, varid,
     $                      rmeanb,(/1,1,1,mi-1/),(/im,jm,kbm1,1/)) )
              end if
              rmeanf = rmean
            end if

            call check( nf90_close(ncid) )

            do k=1,kbm1
              rmean(:,:,k) = rmean(:,:,k)*fsm(:,:)
            end do
            rmean(:,:,kb) = rmean(:,:,kbm1)

          end if

          if (BC%ipl) then
            rmean = rmeanb+fac*(rmeanf-rmeanb)      ! TODO: Not sure if it is correct to intepolate rmean since interpolated field won't correspond to interpolated tclim and sclim.
          end if

          return
!
        case (1) ! Not a boundary condition but climate

          if (mi.ne.rf_clm) then

            rf_clm = mi

            filename = trim(pth_wrk)//trim(pth_bry)//
     $                 trim(pfx_dmn)//"roms_clm.nc"
            call check( nf90_open(filename, NF90_NOWRITE, ncid) )
!
            call check( nf90_inq_varid(ncid, "temp", varid) )
            call check( nf90_get_var(ncid, varid,
     $           tclim(:,:,kbm1:1:-1),(/1,1,1,mi/),(/im,jm,kbm1,1/)) )

            if (BC%ipl) then

              if (mi.ne.1) then
                call check( nf90_get_var(ncid, varid,
     $               tclimb(:,:,kbm1:1:-1),(/1,1,mi-1/),(/im,jm,1/)) )
              else
                call check( nf90_get_var(ncid, varid,
     $               tclimb(:,:,kbm1:1:-1),(/1,1,  12/),(/im,jm,1/)) )
              end if

              tclimf = tclim

            end if
!
            call check( nf90_inq_varid(ncid, "salt", varid) )
            call check( nf90_get_var(ncid, varid,
     $           sclim(:,:,kbm1:1:-1),(/1,1,1,mi/),(/im,jm,kbm1,1/)) )

            if (BC%ipl) then

              if (mi.ne.1) then
                call check( nf90_get_var(ncid, varid,
     $               sclimb(:,:,kbm1:1:-1),(/1,1,mi-1/),(/im,jm,1/)) )
              else
                call check( nf90_get_var(ncid, varid,
     $               sclimb(:,:,kbm1:1:-1),(/1,1,  12/),(/im,jm,1/)) )
              end if

              sclimf = sclim

            end if

            call check( nf90_close(ncid) )
!
            do k=1,kbm1
                tclim(:,:,k) = tclim(:,:,k)*fsm(:,:)
                sclim(:,:,k) = sclim(:,:,k)*fsm(:,:)
            end do

            tclim(:,:,kb) = tclim(:,:,kbm1)
            sclim(:,:,kb) = sclim(:,:,kbm1)

            write(*,*) "Read climate: ", mi

          end if
!     Perform interpolation
          ! FSM is already defined! Feel free to apply it! (Does not yet applied!!!)
          if (BC%ipl) then
            tclim = (tclimb+fac*(tclimf-tclimb))
            sclim = (sclimb+fac*(sclimf-sclimb))
          end if

          return

        case (2) ! elevation
!
          if (mi.ne.rf_el) then
!     If we move to the next month...
            rf_el = mi

            filename = trim(pth_wrk)//trim(pth_bry)//
     $                 trim(pfx_dmn)//"roms_bry.nc"
            call check( nf90_open(filename, NF90_NOWRITE, ncid) )
            call check( nf90_inq_varid(ncid, "zeta_south", varid) )
!     ...read neccessary fields
!     South:
            call check( nf90_get_var(ncid, varid, els,
     $                            (/1,mi/), (/im,1/)) )

            if (BC%ipl) then

              if (mi.ne.1) then
                call check( nf90_get_var(ncid, varid, elsb,
     $                            (/1,mi-1/), (/im,1/)) )
              else
                call check( nf90_get_var(ncid, varid, elsb,
     $                            (/1,12/), (/im,1/)) )
              end if

              elsf = els

            end if

            call check( nf90_close(ncid) )

            write(*,*) "Read el BCs:  ", mi

          end if
!     Perform interpolation
!     South:
          if (BC%ipl) els(:) = (elsb(:)+fac*(elsf(:)-elsb(:)))*fsm(:,1)

          return
C
        case (3) ! u and v
C
          if (mi.ne.rf_uv) then
!        If we move to the next month...
            rf_uv = mi

            filename = trim(pth_wrk)//trim(pth_bry)//
     $                 trim(pfx_dmn)//"roms_bry.nc"
            call check( nf90_open(filename, NF90_NOWRITE, ncid) )
!        ...read neccessary fields
!     South:
            call check( nf90_inq_varid(ncid, "u_south", varid) )

            call check( nf90_get_var(ncid, varid,
     $           ubs(2:im,kbm1:1:-1), (/1,1,mi/), (/imm1,kbm1,1/)) )

            if (BC%ipl) then

              if (mi.ne.1) then
                call check( nf90_get_var(ncid, varid,
     $           ubsb(2:im,kbm1:1:-1),(/1,1,mi-1/),(/imm1,kbm1,1/)) )
              else
                call check( nf90_get_var(ncid, varid,
     $           ubsb(2:im,kbm1:1:-1),(/1,1,12/),(/imm1,kbm1,1/)) )
              end if

              ubsf = ubs

            end if

            call check( nf90_inq_varid(ncid, "v_south", varid) )

            call check( nf90_get_var(ncid, varid,
     $           vbs(:,kbm1:1:-1),(/1,1,mi/),(/im,kbm1,1/)) )

            if (BC%ipl) then

              if (mi.ne.1) then
                call check( nf90_get_var(ncid, varid,
     $           vbsb(:,kbm1:1:-1),(/1,1,mi-1/), (/im,kbm1,1/)) )
              else
                call check( nf90_get_var(ncid, varid,
     $           vbsb(:,kbm1:1:-1),(/1,1,12/),(/im,kbm1,1/)) )
              end if

              vbsf = vbs

            else

!     Calculate vertical means for 2D mode
              do i=1,im
                uabs(i) = 0.
                vabs(i) = 0.
                do k=1,kb
                  uabs(i) = uabs(i)+ubs(i,k)*dz(k)
                  vabs(i) = vabs(i)+vbs(i,k)*dz(k)
                end do
              end do

            end if

            call check( nf90_close(ncid) )

            write(*,*) "Read vel BCs: ", mi

          end if
!     Perform interpolation
!     South:
!          do k=1,kbm1   ! FSM is already defined! Feel free to apply it!
          if (BC%ipl) then

            ubs(:,1:kbm1) = ubsb(:,1:kbm1)
     $                     +fac*(ubsf(:,1:kbm1)-ubsb(:,1:kbm1))
            vbs(:,1:kbm1) = vbsb(:,1:kbm1)
     $                     +fac*(vbsf(:,1:kbm1)-vbsb(:,1:kbm1))
!          do k=1,kbm1              ! This doesn't work, since it completely wipes out u
!            ubs(:,k) = ubs(:,k)*dum(:,1)
!            vbs(:,k) = vbs(:,k)*dvm(:,2) ! index 1, maybe?
!          end do
!          end do

!     Calculate vertical mean BCs for 2D mode (for interpolated conditions; for stepped this code is performed once right after reading)
!
            do i=1,im
              uabs(i) = 0.
              vabs(i) = 0.
              do k=1,kb
                uabs(i) = uabs(i)+ubs(i,k)*dz(k)
                vabs(i) = vabs(i)+vbs(i,k)*dz(k)
              end do
            end do

          end if

          return
C
        case (4) ! temperature and salinity
C
          if (mi.ne.rf_ts) then
!        If we move to the next month...
            rf_ts = mi

            filename = trim(pth_wrk)//trim(pth_bry)//
     $                 trim(pfx_dmn)//"roms_clm.nc"
            call check( nf90_open(filename, NF90_NOWRITE, ncid) )
            call check( nf90_inq_varid(ncid, "temp", varid) )
!    Read climatological BCs
            call check( nf90_get_var(ncid, varid,
     $           tbs(:,kbm1:1:-1),(/ 1, 1,1,mi/), (/im,1,kbm1,1/)) )
            call check( nf90_get_var(ncid, varid,
     $           tbn(:,kbm1:1:-1),(/ 1,jm,1,mi/), (/im,1,kbm1,1/)) )
            call check( nf90_get_var(ncid, varid,
     $           tbw(:,kbm1:1:-1),(/ 1, 1,1,mi/), (/1,jm,kbm1,1/)) )
            call check( nf90_get_var(ncid, varid,
     $           tbe(:,kbm1:1:-1),(/im, 1,1,mi/), (/1,jm,kbm1,1/)) )

            if (BC%ipl) then

              if (mi.ne.1) then
!            south
                call check( nf90_get_var(ncid, varid,
     $           tbsb(:,kbm1:1:-1),(/1,1,1,mi-1/), (/im,1,kbm1,1/)) )
!            north
                call check( nf90_get_var(ncid, varid,
     $           tbnb(:,kbm1:1:-1),(/1,jm,1,mi-1/), (/im,1,kbm1,1/)) )
!            west
                call check( nf90_get_var(ncid, varid,
     $           tbwb(:,kbm1:1:-1),(/1,1,1,mi-1/), (/1,jm,kbm1,1/)) )
!            east
                call check( nf90_get_var(ncid, varid,
     $           tbeb(:,kbm1:1:-1),(/im,1,1,mi-1/), (/1,jm,kbm1,1/)) )
              else
!            south
                call check( nf90_get_var(ncid, varid,
     $           tbsb(:,kbm1:1:-1),(/1,1,1,12/),(/im,1,kbm1,1/)) )
!            north
                call check( nf90_get_var(ncid, varid,
     $           tbnb(:,kbm1:1:-1),(/1,jm,1,12/),(/im,1,kbm1,1/)) )
!            west
                call check( nf90_get_var(ncid, varid,
     $           tbwb(:,kbm1:1:-1),(/1,1,1,12/),(/1,jm,kbm1,1/)) )
!            east
                call check( nf90_get_var(ncid, varid,
     $           tbeb(:,kbm1:1:-1),(/im,1,1,12/),(/1,jm,kbm1,1/)) )
              end if

              tbsf = tbs
              tbnf = tbn
              tbwf = tbw
              tbef = tbe

            end if
!
            call check( nf90_inq_varid(ncid, "salt", varid) )
!    Read climatological salinity BCs
            call check( nf90_get_var(ncid, varid,
     $           sbs(:,kbm1:1:-1),(/ 1, 1,1,mi/), (/im,1,kbm1,1/)) )
            call check( nf90_get_var(ncid, varid,
     $           sbn(:,kbm1:1:-1),(/ 1,jm,1,mi/), (/im,1,kbm1,1/)) )
            call check( nf90_get_var(ncid, varid,
     $           sbw(:,kbm1:1:-1),(/ 1, 1,1,mi/), (/1,jm,kbm1,1/)) )
            call check( nf90_get_var(ncid, varid,
     $           sbe(:,kbm1:1:-1),(/im, 1,1,mi/), (/1,jm,kbm1,1/)) )

            if (BC%ipl) then

              if (mi.ne.1) then
!            south
                call check( nf90_get_var(ncid, varid,
     $           sbsb(:,kbm1:1:-1),(/1, 1,1,mi-1/), (/im,1,kbm1,1/)) )
!            north
                call check( nf90_get_var(ncid, varid,
     $           sbnb(:,kbm1:1:-1),(/1,jm,1,mi-1/), (/im,1,kbm1,1/)) )
!            west
                call check( nf90_get_var(ncid, varid,
     $           sbwb(:,kbm1:1:-1),(/ 1,1,1,mi-1/), (/1,jm,kbm1,1/)) )
!            east
                call check( nf90_get_var(ncid, varid,
     $           sbeb(:,kbm1:1:-1),(/im,1,1,mi-1/), (/1,jm,kbm1,1/)) )
              else
!            south
                call check( nf90_get_var(ncid, varid,
     $           sbsb(:,kbm1:1:-1),(/1, 1,1,12/), (/im,1,kbm1,1/)) )
!            north
                call check( nf90_get_var(ncid, varid,
     $           sbnb(:,kbm1:1:-1),(/1,jm,1,12/), (/im,1,kbm1,1/)) )
!            west
                call check( nf90_get_var(ncid, varid,
     $           sbwb(:,kbm1:1:-1),(/ 1,1,1,12/), (/1,jm,kbm1,1/)) )
!            east
                call check( nf90_get_var(ncid, varid,
     $           sbeb(:,kbm1:1:-1),(/im,1,1,12/), (/1,jm,kbm1,1/)) )
              end if

              sbsf = sbs
              sbnf = sbn
              sbwf = sbw
              sbef = sbe

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
            if (mi.ne.1) then
              call check( nf90_get_var(ncid, varid,
     $         tbsb(:,kbm1:1:-1),(/1,1,mi-1/), (/im,kbm1,1/)) )
              call check( nf90_get_var(ncid, varid,
     $         tbsf(:,kbm1:1:-1),(/1,1,mi/),   (/im,kbm1,1/)) )
            else
              call check( nf90_get_var(ncid, varid,
     $         tbsb(:,kbm1:1:-1),(/1,1,12/),(/im,kbm1,1/)) )
              call check( nf90_get_var(ncid, varid,
     $         tbsf(:,kbm1:1:-1),(/1,1,1/), (/im,kbm1,1/)) )
            end if
            call check( nf90_inq_varid(ncid, "salt_south", varid) )
            if (mi.ne.1) then
              call check( nf90_get_var(ncid, varid,
     $         sbsb(:,kbm1:1:-1),(/1,1,mi-1/), (/im,kbm1,1/)) )
              call check( nf90_get_var(ncid, varid,
     $         sbsf(:,kbm1:1:-1),(/1,1,mi/),   (/im,kbm1,1/)) )
            else
              call check( nf90_get_var(ncid, varid,
     $         sbsb(:,kbm1:1:-1),(/1,1,12/),(/im,kbm1,1/)) )
              call check( nf90_get_var(ncid, varid,
     $         sbsf(:,kbm1:1:-1),(/1,1,1/), (/im,kbm1,1/)) )
            end if

            call check( nf90_close(ncid) )
            end if

            write(*,*) "Read TS BCs:  ", mi

          end if

          if (BC%ipl) then
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
          end if

          return
C
          case (43)   ! Surface temperature and salinity from ICOADS

          if (mi.ne.rf_sts) then
!        If we move to the next month...
            rf_sts = mi

            filename = trim(pth_wrk)//trim(pth_bry)//
     $                 trim(pfx_dmn)//"roms_frc.nc"
            call check( nf90_open(filename, NF90_NOWRITE, ncid) )
!        ...read neccessary fields
!     South:
            call check( nf90_inq_varid(ncid, "SST", varid) )

            call check( nf90_get_var(ncid, varid,
     $           tsurf,(/1,1,mi/), (/im,jm,1/)) )

            if (BC%ipl) then

              if (mi.ne.1) then
                call check( nf90_get_var(ncid, varid,
     $           tsurfb,(/1,1,mi-1/), (/im,jm,1/)) )
              else
                call check( nf90_get_var(ncid, varid,
     $           tsurfb,(/1,1,12/), (/im,jm,1/)) )
              end if

              tsurff = tsurf

            end if
!
            call check( nf90_inq_varid(ncid, "SSS", varid) )

            call check( nf90_get_var(ncid, varid,
     $           ssurf,(/1,1,mi/), (/im,jm,1/)) )

            if (BC%ipl) then

              if (mi.ne.1) then
                call check( nf90_get_var(ncid, varid,
     $           ssurfb,(/1,1,mi-1/), (/im,jm,1/)) )
              else
                call check( nf90_get_var(ncid, varid,
     $           ssurfb,(/1,1,12/), (/im,jm,1/)) )
              end if

              ssurff = ssurf

            end if

            call check( nf90_close(ncid) )

            write(*,*) "Read ST/S BCs:", mi

          end if
!
          if (BC%ipl) then
!     Perform interpolation
!     South:
!          do k=1,kbm1   ! FSM is already defined! Feel free to apply it!
            tsurf = (tsurfb+fac*(tsurff-tsurfb))*fsm
            ssurf = (ssurfb+fac*(ssurff-ssurfb))*fsm
!          end do
          end if

          return
!
          case (45)   ! Temperature and salinity of surface level from WOA

          if (mi.ne.rf_sts) then
!        If we move to the next month...
            rf_sts = mi

            filename = trim(pth_wrk)//trim(pth_bry)//
     $                 trim(pfx_dmn)//"roms_clm.nc"
            call check(nf90_open(filename, NF90_NOWRITE, ncid))
!        ...read neccessary fields
!     South:
            call check( nf90_inq_varid(ncid, "temp", varid) )

            call check( nf90_get_var(ncid, varid, tsurf,
     $           (/1,1,kbm1,mi/), (/im,jm,1,1/)) )

            if (BC%ipl) then

              if (mi.ne.1) then
                call check( nf90_get_var(ncid, varid, tsurfb,
     $           (/1,1,kbm1,mi-1/), (/im,jm,1,1/)) )
              else
                call check( nf90_get_var(ncid, varid, tsurfb,
     $           (/1,1,kbm1,  12/), (/im,jm,1,1/)) )
              end if

              tsurff = tsurf

            end if
!
            call check( nf90_inq_varid(ncid, "salt", varid) )

            call check( nf90_get_var(ncid, varid, ssurf,
     $           (/1,1,kbm1,mi/), (/im,jm,1,1/)) )

            if (BC%ipl) then

              if (mi.ne.1) then
                call check( nf90_get_var(ncid, varid,
     $           ssurfb,(/1,1,kbm1,mi-1/), (/im,jm,1,1/)) )
              else
                call check( nf90_get_var(ncid, varid,
     $           ssurfb,(/1,1,kbm1,  12/), (/im,jm,1,1/)) )
              end if

              ssurff = ssurf

            end if

            call check( nf90_close(ncid) )

            write(*,*) "Read 1l TS BCs:", mi

          end if

          if (BC%ipl) then
!     Perform interpolation
!     South:
!          do k=1,kbm1   ! FSM is already defined! Feel free to apply it!
            tsurf = (tsurfb+fac*(tsurff-tsurfb))*fsm
            ssurf = (ssurfb+fac*(ssurff-ssurfb))*fsm
!          end do
          end if
!
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
      end subroutine bry_roms
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
      integer :: i,j,k, ncid,varid, dist_num,p
      double precision :: trans_tot, trans_dist

      character(len=256) filename

      include 'pomNW.c'

C
      select case (idx)
!
          case (0)            ! Read rmean from pom grid provided file. TODO: Maybe it should be better to move this to case 1 and read all rmean, tclim and sclim at the same time?

          if (mi /= rf_rmn) then

            rf_rmn = mi

            filename = trim(pth_wrk)//trim(pth_grd)
     $                 //trim(pfx_dmn)//"pom_clm.nc"
            call check( nf90_open(filename, NF90_NOWRITE, ncid) )

            call check( nf90_inq_varid(ncid, "Rmean", varid) )
            call check( nf90_get_var(ncid, varid,
     $           rmean,(/1,1,1,mi/),(/im,jm,kb,1/)) )

            if (BC%ipl) then
              if (mi==1) then
                call check( nf90_get_var(ncid, varid,
     $                      rmeanb,(/1,1,1,12/),(/im,jm,kb,1/)) )
              else
                call check( nf90_get_var(ncid, varid,
     $                      rmeanb,(/1,1,1,mi-1/),(/im,jm,kb,1/)) )
              end if
              rmeanf = rmean
            end if

            call check( nf90_close(ncid) )

            do k=1,kbm1
              rmean(:,:,k) = rmean(:,:,k)*fsm(:,:)
            end do
            
            write(*,*) "[-] Read background density: ", mi

          end if

          if (BC%ipl) then
            rmean = rmeanb+fac*(rmeanf-rmeanb)      ! TODO: Not sure if it is correct to intepolate rmean since interpolated field won't correspond to interpolated tmean and smean.
          end if
          

          return
!
        case (1) ! Not a boundary condition but climate

          if (mi /= rf_clm) then

            rf_clm = mi

            filename = trim(pth_wrk)//trim(pth_bry)//
     $                 trim(pfx_dmn)//"pom_clm.nc"
            call check( nf90_open(filename, NF90_NOWRITE, ncid) )
!
            call check( nf90_inq_varid(ncid, "Tclim", varid) )
            call check( nf90_get_var(ncid, varid,
     $           tclim, (/1,1,1,mi/), (/im,jm,kb,1/)) )

            if (BC%ipl) then

              if (mi == 1) then
                call check( nf90_get_var(ncid, varid,
     $               tclimb,(/1,1,1,  12/), (/im,jm,kb,1/)) )
              else
                call check( nf90_get_var(ncid, varid,
     $               tclimb,(/1,1,1,mi-1/), (/im,jm,kb,1/)) )
              end if
              tclimf = tclim

            end if
!
            call check( nf90_inq_varid(ncid, "Sclim", varid) )
            call check( nf90_get_var(ncid, varid,
     $           sclim, (/1,1,1,mi/), (/im,jm,kb,1/)) )

            if (BC%ipl) then

              if (mi == 1) then
                call check( nf90_get_var(ncid, varid,
     $               sclimb, (/1,1,1,  12/), (/im,jm,kb,1/)) )
              else
                call check( nf90_get_var(ncid, varid,
     $               sclimb, (/1,1,1,mi-1/), (/im,jm,kb,1/)) )
              end if
              sclimf = sclim

            end if

            call check( nf90_close(ncid) )
!
            do k=1,kb
              tclim(:,:,k) = tclim(:,:,k)*fsm(:,:)
              sclim(:,:,k) = sclim(:,:,k)*fsm(:,:)
            end do

            write(*,*) "[-] Read climate:   ", mi
!     Update surface T and S
            tsurf = tclim(:,:,1)
            ssurf = sclim(:,:,1)
            write(*,*) "[-] Got surface TS: ", mi

          end if
!     Perform interpolation
          if (BC%ipl) then
            tclim = (tclimb+fac*(tclimf-tclimb))
            sclim = (sclimb+fac*(sclimf-sclimb))
!     Update surface T and S
            tsurf = tclim(:,:,1)
            ssurf = sclim(:,:,1)
          end if

          return

        case (2) ! elevation (TODO: implement elevation)
!
          if (mi.ne.rf_el) then
!     If we move to the next month...
            rf_el = mi
            
            
            filename = trim(pth_wrk)//trim(pth_grd)
     $                 //trim(pfx_dmn)//"pom_bry.nc"
            call check( nf90_open(filename, NF90_NOWRITE, ncid) )

            if (BC%bnd%nth) then
              call check( nf90_inq_varid(ncid, "north.el", varid) )
              call check( nf90_get_var(ncid, varid,
     $                    eln,(/1,mi/),(/im,1/)) )
            end if
            if (BC%bnd%est) then
              call check( nf90_inq_varid(ncid, "east.el",  varid) )
              call check( nf90_get_var(ncid, varid,
     $                    ele,(/1,mi/),(/jm,1/)) )
            end if
            if (BC%bnd%sth) then
              call check( nf90_inq_varid(ncid, "south.el", varid) )
              call check( nf90_get_var(ncid, varid,
     $                    els,(/1,mi/),(/im,1/)) )
            end if
            if (BC%bnd%wst) then
              call check( nf90_inq_varid(ncid, "west.el", varid) )
              call check( nf90_get_var(ncid, varid,
     $                    elw,(/1,mi/),(/im,1/)) )
            end if

            call check( nf90_close(ncid) )

            write(*,*) "[-] Elevation read. ", mi
            
          end if

          
          return
!
        case (3) ! u and v
!
          if (mi.ne.rf_uv) then
!        If we move to the next month...
            rf_uv = mi
            
            write(*,*) "[@] TODO: implement variable BCs for currents."

            filename = trim(pth_wrk)//trim(pth_grd)
     $                 //trim(pfx_dmn)//"pom_bry.nc"
            call check( nf90_open(filename, NF90_NOWRITE, ncid) )

            if (BC%bnd%nth) then
              call check( nf90_inq_varid(ncid, "north.v", varid) )
              call check( nf90_get_var(ncid, varid,
     $                    vbn,(/1,1,mi/),(/im,kb,1/)) )
              vabn = 0.
              do k=1,kbm1
                vabn(:) = vabn(:) + dz(k)*vbn(:,k)
              end do
            end if
            if (BC%bnd%est) then
              call check( nf90_inq_varid(ncid, "east.u",  varid) )
              call check( nf90_get_var(ncid, varid,
     $                    ube,(/1,1,mi/),(/jm,kb,1/)) )
              uabe = 0.
              do k=1,kbm1
                uabe(:) = uabe(:) + dz(k)*ube(:,k)
              end do
            end if
            if (BC%bnd%sth) then
              call check( nf90_inq_varid(ncid, "south.v", varid) )
              call check( nf90_get_var(ncid, varid,
     $                    vbs,(/1,1,mi/),(/im,kb,1/)) )
              vabs = 0.
              do k=1,kbm1
                vabs(:) = vabs(:) + dz(k)*vbs(:,k)
              end do
            end if
            if (BC%bnd%wst) then
              call check( nf90_inq_varid(ncid, "west.u", varid) )
              call check( nf90_get_var(ncid, varid,
     $                    ube,(/1,1,mi/),(/im,kb,1/)) )
              uabe = 0.
              do k=1,kbm1
                uabe(:) = uabe(:) + dz(k)*ube(:,k)
              end do
            end if

            call check( nf90_close(ncid) )
            
            write(*,*) "[-] Read boundary velocities"
            
            if (BC%bnd%vol) call volConserve

          end if

          return
!
        case (4) ! vertical TS
!
!     Get boundary TS from Tclim and Sclim w/o reading the `clm` file for now.
          tbn = tclim(:,jm,:)
          sbn = sclim(:,jm,:)
          tbe = tclim(im,:,:)
          sbe = sclim(im,:,:)
          tbs = tclim(:, 1,:)
          sbs = sclim(:, 1,:)
          tbw = tclim( 1,:,:)
          sbw = sclim( 1,:,:)

          return
!
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
      end subroutine bry
!
      subroutine upd_mnth(tmp, ipl)

        include 'pomNW.c'

        logical,intent( in) :: ipl
        double precision                :: tind, tmp, b, e
!
        tind = mod(tmp,365.)
!
!       If interpolation is enabled we slice months up at their middles.
!
!       Climatological cycle doesn't work for ipl yet.
        if (ipl) then

          if (tind.lt.15) then
            mi = 1
            b = -15.5
            e = 15.
          else
            if (tind.lt.44) then
              mi = 2
              b = 15.
              e = 44.
            else
              if (tind.lt.73.5) then
                mi = 3
                b = 44.
                e = 73.5
              else
                if (tind.lt.104) then
                  mi = 4
                  b = 73.5
                  e = 104.
                else
                  if (tind.lt.134.5) then
                    mi = 5
                    b = 104.
                    e = 134.5
                  else
                    if (tind.lt.165) then
                      mi = 6
                      b = 134.5
                      e = 165.
                    else
                      if (tind.lt.195.5) then
                        mi = 7
                        b = 165.
                        e = 195.5
                      else
                        if (tind.lt.226.5) then
                          mi = 8
                          b = 195.5
                          e = 226.5
                        else
                          if (tind.lt.258) then
                            mi = 9
                            b = 226.5
                            e = 258.
                          else
                            if (tind.lt.288.5) then
                              mi = 10
                              b = 258.
                              e = 288.5
                            else
                              if (tind.lt.319) then
                                mi = 11
                                b = 288.5
                                e = 319.
                              else
                                if (tind.lt.349.5) then
                                  mi = 12
                                  b = 319.
                                  e = 349.5
                                else
                                  mi = 1
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

        else
!
!     If interpolation is disabled we just detect which month the day is in.
!
          if (tind.lt.31) then
            mi = 1
          else
            if (tind.lt.59) then
              mi = 2
            else
              if (tind.lt.90) then
                mi = 3
              else
                if (tind.lt.120) then
                  mi = 4
                else
                  if (tind.lt.151) then
                    mi = 5
                  else
                    if (tind.lt.181) then
                      mi = 6
                    else
                      if (tind.lt.212) then
                        mi = 7
                      else
                        if (tind.lt.243) then
                          mi = 8
                        else
                          if (tind.lt.273) then
                            mi = 9
                          else
                            if (tind.lt.304) then
                              mi = 10
                            else
                              if (tind.lt.334) then
                                mi = 11
                              else
                                if (tind.lt.365) then
                                  mi = 12
                                else
                                  ! How did you get here?
                                  mi = 1
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
        end if
!
        return

      end
!
      subroutine clm_warp()
!
        include 'pomNW.c'
!
!       If climatological cycle is less than 12 we need to warp.    !TODO: make cycling in days. It should be easier and more flexible.
!
        if (clm_cycle<12) then
          mi = modulo(mi-m0, clm_cycle)+m0
        end if
!
      end subroutine
!

      subroutine time2date(time_in, time_off, date)

        character(len=*), intent(in) :: time_off
        double precision                         :: time_in,  tmp
        character*26,     intent(out):: date
        integer  :: YYYY, MM, DD, hh, mi, ss, th, tm
        character :: sign

        read(time_off,
     $       '(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2,1x,a1,i2,1x,i2)')
     $       YYYY, MM, DD, hh, mi, ss, sign, th, tm

        tmp  = time_in - float(floor(time_in))
        ss   = ss + mod(floor(tmp*60.*60.*24.), 60)
        mi   = mi + floor(ss/60.)
        ss   = mod(ss, 60)
        mi   = mi + mod(floor(tmp*60.*24.), 60)
        hh   = hh + floor(mi/60.)
        mi   = mod(mi, 60)
        hh   = hh + mod(floor(tmp*24.), 24)
        DD   = DD + floor(hh/24.)
        hh   = mod(hh, 24)

        tmp  = floor(time_in)
        YYYY = YYYY + floor((DD-1+tmp)/365.)
        DD   = mod((DD+tmp),365.)

        if (DD==0) then
          MM = 12
          DD = 31
        else
          if (DD.le.31) then
            MM = 1
          else
            if (DD.le.59) then
              MM = 2
              DD = DD-31
            else
              if (DD.le.90) then
                MM = 3
                DD = DD-59
              else
                if (DD.le.120) then
                  MM = 4
                  DD = DD-90
                else
                  if (DD.le.151) then
                    MM = 5
                    DD = DD-120
                  else
                    if (DD.le.181) then
                      MM = 6
                      DD = DD-151
                    else
                      if (DD.le.212) then
                        MM = 7
                        DD = DD-181
                      else
                        if (DD.le.243) then
                          MM = 8
                          DD = DD-212
                        else
                          if (DD.le.273) then
                            MM = 9
                            DD = DD-243
                          else
                            if (DD.le.304) then
                              MM = 10
                              DD = DD-273
                            else
                              if (DD.le.334) then
                                MM = 11
                                DD = DD-304
                              else
                                if (DD.le.365) then
                                  MM = 12
                                  DD = DD-334
                                else
                                  ! How did you get here?
                                  MM = 1
                                  DD = DD-365
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
        end if
        !DD = DD+1

        write(date,
     $        '(i4,"-",i2.2,"-",i2.2," ",i2.2,":"
     $         ,i2.2,":",i2.2," ",a1,i2.2,":",i2.2)')
     $        YYYY, MM, DD, hh, mi, ss, sign, th, tm
        !write(*, *) date
        return

      end
      
      subroutine volConserve
        
        implicit none
        
        include 'pomNW.c'
        
        double precision :: wn(im,kb), we(jm,kb),
     $                      ws(im,kb), ww(jm,kb),
     $                      an(im,kb), ae(jm,kb),
     $                      as(im,kb), aw(jm,kb),
     $                      trn, tre, trs, trw  ! total boundary section transport
        double precision :: trans_tot, trans_discr
        integer :: i,j,k,p
        
        write(*,*) "[-] Volume conservation:"
        
! Calculate areas and weights
        do k=1,kbm1
          an(:,k) = dx(:,jm)*(h(:,jm)+h(:,jmm1))*.5*dz(k)
          ae(:,k) = dy(im,:)*(h(im,:)+h(imm1,:))*.5*dz(k)
          as(:,k) = dx(:, 1)*h(:,1)*dz(k)
          aw(:,k) = dy( 1,:)*h(1,:)*dz(k)
        end do
        an(:,kb) = 0.
        ae(:,kb) = 0.
        as(:,kb) = 0.
        aw(:,kb) = 0.
        
        do k=1,kbm1
          wn(:,k) = dvm(:,jm)*an(:,k)*vbn(:,k)
          we(:,k) = dum(im,:)*ae(:,k)*ube(:,k)
          ws(:,k) = dvm(:, 1)*as(:,k)*vbs(:,k)
          ww(:,k) = dum( 1,:)*aw(:,k)*ubw(:,k)
        end do
        wn(:,kb) = 0.
        we(:,kb) = 0.
        ws(:,kb) = 0.
        ww(:,kb) = 0.
        
        trn = sum(wn)
        tre = sum(we)
        trs = sum(ws)
        trw = sum(ww)
        
        trans_tot   = trn+tre+trs+trw
        trans_discr = trs-trn+trw-tre
        write(*,*) "[ ] Transport discrepancy before: ", trans_discr
        
        wn = wn/trans_tot
        we = we/trans_tot
        ws = ws/trans_tot
        ww = ww/trans_tot
        
        vbn(:,1:kbm1) = vbn(:,1:kbm1)
     $                 +wn(:,1:kbm1)*trans_discr/an(:,1:kbm1)
        vbs(:,1:kbm1) = vbs(:,1:kbm1)
     $                 -ws(:,1:kbm1)*trans_discr/as(:,1:kbm1)
        ube(:,1:kbm1) = ube(:,1:kbm1)
     $                 +we(:,1:kbm1)*trans_discr/ae(:,1:kbm1)
        ubw(:,1:kbm1) = ubw(:,1:kbm1)
     $                 -ww(:,1:kbm1)*trans_discr/aw(:,1:kbm1)
        vabn = 0.
        vabs = 0.
        uabw = 0.
        uabe = 0.
        do k=1,kbm1
          vabn = vabn + vbn(:,k)*dz(k)
          vabs = vabs + vbs(:,k)*dz(k)
          uabw = uabw + ubw(:,k)*dz(k)
          uabe = uabe + ube(:,k)*dz(k)
        end do
        
! Check transport discrepancy
        do k=1,kbm1
          wn(:,k) = dvm(:,jm)*an(:,k)*vbn(:,k)
          we(:,k) = dum(im,:)*ae(:,k)*ube(:,k)
          ws(:,k) = dvm(:, 1)*as(:,k)*vbs(:,k)
          ww(:,k) = dum( 1,:)*aw(:,k)*ubw(:,k)
        end do
        wn(:,kb) = 0.
        we(:,kb) = 0.
        ws(:,kb) = 0.
        ww(:,kb) = 0.
        
        trn = sum(wn)
        tre = sum(we)
        trs = sum(ws)
        trw = sum(ww)
        
        trans_tot   = trn+tre+trs+trw
        trans_discr = trs-trn+trw-tre
        write(*,*) "[ ] Transport discrepancy after: ", trans_discr
              
      end subroutine

      subroutine ncclose(ncid)
        use netcdf
        integer, intent(in) :: ncid
        integer varid
        
        call check( nf90_redef(ncid) )
        
        call check( nf90_inq_varid(ncid, "T", varid) )
        call check( nf90_put_att(ncid, varid, "_FillValue", 0.) )
        call check( nf90_inq_varid(ncid, "S", varid) )
        call check( nf90_put_att(ncid, varid, "_FillValue", 0.) )
        call check( nf90_inq_varid(ncid, "R", varid) )
        call check( nf90_put_att(ncid, varid, "_FillValue", 0.) )
!        call check( nf90_inq_varid(ncid, "Rmean", varid) )
!        call check( nf90_put_att(ncid, varid, "_FillValue", 0.) )
        
        call check( nf90_enddef(ncid) )
        
        call check( nf90_close(ncid) )
        
        contains
          subroutine check(status)
            integer, intent ( in) :: status
!            if (DBG) write(*,*) status
            if(status /= nf90_noerr) then
              stop "Stopped"
            end if
          end subroutine check
          
      end subroutine ncclose
!
!      include 'pom2k.n'                                       ! *netCDF*
C
C     End of source code
C
C-----------------------------------------------------------------------