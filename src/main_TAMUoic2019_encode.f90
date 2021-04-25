!+
! TAMUoic2019 
! 
!   Encode a single-scattering property databse of oriented ice crystals
!   based on the TAMUoic2019 data kernel
!
!   History: 
!     Masa Saito  07/10/2019 Released as v1.0.0.
!     Masa Saito  04/24/2021 Encoder and kernel technique implemented.
!                            Released as v1.1.0
! 
!
!   Reference:
!   Saito, M., and P. Yang, Oriented ice crystals: A single-scattering property
!        database for applications to lidar and optical phenomenon simulations, 
!        J. Atmos. Sci., 76, 2635–2652.
!
!-
program tamuoic2019encode

  use globals
  use hparx, only : open_dir, gridIdx_loc, open_seq, getCmdArgs, err_read
  implicit none

  ! TAMUoic2019 binary kernel table
  integer, parameter :: IUI        = 101   ! file ID
  integer, parameter :: Tab_nsiz   = 165   ! # of size bin
  integer, parameter :: Tab_nwav   = 3     ! # of wavelength bin
  integer, parameter :: Tab_nhab   = 2     ! # of ice crystal habit bin
  integer, parameter :: Tab_npang  = 498   ! # of scattering polar angle bin
  integer, parameter :: Tab_naang  = 180   ! # of scattering azimuth angle bin
  integer, parameter :: Tab_niang  = 94    ! # of incident polar angle bin
  integer, parameter :: Tab_nmet   = 6     ! # of meta information
  integer, save      :: Tab_npol           ! # of polarization elements
  integer            :: Tab_recls, Tab_reclp
  real(R_)           :: Tab_siz(Tab_nsiz),   Tab_wav(Tab_nwav)
  character(len=14)  :: Tab_hab(Tab_nhab)
  real(R_)           :: Tab_pang(Tab_npang), Tab_aang(Tab_naang), Tab_iang(Tab_niang) 

  ! TAMUoic2019 ASCII database
  integer, parameter    :: IUO       = 701 ! file ID
  character(len=256)    :: Oic_outpath     ! output dust database path
  integer, save         :: Oic_nsiz        ! # of size bin 
  integer, save         :: Oic_nwav        ! # of wavelength
  integer, parameter    :: Oic_npang = 498 ! # of scattering polar angle bin
  integer, parameter    :: Oic_naang = 180 ! # of scattering azimuth angle bin
  integer, save         :: Oic_niang       ! # of incident polar angle bin
  integer, save         :: Oic_npol        ! (Full phase matrix elements: P11, ..., P44)
  real(R_), allocatable :: Oic_siz(:)      ! (nsiz) maximum dimension
  real(R_), allocatable :: Oic_wav(:)      ! (nwav) wavelength
  real(R_), allocatable :: Oic_pang(:)     ! (npang) scattering polar angle
  real(R_), allocatable :: Oic_aang(:)     ! (naang) scattering azimuth angle
  real(R_), allocatable :: Oic_iang(:)     ! (niang) incident polar angle
  character(len=14)     :: Oic_hab

  ! TAMUoic2019 Core Process
  call TAMUoic2019__init()
  call TAMUoic2019__encode()
  call TAMUoic2019__final()

stop
contains


   !+
   ! TAMUoic2019: Read namelist and initialization
   !-
   subroutine TAMUoic2019__init()
     
      integer            :: ios, ivar, narg, mwarn
      character(len=256) :: argv(10) ! argument to be red
      character(len=256) :: argmsg   ! argment message
      character(len=256), save :: outpath            ! output path of horizontally oriented phase matrix
      character(len=256)       :: TAMUoic2019path    ! TAMUoic2019 binary table path
      integer,  save           :: nsiz, nwav, niang  ! # of size/wavelength bins
      real(R_), save           :: siz(Tab_nsiz), iang(Tab_niang), wav(Tab_nwav)
      character(len=14), save  :: hab = 'OrientedPlate '
      character(len=6),  save  :: pol = 'Vector'
      namelist /TAMUoic2019/TAMUoic2019path, outpath, &
                            nsiz, nwav, niang, siz, wav, iang, hab, pol  

      ! Read namelist
      argmsg = ' Usage: tamuoic2019encode nmlst > output '  
      narg = 1
      call getCmdArgs(narg, argv, argmsg)
      open(IUI, file = trim(adjustl(argv(1))), status = 'old')
      read(IUI, nml=TAMUoic2019, iostat=ios)
      call err_read(ios, IUI, 'Main: Error in reading the namelist data.')
      close(IUI)

 
      ! Initialization
      Oic_outpath = outpath
      Oic_nsiz    = nsiz
      Oic_nwav    = nwav
      Oic_niang   = niang
      Oic_hab     = hab
      if      (pol == 'Scaler') then
         Tab_npol = 1
      else if (pol == 'Vector') then
         Tab_npol = 16
      else
         stop 'ERROR: "pol" should be either "Scaler" or "Vector".' 
      endif
      Oic_npol = Tab_npol
      call TAMUoic2019__setDim()
      allocate( Oic_siz(Oic_nsiz), Oic_wav(Oic_nwav), Oic_iang(Oic_niang) )
      allocate( Oic_pang(Oic_npang), Oic_aang(Oic_naang)                  )
      Oic_pang(:) = Tab_pang(:)
      Oic_aang(:) = Tab_aang(:)
      Oic_iang(:) = iang(1:niang)
      Oic_siz(:)  = siz(1:nsiz)
      Oic_wav(:)  = wav(1:nwav)


      ! Inquire and load database
      inquire (iolength=Tab_recls) (1.0_R4_, ivar=1, Tab_nmet)
      inquire (iolength=Tab_reclp) (1.0_R4_, ivar=1, Tab_npang*Tab_naang)
      call open_dir(IUI,    trim(TAMUoic2019path)//'/TAMUoic2019isca.bin', Tab_recls, 'old')
      call open_dir(IUI+1,  trim(TAMUoic2019path)//'/TAMUoic2019P11.bin',  Tab_reclp, 'old')
      if (Oic_npol > 1) then
         call open_dir(IUI+2,  trim(TAMUoic2019path)//'/TAMUoic2019P12.bin', Tab_reclp, 'old')
         call open_dir(IUI+3,  trim(TAMUoic2019path)//'/TAMUoic2019P13.bin', Tab_reclp, 'old')
         call open_dir(IUI+4,  trim(TAMUoic2019path)//'/TAMUoic2019P14.bin', Tab_reclp, 'old')
         call open_dir(IUI+5,  trim(TAMUoic2019path)//'/TAMUoic2019P21.bin', Tab_reclp, 'old')
         call open_dir(IUI+6,  trim(TAMUoic2019path)//'/TAMUoic2019P22.bin', Tab_reclp, 'old')
         call open_dir(IUI+7,  trim(TAMUoic2019path)//'/TAMUoic2019P23.bin', Tab_reclp, 'old')
         call open_dir(IUI+8,  trim(TAMUoic2019path)//'/TAMUoic2019P24.bin', Tab_reclp, 'old')
         call open_dir(IUI+9,  trim(TAMUoic2019path)//'/TAMUoic2019P31.bin', Tab_reclp, 'old')
         call open_dir(IUI+10, trim(TAMUoic2019path)//'/TAMUoic2019P32.bin', Tab_reclp, 'old')
         call open_dir(IUI+11, trim(TAMUoic2019path)//'/TAMUoic2019P33.bin', Tab_reclp, 'old')
         call open_dir(IUI+12, trim(TAMUoic2019path)//'/TAMUoic2019P34.bin', Tab_reclp, 'old')
         call open_dir(IUI+13, trim(TAMUoic2019path)//'/TAMUoic2019P41.bin', Tab_reclp, 'old')
         call open_dir(IUI+14, trim(TAMUoic2019path)//'/TAMUoic2019P42.bin', Tab_reclp, 'old')
         call open_dir(IUI+15, trim(TAMUoic2019path)//'/TAMUoic2019P43.bin', Tab_reclp, 'old')
         call open_dir(IUI+16, trim(TAMUoic2019path)//'/TAMUoic2019P44.bin', Tab_reclp, 'old')
      endif


      ! Warning
      mwarn = 0
      if (maxval(Oic_siz) > maxval(Tab_siz) .or. minval(Oic_siz) < minval(Tab_siz)) then
         print*,  ' ERROR: A specified maximum diameter range', minval(Oic_siz), '-', maxval(Oic_siz), &
         'should be within 50-10,000 (µm).'
         mwarn = 1 
      end if
      if (maxval(Oic_wav) > maxval(Tab_wav) .or. minval(Oic_wav) < minval(Tab_wav)) then
         print*,  ' ERROR: A specified wavelength range', minval(Oic_wav), '-', maxval(Oic_wav), &
         'should be within 0.355-1.064 (µm).'
         mwarn = 1 
      end if
      if (maxval(Oic_iang) > maxval(Tab_iang) .or. minval(Oic_iang) < minval(Tab_iang)) then
         print*,  ' ERROR: A specified incident polar angle range', minval(Oic_iang), '-', maxval(Oic_iang), &
         'should be within 0-90 (°).'
         mwarn = 1 
      end if
      if (mwarn == 1) stop ' TAMUoic2019 encoder halted. Please modify the input file.'

   end subroutine TAMUoic2019__init


   !+
   ! TAMUoic2019: Finalize the database
   !-
   subroutine TAMUoic2019__final()

      integer :: ipol
      deallocate( Oic_siz, Oic_wav, Oic_pang, Oic_aang, Oic_iang )
      do ipol = 0, Oic_npol
         close(IUI+ipol)
      end do

   end subroutine TAMUoic2019__final


   !+
   ! TAMUoic2019: Encode an ASCII database based on user-specified configuration 
   !-
   subroutine TAMUoic2019__encode()
     
      integer, parameter  :: NBIN = 8 ! 2**3
      integer  :: iwav,       isiz,       iiang,       ihab
      integer  :: jwav(NBIN), jsiz(NBIN), jiang(NBIN)   !, jhab(NBIN)
      integer  :: idxwav,     idxsiz,     idxiang !, idxhab
      real(R_) :: swav,       ssiz,       siang   !,   shab
      integer  :: ipol,   ibin
      integer  :: irec, ioss, iosp
      real(R_), allocatable :: wrk1s(:,:), wrk1p(:,:,:)     ! wrk for isca/p??.dat.    
      real(R_), allocatable :: wrk2s(:),   wrk2p(:,:)       ! wrk for isca/p??.dat.    
      real(R4_)             :: TAMUoic2019STabBuf(Tab_nmet)            ! buffer for isca.dat kernel table
      real(R4_)             :: TAMUoic2019PTabBuf(Tab_npang*Tab_naang) ! buffer for PMat.dat kernel table 

      ! Set properties
      allocate( wrk1s(Tab_nmet,NBIN), wrk1p(Tab_npang*Tab_naang,Tab_npol,NBIN) )
      allocate( wrk2s(Tab_nmet),      wrk2p(Tab_npang*Tab_naang,Tab_npol)      )
      call TAMUoic2019_database_open()
      if       (Oic_hab == 'OrientedPlate ') then
         ihab = 1 ! HOP
      else if  (Oic_hab == 'OrientedColumn') then 
         ihab = 2 ! HOC
      else
         stop 'ERROR: Oic_hab should be either OrientedPlate or OrientedColumn'
      end if
  
      ! Initialize binary coefficient
      do ibin = 1, NBIN
         jwav(ibin)  = int(mod((ibin-1)/4, 2)) 
         jiang(ibin) = int(mod((ibin-1)/2, 2)) 
         jsiz(ibin)  = int(mod( ibin-1,    2))
      end do 

      ! Interpolation
      do iwav = 1, Oic_nwav  ! wavelength
         call gridIdx_loc(Tab_wav(:), Oic_wav(iwav), idxwav, swav, 1)
         swav = max(0.0_R_, min(1.0_R_, swav))
         do iiang = 1, Oic_niang ! incident polar angle
            call gridIdx_loc(Tab_iang(:), Oic_iang(iiang), idxiang, siang, 1)
            siang = max(0.0_R_, min(1.0_R_, siang))
            do isiz = 1, Oic_nsiz ! maximum diameter
               call gridIdx_loc(Tab_siz(:), Oic_siz(isiz), idxsiz, ssiz, 1)
               ssiz = max(0.0_R_, min(1.0_R_, ssiz))

               ! Read data kernel
               do ibin = 1, NBIN
                  irec = idxsiz+jsiz(ibin) + Tab_nsiz*(idxiang-1+jiang(ibin) + &
                          Tab_niang*(idxwav-1+jwav(ibin) + Tab_nwav*(ihab-1)))
                  read(IUI, rec=irec, iostat=ioss) TAMUoic2019STabBuf(:) 
                  if (ioss /= 0) stop 'ERROR: reading isca.dat kernel table'
                  wrk1s(:,ibin) = TAMUoic2019STabBuf(:)
                  do ipol = 1, Oic_npol
                     read(IUI+ipol, rec=irec, iostat=iosp) TAMUoic2019PTabBuf(:) 
                     if (iosp /= 0) stop 'ERROR: reading PMat.dat kernel table'
                     wrk1p(:,ipol,ibin) = TAMUoic2019PTabBuf(:)
                  end do
               end do

               ! Get scattering properties for a single case
               wrk2s(:)   = (1.0_R_-ssiz)*(1.0_R_-siang)*(1.0_R_-swav) * wrk1s(:,1) & 
                          +         ssiz *(1.0_R_-siang)*(1.0_R_-swav) * wrk1s(:,2) &
                          + (1.0_R_-ssiz)*        siang *(1.0_R_-swav) * wrk1s(:,3) &
                          +         ssiz *        siang *(1.0_R_-swav) * wrk1s(:,4) &
                          + (1.0_R_-ssiz)*(1.0_R_-siang)*        swav  * wrk1s(:,5) &
                          +         ssiz *(1.0_R_-siang)*        swav  * wrk1s(:,6) &
                          + (1.0_R_-ssiz)*        siang *        swav  * wrk1s(:,7) &
                          +         ssiz *        siang *        swav  * wrk1s(:,8)
               wrk2p(:,:) = (1.0_R_-ssiz)*(1.0_R_-siang)*(1.0_R_-swav) * wrk1p(:,:,1) &
                          +         ssiz *(1.0_R_-siang)*(1.0_R_-swav) * wrk1p(:,:,2) &
                          + (1.0_R_-ssiz)*        siang *(1.0_R_-swav) * wrk1p(:,:,3) & 
                          +         ssiz *        siang *(1.0_R_-swav) * wrk1p(:,:,4) &
                          + (1.0_R_-ssiz)*(1.0_R_-siang)*        swav  * wrk1p(:,:,5) &
                          +         ssiz *(1.0_R_-siang)*        swav  * wrk1p(:,:,6) &
                          + (1.0_R_-ssiz)*        siang *        swav  * wrk1p(:,:,7) &
                          +         ssiz *        siang *        swav  * wrk1p(:,:,8) 

               ! Write 
               call TAMUoic2019_convert(wrk2s, wrk2p)
               call TAMUoic2019_write_database(isiz, iiang, iwav, wrk2s, wrk2p)
            end do
         end do
      end do
      call TAMUoic2019_database_close()
      deallocate(wrk1s, wrk1p, wrk2s, wrk2p)

   end subroutine TAMUoic2019__encode


   !+ 
   ! TAMUoic2019: Conversion from the kernel to single-scattering properties
   !-
   subroutine TAMUoic2019_convert(wrks, wrkp)
 
      real(R_), intent(inout)  :: wrks(:)     ! (nvar=6)           
      real(R_), intent(inout)  :: wrkp(:,:)   ! (npang*naang,npol) 
      real(R_) :: bufs(size(wrks, 1))
      integer  :: ipol
      !//NOTE 'wrks' includes (Vol, Ageo, Aproj, Cext, Csca, g*Csca) 
      !       'wrkp' includes (Csca*P??; P11, P12, P13, ..., P44)

      ! copy 
      bufs(:)   = wrks(:)

      ! Convert single-scattering properties
      wrks(4) = bufs(4)/bufs(2)  ! Qext = Cext/Ageo
      wrks(5) = bufs(5)/bufs(4)  ! SSA  = Csca/Cext
      wrks(6) = bufs(6)/bufs(5)  ! g    = g*Csca/Csca 

      ! Convert phase matrix and normalization (Csca*P?? in the database)
      if (Oic_npol > 1) then
         do ipol = 2, Oic_npol ! Csca*P?? --> P??/P11
            wrkp(:,ipol) = wrkp(:,ipol)/wrkp(:,1)
         end do 
      end if
      wrkp(:,1) = wrkp(:,1)/bufs(5) ! Csca*P11 --> P11
      !//Note (04/20/2021): TAMUoic2019 is now using a kernel technique.
      !
      ! Reference: Twomey, S. (1977). Introduction to the mathematics of
      ! inversion in remote sensing and indirect measurements. New York:
      ! Elsevier.

   end subroutine TAMUoic2019_convert


   !+ 
   ! TAMUoic2019: Output a database 
   !-
   subroutine TAMUoic2019_write_database(isiz, iiang, iwav, isca, pmat)
 
      integer, intent(in)  :: isiz, iiang, iwav
      real(R_), intent(in) :: isca(:)     ! (nvar=8) This will be (nvar=6)         
      real(R_), intent(in) :: pmat(:,:)   ! (npang*naang,npol) 
      integer :: iaang, ipol

      ! Write 
      write(IUO, '(9es16.7e2)') Oic_wav(iwav), Oic_siz(isiz), Oic_iang(iiang), isca(:)
      do ipol = 1, Oic_npol 
         do iaang = 1, Oic_naang
            write(IUO+ipol, '(1000es16.7e2)') pmat(1+(iaang-1)*Oic_npang:iaang*Oic_npang,ipol) 
         end do
      end do

   end subroutine TAMUoic2019_write_database


   !+
   ! TAMUoic2019: Open a database
   !-
   subroutine TAMUoic2019_database_open()

     character(len=8) :: fname(0:Oic_npol)  ! Filename of input file
     integer  :: ipol

     ! Set filenames
     fname(0) = "isca.dat"
     fname(1) = "P11.dat"
     if (Oic_npol > 1) then
        fname(2)  = "P12.dat"
        fname(3)  = "P13.dat"
        fname(4)  = "P14.dat"
        fname(5)  = "P21.dat"
        fname(6)  = "P22.dat"
        fname(7)  = "P23.dat"
        fname(8)  = "P24.dat"
        fname(9)  = "P31.dat"
        fname(10) = "P32.dat"
        fname(11) = "P33.dat"
        fname(12) = "P34.dat"
        fname(13) = "P41.dat"
        fname(14) = "P42.dat"
        fname(15) = "P43.dat"
        fname(16) = "P44.dat"
     end if

     ! Open the isca.dat and P??.dat
     call open_seq(IUO, trim(Oic_outpath) //'/'// trim(fname(0)), 'unknown')
     do ipol = 1, Oic_npol
        call open_seq(IUO+ipol, trim(Oic_outpath) //'/'// trim(fname(ipol)), 'unknown')
        write(IUO+ipol, '(1000es16.7e2)') Oic_pang(1:Oic_npang)
     end do

   end subroutine TAMUoic2019_database_open


   !+
   ! TAMUoic2019: Close a database
   !-
   subroutine TAMUoic2019_database_close()

     integer  :: ipol

     ! close files
     do ipol = 0, Oic_npol
        close(IUO+ipol)
     end do

   end subroutine TAMUoic2019_database_close


   !+ 
   ! TAMUoic2019: Set dimension in the kernel table
   !-
   subroutine TAMUoic2019__setDim()

      integer  :: idx

      ! Tab bin: Maximum diameter, wavelength, and particle habit
      do idx = 1, Tab_nsiz
         if (idx <= 21) then
            Tab_siz(idx) = 5.0000000E+01_R_ + 2.5000000E+00_R_ * real(idx-1)
         elseif(idx <= 93) then
            Tab_siz(idx) = 1.0000000E+02_R_ + 1.2500000E+01_R_ * real(idx-21)
         else
            Tab_siz(idx) = 1.0000000E+03_R_ + 1.2500000E+02_R_ * real(idx-93)
         endif
      end do
      Tab_wav(1:Tab_nwav) = (/ 3.5500000E-01_R_, 5.3200000E-01_R_, 1.0640000E+00_R_ /)
      Tab_hab(1:Tab_nhab) = (/ 'OrientedPlate ', 'OrientedColumn' /)

      ! Tab bin: Scattering polar angle
      do idx = 1, Tab_npang
         if    (idx <= 201) then
            Tab_pang(idx) =                    1.0000000E-02_R_ * real(idx-1)
         elseif(idx <= 261) then
            Tab_pang(idx) = 2.0000000E+00_R_ + 5.0000000E-02_R_ * real(idx-201)
         elseif(idx <= 311) then
            Tab_pang(idx) = 5.0000000E+00_R_ + 1.0000000E-01_R_ * real(idx-261)
         elseif(idx <= 321) then
            Tab_pang(idx) = 1.0000000E+01_R_ + 5.0000000E-01_R_ * real(idx-311)
         elseif(idx <= 482) then
            Tab_pang(idx) = 1.5000000E+01_R_ + 1.0000000E+00_R_ * real(idx-321)
         else
            Tab_pang(idx) = 1.7600000E+02_R_ + 2.5000000E-01_R_ * real(idx-482)
         endif
      end do

      ! Tab bin: Scattering azimuthal angle
      do idx = 1, Tab_naang
         Tab_aang(idx) = 2.0000000E+00_R_ * real(idx-1) 
      end do

      ! Tab bin: Incident polar angle
      do idx = 1, Tab_niang
         if    (idx <= 31) then
            Tab_iang(idx) =                    1.0000000E-02_R_ * real(idx-1)
         elseif(idx <= 58) then
            Tab_iang(idx) = 3.0000000E-01_R_ + 1.0000000E-01_R_ * real(idx-31)
         elseif(idx <= 72) then
            Tab_iang(idx) = 3.0000000E+00_R_ + 5.0000000E-01_R_ * real(idx-58)
         elseif(idx <= 82) then
            Tab_iang(idx) = 1.0000000E+01_R_ + 2.0000000E+00_R_ * real(idx-72)
         else
            Tab_iang(idx) = 3.0000000E+01_R_ + 5.0000000E+00_R_ * real(idx-82)
         endif
      end do

   end subroutine TAMUoic2019__setDim

end program tamuoic2019encode

