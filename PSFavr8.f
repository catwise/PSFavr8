c PSFavr8: average PSFs for cryo/postcryo and rotate
c
c            Supports Option 0 and Option 1
c
c            Focal-Plane PSF file names:
c
c            w1-psf-cryo-wpro.fits
c            w1-psf-postcryo-wpro.fits
c            w1-psfunc-cryo-wpro.fits
c            w1-psfunc-postcryo-wpro.fits
c            w2-psf-cryo-wpro.fits
c            w2-psf-postcryo-wpro.fits
c            w2-psfunc-cryo-wpro.fits
c            w2-psfunc-postcryo-wpro.fits
c
c            WPhotpmc PSF file names (in Asce and Desc tile subdirs):
c
c            unWISE-w1-psf-wpro-01x01-01x01.fits
c            unWISE-w1-psfunc-wpro-01x01-01x01.fits
c            unWISE-w2-psf-wpro-01x01-01x01.fits
c            unWISE-w2-psfunc-wpro-01x01-01x01.fits
c            unWISE-w3-psf-wpro-01x01-01x01.fits        (if generated)
c            unWISE-w3-psfunc-wpro-01x01-01x01.fits     (if generated)
c            unWISE-w4-psf-wpro-01x01-01x01.fits        (if generated)
c            unWISE-w4-psfunc-wpro-01x01-01x01.fits     (if generated)
c
c            ICORE PSF file names (in tile directory):
c
c            unwise-w1-psf-awaic-01x01-01x01.fits
c            unwise-w2-psf-awaic-01x01-01x01.fits
c
c vsn 1.0  B70725 initial version, cloned from PSFavr8.f
c vsn 1.1  B70804 fixed cryo/postcryo averging for case of no cryo or
c                 postcryo for given asce or desc.
c vsn 1.2  B70830 changes error-exit stop statements to call exit(64) 
c vsn 1.3  B70914 cloned from PSFavr8sd; added option 0 support
c vsn 1.4  B70916 added optional scale factors for psfuncs
c vsn 1.5  B71004 made -sc1 and -sp1 defaults 1.65
c vsn 1.6  B71012 changed psfunc scaling from input to output; no more
c                 -sc1 -sc2 -sp1 -sp2, now s01 s02 s11 s22 for option-0
c                 W1 and W2 resp., and option-1 W1 and W2 resp.
c vsn 1.6  B71015 added comments about -h and -mh to tutorial
c vsn 1.6  B71020 added handling for PA zero crossings by allowing up
c                 3600 histogram cells and jumping past those with
c                 zero counts
c vsn 1.6  B71025 added missing CTYPE# lines to awaic headers
c vsn 1.6  B71030 changed CTYPE# from "sin" to "tan" (WISE to unWISE)
c vsn 1.6  B71122 changed MJD0 from start of hibernation to end of
c                 4-band cryo
c vsn 1.6  B71128 made MJD0 a command-line parameter
c vsn 1.6  B71207 eliminates "-w" because it is not needed; left in
c                 Slash variable for convenience, '/' default
c vsn 1.61 B80111 added psfunc scale factor to FITS header
c vsn 1.62 B80128 set default psfunc scale factors back to 1.0
c          B80226 changed "end of cryo" to "end of cryo PSF"
c vsn 1.63 B80306 delete FITS files with output names that already exist
c
c=======================================================================
c
      Integer*4     IArgC, LNBlnk, FileID, nOSc, nEpoch,
     +              NPlanes, NRows, NCols, I, J, N, IStat,
     +              ImRPixR, NArgs, NArg, NPix, status, n1ang, n2ang,
     +              naxes(2), nOvrSamp, nAngsw1ac, nAngsw1dc, nAngsw2ac,
     +              nAngsw2dc, nAngsw1ap, nAngsw1dp, nAngsw2ap,
     +              nAngsw2dp, ndir, nHmax, nW1acryo, nW1apostcryo,
     +              nW2acryo, nW2apostcryo, nW1dcryo, nW1dpostcryo,
     +              nW2dcryo, nW2dpostcryo
      Character*200 InFNam, OutFNam, CanPath,
     +              OutPath, w1paNam, w2paNam
      Character*80  hdrline(100), record
      Character*65  EpochStr
      Character*11  Vsn, TmpStr, Flag
      Character*9   TmpStr9
      Character*8   cdate, ctime, TileNam
      Character*2   TmpStr2
      Character*1   Slash
      Logical*4     dbg, GotIN, GotOut, w3, w4, da, ge,
     +              GotAng1, GotAng2
c                                                            ! Focal-Plane PSFs
      Real*4        w1psfc(641,641), w1psfuncc(641,641),     ! W1 cryo
     +              w1psfp(641,641), w1psfuncp(641,641),     ! W1 postcryo
     +              w2psfc(641,641), w2psfuncc(641,641),     ! W2 cryo
     +              w2psfp(641,641), w2psfuncp(641,641),     ! W2 postcryo
     +              wgt(641,641)
c                                                            ! Output wpro PSFs
      Real*4        w1psfca(641,641), w1psfuncca(641,641),   ! W1 cryo asce
     +              w1psfcd(641,641), w1psfunccd(641,641),   ! W1 cryo desc
     +              w2psfca(641,641), w2psfuncca(641,641),   ! W2 cryo asce
     +              w2psfcd(641,641), w2psfunccd(641,641),   ! W2 cryo desc
     +              w1psfpa(641,641), w1psfuncpa(641,641),   ! W1 postcryo asce
     +              w1psfpd(641,641), w1psfuncpd(641,641),   ! W1 postcryo desc
     +              w2psfpa(641,641), w2psfuncpa(641,641),   ! W2 postcryo asce
     +              w2psfpd(641,641), w2psfuncpd(641,641)    ! W2 postcryo desc
c
c NOTE: W1 cryo & postcryo will be epoch-averaged and only one output
c       unless W3 = T, in which case W1 cryo will be output as W3, and
c       W1 postcryo will be output as W1 with no epoch averaging.
c     
c       W2 cryo & postcryo will be epoch-averaged and only one output
c       unless W4 = T, in which case W2 cryo will be output as W4, and
c       W2 postcryo will be output as W2 with no epoch averaging.
c
c       When cryo and post-cryo are averaged, they are averaged into the
c       cryo array, and that is used for output.
c
c                                                            ! Output awaic PSFs
      Real*4        w1psfawaic(27,27),    w2psfawaic(27,27)  ! W1 & W2
c                                                            ! Option 0 PSFs
      Real*4       w1psfopt0(641,641),    w2psfopt0(641,641),! W1 & W2
     +          w1psfuncopt0(641,641), w2psfuncopt0(641,641)
c
      Real*4, allocatable :: w1angsac(:), w2angsac(:),       ! ascending  cryo
     +                       w1angsdc(:), w2angsdc(:),       ! descending cryo
     +                       w1angsap(:), w2angsap(:),       ! ascending  postcryo
     +                       w1angsdp(:), w2angsdp(:)        ! descending postcryo
      Real*8        Ang1acMin, Ang1acMax, Ang2acMin, PixSum,
     +              Ang2acMax, Ang1dcMin, Ang1dcMax, Ang2dcMin,
     +              Ang2dcMax, Ang1apMin, Ang1apMax, Ang2apMin,
     +              Ang2apMax, Ang1dpMin, Ang1dpMax, Ang2dpMin,
     +              Ang2dpMax, wgt1, wgt2, MJD, MJD0, w1ang, w2ang,
     +              f1ac, f1dc, f1ap, f1dp, f2ac, f2dc, f2ap, f2dp,
     +              MJD1, MJD2
      real*4        dH1ac, dH1dc, Dh2ac, dH2dc,
     +              dH1ap, dH1dp, Dh2ap, dH2dp,
     +              nullval, dH, xlo, xhi, ylo, yhi,
     +              s01, s02, s11, s12
      logical*4 anynull
c
      common /vdt/ cdate,ctime,vsn
c
      Data Vsn/'1.63 B80306'/, nOvrSamp/11/, nHmax/3600/, dH/0.1/,
     +     GotIn,GotOut,GotAng1,GotAng2/4*.false./, da/.false./,
     +     dbg/.false./, Slash/'/'/, TileNam/'NotGiven'/,
     +     Ang1acMin,Ang2acMin,Ang1dcMin,Ang2dcMin/4*99999.9/,
     +     Ang1acMax,Ang2acMax,Ang1dcMax,Ang2dcMax/4*-99999.9/,
     +     Ang1apMin,Ang2apMin,Ang1dpMin,Ang2dpMin/4*99999.9/,
     +     Ang1apMax,Ang2apMax,Ang1dpMax,Ang2dpMax/4*-99999.9/,
     +     w3/.false./, w4/.false./, ge/.false./, MJD2/-9.9/,
     +     nW1acryo,nW1apostcryo,nW2acryo,nW2apostcryo/4*0/,
     +     nW1dcryo,nW1dpostcryo,nW2dcryo,nW2dpostcryo/4*0/,
     +     nEpoch/-1/, s01/1.0/, s02/1.0/, s11/1.0/, s12/1.0/
c     
c     Data MJD0/55469.278/ ! JD = 2455468.5, September 29, 2010
c                          !                 end of 3-band cryo
c     Data MJD0/55994.0/   ! JD = 2455468.5, February 1, 2011
c                          !                 start of hibernation
c     Data MJD0/55414.932/ ! JD = 2455414.441, August 6, 2010
c                          !                   end of 4-band cryo
      Data MJD0/55480.0/   ! JD = 2455480.0, October 11, 2010
c                                            Peter's Preference
c                                            end of cryo PSF
c
c=======================================================================
c
      NArgs = IArgC()
      If (NArgs .lt. 1) then
        print *,'PSFavr8 vsn ', Vsn
        print *
        print *,'Usage: PSFavr8 <flags specifications>'
        print *
        print *,'Where the REQUIRED flags and specifications are:'
        print *,'    -i   pathname for the focal-plane PSFs'
        print *,'    -o   directory name for output PSF files'
        print *,'         (must have "Asce" and "Desc" subdirectories)'
        print *,'    -a1  W1 filename for rotation angles'
        print *,'    -a2  W2 filename for rotation angles'
        print *
        print *,'The OPTIONAL flags and specifications are:'
        print *,'    -t   tile name (e.g., 1124p045)'
        print *,'    -s01 scale factor for W1 option-0 psfunc (1.0)'
        print *,'    -s02 scale factor for W2 option-0 psfunc (1.0)'
        print *,'    -s11 scale factor for W1 option-1 psfunc (1.0)'
        print *,'    -s12 scale factor for W2 option-1 psfunc (1.0)'
        print *,'    -w3  output W1 cryo PSF as W3'
        print *,'         (no cryo/postcryo averaging)'
        print *,'    -w4  output W2 cryo PSF as W4'
        print *,'         (no cryo/postcryo averaging)'
        print *,'    -os  oversampling factor for interpolation (11)'
        print *,'    -e   generate a tile epoch text file'
        print *,'    -h   histogram resolution for angles (0.1 deg)'
        print *,'    -mh  maximum number of histogram cells (3600)'
        print *,'    -d   turn on debug prints'
        print *,'    -da  dump angle histograms to stdout'
        print *,'    -m   MJD separating cryo and post-cryo (55480.0)'
c       print *,'    -w   testing on a Windows machine' ! don't need this after all
        print *
        print *,'If the -h and -mh values cannot span the angle range,'
        print *,'the histogram resolution (-h) will be increased as'
        print *,'needed to do so.'
        stop
      end if
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      NArg = 0
5     NArg = NArg + 1
      call GetArg(NArg,Flag)
      call UpCase(Flag)
c                                      ! input path for focal-plane PSFs
      If (Flag .eq. '-I') then
        call NextNarg(NArg,Nargs)
        call GetArg(NArg,CanPath)
        if (CanPath(LNBlnk(CanPath):LNBlnk(CanPath)) .ne. Slash)
     +  CanPath = CanPath(1:LNBlnk(CanPath))//Slash  
        GotIn = .true.
c                                      ! Turn debug prints on
      else if (Flag .eq. '-D') then
        dbg = .true.
        print *,'Debug prints enabled'
c                                      ! Dump angle histograms
      else if (Flag .eq. '-DA') then
        da = .true.
        print *,'Angle histograms will be dumped to sdtout'
c                                  
      else if (Flag .eq. '-E') then    ! Generate tile epoch text file
        ge = .true.
        print *,'Tile epoch text file will be generated'
c                                      ! output filename stem
      else if (Flag .eq. '-O') then
        call NextNarg(NArg,Nargs)
        call GetArg(NArg,OutPath)
        if (OutPath(LNBlnk(OutPath):LNBlnk(OutPath)) .ne. Slash)
     +  OutPath = OutPath(1:LNBlnk(OutPath))//Slash 
        GotOut = .true.
c                                      ! W1 rotation angles
      else if (Flag .eq. '-T') then
        call NextNarg(NArg,Nargs)
        call GetArg(NArg,TileNam)
        if (dbg) print *,'TileNam = ', TileNam
c                                   
      else if (Flag .eq. '-A1') then   ! W1 rotation angles
        call NextNarg(NArg,Nargs)
        call GetArg(NArg,w1paNam)
        if (dbg) print *,'w1paNam = ', w1paNam(1:lnblnk(w1paNam))
        GotAng1 = .true.
        if (Access(w1paNam(1:lnblnk(w1paNam)),' ') .ne. 0) then
          print *,'File not found: ',w1paNam(1:lnblnk(w1paNam))
          call exit(64)
         end if
c                                      ! W2 rotation angles
      else if (Flag .eq. '-A2') then
        call NextNarg(NArg,Nargs)
        call GetArg(NArg,w2paNam)
        if (dbg) print *,'w2paNam = ', w2paNam(1:lnblnk(w2paNam))
        GotAng2 = .true.
        if (Access(w2paNam(1:lnblnk(w2paNam)),' ') .ne. 0) then
          print *,'File not found: ',w2paNam(1:lnblnk(w2paNam))
          call exit(64)
         end if
c                                      ! Oversampling factor
      else if (Flag .eq. '-OS') then
        call NextNarg(NArg,Nargs)
        call GetArg(NArg,TmpStr)
        read (TmpStr, *, err=3000) nOvrSamp
        if (dbg) print *,'Oversampling factor = ',nOvrSamp
        if (mod(nOvrSamp,2) .eq. 0) then
          nOvrSamp = nOvrSamp + 1
          print *,'Oversampling factor must be odd; changed to ',nOvrSamp
        end if        
c                                      ! W1 option-0 psfunc scale factor
      else if (Flag .eq. '-S01') then
        call NextNarg(NArg,Nargs)
        call GetArg(NArg,TmpStr)
        read (TmpStr, *, err=3000) s01
        if (dbg) print *,'W1 option-0 psfunc scale factor = ',s01
c                                    
      else if (Flag .eq. '-S02') then  ! W2 option-0 psfunc scale factor
        call NextNarg(NArg,Nargs)
        call GetArg(NArg,TmpStr)
        read (TmpStr, *, err=3000) s02
        if (dbg) print *,'W2 option-0 psfunc scale factor = ',s02
c                                    
      else if (Flag .eq. '-S11') then  ! W1 option-1 psfunc scale factor
        call NextNarg(NArg,Nargs)
        call GetArg(NArg,TmpStr)
        read (TmpStr, *, err=3000) s11
        if (dbg) print *,'W1 option-1 psfunc scale factor = ',s11
c                                    
      else if (Flag .eq. '-S12') then  ! W2 option-1 psfunc scale factor
        call NextNarg(NArg,Nargs)
        call GetArg(NArg,TmpStr)
        read (TmpStr, *, err=3000) s12
        if (dbg) print *,'W2 option-1 psfunc scale factor = ',s12
c                                  
      else if (Flag .eq. '-H') then    ! Histogram resolution
        call NextNarg(NArg,Nargs)
        call GetArg(NArg,TmpStr)
        read (TmpStr, *, err=3000) dH
        if (dbg) print *,'Histogram resolution = ',dH
        if (dH .le. 0.0) then
          print *,'Histogram resolution must be > 0'
          call exit(64)
        end if        
c                                  
      else if (Flag .eq. '-M') then    ! MJD0
        call NextNarg(NArg,Nargs)
        call GetArg(NArg,TmpStr)
        read (TmpStr, *, err=3000) MJD0
        if (dbg) print *,'MJD0 = ',MJD0
c                                   
      else if (Flag .eq. '-MH') then   ! Max # Histogram cells
        call NextNarg(NArg,Nargs)
        call GetArg(NArg,TmpStr)
        read (TmpStr, *, err=3000) nHmax
        if (dbg) print *,'Max # Histogram cells = ',nHmax
        if (nHmax .le. 0) then
          print *,'Max # Histogram cells must be > 0'
          call exit(64)
        end if        
c                                      ! W1 cryo as W3
      else if (Flag .eq. '-W3') then
        w3 = .true.
        if (dbg) print *,'w3 = T'
c                                      ! W1 cryo as W3
      else if (Flag .eq. '-W4') then
        w4 = .true.
        if (dbg) print *,'w4 = T'
c                                      ! Windows directory backslashes
c     else if (Flag .eq. '-W') then    ! Windows DOS wil accept forward
c       Slash = '\'                    ! slash in file names
c       if (dbg) print *,'Window test; "slash" character = "\"'
      end if
c 
      If (NArg .lt. NArgs) Go to 5
      call signon('PSFavr8')
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c                                      ! Verify input
      If (.not.GotIn) then
        print *,'ERROR: Input path to focal-plane PSFs not specified'
        call exit(64)
      end if
c
      If (.not.GotOut) then
        print *,'ERROR: Output filename stem not specified'
        call exit(64)
      end if
c
      If (.not.GotAng1) then
        print *,'ERROR: W1 Rotation angle filename not specified'
        call exit(64)
      end if
c
      If (.not.GotAng2) then
        print *,'ERROR: W2 Rotation angle filename not specified'
        call exit(64)
      end if
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c                                      ! Set up tile epoch text file
      if (ge) open (12, file = TileNam//'-epochs.txt')
c      
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c                                      ! Set up rotation angles for W1
      open (10, file = w1paNam)
      n1ang = 0
10    read (10,*, end=15,err=3001) ndir, w1ang, MJD
      if (ge) then
        if (MJD-MJD2 .gt. 9.9) then    ! NOTE: this will probably not
          if (nEpoch .ge. 0) then      !       work near a pole
            write(TmpStr9,'(F9.2)') MJD2
            EpochStr = EpochStr(1:lnblnk(EpochStr))//Tmpstr9
            write (12,*) EpochStr
          end if
          nEpoch = nEpoch + 1
          write(TmpStr2,'(I2)') nEpoch
          write(TmpStr9,'(F9.2)') MJD
          if (ndir .eq. 1) then
            EpochStr = 'Epoch '//TmpStr2//':  Ascending'
          else
            EpochStr = 'Epoch '//TmpStr2//': Descending'
          end if
          if (MJD .lt. MJD0) then
            EpochStr = EpochStr(1:lnblnk(EpochStr))
     +                 //' Cryo,      MJD Range ='//TmpStr9//' to' 
          else
            EpochStr = EpochStr(1:lnblnk(EpochStr))
     +                 //' Post-Cryo, MJD Range ='//TmpStr9//' to'
          end if
        end if
        MJD2 = MJD
      end if
      n1ang = n1ang + 1
      if (ndir .eq. 1) then                           ! 1 --> ascending scan
        if (MJD .lt. MJD0) then                       ! cryo
          nW1acryo = nW1acryo + 1
          if (w1ang .gt. Ang1acMax) Ang1acMax = w1ang
          if (w1ang .lt. Ang1acMin) Ang1acMin = w1ang
        else                                          ! postcryo
          nW1apostcryo = nW1apostcryo + 1
          if (w1ang .gt. Ang1apMax) Ang1apMax = w1ang
          if (w1ang .lt. Ang1apMin) Ang1apMin = w1ang
        end if
      end if
      if (ndir .eq. 0) then                           ! 0 --> descending scan
        if (MJD .lt. MJD0) then                       ! cryo
          nW1dcryo = nW1dcryo + 1
          if (w1ang .gt. Ang1dcMax) Ang1dcMax = w1ang
          if (w1ang .lt. Ang1dcMin) Ang1dcMin = w1ang
        else                                          ! postcryo
          nW1dpostcryo = nW1dpostcryo + 1
          if (w1ang .gt. Ang1dpMax) Ang1dpMax = w1ang
          if (w1ang .lt. Ang1dpMin) Ang1dpMin = w1ang
        end if
      end if
      go to 10
c      
15    rewind(10)
c                                      ! Gen Epoch File if requested
      if (ge) then
        if (nEpoch .ge. 0) then
          write(TmpStr9,'(F9.2)') MJD2
          EpochStr = EpochStr(1:lnblnk(EpochStr))//Tmpstr9
          write (12,*) EpochStr
        end if
      end if 
c                                      ! Allocate W1 angle arrays
      dH1ac = dH
      nAngsw1ac = NInt((Ang1acMax - Ang1acMin)/dH1ac)
      if (nAngsw1ac .lt. 1) nAngsw1ac = 1
      if (nAngsw1ac .gt. nHmax) then
        nAngsw1ac = nHmax
        print *,'WARNING: adjusting angle histogram for W1 ac'
      end if     
      dH1ac = (Ang1acMax - Ang1acMin)/float(nAngsw1ac)
      allocate(w1angsac(nAngsw1ac))
      if (.not.allocated(w1angsac)) then
        print *,'ERROR: allocation of w1angsac failed'
        print *,'       no. elements =',nAngsw1ac
        call exit(64)
      end if
c      
      dH1ap = dH
      nAngsw1ap = NInt((Ang1apMax - Ang1apMin)/dH1ap)
      if (nAngsw1ap .lt. 1) nAngsw1ap = 1
      if (nAngsw1ap .gt. nHmax) then
        nAngsw1ap = nHmax
        print *,'WARNING: adjusting angle histogram for W1 ap'
      end if     
      dH1ap = (Ang1apMax - Ang1apMin)/float(nAngsw1ap)
      allocate(w1angsap(nAngsw1ap))
      if (.not.allocated(w1angsap)) then
        print *,'ERROR: allocation of w1angsap failed'
        print *,'       no. elements =',nAngsw1ap
        call exit(64)
      end if
c      
      dH1dc = dH
      nAngsw1dc = NInt((Ang1dcMax - Ang1dcMin)/dH1dc)
      if (nAngsw1dc .lt. 1) nAngsw1dc = 1
      if (nAngsw1dc .gt. nHmax) then
        nAngsw1dc = nHmax
        print *,'WARNING: adjusting angle histogram for W1 dc'
      end if     
      dH1dc = (Ang1dcMax - Ang1dcMin)/float(nAngsw1dc)
      allocate(w1angsdc(nAngsw1dc))
      if (.not.allocated(w1angsdc)) then
        print *,'ERROR: allocation of w1angsdc failed'
        print *,'       no. elements =',nAngsw1dc
        call exit(64)
      end if
c      
      dH1dp = dH
      nAngsw1dp = NInt((Ang1dpMax - Ang1dpMin)/dH1dp)
      if (nAngsw1dp .lt. 1) nAngsw1dp = 1
      if (nAngsw1dp .gt. nHmax) then
        nAngsw1dp = nHmax
        print *,'WARNING: adjusting angle histogram for W1 dp'
      end if     
      dH1dp = (Ang1dpMax - Ang1dpMin)/float(nAngsw1dp)
      allocate(w1angsdp(nAngsw1dp))
      if (.not.allocated(w1angsdp)) then
        print *,'ERROR: allocation of w1angsdp failed'
        print *,'       no. elements =',nAngsw1dp
        call exit(64)
      end if
c
      w1angsac = 0.0
      w1angsap = 0.0
      w1angsdc = 0.0
      w1angsdp = 0.0
c     
      do 20 I = 1, n1ang
        read(10, *, err=3001) ndir, w1ang, MJD
        if (ndir .eq. 1) then
          if (MJD .lt. MJD0) then
            n = NInt((w1ang - Ang1acMin)/dH1ac + 0.5)
            if (n .lt. 1) n = 1
            if (n .gt. nAngsw1ac) n = nAngsw1ac
            w1angsac(n) = w1angsac(n) + 1.0
          else
            n = NInt((w1ang - Ang1apMin)/dH1ap + 0.5)
            if (n .lt. 1) n = 1
            if (n .gt. nAngsw1ap) n = nAngsw1ap
            w1angsap(n) = w1angsap(n) + 1.0
          end if
        end if
        if (ndir .eq. 0) then
          if (MJD .lt. MJD0) then
            n = NInt((w1ang - Ang1dcMin)/dH1dc + 0.5)
            if (n .lt. 1) n = 1
            if (n .gt. nAngsw1dc) n = nAngsw1dc
            w1angsdc(n) = w1angsdc(n) + 1.0
          else
            n = NInt((w1ang - Ang1dpMin)/dH1dp + 0.5)
            if (n .lt. 1) n = 1
            if (n .gt. nAngsw1dp) n = nAngsw1dp
            w1angsdp(n) = w1angsdp(n) + 1.0
          end if
        end if
20    continue
      close(10)
      call ChkPAX(w1angsac,nAngsw1ac,dH1ac)      
      call ChkPAX(w1angsdc,nAngsw1dc,dH1dc)      
      call ChkPAX(w1angsap,nAngsw1ap,dH1ap)      
      call ChkPAX(w1angsdp,nAngsw1dp,dH1dp)      
c
      if (da) then
        print *
        print *,'W1 Cryo Ascending angle histogram (bin# Angle Count):'
        do 22 n = 1, nAngsw1ac
          w1ang = (float(n) - 0.5)* dH1ac + Ang1acMin
          if ((w1angsac(n) .gt. 0.0) .or. (n .eq. 1))
     +    write (6, '(I4,2F10.2)') n, w1ang, w1angsac(n)
22      continue
        print *
        print *,'W1 Cryo Descending angle histogram (bin# Angle Count):'
        do 24 n = 1, nAngsw1dc
          w1ang = (float(n) - 0.5)* dH1dc + Ang1dcMin
          if ((w1angsdc(n) .gt. 0.0) .or. (n .eq. 1))
     +    write (6, '(I4,2F10.2)') n, w1ang, w1angsdc(n)
24      continue
        print *
        print *,'W1 Post-Cryo Ascending angle histogram (bin# Angle Count):'
        do 26 n = 1, nAngsw1ap
          w1ang = (float(n) - 0.5)* dH1ap + Ang1apMin
          if ((w1angsap(n) .gt. 0.0) .or. (n .eq. 1))
     +    write (6, '(I4,2F10.2)') n, w1ang, w1angsap(n)
26      continue
        print *
        print *,'W1 Post-Cryo Descending angle histogram (bin# Angle Count):'
        do 28 n = 1, nAngsw1dp
          w1ang = (float(n) - 0.5)* dH1dp + Ang1dpMin
          if ((w1angsdp(n) .gt. 0.0) .or. (n .eq. 1))
     +    write (6, '(I4,2F10.2)') n, w1ang, w1angsdp(n)
28      continue
      end if
c
      if (n1ang .gt. 0) then
        f1ac = dfloat(nW1acryo)/dfloat(n1ang)
        f1dc = dfloat(nW1dcryo)/dfloat(n1ang)
        f1ap = dfloat(nW1apostcryo)/dfloat(n1ang)
        f1dp = dfloat(nW1dpostcryo)/dfloat(n1ang)
      else
        f1ac = 0.0d0
        f1dc = 0.0d0
        f1ap = 0.0d0
        f1dp = 0.0d0
      end if
      if (dbg) then
        print *,'n1ang: ',n1ang
        print *,'nW1acryo, f1ac:     ', nW1acryo, f1ac
        print *,'nW1dcryo, f1dc:     ', nW1dcryo, f1dc
        print *,'nW1apostcryo, f1ap: ', nW1apostcryo, f1ap
        print *,'nW1dpostcryo, f1dp: ', nW1dpostcryo, f1dp
      end if
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c                                      ! Set up rotation angles for W2
      open (10, file = w2paNam)
      n2ang = 0
30    read (10,*, end=35,err=3002) ndir, w2ang, MJD
      n2ang = n2ang + 1
      if (ndir .eq. 1) then                           ! 1 --> ascending scan
        if (MJD .lt. MJD0) then                       ! cryo
          nW2acryo = nW2acryo + 1
          if (w2ang .gt. Ang2acMax) Ang2acMax = w2ang
          if (w2ang .lt. Ang2acMin) Ang2acMin = w2ang
        else                                          ! postcryo
          nW2apostcryo = nW2apostcryo + 1
          if (w2ang .gt. Ang2apMax) Ang2apMax = w2ang
          if (w2ang .lt. Ang2apMin) Ang2apMin = w2ang
        end if
      end if
      if (ndir .eq. 0) then                           ! 0 --> descending scan
        if (MJD .lt. MJD0) then                       ! cryo
          nW2dcryo = nW2dcryo + 1
          if (w2ang .gt. Ang2dcMax) Ang2dcMax = w2ang
          if (w2ang .lt. Ang2dcMin) Ang2dcMin = w2ang
        else                                          ! postcryo
          nW2dpostcryo = nW2dpostcryo + 1
          if (w2ang .gt. Ang2dpMax) Ang2dpMax = w2ang
          if (w2ang .lt. Ang2dpMin) Ang2dpMin = w2ang
        end if
      end if
      go to 30
c      
35    rewind(10)
      dH2ac = dH
      nAngsw2ac = NInt((Ang2acMax - Ang2acMin)/dH2ac)
      if (nAngsw2ac .lt. 1) nAngsw2ac = 1
      if (nAngsw2ac .gt. nHmax) then
        nAngsw2ac = nHmax
        print *,'WARNING: adjusting angle histogram for W2 ac'
      end if     
      dH2ac = (Ang2acMax - Ang2acMin)/float(nAngsw2ac)
      allocate(w2angsac(nAngsw2ac))
      if (.not.allocated(w2angsac)) then
        print *,'ERROR: allocation of w2angsac failed'
        print *,'       no. elements =',nAngsw2ac
        call exit(64)
      end if
c      
      dH2ap = dH
      nAngsw2ap = NInt((Ang2apMax - Ang2apMin)/dH2ap)
      if (nAngsw2ap .lt. 1) nAngsw2ap = 1
      if (nAngsw2ap .gt. nHmax) then
        nAngsw2ap = nHmax
        print *,'WARNING: adjusting angle histogram for W2 ap'
      end if     
      dH2ap = (Ang2apMax - Ang2apMin)/float(nAngsw2ap)
      allocate(w2angsap(nAngsw2ap))
      if (.not.allocated(w2angsap)) then
        print *,'ERROR: allocation of w2angsap failed'
        print *,'       no. elements =',nAngsw2ap
        call exit(64)
      end if
c      
      dH2dc = dH
      nAngsw2dc = NInt((Ang2dcMax - Ang2dcMin)/dH2dc)
      if (nAngsw2dc .lt. 1) nAngsw2dc = 1
      if (nAngsw2dc .gt. nHmax) then
        nAngsw2dc = nHmax
        print *,'WARNING: adjusting angle histogram for W2 dc'
      end if     
      dH2dc = (Ang2dcMax - Ang2dcMin)/float(nAngsw2dc)
      allocate(w2angsdc(nAngsw2dc))
      if (.not.allocated(w2angsdc)) then
        print *,'ERROR: allocation of w2angsdc failed'
        print *,'       no. elements =',nAngsw2dc
        call exit(64)
      end if
c      
      dH2dp = dH
      nAngsw2dp = NInt((Ang2dpMax - Ang2dpMin)/dH2dp)
      if (nAngsw2dp .lt. 1) nAngsw2dp = 1
      if (nAngsw2dp .gt. nHmax) then
        nAngsw2dp = nHmax
        print *,'WARNING: adjusting angle histogram for W2 dp'
      end if     
      dH2dp = (Ang2dpMax - Ang2dpMin)/float(nAngsw2dp)
      allocate(w2angsdp(nAngsw2dp))
      if (.not.allocated(w2angsdp)) then
        print *,'ERROR: allocation of w2angsdp failed'
        print *,'       no. elements =',nAngsw2dp
        call exit(64)
      end if
c
      w2angsac = 0.0
      w2angsap = 0.0
      w2angsdc = 0.0
      w2angsdp = 0.0
c     
      do 40 I = 1, n2ang
        read(10, *, err=3002) ndir, w2ang, MJD
        if (ndir .eq. 1) then
          if (MJD .lt. MJD0) then
            n = NInt((w2ang - Ang2acMin)/dH2ac + 0.5)
            if (n .lt. 1) n = 1
            if (n .gt. nAngsw2ac) n = nAngsw2ac
            w2angsac(n) = w2angsac(n) + 1.0
          else
            n = NInt((w2ang - Ang2apMin)/dH2ap + 0.5)
            if (n .lt. 1) n = 1
            if (n .gt. nAngsw2ap) n = nAngsw2ap
            w2angsap(n) = w2angsap(n) + 1.0
          end if
        end if
        if (ndir .eq. 0) then
          if (MJD .lt. MJD0) then
            n = NInt((w2ang - Ang2dcMin)/dH2dc + 0.5)
            if (n .lt. 1) n = 1
            if (n .gt. nAngsw2dc) n = nAngsw2dc
            w2angsdc(n) = w2angsdc(n) + 1.0
          else
            n = NInt((w2ang - Ang2dpMin)/dH2dp + 0.5)
            if (n .lt. 1) n = 1
            if (n .gt. nAngsw2dp) n = nAngsw2dp
            w2angsdp(n) = w2angsdp(n) + 1.0
          end if
        end if
40    continue
      close(10)
      call ChkPAX(w2angsac,nAngsw2ac,dH2ac)      
      call ChkPAX(w2angsdc,nAngsw2dc,dH2dc)      
      call ChkPAX(w2angsap,nAngsw2ap,dH2ap)      
      call ChkPAX(w2angsdp,nAngsw2dp,dH2dp)      
c      
      if (da) then
        print *
        print *,'W2 Cryo Ascending angle histogram (bin# Angle Count):'
        do 42 n = 1, nAngsw2ac
          w2ang = (float(n) - 0.5)* dH2ac + Ang2acMin
          if ((w2angsac(n) .gt. 0.0) .or. (n .eq. 1))
     +    write (6, '(I4,2F10.2)') n, w2ang, w2angsac(n)
42      continue
        print *
        print *,'W2 Cryo Descending angle histogram (bin# Angle Count):'
        do 44 n = 1, nAngsw2dc
          w2ang = (float(n) - 0.5)* dH2dc + Ang2dcMin
          if ((w2angsdc(n) .gt. 0.0) .or. (n .eq. 1))
     +    write (6, '(I4,2F10.2)') n, w2ang, w2angsdc(n)
44      continue
        print *
        print *,'W2 Post-Cryo Ascending angle histogram (bin# Angle Count):'
        do 46 n = 1, nAngsw2ap
          w2ang = (float(n) - 0.5)* dH2ap + Ang2apMin
          if ((w2angsap(n) .gt. 0.0) .or. (n .eq. 1))
     +    write (6, '(I4,2F10.2)') n, w2ang, w2angsap(n)
46      continue
        print *
        print *,'W2 Post-Cryo Descending angle histogram (bin# Angle Count):'
        do 48 n = 1, nAngsw2dp
          w2ang = (float(n) - 0.5)* dH2dp + Ang2dpMin
          if ((w2angsdp(n) .gt. 0.0) .or. (n .eq. 1))
     +    write (6, '(I4,2F10.2)') n, w2ang, w2angsdp(n)
48      continue
      end if
c
      if (n2ang .gt. 0) then
        f2ac = dfloat(nW2acryo)/dfloat(n2ang)
        f2dc = dfloat(nW2dcryo)/dfloat(n2ang)
        f2ap = dfloat(nW2apostcryo)/dfloat(n2ang)
        f2dp = dfloat(nW2dpostcryo)/dfloat(n2ang)
      else
        f2ac = 0.0d0
        f2dc = 0.0d0
        f2ap = 0.0d0
        f2dp = 0.0d0
      end if
      if (dbg) then
        print *,'n2ang: ',n2ang
        print *,'nW2acryo, f2ac:     ', nW2acryo, f2ac
        print *,'nW2dcryo, f2dc:     ', nW2dcryo, f2dc
        print *,'nW2apostcryo, f2ap: ', nW2apostcryo, f2ap
        print *,'nW2dpostcryo, f2dp: ', nW2dpostcryo, f2dp
      end if
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c                                      ! Read in all 8 PSF files
c
      InFnam = CanPath(1:LNBlnk(CanPath))//'w1-psf-cryo-wpro.fits'
c
      Call GetNAX(InFNam,NCols,NRows,NPlanes,1,xlo,xhi,ylo,yhi,FileID)
      if (dbg) print *,'GetNAX returned NCols,NRows,NPlanes,FileID:',
     +                  NCols,NRows,NPlanes,FileID ! dbg
      if (dbg) print *,InFNam(1:lnblnk(InFNam))
      if (dbg) print *,'xlo,xhi,ylo,yhi: ', xlo, xhi, ylo, yhi
      if ((Ncols .ne. 641) .or. (NRows .ne. 641) .or. (NPlanes .ne.1))
     +then
        print *,'ERROR: NCols, NRows, NPlanes = ',
     +           Ncols,', ',NRows,', ',NPlanes
        print *,'       in file: ',InFNam(1:lnblnk(InFNam))
        print *,'       should be 641 by 641'
        call exit(64)
      end if
c
      NPix = NCols*NRows
      naxes(1) = NCols
      naxes(2) = NRows
c
      status = 0
      call ftgpve(FileID,1,1,NPix,nullval,w1psfc,anynull,status)
      if (dbg) print *,'ftgpve returned ',NPix,' pixels'              ! dbg
c
C  The FITS file must always be closed before exiting the program. 
C  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(FileID, status)
      call ftfiou(FileID, status)
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     
      InFnam = CanPath(1:LNBlnk(CanPath))//'w1-psfunc-cryo-wpro.fits'
c
      Call GetNAX(InFNam,NCols,NRows,NPlanes,1,xlo,xhi,ylo,yhi,FileID)
      if (dbg) print *,'GetNAX returned NCols,NRows,NPlanes,FileID:',
     +                  NCols,NRows,NPlanes,FileID ! dbg
      if (dbg) print *,InFNam(1:lnblnk(InFNam))
      if (dbg) print *,'xlo,xhi,ylo,yhi: ', xlo, xhi, ylo, yhi
      if ((Ncols .ne. 641) .or. (NRows .ne. 641) .or. (NPlanes .ne.1))
     +then
        print *,'ERROR: NCols, NRows, NPlanes = ',
     +           Ncols,', ',NRows,', ',NPlanes
        print *,'       in file: ',InFNam(1:lnblnk(InFNam))
        print *,'       should be 641 by 641'
        call exit(64)
      end if
c
      status = 0
      call ftgpve(FileID,1,1,NPix,nullval,w1psfuncc,anynull,status)
      if (dbg) print *,'ftgpve returned ',NPix,' pixels'              ! dbg
c
C  The FITS file must always be closed before exiting the program. 
C  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(FileID, status)
      call ftfiou(FileID, status)
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     
      InFnam = CanPath(1:LNBlnk(CanPath))//'w1-psf-postcryo-wpro.fits'
c
      Call GetNAX(InFNam,NCols,NRows,NPlanes,1,xlo,xhi,ylo,yhi,FileID)
      if (dbg) print *,'GetNAX returned NCols,NRows,NPlanes,FileID:',
     +                  NCols,NRows,NPlanes,FileID ! dbg
      if (dbg) print *,InFNam(1:lnblnk(InFNam))
      if (dbg) print *,'xlo,xhi,ylo,yhi: ', xlo, xhi, ylo, yhi
      if ((Ncols .ne. 641) .or. (NRows .ne. 641) .or. (NPlanes .ne.1))
     +then
        print *,'ERROR: NCols, NRows, NPlanes = ',
     +           Ncols,', ',NRows,', ',NPlanes
        print *,'       in file: ',InFNam(1:lnblnk(InFNam))
        print *,'       should be 641 by 641'
        call exit(64)
      end if
c
      status = 0
      call ftgpve(FileID,1,1,NPix,nullval,w1psfp,anynull,status)
      if (dbg) print *,'ftgpve returned ',NPix,' pixels'              ! dbg
c
C  The FITS file must always be closed before exiting the program. 
C  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(FileID, status)
      call ftfiou(FileID, status)
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     
      InFnam = CanPath(1:LNBlnk(CanPath))//'w1-psfunc-postcryo-wpro.fits'
c
      Call GetNAX(InFNam,NCols,NRows,NPlanes,1,xlo,xhi,ylo,yhi,FileID)
      if (dbg) print *,'GetNAX returned NCols,NRows,NPlanes,FileID:',
     +                  NCols,NRows,NPlanes,FileID ! dbg
      if (dbg) print *,InFNam(1:lnblnk(InFNam))
      if (dbg) print *,'xlo,xhi,ylo,yhi: ', xlo, xhi, ylo, yhi
      if ((Ncols .ne. 641) .or. (NRows .ne. 641) .or. (NPlanes .ne.1))
     +then
        print *,'ERROR: NCols, NRows, NPlanes = ',
     +           Ncols,', ',NRows,', ',NPlanes
        print *,'       in file: ',InFNam(1:lnblnk(InFNam))
        print *,'       should be 641 by 641'
        call exit(64)
      end if
c
      status = 0
      call ftgpve(FileID,1,1,NPix,nullval,w1psfuncp,anynull,status)
      if (dbg) print *,'ftgpve returned ',NPix,' pixels'              ! dbg
c
C  The FITS file must always be closed before exiting the program. 
C  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(FileID, status)
      call ftfiou(FileID, status)
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      InFnam = CanPath(1:LNBlnk(CanPath))//'w2-psf-cryo-wpro.fits'
c
      Call GetNAX(InFNam,NCols,NRows,NPlanes,1,xlo,xhi,ylo,yhi,FileID)
      if (dbg) print *,'GetNAX returned NCols,NRows,NPlanes,FileID:',
     +                  NCols,NRows,NPlanes,FileID ! dbg
      if (dbg) print *,InFNam(1:lnblnk(InFNam))
      if (dbg) print *,'xlo,xhi,ylo,yhi: ', xlo, xhi, ylo, yhi
      if ((Ncols .ne. 641) .or. (NRows .ne. 641) .or. (NPlanes .ne.1))
     +then
        print *,'ERROR: NCols, NRows, NPlanes = ',
     +           Ncols,', ',NRows,', ',NPlanes
        print *,'       in file: ',InFNam(1:lnblnk(InFNam))
        print *,'       should be 641 by 641'
        call exit(64)
      end if
c
      status = 0
      call ftgpve(FileID,1,1,NPix,nullval,w2psfc,anynull,status)
      if (dbg) print *,'ftgpve returned ',NPix,' pixels'              ! dbg
c
C  The FITS file must always be closed before exiting the program. 
C  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(FileID, status)
      call ftfiou(FileID, status)
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     
      InFnam = CanPath(1:LNBlnk(CanPath))//'w2-psfunc-cryo-wpro.fits'
c
      Call GetNAX(InFNam,NCols,NRows,NPlanes,1,xlo,xhi,ylo,yhi,FileID)
      if (dbg) print *,'GetNAX returned NCols,NRows,NPlanes,FileID:',
     +                  NCols,NRows,NPlanes,FileID ! dbg
      if (dbg) print *,InFNam(1:lnblnk(InFNam))
      if (dbg) print *,'xlo,xhi,ylo,yhi: ', xlo, xhi, ylo, yhi
      if ((Ncols .ne. 641) .or. (NRows .ne. 641) .or. (NPlanes .ne.1))
     +then
        print *,'ERROR: NCols, NRows, NPlanes = ',
     +           Ncols,', ',NRows,', ',NPlanes
        print *,'       in file: ',InFNam(1:lnblnk(InFNam))
        print *,'       should be 641 by 641'
        call exit(64)
      end if
c
      status = 0
      call ftgpve(FileID,1,1,NPix,nullval,w2psfuncc,anynull,status)
      if (dbg) print *,'ftgpve returned ',NPix,' pixels'              ! dbg
c
C  The FITS file must always be closed before exiting the program. 
C  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(FileID, status)
      call ftfiou(FileID, status)
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     
      InFnam = CanPath(1:LNBlnk(CanPath))//'w2-psf-postcryo-wpro.fits'
c
      Call GetNAX(InFNam,NCols,NRows,NPlanes,1,xlo,xhi,ylo,yhi,FileID)
      if (dbg) print *,'GetNAX returned NCols,NRows,NPlanes,FileID:',
     +                  NCols,NRows,NPlanes,FileID ! dbg
      if (dbg) print *,InFNam(1:lnblnk(InFNam))
      if (dbg) print *,'xlo,xhi,ylo,yhi: ', xlo, xhi, ylo, yhi
      if ((Ncols .ne. 641) .or. (NRows .ne. 641) .or. (NPlanes .ne.1))
     +then
        print *,'ERROR: NCols, NRows, NPlanes = ',
     +           Ncols,', ',NRows,', ',NPlanes
        print *,'       in file: ',InFNam(1:lnblnk(InFNam))
        print *,'       should be 641 by 641'
        call exit(64)
      end if
c
      status = 0
      call ftgpve(FileID,1,1,NPix,nullval,w2psfp,anynull,status)
      if (dbg) print *,'ftgpve returned ',NPix,' pixels'              ! dbg
c
C  The FITS file must always be closed before exiting the program. 
C  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(FileID, status)
      call ftfiou(FileID, status)
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     
      InFnam = CanPath(1:LNBlnk(CanPath))//'w2-psfunc-postcryo-wpro.fits'
c
      Call GetNAX(InFNam,NCols,NRows,NPlanes,1,xlo,xhi,ylo,yhi,FileID)
      if (dbg) print *,'GetNAX returned NCols,NRows,NPlanes,FileID:',
     +                  NCols,NRows,NPlanes,FileID ! dbg
      if (dbg) print *,InFNam(1:lnblnk(InFNam))
      if (dbg) print *,'xlo,xhi,ylo,yhi: ', xlo, xhi, ylo, yhi
      if ((Ncols .ne. 641) .or. (NRows .ne. 641) .or. (NPlanes .ne.1))
     +then
        print *,'ERROR: NCols, NRows, NPlanes = ',
     +           Ncols,', ',NRows,', ',NPlanes
        print *,'       in file: ',InFNam(1:lnblnk(InFNam))
        print *,'       should be 641 by 641'
        call exit(64)
      end if
c
      status = 0
      call ftgpve(FileID,1,1,NPix,nullval,w2psfuncp,anynull,status)
      if (dbg) print *,'ftgpve returned ',NPix,' pixels'              ! dbg
c
C  The FITS file must always be closed before exiting the program. 
C  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(FileID, status)
      call ftfiou(FileID, status)
c      
      w1psfopt0    = 0.0
      w2psfopt0    = 0.0
      w1psfuncopt0 = 0.0
      w2psfuncopt0 = 0.0
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c                                      ! W1 ascending
c
      call Rot8(w1psfc,w1psfuncc,nAngsw1ac,dH1ac,w1angsac,nOvrSamp,
     +          Ang1acMin,w1psfca,w1psfuncca,wgt)
      call Rot8(w1psfp,w1psfuncp,nAngsw1ap,dH1ap,w1angsap,nOvrSamp,
     +          Ang1apMin,w1psfpa,w1psfuncpa,wgt)
      w1psfopt0    = w1psfopt0    + f1ac*w1psfca    + f1ap*w1psfpa
      w1psfuncopt0 = w1psfuncopt0 + f1ac*w1psfuncca + f1ap*w1psfuncpa
c
      hdrline(1)  = 'EXTEND  =                    T /'
     +              //' TAPE MAY HAVE STANDARD FITS EXTENSIONS'         
      hdrline(2)  = 'CDELT2  =          9.54861E-05 /'
      hdrline(3)  = 'CRPIX1  =              321.000 /'
      hdrline(4)  = 'CRPIX2  =              321.000 /'
      hdrline(5)  = 'XLO     =                1.000 /'
      hdrline(6)  = 'YLO     =                1.000 /'
      hdrline(7)  = 'XHI     =             2048.000 /'
      hdrline(8)  = 'YHI     =             2048.000 /'
      write(hdrline(9), '(''NOVRSAMP=           '',I10,'' /'')')
     +        nOvrSamp
      hdrline(12) = 'TILE    = '''//TileNam//''''
      hdrline(14) = 'PSFTYPE = ''W1'''
c         
      OutFNam = OutPath(1:LNBlnk(OutPath))
     +          //'Asce'//Slash//'unWISE-w1-psf-wpro-01x01-01x01.fits'
c
      if (W3) then                     ! W1 cryo --> W3
        OutFNam = OutPath(1:LNBlnk(OutPath))
     +            //'Asce'//Slash//'unWISE-w3-psf-wpro-01x01-01x01.fits'
        hdrline(14) = 'PSFTYPE = ''W1 CRYO LABELED W3'''
      else if ((nW1acryo .gt. 0) .and. (nW1apostcryo .gt. 0)) then
        wgt1 = dfloat(nW1acryo)/dfloat(nW1acryo + nW1apostcryo)
        wgt2 = 1.0d0 - wgt1
        w1psfuncca = sqrt((wgt1*w1psfuncca)**2 + (wgt2*w1psfuncpa)**2
     +                + (w1psfca-w1psfpa)**2)
        w1psfca = wgt1*w1psfca + wgt2*w1psfpa
        PixSum = 0.0d0
        do 70 j = 1, 641
          do 60 i = 1, 641
            PixSum = PixSum + w1psfca(i,j)
60        continue
70      continue        
        do 90 j = 1, 641
          do 80 i = 1, 641
            w1psfca(i,j) = w1psfca(i,j)/PixSum
80        continue
90      continue        
        if (Ang1apMin .lt. Ang1acMin) Ang1acMin = Ang1apMin
        if (Ang1apMax .gt. Ang1acMax) Ang1acMax = Ang1apMax
        hdrline(14) = 'PSFTYPE = ''W1 CRYO+POSTCRYO'''
      else if (nW1apostcryo .gt. 0) then
        w1psfca    = w1psfpa
        w1psfuncca = w1psfuncpa
        Ang1acMin = Ang1apMin
        Ang1acMax = Ang1apMax
        hdrline(14) = 'PSFTYPE = ''W1 POSTCRYO'''
      end if
c      
      write(hdrline(10),'(''ECLANGMN=           '',f10.4,'' /'')')
     +      Ang1acMin
      write(hdrline(11),'(''ECLANGMX=           '',f10.4,'' /'')')
     +      Ang1acMax
      hdrline(13) = 'SCANDIR = ''ASCENDING'''
      write(hdrline(15),'(''FPSFUNC =           '',f10.4,'' /'')') s11
c
      call wrfitsr(2,naxes,NPix,w1psfca,OutFNam,status,15,hdrline)
      if (status .ne.0) then
        print *, 'ERROR: output status = ',status,' on file:'
        print *, OutFNAM(1:lnblnk(OutFNam))
      else        
        print *, OutFNam(1:lnblnk(OutFNam)),' generated'
      end if
c
      if (W3) then
        OutFNam = OutPath(1:LNBlnk(OutPath))
     +            //'Asce'//Slash//'unWISE-w3-psfunc-wpro-01x01-01x01.fits'
        hdrline(14) = 'PSFTYPE = ''W1 CRYO LABELED W3'''
      else
        OutFNam = OutPath(1:LNBlnk(OutPath))
     +            //'Asce'//Slash//'unWISE-w1-psfunc-wpro-01x01-01x01.fits'
        hdrline(14) = 'PSFTYPE = ''W1 CRYO+POSTCRYO'''
      end if
c
      w1psfuncca = s11*w1psfuncca
      call wrfitsr(2,naxes,NPix,w1psfuncca,OutFNam,status,15,hdrline)
      if (status .ne.0) then
        print *,'ERROR: output status = ',status,' on file:'
        print *,OutFNAM(1:lnblnk(OutFNam))
      else        
        print *, OutFNam(1:lnblnk(OutFNam)),' generated'
      end if
c      
      if (W3) then                     ! W1 post-cryo as W1
        OutFNam = OutPath(1:LNBlnk(OutPath))
     +            //'Asce'//Slash//'unWISE-w1-psf-wpro-01x01-01x01.fits'
        write(hdrline(10),'(''ECLANGMN=           '',f10.4,'' /'')')
     +        Ang1apMin
        write(hdrline(11),'(''ECLANGMX=           '',f10.4,'' /'')')
     +        Ang1apMax
        hdrline(14) = 'PSFTYPE = ''W1 POSTCRYO'''
        call wrfitsr(2,naxes,NPix,w1psfpa,OutFNam,status,15,hdrline)
        if (status .ne.0) then
          print *,'ERROR: output status = ',status,' on file:'
          print *,OutFNAM(1:lnblnk(OutFNam))
        else        
          print *, OutFNam(1:lnblnk(OutFNam)),' generated'
        end if
        OutFNam = OutPath(1:LNBlnk(OutPath))
     +            //'Asce'//Slash//'unWISE-w1-psfunc-wpro-01x01-01x01.fits'
        w1psfuncpa = s11*w1psfuncpa
        call wrfitsr(2,naxes,NPix,w1psfuncpa,OutFNam,status,15,hdrline)
        if (status .ne.0) then
          print *,'ERROR: output status = ',status,' on file:'
          print *,OutFNAM(1:lnblnk(OutFNam))
        else        
          print *, OutFNam(1:lnblnk(OutFNam)),' generated'
        end if
      end if
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c                                      ! debug check on counters
      if (dbg) then      
        OutFNam = OutPath(1:LNBlnk(OutPath))//'w1a-counters.fits'
        call wrfitsr(2,naxes,NPix,wgt,OutFNam,status,14,hdrline)
        if (status .ne.0) then
          print *,'ERROR: output status = ',status,' on file:'
          print *,OutFNAM(1:lnblnk(OutFNam))
        else        
          print *, OutFNam(1:lnblnk(OutFNam)),' generated'
        end if
      end if
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c                                      ! W1 descending
c
      call Rot8(w1psfc,w1psfuncc,nAngsw1dc,dH1dc,w1angsdc,nOvrSamp,
     +          Ang1dcMin,w1psfcd,w1psfunccd,wgt)
      call Rot8(w1psfp,w1psfuncp,nAngsw1dp,dH1dp,w1angsdp,nOvrSamp,
     +          Ang1dpMin,w1psfpd,w1psfuncpd,wgt)
      w1psfopt0    = w1psfopt0    + f1dc*w1psfcd    + f1dp*w1psfpd
      w1psfuncopt0 = w1psfuncopt0 + f1dc*w1psfunccd + f1dp*w1psfuncpd
c         
      OutFNam = OutPath(1:LNBlnk(OutPath))
     +          //'Desc'//Slash//'unWISE-w1-psf-wpro-01x01-01x01.fits'
c
      hdrline(14) = 'PSFTYPE = ''W1'''
      if (W3) then                     ! W1 cryo --> W3
        OutFNam = OutPath(1:LNBlnk(OutPath))
     +            //'Desc'//Slash//'unWISE-w3-psf-wpro-01x01-01x01.fits'
        hdrline(14) = 'PSFTYPE = ''W1 CRYO LABELED W3'''
      else if ((nW1dcryo .gt. 0) .and. (nW1dpostcryo .gt. 0)) then
        wgt1 = dfloat(nW1dcryo)/dfloat(nW1dcryo + nW1dpostcryo)
        wgt2 = 1.0d0 - wgt1
        w1psfunccd = sqrt((wgt1*w1psfunccd)**2 + (wgt2*w1psfuncpd)**2
     +                + (w1psfcd-w1psfpd)**2)
        w1psfcd = wgt1*w1psfcd + wgt2*w1psfpd
        PixSum = 0.0d0
        do 170 j = 1, 641
          do 160 i = 1, 641
            PixSum = PixSum + w1psfcd(i,j)
160       continue
170     continue        
        do 190 j = 1, 641
          do 180 i = 1, 641
            w1psfcd(i,j) = w1psfcd(i,j)/PixSum
180       continue
190     continue        
        if (Ang1dpMin .lt. Ang1dcMin) Ang1dcMin = Ang1dpMin
        if (Ang1dpMax .gt. Ang1dcMax) Ang1dcMax = Ang1dpMax
        hdrline(14) = 'PSFTYPE = ''W1 CRYO+POSTCRYO'''
      else if (nW1dpostcryo .gt. 0) then
        w1psfcd    = w1psfpd
        w1psfunccd = w1psfuncpd
        Ang1dcMin = Ang1dpMin
        Ang1dcMax = Ang1dpMax
        hdrline(14) = 'PSFTYPE = ''W1 POSTCRYO'''
      end if
c      
      write(hdrline(10),'(''ECLANGMN=           '',f10.4,'' /'')')
     +      Ang1dcMin
      write(hdrline(11),'(''ECLANGMX=           '',f10.4,'' /'')')
     +      Ang1dcMax
      hdrline(13) = 'SCANDIR = ''DESCENDING'''
      write(hdrline(15),'(''FPSFUNC =           '',f10.4,'' /'')') s11
      call wrfitsr(2,naxes,NPix,w1psfcd,OutFNam,status,15,hdrline)
      if (status .ne.0) then
        print *,'ERROR: output status = ',status,' on file:'
        print *,OutFNAM(1:lnblnk(OutFNam))
      else        
        print *, OutFNam(1:lnblnk(OutFNam)),' generated'
      end if
c
      if (W3) then
        OutFNam = OutPath(1:LNBlnk(OutPath))
     +            //'Desc'//Slash//'unWISE-w3-psfunc-wpro-01x01-01x01.fits'
      else
        OutFNam = OutPath(1:LNBlnk(OutPath))
     +            //'Desc'//Slash//'unWISE-w1-psfunc-wpro-01x01-01x01.fits'
      end if
c
      w1psfunccd = s11*w1psfunccd
      call wrfitsr(2,naxes,NPix,w1psfunccd,OutFNam,status,15,hdrline)
      if (status .ne.0) then
        print *,'ERROR: output status = ',status,' on file:'
        print *,OutFNAM(1:lnblnk(OutFNam))
      else        
        print *, OutFNam(1:lnblnk(OutFNam)),' generated'
      end if
c      
      if (W3) then                     ! W1 post-cryo as W1
        OutFNam = OutPath(1:LNBlnk(OutPath))
     +            //'Desc'//Slash//'unWISE-w1-psf-wpro-01x01-01x01.fits'
        write(hdrline(10),'(''ECLANGMN=           '',f10.4,'' /'')')
     +        Ang1dpMin
        write(hdrline(11),'(''ECLANGMX=           '',f10.4,'' /'')')
     +        Ang1dpMax
        hdrline(14) = 'PSFTYPE = ''W1 POSTCRYO'''
        call wrfitsr(2,naxes,NPix,w1psfpd,OutFNam,status,15,hdrline)
        if (status .ne.0) then
          print *,'ERROR: output status = ',status,' on file:'
          print *,OutFNAM(1:lnblnk(OutFNam))
        else        
          print *, OutFNam(1:lnblnk(OutFNam)),' generated'
        end if
        OutFNam = OutPath(1:LNBlnk(OutPath))
     +            //'Desc'//Slash//'unWISE-w1-psfunc-wpro-01x01-01x01.fits'
        w1psfuncpd = s11*w1psfuncpd
        call wrfitsr(2,naxes,NPix,w1psfuncpd,OutFNam,status,15,hdrline)
        if (status .ne.0) then
          print *,'ERROR: output status = ',status,' on file:'
          print *,OutFNAM(1:lnblnk(OutFNam))
        else        
          print *, OutFNam(1:lnblnk(OutFNam)),' generated'
        end if
      end if
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c                                      ! W2 ascending
c
      call Rot8(w2psfc,w2psfuncc,nAngsw2ac,dH2ac,w2angsac,nOvrSamp,
     +          Ang2acMin,w2psfca,w2psfuncca,wgt)
      call Rot8(w2psfp,w2psfuncp,nAngsw2ap,dH2ap,w2angsap,nOvrSamp,
     +          Ang2apMin,w2psfpa,w2psfuncpa,wgt)
      w2psfopt0    = w2psfopt0    + f2ac*w2psfca    + f2ap*w2psfpa
      w2psfuncopt0 = w2psfuncopt0 + f2ac*w2psfuncca + f2ap*w2psfuncpa
c      
      OutFNam = OutPath(1:LNBlnk(OutPath))
     +          //'Asce'//Slash//'unWISE-w2-psf-wpro-01x01-01x01.fits'
c
      hdrline(14) = 'PSFTYPE = ''W2'''
      if (W4) then                     ! W2 cryo as W4
        OutFNam = OutPath(1:LNBlnk(OutPath))
     +            //'Asce'//Slash//'unWISE-w4-psf-wpro-01x01-01x01.fits'
        hdrline(14) = 'PSFTYPE = ''W2 CRYO LABELED W4'''
      else if ((nW2acryo .gt. 0) .and. (nW2apostcryo .gt. 0)) then
        wgt1 = dfloat(nW2acryo)/dfloat(nW2acryo + nW2apostcryo)
        wgt2 = 1.0d0 - wgt1
        w2psfuncca = sqrt((wgt1*w2psfuncca)**2 + (wgt2*w2psfuncpa)**2
     +                + (w2psfca-w2psfpa)**2)
        w2psfca = wgt1*w2psfca + wgt2*w2psfpa
        PixSum = 0.0d0
        hdrline(14) = 'PSFTYPE = ''W2'''
        do 270 j = 1, 641
          do 260 i = 1, 641
            PixSum = PixSum + w2psfca(i,j)
260       continue
270     continue        
        if (dbg) print *,'W2 PSF Asce pixel sum before renorm: ',PixSum
        do 290 j = 1, 641
          do 280 i = 1, 641
            w2psfca(i,j) = w2psfca(i,j)/PixSum
280       continue
290     continue        
        if (dbg) then
          PixSum = 0.0d0
          do 310 j = 1, 641
            do 300 i = 1, 641
              PixSum = PixSum + w2psfca(i,j)
300         continue
310       continue        
          print *,'W2 PSF Asce pixel sum after renorm:  ',PixSum
        end if
        if (Ang2apMin .lt. Ang2acMin) Ang2acMin = Ang2apMin
        if (Ang2apMax .gt. Ang2acMax) Ang2acMax = Ang2apMax
         hdrline(14) = 'PSFTYPE = ''W2 CRYO+POSTCRYO'''
      else if (nW2apostcryo .gt. 0) then
        w2psfca    = w2psfpa
        w2psfuncca = w2psfuncpa
        Ang2acMin = Ang2apMin
        Ang2acMax = Ang2apMax
        hdrline(14) = 'PSFTYPE = ''W2 POSTCRYO'''
       end if
      write(hdrline(10),'(''ECLANGMN=           '',f10.4,'' /'')')
     +      Ang2acMin
      write(hdrline(11),'(''ECLANGMX=           '',f10.4,'' /'')')
     +      Ang2acMax
c
      hdrline(13) = 'SCANDIR = ''ASCENDING'''
      write(hdrline(15),'(''FPSFUNC =           '',f10.4,'' /'')') s12
      call wrfitsr(2,naxes,NPix,w2psfca,OutFNam,status,15,hdrline)
      if (status .ne.0) then
        print *,'ERROR: output status = ',status,' on file:'
        print *,OutFNAM(1:lnblnk(OutFNam))
      else        
        print *, OutFNam(1:lnblnk(OutFNam)),' generated'
      end if
c
      if (W4) then
        OutFNam = OutPath(1:LNBlnk(OutPath))
     +            //'Asce'//Slash//'unWISE-w4-psfunc-wpro-01x01-01x01.fits'
      else
        OutFNam = OutPath(1:LNBlnk(OutPath))
     +            //'Asce'//Slash//'unWISE-w2-psfunc-wpro-01x01-01x01.fits'
      end if
c
      w2psfuncca = s12*w2psfuncca
      call wrfitsr(2,naxes,NPix,w2psfuncca,OutFNam,status,15,hdrline)
      if (status .ne.0) then
        print *,'ERROR: output status = ',status,' on file:'
        print *,OutFNAM(1:lnblnk(OutFNam))
      else        
        print *, OutFNam(1:lnblnk(OutFNam)),' generated'
      end if
c      
      if (W4) then                     ! W2 post-cryo as W2
        OutFNam = OutPath(1:LNBlnk(OutPath))
     +            //'Asce'//Slash//'unWISE-w2-psf-wpro-01x01-01x01.fits'
        write(hdrline(10),'(''ECLANGMN=           '',f10.4,'' /'')')
     +        Ang2apMin
        write(hdrline(11),'(''ECLANGMX=           '',f10.4,'' /'')')
     +        Ang2apMax
        hdrline(14) = 'PSFTYPE = ''W2 POSTCRYO'''
        call wrfitsr(2,naxes,NPix,w2psfpa,OutFNam,status,15,hdrline)
        if (status .ne.0) then
          print *,'ERROR: output status = ',status,' on file:'
          print *,OutFNAM(1:lnblnk(OutFNam))
        else        
          print *, OutFNam(1:lnblnk(OutFNam)),' generated'
        end if
        OutFNam = OutPath(1:LNBlnk(OutPath))
     +            //'Asce'//Slash//'unWISE-w2-psfunc-wpro-01x01-01x01.fits'
        w2psfuncpa = s12*w2psfuncpa
        call wrfitsr(2,naxes,NPix,w2psfuncpa,OutFNam,status,15,hdrline)
        if (status .ne.0) then
          print *,'ERROR: output status = ',status,' on file:'
          print *,OutFNAM(1:lnblnk(OutFNam))
        else        
          print *, OutFNam(1:lnblnk(OutFNam)),' generated'
        end if
      end if
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c                                      ! W2 descending
c
      call Rot8(w2psfc,w2psfuncc,nAngsw2dc,dH2dc,w2angsdc,nOvrSamp,
     +          Ang2dcMin,w2psfcd,w2psfunccd,wgt)
      call Rot8(w2psfp,w2psfuncp,nAngsw2dp,dH2dp,w2angsdp,nOvrSamp,
     +          Ang2dpMin,w2psfpd,w2psfuncpd,wgt)
      w2psfopt0    = w2psfopt0    + f2dc*w2psfcd    + f2dp*w2psfpd
      w2psfuncopt0 = w2psfuncopt0 + f2dc*w2psfunccd + f2dp*w2psfuncpd
c      
      OutFNam = OutPath(1:LNBlnk(OutPath))
     +          //'Desc'//Slash//'unWISE-w2-psf-wpro-01x01-01x01.fits'
c
      hdrline(14) = 'PSFTYPE = ''W2'''
      if (W4) then
        OutFNam = OutPath(1:LNBlnk(OutPath))
     +            //'Desc'//Slash//'unWISE-w4-psf-wpro-01x01-01x01.fits'
        hdrline(14) = 'PSFTYPE = ''W2 CRYO LABELED W4'''
      else if ((nW2dcryo .gt. 0) .and. (nW2dpostcryo .gt. 0)) then
        wgt1 = dfloat(nW2dcryo)/dfloat(nW2dcryo + nW2dpostcryo)
        wgt2 = 1.0d0 - wgt1
        w2psfunccd = sqrt((wgt1*w2psfunccd)**2 + (wgt2*w2psfuncpd)**2
     +                + (w2psfcd-w2psfpd)**2)
        w2psfcd = wgt1*w2psfcd + wgt2*w2psfpd
        PixSum = 0.0d0
        do 370 j = 1, 641
          do 360 i = 1, 641
            PixSum = PixSum + w2psfcd(i,j)
360       continue
370     continue        
        do 390 j = 1, 641
          do 380 i = 1, 641
            w2psfcd(i,j) = w2psfcd(i,j)/PixSum
380       continue
390     continue        
        if (Ang2dpMin .lt. Ang2dcMin) Ang2dcMin = Ang2dpMin
        if (Ang2dpMax .gt. Ang2dcMax) Ang2dcMax = Ang2dpMax
        hdrline(14) = 'PSFTYPE = ''W2 CRYO+POSTCRYO'''
       else if (nW2dpostcryo .gt. 0) then
        w2psfcd    = w2psfpd
        w2psfunccd = w2psfuncpd
        Ang2dcMin = Ang2dpMin
        Ang2dcMax = Ang2dpMax
        hdrline(14) = 'PSFTYPE = ''W2 POSTCRYO'''
       end if
      write(hdrline(10),'(''ECLANGMN=           '',f10.4,'' /'')')
     +      Ang2dcMin
      write(hdrline(11),'(''ECLANGMX=           '',f10.4,'' /'')')
     +      Ang2dcMax
c
      hdrline(13) = 'SCANDIR = ''DESCENDING'''
      write(hdrline(15),'(''FPSFUNC =           '',f10.4,'' /'')') s12
      call wrfitsr(2,naxes,NPix,w2psfcd,OutFNam,status,15,hdrline)
      if (status .ne.0) then
        print *,'ERROR: output status = ',status,' on file:'
        print *,OutFNAM(1:lnblnk(OutFNam))
      else        
        print *, OutFNam(1:lnblnk(OutFNam)),' generated'
      end if
c
      if (W4) then
        OutFNam = OutPath(1:LNBlnk(OutPath))
     +            //'Desc'//Slash//'unWISE-w4-psfunc-wpro-01x01-01x01.fits'
      else
        OutFNam = OutPath(1:LNBlnk(OutPath))
     +            //'Desc'//Slash//'unWISE-w2-psfunc-wpro-01x01-01x01.fits'
      end if
c
      w2psfunccd = s12*w2psfunccd
      call wrfitsr(2,naxes,NPix,w2psfunccd,OutFNam,status,15,hdrline)
      if (status .ne.0) then
        print *,'ERROR: output status = ',status,' on file:'
        print *,OutFNAM(1:lnblnk(OutFNam))
      else        
        print *, OutFNam(1:lnblnk(OutFNam)),' generated'
      end if   
c      
      if (W4) then                     ! W2 post-cryo as W2
        OutFNam = OutPath(1:LNBlnk(OutPath))
     +            //'Desc'//Slash//'unWISE-w2-psf-wpro-01x01-01x01.fits'
        write(hdrline(10),'(''ECLANGMN=           '',f10.4,'' /'')')
     +        Ang2dpMin
        write(hdrline(11),'(''ECLANGMX=           '',f10.4,'' /'')')
     +        Ang2dpMax
        hdrline(14) = 'PSFTYPE = ''W2 POSTCRYO'''
        write(hdrline(15),'(''FPSFUNC =           '',f10.4,'' /'')') s12
        call wrfitsr(2,naxes,NPix,w2psfpd,OutFNam,status,15,hdrline)
        if (status .ne.0) then
          print *,'ERROR: output status = ',status,' on file:'
          print *,OutFNAM(1:lnblnk(OutFNam))
        else        
          print *, OutFNam(1:lnblnk(OutFNam)),' generated'
        end if
        OutFNam = OutPath(1:LNBlnk(OutPath))
     +            //'Desc'//Slash//'unWISE-w2-psfunc-wpro-01x01-01x01.fits'
        w2psfuncpd = s12*w2psfuncpd
        call wrfitsr(2,naxes,NPix,w2psfuncpd,OutFNam,status,15,hdrline)
        if (status .ne.0) then
          print *,'ERROR: output status = ',status,' on file:'
          print *,OutFNAM(1:lnblnk(OutFNam))
        else        
          print *, OutFNam(1:lnblnk(OutFNam)),' generated'
        end if
      end if
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c                                      ! Option 0 PSFs
      OutFNam = OutPath(1:LNBlnk(OutPath))
     +          //'unWISE-w1-psf-wpro-01x01-01x01.fits'
      if (Ang1apMin .lt. Ang1acMin) Ang1acMin = Ang1apMin
      if (Ang1apMax .gt. Ang1acMax) Ang1acMax = Ang1apMax
      if ((abs(Ang1acMin) .le. 360.0) .and. (abs(Ang1acMax) .le. 360.0)) then
        write(hdrline(10),'(''ECANGAMN=           '',f10.4,'' /'')')
     +        Ang1acMin
        write(hdrline(11),'(''ECANGAMX=           '',f10.4,'' /'')')
     +        Ang1acMax
      else
        hdrline(10) = 'ECANGAMN= null'
        hdrline(11) = 'ECANGAMX= null'
      end if
      if (Ang1dpMin .lt. Ang1dcMin) Ang1dcMin = Ang1dpMin
      if (Ang1dpMax .gt. Ang1dcMax) Ang1dcMax = Ang1dpMax
      if ((abs(Ang1dpMin) .le. 360.0) .and. (abs(Ang1dpMax) .le. 360.0)) then
        write(hdrline(12),'(''ECANGDMN=           '',f10.4,'' /'')')
     +        Ang1dpMin
        write(hdrline(13),'(''ECANGDMX=           '',f10.4,'' /'')')
     +        Ang1dpMax
      else
        hdrline(12) = 'ECANGDMN= null'
        hdrline(13) = 'ECANGDMX= null'
      end if
      hdrline(14) = 'PSFTYPE = ''W1 OPTION 0'''
      hdrline(15) = 'TILE    = '''//TileNam//''''
      hdrline(16) = 'SCANDIR = ''ASCENDING AND DESCENDING AVERAGED'''
      write(hdrline(17),'(''FPSFUNC =           '',f10.4,'' /'')') s01
      call wrfitsr(2,naxes,NPix,w1psfopt0,OutFNam,status,17,hdrline)
      if (status .ne.0) then
        print *,'ERROR: output status = ',status,' on file:'
        print *,OutFNAM(1:lnblnk(OutFNam))
      else        
        print *, OutFNam(1:lnblnk(OutFNam)),' generated'
      end if
      OutFNam = OutPath(1:LNBlnk(OutPath))
     +          //'unWISE-w1-psfunc-wpro-01x01-01x01.fits'
      w1psfuncopt0 = s01*w1psfuncopt0
      call wrfitsr(2,naxes,NPix,w1psfuncopt0,OutFNam,status,17,hdrline)
      if (status .ne.0) then
        print *,'ERROR: output status = ',status,' on file:'
        print *,OutFNAM(1:lnblnk(OutFNam))
      else        
        print *, OutFNam(1:lnblnk(OutFNam)),' generated'
      end if
c
      OutFNam = OutPath(1:LNBlnk(OutPath))
     +          //'unWISE-w2-psf-wpro-01x01-01x01.fits'
      if (Ang2apMin .lt. Ang2acMin) Ang2acMin = Ang2apMin
      if (Ang2apMax .gt. Ang2acMax) Ang2acMax = Ang2apMax
      if ((abs(Ang2acMin) .le. 360.0) .and. (abs(Ang2acMax) .le. 360.0)) then
        write(hdrline(10),'(''ECANGAMN=           '',f10.4,'' /'')')
     +        Ang2acMin
        write(hdrline(11),'(''ECANGAMX=           '',f10.4,'' /'')')
     +        Ang2acMax
      else
        hdrline(10) = 'ECANGAMN= null'
        hdrline(11) = 'ECANGAMX= null'
      end if
      if (Ang2dpMin .lt. Ang2dcMin) Ang2dcMin = Ang2dpMin
      if (Ang2dpMax .gt. Ang2dcMax) Ang2dcMax = Ang2dpMax
      if ((abs(Ang2dcMin) .le. 360.0) .and. (abs(Ang2dcMax) .le. 360.0)) then
        write(hdrline(12),'(''ECANGDMN=           '',f10.4,'' /'')')
     +        Ang2dcMin
        write(hdrline(13),'(''ECANGDMX=           '',f10.4,'' /'')')
     +        Ang2dcMax
      else
        hdrline(12) = 'ECANGDMN= null'
        hdrline(13) = 'ECANGDMX= null'
      end if
      hdrline(14) = 'PSFTYPE = ''W2 OPTION 0'''
      write(hdrline(17),'(''FPSFUNC =           '',f10.4,'' /'')') s02
      call wrfitsr(2,naxes,NPix,w2psfopt0,OutFNam,status,17,hdrline)
      if (status .ne.0) then
        print *,'ERROR: output status = ',status,' on file:'
        print *,OutFNAM(1:lnblnk(OutFNam))
      else        
        print *, OutFNam(1:lnblnk(OutFNam)),' generated'
      end if
      OutFNam = OutPath(1:LNBlnk(OutPath))
     +          //'unWISE-w2-psfunc-wpro-01x01-01x01.fits'
      w2psfuncopt0 = s02*w2psfuncopt0
      call wrfitsr(2,naxes,NPix,w2psfuncopt0,OutFNam,status,17,hdrline)
      if (status .ne.0) then
        print *,'ERROR: output status = ',status,' on file:'
        print *,OutFNAM(1:lnblnk(OutFNam))
      else        
        print *, OutFNam(1:lnblnk(OutFNam)),' generated'
      end if
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c                                      ! AWAIC PSFs
c     CRPIX1  =                   14
c     CRPIX2  =                   14
c     CTYPE1  = 'RA---TAN'          
c     CTYPE2  = 'DEC--TAN'          
c     CDELT1  =       -0.00019097222    
c     CDELT2  =        0.00019097222    ! 2*9.54861E-05
c     FPALOCX =                  508    
c     FPALOCY =                  508    
c
      call DownSamp(w1psfopt0,nOvrSamp,w1psfawaic)
      call DownSamp(w2psfopt0,nOvrSamp,w2psfawaic)
c     
      hdrline(1)  = 'EXTEND  =                    T /'
     +              //' TAPE MAY HAVE STANDARD FITS EXTENSIONS'         
      hdrline(2)  = 'CRPIX1  =               14.000 /'
      hdrline(3)  = 'CRPIX2  =               14.000 /'
      hdrline(4)  = 'CTYPE1  = ''RA---TAN'''
      hdrline(5)  = 'CTYPE2  = ''DEC--TAN'''
      hdrline(6)  = 'CDELT1  =       -0.00019097222 /'
      hdrline(7)  = 'CDELT2  =        0.00019097222 /'
      hdrline(8)  = 'FPALOCX =                  508 /'
      hdrline(9)  = 'FPALOCY =                  508 /'
      write(hdrline(10), '(''NOVRSAMP=           '',I10,'' /'')')
     +      nOvrSamp
      hdrline(11) = 'TILE    = '''//TileNam//''''
      naxes(1) = 27
      naxes(2) = 27
      NPix     = 729
      OutFNam  = OutPath(1:LNBlnk(OutPath))
     +           //'unwise-w1-psf-awaic-01x01-01x01.fits'
      call wrfitsr(2,naxes,NPix,w1psfawaic,OutFNam,status,11,hdrline)
      if (status .ne.0) then
        print *,'ERROR: output status = ',status,' on file:'
        print *,OutFNAM(1:lnblnk(OutFNam))
      else        
        print *, OutFNam(1:lnblnk(OutFNam)),' generated'
      end if
      OutFNam  = OutPath(1:LNBlnk(OutPath))
     +           //'unwise-w2-psf-awaic-01x01-01x01.fits'
      call wrfitsr(2,naxes,NPix,w2psfawaic,OutFNam,status,11,hdrline)
      if (status .ne.0) then
        print *,'ERROR: output status = ',status,' on file:'
        print *,OutFNAM(1:lnblnk(OutFNam))
      else        
        print *, OutFNam(1:lnblnk(OutFNam)),' generated'
      end if
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     
      call signoff('PSFavr8')
c
      Stop
c
3000  print *,'ERROR in numeric string: ',TmpStr
      call exit(64)      
c
3001  print *,'ERROR in numeric string in w1pa file'
      call exit(64)      
c
3002  print *,'ERROR in numeric string in w2pa file'
      call exit(64)      
c
      end
c
c=======================================================================
c
      Subroutine GetNAX(FilNam,NAXIS1,NAXIS2,NAXIS3,IDOp,
     +                  xlo,xhi,ylo,yhi,unit)
c-----------------------------------------------------------------------
c
c    Gets the dimensions of the file named FilNam; if IDOp = 0, the
c    file is closed, otherwise it's left open with handle unit
c
c-----------------------------------------------------------------------
c
      Character*200 FilNam
      character*80 comment
      Integer*4     NAXIS, NAXIS1, NAXIS2, NAXIS3, IStat, ImOpen,
     +              ImRKeyI, unit, IDOp, LNBlnk, ImClose
      integer status,readwrite,blocksize,naxes(3)
      Real*4        xlo, xhi, ylo, yhi
      logical*4 dbg
      data dbg/.false./
c
c-----------------------------------------------------------------------
c
C  The STATUS parameter must always be initialized.
      status=0
C  Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)
      if (dbg) print *,'GetNAX: unit, status:',unit,status ! dbg
C  Open the FITS file 
      readwrite=0
      call ftopen(unit,FilNam,readwrite,blocksize,status)
      if (status /= 0) then
          write(6,'(a)') 'GetNAX: Could not read '//trim(FilNam)
          istat = 3
          return
      endif
c
C  Determine the size of the image.
      call ftgknj(unit,'NAXIS',1,2,naxes,NAXIS,status)
      if (dbg) print *,'GetNAX: naxes, naxis, status:',naxes, naxis, status ! dbg
c
C  Check that it found both NAXIS1 and NAXIS2 keywords.
      if (NAXIS .lt. 2)then
          print *,'GetNAX: Failed to read the NAXISn keywords.'
        istat = 4
          return
       end if
c
      NAXIS1 = naxes(1)
      NAXIS2 = naxes(2)
c
      If (NAXIS .gt. 2) then
        NAXIS3 = naxes(3)
      Else
        NAXIS3 = 1
      End If
c
      status = 0
      if (dbg) print *,'about to call ftgkye(unit, ''XLO'', xlo, comment, status)'
      call ftgkye(unit, 'XLO', xlo, comment, status)
      if (dbg) print *,'XLO: ',xlo      
      if (status .ne. 0) then 
        print *,'ERROR reading XLO, XHI, YLO, YHI; status = ',status
        call exit(64)
      end if
      call ftgkye(unit, 'XHI', xhi, comment, status) 
      if (dbg) print *,'XHI: ',xhi     
      if (status .ne. 0) then 
        print *,'ERROR reading XLO, XHI, YLO, YHI; status = ',status
        call exit(64)
      end if
      call ftgkye(unit, 'YLO', ylo, comment, status) 
      if (dbg) print *,'YLO: ',ylo      
      if (status .ne. 0) then 
        print *,'ERROR reading XLO, XHI, YLO, YHI; status = ',status
        call exit(64)
      end if
      call ftgkye(unit, 'YHI', yhi, comment, status)
      if (dbg) print *,'YHI: ',yhi      
      if (status .ne. 0) then 
        print *,'ERROR reading XLO, XHI, YLO, YHI; status = ',status
        call exit(64)
      end if
c
      If (IDOp .eq. 0) then
        call ftclos(unit, status)
        call ftfiou(unit, status)
      end if
c
      return
c
      end
c
c=======================================================================
c
      Subroutine NextNarg(NArg,NArgs)
c
      integer NArg, NArgs
c
c-----------------------------------------------------------------------
c
      if (NArg .lt. NArgs) then
        NArg = NArg + 1
        return
      else
        print *,'ERROR: expected another argument but none found'
        call exit(64)
      end if
      return
      end
c
c=======================================================================
c
      subroutine upcase(string)
      character*(*) string
      integer*4 j, lnblnk
c
      do 10 j = 1,lnblnk(string)
         if(string(j:j) .ge. "a" .and. string(j:j) .le. "z") then
            string(j:j) = achar(iachar(string(j:j)) - 32)
         end if
10    continue
      return
      end
c
c=======================================================================
c
      subroutine wrfitsr(naxis,naxes,lsize,array,fout,status,nhdrline,hdrline)
c-----------------------------------------------------------------------
c
c  General-purpose real*4 FITS output subroutine
c
c    NOTE (JWF B60826): this version of FITSIO *MUST* have the first few
c                       header parameters in a certain order or else it 
c                       will fail on a read attempt; the following has worked:
c
c SIMPLE  =                    T / Written by IDL:  Fri Mar 25 12:38:25 2011      
c BITPIX  =                  -32 / Number of bits per data pixel                  
c NAXIS   =                    2 / Number of data axes                            
c NAXIS1  =                  641 /Number of positions along axis 1                
c NAXIS2  =                  641 /Number of positions along axis 2                
c  ...
c END                                                                             
c
c-----------------------------------------------------------------------
c
      integer*4     status,blocksize,outunit,naxis,naxes(naxis),bitpix,
     +              fpixel, group, nelements, nhdrline, n, lsize,
     +              Access, system
      real(4)       array(lsize)
      character*(*) fout
      character(80) hdrline(nhdrline)
      logical*4     simple, extend
c
      data          blocksize/1/, bitpix/-32/, simple/.true./,
     +              extend/.false./, fpixel/1/, group/1/
c
c-----------------------------------------------------------------------
c
c                                      ! delete existing file with this
c                                      ! same name, if any
c
      if (Access(fout(1:lnblnk(fout)),' ') .eq. 0)
     +    n = system('rm '//fout(1:lnblnk(fout)))
c
c  The STATUS parameter must be initialized before using FITSIO.  A
c  positive value of STATUS is returned whenever a serious error occurs.
c  FITSIO uses an `inherited status' convention, which means that if a
c  subroutine is called with a positive input value of STATUS, then the
c  subroutine will exit immediately, preserving the status value. For 
c  simplicity, this program only checks the status value at the end of 
c  the program, but it is usually better practice to check the status 
c  value more frequently.
c     print *,'wrfitsr: file name = ',fout(1:lnblnk(fout))    ! dbg
c     do 10 n = 1, nhdrline                                   ! dbg
c       print *,'wrfitsr: header line ',n                     ! dbg
c       print *,hdrline(n)(1:lnblnk(hdrline(n)))              ! dbg
c10   continue                                                ! dbg
c
C  Get  unused Logical Unit Numbers to use to open the FITS files.
      status = 0 
      call ftgiou(outunit,status)
      if (status .ne. 0) return
c
      call ftinit(outunit,fout,blocksize,status)
      if (status .ne. 0) return
c
C  Initialize parameters about the FITS image.
C  BITPIX = 16 means that the image pixels will consist of 16-bit
C  integers.  The size of the image is given by the NAXES values. 
C  The EXTEND = TRUE parameter indicates that the FITS file
C  may contain extensions following the primary array.
c
C  Write the required header keywords to the file
      call ftphpr(outunit,simple,bitpix,naxis,naxes,0,1,extend,status)
      if (status .ne. 0) return
c
      if (nhdrline .gt. 0) then
        do 100 n = 1, nhdrline
          status = 0
          call ftprec(outunit,hdrline(n),status)
          if (status .ne. 0) then
            print *, 'wrfitsr WARNING - failed to place header line:'
            print *,hdrline(n)(1:lnblnk(hdrline(n)))
          end if
100     continue
      end if      
c
C  Write the array to the FITS file.
C  The last letter of the subroutine name defines the datatype of the
C  array argument; in this case the 'J' indicates that the array has an
C  integer*4 datatype. ('I' = I*2, 'E' = Real*4, 'D' = Real*8).
C  The 2D array is treated as a single 1-D array with NAXIS1 * NAXIS2
C  total number of pixels.  GROUP is seldom used parameter that should
C  almost always be set = 1.
      nelements = naxes(1)*naxes(2)
      if (naxis .eq. 3) nelements = nelements*naxes(3)
      status = 0
      call ftppre(outunit,group,fpixel,nelements,array,status)
      if (status .ne. 0) return
c
C  The FITS file must always be closed before exiting the program. 
C  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(outunit, status)
      call ftfiou(outunit, status)
c
      return
      end
c      
c=======================================================================
c
      subroutine SignOn(pgmnam)
c
c *** signon- routine which provides sign-on and sign-off messages
c             (orig by John Fowler- mod by Howard McCallon-041214-SIRTF)
c
c     inputs:  pgmnam = program name                                 [call arg]
c
c     outputs: message to stdout
c
      character*(*) pgmnam
      character vsn*11,cdate*8,ctime*8,Fmt*11,FLen*4
      integer*4 onoff,jdate(3),jtime(3),lnblnk
      real*4    dummyt,second(2),etime
c
      common /vdt/ cdate,ctime,vsn
c##
      onoff = 1
c
c         i. obtain date
c
100   cdate = '00-00-00'
      call idate(jdate)    ! Linux call
c
      jdate(3) = mod(jdate(3), 100)
      write(cdate(1:2), '(i2)') jdate(2)
      write(cdate(4:5), '(i2)') jdate(1)
      write(cdate(7:8), '(i2)') jdate(3)
c
      if(cdate(4:4) .eq. ' ') cdate(4:4) = '0'
      if(cdate(7:7) .eq. ' ') cdate(7:7) = '0'
c
c         ii. obtain time
c
      ctime = '00:00:00'
      call itime(jtime)
      write(ctime(1:2), '(i2)') jtime(1)
      write(ctime(4:5), '(i2)') jtime(2)
      write(ctime(7:8), '(i2)') jtime(3)
c
      if(ctime(4:4) .eq. ' ') ctime(4:4) = '0'
      if(ctime(7:7) .eq. ' ') ctime(7:7) = '0'
c
c         iii. set up format for pgmnam
c
      write(Flen,'(I4)') lnblnk(pgmnam)
      Fmt = '(A'//Flen//'$)'
c
c         iv. write out results
c
      write(*,Fmt) pgmnam
      if(onoff .eq. 1) then                      ! sign on
        write(*,301) vsn,cdate,ctime
      else                                       ! sign off
        dummyt = etime(second)
        write(*,302) vsn,cdate,ctime,second
      endif
  301 format(' version: ',a11,' - execution begun on ',a8,' at ',a8)
  302 format(' version: ',a11,' - execution ended on ',a8,' at ',a8
     *    /1x,f9.2,' cpu seconds used;',f8.2,' system seconds used.')
c
      return
c
      entry SignOff(pgmnam)
      OnOff = 2
      go to 100
c
      end
c      
c=======================================================================
c
      subroutine Rot8(psf,psfunc,nAngs,dH,AngHist,nOvrSamp,AngMin,
     +                psfout,psfuncout,wgt)
c     
      real*4    wgt(641,641),    sqr(641,641),    var(641,641),
     +          psf(641,641),    psfunc(641,641), AngHist(nAngs),
     +          psfout(641,641), psfuncout(641,641),
     +          OvrScale, ang, dH
      Real*8    SumPix, AngMin, d2r, ca, sa, dx, dy, rx, ry
      integer*4 nAngs, nOsc, nOvrSamp, n, i, j, ii, jj, ir, jr
c      
      Data d2r/1.745329252d-2/
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      psfout    = 0.0
      psfuncout = 0.0
      wgt       = 0.0
      if (nAngs .lt. 1) return
      nOSc      = nOvrSamp/2 + 1
      OvrScale  = 1.0d0/dfloat(nOvrSamp)
      sqr       = 0.0
c      
      do 200 n = 1, nAngs
        if (AngHist(n) .le. 0.0) go to 200
        ang = (float(n) - 0.5)* dH + AngMin
        ca = dcos(-d2r*ang)                      ! rotate through -PA
        sa = dsin(-d2r*ang)
        do 180 j = 1, 641
          do 160 i = 1, 641
            do 140 jj = 1, nOvrSamp
              dy = j - 321 + OvrScale*(jj - nOSc)
              do 120 ii = 1, nOvrSamp
                dx = i - 321 + OvrScale*(ii - nOSc)
                rx = ca*dx + sa*dy
                ry = ca*dy - sa*dx
                ir = NInt(321.0 + rx)
                jr = NInt(321.0 + ry)
                if   ((ir .gt. 0) .and. (ir .lt. 642)
     +          .and. (jr .gt. 0) .and. (jr .lt. 642)) then
                  psfout(ir,jr)    = psfout(ir,jr)
     +                             + AngHist(n)*psf(i,j)
                  sqr(ir,jr)       = sqr(ir,jr)
     +                             + AngHist(n)*psf(i,j)**2
                  psfuncout(ir,jr) = psfuncout(ir,jr)
     +                             + AngHist(n)*psfunc(i,j)
                  wgt(ir,jr) = wgt(ir,jr) + AngHist(n)
                end if
120           continue
140         continue
160       continue
180     continue
200   continue
c
      SumPix = 0.0d0
      do 500 j = 1, 641
        do 400 i = 1, 641
          if (wgt(i,j) .gt. 0.0) then
            psfout(i,j)    = psfout(i,j)/wgt(i,j)
            psfuncout(i,j) = psfuncout(i,j)/wgt(i,j)
            sqr(i,j)       = sqr(i,j)/wgt(i,j)
          else
            psfout(i,j)    = 1.0e-25
            psfuncout(i,j) = 1.0e-25
            sqr(i,j)       = 0.0
          end if
          SumPix = SumPix + psfout(i,j)
400     continue
500   continue
c                                      ! normalize PSFs, not PSFuncs
      psfout = psfout/SumPix
c
      var = sqr - psfout**2
      do 700 j = 1, 641
        do 600 i = 1, 641
          if ((wgt(i,j) .gt. 0.0) .and. (var(i,j) .gt. 0.0))
     +         psfuncout(i,j) = sqrt(psfuncout(i,j)**2 + (var(i,j)))
600      continue
700    continue        
c     
      return
c      
      end
c      
c=======================================================================
c
      subroutine DownSamp(psf,nOvrSamp,psfout)
c     
      real*4    psf(641,641), psfout(27,27), wgt(27,27), OvrScale
      Real*8    SumPix, dx, dy
      integer*4 nOsc, nOvrSamp, i, j, ii, jj, ir, jr
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      nOSc      = nOvrSamp/2 + 1
      OvrScale  = 1.0d0/dfloat(nOvrSamp)
      wgt       = 0.0
      psfout    = 0.0
c      
      do 180 j = 1, 641
        do 160 i = 1, 641
          do 140 jj = 1, nOvrSamp
            dy = j - 321 + OvrScale*(jj - nOSc)
            do 120 ii = 1, nOvrSamp
              dx = i - 321 + OvrScale*(ii - nOSc)
              ir = NInt(14.0 + dx/2.0)
              jr = NInt(14.0 + dy/2.0)
              if   ((ir .gt. 0) .and. (ir .lt. 28)
     +        .and. (jr .gt. 0) .and. (jr .lt. 28)) then
                psfout(ir,jr)    = psfout(ir,jr) + psf(i,j)
                wgt(ir,jr) = wgt(ir,jr) + 1.0
              end if
120         continue
140       continue
160     continue
180   continue
c
      SumPix = 0.0d0
      do 500 j = 1, 27
        do 400 i = 1, 27
          if (wgt(i,j) .gt. 0.0) then
            psfout(i,j)    = psfout(i,j)/wgt(i,j)
          else
            psfout(i,j)    = 1.0e-25
          end if
          SumPix = SumPix + psfout(i,j)
400     continue
500   continue
c
      psfout = psfout/SumPix
c
      return
c
      end
c      
c=======================================================================
c
c     Check for embedded zeroes in Angs array, indicating a PA zero
c     crossing; allow if found, because those histogram elements will
c     be skipped in the rotation subroutine; otherwise collapse array
c     to less than 100 elements.
c    
      subroutine ChkPAX(Angs,nAngs,dH)
c
      integer*4 nAngs, n, nZero, nCompress, j, nOff, nMax
      real*4    Angs(nAngs), dH, Ang
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      if (nAngs .lt. 100) return
      nZero = 0
      nMax  = nAngs
c     print *,'ChkPAX: nAngs =', nAngs   ! dbg
c     
      do 10 n = 1, nAngs
        if (Angs(n) .le. 0.0) nZero = nZero + 1
10    continue
c
c     print *,'ChkPAX: nZero =', nZero   ! dbg
      if (nZero .gt. nAngs/2) return
c 
      nCompress = nAngs/50
      if (nCompress .lt. 2) return
      print *,'WARNING: compressing angle histogram by a factor of',
     +         nCompress
      print *,'         original no. of cells:  ', nAngs     
c      
      nOff = 0    
      do 50 n = 1, nAngs
        Ang = 0.0
        do 30 j = 1, nCompress
          nOff = nOff + 1
          if (nOff .gt. nAngs) go to 40
          Ang = Ang + Angs(nOff)
30      continue
40      Angs(n) = Ang
        nMax    = n
        if (nOff .gt. nAngs) go to 100
50    continue
c 
100   nAngs = nMax
      print *,'         compressed no. of cells:', nAngs     
      dH = float(nCompress)*dH
c    
      return
      end
