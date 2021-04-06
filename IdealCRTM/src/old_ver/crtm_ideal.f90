PROGRAM crtm_ideal
!
! May 13, 2020: No RH dependency,
!
USE netcdf
USE CRTM_Module
USE CRTM_Model_Profiles
USE CRTM_SpcCoeff        , ONLY: SC
USE SpcCoeff_Define      , ONLY: SpcCoeff_type, SpcCoeff_Inspect
IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME   = 'crtm_ideal'
  CHARACTER(*), PARAMETER :: PROGRAM_VERSION_ID = &
    '$Id: check_crtm.fpp 60152 2015-08-13 19:19:13Z paul.vandelst@noaa.gov $'

  CHARACTER(*), PARAMETER :: ENDIAN_TYPE='big_endian'
  CHARACTER(*), PARAMETER :: COEFFICIENT_PATH='./coefficients/'

  ! Directory location of results for comparison [NOT USED YET]
  CHARACTER(*), PARAMETER :: RESULTS_PATH = './results/'

  ! Profile dimensions
  INTEGER, PARAMETER :: N_PROFILES  = 1
  INTEGER, PARAMETER :: N_LAYERS    = 100
  INTEGER, PARAMETER :: N_ABSORBERS = 3
  INTEGER, PARAMETER :: Absorber_IDs(N_ABSORBERS) = &
                        (/H2O_ID, &
                           O3_ID, & 
                          CO2_ID/)
  INTEGER, PARAMETER :: Absorber_units(N_ABSORBERS) = &
                        (/MASS_MIXING_RATIO_UNITS, &
                        VOLUME_MIXING_RATIO_UNITS, &
                        VOLUME_MIXING_RATIO_UNITS/)
  INTEGER, PARAMETER :: N_CLOUDS    = 0
  INTEGER :: N_AEROSOLS

  ! Sensor information
  INTEGER     , PARAMETER :: N_SENSORS = 2
  CHARACTER(*), PARAMETER :: SENSOR_ID(N_SENSORS) = (/'iasi_metop-a   ', &
                                                      'v.modis_aqua   '/)

  ! Some pretend geometry angles. The scan angle is based
  ! on the default Re (earth radius) and h (satellite height)
  REAL(fp), PARAMETER :: ZENITH_ANGLE = 0.0_fp !30.0_fp
  REAL(fp), PARAMETER :: SCAN_ANGLE   = 0.0_fp !26.37293341421_fp
  REAL(fp), PARAMETER :: secant_term  = 1./COS(ZENITH_ANGLE)

  ! ---------
  ! Variables
  ! ---------
  CHARACTER(256) :: message, version
  INTEGER :: err_stat, alloc_stat
  INTEGER :: n_channels
  INTEGER :: i, m, n, ii, k
  character(len=20),dimension(1) :: sensorlist
  CHARACTER(len=8) :: absorbername
  REAL(fp), ALLOCATABLE :: work1d(:),work2d(:,:),replacetmp(:),ptau5(:,:)
  REAL(fp) :: maxtmp, totalaerconc, tmptotalod
  INTEGER :: slv,elv,newslv,newelv,maxlvidx
  ! For 20-80 autogen, 
  REAL(fp), PARAMETER :: user_pi=4*atan(1.d0)
  INTEGER :: nblvs, kk, nlvs_process, st_proc_lv, nb
  REAL(fp) :: abovemode, belowmode, alvsmean, blvsmean, modeinc
  REAL(fp) :: degint_bl
  REAL(fp), ALLOCATABLE :: weighting(:)
  REAL(fp) :: deg,gentotal
  REAL(fp) :: decaying_ratio,residual
  
  ! Namelist
  LOGICAL :: lextprofile
  INTEGER :: naers
  CHARACTER(len=256) :: aerprof
  REAL(fp) :: scalef
  INTEGER :: lumptolvs,puttolv,shiftthick,shiftlvs
  CHARACTER(len=32) :: varname
  CHARACTER(len=16) :: binsname(5)
  REAL(fp) :: totalconc
  INTEGER :: thicknesslvs, atlv
  INTEGER :: genmethod
  INTEGER :: modelv
  REAL(fp) :: bmoderatio,sfc_mode_ratio
  REAL(fp) :: binspar(5)
  REAL(fp) :: landcover,lai
  INTEGER :: landtype

  namelist/main/lextprofile,naers
  namelist/extprof/aerprof,scalef,lumptolvs,puttolv,shiftthick,shiftlvs
  namelist/autogen/varname,binsname,binspar,genmethod,totalconc,thicknesslvs,atlv,modelv,bmoderatio,sfc_mode_ratio
  namelist/sfcsetup/landcover,landtype,lai
  ! Aerosol Profile
  CHARACTER(len=10) :: aername
  INTEGER :: aerlvs
  REAL(fp), DIMENSION(:), ALLOCATABLE :: aerconc, aerpres
    

  ! ============================================================================
  ! STEP 3. **** DEFINE THE CRTM INTERFACE STRUCTURES ****
  !
  ! 3a. Define the "non-demoninational" arguments
  ! ---------------------------------------------
  TYPE(CRTM_ChannelInfo_type),DIMENSION(1)             :: chinfo
  TYPE(CRTM_Geometry_type),DIMENSION(1)                :: geo

  ! 3b. Define the FORWARD variables
  ! --------------------------------
  TYPE(CRTM_Atmosphere_type),DIMENSION(1)              :: atm
  TYPE(CRTM_Surface_type),DIMENSION(1)                 :: sfc
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts(:,:)


  ! 3c. Define the K-MATRIX variables
  ! ---------------------------------
  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: atm_K(:,:)
  TYPE(CRTM_Surface_type)   , ALLOCATABLE :: sfc_K(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts_K(:,:)

  !
  !TYPE(SpcCoeff_type), ALLOCATABLE :: SC(:)

  ! Define NetCDF output variables
  INTEGER :: ncid, levid, lev_varid, varid, &
             chnlid, chnl_varid, schnlid, waveid
  INTEGER :: presid, tempid, preskid, tempkid
  INTEGER :: dimid2d(2)
  INTEGER, ALLOCATABLE :: levidx(:), chlidx(:)
  CHARACTER(len=4) :: crtmfunc

  SENSORLOOP: DO m=1,N_SENSORS
  sensorlist(1)=TRIM(SENSOR_ID(m))
  WRITE( *,'(/5x,"Initializing the CRTM...",a)' ) TRIM(sensorlist(1))
  err_stat = CRTM_Init( sensorlist, &
                        chinfo, &
                        File_Path=COEFFICIENT_PATH, &
                        Quiet=.TRUE.)
  IF ( err_stat /= SUCCESS ) THEN
    message = 'Error initializing CRTM for '//sensorlist(1)
    CALL Display_Message( PROGRAM_NAME, message, FAILURE )
    STOP
  END IF

  ! Read in namelist
  OPEN(10,file='idealcrtm.nml')
  READ(10,main)
  IF (lextprofile) THEN
     READ(10,extprof)
     WRITE(6,*) 'Aerosol profile is ',trim(aerprof)
  ELSE
     READ(10,autogen)
  END IF
  READ(10,sfcsetup)
  CLOSE(10)
  N_AEROSOLS=naers

  ! ==========================================================================
  ! STEP 5. **** ALLOCATE STRUCTURE ARRAYS ****
  !
  ! 5a. Determine the number of channels
  !     for the current sensor
  ! ------------------------------------
  n_channels = CRTM_ChannelInfo_n_Channels(chinfo(1))
  WRITE( *,'(/5x,"Processing a total of ",i0," channels...")' ) n_channels
  WRITE( *,'(7x,i0," from ",a)' ) &
  CRTM_ChannelInfo_n_Channels(chinfo(1)), TRIM(sensorlist(1))

  ! 5b. Allocate the ARRAYS
  ! -----------------------
  ALLOCATE( rts( n_channels, 1 ), &
            atm_K( n_channels, 1 ), &
            sfc_K( n_channels, 1 ), &
            rts_K( n_channels, 1 ), &
            chlidx( n_channels ), &
            levidx( N_LAYERS ), &
            ptau5 ( N_LAYERS,n_channels ), &
            work2d( N_LAYERS,n_channels ), &
            STAT = alloc_stat )
  IF ( alloc_stat /= 0 ) THEN
    message = 'Error allocating structure arrays'
    CALL Display_Message( PROGRAM_NAME, message, FAILURE )
    STOP
  END IF

  ! 5c. Allocate the STRUCTURE INTERNALS
  !     NOTE: Only the Atmosphere structures
  !           are allocated in this example
  ! ----------------------------------------
  CALL CRTM_Atmosphere_Create( atm, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  IF ( ANY(.NOT. CRTM_Atmosphere_Associated(atm)) ) THEN
    message = 'Error allocating CRTM Forward Atmosphere structure'
    CALL Display_Message( PROGRAM_NAME, message, FAILURE )
    STOP
  END IF

  CALL CRTM_RTSolution_Create( rts, N_LAYERS )
  CALL CRTM_RTSolution_Create( rts_K, N_LAYERS )

  DO ii=1,n_Channels
     atm_K(ii,1)=atm(1)
     sfc_K(ii,1)=sfc(1)
  ENDDO
  ! ==========================================================================
  ! STEP 6. **** ASSIGN INPUT DATA ****
  !
  ! 6a. Atmosphere and Surface input
  !     NOTE: that this is the hard part (in my opinion :o). The mechanism by
  !     by which the atmosphere and surface data are loaded in to their
  !     respective structures below was done purely to keep the step-by-step
  !     instructions in this program relatively "clean".
  ! ------------------------------------------------------------------------
  CALL CRTM_Get_Model_Profile(Absorber_IDs, atm(1)%Level_Pressure(:),  &
                              atm(1)%Temperature(:), atm(1)%Absorber(:,:) )
  atm(1)%Absorber_ID(:) = Absorber_IDs
  atm(1)%Absorber_Units(:) = Absorber_units
!  atm(1)%Absorber(:,:) = 0.
!  atm(1)%Absorber(:,2) = atm(1)%Absorber(:,2)*0.6
  
  ! Interpolation for model layers pressure
  DO i=1,N_LAYERS
     atm(1)%Pressure(i)=exp(0.5*(log(atm(1)%Level_Pressure(i-1))+log(atm(1)%Level_Pressure(i))))
  ENDDO

  ! Read in the aerosol profile and interpolate to model layers
  IF (N_AEROSOLS .NE. 0 ) THEN
     IF ( lextprofile ) THEN
     WRITE(*,*) 
     OPEN(11, file=aerprof, form='formatted')
     READ(11, *) aername, aerlvs
     ALLOCATE(aerconc(aerlvs), aerpres(aerlvs))
     DO i=1,aerlvs
        READ(11, *) aerpres(i), aerconc(i)
     ENDDO
     CLOSE(11)
     SELECT CASE (trim(aername))
     CASE('DU001')
         atm(1)%Aerosol(1)%Type = DUST_AEROSOL
         atm(1)%Aerosol(1)%Effective_Radius = 0.55
     CASE('DU002')
         atm(1)%Aerosol(1)%Type = DUST_AEROSOL
         atm(1)%Aerosol(1)%Effective_Radius = 1.4
     CASE('DU003')
         atm(1)%Aerosol(1)%Type = DUST_AEROSOL
         atm(1)%Aerosol(1)%Effective_Radius = 2.4
     CASE('DU004')
         atm(1)%Aerosol(1)%Type = DUST_AEROSOL
         atm(1)%Aerosol(1)%Effective_Radius = 4.5
     CASE('DU005')
         atm(1)%Aerosol(1)%Type = DUST_AEROSOL
         atm(1)%Aerosol(1)%Effective_Radius = 8.0
     END SELECT
 
  ! Interpolation
     DO i=1,N_LAYERS
        loop_aerlvs: DO ii=1,aerlvs
           
           IF (atm(1)%Pressure(i) .GT. aerpres(ii) .AND. &
               atm(1)%Pressure(i) .LE. aerpres(ii+1) ) THEN
               atm(1)%Aerosol(1)%Concentration(i)=aerconc(ii+1)*((log(atm(1)%Pressure(i)) - log(aerpres(ii))) / &
                                                                 (log(aerpres(ii+1))      - log(aerpres(ii))))+ &
                                                  aerconc(ii)*((log(aerpres(ii+1)) - log(atm(1)%Pressure(i))) / &
                                                               (log(aerpres(ii+1)) - log(aerpres(ii))))
               EXIT loop_aerlvs  
           ELSE IF (atm(1)%Pressure(i) .GT. aerpres(aerlvs) ) THEN
               atm(1)%Aerosol(1)%Concentration(i)=aerconc(aerlvs)
               EXIT loop_aerlvs  
           ELSE IF (atm(1)%Pressure(i) .LT. aerpres(1) ) THEN
               atm(1)%Aerosol(1)%Concentration(i)=aerconc(1)
               EXIT loop_aerlvs  
           ENDIF
        ENDDO loop_aerlvs
     ENDDO
     atm(1)%Aerosol(1)%Concentration(:)=atm(1)%Aerosol(1)%Concentration(:)*scalef
     IF (lumptolvs .gt.0 ) THEN
        totalaerconc=sum(atm(1)%Aerosol(1)%Concentration(:))
        atm(1)%Aerosol(1)%Concentration(:)=0.
        atm(1)%Aerosol(1)%Concentration(puttolv:puttolv+lumptolvs-1)=totalaerconc/lumptolvs
     END IF
     IF (shiftthick .gt. 0) THEN
        maxtmp=-9999
        do i=1,N_LAYERS
           if ( atm(1)%Aerosol(1)%Concentration(i) > maxtmp ) then
              maxlvidx=i
              maxtmp=atm(1)%Aerosol(1)%Concentration(i)
           end if
        end do
        IF (mod(shiftthick,2) .ne. 0) THEN
           slv=maxlvidx-int(shiftthick/2.)
        ELSE
           slv=maxlvidx-int(shiftthick/2.)+1
        END IF
        elv=maxlvidx+int(shiftthick/2.)
        allocate(replacetmp(shiftthick))
        replacetmp(:)=atm(1)%Aerosol(1)%Concentration(slv:elv)
        newslv=slv+shiftlvs
        newelv=elv+shiftlvs
        atm(1)%Aerosol(1)%Concentration(slv:elv)=atm(1)%Aerosol(1)%Concentration(newslv:newelv)
        atm(1)%Aerosol(1)%Concentration(newslv:newelv)=replacetmp(:)
        deallocate(replacetmp)
     END IF
     DEALLOCATE(aerconc, aerpres)
     ! End of if read in pregenerated profile
     ELSE
     ! Start autogen
     SELECT CASE (trim(varname))
     CASE('gocart_dust')
         if (naers .ne. 5) then
            write(6,*) '!!!Error!!! naers should be 5 for gocart_dust' 
            stop
         end if
         do i=1,naers
            atm(1)%Aerosol(i)%Type = DUST_AEROSOL
         end do
         atm(1)%Aerosol(1)%Effective_Radius = 0.55
         atm(1)%Aerosol(2)%Effective_Radius = 1.4
         atm(1)%Aerosol(3)%Effective_Radius = 2.4
         atm(1)%Aerosol(4)%Effective_Radius = 4.5
         atm(1)%Aerosol(5)%Effective_Radius = 8.0
     END SELECT
 
        IF (genmethod .eq. 1) THEN
           ! boxcar 
           WRITE(*,*) "Generate boxcar artificial profile"
           do nb=1,N_AEROSOLS
              atm(1)%Aerosol(nb)%Concentration(atlv:atlv+thicknesslvs-1)=totalconc/thicknesslvs*binspar(nb)
           end do
        ELSE IF (genmethod .eq. 2) THEN
           ! 20-80: use totalconc, modelv, extendlvs, sfc_mode_ratio 
           WRITE(*,*) "Generate artificial profile"
           do nb=1,naers
              atm(1)%Aerosol(nb)%Concentration(:)=0.
           end do
           abovemode=totalconc*(1.-bmoderatio)
           belowmode=totalconc*bmoderatio
           nblvs=N_LAYERS-modelv+1  ! number of levels below mode (including mode) 
           blvsmean=belowmode/nblvs

           !write(6,*) bmoderatio,abovemode,belowmode
           
           atm(1)%Aerosol(naers)%Concentration(modelv:N_LAYERS)=blvsmean

           !do k=1,N_LAYERS
           !   write(6,*) atm(1)%Aerosol(1)%Concentration(k)
           !end do

           degint_bl=180.0_fp/(nblvs-1)
           
           allocate(weighting(nblvs))
           weighting=0.
           ! below the mode (included)
           do k=1,nblvs
              kk=N_LAYERS-k+1
              if (kk .eq. N_LAYERS) then
                 deg=180.
                 weighting(k)=cos(deg*pi/180.0_fp)
              else if (kk .lt. N_LAYERS .and. kk .ge. modelv) then
                 deg=deg-degint_bl
                 weighting(k)=cos(deg*pi/180.0_fp)
              end if
              !write(6,*) kk,weighting(k)
           end do
              
           modeinc=blvsmean*((1.-sfc_mode_ratio)/(1.+sfc_mode_ratio))
           !write(6,*) 'increment below mode',modeinc

           do k=1,nblvs
              kk=N_LAYERS-k+1
              atm(1)%Aerosol(naers)%Concentration(kk)=atm(1)%Aerosol(naers)%Concentration(kk)+modeinc*weighting(k)
              !write(6,*) kk,atm(1)%Aerosol(1)%Concentration(kk)
           end do 

           !write(6,*) 'Finish the profile below mode'

           decaying_ratio=atm(1)%Aerosol(naers)%Concentration(modelv)*0.8/abovemode

           residual=abovemode
           do k=1,modelv-1
              kk=(modelv-1)-k+1
              if (kk .eq. modelv-1) then ! 0.8 is the ratio of modelv-1 and modelv
                 atm(1)%Aerosol(naers)%Concentration(kk)=atm(1)%Aerosol(naers)%Concentration(modelv)*0.8
              else if (kk .lt. modelv-1) then
                 atm(1)%Aerosol(naers)%Concentration(kk)=residual*decaying_ratio
              else
                 write(6,*) 'Error: kk is larger than modelv'
                 stop
              end if
              residual=residual-atm(1)%Aerosol(naers)%Concentration(kk)
              !write(6,*) kk,atm(1)%Aerosol(1)%Concentration(kk)
           end do
           
           write(*,*) 'total concentration: ',totalconc,sum(atm(1)%Aerosol(naers)%Concentration(:))

           ! spread to different bins
           do nb=1,naers
              atm(1)%Aerosol(nb)%Concentration(:)=atm(1)%Aerosol(naers)%Concentration(:)*binspar(nb)
           end do

           !do k=1,nblvs+extendlvs+1
           !   write(*,*) degs(k),cosvalue(k)
           !end do
           deallocate(weighting)
        ELSE IF (genmethod .eq. 99) THEN
           do nb=1,naers
              atm(1)%Aerosol(nb)%Concentration(:)=0_fp
           end do
        END IF ! genmethod
     END IF
  ENDIF

  CALL Load_Sfc_Data()
 

  ! 6b. Geometry input
  ! ------------------
  ! All profiles are given the same value
  !  The Sensor_Scan_Angle is optional.
  CALL CRTM_Geometry_SetValue( geo, &
                               Sensor_Zenith_Angle = ZENITH_ANGLE, &
                               Sensor_Scan_Angle   = SCAN_ANGLE )
                              !Source_Zenith_Angle = 0.  
                              !Source_Azimuth_Angle= 180.  
  ! ==========================================================================
  ! STEP 7. **** INITIALIZE THE K-MATRIX ARGUMENTS ****
  !
  ! 7a. Zero the K-matrix OUTPUT structures
  ! ---------------------------------------
  CALL CRTM_Atmosphere_Zero( atm_K )
  CALL CRTM_Surface_Zero( sfc_K )

  ! 7b. Inintialize the K-matrix INPUT so
  !     that the results are dTb/dx
  ! -------------------------------------
  DO i=1,n_channels
     rts_K(i,1)%Radiance               = ZERO
     rts_K(i,1)%Brightness_Temperature = ONE
  END DO
  ! ==========================================================================
  ! ==========================================================================
  ! STEP 8. **** CALL THE CRTM FUNCTIONS FOR THE CURRENT SENSOR ****
  !
  !write(6,*) atm(1)%Absorber_ID(:)
  !write(6,*) atm(1)%Absorber_Units(:)
  !write(6,*) atm(1)%Absorber(:,2)

  ! 8a. The K-matrix model
  ! ----------------------
  IF ( sensorlist(1) .ne. 'v.modis_aqua' ) THEN
     err_stat = CRTM_K_Matrix( atm        , &  ! FORWARD  Input
                               sfc        , &  ! FORWARD  Input
                               rts_K      , &  ! K-MATRIX Input
                               geo        , &  ! Input
                               chinfo     , &  ! Input
                               atm_K      , &  ! K-MATRIX Output
                               sfc_K      , &  ! K-MATRIX Output
                               rts          )  ! FORWARD  Output
     crtmfunc='kmod'
  ELSE
     err_stat = CRTM_AOD_K( atm        , &  ! FORWARD  Input
                            rts_K      , &  ! K-MATRIX Input
                            chinfo     , &  ! Input
                            rts        , &  ! FORWARD  Output
                            atm_K        )  ! K-MATRIX Output
     crtmfunc='aodk'
  END IF
     

  IF ( err_stat /= SUCCESS ) THEN
    message = 'Error calling CRTM K-Matrix Model for '//TRIM(SENSOR_ID(n))
    CALL Display_Message( PROGRAM_NAME, message, FAILURE )
    STOP
  END IF
 
  do i=1,n_channels
     tmptotalod=0.
     do k=1,N_LAYERS
        tmptotalod=tmptotalod+rts(i,1)%Layer_Optical_Depth(k)
        ptau5(k,i)=exp(-min(limit_exp,tmptotalod*secant_term))
     end do
  end do
 
  do i=1,N_LAYERS
     levidx(i)=i
  enddo

  do i=1,n_channels
     chlidx(i)=i
  enddo
  !write(6,*) 'nlayers',atm(1)%n_Layers
  !write(6,*) 'level indice',levidx
  !write(6,*) 'channel indice',chlidx

  ! Write out NetCDF for FORWARD input and output
  call check_nc( nf90_create(trim(crtmfunc)//'_inspection.nc',NF90_CLOBBER, ncid) )
  call check_nc( nf90_def_dim(ncid,'layers',atm(1)%n_Layers,levid) ) !atm%n_Layers
  call check_nc( nf90_def_dim(ncid,'wavenumber',n_channels,chnlid) ) !atm%n_Layers
  !call check_nc( nf90_def_dim(ncid,'channels',n_channels,chnlid) ) !atm%n_Layers
  dimid2d=(/levid,chnlid/)

  call check_nc( nf90_def_var(ncid,'layers',NF90_INT,levid,lev_varid) )
  call check_nc( nf90_put_att(ncid,lev_varid,"long_name","model layers") )

  call check_nc( nf90_def_var(ncid,'channels',NF90_INT,chnlid,chnl_varid) )
  call check_nc( nf90_put_att(ncid,chnl_varid,"long_name","sensor channels") )

  call check_nc( nf90_def_var(ncid,'wavenumber',NF90_REAL,chnlid,waveid) )
  call check_nc( nf90_put_att(ncid,waveid,"long_name","Wavenumber") )
  call check_nc( nf90_put_att(ncid,waveid,"units","cm^-1") )
 
  call check_nc( nf90_def_var(ncid,'sensor_channels_idx',NF90_INT,chnlid,schnlid) )
  call check_nc( nf90_put_att(ncid,schnlid,"long_name","Sensor Channels Index") )

  call check_nc( nf90_def_var(ncid,'pres',NF90_REAL,levid,presid) )
  call check_nc( nf90_put_att(ncid,presid,"long_name","Layer pressure") )
  call check_nc( nf90_put_att(ncid,presid,"units","hPa") )

  call check_nc( nf90_def_var(ncid,'pres_K',NF90_REAL,dimid2d,preskid) )
  call check_nc( nf90_put_att(ncid,preskid,"long_name","Jacobian for Layer pressure") )
  call check_nc( nf90_put_att(ncid,preskid,"units","K hPa^-1") )

  call check_nc( nf90_def_var(ncid,'temp',NF90_REAL,levid,tempid) )
  call check_nc( nf90_put_att(ncid,tempid,"long_name","Layer temperature") )
  call check_nc( nf90_put_att(ncid,tempid,"units","K") )

  call check_nc( nf90_def_var(ncid,'temp_K',NF90_REAL,dimid2d,tempkid) )
  call check_nc( nf90_put_att(ncid,tempkid,"long_name","Jacobian for Layer temperature") )
  call check_nc( nf90_put_att(ncid,tempkid,"units","K K^-1") )

  !--- Global Attributes
  call check_nc( nf90_put_att(ncid,NF90_GLOBAL,'Sensor_ID',trim(SC(1)%Sensor_Id)) )
  call check_nc( nf90_put_att(ncid,NF90_GLOBAL,'N_AEROSOLS',N_AEROSOLS) )
  call check_nc( nf90_put_att(ncid,NF90_GLOBAL,'RT_Algorithm_Name',rts(1,1)%RT_Algorithm_Name) )
  IF (lextprofile) THEN
     call check_nc( nf90_put_att(ncid,NF90_GLOBAL,'AerosolFile',aerprof) )
     call check_nc( nf90_put_att(ncid,NF90_GLOBAL,'Scale_Factor',scalef) )
     call check_nc( nf90_put_att(ncid,NF90_GLOBAL,'Lump_to_levels',lumptolvs) )
     if (lumptolvs .gt. 0) then
        call check_nc( nf90_put_att(ncid,NF90_GlOBAL,'Lump_to',puttolv) )
     end if
     if (shiftthick .ne. 0) then
        call check_nc( nf90_put_att(ncid,NF90_GlOBAL,'Shift_Thickness',shiftthick) )
        call check_nc( nf90_put_att(ncid,NF90_GlOBAL,'Shift_Levels',shiftlvs) )
     end if
  ELSE
     call check_nc( nf90_put_att(ncid,NF90_GLOBAL,'Total_Concentration',totalconc) )
     call check_nc( nf90_put_att(ncid,NF90_GLOBAL,'Gen_Method',genmethod) )
     if (genmethod .eq. 1) then
        call check_nc( nf90_put_att(ncid,NF90_GLOBAL,'Thickness_Levels',thicknesslvs) )
        call check_nc( nf90_put_att(ncid,NF90_GLOBAL,'At_Level',atlv) )
     else if (genmethod .eq. 2) then
        call check_nc( nf90_put_att(ncid,NF90_GLOBAL,'Mode_at_Lv',modelv) )
        call check_nc( nf90_put_att(ncid,NF90_GLOBAL,'DecayingRatio_above_Mode',decaying_ratio) )
        call check_nc( nf90_put_att(ncid,NF90_GLOBAL,'Fraction_below_Mode',bmoderatio) )
        call check_nc( nf90_put_att(ncid,NF90_GLOBAL,'Sfc_Mode_Ratio',sfc_mode_ratio) )
     end if
  END IF
  IF (crtmfunc .eq. 'aodk') THEN
     call check_nc( nf90_put_att(ncid,NF90_GLOBAL,'AOD550',sum(rts(4,1)%Layer_Optical_Depth(:))) )
  END IF
  call check_nc( nf90_put_att(ncid,NF90_GLOBAL,'Land_Type',sfc(1)%Land_Type) )
  call check_nc( nf90_put_att(ncid,NF90_GLOBAL,'Land_Coverage',sfc(1)%Land_Coverage) )
  call check_nc( nf90_put_att(ncid,NF90_GLOBAL,'Land_Temperature',sfc(1)%Land_Temperature) )
  call check_nc( nf90_put_att(ncid,NF90_GLOBAL,'Water_Coverage',sfc(1)%Water_Coverage) )
  call check_nc( nf90_put_att(ncid,NF90_GLOBAL,'Water_Temperature',sfc(1)%Water_Temperature) )
  call check_nc( nf90_put_att(ncid,NF90_GLOBAL,'Snow_Coverage',sfc(1)%Snow_Coverage) )
  call check_nc( nf90_put_att(ncid,NF90_GLOBAL,'Snow_Temperature',sfc(1)%Snow_Temperature) )
  call check_nc( nf90_put_att(ncid,NF90_GLOBAL,'Ice_Coverage',sfc(1)%Ice_Coverage) )
  call check_nc( nf90_put_att(ncid,NF90_GLOBAL,'Ice_Temperature',sfc(1)%Ice_Temperature) )
  call check_nc( nf90_put_att(ncid,NF90_GLOBAL,'Leaf_Area_Index',sfc(1)%Lai) )

  call check_nc( nf90_enddef(ncid) )

  call check_nc( nf90_put_var(ncid,lev_varid,levidx) )
  call check_nc( nf90_put_var(ncid,chnl_varid,chlidx) )
  call check_nc( nf90_put_var(ncid,waveid,SC(1)%Wavenumber))
  call check_nc( nf90_put_var(ncid,schnlid,SC(1)%Sensor_Channel))

  call check_nc( nf90_put_var(ncid,presid,atm(1)%Pressure(:)) )
  call check_nc( nf90_put_var(ncid,tempid,atm(1)%Temperature(:)) ) 

  work2d(:,:)=0.
  do ii = 1, n_channels
     work2d(:,ii)=atm_k(ii,1)%Pressure(:)
     call check_nc( nf90_put_var(ncid,preskid,work2d) )
  end do
  work2d(:,:)=0.
  do ii = 1, n_channels
     work2d(:,ii)=atm_k(ii,1)%Temperature(:)
     call check_nc( nf90_put_var(ncid,tempkid,work2d) )
  end do
! Absorbers 
  if ( atm(1)%n_absorbers .ne. 0 ) then
     do i = 1, atm(1)%n_absorbers
        absorbername=trim(ABSORBER_ID_NAME(atm(1)%Absorber_Id(i)))
        work2d(:,:)=0.
        call check_nc( nf90_redef(ncid) )
        call check_nc( nf90_def_var(ncid,trim(absorbername),NF90_REAL,levid,varid) )
        call check_nc( nf90_put_att(ncid,varid,"units",ABSORBER_UNITS_NAME(atm(1)%Absorber_Units(i))) )
        call check_nc( nf90_enddef(ncid) )
        call check_nc( nf90_put_var(ncid,varid,atm(1)%Absorber(:,i)) )
        
        call check_nc( nf90_redef(ncid) )
        call check_nc( nf90_def_var(ncid,trim(absorbername)//'_K',NF90_REAL,dimid2d,varid) )
        call check_nc( nf90_put_att(ncid,varid,"units",'K ('//trim(ABSORBER_UNITS_NAME(atm(1)%Absorber_Units(i)))//')^-1') )
        call check_nc( nf90_put_att(ncid,varid,"long_name","Jacobian for "//trim(absorbername)) )
        call check_nc( nf90_enddef(ncid) )
        do ii = 1, n_channels
           work2d(:,ii)=atm_k(ii,1)%Absorber(:,i)
        end do 
        call check_nc( nf90_put_var(ncid,varid,work2d) ) 
     end do
  end if 
! Aerosol relative output
  if ( atm(1)%n_aerosols .ne. 0 ) then
     do i = 1, atm(1)%n_aerosols
        aername=binsname(i)
        call check_nc( nf90_redef(ncid) )
           call check_nc( nf90_def_var(ncid,trim(aername)//'_Re',NF90_REAL,levid,varid) )
           call check_nc( nf90_put_att(ncid,varid,"units","µm") )
        call check_nc( nf90_enddef(ncid) )

        call check_nc( nf90_put_var(ncid,varid,atm(1)%aerosol(i)%Effective_Radius(:)) )

        call check_nc( nf90_redef(ncid) )
           call check_nc( nf90_def_var(ncid,trim(aername)//'_Conc',NF90_REAL,levid,varid) )
           call check_nc( nf90_put_att(ncid,varid,"units","kg m^-2") )
        call check_nc( nf90_enddef(ncid) )

        call check_nc( nf90_put_var(ncid,varid,atm(1)%aerosol(i)%Concentration(:)) )

        call check_nc( nf90_redef(ncid) )
        call check_nc( nf90_def_var(ncid,trim(aername)//'_Re_K',NF90_REAL,dimid2d,varid) )
        call check_nc( nf90_put_att(ncid,varid,"units",'K µm^-1') )
        call check_nc( nf90_put_att(ncid,varid,"long_name","Jacobian for "//trim(aername)//" Re") )
        call check_nc( nf90_enddef(ncid) )
        work2d(:,:)=0.
        do ii = 1, n_channels
           work2d(:,ii)=atm_k(ii,1)%Aerosol(i)%Effective_Radius(:)
        end do 
        call check_nc( nf90_put_var(ncid,varid,work2d) ) 

        call check_nc( nf90_redef(ncid) )
        call check_nc( nf90_def_var(ncid,trim(aername)//'_Conc_K',NF90_REAL,dimid2d,varid) )
        call check_nc( nf90_put_att(ncid,varid,"units",'K (kg/m^-2)^-1') )
        call check_nc( nf90_put_att(ncid,varid,"long_name","Jacobian for "//trim(aername)//" Conc") )
        call check_nc( nf90_enddef(ncid) )
        work2d(:,:)=0.
        do ii = 1, n_channels
           work2d(:,ii)=atm_k(ii,1)%Aerosol(i)%Concentration(:)
        end do 
        call check_nc( nf90_put_var(ncid,varid,work2d) ) 
     end do
  end if 
! RT Solution output
  allocate(work1d(n_channels))
  work1d(:)=rts(:,1)%SOD
  call nc_write_real_var1d(ncid,"SOD",chnlid,work1d, &
                    varlongname="Scattering Optical Depth")
  work1d(:)=rts(:,1)%Surface_Emissivity
  call nc_write_real_var1d(ncid,"Surface_Emissivity",chnlid,work1d, &
                    varlongname="Surface Emissivity")
  work1d(:)=rts(:,1)%Surface_Reflectivity
  call nc_write_real_var1d(ncid,"Surface_Reflectivity",chnlid,work1d, &
                    varlongname="Surface Reflectivity")
  work1d(:)=rts(:,1)%Surface_Planck_Radiance
  call nc_write_real_var1d(ncid,"Surface_Planck_Radiance",chnlid,work1d, &
                    varlongname="Surface Planck Radiance")
  work1d(:)=rts(:,1)%Radiance
  call nc_write_real_var1d(ncid,"Radiance",chnlid,work1d, &
                    varlongname="Radiance")
  work1d(:)=rts(:,1)%Brightness_Temperature
  call nc_write_real_var1d(ncid,"BT",chnlid,work1d, &
                    varlongname="Brightness Temperature", varunits="K")
  do ii = 1, n_channels
     work2d(:,ii)=rts(ii,1)%Layer_Optical_Depth(:)
  end do
  call nc_write_real_var2d(ncid,"OD",dimid2d,work2d, &
                    varlongname="Layer Optical Depth")

  work1d(:)=rts_k(:,1)%SOD
  call nc_write_real_var1d(ncid,"SOD_K",chnlid,work1d, &
                    varlongname="Jacobian for Scattering Optical Depth")
  work1d(:)=rts_k(:,1)%Surface_Emissivity
  call nc_write_real_var1d(ncid,"Surface_Emissivity_K",chnlid,work1d, &
                    varlongname="Jacobian for Surface Emissivity")
  work1d(:)=rts_k(:,1)%Surface_Reflectivity
  call nc_write_real_var1d(ncid,"Surface_Reflectivity_K",chnlid,work1d, &
                    varlongname="Jacobian for Surface Reflectivity")
  work1d(:)=rts_k(:,1)%Radiance
  call nc_write_real_var1d(ncid,"Radiance_K",chnlid,work1d, &
                    varlongname="Jacobian for Radiance")
  work1d(:)=rts_k(:,1)%Brightness_Temperature
  call nc_write_real_var1d(ncid,"BT_K",chnlid,work1d, &
                    varlongname="Jacobian for Brightness Temperature", varunits="K")
  do ii = 1, n_channels
     work2d(:,ii)=rts_k(ii,1)%Layer_Optical_Depth(:)
  end do
  call nc_write_real_var2d(ncid,"OD_K",dimid2d,work2d, &
                    varlongname="Jacobian for Layer Optical Depth", varunits="K")
  ! Surface sensitivity
  if ( sfc(1)%Water_Coverage .ne. 0 ) then
     work1d(:)=sfc_K(:,1)%Water_Temperature
     call nc_write_real_var1d(ncid,"WaterT_K",chnlid,work1d, &
                    varlongname="Jacobian for Water Temperature", varunits="K K^-1")
  end if
  if ( sfc(1)%Land_Coverage .ne. 0 ) then
     work1d(:)=sfc_K(:,1)%Land_Temperature
     call nc_write_real_var1d(ncid,"LandT_K",chnlid,work1d, &
                    varlongname="Jacobian for Land Temperature", varunits="K K^-1")
  end if
  ! Transmittance
  call nc_write_real_var2d(ncid,"ptau",dimid2d,ptau5, &
                    varlongname="Transmittance")

  call check_nc(nf90_close(ncid))
   
  ! ==========================================================================
  ! STEP 9. **** CLEAN UP ****
  !
  ! 9a. Deallocate the structures
  ! -----------------------------
  CALL CRTM_Atmosphere_Destroy(atm_K)
  CALL CRTM_Atmosphere_Destroy(atm)

  ! 9b. Deallocate the arrays
  ! -------------------------
  DEALLOCATE(rts, rts_K, sfc_k, atm_k, chlidx, ptau5, levidx, work1d, work2d, STAT = alloc_stat)

  ! ==========================================================================
  ! 10. **** DESTROY THE CRTM ****
  !
  WRITE( *, '( /5x, "Destroying the CRTM..." )' )
  err_stat = CRTM_Destroy( chinfo )
  IF ( err_stat /= SUCCESS ) THEN
    message = 'Error destroying CRTM'
    CALL Display_Message( PROGRAM_NAME, message, FAILURE )
    STOP
  END IF

  END DO SENSORLOOP

CONTAINS
 SUBROUTINE nc_write_real_var1d(inncid,varname,dimid,vararray,varlongname,varunits)
 IMPLICIT NONE
 INTEGER, INTENT(in) :: inncid, dimid
 CHARACTER(*), INTENT(in) :: varname
 CHARACTER(*), INTENT(in), OPTIONAL :: varlongname, varunits
 REAL(fp), INTENT(in) :: vararray(:)
 
 INTEGER :: varid

   call check_nc( nf90_redef(inncid) )
   call check_nc( nf90_def_var(inncid,varname,NF90_REAL,dimid,varid) )
   if ( present(varlongname) ) then
      call check_nc( nf90_put_att(inncid,varid,"long_name",trim(varlongname)) )
   end if
   if ( present(varunits) ) then
      call check_nc( nf90_put_att(inncid,varid,"units",trim(varunits)) )
   end if
   call check_nc( nf90_enddef(inncid) )
   
   call check_nc( nf90_put_var(inncid,varid,vararray) )

 END SUBROUTINE nc_write_real_var1d
 SUBROUTINE nc_write_real_var2d(inncid,varname,dimid,vararray,varlongname,varunits)
 IMPLICIT NONE
 INTEGER, INTENT(in) :: inncid, dimid(2)
 CHARACTER(*), INTENT(in) :: varname
 CHARACTER(*), INTENT(in), OPTIONAL :: varlongname, varunits
 REAL(fp), INTENT(in) :: vararray(:,:)

 INTEGER :: varid

   call check_nc( nf90_redef(inncid) )
   call check_nc( nf90_def_var(inncid,varname,NF90_REAL,dimid,varid) )
   if ( present(varlongname) ) then
      call check_nc( nf90_put_att(inncid,varid,"long_name",trim(varlongname)) )
   end if
   if ( present(varunits) ) then
      call check_nc( nf90_put_att(inncid,varid,"units",trim(varunits)) )
   end if
   call check_nc( nf90_enddef(inncid) )

   call check_nc( nf90_put_var(inncid,varid,vararray) )

 END SUBROUTINE nc_write_real_var2d

! SUBROUTINE Get_RH()
!   allocate(tsen(nsig))
!   do k=1,nsig
!      tv=tmp(k)*(1.+fv*spfh(k))
!      tsen(k)=tv/(1.+fv*max(0.,spfh(k)))
!   end do
!!
!! Generate qsat, code comes from GSI, genqsat.f90
!!
!   mint=340.
!   lmint=1
!   do k=1,nsig
!      if((prsl(k) < 30000.0 .and.  &
!          prsl(k) > 2000.0) .and.  &
!          tsen(k) < mint)then
!          lmint=k
!          mint=tsen(k)
!      end if
!   end do
!    write(6,*) lmint,mint
!   tdry = mint
!   tr = ttp/tdry
!   if (tdry >= ttp .or. .not. ice) then
!      estmax = psat * (tr**xa) * exp(xb*(1.0-tr))
!   elseif (tdry < tmix) then
!      estmax = psat * (tr**xai) * exp(xbi*(1.0-tr))
!   else
!      w  = (tdry - tmix) / (ttp - tmix)
!      estmax =  w * psat * (tr**xa) * exp(xb*(1.0-tr)) &
!              + (1.0-w) * psat * (tr**xai) * exp(xbi*(1.0-tr))
!   endif
!   allocate(qsat(nsig),rh(nsig))
!   do k=1,nsig
!      tdry = tsen(k)
!      tr = ttp/tdry
!      if (tdry >= ttp .or. .not. ice) then
!         es = psat * (tr**xa) * exp(xb*(1.0-tr))
!      elseif (tdry < tmix) then
!         es = psat * (tr**xai) * exp(xbi*(1.0-tr))
!      else
!         !esw = psat * (tr**xa) * exp(xb*(1.0-tr))
!         !esi = psat * (tr**xai) * exp(xbi*(1.0-tr))
!         w  = (tdry - tmix) / (ttp - tmix)
!         es =  w * psat * (tr**xa) * exp(xb*(1.0-tr)) &
!                  + (1.0-w) * psat * (tr**xai) * exp(xbi*(1.0-tr))
!      end if
!      pw = prsl(k)
!      esmax = es
!      if(lmint > k )then
!         esmax=0.1*pw
!         esmax=min(esmax,estmax)
!      end if
!      es2=min(es,esmax)
!      qsat(k) = eps * es2 / (pw - omeps * es2)
!      rh(k)=spfh(k)/qsat(k)
!      !write(6,*) spfh(k),qsat(k),rh(k)
!   end do
! END SUBROUTINE Get_RH
   function GOCART_Aerosol_size( nbin, itype,  & ! Input
                                       lrh ) & ! Input in 0-1
                           result( R_eff  )   ! in micrometer
   use crtm_aerosolcoeff, only: AeroC, CRTM_AerosolCoeff_Load
   use crtm_module, only:SULFATE_AEROSOL,BLACK_CARBON_AEROSOL,ORGANIC_CARBON_AEROSOL,&
       DUST_AEROSOL,SEASALT_SSAM_AEROSOL,SEASALT_SSCM1_AEROSOL,SEASALT_SSCM2_AEROSOL,SEASALT_SSCM3_AEROSOL
   implicit none
!
!   modified from a function provided by Quanhua Liu
!
   integer,intent(in) :: nbin, itype
   real(8),intent(in) :: lrh

   integer :: j1,j2,k
   integer :: errstat
   real(8) :: h1
   real(8) :: R_eff

   errstat=CRTM_AerosolCoeff_Load("./coefficients/AerosolCoeff.bin")

   if ( itype==DUST_AEROSOL ) then
      if (nbin==1) then
           R_eff = 0.55
      else if (nbin==2) then
           R_eff = 1.4
      else if (nbin==3) then
           R_eff = 2.4
      else if (nbin==4) then
           R_eff = 4.5
      else if (nbin==5) then
           R_eff = 8.0
      end if
      return
   else if ( itype==BLACK_CARBON_AEROSOL .and. nbin==1 ) then
      R_eff = AeroC%Reff(1,itype )
      return
   else if ( itype==ORGANIC_CARBON_AEROSOL .and. nbin==1 ) then
      R_eff = AeroC%Reff(1,itype )
      return
   endif

  j2 = 0

  if ( lrh < AeroC%RH(1) ) then
     j1 = 1
  else if ( lrh > AeroC%RH(AeroC%n_RH) ) then
     j1 = AeroC%n_RH
  else
     do k = 1, AeroC%n_RH-1
        if ( lrh <= AeroC%RH(k+1) .and. lrh > AeroC%RH(k) ) then
           j1 = k
           j2 = k+1
           h1 = (lrh-AeroC%RH(k))/(AeroC%RH(k+1)-AeroC%RH(k))
           exit
        endif
     enddo
  endif
  if ( j2 == 0 ) then
     R_eff = AeroC%Reff(j1,itype )
  else
     R_eff = (1.0-h1)*AeroC%Reff(j1,itype ) + h1*AeroC%Reff(j2,itype )
  endif

  return
  end function GOCART_Aerosol_size
  SUBROUTINE Load_Sfc_Data()

    ! 4a.0 Surface type definitions for default SfcOptics definitions
    !      For IR and VIS, this is the NPOESS reflectivities.
    ! ---------------------------------------------------------------
    INTEGER, PARAMETER :: TUNDRA_SURFACE_TYPE         = 10  ! NPOESS Land surface type for IR/VIS Land SfcOptics
    INTEGER, PARAMETER :: SCRUB_SURFACE_TYPE          =  7  ! NPOESS Land surface type for IR/VIS Land SfcOptics
    INTEGER, PARAMETER :: COARSE_SOIL_TYPE            =  1  ! Soil type for MW land SfcOptics
    INTEGER, PARAMETER :: GROUNDCOVER_VEGETATION_TYPE =  7  ! Vegetation type for MW Land SfcOptics
    INTEGER, PARAMETER :: BARE_SOIL_VEGETATION_TYPE   = 11  ! Vegetation type for MW Land SfcOptics
    INTEGER, PARAMETER :: SEA_WATER_TYPE              =  1  ! Water type for all SfcOptics
    INTEGER, PARAMETER :: FRESH_SNOW_TYPE             =  2  ! NPOESS Snow type for IR/VIS SfcOptics
    INTEGER, PARAMETER :: FRESH_ICE_TYPE              =  1  ! NPOESS Ice type for IR/VIS SfcOptics
    REAL(fp) :: watercover

    ! 4a.1 Profile #1
    ! ---------------
    ! ...Land surface characteristics
    watercover=1.0_fp-landcover
    !write(6,*) 'land,water=',landcover,watercover
    sfc(1)%Land_Coverage     = landcover
    sfc(1)%Land_Type         = landtype
    sfc(1)%Land_Temperature  = 300.0_fp
    sfc(1)%Lai               = lai
    sfc(1)%Soil_Type         = COARSE_SOIL_TYPE
    sfc(1)%Vegetation_Type   = GROUNDCOVER_VEGETATION_TYPE
    ! ...Water surface characteristics
    sfc(1)%Water_Coverage    = watercover
    sfc(1)%Water_Type        = SEA_WATER_TYPE
    sfc(1)%Water_Temperature = 300.0_fp
    ! ...Snow coverage characteristics
    sfc(1)%Snow_Coverage    = 0.0_fp
    sfc(1)%Snow_Type        = FRESH_SNOW_TYPE
    sfc(1)%Snow_Temperature = 265.0_fp
    ! ...Ice surface characteristics
    sfc(1)%Ice_Coverage    = 0.0_fp
    sfc(1)%Ice_Type        = FRESH_ICE_TYPE
    sfc(1)%Ice_Temperature = 269.0_fp

  END SUBROUTINE Load_Sfc_Data

  subroutine check_nc(status)
    integer, intent(in) :: status

    if(status /= nf90_noerr) then
      print *, 'netCDF error!:'
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check_nc

END PROGRAM crtm_ideal
