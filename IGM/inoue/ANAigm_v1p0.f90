! ANAigm.f90
!
! Analytic mean IGM transmission curve based on Inoue, Shimizu, Iwata (2014).
! The analytic approximation is valid if the calculation wavelength is larger than 
! 912 A (Hydrogen Lyman limit) in the observer's rest-frame. Be careful for low-z objects.
!
! 2013 December 28 by AKI
! 2014 February 13 by AKI: "LLS" => "DLA" as in the paper.
! 2014 February 14 by AKI: A typo in tLCLAF(zS,lobs) was corrected.
!----------------------------------------------------------------------
program ANAigm

  implicit none
  integer, parameter :: NzS = 1
  real(8), dimension(40) :: lam1  ! Wavelength of Lyman series lines [A]
  real(8), dimension(40,3) :: ALAF  ! Coefficients for Lyman series optical depth by the LAF comp.
  real(8), dimension(40,2) :: ADLA  ! Coefficients for Lyman series optical depth by the DLA comp.
  real(8) :: lobs, tau, lrest, Zs
  real(8) :: dlest = 1d0  ! A
  ! real(8) :: zS = 3d0
  integer :: i, j, k, GETZ, OpenStatus, InputStatus
  CHARACTER(LEN=100) :: option1
  CHARACTER(LEN=200) :: idlpath, dlainput, lafinput

  !grab idl path 
  CALL GETENV("MIDL",idlpath)
  dlainput = trim(idlpath) // "IGM/inoue/DLAcoeff.txt"
  lafinput = trim(idlpath) // "IGM/inoue/LAFcoeff.txt"
  

  ! Grab the input 
  GETZ = IARGC()
  if (GETZ < 1) then 
     zS=3d0
  else
     CALL GETARG(1,option1)  
     read(option1, *) zS
  endif

     ! Coefficients input
     open (15, file = trim(lafinput), status = "old", iostat = OpenStatus)
     if (OpenStatus > 0) stop "*** Cannot open the file! ***"
     do
        read (15, *, iostat = InputStatus) i, lam1(i), ALAF(i,1), ALAF(i,2), ALAF(i,3)
        if (InputStatus > 0) stop "*** Input error! ***"
        if (InputStatus < 0) exit ! End of File
     end do
     close(15)
     open (15, file = trim(dlainput), status = "old", iostat = OpenStatus)
     if (OpenStatus > 0) stop "*** Cannot open the file! ***"
     do
        read (15, *, iostat = InputStatus) i, lam1(i), ADLA(i,1), ADLA(i,2)
        if (InputStatus > 0) stop "*** Input error! ***"
        if (InputStatus < 0) exit ! End of File
     end do
     close(15)
     
     
     do i = 0, 1300
        lrest = 1d2 + real(i) * dlest
        
        write (*,'(f7.1)',advance='no') lrest
        
        do k = 1, NzS
           
           ! zS = real(k) * dzS
           
           lobs = lrest * (1d0 + zS)
           
           tau = tLSLAF(zS,lobs,lam1,ALAF) + tLSDLA(zS,lobs,lam1,ADLA) &
                + tLCLAF(zS,lobs) + tLCDLA(zS,lobs)
           
           !        tau = tLSLAF(zS,lobs,lam1,ALAF)
           
           !        tau = tLCLAF(zS,lobs) 
           
           !        tau = tLSDLA(zS,lobs,lam1,ADLA)
           
           !        tau = tLCDLA(zS,lobs) 
           
           write (*,'(es12.5)',advance='no') exp(-tau)
        end do
        
        write (*,*)
     end do
     
   contains
     
     ! Lyman series optical depth by the DLA component in Inoue, Shimizu, Iwata (2014)
     ! Input: Source redshift (zS), observers' frame wavelength in A (lobs),
     ! Lyman series wavelengths in A (lam1), and coefficients (A).
     function tLSDLA(zS,lobs,lam1,A)
       
       real(8) :: tLSDLA
       real(8), intent(in) :: zS
       real(8), intent(in) :: lobs
       real(8), intent(in), dimension(40) :: lam1  ! Wavelength of Lyman series lines [A]
       real(8), intent(in), dimension(40,2) :: A  ! Coefficients for Lyman series optical depth
       real(8), parameter :: z1DLA = 2d0
       integer :: j
       
       tLSDLA = 0d0  ! Initialize
       do j = 2, 40
          if ((lobs < lam1(j)*(1d0+zS)) .and. (lobs > lam1(j))) then
             if (lobs < lam1(j)*(1d0+z1DLA)) then
                tLSDLA = tLSDLA + A(j,1) * (lobs / lam1(j)) ** 2d0
             else
                tLSDLA = tLSDLA + A(j,2) * (lobs / lam1(j)) ** 3d0
             end if
          end if
       end do
       
     end function tLSDLA
     
     ! Lyman series optical depth by the LAF component in Inoue, Shimizu, Iwata (2014)
     ! Input: Source redshift (zS), observers' frame wavelength in A (lobs),
     ! Lyman series wavelengths in A (lam1), and coefficients (A).
     function tLSLAF(zS,lobs,lam1,A)
       
       real(8) :: tLSLAF
       real(8), intent(in) :: zS
       real(8), intent(in) :: lobs
       real(8), intent(in), dimension(40) :: lam1  ! Wavelength of Lyman series lines [A]
       real(8), intent(in), dimension(40,3) :: A  ! Coefficients for Lyman series optical depth
       real(8), parameter :: z1LAF = 1.2d0
       real(8), parameter :: z2LAF = 4.7d0
       integer :: j
       
       tLSLAF = 0d0  ! Initialize
       do j = 2, 40
          if ((lobs < lam1(j)*(1d0+zS)) .and. (lobs > lam1(j))) then
             if (lobs < lam1(j)*(1d0+z1LAF)) then
                tLSLAF = tLSLAF + A(j,1) * (lobs / lam1(j)) ** 1.2d0
             else if (lobs < lam1(j)*(1d0+z2LAF)) then
                tLSLAF = tLSLAF + A(j,2) * (lobs / lam1(j)) ** 3.7d0
             else
                tLSLAF = tLSLAF + A(j,3) * (lobs / lam1(j)) ** 5.5d0
             end if
          end if
       end do
       
     end function tLSLAF
     
     ! Lyman continuum optical depth by the DLA component in Inoue, Shimizu, Iwata (2014)
     ! Input: Source redshift (zS) and observers' frame wavelength in A (lobs)
     ! NOTE: This function is valid only for lobs > lamL.
     function tLCDLA(zS,lobs)
       
       real(8) :: tLCDLA
       real(8), intent(in) :: zS
       real(8), intent(in) :: lobs
       real(8), parameter :: z1DLA = 2d0
       real(8), parameter :: lamL = 911.8d0  ! Lyman limit wavelength [A]
       
       if (lobs > lamL*(1d0+zS)) then
          tLCDLA = 0d0
       else if (zS < z1DLA) then
          tLCDLA = 0.2113d0 * (1d0+zS)**2d0 - 0.07661d0 * (1d0+zS)**2.3d0 * (lobs/lamL)**(-3d-1) &
               - 0.1347d0 * (lobs/lamL)**2d0
       else 
          if (lobs > lamL*(1d0+z1DLA)) then
             tLCDLA = 0.04696d0 * (1d0+zS)**3d0 &
                  - 0.01779d0 * (1d0+zS)**3.3d0 * (lobs/lamL)**(-3d-1) &
                  - 0.02916 * (lobs/lamL)**3d0
          else
             tLCDLA = 0.6340d0 + 0.04696d0 * (1d0+zS)**3d0 &
                  - 0.01779d0 * (1d0+zS)**3.3d0 * (lobs/lamL)**(-3d-1) &
                  - 0.1347d0 * (lobs/lamL)**2d0 - 0.2905 * (lobs/lamL)**(-3d-1)
          end if
       end if
       
     end function tLCDLA
     
     ! Lyman continuum optical depth by the LAF component in Inoue, Shimizu, Iwata (2014)
     ! Input: Source redshift (zS) and observers' frame wavelength in A (lobs)
     ! NOTE: This function is valid only for lobs > lamL.
     function tLCLAF(zS,lobs)
       
       real(8) :: tLCLAF
       real(8), intent(in) :: zS
       real(8), intent(in) :: lobs
       real(8), parameter :: z1LAF = 1.2d0
       real(8), parameter :: z2LAF = 4.7d0
       real(8), parameter :: lamL = 911.8d0  ! Lyman limit wavelength [A]
       
       if (lobs > lamL*(1d0+zS)) then
          tLCLAF = 0d0
       else if (zS < z1LAF) then
          tLCLAF = 0.3248d0 * ((lobs/lamL)**1.2d0 - (1d0+zS)**(-9d-1) * (lobs/lamL)**2.1d0)
       else if (zS < z2LAF) then
          if (lobs > lamL*(1d0+z1LAF)) then
             tLCLAF = 2.545d-2 * ((1d0+zS)**1.6d0 * (lobs/lamL)**2.1d0 - (lobs/lamL)**3.7d0)
          else
             tLCLAF = 2.545d-2 * (1d0+zS)**1.6d0 * (lobs/lamL)**2.1d0 &
                  + 0.3248d0 * (lobs/lamL)**1.2d0 - 0.2496d0 * (lobs/lamL)**2.1d0
          end if
       else
          if (lobs > lamL*(1d0+z2LAF)) then
             tLCLAF = 5.221d-4 * ((1d0+zS)**3.4d0 * (lobs/lamL)**2.1d0 - (lobs/lamL)**5.5d0)
          else if (lobs > lamL*(1d0+z1LAF)) then
             tLCLAF = 5.221d-4 * (1d0+zS)**3.4d0 * (lobs/lamL)**2.1d0 &
                  + 0.2182d0 * (lobs/lamL)**2.1d0 - 2.545d-2 * (lobs/lamL)**3.7d0
          else
             tLCLAF = 5.221d-4 * (1d0+zS)**3.4d0 * (lobs/lamL)**2.1d0 &
                  + 0.3248d0 * (lobs/lamL)**1.2d0 - 3.140d-2 * (lobs/lamL)**2.1d0
          end if
       end if
       
     end function tLCLAF
    
end program ANAigm

