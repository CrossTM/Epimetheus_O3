!********************************************************************************************************************************
  SUBROUTINE READ_ARGUMENT(INPUT_FILE,OUTPUT_FILE)
 !********************************************************************************************************************************
 ! READ INPUT AND OUTPUT FILE NAME GIVEN IN ARGUMENT
 ! IF THERE IS NOT OUTPUT FILE NAME CREATE ONE USING INPUT FILE NAME MORE ".OUT"
 !
 ! INPUT_FILE    : STRING, NAME OF THE INPUT  FILE   (OUTPUT)
 ! OUTPUT_FILE   : STRING, NAME OF THE OUTPUT FILE   (OUTPUT)
 !
 ! WARNING: HAVE TO CHECK COMPATIBILITY WITH MPI FOR TOTAL NUMBER OF ARGUMENT (MAXIMUM ALLOWED BY THIS SUBROUTINE: 2 )
 !
 ! BRUNO MADEBENE                                                                                    LAST MODIFICATION :  5/06/2006
 ! DAVID BENOIT: FEW COSMETIC CHANGES AND CHANGED .OUT to .out                                       LAST MODIFICATION : 14/09/2006
 ! DAVID BENOIT: REMOVED AUTOMATIC GENERATION OF OUTPUT FILE, GOES TO STDOUT INSTEAD                 LAST MODIFICATION : 24/11/2006
 !*********************************************************************************************************************************
 IMPLICIT NONE

 CHARACTER, INTENT(OUT)  :: INPUT_FILE*(*),OUTPUT_FILE*(*)   ! NAME OF THE INPUT AND OUTPUT FILE

 INTEGER                 :: NB_ARGUMENT                      ! NUMBER OF ARGUMENT ON THE COMMANDE LINE
 INTEGER                 :: IARGC

 ! FIND NUMBER OF ARGUMENT GIVEN BY USER
 NB_ARGUMENT= IARGC()
 ! CHECK WHAT TO DO DEPENDING ON THE NUMBER OF ARGUMENT
 SELECT CASE (NB_ARGUMENT)
   CASE(0)
   ! NO ARGUMENT, PROBLEM: NO INPUT FILE, PRINT ERROR MESSAGE ON SCREEN AND STOP PROGRAM
     WRITE(*,*)
     WRITE(*,*) ' ERROR: NO INPUT FILE'
     WRITE(*,*) ' EXITING PROGAM'
     WRITE(*,*)
     STOP
   CASE(1)
   ! ONLY ONE ARGUMENT, ASSUMED TO BE INPUT FILE
     CALL GETARG(1,INPUT_FILE)
     WRITE(*,*)
     WRITE(*,*) ' INPUT  FILE        : ',TRIM(ADJUSTL(INPUT_FILE))
     WRITE(*,*)
   CASE(2)
   ! TWO ARGUMENT, FIRST IS ASSUME TO BE INPUT_FILE NAME, SECOND, OUTPUT FILE NAME
     CALL GETARG(1,INPUT_FILE)
     CALL GETARG(2,OUTPUT_FILE)
     WRITE(*,*)
     WRITE(*,*) ' INPUT  FILE        : ',TRIM(ADJUSTL(INPUT_FILE))
     WRITE(*,*) ' OUTPUT FILE        : ',TRIM(ADJUSTL(OUTPUT_FILE))
     WRITE(*,*)
   ! MORE THAN TWO ARGUMENT: LET IT SLIP, GRABBING THE INPUT FILE NAME AND ASSUMING THAT THE REST IS OUTPUT REDIRECTION
   CASE DEFAULT
     CALL GETARG(1,INPUT_FILE)
     WRITE(*,*)
     WRITE(*,*) ' INPUT  FILE        : ',TRIM(ADJUSTL(INPUT_FILE))
     WRITE(*,*)
 END SELECT

 END SUBROUTINE READ_ARGUMENT
 !********************************************************************************************************************************


program ozonepes
implicit none

! NEEDS:  r1/angs, r2/angs, theta/rad
double precision :: r1, r2, theta,v
CHARACTER(256) :: inputfile,outputfile

! get filename
CALL READ_ARGUMENT(inputfile, outputfile)

open(unit=1,form='formatted',file=TRIM(ADJUSTL(inputfile)))
read(1,*) r1, r2, theta
!write(*,*) r1,r2,theta
close(unit=1)

! theta is already in radians...
! CONVERT THETA IN RADIANS...
!theta=theta*dacos(-1.D0)/180.D0

write(*,*) "r1=",r1, "r2=",r2, "theta=",theta

!distances have to be in ANGS for this to work
call poten(v,r1,r2,theta)

write(*,*) " energy/cm-1 ",v

! Save energy to file in hartree
open(unit=2,form='formatted',file='ENERGY')
write(2,*) v/219474.624d0
close(unit=2)

end program ozonepes


