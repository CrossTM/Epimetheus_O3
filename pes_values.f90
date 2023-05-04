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

! NEEDS:  r1/bohrs, r2/bohrs, cos( theta )
double precision :: r1, r2, cos_theta,theta,v
double precision :: A1(3),A2(3),A3(3)
double precision :: v1(3),v2(3)
double precision :: n1,n2
integer :: natoms
CHARACTER(256) :: inputfile,outputfile,title
CHARACTER(3) :: atomname

! get filename
CALL READ_ARGUMENT(inputfile, outputfile)

!inputfile='test.xyz'

open(unit=1,form='formatted',file=TRIM(ADJUSTL(inputfile)))
read(1,*) natoms
write(*,*) "reading", natoms,"atoms"

!check that we have 3 atoms
if (natoms .NE. 3) STOP

read(1,*) title

! We know that there are 3 atoms so a completely unrolled loop will do...
! ASSUMING ANGTROMS
read(1,*) atomname,A1(1),A1(2),A1(3) 
read(1,*) atomname,A2(1),A2(2),A2(3) 
read(1,*) atomname,A3(1),A3(2),A3(3) 

close(unit=1)

write(*,*) 'molecule name: ',TRIM(ADJUSTL(title))
write(*,*) A1
write(*,*) A2
write(*,*) A3

! computing distance vectors

v1(1)=(A1(1)-A2(1))
v1(2)=(A1(2)-A2(2))
v1(3)=(A1(3)-A2(3))

v2(1)=(A3(1)-A2(1))
v2(2)=(A3(2)-A2(2))
v2(3)=(A3(3)-A2(3))

!norm
n1=dsqrt(v1(1)**2+v1(2)**2+v1(3)**2)
n2=dsqrt(v2(1)**2+v2(2)**2+v2(3)**2)

!compute distance
!r1=n1/0.529177208D0
!r2=n2/0.529177208D0
!distances have to be in ANGS for this to work
r1=n1
r2=n2

!compute cosine of angle
cos_theta=dot_product(v1,v2)/(n1*n2)
theta=dacos(cos_theta)
! APPARENTLY THETA STAYS IN RADIANS...
write(*,*) 'angle in degrees',theta*180./acos(-1.)
!r1=0.730
!r2=0.730
!cos_theta=0.

write(*,*) "r1=",r1, "r2=",r2, "cos(t)=",cos_theta,"theta=",theta

call poten(v,r1,r2,theta)

write(*,*) TRIM(ADJUSTL(title))," energy/cm-1 ",v

! Save energy to file in hartree
open(unit=2,form='formatted',file='ENERGY')
write(2,*) v/219474.624d0
close(unit=2)

end program ozonepes


