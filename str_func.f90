CHARACTER(LEN=80) FUNCTION CenterStr(str)

! Centers the input string within an 80 character output string.

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  CHARACTER(LEN=80) :: str   ! the string input
  INTEGER :: length          ! length of the string without trailing blanks
  INTEGER :: midpoint        ! the midpoint in the string
  INTEGER :: i

  IF (TEST3) PRINT *,"CenterStr" !TEST3

  str=TRIM(str)
  length=LEN_TRIM(str)
  midpoint=(length/2)+1

  CenterStr=""

  DO i=1,length
    CenterStr(40-midpoint+i:40-midpoint+i)=str(i:i)
  END DO

  CenterStr(1:1) = '|'
  CenterStr(80:80) = '|'

END FUNCTION CenterStr
SUBROUTINE ComplStr(str)
!
! Returns the DNA-complement of the section of a string
!
!
!                 A .................... = T
!                 C .................... = G
!                 G .................... = C
!                 T .................... = A
!                 M = A or C ........... = K
!                 R = A or G ........... = Y
!                 W = A or T ........... = W
!                 S = C or G ........... = S
!                 Y = C or T ........... = R
!                 K = G or T ........... = M
!                 V = A or C or G ...... = B
!                 H = A or C or T ...... = D
!                 D = A or G or T ...... = H
!                 B = C or G or T ...... = V
!                 N = A or C or G or T . = N

  USE dnaworks_test
  IMPLICIT NONE

  CHARACTER(LEN=*) :: str
  INTEGER :: i

  IF (TEST3) PRINT *,'ComplStr'

  DO i=1,LEN_TRIM(str)
    SELECT CASE(str(i:i))
      CASE('A')
        str(i:i)="T"
      CASE('T')
        str(i:i)="A"
      CASE('G')
        str(i:i)="C"
      CASE('C')
        str(i:i)="G"
      CASE('M')
        str(i:i)="K"
      CASE('R')
        str(i:i)="Y"
      CASE('W')
        str(i:i)="W"
      CASE('S')
        str(i:i)="S"
      CASE('Y')
        str(i:i)="R"
      CASE('K')
        str(i:i)="M"
      CASE('V')
        str(i:i)="B"
      CASE('H')
        str(i:i)="D"
      CASE('D')
        str(i:i)="H"
      CASE('B')
        str(i:i)="V"
      CASE('N')
        str(i:i)="N"
    END SELECT
  END DO

END SUBROUTINE ComplStr
INTEGER FUNCTION NT2Int(nt)
!
! Converts a nt into an integer.

  USE dnaworks_test
  IMPLICIT NONE

  CHARACTER(LEN=1) :: nt

  IF (TEST3) PRINT *,"NT2Int" !TEST3

  SELECT CASE(nt(1:1)) ! convert sequence to num representation
    CASE('A')
      NT2Int=-1
    CASE('T')
      NT2Int=1
    CASE('C')
      NT2Int=-3
    CASE('G')
      NT2Int=3
  END SELECT

END FUNCTION NT2Int
SUBROUTINE RevComplStr(str)
!
! Returns the reverse complement of a string

  USE dnaworks_test
  IMPLICIT NONE

  CHARACTER(LEN=*) :: str

  IF (TEST3) PRINT *,'RevComplStr'

  CALL RevStr(str)
  CALL ComplStr(str)

END SUBROUTINE RevComplStr
SUBROUTINE RevStr(str)
!
! Returns the reverse of a string

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  CHARACTER(LEN=*) :: str
  INTEGER :: i,j

  IF (TEST3) PRINT *,'RevStr'

  DO i=1,LEN_TRIM(str)
    j=(LEN_TRIM(str))-i+1
    SCRATCH(i:i)=str(j:j)
  END DO

  str=SCRATCH(1:(LEN_TRIM(str)))

END SUBROUTINE RevStr
INTEGER FUNCTION StrToInt(str)
!
! Converts a string into an integer.  All spaces are assigned zero.

  USE dnaworks_test
  IMPLICIT NONE

  CHARACTER(LEN=*) :: str
  INTEGER :: i,strlen,val,a,j

  IF (TEST3) PRINT *,"StrToInt"

! initial values

  j=1
  StrToInt=0

! find length of string

  strlen=LEN(str)

! convert all non-numerical characters to spaces, except for '-' sign

  DO i=1,strlen
    a=IACHAR(str(i:i))
    IF ((a.ge.48.and.a.le.57).or.(str(i:i).eq."-")) THEN
      str(i:i)=str(i:i)
    ELSE
      str(i:i)=" "
    END IF
  END DO

! shift string to right

  str=ADJUSTR(str)

! convert to integer,going from right to left

  DO i=strlen,1,-1
    a=IACHAR(str(i:i))
    IF (a.ge.48.and.a.le.57) THEN
      val=a-48
    ELSE
      val=0
    END IF
    StrToInt=(val*j)+StrToInt
    j=j*10
  END DO

! find sign

  IF ((INDEX(str,'-')).gt.0) StrToInt=StrToInt*(-1)

END FUNCTION StrToInt
REAL FUNCTION StrToReal(str)
!
! Converts a string into an real number.  Blanks are ignored.

  USE dnaworks_test
  IMPLICIT NONE

  CHARACTER(LEN=*) :: str
  INTEGER :: i,strlen,a
  INTEGER :: de
  INTEGER,EXTERNAL :: StrToInt
  REAL :: front_real,back_real,ex_real,de_real,si
  CHARACTER(LEN=1) :: b

  IF (TEST3) PRINT *,"StrToReal" !TEST3

! initial values

  StrToReal=0
  front_real=0
  back_real=0
  si=1

! convert all useless characters to spaces

  DO i=1,LEN(str)
    b=str(i:i)
    a=IACHAR(b)
    IF ((a.ge.48.and.a.le.57).or.(b.eq."-").or.(b.eq.'.').or.&
        (b.eq.'e').or.(b.eq.'E')) THEN
      str(i:i)=str(i:i)
    ELSE
      str(i:i)=" "
    END IF
  END DO

! shift string to left

  str=ADJUSTL(str)
  strlen=LEN_TRIM(str)

! find sign

  IF (str(1:1).eq.'-') si=-1

! find exponents, if there

  IF ((INDEX(str,'e')).gt.0) THEN
    ex_real = (REAL(10))**(REAL(StrToInt(str(((INDEX(str,'e'))+1):strlen))))
    strlen = (INDEX(str,'e'))-1
  ELSE IF ((INDEX(str,'E')).gt.0) THEN

    ex_real = (REAL(10))**(REAL(StrToInt(str(((INDEX(str,'E'))+1):strlen))))
    strlen = (INDEX(str,'E'))-1
  ELSE
    ex_real=1
  END IF

! find decimal position

  de=INDEX(str,'.')

! find values

  IF (de.gt.0) THEN
    front_real = REAL(ABS(StrToInt(str(1:(de-1)))))
    de_real=(REAL(1))/((REAL(10))**(strlen-de))
    back_real = (REAL(ABS(StrToInt(str((de+1):strlen)))))*de_real
  ELSE
    front_real = REAL(ABS(StrToInt(str(1:strlen))))
    back_real=0
  END IF

  StrToReal=(front_real+back_real)*si*ex_real

END FUNCTION StrToReal
SUBROUTINE ToLowerCase(str)

! Converts a string to lower case

  USE dnaworks_test
  IMPLICIT NONE

  CHARACTER(LEN=*) :: str
  INTEGER :: i,j

  IF (TEST3) PRINT *,'ToLowerCase'

  DO i=1,LEN_TRIM(str)
    j = ICHAR(str(i:i))
    IF (65.LE.j.AND.j.LE.90) str(i:i)=CHAR(j+32)
  END DO

END SUBROUTINE ToLowerCase
SUBROUTINE ToUpperCase(str)

! Converts a string to upper case

  USE dnaworks_test
  IMPLICIT NONE

  CHARACTER(LEN=*) :: str
  INTEGER :: i,j

  IF (TEST3) PRINT *,'ToUpperCase'

  DO i=1,LEN_TRIM(str)
    j = ICHAR(str(i:i))
    IF (97.LE.j.AND.j.LE.122) str(i:i)=CHAR(j-32)
  END DO

END SUBROUTINE ToUpperCase
