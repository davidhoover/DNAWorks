CHARACTER(LEN=10) FUNCTION CurrentDate()
!
! Returns the date as DD/MM/YYYY in a 10 character string

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: values(8)
  CHARACTER(LEN=8) :: date
  CHARACTER(LEN=10) :: time
  CHARACTER(LEN=5) :: zone
  CHARACTER(LEN=2) :: MM,DD
  CHARACTER(LEN=4) :: YYYY

  IF (TEST3) PRINT *,"CurrentDate" !TEST3

  CALL DATE_AND_TIME(date,time,zone,values)

! Capture months

  IF (values(2).GE.10) THEN
    WRITE(UNIT=MM,FMT="(i2)") values(2)
  ELSE
    WRITE(UNIT=MM,FMT="('0',i1)") values(2)
  END IF

! Capture days

  IF (values(3).GE.10) THEN
    WRITE(UNIT=DD,FMT="(i2)") values(3)
 ELSE
    WRITE(UNIT=DD,FMT="('0',i1)") values(3)
  END IF

! Capture year

  WRITE(UNIT=YYYY,FMT="(i4)") values(1)

  CurrentDate = MM//'/'//DD//'/'//YYYY

END FUNCTION CurrentDate
CHARACTER(LEN=6) FUNCTION CurrentDateNice()
!
! Returns the date as YYMMDD in an 6 character string

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: values(8)
  CHARACTER(LEN=8) :: date
  CHARACTER(LEN=10) :: time
  CHARACTER(LEN=5) :: zone
  CHARACTER(LEN=2) :: MM,DD
  CHARACTER(LEN=4) :: YYYY

  IF (TEST3) PRINT *,"CurrentDateNice" !TEST3

  CALL DATE_AND_TIME(date,time,zone,values)

! Capture months

  IF (values(2).GE.10) THEN
    WRITE(UNIT=MM,FMT="(i2)") values(2)
  ELSE
    WRITE(UNIT=MM,FMT="('0',i1)") values(2)
  END IF

! Capture days

  IF (values(3).GE.10) THEN
    WRITE(UNIT=DD,FMT="(i2)") values(3)
 ELSE
    WRITE(UNIT=DD,FMT="('0',i1)") values(3)
  END IF

! Capture year

  WRITE(UNIT=YYYY,FMT="(i4)") values(1)

  CurrentDateNice = YYYY(3:4)//MM//DD

END FUNCTION CurrentDateNice
CHARACTER(LEN=8) FUNCTION CurrentTime()
!
! Returns the time as HH:MM:SS in a 8 character string

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: values(8)
  CHARACTER(LEN=8) :: date
  CHARACTER(LEN=10) :: time
  CHARACTER(LEN=5) :: zone
  CHARACTER(LEN=2) :: HH,MM,SS

  IF (TEST3) PRINT *,"CurrentTime" !TEST3

  CALL DATE_AND_TIME(date,time,zone,values)

! Capture hours

  IF (values(5).GE.10) THEN
    WRITE(UNIT=HH,FMT="(i2)") values(5)
  ELSE IF (values(5).EQ.0) THEN
    WRITE(UNIT=HH,FMT="('00')")
  ELSE
    WRITE(UNIT=HH,FMT="('0',i1)") values(5)
  END IF

! Capture minutes

  IF (values(6).GE.10) THEN
    WRITE(UNIT=MM,FMT="(i2)") values(6)
  ELSE IF (values(6).EQ.0) THEN
    WRITE(UNIT=MM,FMT="('00')")
  ELSE
    WRITE(UNIT=MM,FMT="('0',i1)") values(6)
  END IF

! Capture seconds

  IF (values(7).GE.10) THEN
    WRITE(UNIT=SS,FMT="(i2)") values(7)
  ELSE IF (values(7).EQ.0) THEN
    WRITE(UNIT=SS,FMT="('00')")
  ELSE
    WRITE(UNIT=SS,FMT="('0',i1)") values(7)
  END IF

  CurrentTime = HH//':'//MM//':'//SS

END FUNCTION CurrentTime
INTEGER FUNCTION CurrentTimeSeconds()
!
! Returns the time in seconds since Wed Dec 31 19:00:00 1969, adjusting for
! leap years.  However, there is NO adjustment for daylight saving time

  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: values(8),x
  CHARACTER(LEN=8) :: date
  CHARACTER(LEN=10) :: time
  CHARACTER(LEN=5) :: zone
  LOGICAL :: leap_year

  IF (TEST3) PRINT *,"CurrentTimeSeconds" !TEST3
!  PRINT *,"CurrentTimeSeconds" !TEST3

  leap_year=.FALSE.

! x = number of seconds from 12/31/1969 19:00:00 to 1/1/2001 00:00:00

  x=978325200

! find current date and time

  CALL DATE_AND_TIME(date,time,zone,values)

! is this year a leap year?

  IF (MOD(values(1),4).eq.0) leap_year=.TRUE.

! find yearly sums

  x=x+((values(1)-2001)*(365*86400))+(((values(1)-2001)/4)*86400)

! find monthly sums

  SELECT CASE(values(2))
    CASE(2)
      x=x+(31*86400)
    CASE(3)
      x=x+(59*86400)
    CASE(4)
      x=x+(90*86400)
    CASE(5)
      x=x+(120*86400)
    CASE(6)
      x=x+(151*86400)
    CASE(7)
      x=x+(181*86400)
    CASE(8)
      x=x+(212*86400)
    CASE(9)
      x=x+(243*86400)
    CASE(10)
      x=x+(273*86400)
    CASE(11)
      x=x+(304*86400)
    CASE(12)
      x=x+(334*86400)
  END SELECT

! correct for leap day in February

  IF ((leap_year).and.(values(2).gt.2)) x=x+86400

! are we in daylight savings time?

!  spring=7-(MOD((2800+values(1)-2474+(INT(values(1)/4))),7))
!  fall=spring+20
!  IF (spring.le.4) fall=fall+7
!
!  IF ((values(2).ge.4).and.(values(2).le.10).and.&
!     &(values(3).ge.spring).and.(values(3).le.fall).and.&
!     &(values(5).ge.2)) x=x-3600

! find final sum

    CurrentTimeSeconds=x+((values(3)-1)*86400)+(values(5)*3600)+(values(6)*60)+values(7)
!  PRINT *,CurrentTimeSeconds !TEST3

END FUNCTION CurrentTimeSeconds
