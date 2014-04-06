SUBROUTINE Send_Email

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER,EXTERNAL :: CurrentTimeSeconds
  INTEGER :: start
  CHARACTER(LEN=1000) :: text
  CHARACTER(LEN=500) :: text1,text2

  IF (TEST0) PRINT *,"Send_Email" !TEST0

  WRITE(text1,FMT="(a80,' -s ""DNAWorks Output - ',a80)") MAILPATH,jobname
  WRITE(text2,FMT="('"" <',a,' ',a80)") outputfile,email

  text=text1(1:LEN_TRIM(text1))//text2(1:LEN_TRIM(text2))

!  PRINT *,text

  CALL SYSTEM(text)

! The following is a waste of time.  It should take about 10 seconds to
! go through the loop. This should give the program enough time to send out
!  an email.

  start=CurrentTimeSeconds()
  DO WHILE (CurrentTimeSeconds()-start.LT.10)
  END DO

END SUBROUTINE Send_Email
