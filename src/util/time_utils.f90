module time_utils

  use nrtype
  implicit none

  public :: date_extractor
  public :: juldayss
  public :: caldatss

  contains

    subroutine date_extractor(refDate,iy,im,id,ih)
    ! used to extract the date from a units string
    ! (based on a routine written by David Rupp)
    IMPLICIT NONE
    ! input
    character(*),intent(in)             :: refDate     ! reference time string
    ! output
    integer(i4b),intent(out)            :: iy,im,id,ih ! reference time (year,month,day,hour)
    ! internal
    character(len=64)                   :: refd        ! temporary time and units string
    character(len=4)                    :: cyyyy       ! char year extracted from refDate
    character(len=2)                    :: cmm,cdd,chh ! char month and day and hour extracted from refDate
    integer(i4b)                        :: posit       ! position in refDate string

    ! strip out time units, if they exist (seconds since , days since , hours since )
    REFD = TRIM(REFDATE)
    POSIT = INDEX(REFDATE, 'since')
    IF (POSIT.GT.0) REFD = REFD(POSIT+6:len(refD))  ! +6 because 'since' has 5 characters
    ! get the year
    POSIT = INDEX(REFD, '-')     ! up to -
    CYYYY = REFD(1:POSIT-1)
    ! get the month
    REFD = REFD(POSIT+1:len(refD))
    POSIT = INDEX(REFD, '-')     ! up to -
    CMM = REFD(1:POSIT-1)
    ! get the day
    REFD = REFD(POSIT+1:len(refD))
    POSIT = INDEX(REFD, ' ')     ! up to space
    CDD = REFD(1:POSIT-1)
    ! get the hour
    REFD = REFD(POSIT+1:len(refD))
    POSIT = INDEX(REFD, ':')    ! up to :
    IF (POSIT.GT.0) THEN
     CHH = REFD(1:POSIT-1)
    ELSE
     CHH = '00'
    ENDIF
    ! convert to integers
    READ(CYYYY,'(i4)') IY
    READ(CMM,'(i2)') IM
    READ(CDD,'(i2)') ID
    READ(CHH,'(i2)') IH

    end subroutine date_extractor

    SUBROUTINE juldayss(iy,im,id,ih,            &  ! input
                        juldayFrac,ierr,message)   ! output
    !-----------------------------------------------------------------------
    ! compute the julian day
    !-----------------------------------------------------------------------
    !       Based on code from "Numerical Recipes in Fortran-77
    !
    !       Ref:    Press, W.H., S.A. Teukolsky, W.T. Vetterling, and B.P.
    !                 Flannery, 1992:  Numerical Recipes in Fortran 77:
    !                 The Art of Scientific Computing (2nd Ed.)  Cambridge
    !                 University Press, 933pp.
    !-----------------------------------------------------------------------
    !   Modified by David Rupp 2006-March-07 to account for hours with
    !   Output julian time in units of seconds from date
    !-----------------------------------------------------------------------
    IMPLICIT NONE
    ! input
    integer(i4b),intent(in)        :: iy,im,id,ih          ! year, month, day, hour
    ! output
    real(wp),intent(out)           :: juldayFrac           ! julian day (fraction of days)
    integer(i4b), intent(out)      :: ierr                 ! error code
    character(*), intent(out)      :: message              ! error message
    ! internal
    integer(i4b)                   :: julday               ! julian day (whole days)
    integer(i4b),parameter         :: IGREG=15+31*(10+12*1582) !IGREG = 588829
    integer(i4b)                   :: ja,jm,jy             ! temporary
    !-----------------------------------------------------------------------
    ! initialize error control
    !-----------------------------------------------------------------------
    ierr=0; message='juldayss/'
    !-----------------------------------------------------------------------
    jy=iy
    if (jy.eq.0)then; ierr=20; message=trim(message)//'there is no year zero'; return; endif
    if (jy.lt.0) jy=jy+1
    if (im.gt.2) then
     jm=im+1
    else
     jy=jy-1
     jm=im+13
    endif
    julday=int(365.25*jy)+int(30.6001*jm)+id+1720995
    if (id+31*(im+12*iy).ge.IGREG) then
     ja=int(0.01*jy)
     julday=julday+2-ja+int(0.25*ja)
    endif
    juldayFrac = real(julday, KIND(sp) ) &
               + real(ih,  KIND(sp) )/24._wp
    end SUBROUTINE juldayss

    SUBROUTINE caldatss(juliandd,iyyy,im,id,ih,imin,asec)
    !-----------------------------------------------------------------------
    !       Based on code from "Numerical Recipes in Fortran-77
    !
    !       Ref:    Press, W.H., S.A. Teukolsky, W.T. Vetterling, and B.P.
    !                 Flannery, 1992:  Numerical Recipes in Fortran 77:
    !                 The Art of Scientific Computing (2nd Ed.)  Cambridge
    !                 University Press, 933pp.
    !-----------------------------------------------------------------------
    !       Modified by David Rupp 2006-March-07 to account for hours
    !       Output is yyyy, mm, dd, hh
    !-----------------------------------------------------------------------
    IMPLICIT NONE
    ! input
    real(wp),intent(in)       :: juliandd   ! julian day in days since the beginning of time
    ! output
    integer(i4b),intent(out)  :: iyyy,im,id ! day of the year
    integer(i4b),intent(out)  :: ih,imin    ! time of day
    real(wp),intent(out)      :: asec       ! second
    ! internal
    integer(i4b)              :: julian     ! julian day
    integer(i4b),parameter    :: IGREG=2299161
    real(wp)                  :: hours, minutes
    integer(i4b)              :: ja,jalpha,jb,jc,jd,je
    ! get the julian day
    julian = int(juliandd)
    ! gets the hours, (remaining decimal)*24
    hours = (juliandd-julian)*24
    ih = int(hours)  ! convert to an integer
    ! get the minutes, (remaining decimal)*60
    minutes = (hours-ih)*60
    imin = int(minutes)
    ! get the seconds (keep as a decimal
    asec = (minutes-imin)*60
    ! uses the integer julian from above (below original num rec)
    if(julian.ge.IGREG)then
      jalpha=int(((julian-1867216)-0.25)/36524.25)
      ja=julian+1+jalpha-int(0.25*jalpha)
    else if(julian.lt.0)then
      ja=julian+36525*(1-julian/36525)
    else
      ja=julian
    endif
    jb=ja+1524
    jc=int(6680.+((jb-2439870)-122.1)/365.25)
    jd=365*jc+int(0.25*jc)
    je=int((jb-jd)/30.6001)
    id=jb-jd-int(30.6001*je)
    im=je-1
    if(im.gt.12)im=im-12
    iyyy=jc-4715
    if(im.gt.2)iyyy=iyyy-1
    if(iyyy.le.0)iyyy=iyyy-1
    if(julian.lt.0)iyyy=iyyy-100*(1-julian/36525)
    END SUBROUTINE caldatss

end module time_utils
