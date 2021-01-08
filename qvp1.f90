 subroutine qvp(zdr,ranger,itimes,elev_deg,numzdrgates,nazm_zdrs,h_radar,xlat,xlon, &
                year,mon,day,time,depth)
!
! basic subroutine to calculate quasi-vertical profile (QVP) of ZDR observations
! from a time series of a single multiple elevation angle scan.  The QVP is used 
! to estimate convective boundary layer (CBL) depth.  The best elevation angle 
! to use is ~4.5 degrees.  
!
! Based upon the following:
!
! Original Matthew Kumjian MATLAB code
! 16-18 October 2012; Boulder CO
!
! Modified for use with QVPs (Kumjian and Lombardo 2017, MWR) in Feb 2015,
! State College, PA -- mk
!
! Code to identify CBL top from John Banghoff, MS thesis, Penn State 2019
!
! Smoothing and sunrise/sunset added by David Stensrud, 2019-2020
!
! Commented for use in EMC -- 7/6/2020 - mk, ds.  Converted to Fortran90 9/2020 ds
!
! It is assumed that the ZDR observations provided begin no later than 0900 UTC and extend 
! through at least a few hours after sunrise in order to capture the daytime CBL rise. 
! Program could be run every hour to update the time series throughout the daytime.
! No need to run the program at night.  
!
!  Input:
!  zdr           ZDR values in 2d array (times, azimuth, radial) 
!  ranger        radial range distances from radar (m) in 1d array (radial)
!  itimes        number of radar scans provided (should be sequential in time)
!  elev_deg      radar elevation angle of ZDR observations in 1d array (itimes)
!  numzdrgates   number of observations along radial in sweep in 1d array (itimes)
!  nazm_zdrs     number of azimuth angles in sweep in 1d array (itimes)
!  h_radar       radar height above the ground (m)
!  xlat          latitude of radar location 
!  xlon          longitude of radar location
!  year          year of radar observation
!  mon           month of radar observation
!  day           day of radar observation
!  time          utc time of radar observation in 1d array (itimes)
!                format is hour.fraction of hour.  So a value of 6.5 = 6:30
!
!  Output:
!
! Assumed that there are no more than 720 times (elevation scans) in the observations for a given day.
! and that there are no more than 2000 gates along a radial.
!
!  depth         array of CBL depths (m) in 1d array (itimes)
!
  real :: zdr(itimes,nazm_zdr_max,numzdrgate_max)
  real :: elev_deg(itimes)
  real :: numzdrgates(itimes)
  real :: nazm_zdrs(itimes)
  real :: time(itimes)
  real :: depth(itimes)
  real :: ranger(numzdrgate_max)
  real :: height(720,2000)
  real :: qvp_zdr(720,2000)
  real :: slice(7)
  integer :: mond(12)
!
! define constants
!
!  a           earth equatorial radius in m
!  b           earth polar radius in m
!  ke          constant for 4/3 equivalent earth radius assumption
!  pi          3.1415...
!  mond        days in each month of the year, non-leap year months
!  solar       solar constant (W m-2)
!  dtr         conversion factor degrees to radians
!
  ke = 4.0 / 3.0
  a = 6378137.
  b = 6356752.3 
  pi = 4.0*atan(1.0)
  mond = (/31,28,31,30,31,30,31,31,30,31,30,31/)
  solar=1368.
  dtr=pi/180.
!
! print out input values to check
!
  print *,' QVP code: xlat and xlon ', xlat, xlon
!
! loop over all data times
!
 do it=1,itimes
!
! calculate height of radar beam above ground for all ranges
! for this particular scan
!
  elrad = pi*elev_deg(it)/180.

  iend=numzdrgates(it)
  do i=1,iend
   height(it,i)=sqrt(ranger(i)*ranger(i) + (ke*a)*(ke*a) + 2.*ranger(i)*ke*a*sin(elrad))  &
                  -ke*a + h_radar
  enddo
!
! test height values in code
!
  if(it == 1) then
    dh = height(it,40)-height(it,39)
    print *,' QVP code:  Change in height ',dh
  end if

!
! calculate mean zdr for each range gate (average over all azimuth angles)
!      no quality control is done here - assumed that the data are qc'ed before 
!      passing to subroutine
!
  iend=numzdrgates(it)
  kend=nazm_zdrs(it)
  do i=1,iend
   do k=1,kend
    qvp_zdr(it,i) = qvp_zdr(it,i)+zdr(it,k,i)
   enddo
    qvp_zdr(it,i)=qvp_zdr(it,i)/nazm_zdrs(it)
  enddo
!
! end itimes loop - have created a time-height profile of mean zdr from all available observation
! times.
!
 enddo
! 
! We have created the time-height array of ZDR values that constitute the QVP.  Now we smooth
! the array to remove some of the noise. The smoother is a 15-point version that worked well for
! a variety of test cases, but is imperfect as it updates the array as the smoothing is ongoing.
! This is done to save memory but could be updated as needed.  Not sure it makes much difference
! to the result.
!
! Assumption:  number of range gates is the same for all scans
!              if not true, then find minimum number of range gates for all scans and use it in
!              the loop below.  
!
  lend=itimes-2
  mend=numzdrgates(1)-2
  do l=3,lend
    do m = 3,mend
        part0 = qvp_zdr(l-2,m-1)+qvp_zdr(l-2,m)+qvp_zdr(l-2,m+1);
        part1 = qvp_zdr(l-1,m-1)+qvp_zdr(l-1,m)+qvp_zdr(l-1,m+1);
        part2 = qvp_zdr(l,m-1)+qvp_zdr(l,m)+qvp_zdr(l,m+1);
        part3 = qvp_zdr(l+1,m-1)+qvp_zdr(l+1,m)+qvp_zdr(l+1,m+1);
        part4 = qvp_zdr(l+2,m-1)+qvp_zdr(l+2,m)+qvp_zdr(l+2,m+1);
        qvp_zdr(l,m) = (1./15.)*(part0+part1+part2+part3+part4);
    enddo
  enddo
!
! Next step is to calculate the CBL depth as a function of time, which is indicated by the minimum value
! of ZDR in the column at each observation time.  One then connects the height of the ZDR minimum values
! for each observation time to create the CBL height as a function of time.  
!
! To improve algorithm efficiency and accuracy, it is important to know the latitude, longitude, and UTC
! time of the observations.  This allows one to determine is the observation is after sunrise and before
! sunset.  
!
! Input:  year (4 digit integer)
!         month (2 digit integer)
!         day (2 digit integer)
!
! Check for leap year
!
  wtest = year/4.
  ytest = nint(wtest)
  ztest=abs(wtest-ytest)
  leap=0
  if(ztest < 0.25) then
      leap=1
  end if
!
  julian=0
  if (mon >= 2) then
    do k=1,mon-1
     julian=julian+mond(k)
    enddo
  end if

  if (mon > 2) then
      julian=julian+leap
  end if
!
! now add the day of the month to the total
!
  julian=julian+day
!
! Determine UTC time of sunrise and sunset.  Start looping over times
! before sunrise, so use 5 UTC.  Change as needed.  Time of sunset may be
! after 24 UTC in the warm season. 
!  variable stime is time of sunset in hours (so 6.5 would be 6:30 am)
    
    sunrise=-10.
    sunset=10.
    
    delta=23.45*cos( (2.*pi*(1.*julian-173.)) / 365.25)
    
    do k=1,288
        xh=( ( 5. + (k-1)*5./60. - 12.0)*pi )/12. + xlon*pi/180.
        stime = 5.+(k-1)*5./60.
        rnside=sin(xlat*dtr)*sin(delta*dtr)+cos(xlat*dtr)*cos(delta*dtr)*cos(xh)
        qs=solar*rnside
        if (qs > 0.0 .and. sunrise < 0.0) then 
            sunrise=stime
            sunset=-10.
        end if
        if (qs < 0.0 .and. sunset < 0.0) then
            sunset=stime
        end if
    enddo
!
! verify that sunrise and sunset are working
!
  print *,' QVP code: Sunrise and Sunset in UTC for this radar are: '
  print *,sunrise,'   ',sunset
  print *,' QVP code: Example time format: '
  print *,time(itimes)

! Estimate CBL depth by searching for the minimum value in ZDR in
! vertical profile, but constrained to be near the ground at sunrise
! and moving upward or constant as time increases with a maximum
! vertical change in CBL height determined by time interval between
! radar scans. The array "depth" holds the CBL depth profile as function
! of time in units of array index value, not meters.
    
    kstart=10
     
    depth(1) = kstart
    depth(2) = kstart

!
! loop over observation times
!     note that lchg is the max change in height allowed from one observation time
!     to the next observation time, in array index values.  lchg is larger
!     after sunrise than during night
!
    kend=itimes-1
    do k=3,kend
         
!
! after sunrise then
!
    if (time(k) >= sunrise) then 
        lchg = nint(200*(time(k)-time(k-1)))
!
! define start (ist) and end (ien) indices for finding minimum value
!
        ist=depth(k-1)-lchg
        if (ist < 1) then 
            ist = 1
        end if
        ien=depth(k-1)+lchg

        zdrmax=10.
        depth(k) = depth(k-1)
        do l=ist,ien
            if (qvp_zdr(k,l) < zdrmax) then
                zdrmax = qvp_zdr(k,l)
                depth(k) = l
            end if
        enddo
!
! before sunrise then
!            
    else
       depth(k)=kstart
    end if

    if (depth(k) < kstart) then
        depth(k) = kstart
    end if

    enddo

!
! smooth CBL height evolution
! optional smoothing of the CBL height over time
!
    iend=itimes-5    
    do i = 5,iend
        slice(1)=depth(i-3)
        slice(2)=depth(i-2)
        slice(3)=depth(i-1)
        slice(4)=depth(i)
        slice(5)=depth(i+1)
        slice(6)=depth(i+2)
        slice(7)=depth(i+3)
        call mdian(slice,7,xmed)
        depth(i)=xmed
    enddo

end subroutine qvp
!***************************************************************
!*  Calculate the median value of an array with the Heapsort   *
!*  method                                                     *
!* ----------------------------------------------------------- *
!* REFERENCE:                                                  *
!*      "NUMERICAL RECIPES By W.H. Press, B.P. Flannery,       *
!*       S.A. Teukolsky and W.T. Vetterling, Cambridge         *
!*       University Press, 1986" [BIBLI 08].                   *
!*******************************************************
!* Given an array X of N numbers, returns their median *
!* value XMED. The array X is modified and returned    *
!* sorted in ascending order.                          *
!*******************************************************
SUBROUTINE MDIAN(X,N,XMED)
  real X(N)
  call hpsort(N,X)
  N2=N/2
  if (2*N2.eq.N) then
    XMED = 0.5*(X(N2)+X(N2+1))
  else
    XMED = X(N2+1)
  endif
  return
END

!*****************************************************
!*  Sorts an array RA of length N in ascending order *
!*                by the Heapsort method             *
!* ------------------------------------------------- *
!* INPUTS:                                           *
!*	    N	  size of table RA                   *
!*      RA	  table to be sorted                 *
!* OUTPUT:                                           *
!*	    RA    table sorted in ascending order    *
!*                                                   *
!* NOTE: The Heapsort method is a N Log2 N routine,  *
!*       and can be used for very large arrays.      *
!*****************************************************         
SUBROUTINE HPSORT(N,RA)
  real RA(N)
  L=N/2+1
  IR=N
  !The index L will be decremented from its initial value during the
  !"hiring" (heap creation) phase. Once it reaches 1, the index IR 
  !will be decremented from its initial value down to 1 during the
  !"retirement-and-promotion" (heap selection) phase.
10 continue
  if(L > 1)then
    L=L-1
    RRA=RA(L)
  else
    RRA=RA(IR)
    RA(IR)=RA(1)
    IR=IR-1
    if(IR.eq.1)then
      RA(1)=RRA
      return
    end if
  end if
  I=L
  J=L+L
20 if(J.le.IR)then
  if(J < IR)then
    if(RA(J) < RA(J+1))  J=J+1
  end if
  if(RRA < RA(J))then
    RA(I)=RA(J)
    I=J; J=J+J
  else
    J=IR+1
  end if
  goto 20
  end if
  RA(I)=RRA
  goto 10
END
