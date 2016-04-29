module Date_Utility
  implicit none
  private
  public :: Number_of_Days, Days_since_to_Date, Days_since_to_DateTime, Days_to_TimeStamp, &
            Date_Parse, Days_to_Date, Days_to_DateTime, Number_of_Days_with_Time, Days_in_Between, &
            DateTime_Build
  integer*1, parameter :: N_MONTHS = 12
  integer*2, parameter :: MONTH_DAYS(N_MONTHS) = (/ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 /)
  
  contains
!___________________________________________________________________
!
! Check for leap year
!
    logical function Is_Leap_Year(Year)
    
      integer*2, intent(in) :: Year
      Is_Leap_Year = .false.
      if ( ( mod( Year, 4 ) == 0 .and. mod( Year, 100 ) /= 0 ) .or. mod( Year, 400 ) == 0 ) Is_Leap_Year = .true.
      
    end function Is_Leap_Year
!___________________________________________________________________
!
! Get the number of a day in a year
!
    function Date_to_Day_of_Year( Day_of_Month, Month, Year ) result(Day_Of_Year)
      
      integer*1, intent(in) :: Day_of_Month
      integer*1, intent(in) :: Month
      integer*2, intent(in) :: Year
      integer*2 :: Day_of_Year
      integer*2 :: Days_of_Month(N_MONTHS)

      Day_Of_Year = -1
! Compute days per month
      Days_of_Month = MONTH_DAYS
! Error checking
      if ( Year < 1 ) return
      if ( Month < 1 .or. Month > N_MONTHS ) return
      if ( Day_of_Month > Days_of_Month(Month) ) return
! Compute day of year
      Day_of_Year = Days_of_Month(Month) + Day_of_Month
      if ( Is_Leap_Year(Year) ) Day_of_Year = Day_of_Year + 1
      
    end function Date_to_Day_of_Year
!___________________________________________________________________
!
! Parse date string 'YYYY-MM-DD'
!
    function Date_Parse(time_str) result(Date)
        
      character(len=*), intent(in) :: time_str
      integer*2, dimension(3) :: Date
        
      read(time_str, '(i4,1x,i2,1x,i2)') Date
        
    end function Date_Parse
!___________________________________________________________________
!
! Parse datetime string 'YYYY-MM-DD hh:ii:ss'
!
    function DateTime_Parse(time_str) result(DateTime)
        
      character(len=*), intent(in) :: time_str
      integer*2, dimension(6) :: DateTime
        
      read(time_str, '(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)') DateTime
        
    end function DateTime_Parse
!___________________________________________________________________
!
! Parse datetime string 'YYYY-MM-DD hh:ii:ss <+/->HH:II'
!
    function TimeStamp_Parse(time_str) result(DateTime)
        
      character(len=*), intent(in) :: time_str
      integer*2, dimension(8) :: DateTime
      character :: sign
        
      read(time_str, '(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2,1x,a1,i2,1x,i2)') DateTime(1:6), sign, DateTime(7:8)
      if (sign=='-') then
        if (DateTime(7)/=0) then
          DateTime(7) = -DateTime(7)
        else
          DateTime(8) = -DateTime(8)
        end if
      end if
        
    end function TimeStamp_Parse
!___________________________________________________________________
!
! Build datetime string (with +00:00 timezone)
!
    character(len=26) function TimeStamp_Build(DateTime)
        
      integer*2, dimension(8) :: DateTime
      character :: sign
        
      sign = '+'
      if (DateTime(7)<0.or.DateTime(8)<0) sign = '-'
      write(DateTime_Build,                                         &
          '(i4,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a,a,i0.2,a,i0.2)')&
                            DateTime(1), "-",                       &
                            DateTime(2), "-",                       &
                            DateTime(3), " ",                       &
                            DateTime(4), ":",                       &
                            DateTime(5), ":",                       &
                            DateTime(6), " ", sign,                 &
                            DateTime(7), ":", DateTime(8)
        
    end function TimeStamp_Build
!___________________________________________________________________
!
! Build datetime string from days (with +00:00 timezone)
!
    character(len=26) function Days_to_TimeStamp(Days)
        
      double precision, intent(in) :: Days
        
      Days_to_TimeStamp = TimeStamp_Build(Days_to_TimeStamp(Days))
        
    end function Days_to_TimeStamp
!___________________________________________________________________
!    
    integer function Number_of_Days(time_start)
    
      character(len=*), intent(in) :: time_start
      integer*2 :: Year
      integer*1 :: Month, Day
      integer*2 :: Days_of_Month(N_MONTHS)
      integer*2 :: Date(3)
      
      Days_of_Month = MONTH_DAYS
      
      Date = Date_Parse(time_start)
      Year  = Date(1)
      Month = Date(2)
      Day   = Date(3)
      
      Number_of_Days = int(Year*365.            &
                          + ceiling(Year/4.)    &
                          - ceiling(Year/100.)  &
                          + ceiling(Year/400.)  &
                          + Days_of_Month(Month) + Day)
                          
      if ( Month>2 .and. Is_Leap_Year(Year) ) Number_of_Days = Number_of_Days + 1
      
    end function Number_of_Days
!
    double precision function Number_of_Days_with_Time(time_start)
    
      character(len=*), intent(in) :: time_start
      integer*1 :: Hour, Minute, Seconds
      integer   :: Days
      integer*2 :: DateTime(6)
      double precision :: Time
      
      Days = Number_of_Days(time_start)
      
      DateTime = DateTime_Parse(time_start)
      Hour    = DateTime(4)
      Minute  = DateTime(5)
      Seconds = DateTime(6)
      
      Time = (Hour + Minute/60. + Seconds/3600.)/24.
      
      Number_of_Days_with_Time = real(Days, 8) + Time
      
    end function Number_of_Days_with_Time
!
!    
    function Days_to_Date(Days) result(Date)
      
      double precision, intent(in) :: Days
      integer*2 :: Date(3)
      integer*2 :: Days_of_Month(N_MONTHS)
      integer*1 :: leap
      
      Days_of_Month = MONTH_DAYS
      
      Date(1) = int(Days/365.24)  ! 365.24 - mean year length in days (leap years are taken into account)
      Date(3) = Days-(Date(1)*365.+ceiling(Date(1)/4.)-ceiling(Date(1)/100.)+ceiling(Date(1)/400.))  ! Tmp storage
      ! correction variable to fix month detection for leap years
      leap = 0
      if (Is_Leap_Year(Date(1))) leap = 1
      Date(2) = minloc(Days_of_Month, 1, Days_of_Month >= Date(3)-leap)-1
      
      Date(3) = floor(float( Date(3) - Days_of_Month(Date(2)) ))
      if ( Date(2)>2 .and. Is_Leap_Year(Date(1)) ) Date(3) = Date(3) - 1
      
    end function Days_to_Date
!
!
    function Days_to_DateTime(Days) result(DateTime)
      double precision, intent(in) :: Days
      double precision :: Time
      integer*2 :: DateTime(6)
      
      DateTime(1:3) = Days_to_Date(Days)
      
      Time = Days - floor(Days)
      
      DateTime(4) = floor(Time*24.)
      DateTime(5) = floor((Time*24.-DateTime(4))*60.)
      DateTime(6) = ((Time*24.-DateTime(4))*60.-DateTime(5))*60.+.5 ! 0.5 is to compensate for integer rounding. It makes 60 seconds possible, though.
      
    end function Days_to_DateTime
!
!
    function Days_to_TimeStamp(Days, TZHour, TZMinute) result(TimeStamp)
      double precision, intent(in) :: Days
      integer*1,        intent(in) :: TZHour, TZMinute
      double precision :: Time
      integer*2 :: TimeStamp(8)
      
      if ((TZHour<0.or.TZHour>12).or.(TZMinute<0.or.TZMinute>60)) then
        write(*,*) "[date_utility] /invalid time/"
        write(*,*) "    Error passing timezone parameters to function Days_to_TimeStamp."
        write(*,*) ""
      end if
      DateTime(1:6) = Days_to_DateTime(Days)
      
      Time = Days - floor(Days)
      
      DateTime(4) = floor(Time*24.)
      DateTime(5) = floor((Time*24.-DateTime(4))*60.)
      DateTime(6) = ((Time*24.-DateTime(4))*60.-DateTime(5))*60.+.5 ! 0.5 is to compensate for integer rounding. It makes 60 seconds possible, though.
      
    end function Days_to_TimeStamp
!
!
    function Days_since_to_Date(time_start, Days) result(Date)
      character(len=*), intent(in) :: time_start
      double precision, intent(in) :: Days
      double precision :: Days_of_Start
      integer*2 :: Date(3)
      
      Days_of_Start = Number_of_Days(time_start) + Days
      Date = Days_to_Date(Days_of_Start)      
      
    end function Days_since_to_Date
!
    function Days_since_to_DateTime(time_start, Time) result(DateTime)
      character(len=*), intent(in) :: time_start
      double precision, intent(in) :: Time
      double precision :: Time_of_Start
      integer*2 :: DateTime(6)
      
      Time_of_Start = Number_of_Days_with_Time(time_start) + Time
      DateTime = Days_to_DateTime(Time_of_Start)
      
    end function Days_since_to_DateTime
    
    integer function Days_in_Between(time_start, time_end)
      character(len=*), intent(in) :: time_start, time_end
      
      Days_in_Between = abs(Number_of_Days(time_end) - Number_of_Days(time_start))
      
    end function Days_in_Between
    
!    integer function Days_in_Between_I_S(time_start, time_end)
!      double precision, intent(in) :: time_start
!      character(len=*), intent(in) :: time_end
!      
!      Days_in_Between = abs(Number_of_Days(time_end) - time_start)
!      
!    end function Days_in_Between_I_S

end module Date_Utility