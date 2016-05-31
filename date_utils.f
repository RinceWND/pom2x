module Date_Utility
  implicit none
  private
  public :: Date_since, Days_in_between, get_Month, get_Month_Range, Number_of_Days, get_Day_of_Year, Days_to_TimeStamp
  public :: T_TimeStamp, T_Zone
  integer*1, parameter :: N_MONTHS = 13
  integer*2, parameter :: MONTH_DAYS(N_MONTHS) = (/ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365 /)
!  integer*2, parameter :: YEAR_OFFSET = 1900
  type T_Date
  sequence
    integer*2 :: year  = 2000
    integer*1 :: month = 1
    integer*1 :: day   = 1
  end type T_Date
  
  type T_TimeShort
  sequence
    integer*1 :: hour   = 0
    integer*1 :: minute = 0
  end type T_TimeShort
  
  type T_Time
  sequence
    integer*1        :: hour    = 0
    integer*1        :: minute  = 0
    double precision :: seconds = 0.d0
  end type T_Time
  
  type T_Zone
  sequence
    logical   :: negative = .false.
    integer*1 :: hour     = 0
    integer*1 :: minute   = 0
  end type T_Zone
  
  type T_DateTime
  sequence
    type (T_Date) :: date
    type (T_Time) :: time
  end type T_DateTime
  
  type T_DateTimeShort
  sequence
    type (T_Date)      :: date
    type (T_TimeShort) :: time
  end type T_DateTimeShort
  
  type T_TimeStampShort
  sequence
    type (T_Date)      :: date
    type (T_TimeShort) :: time
    type (T_Zone)      :: zone
  end type T_TimeStampShort
  
  type T_TimeStamp
  sequence
    type (T_Date) :: date
    type (T_Time) :: time
    type (T_Zone) :: zone
  end type T_TimeStamp
  
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
! Parse datetime string 'YYYY-MM-DD hh:ii:ss {+|-}HH:II'
!
    double precision function get_Day_of_Year(Time)
      
      double precision, intent(in) :: Time
      type(T_TimeStamp)            :: TimeStamp
      type(T_Zone)                 :: Zone
      double precision             :: TimeOfTheDay
      
      TimeStamp = Days_to_TimeStamp(Time, Zone)
      TimeOfTheDay = Time - floor(Time)
      
      get_Day_Of_Year = -1.
! Error checking
      if ( TimeStamp%date%year < 1 ) return
      if ( TimeStamp%date%month < 1 .or. TimeStamp%date%month > N_MONTHS-1 ) return
! Compute day of year
      get_Day_of_Year = MONTH_DAYS(TimeStamp%date%month) + TimeStamp%date%day
      if ( Is_Leap_Year(TimeStamp%date%year) .and. TimeStamp%date%month > 2) get_Day_of_Year = get_Day_of_Year + 1
      get_Day_of_Year = get_Day_of_Year + TimeOfTheDay
    
    end function get_Day_of_Year
!___________________________________________________________________
!
! Parse datetime string 'YYYY-MM-DD hh:ii:ss {+|-}HH:II'
!
    function TimeStamp_Parse(TimeStampString) result(TimeStamp)
        
      character(len=*), intent(in) :: TimeStampString
      type (T_TimeStamp) :: TimeStamp
      integer*1 :: offset
      integer   :: sec
      character :: sign
      
! Error checks
      if (len(TimeStampString)<10) then
        write(*,*) "[date_utility] /invalid timestamp format/"
        write(*,*) "    The TimeStamp string is too short."
        write(*,*) ""
        return
      end if
      
      read(TimeStampString, '(i4,1x,i2,1x,i2)') TimeStamp%date%year,  &
                                                TimeStamp%date%month, &
                                                TimeStamp%date%day
      offset = 18
      
      read(TimeStampString(offset-6:offset-2), '(i2,1x,i2)') TimeStamp%time%hour, TimeStamp%time%minute
! Detect a format possibility
      if (TimeStampString(offset-1:offset-1).eq.":") then
        read(TimeStampString(offset:offset+1), '(i2)') sec
        TimeStamp%time%seconds = real(sec, 8)
        offset = offset+3
      end if
      read(TimeStampString(offset:offset+5), '(a1,i2,1x,i2)') sign, TimeStamp%zone%hour, TimeStamp%zone%minute
      if (sign=='-') TimeStamp%zone%negative = .true.
        
    end function TimeStamp_Parse
!___________________________________________________________________
!
! Build datetime string (with +00:00 timezone)
!
    character(len=26) function TimeStamp_Build(TimeStamp, mode)
        
      type (T_TimeStamp), intent(in) :: TimeStamp
      character*1 :: mode
      character*1 :: sign
        
      sign = "+"
      if ( TimeStamp%zone%negative ) sign = "-"
      if ( mode=="s" ) then
        write(TimeStamp_Build,                                      &
         '(i4,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a,a1,i0.2,a,i0.2,a3)')    &
                            TimeStamp%date%year,    "-",            &
                            TimeStamp%date%month,   "-",            &
                            TimeStamp%date%day,     " ",            &
                            TimeStamp%time%hour,    ":",            &
                            TimeStamp%time%minute,  " ",            &
                            sign, TimeStamp%zone%hour, ":",         &
                                  TimeStamp%zone%minute, "   "
      else
        write(TimeStamp_Build,                                      &
         '(i4,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a,a1,i0.2,a,i0.2)')&
                            TimeStamp%date%year,    "-",            &
                            TimeStamp%date%month,   "-",            &
                            TimeStamp%date%day,     " ",            &
                            TimeStamp%time%hour,    ":",            &
                            TimeStamp%time%minute,  ":",            &
                            int(TimeStamp%time%seconds), " ",       &
                            sign, TimeStamp%zone%hour, ":",         &
                                  TimeStamp%zone%minute
      end if
        
    end function TimeStamp_Build
!___________________________________________________________________
!
! Build datetime string from days
!
    character(len=26) function TimeStamp(Days, TimeZone)
        
      double precision, intent(in) :: Days
      type (T_Zone) TimeZone
        
      TimeStamp = TimeStamp_Build(Days_to_TimeStamp(Days, TimeZone), "s")
        
    end function TimeStamp
!___________________________________________________________________
!
! Get number of days from TimeStamp
!
    double precision function days_number(TimeStamp)
    
      type (T_TimeStamp) :: TimeStamp
      
      days_number = int(TimeStamp%date%year*365.           &
                      + ceiling(TimeStamp%date%year/4.)    &
                      - ceiling(TimeStamp%date%year/100.)  &
                      + ceiling(TimeStamp%date%year/400.)  &
                      + MONTH_DAYS(TimeStamp%date%month) + TimeStamp%date%day)
                          
      if ( TimeStamp%date%month>2 .and. Is_Leap_Year(TimeStamp%date%year) ) days_number = days_number + 1
      
      days_number = days_number + time_to_days(TimeStamp%time)
      
    end function days_number
!___________________________________________________________________
!
!
!
    double precision function Number_of_Days(TimeStampString)
    
      character(len=*), intent(in) :: TimeStampString
      type (T_TimeStamp) :: TimeStamp
      double precision :: Days, Time
      
      TimeStamp = TimeStamp_Parse(TimeStampString)

      Number_of_Days = days_number(TimeStamp)
      
    end function Number_of_Days
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
    double precision function time_to_seconds(Time)
      type (T_Time), intent(in) :: Time
      
      time_to_seconds = Time%seconds + 60.*(Time%minute + Time%hour*60.)
      
    end function time_to_seconds
    
    double precision function time_to_days(Time)
      type (T_Time), intent(in) :: Time
      
      time_to_days = ((Time%seconds/60. + Time%minute)/60. + Time%hour)/24.
      
    end function time_to_days
!
    function Days_to_TimeStamp(Days, TimeZone) result(TimeStamp)
      double precision, intent(in) :: Days
      type (T_Zone), intent(in) :: TimeZone
      type (T_TimeStamp) :: TimeStamp
      double precision :: Time
      integer*1 :: leap
      
      TimeStamp%date%year = int(Days/365.24)  ! 365.24 - mean year length in days (leap years are taken into account)
      Time = floor(Days-(TimeStamp%date%year*365.   &
                  +ceiling(TimeStamp%date%year/4.)  &
                  -ceiling(TimeStamp%date%year/100.)&
                  +ceiling(TimeStamp%date%year/400.)))  ! Tmp storage
      ! correction variable to fix month detection for leap years
      leap = 0
      if (Is_Leap_Year(TimeStamp%date%year) .and. Time > 59.) leap = 1
      TimeStamp%date%month = minloc(MONTH_DAYS, 1, MONTH_DAYS >= Time-leap)-1
      
      TimeStamp%date%day = floor(Time - float( MONTH_DAYS(TimeStamp%date%month) ))
      if ( TimeStamp%date%month>2 .and. Is_Leap_Year(TimeStamp%date%year) ) TimeStamp%date%day = TimeStamp%date%day - 1
      
      TimeStamp%time%seconds = Days - floor(Days)
      TimeStamp%time%hour    = floor(TimeStamp%time%seconds*24.)
      TimeStamp%time%minute  = floor((TimeStamp%time%seconds*24.-TimeStamp%time%hour)*60.)
      TimeStamp%time%seconds = ((TimeStamp%time%seconds*24.-TimeStamp%time%hour)*60.-TimeStamp%time%minute)*60.
      
      TimeStamp%zone = TimeZone
      
    end function Days_to_TimeStamp
!
!
    function Days_since_to_TimeStamp(TimeStartString, Time) result(TimeStamp)
      character(len=*), intent(in) :: TimeStartString
      double precision, intent(in) :: Time
      double precision :: NewTime
      type (T_TimeStamp) :: TimeStamp
      
      TimeStamp = TimeStamp_Parse(TimeStartString)
      NewTime = days_number(TimeStamp)+Time
      TimeStamp = Days_to_TimeStamp(NewTime, TimeStamp%zone)
      
    end function Days_since_to_TimeStamp
!
! [ Exported function ]
! 
!  Returns timestamp string of date plus specified amount of days.
!
    character(len=26) function Date_since(TimeStartString, Time) result(TimeEndString)
      character(len=*), intent(in) :: TimeStartString
      double precision, intent(in) :: Time
      type (T_TimeStamp) :: TimeStamp
      
      TimeStamp = Days_since_to_TimeStamp(TimeStartString, Time)
      TimeEndString = TimeStamp_Build(TimeStamp, "l")
    end function Date_since
!    
    double precision function Days_in_Between(time_start, time_end)
      character(len=*), intent(in) :: time_start, time_end
      type (T_TimeStamp) :: TimeStart, TimeEnd
      
      TimeStart = TimeStamp_Parse(time_start)
      TimeEnd   = TimeStamp_Parse(time_end)
      
      Days_in_Between = Number_of_Days(time_end) - Number_of_Days(time_start)
      
    end function Days_in_Between
    
    integer function get_Month(Time)
      
      double precision, intent(in) :: Time
      type (T_TimeStamp)           :: TimeStamp
      type (T_Zone)                :: Zone
      
      TimeStamp = Days_to_TimeStamp(Time, Zone)
      
      get_Month = TimeStamp%date%month
      
    end function get_Month
!
!
!
    function get_Month_Range(Time) result(Range)
      
      double precision, intent(in) :: Time
      integer, dimension(2)        :: Range
      type (T_TimeStamp)           :: TimeStamp
      type (T_Zone)                :: Zone
      
      TimeStamp = Days_to_TimeStamp(Time, Zone)
      
      Range(1) = MONTH_DAYS(TimeStamp%date%month)
      Range(2) = MONTH_DAYS(TimeStamp%date%month+1)
      
      if (Is_Leap_Year(TimeStamp%date%year)) then
        if (TimeStamp%date%month >= 2) then
          Range(2) = Range(2)+1
          if (TimeStamp%date%month > 2) Range(1) = Range(1)+1
        end if
      end if
      
    end function
    
!    integer function Days_in_Between_I_S(time_start, time_end)
!      double precision, intent(in) :: time_start
!      character(len=*), intent(in) :: time_end
!      
!      Days_in_Between = abs(Number_of_Days(time_end) - time_start)
!      
!    end function Days_in_Between_I_S

end module Date_Utility