!     
! File:   ini2k.f
! Author: pkharlamov
!
! Created on 13 Май 2015 г., 13:18
!
module ini_get

    implicit none

    private
    !public :: readIni

    !character*256 filepath

    contains
    
      subroutine readIni
        implicit none
        print *, "read...inrffn"
        call readLine("")
        call readLine(" whining ")
      end subroutine readIni

      subroutine readLine(pesky)
        implicit none
        character(len=*) pesky
        print *, "[ lining ]"//pesky
      end subroutine readLine

end module ini_get