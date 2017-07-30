MODULE copyright_module
  CHARACTER (LEN=*), PARAMETER :: copyright_text = '(&
       &"     ___________________________________________________________________________"/&
       &"                                                                                "/&
       &"                 MX Toolbox (MXT)                                               "/&
       &"                 _______________________________________________________________"/&
       &"                 version 1.0 / January 2017                                     "/&
       &"                                                                                "/&
       &"                 Copyright (c) Rangana Warshamanage, PSI, Switzerland, 2017.    "/&
       &"                                                                                "/&
       &"                                                                                "/&
       &"                 This toolbox contains useful tools to analyze diffraction data "/&
       &"                 derived mainly from XDS and XSCALE. MXT is a byproduct of      "/&
       &"                 data processing and selection pipeline developed at SLS/PSI.   "/&
       &"     ___________________________________________________________________________"/&
       &) '
END MODULE copyright_module
MODULE menu_module
  CHARACTER (LEN=*), PARAMETER :: menu_text = '(&
       &"                                                        "/&
       &"     _____________________MAIN MENU_____________________"/&
       &"                                                        "/&
       &"     1 - Data set correlation coefficient calculation   "/&
       &"     2 - Data frame correlation coefficient calculation "/&
       &"     3 - I/Sig calculation                              "/&
       &"     4 - CC(1/2) calculation                            "/&
       &"     5 - Chi-squared test                               "/&
       &"     6 - Dataset remover                                "/&
       &"     7 - Exit                                           "/&
       &)'
END MODULE menu_module
MODULE parameter_module
  character*120 :: hklfile,arg_name
  integer       :: iarg,narg,choice
  !real          :: temp
  !logical       :: rescutoff
END MODULE parameter_module

MODULE global_module
  use copyright_module
  use menu_module
  use subroutines
  use parameter_module
END MODULE global_module


PROGRAM cc_datasets
  use global_module

  IMPLICIT NONE


  WRITE(*,copyright_text)
  !WRITE(*,menu_text)

  rescutoff = .false.

  narg = IARGC()
  IF(narg < 1)THEN
     WRITE(*,menu_text)
     PRINT*, 'usage: ./mxt xscalefile [-rc]'
     PRINT*, 'xscalefile = your corresponding XSCALE.ahkl'
     PRINT*, '-rc [optional] = include resolution cutoffs'
     PRINT*,
     PRINT*,
     STOP
  ELSE
     CALL GETARG(1,hklfile)
  END IF
  IF(narg > 1)THEN
     do iarg=2, narg
        CALL GETARG(iarg,arg_name)
        if(adjustl(arg_name)=='-rc')then
           rescutoff = .true.
        else
           print*,'Unknown argument'
           stop
        end if
     end do
  END IF

10 WRITE(*,menu_text)
  print*, 'Enter your choice:'
  read(*,*) choice
  select case(choice)
  case (1)
     print*, 'Correlation Coefficient data set calculation'
     if(rescutoff) call getcutoff(resol_ll,resol_ul)
     call cc_set(hklfile)
     go to 10
  case (2)
     print*, 'Correlation Coefficient frame calculation'
     if(rescutoff) call getcutoff(resol_ll,resol_ul)
     call cc_frames(hklfile)
     go to 10
  case (3)
     print*, 'Data set I/sigma calculation'
     if(rescutoff) call getcutoff(resol_ll,resol_ul)
     call xscale2stat(hklfile)
     go to 10
  case (4)
     print*, 'Half correlations calculation'
     if(rescutoff) call getcutoff(resol_ll,resol_ul)
     call cchalf_xscale(hklfile)
     go to 10
  case (5)
     print*, 'Chi-squared calculation'
     go to 10
  case (6)
     print*, 'Dataset remover (inputs: ccd_bin#.txt, xscale.inp)'
     call cc_set_remover
     go to 10
  case (7)
     print*, 'Exiting toolbox'
     go to 20
  case default
     go to 10
  end select
20 stop

CONTAINS

  subroutine getcutoff(lowres,highres)
    implicit none
    real          :: temp
    real,intent(out) :: lowres,highres
    
    write(*,*) 'Enter your resolution cutoffs: low high'
    read(*,*) lowres, highres
    if(lowres < highres)then
       temp = lowres
       lowres = highres
       highres = temp
    end if
    print*, 'Resolution: low,high:',lowres,' ',highres
    !stop
  end subroutine getcutoff
  

END PROGRAM cc_datasets
