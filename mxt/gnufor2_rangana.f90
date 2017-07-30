!I'm grateful to the author of this kit. Unfortunately I cannot remember where I 
!downloaded this from. I have changed the code as I needed.

!***********************************************************************************
module gnufor2_rangana
  implicit none

  public :: plot_1
  
  !***********************************************************************************
  ! these are default parameters which control linewidth, colors and terminal
  !***********************************************************************************
  character(len=3), parameter	:: default_linewidth='1'
  character(len=100), parameter	:: default_color1='blue'
  character(len=100), parameter	:: default_color2='dark-green'
  character(len=100), parameter	:: default_color3='orange-red'
  character(len=100), parameter	:: default_color4='dark-salmon'
  character(len=100), parameter	:: default_terminal='X11'
  character(len=100), parameter	:: default_palette='CMY'
  !***********************************************************************************

  !***********************************************************************************
  !***********************************************************************************
  !***********************************************************************************
contains
  !***********************************************************************************
  !***********************************************************************************
  !***********************************************************************************
  function my_date_and_time() result(f_result)
    !***********************************************************************************
    ! this function creates a string with current date and time
    ! it is a default method to name output files
    !***********************************************************************************
    implicit none
    character(len=8)  :: date
    character(len=10) :: time
    character(len=33) :: f_result
    !***********************************************************************************
    call date_and_time(date,time)
    f_result= 'date_'//date(7:8)//'-'//date(5:6)//'-'//date(1:4)//'_time_'//time(1:2)//':'//time(3:4)//':'//time(5:10)
    !***********************************************************************************
  end function my_date_and_time
  !***********************************************************************************
  !***********************************************************************************
  !***********************************************************************************
  function output_terminal(terminal) result(f_result)
    !***********************************************************************************
    implicit none
    character(len=*),intent(in)	:: terminal
    integer, parameter		:: Nc=35
    character(len=Nc)		:: f_result
    !***********************************************************************************
    select case(terminal)
    case('ps')
       f_result='postscript landscape color'
    case default
       f_result=terminal
    end select
    !***********************************************************************************
  end function output_terminal
  !***********************************************************************************

  !***********************************************************************************
  subroutine plot_1(x1,y1,style,pause,color1,terminal,filename,polar,persist,input,linewidth)
    !***********************************************************************************
    ! this subroutine plots a two-dimensional graph
    !***********************************************************************************
    implicit none
    integer, intent(in)	        :: x1(:)
    real(kind=4), intent(in)	::y1(:)
    real(kind=4), optional	:: pause,linewidth
    character(len=*),optional	:: style, color1, terminal, filename, polar, persist, input
    integer 			:: i, ierror, ios, file_unit, Nx1
    character(len=100)		:: data_file_name, command_file_name, my_linewidth
    integer, parameter		:: Nc=20
    character(len=Nc)		:: my_line_type1, my_color1, my_range, my_pause, my_persist
    !***********************************************************************************
    if (present(input)) then
       data_file_name='data_file_'//input//'.txt'
       command_file_name='command_file_'//input//'.txt'		
    else
       data_file_name='data_file.txt'
       command_file_name='command_file.txt'
    end if
    !***********************************************************************************
    Nx1=size(x1)
    if ((size(x1).ne.size(y1))) then
       print *,'subroutine plot ERROR: size(x) is not equal to size(y)'
       stop
    end if
    if (present(style).and.(len(style).ne.3)) then
       print *,'subroutine plot ERROR: argument "style" has wrong number of characters'
       stop
    end if
    !***********************************************************************************
    ierror=0	
    call get_unit(file_unit)	
    if (file_unit==0) then
       ierror=1
       print *,'write_vector_data - fatal error! Could not get a free FORTRAN unit.'
       stop
    end if
    open (unit=file_unit, file=data_file_name, status='replace', iostat=ios)	
    if (ios/=0) then
       ierror=2
       print *,'write_vector_data - fatal error! Could not open the terminal data file.'
       stop
    end if
    !***********************************************************************************	
    ! here we write the date to the data_file - the gnuplot will read this data later
    !***********************************************************************************	
    do i=1,Nx1
       write (file_unit,'(I5,F8.3)') x1(i), y1(i)
    end do
    !***********************************************************************************	
    close (unit=file_unit)
    !***********************************************************************************
    ierror = 0
    call get_unit(file_unit)
    if (file_unit==0) then
       ierror=1
       print *,'write_vector_data - fatal error! Could not get a free FORTRAN unit.'
       stop
    end if
    open (unit=file_unit, file=command_file_name, status='replace',	iostat=ios)
    if (ios/=0) then
       ierror=2
       print *,'write_vector_data - fatal error! Could not open the terminal command file.'
       stop
    end if
    !***********************************************************************************
    ! here we write the commands to the commands file which gnuplot will execute
    !***********************************************************************************
    my_line_type1='lines'
    if (present(style)) then
       if ((style(3:3)=='-')) then
          my_line_type1='linespoints'
       else
          my_line_type1='points'
       end if
    end if
    if (present(linewidth)) then
       write (	my_linewidth,'(e9.3)') linewidth
    else
       my_linewidth=trim(default_linewidth)
    end if
    if (present(color1)) then
       my_color1='"'//trim(color1)//'"'
    else
       my_color1='"'//trim(default_color1)//'"'
    end if
    !***********************************************************************************
    my_persist='persist '
    if (present(persist).and.(persist=='no')) my_persist=' '
    if (present(terminal)) then
       write ( file_unit, '(a)' ) 'set terminal '// trim(output_terminal(terminal)) 
       if (present(filename)) then
          write ( file_unit, '(a)' ) 'set output "'// trim(filename) //'"' 
       else
          write ( file_unit, '(a)' ) 'set output "'//my_date_and_time()//'"' 
       end if
    else
       write ( file_unit, '(a)' ) 'set terminal ' // trim(default_terminal) // ' ' &
            & //trim(my_persist) //' title  "Gnuplot"' 
    end if
    !***********************************************************************************
    write ( file_unit, '(a)' ) 'unset key'
    if (present(polar).and.(polar=='yes')) then
       write (my_range,'(f8.3)') maxval(abs(y1))
       write ( file_unit, '(a)' ) 'set xrange [-'//trim(my_range)//':'//trim(my_range)//']'
       write ( file_unit, '(a)' ) 'set yrange [-'//trim(my_range)//':'//trim(my_range)//']'
       write ( file_unit, '(a)' ) 'set size square'
       write ( file_unit, '(a)' ) 'set polar'
       write ( file_unit, '(a)' ) 'set grid polar'
    else
       write ( file_unit, '(a)' ) 'set grid'
    end if
    !***********************************************************************************
    if (present(style)) then
       write ( file_unit, '(a,i2,a)' ) 'plot "' // trim (data_file_name) &
            &//'" using 1:2 with ' // trim(my_line_type1) // ' pointtype ' // &
            & style(1:2) // ' linecolor rgb ' // trim(my_color1) // ' linewidth '// trim(my_linewidth) 
    else 
       write ( file_unit, '(a,i2,a)' ) 'plot "' // trim (data_file_name) &
            & //'" using 1:2 with ' // trim(my_line_type1)  // ' linecolor rgb '&
            & // trim(my_color1) // ' linewidth '// trim(my_linewidth)  
    end if
    !***********************************************************************************
    if (present(pause)) then
       if (pause<0.0) then
          write ( file_unit, '(a)' ) 'pause -1 "press RETURN to continue"'
       else 
          write (	my_pause,'(e9.3)') pause
          write ( file_unit, '(a)' ) 'pause ' // trim(my_pause) 
       end if
    else
       write ( file_unit, '(a)' ) 'pause 0'
    end if
    !***********************************************************************************
    write ( file_unit, '(a)' ) 'q'
    close ( unit = file_unit )
    !***********************************************************************************
    call run_gnuplot (command_file_name) 
    !***********************************************************************************
  end subroutine plot_1
  !***********************************************************************************
  !***********************************************************************************

  !***********************************************************************************
  subroutine run_gnuplot(command_file_name)
    !***********************************************************************************
    implicit none
    character (len = 100) command
    character (len = *) command_file_name
    integer status
    integer system
    !***********************************************************************************
    !  Issue a command to the system that will startup GNUPLOT, using
    !  the file we just wrote as input.
    !***********************************************************************************
    write (command, *) 'gnuplot ' // trim (command_file_name)		
    status=system(trim(command))	
    if (status.ne.0) then
       print *,'RUN_GNUPLOT - Fatal error!'
       stop
    end if
    return
    !***********************************************************************************
  end subroutine run_gnuplot
  !***********************************************************************************
  !***********************************************************************************
  !***********************************************************************************
  subroutine get_unit(iunit)
    !***********************************************************************************
    implicit none
    integer i
    integer ios
    integer iunit
    logical lopen
    !***********************************************************************************	
    iunit=0
    do i=1,99
       if (i/= 5 .and. i/=6) then	
          inquire (unit=i, opened=lopen, iostat=ios)
          if (ios==0) then
             if (.not.lopen) then
                iunit=i
                return
             end if
          end if

       end if
    end do
    return
  end subroutine get_unit
  !***********************************************************************************
end module gnufor2_rangana
