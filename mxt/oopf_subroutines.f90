MODULE subroutines

  USE symops_mod
  USE resolution_mod
  USE gnufor2_rangana

  IMPLICIT NONE

  REAL, PUBLIC :: resol_ll,resol_ul
  LOGICAL, PUBLIC :: rescutoff

  TYPE, PRIVATE :: dataset
     integer iiset
     TYPE(dataset),POINTER :: nextset
  END TYPE dataset

  TYPE, PRIVATE :: dataframe
     integer iframe
     TYPE(dataframe),POINTER :: nextframe
  END TYPE dataframe

  TYPE, PRIVATE :: dataint
     real iintensity,sig
     TYPE(dataint),POINTER :: nextint
  END TYPE dataint

  TYPE, PRIVATE ::uniq
     integer uindex,ufrq
     real    hkl_resol
     type(dataset) sn
     type(dataint) inte
     type(dataframe) fn
     type(uniq), POINTER :: nextuniq
  END TYPE uniq

  type(uniq),      pointer :: ob,firstuniq
  type(dataset),   pointer :: sn, firstsn
  type(dataframe), pointer :: fn, firstfn
  type(dataint),   pointer :: inte, firstint

  character(len=80)   :: line
  character(len=120)  :: filename_template
  character*2         :: strng
  integer             :: unit,ier,ncheck,ispgr,iset,nset,stat_ar
  integer             :: ih,ik,il,jh,jk,jl,ihkl,oldhkl,iop
  real                :: fsq,sigma,dum4(4)
  real                :: start,finish,s2,resol,lowest_resol,highest_resol
  logical             :: last,xscale

  integer, dimension(:), allocatable :: uindex, setnum,frmnum
  real,    dimension(:), allocatable :: intensity,sig_ar,hkl_resol
  real,    dimension(:), allocatable :: x,y


CONTAINS

  subroutine cc_set(xscalefile)

1000 FORMAT(A)
1002 FORMAT(3I6,2G11.3,4F8.1,I3) !xscale.ahkl

    character (len=120),intent(in) :: xscalefile
    character*255 :: ccdfile
    integer       :: ni,nj,nbin,nobs,nuniq
    integer       :: i,j

    real,    dimension(:),   allocatable :: bin_arr
    real,    dimension(:,:), allocatable :: sumI
    integer, dimension(:,:), allocatable :: nfreq

    nuniq=0   ! Number of unique reflections
    nobs=0    ! Number of observations
    last=.false.
    ihkl=0;i = 0

    if(rescutoff) print*, 'Resolution cutoffs have enabled!'

    do i=1,len_trim(xscalefile)
       if(xscalefile(i:i)=='.') exit
    end do
    read(xscalefile(1:i-1),1000) filename_template

    call xscaleopen(xscalefile)

    allocate(ob); nullify(ob%nextuniq); firstuniq => ob
    allocate(sn);nullify(sn%nextset); firstsn => sn
    allocate(inte); nullify(inte%nextint);firstint => inte

    i = 0
    do
       read(unit,1000) line
       if(line(1:12) == '!END_OF_DATA') last=.true.
       if(.not.last)then
          if(xscale)then
             read(line,*) ih,ik,il,fsq,sigma,dum4,iset
          else
             print*,'This file is not a xscale.ahkl type file.'
             stop
          end if
          !if(sigma<0.) cycle
          if((fsq/sigma)<0.1)then
             write(30,1002) ih,ik,il,fsq,sigma,dum4,iset
             cycle
          end if
          call getres(ih,ik,il,s2,resol)

          if(rescutoff)then
             if(resol>resol_ll.or.resol<resol_ul) cycle !this removes not matching reflections
          end if

          nobs=nobs+1
          call asuput(ih,ik,il,jh,jk,jl,ihkl,iop,op,size(op,3)/ncenop)
          write(31,*) ih,ik,il,ihkl,resol
          if(nobs==1) oldhkl=ihkl
          i = i + 1
          if(i==1)then
             lowest_resol = resol
             highest_resol = resol
          else
             if(lowest_resol < resol) lowest_resol = resol
             if(highest_resol > resol) highest_resol = resol
          end if
          if(i==1) nuniq=1
          !if(i==1) write(31,*) ih,ik,il,fsq,sigma,fsq/sigma,iset
          if(ihkl/=oldhkl.or.last) nuniq = nuniq + 1
          !if(ihkl/=oldhkl.or.last) write(31,*) ih,ik,il,fsq,sigma,fsq/sigma,iset
          ob%uindex    = nuniq
          ob%sn%iiset    = iset
          ob%inte%iintensity = fsq
          ob%inte%sig = sigma
          ob%hkl_resol = resol
          allocate(ob%nextuniq); ob => ob%nextuniq;nullify(ob%nextuniq)
          allocate(sn%nextset); sn => sn%nextset;nullify(sn%nextset)
          allocate(inte%nextint); inte => inte%nextint;nullify(inte%nextint)
       end if
       if(last) exit
       oldhkl=ihkl
       !print*, ih,ik,il,resol,lowest_resol,highest_resol
    end do
    xscale = .false.
    print*, 'Reading complete'
    print*,'# obs (excluding misfits), # unique =',nobs,nuniq
    print*, 'lowest resolution=', lowest_resol
    print*, 'highest resolution', highest_resol
    print*,
    print*, 'Enter the number of resolution bins you want:'
    read(*,*) nbin
    call setbin(nbin,lowest_resol,highest_resol,2,bin_arr)
!!$    do i=1,nbin
!!$       print*, bin_arr(i)
!!$    end do
    !print*,'Allocating arrays for data...'
    if(nobs>0)then
       allocate(uindex(nobs),stat=stat_ar)
       if(stat_ar<0)print*,'Allocation for uindex - problem'
       allocate(setnum(nobs),stat=stat_ar)
       if(stat_ar<0)print*,'Allocation for setnum - problem'
       allocate(intensity(nobs),stat=stat_ar)
       if(stat_ar<0)print*,'Allocation for intensity - problem'
       allocate(hkl_resol(nobs),stat=stat_ar)
       if(stat_ar<0)print*,'Allocation for resolution - problem'
       allocate(sig_ar(nobs),stat=stat_ar)
       if(stat_ar<0)print*,'Allocation for sigmas - problem'
       allocate(sumI(nuniq,nset),stat=stat_ar)
       if(stat_ar<0)print*,'Allocation for sumI - problem'
       allocate(nfreq(nuniq,nset),stat=stat_ar)
       if(stat_ar<0)print*,'Allocation for nfreq - problem'
    end if
    !print*,'Allocation...done'

    !Initializing arrays
    do i = 1, nobs
       uindex(i) = 0
       setnum(i) = 0
       intensity(i) = 0.
       sig_ar(i) = 0.
       hkl_resol(i) = 0.
    end do
    do j = 1, nuniq
       do i = 1, nset
          sumI(j,i) = 0.
          nfreq(j,i) = 0
       end do
    end do

    !Transfering data to arrays
    i = 0
    ob => firstuniq
    sn => firstsn
    inte => firstint
    do while(associated(ob%nextuniq))
       i = i+1
       uindex(i) = ob%uindex
       setnum(i) = ob%sn%iiset
       hkl_resol(i) = ob%hkl_resol
       intensity(i) = ob%inte%iintensity
       sig_ar(i) = ob%inte%sig
       ob => ob%nextuniq
       sn => sn%nextset
       inte => inte%nextint
       !print*, uindex(i),setnum(i),hkl_resol(i),intensity(i)
    end do

    do i=1,nbin
       print*, bin_arr(i)
    end do
    do i=1,size(bin_arr)
       print*, i,bin_arr(i)
       write(strng,'(I2)') i
       ccdfile = trim(adjustl(filename_template))//'ccd_bin'//trim(adjustl(strng))//'.txt'
       do j=1,nobs
          if(i==1)then
             if(hkl_resol(j) >= bin_arr(i))then
                write(40,*) uindex(j),setnum(j),intensity(j),sig_ar(j),intensity(j)/sig_ar(j),hkl_resol(j) 
                sumI(uindex(j),setnum(j)) = sumI(uindex(j),setnum(j)) + intensity(j)
                nfreq(uindex(j),setnum(j)) = nfreq(uindex(j),setnum(j)) + 1
                !print*, sumI(uindex(j),setnum(j)),nfreq(uindex(j),setnum(j))
                write(45,*) uindex(j),setnum(j),sumI(uindex(j),setnum(j)),nfreq(uindex(j),setnum(j))
             endif
          elseif(i > 1)then
             if(hkl_resol(j) >= bin_arr(i) .and. hkl_resol(j) < bin_arr(i-1))then
                sumI(uindex(j),setnum(j)) = sumI(uindex(j),setnum(j)) + intensity(j)
                nfreq(uindex(j),setnum(j)) = nfreq(uindex(j),setnum(j)) + 1
             endif
          endif
       end do
       print*, 'calculating CC for ', bin_arr(i),i
       call cc_calc(nset,nuniq,sumI,nfreq,bin_arr(i),ccdfile)
       do nj=1,nuniq
          do ni=1,nset
             sumI(nj,ni) = 0.
             nfreq(nj,ni) = 0
          end do
       end do
    end do
    deallocate(uindex)!;print*,'uindex - deallocated'
    deallocate(setnum)!;print*,'setnum - deallocated'
    deallocate(intensity)!;print*,'intensity - deallocated'
    deallocate(hkl_resol)!;print*,'resolution - deallocated'
    deallocate(sig_ar)!;print*,'sigmas - deallocated'
    deallocate(sumI)!;print*,'sumI - deallocated'
    deallocate(nfreq)!;print*,'nfreq - deallocated'
    close(unit);close(40);close(30)
  end subroutine cc_set



  subroutine cc_frames(xscalefile)

    integer,parameter :: maxref=50000
    character(len=120),intent(in) :: xscalefile
    character(len=255) :: nframeset,ccframe,ccframeavg
    logical  :: present
    integer  :: i,j,ii,jj,is,k,l,nx,ny,nuniq,nobs,npairs
    integer  :: oldset,max_ii,olduindx,oldframe,frm,nfrms,ncc
    real     :: sumCC
    real     :: sumxx,sumyy,sumX,sumY,sumXY,sumX2,sumY2,AA,BB,CC,DD,CC2
    integer  :: ifq(maxref)
    integer, dimension(:),   allocatable :: set
    integer, dimension(:,:), allocatable :: dset,dframe,np
    real,    dimension(:,:), allocatable :: dint,CorC,SE

1000 FORMAT(A)

    npairs = 0;nx = 0;ny = 0;nuniq = 0;nobs=0
    sumX = 0.;sumY = 0.;sumXY = 0.;sumX2 = 0.;sumY2 = 0.
    AA = 0.;BB = 0.;CC = 0.;DD = 0.
    ihkl=0;i = 0;ii = 1;jj = 1
    max_ii = 1 ! This variable carrys the number of times a given uniq reflection be found in the dataset.
    last=.false.;present = .false.

    if(rescutoff) print*, 'Resolution cutoffs have enabled!'

    do i=1,len_trim(xscalefile)
       if(xscalefile(i:i)=='.') exit
    end do
    read(xscalefile(1:i-1),1000) filename_template

    !Declaration of output files
    nframeset = trim(adjustl(filename_template))//'_setandframes.txt'
    ccframe   = trim(adjustl(filename_template))//'_cc_frames.txt'
    ccframeavg = trim(adjustl(filename_template))//'_cc_averaged.txt'
    open(24,file=trim(adjustl(nframeset)),status='replace')
    open(28,file=trim(adjustl(ccframe)),status='replace')
    open(29,file=trim(adjustl(ccframeavg)),status='replace')
    write(24,1000)'set_number    frame_number'
    write(28,1000)'set#  frame#  cc_frame  err(cc)  No.of_pairs'
    write(29,1000)'frame#  avg_cc_frame'

    call xscaleopen(xscalefile)
    allocate(set(nset))
    allocate(ob); nullify(ob%nextuniq); firstuniq => ob
    allocate(sn);nullify(sn%nextset); firstsn => sn
    allocate(fn); nullify(fn%nextframe);firstfn => fn
    allocate(inte); nullify(inte%nextint);firstint => inte

    do
       read(unit,1000) line
       if(line(1:12) == '!END_OF_DATA') last=.true.
       if(.not.last)then
          if(xscale)then
             read(line,*)ih,ik,il,fsq,sigma,dum4(1),dum4(2),dum4(3),dum4(4),iset
          else
             !read(line,*)ih,ik,il,fsq,sig,dum4,peak,corr,psi
             print*,'This file is not a xscale.ahkl type file.'
             stop
          end if
          if((fsq/sigma)<0.) cycle
          If(abs(dum4(3)-nint(dum4(3)))<0.05) cycle
          frm = INT(dum4(3))+1
          call getres(ih,ik,il,s2,resol)

          if(rescutoff)then
             if(resol>resol_ll.or.resol<resol_ul) cycle !this removes not matching reflections
          end if

          nobs=nobs+1
          call asuput(ih,ik,il,jh,jk,jl,ihkl,iop,op,size(op,3)/ncenop)
          if(nobs==1)then
             oldhkl=ihkl
             nuniq = 1
             oldset = iset
             lowest_resol = resol
             highest_resol = resol
          else
             if(lowest_resol < resol) lowest_resol = resol
             if(highest_resol > resol) highest_resol = resol
          end if
          if(ihkl == oldhkl)then
             if(nobs == 1) ii = 0
             i = nuniq
             ii = ii + jj
             if(ii > max_ii) max_ii = ii
             ifq(i) = ii ! Choice for better performance
             ob%uindex = nuniq
             ob%sn%iiset = iset
             ob%fn%iframe = frm
             ob%inte%iintensity = fsq
             !WRITE(22,*) i,ii,ob%uindex,ob%setnum%iiset,ob%framenum%iframe,ob%intensity%iintensity !For testing
          else ! IF (ihkl /= oldhkl) THEN
             nuniq = nuniq+1
             j = nuniq
             jj = 1
             ifq(j) = jj ! Choice for better performance         
             ob%uindex = nuniq
             ob%sn%iiset = iset
             ob%fn%iframe = frm
             ob%inte%iintensity = fsq
             !WRITE(22,*) j,jj,ob%uindex,ob%setnum%iiset,ob%framenum%iframe,ob%intensity%iintensity !For testing
             ii = 1
          end if
          allocate(ob%nextuniq);  ob => ob%nextuniq;  nullify(ob%nextuniq)
          allocate(sn%nextset);   sn => sn%nextset;   nullify(sn%nextset)
          allocate(fn%nextframe); fn => fn%nextframe; nullify(fn%nextframe)
          allocate(inte%nextint);  inte => inte%nextint; nullify(inte%nextint)
       end if
       if(last) exit
       oldhkl=ihkl
       oldset = iset
    end do
    xscale = .false.
    print*, 'Reading complete.'
    !print*,'max_ii=',max_ii !For testing
    print*,'# obs (excluding misfits), # unique =',nobs,nuniq
    print*, 'lowest resolution=', lowest_resol
    print*, 'highest resolution', highest_resol

    !print*,'Allocating arrays for data...'
    if(nobs>0)then
       allocate(uindex(nobs),stat=stat_ar)
       if(stat_ar<0)print*,'Allocation for uindex - problem'
       allocate(setnum(nobs),stat=stat_ar)
       if(stat_ar<0)print*,'Allocation for setnum - problem'
       allocate(frmnum(nobs),stat=stat_ar)
       if(stat_ar<0)print*,'Allocation for frmnum - problem'
       !allocate(hkl_resol(nobs))
       !allocate(sig_ar(nobs))
       allocate(x(nobs),stat=stat_ar)
       if(stat_ar<0)print*,'Allocation for x - problem'
       allocate(y(nobs),stat=stat_ar)
       if(stat_ar<0)print*,'Allocation for y - problem'
       allocate(dset(nobs,max_ii),stat=stat_ar)
       if(stat_ar<0)print*,'Allocation for dset - problem'
       allocate(dframe(nobs,max_ii),stat=stat_ar)
       if(stat_ar<0)print*,'Allocation for dframe - problem'
       allocate(dint(nobs,max_ii),stat=stat_ar)
       if(stat_ar<0)print*,'Allocation for dint - problem'
    end if
    !print*,'Allocation...done'

    !Transfering data to arrays...
    i = 0; j = 0
    ob => firstuniq
    sn => firstsn
    fn => firstfn
    inte => firstint
    do while(associated(ob%nextuniq))
       i = i + 1
       setnum(i) = ob%sn%iiset
       frmnum(i) = ob%fn%iframe
       if(i==1) olduindx = ob%uindex
       if(olduindx == ob%uindex)then
          j = j + 1
          k = ob%uindex
          uindex(k) = k
          dset(k,j)   = ob%sn%iiset
          dframe(k,j) = ob%fn%iframe
          dint(k,j)   = ob%inte%iintensity
       elseif(olduindx /= ob%uindex)then
          j = 1
          k = ob%uindex
          uindex(k)   = k
          dset(k,j)   = ob%sn%iiset
          dframe(k,j) = ob%fn%iframe
          dint(k,j)   = ob%inte%iintensity
       end if
       olduindx = ob%uindex
       ob => ob%nextuniq
       sn => sn%nextset
       fn => fn%nextframe
       inte => inte%nextint
    end do

    oldframe = 0
    do is = 1, nset
       do j = 1, nobs
          if(is == setnum(j))then
             if(oldframe < frmnum(j)) oldframe = frmnum(j)
          end if
       end do
       set(is) = oldframe
       write(24,*) is, set(is) !Dataset and #frames
       oldframe = 0
    end do
    nfrms = maxval(set); print*,'Highest frame-number of all datasets=', nfrms

    allocate(np(nfrms,nset)); allocate(CorC(nfrms,nset));allocate(SE(nfrms,nset))
    do i = 1, nset
       do j = 1, nfrms
          CorC(j,i) = 0.
          SE(j,i)   = 0.
          np(j,i)   = 0
       end do
    end do

    !REMARK: It is not easy to include CC calculation in separate subroutine.
    !Many dependancies and hence decided to keep the CC calculation here.

    print*, 'CALCULATING. PLEASE WAIT...'
    call cpu_time(start)
    do i = 1, nset
       do j = 1, set(i)  !Maximum number of frames in ith dataset
          do k = 1, nuniq
             do l = 1, ifq(k)  !Occurence of kth unique reflection in all data
                if((dset(k,l)==i .and. dframe(k,l)==j) .and. (dint(k,l) /= 0))then
                   nx = nx + 1
                   sumxx = sumxx + dint(k,l)
                   present = .true.
                end if
             end do

             if(present)then
                do l = 1, ifq(k)  !Occurence of kth unique reflection in all data
                   if((.not.(dset(k,l)==i .and. dframe(k,l)==j)).and.(dint(k,l) /= 0))then 
                      ny = ny + 1
                      sumyy = sumyy + dint(k,l)
                   end if
                end do
                if(ny /= 0)then
                   x(k) = sumxx/nx
                   y(k) = sumyy/ny
                   npairs = npairs + 1

!!$                 WRITE(25,*) nx,ny,npairs,i,j,k,xx(k),yy(k) !For testing
!!$                 WRITE(26,*) nx,sumxx,xx(k),ny,sumyy,yy(k)  !For testing

                   sumX = sumX + x(k)         ! SUM(Xi)
                   sumY = sumY + y(k)         ! SUM(Yi)
                   sumXY = sumXY + x(k)*y(k) ! SUM(XiYi)

                   sumX2 = sumX2 + x(k)**2    ! SUM(X^2)
                   sumY2 = sumY2 + y(k)**2    ! SUM(Y^2)
                   ! PRINT*, sumX,sumY,sumXY,sumX2,sumY2
                end if
             end if
             present = .false.
             nx = 0
             ny = 0
             sumxx = 0.0
             sumyy = 0.0
          end do
          !Here do the CC calculation.
          !PRINT*, npairs
          if(sumX /= 0 .or. sumY /= 0)then
             AA = sumXY - (sumX * sumY)/npairs
             BB = sumX2 - ((sumX)**2)/npairs
             CC = sumY2 - ((sumY)**2)/npairs
             DD = sqrt(BB*CC)
             CorC(j,i) = AA/DD
             np(j,i) = npairs
!!$           WRITE(27,*) sumX,sumY,sumXY,sumX2,sumY2,AA,DD,CorC(j,i),np(j,i),i,j !For testing
             npairs = 0
             sumX = 0.
             sumY = 0.
             sumXY = 0.
             sumX2 = 0.
             sumY2 = 0.
          end if
          !Standard error in CC
          CC2 = (CorC(j,i)**2)
          if(np(j,i)>2)then
             SE(j,i) = sqrt((1-CC2)/(np(j,i)-2))
          else
             SE(j,i) = 0.
          end if
          !if(np(j,i) > 3) write(6,'(2I5,2F8.3,I5)') i,j,CorC(j,i),SE(j,i),np(j,i) !Screen output
       end do
    end do
    call cpu_time(finish)
    print*, 'CPU time =',finish-start,'  seconds'

    print*, 'PRINTING CC VALUES TO ARRAYS...'
    do i = 1, nset
       do j = 1, nfrms
          if(np(j,i) < 3)then  !Filter out not-defined situations for the CC calculation
             CorC(j,i) = 0.0
             SE(j,i)   = 0.0
             np(j,i)   = 0
          end if
          write(28,'(2I5,2F8.3,I5)') i,j,CorC(j,i),SE(j,i),np(j,i) !Writing to cc_frames.txt
       end do
    end do

    !Here calculates the averaged CC for each frame aveaged over all datasets.
    sumCC = 0.
    ncc = 0  
    do j = 1, nfrms
       do i = 1, nset
          if(CorC(j,i) /= 0)then
             ncc = ncc + 1
             sumCC = sumCC + CorC(j,i)
          end if
       end do
       if(ncc /= 0)then
          write(29,'(I5,F8.3)') j,sumCC/ncc !Writing to cc_averaged.txt
          sumCC = 0.
          ncc = 0
       end if
    end do

    deallocate(uindex)!;print*,'uindex - deallocated'
    deallocate(setnum)! ;print*,'setnum - deallocated'
    deallocate(frmnum)! ;print*,'frmnum - deallocated'
    !deallocate(hkl_resol);print*,'hkl_resol - deallocated'
    !deallocate(sig_ar);print*,'sig_ar - deallocated'
    deallocate(x)!;print*,'x - deallocated'
    deallocate(y)!;print*,'y - deallocated'
    deallocate(dset)!;print*,'dset - deallocated'
    deallocate(dframe)!;print*,'dframe - deallocated'
    deallocate(dint)!;print*,'dint - deallocated'

    print*, 'EVERYTHING IS SUCCESSFUL.' 
    print*, 'Please look at cc_frames.txt AND cc_averaged.txt in your working directory. Good luck!'
    close(unit);close(24);close(28);close(29)

  end subroutine cc_frames



  subroutine xscale2stat(xscalefile)

    character (len=120),intent(in) :: xscalefile
    character*255 :: statfile
    integer       :: ni,nj,nbin,nobs,nuniq
    integer       :: i,j,peak,corr,ii,sumNfreq
    real          :: psi,tempsum

    real,    dimension(:),   allocatable :: bin_arr
    real,    dimension(:,:), allocatable :: sumI,sumSig2
    integer, dimension(:,:), allocatable :: nfreq

1000 FORMAT(A)
1001 FORMAT(3I6,2E11.3,2F8.2,F9.2,F10.6,2I4,F8.3)

    nuniq=0   ! Number of unique reflections
    nobs=0    ! Number of observations
    last=.false.
    ihkl=0;i = 0
    ii=0

    if(rescutoff) print*, 'Resolution cutoffs have enabled!'

    do i=1,len_trim(xscalefile)
       if(xscalefile(i:i)=='.') exit
    end do
    read(xscalefile(1:i-1),1000) filename_template

    call xscaleopen(xscalefile)

    allocate(ob); nullify(ob%nextuniq); firstuniq => ob
    allocate(sn);nullify(sn%nextset); firstsn => sn
    allocate(inte); nullify(inte%nextint);firstint => inte

    i = 0
    do
       read(unit,1000) line
       if(line(1:12) == '!END_OF_DATA') last=.true.
       if(.not.last)then
          if(xscale)then
             read(line,*)ih,ik,il,fsq,sigma,dum4,iset
          else
             read(line,1001)ih,ik,il,fsq,sigma,dum4,peak,corr,psi
          end if
          if(sigma<0.)cycle    
          call getres(ih,ik,il,s2,resol)

          if(rescutoff)then
             if(resol>resol_ll.or.resol<resol_ul) cycle !this removes not matching reflections
          end if

          nobs=nobs+1
          call asuput(ih,ik,il,jh,jk,jl,ihkl,iop,op,size(op,3)/ncenop)
          if(nobs==1) oldhkl=ihkl
          i = i + 1
          if(i==1)then
             nuniq = 1
             lowest_resol = resol
             highest_resol = resol
          else
             if(lowest_resol < resol) lowest_resol = resol
             if(highest_resol > resol) highest_resol = resol
          end if
          if(ihkl /= oldhkl.or.last) nuniq = nuniq+1
          ob%uindex    = nuniq
          ob%sn%iiset    = iset
          ob%inte%iintensity = fsq
          ob%inte%sig = sigma
          ob%hkl_resol = resol
          allocate(ob%nextuniq); ob => ob%nextuniq;nullify(ob%nextuniq)
          allocate(sn%nextset); sn => sn%nextset;nullify(sn%nextset)
          allocate(inte%nextint); inte => inte%nextint;nullify(inte%nextint)
       end if
       if(last) exit
       oldhkl=ihkl
    end do
    print*, 'Reading complete'
    print*,'# obs (excluding misfits), # unique =',nobs,nuniq
    print*, 'lowest resolution=', lowest_resol
    print*, 'highest resolution', highest_resol
    print*, 'Enter the number of resolution bins you want:'
    read(*,*) nbin
    call setbin(nbin,lowest_resol,highest_resol,2,bin_arr)
    print*,'Allocating arrays for data...'
    if(nobs>0)then
       allocate(uindex(nobs))
       allocate(setnum(nobs))
       allocate(intensity(nobs))
       allocate(hkl_resol(nobs))
       allocate(sig_ar(nobs))
       allocate(sumI(nuniq,nset)) !these 2D arrays can be huge depending on nobs.
       allocate(nfreq(nuniq,nset))!therefore, that may lead to memory failures.
       allocate(sumSig2(nuniq,nset))
    end if
    print*,'Allocation...done'

    !Initializing arrays
    do i = 1, nobs
       uindex(i) = 0
       setnum(i) = 0
       intensity(i) = 0.
       sig_ar(i) = 0.
       hkl_resol(i) = 0.
    end do
    do j = 1, nuniq
       do i = 1, nset
          sumI(j,i) = 0.
          nfreq(j,i) = 0
          sumSig2(j,i) = 0.
       end do
    end do

    !Transfering data to arrays
    i = 0
    ob => firstuniq
    sn => firstsn
    inte => firstint
    do while(associated(ob%nextuniq))
       i = i+1
       uindex(i) = ob%uindex
       setnum(i) = ob%sn%iiset
       hkl_resol(i) = ob%hkl_resol
       intensity(i) = ob%inte%iintensity
       sig_ar(i) = ob%inte%sig
       ob => ob%nextuniq
       sn => sn%nextset
       inte => inte%nextint
       !print*, uindex(i),setnum(i),hkl_resol(i),intensity(i)
    end do
    print*,i,nobs
    do i=1,nbin !testing
       print*, bin_arr(i) !testing
    end do !testing

    do i=1,size(bin_arr)
       print*, i,bin_arr(i)
       write(strng,'(I2)') i
       statfile = trim(adjustl(filename_template))//'_stat_bin'//trim(adjustl(strng))//'.txt'
       ii = 0
       do j=1,nobs
          if(i==1)then
             if(hkl_resol(j) >= bin_arr(i))then
                ii = ii + 1
                !write(41,*) uindex(j),setnum(j),intensity(j),sig_ar(j),intensity(j)/sig_ar(j),hkl_resol(j) 
                sumI(uindex(j),setnum(j)) = sumI(uindex(j),setnum(j)) + intensity(j)
                nfreq(uindex(j),setnum(j)) = nfreq(uindex(j),setnum(j)) + 1
                sumSig2(uindex(j),setnum(j)) = sumSig2(uindex(j),setnum(j)) + sig_ar(j)**2
                !print*, sumI(uindex(j),setnum(j)),nfreq(uindex(j),setnum(j))
             endif
          elseif(i > 1)then
             if(hkl_resol(j) >= bin_arr(i) .and. hkl_resol(j) < bin_arr(i-1))then
                ii = ii + 1
                sumI(uindex(j),setnum(j)) = sumI(uindex(j),setnum(j)) + intensity(j)
                nfreq(uindex(j),setnum(j)) = nfreq(uindex(j),setnum(j)) + 1
                sumSig2(uindex(j),setnum(j)) = sumSig2(uindex(j),setnum(j)) + sig_ar(j)**2
             endif
          endif
       end do
       tempsum = 0.
       do j=1,nobs
          if(j==1)then
             write(41,*) j,uindex(j),setnum(j),sumI(uindex(j),setnum(j)),nfreq(uindex(j),setnum(j)),&
                  sumI(uindex(j),setnum(j))/nfreq(uindex(j),setnum(j))
             tempsum = sumI(uindex(j),setnum(j))
          else
             if(sumI(uindex(j),setnum(j))/=tempsum)then
                write(41,*) j,uindex(j),setnum(j),sumI(uindex(j),setnum(j)),nfreq(uindex(j),setnum(j)),&
                     sumI(uindex(j),setnum(j))/nfreq(uindex(j),setnum(j))
                tempsum = sumI(uindex(j),setnum(j))
             else
                cycle
             end if
          end if
       end do
       !print*,sumNfreq
       print*, 'ii=',ii
       !stop
       print*,'Please wait...'
       !print*, 'calculating CC for ', bin_arr(i),i
       !call cc_calc(nset,nuniq,sumI,nfreq,bin_arr(i),ccdfile)
       call calc_ivs(nset,nuniq,sumI,nfreq,sumSig2,bin_arr(i),statfile)
       print*,'Successfully returned!'
       do nj=1,nuniq
          do ni=1,nset
             sumI(nj,ni) = 0.
             nfreq(nj,ni) = 0
             sumSig2(nj,ni) = 0.
          end do
       end do
    end do
    deallocate(uindex)   !;print*,'uindex - deallocated'
    deallocate(setnum)   !;print*,'setnum - deallocated'
    deallocate(intensity)!;print*,'intensity - deallocated'
    deallocate(hkl_resol)!;print*,'resolution - deallocated'
    deallocate(sig_ar)   !;print*,'sigmas - deallocated'
    deallocate(sumI)     !;print*,'sumI - deallocated'
    deallocate(nfreq)    !;print*,'nfreq - deallocated'
    deallocate(sumSig2)
    close(unit)
    print*,'Done'
  end subroutine xscale2stat


  subroutine cc_calc(nset,nuniq,sumI,nfreq,bin_resol,ccdfile_name)

    integer, intent(in) :: nset,nuniq
    real, intent(in)    :: sumI(:,:),bin_resol
    integer, intent(in) :: nfreq(:,:)
    integer             :: n,i,j,ii,nx,ny,npairs
    real                :: CorC,SE
    character*80        :: ccdfile_name
    real                :: sumxx,sumyy,sumX,sumY,sumXY,sumX2,sumY2,AA,BB,CC,DD,CC2
    integer,dimension(:), allocatable :: dataset_i
    real,   dimension(:), allocatable :: cc_i

    ! I assume that I want to calculate the CC between set_ii and the others. In this case
    ! X values are calculated using set ii and Y values are calculated using all other sets.
    ! Select the suitable subset of unique reflections from X and Y populations

    OPEN(UNIT=22,FILE=trim(adjustl(ccdfile_name)),status='replace')
    write(22,*) 'Resolution= ', bin_resol
    write(22,*)'set#   cc   err(cc)   No.of_pairs'

    n = 0
    npairs = 0  ! This variable carries the number of pairs between two subsets (X and Y)
    SE = 0.
    nx = 0; ny = 0

    allocate(x(nuniq));allocate(y(nuniq))
    if(allocated(dataset_i)) deallocate(dataset_i)
    if(allocated(cc_i)) deallocate(cc_i)
    allocate(dataset_i(nset));allocate(cc_i(nset))

    !print*, nset, nuniq

    PRINT*, 'CALCULATING. PLEASE WAIT...'
    CALL cpu_time(start)
    DO ii = 1, nset
       DO j = 1, nuniq
          IF (nfreq(j,ii) > 0) THEN
             nx = nx +1  ! Not used
             DO i = 1, nset
                IF ((i /= ii).AND.(nfreq(j,i) > 0)) THEN
                   ny = ny + nfreq(j,i)  
                   sumyy = sumyy + sumI(j,i)    
                END IF
             END DO
             IF(ny /= 0)THEN
                npairs = npairs + 1
                x(npairs) = sumI(j,ii)/nfreq(j,ii) 
                y(npairs) = sumyy/ny ! Kay method of averaging
                !WRITE(25,'(2I5,2ES15.5,I10)') ii,j,x(npairs),y(npairs),npairs ! Fort testing
             END IF
          END IF
          ny = 0
          sumyy = 0.
       END DO  ! Loop over all unique reflections
       ! print*, npairs

       ! Calculation of correlation coefficeint CC between two sets
       ! CorrelationCoefficient = [{SUM(XY)-(SUMX.SUMY)/n}/sqrt[{SUMX^2-(SUMX)^2/n}{SUMY^2-(SUMY)^2/n}]
       sumX = 0.
       sumY = 0.
       sumXY = 0.
       sumX2 = 0.
       sumY2 = 0.
       AA = 0.
       BB = 0.
       CC = 0.
       DD = 0.
       CorC = 0. ! Correlation Coefficeint
       DO i = 1, npairs
          sumX = sumX + x(i)        ! SUM(Xi)
          sumY = sumY + y(i)        ! SUM(Yi)
          sumXY = sumXY + x(i)*y(i) ! SUM(XiYi)

          sumX2 = sumX2 + x(i)**2   ! SUM(X^2)
          sumY2 = sumY2 + y(i)**2   ! SUM(Y^2)
       END DO
       AA = sumXY - (sumX * sumY)/npairs
       BB = sumX2 - ((sumX)**2)/npairs
       CC = sumY2 - ((sumY)**2)/npairs
       DD = sqrt(BB*CC)
       CorC = AA/DD
       CC2 = (CorC**2)
       SE = sqrt((1-CC2)/(npairs-2))
       dataset_i(ii) = ii
       cc_i(ii) = CorC
       WRITE(22,'(I5,2F8.3,I10)') ii,CorC,SE,npairs ! Writing to cc_datasets.txt
       npairs = 0
    END DO  ! Loop over all datasets
    CALL cpu_time(finish)
    print*, 'CPU time for CC calculation=',finish-start,'  seconds'
    deallocate(x);deallocate(y)
    close(22)
    call plot_1(dataset_i,cc_i)
  end subroutine cc_calc


  subroutine calc_ivs(nset,nuniq,sumI,nfreq,sumSig2,bin_resol,statfile_name)
    implicit none
    integer, intent(in) :: nset,nuniq
    real, intent(in) :: sumI(:,:),sumSig2(:,:),bin_resol
    integer, intent(in) :: nfreq(:,:)
    integer :: i,j,nx,k
    real :: start,finish
    real :: SumBinI,avg_binI,avg_binSig,sumIS,avg_IS,diff_sqd_sum,std_dev_IS
    real, allocatable :: avg_Ihkl(:,:),avg_Sighkl(:,:),avg_I_ovr_Sig(:,:),diff_IS(:)
    character*255 :: statfile_name


    open(22,file=trim(adjustl(statfile_name)),status='replace')
    WRITE(22,*) 'resolution-bin   meanI/sig   Std.Dev.   #uniq-refl.   set #'

    allocate(avg_Ihkl(nuniq,nset))
    allocate(avg_Sighkl(nuniq,nset))
    allocate(avg_I_ovr_Sig(nuniq,nset))
    allocate(diff_IS(nuniq))

    sumIS = 0.
    diff_sqd_sum = 0.
    nx = 0

    call cpu_time(start)
    do i = 1, nset
       nx = 0
       sumIS = 0.
       diff_sqd_sum = 0.
       do j = 1, nuniq
          if(nfreq(j,i) > 0)then
             nx = nx +1  
             avg_Ihkl(nx,i) = sumI(j,i)/nfreq(j,i)
             avg_Sighkl(nx,i) = sqrt(sumSig2(j,i))
             avg_I_ovr_Sig(nx,i) = avg_Ihkl(nx,i)/avg_Sighkl(nx,i)
             write(40,*) avg_Ihkl(nx,i),avg_Sighkl(nx,i),avg_I_ovr_Sig(nx,i) 
          end if
       end do  ! Loop over all unique reflections
       !print*,'nuniq=',nuniq,'nx=',nx
       do k=1,nx
          sumIS = sumIS + avg_I_ovr_Sig(k,i)
       end do
       ! print*, 'Bin=',bin_resol,'I/sig=',sumIS/nx,'#uniq refl.',nx,'set=',i
       ! write(22,*) bin_resol,sumIS/nx,nx,i
       ! calculating std. dev
       avg_IS = sumIS/nx
       do k=1,nx
          diff_IS(k) = avg_I_ovr_Sig(k,i) - avg_IS
          diff_sqd_sum = diff_sqd_sum + diff_IS(k)**2 
       end do
       std_dev_IS = sqrt(diff_sqd_sum/nx)      
       !write(22,'(F5.2,2F8.2,2I8)') bin_resol,sumIS/nx,std_dev_IS,nx,i
       write(22,*) bin_resol,sumIS/nx,std_dev_IS,nx,i
    end do  ! Loop over all datasets
    call cpu_time(finish)
    print*, ' CPU time for calculation=',finish-start,'  seconds'
    deallocate(avg_Ihkl);deallocate(avg_Sighkl)
    deallocate(avg_I_ovr_Sig);deallocate(diff_IS)
    close(22);close(40)

  end subroutine calc_ivs


  subroutine xscaleopen(xscalefile)

    character (len=120),intent(in) :: xscalefile

1000 FORMAT(A)

    xscale = .false.

    unit=10
    OPEN(UNIT=unit,FILE=trim(adjustl(xscalefile)),STATUS='old',ACTION='read',iostat=ier)
    if(ier/=0)then
       print*,'ier=',ier
       stop 'could not open'
    end if
    ncheck=0
    do
       read(unit,1000,iostat=ier) line
       if(ier/=0)then
          print*,'ier=',ier
          print*,'line=',line(:LEN_TRIM(line))
          stop 'read error in xds_file'
       end if
       if(line(1:20)=='!SPACE_GROUP_NUMBER=')then
          read(line(21:),*) ispgr
          call getops(ispgr,op)
          write(*,'(a)'),line(:len_trim(line))
       else if(line(1:21).eq.'!UNIT_CELL_CONSTANTS=')then 
          read(line(22:),*) a,b,c,alf,bet,gam
          write(*,'(a)'),line(:len_trim(line))
       else if(line(1:9) == '!ITEM_H=1')then
          ncheck=ncheck+1
       else if(line(1:9) == '!ITEM_K=2')then
          ncheck=ncheck+1
       else if(line(1:9) == '!ITEM_L=3')then
          ncheck=ncheck+1
       else if(line(1:12) == '!ITEM_IOBS=4')then
          ncheck=ncheck+1
       else if(line(1:19) == '!ITEM_SIGMA(IOBS)=5')then
          ncheck=ncheck+1
       else if(line(1:13) == '!ITEM_ISET=10')then
          ncheck=ncheck+1
          xscale = .true.
       else if(line(1:7) == '! ISET=')then
          read(line(8:),*) iset
          nset=max(nset,iset)
       else if(line(1:14) == '!END_OF_HEADER')then
          exit
       end if
    end do
    if(ncheck == 6 .and. xscale) print*, 'this is an xscale input'
    if(ncheck == 5 .and. (.not.xscale)) print*, 'this is an xds_ascii input'
    if(ncheck /= 5 .and. (.not.xscale)) stop 'header items not found'
    if(.not.xscale)then
       iset = 1;nset = 1
    end if

    !if(ncheck /= 6) stop 'header items not found'
    print*,'nset=',nset
  end subroutine xscaleopen


  subroutine cchalf_xscale(xscalefile)

    integer,parameter :: maxref=50000
    character (len=120),intent(in) :: xscalefile
    character*255 :: cchfile
    integer       :: ni,nj,nbin,nobs,nuniq,nx,peak,corr
    integer       :: i,j,k,l,m,n,ii,jj,max_ii,olduindx,npairs
    real          :: sumx,sumy,CorC,SE,psi
    integer       :: ifq(maxref)
    real,    dimension(:),allocatable :: bin_arr,xx,yy
    real,    dimension(:,:), allocatable :: dint,resol_ar


1000 FORMAT(A)
1001 FORMAT(3I6,2E11.3,2F8.2,F9.2,F10.6,2I4,F8.3)

    ii = 1;jj = 1
    max_ii = 1
    nuniq=0   ! Number of unique reflections
    nobs=0    ! Number of observations
    last=.false.
    ihkl=0;i=0

    if(rescutoff) print*, 'Resolution cutoffs have enabled!'

    do i=1,len_trim(xscalefile)
       if(xscalefile(i:i)=='.') exit
    end do
    read(xscalefile(1:i-1),1000) filename_template

    call xscaleopen(xscalefile)

    allocate(ob); nullify(ob%nextuniq); firstuniq => ob
    allocate(inte); nullify(inte%nextint);firstint => inte

    i = 0
    do
       read(unit,1000) line
       if(line(1:12) == '!END_OF_DATA') last=.true.
       if(.not.last)then
          if(xscale)then
             read(line,*) ih,ik,il,fsq,sigma,dum4,iset
          else
             !print*,'This file is not a xscale.ahkl type file.'
             read(line,1001)ih,ik,il,fsq,sigma,dum4,peak,corr,psi
             !stop
          end if
          !if(sigma<0.) cycle
          !if((fsq/sigma)<0.) cycle
          call getres(ih,ik,il,s2,resol)

          if(rescutoff)then
             if(resol>resol_ll.or.resol<resol_ul) cycle !this removes not matching reflections
          end if

          nobs=nobs+1
          call asuput(ih,ik,il,jh,jk,jl,ihkl,iop,op,size(op,3)/ncenop)
          if(nobs==1)then
             oldhkl=ihkl
             nuniq = 1
             lowest_resol = resol
             highest_resol = resol
          else
             if(lowest_resol < resol) lowest_resol = resol
             if(highest_resol > resol) highest_resol = resol
          end if
          if(ihkl == oldhkl)then
             if(nobs == 1) ii = 0
             i = nuniq
             ii = ii + jj
             if(ii > max_ii) max_ii = ii
             ifq(i) = ii ! Choice for better performance
             ob%uindex = nuniq
             ob%hkl_resol = resol
             ob%inte%iintensity = fsq
             !WRITE(22,*) i,ii,ob%uindex,ob%setnum%iiset,ob%framenum%iframe,ob%intensity%iintensity !For testing
          else ! IF (ihkl /= oldhkl) THEN
             nuniq = nuniq+1
             j = nuniq
             jj = 1
             ifq(j) = jj ! Choice for better performance         
             ob%uindex = nuniq
             ob%hkl_resol = resol
             ob%inte%iintensity = fsq
             !WRITE(22,*) j,jj,ob%uindex,ob%setnum%iiset,ob%framenum%iframe,ob%intensity%iintensity !For testing
             ii = 1
          end if
          allocate(ob%nextuniq);  ob => ob%nextuniq;  nullify(ob%nextuniq)
          allocate(inte%nextint);  inte => inte%nextint; nullify(inte%nextint)
       end if
       if(last) exit
       oldhkl=ihkl
       !print*, ih,ik,il,resol,lowest_resol,highest_resol
    end do
    xscale = .false.
    print*, 'Reading complete'
    print*,'# obs (excluding misfits), # unique =',nobs,nuniq
    print*, 'lowest resolution=', lowest_resol
    print*, 'highest resolution', highest_resol
    print*,
    print*, 'Enter the number of resolution bins you want:'
    read(*,*) nbin
    call setbin(nbin,lowest_resol,highest_resol,2,bin_arr)
    print*,'Allocating arrays for data...'
    if(nobs>0)then
       allocate(uindex(nuniq))!;print*,'Allocation for uindex - done'
       allocate(dint(nuniq,max_ii))
       allocate(resol_ar(nuniq,max_ii))
    end if
    print*,'Allocation...done'

    !Transfering data to arrays...
    i = 0; j = 0
    ob => firstuniq
    inte => firstint
    do while(associated(ob%nextuniq))
       i = i + 1
       if(i==1) olduindx = ob%uindex
       if(olduindx == ob%uindex)then
          j = j + 1
          k = ob%uindex
          uindex(k) = k
          resol_ar(k,j) = ob%hkl_resol
          dint(k,j)   = ob%inte%iintensity
       elseif(olduindx /= ob%uindex)then
          j = 1
          k = ob%uindex
          uindex(k)   = k
          resol_ar(k,j) = ob%hkl_resol
          dint(k,j)   = ob%inte%iintensity
       end if
       olduindx = ob%uindex
       ob => ob%nextuniq
       inte => inte%nextint
    end do
    print*, 'Transfering to arrays - Done'
    cchfile = trim(adjustl(filename_template))//'cch_bin.txt'
    open(20,file=trim(adjustl(cchfile)),status='replace')
    write(20,*)'Resolution   Corr    Err(CC)   Npairs'
    do i=1,size(bin_arr)
       !print*, i,bin_arr(i)
       write(strng,'(I2)') i
       allocate(x(nuniq));allocate(y(nuniq))
       do k=1,nuniq
          x(k) = 0.
          y(k) = 0.
       end do
       nx = 0
       do k = 1, nuniq
          m=0;n=0
          sumx=0.;sumy=0.
          do l = 1, ifq(k)  !Occurence of kth unique reflection in all data
             if(i==1)then
                if(resol_ar(k,l) >= bin_arr(i))then
                   if(mod(l,2)==0)then
                      n = n + 1
                      sumx = sumx + dint(k,l)
                      write(49,*) 'k=',k,'x=',dint(k,l)
                   else
                      m = m + 1
                      sumy = sumy + dint(k,l)
                      write(49,*) 'k=',k,'y=',dint(k,l)
                   end if
                end if
             else
                if(resol_ar(k,l) >= bin_arr(i) .and. resol_ar(k,l) < bin_arr(i-1))then
                   if(mod(l,2)==0)then
                      n = n + 1
                      sumx = sumx + dint(k,l)
                   else
                      m = m + 1
                      sumy = sumy + dint(k,l)
                   end if
                end if
             end if
          end do
          if(n>0)then
             nx = nx + 1
             x(nx) = sumx/n
             y(nx) = sumy/m
          end if
       end do
       allocate(xx(nx));allocate(yy(nx))
       do j=1,nx
          write(50,*) j,x(j),y(j)
          xx(j) = x(j)
          yy(j) = y(j)
       end do
       write(50,*)'calculating CChalf for ', bin_arr(i),i
       call cch_calc(bin_arr(i),xx,yy,CorC,SE,npairs)
       write(20,'(3F8.3,I10)') bin_arr(i),CorC,SE,npairs
       deallocate(x);deallocate(y)
       deallocate(xx);deallocate(yy)
    end do
    deallocate(uindex)   !;print*,'uindex - deallocated'
    deallocate(resol_ar)!;print*,'resolution - deallocated'
    close(unit);close(20)
    print*,'Done'
  end subroutine cchalf_xscale

  subroutine cch_calc(bin,x,y,CorC,SE,npairs)
    implicit none
    real,intent(in) ::bin
    real,dimension(:),intent(in) :: x,y
    real :: sumyy,sumX,sumY,sumXY,sumX2,sumY2,AA,BB,CC,DD,CC2
    real,intent(out) :: CorC,SE
    integer :: i
    integer,intent(out) :: npairs

    sumX = 0.
    sumY = 0.
    sumXY = 0.
    sumX2 = 0.
    sumY2 = 0.
    AA = 0.
    BB = 0.
    CC = 0.
    DD = 0.
    CorC = 0.
    npairs = 0
    do i=1,size(x)
       npairs = npairs + 1
       sumX = sumX + x(i)
       sumY = sumY + y(i)
       sumXY = sumXY + x(i)*y(i)
       sumX2 = sumX2 + x(i)**2
       sumY2 = sumY2 + y(i)**2
    end do
    AA = sumXY - (sumX * sumY)/npairs
    BB = sumX2 - ((sumX)**2)/npairs
    CC = sumY2 - ((sumY)**2)/npairs
    DD = sqrt(BB*CC)
    CorC = AA/DD
    CC2 = (CorC**2)
    SE = sqrt((1-CC2)/(npairs-2))
    WRITE(6,'(3F8.3,I10)') bin,CorC,SE,npairs
  end subroutine cch_calc


  subroutine cc_set_remover

  IMPLICIT NONE
  CHARACTER (LEN=132) :: line=''
  INTEGER :: unit,ier,ier1,ier2,iset,np,i,j,nset,lstln
  INTEGER :: narg,cnt
  REAL cutoff,cci,err
  CHARACTER*120 :: ccdfile,statfile,inpfile,inpfile2
  CHARACTER*1 :: choice
  LOGICAL inpline_occured,ccset,ivs
  REAL resi,ivsi,stdi
  INTEGER uniqi

  TYPE :: cc
     integer i_set !set#
     real    cc_i  !cc set
     real    err_i !err in cc
     integer np_i  !#common pairs
     real    res_i !resolution
     real    ivs_i !mean(I/Sig)
     real    std_i !std ivs
     integer uniq_i!#uniq
     type(cc),pointer :: nextcc_i
  END type cc

  TYPE(cc), POINTER :: ob,firstcc

  INTEGER, DIMENSION(:),ALLOCATABLE :: i_set
  REAL   , DIMENSION(:),ALLOCATABLE :: cc_i,ivs_i

1000 FORMAT(A)

  ccset = .false.
  ivs = .false.

  narg = IARGC()
!!$    IF(narg < 1)THEN
!!$       PRINT*, 'usage: ./foo ccdfile inpfile'
!!$       PRINT*, 'ccdfile = your cc_dataset.txt'
!!$       PRINT*, 'inpfile = your corresponding XSCALE.INP'
!!$       STOP
!!$    ELSE
!!$       CALL GETARG(1,ccdfile)
!!$       CALL GETARG(2,inpfile)
!!$    END IF

20 WRITE(*,*) 'Dataset removing is based on either CC-dataset or mean(I/Sig) of each dataset'
  WRITE(*,*) 'Please make your choice, [C] cc-dataset  [S] mean(I/Sig):'
  READ(*,*) choice
  IF(choice=='C'.or.choice=='c')THEN
     ccset = .true.
     ivs = .false.
     WRITE(*,*) 'Your choice - CC data set'
     WRITE(6,*) 'Enter the name of your ccd file:'
     READ(5,*) ccdfile
  END IF
  IF(choice=='S'.or.choice=='s')THEN
     ivs = .true.
     ccset = .false.
     WRITE(*,*) 'Your choice - mean(I/Sig)'
     WRITE(*,*) 'Enter the name of your stat file:'
     READ(5,*) statfile
  END IF
  !IF(choice/='C'.or.choice/='c'.or.choice/='S'.or.choice/='s') go to 20
  if(.not.(ccset.or.ivs)) go to 20

  WRITE(6,*) 'Enter the name of corresponding XSCALE.INP file:'
  READ(5,*) inpfile

  IF(trim(adjustl(inpfile))=='XSCALE_newccd.INP')THEN
     inpfile2 = trim(adjustl(inpfile))//"_old"
     CALL SYSTEM("mv "//trim(adjustl(inpfile))//" "//trim(adjustl(inpfile2)))
  ELSE
     inpfile2 = inpfile
  END IF

  IF(ccset)THEN
     OPEN(UNIT=10,FILE=trim(adjustl(ccdfile)),STATUS='old',ACTION='READ',iostat=ier1)
  END IF
  IF(ivs)THEN
     OPEN(UNIT=10,FILE=trim(adjustl(statfile)),STATUS='old',ACTION='READ',iostat=ier1)
  END IF

  OPEN(UNIT=11,FILE=trim(adjustl(inpfile2)),STATUS='old',ACTION='READ',iostat=ier2)

  IF (ier1/=0 .OR. ier2/=0) THEN
     PRINT*,'ier1=',ier1
     PRINT*,'ier2=',ier2
     STOP 'could not open file'
  END IF

  OPEN(UNIT=20,FILE='XSCALE_newccd.INP',STATUS='unknown')
  WRITE(20,1000) 'OUTPUT_FILE=XSCALE_newccd.ahkl'

  DO
     READ(11,*) line
     IF(trim(adjustl(line(1:11)))=='OUTPUT_FILE') CYCLE
     IF(trim(adjustl(line(1:10)))=='INPUT_FILE')THEN
        BACKSPACE(11); EXIT
     END IF
     WRITE(20,1000) trim(adjustl(line))
  END DO

  IF(ccset)THEN
     DO i = 1, 2
        READ(10,*) line
     END DO
  END IF
  IF(ivs)THEN
     DO i = 1, 1
        READ(10,*) line
     END DO
  END IF

  ALLOCATE(ob); NULLIFY(ob%nextcc_i); firstcc => ob
  nset = 0

  DO
     IF(ccset) READ(10,*,iostat=ier) iset, cci, err, np
     IF(ivs)   READ(10,*,iostat=ier) resi,ivsi,stdi,uniqi,iset
     IF(ier > 0) THEN
        STOP 'Something wrong in ccd/stat file!'
     ELSE IF(ier == 0 .and. ccset)THEN
        nset = nset + 1
        ob%i_set = iset
        ob%cc_i  = cci
        ob%err_i = err
        ob%np_i  = np
        ALLOCATE(ob%nextcc_i); ob => ob%nextcc_i; NULLIFY(ob%nextcc_i)
     ELSE IF(ier == 0 .and. ivs)THEN
        nset = nset + 1
        ob%res_i = resi
        ob%ivs_i = ivsi
        ob%std_i = stdi
        ob%uniq_i = uniqi
        ob%i_set = iset
        ALLOCATE(ob%nextcc_i); ob => ob%nextcc_i; NULLIFY(ob%nextcc_i)
     ELSE IF(ier < 0)THEN
        EXIT !EOF occured
     END IF
  END DO

  print*, 'Number of data sets = ',nset
  ALLOCATE(i_set(nset))
  IF(ccset)THEN
     IF(.NOT.ALLOCATED(cc_i)) ALLOCATE(cc_i(nset))
  END IF
  
  IF(ivs)THEN
     IF(.NOT.ALLOCATED(ivs_i)) ALLOCATE(ivs_i(nset))
  END IF
  
  i = 0
  ob => firstcc
  DO WHILE(ASSOCIATED(ob%nextcc_i))
     i = i + 1
     i_set(i) = ob%i_set
     if(ccset) cc_i(i)  = ob%cc_i
     if(ivs)   ivs_i(i) = ob%ivs_i
     ob => ob%nextcc_i
  END DO

  IF(ccset) PRINT*, 'What is your CC cutoff?'
  IF(ivs) PRINT*, 'What is your I/Sig cutoff?'
  READ(*,*) cutoff

  inpline_occured = .false.

  cnt = 0
  DO i = 1, nset
     inpline_occured = .false.
     IF((ccset.and.(cc_i(i) > cutoff)).OR.(ivs.and.(ivs_i(i) > cutoff)))THEN
        ! IF(ivs.and.(ivs_i(i) > cutoff)) go to 21

        cnt = cnt + 1
        DO
           READ(11,1000,iostat=lstln) line
           IF((trim(adjustl(line(1:10)))=='INPUT_FILE') .AND. (.not.inpline_occured))THEN
              WRITE(20,1000) trim(adjustl(line))
              inpline_occured = .true.
              cycle
           END IF
           IF(inpline_occured)THEN
              IF(trim(adjustl(line(1:10))) /= 'INPUT_FILE')THEN
                 IF(lstln /= 0)THEN
                    EXIT
                 ELSE
                    WRITE(20,1000) trim(adjustl(line))
                    CYCLE
                 END IF
              END IF
              IF(trim(adjustl(line(1:10))) == 'INPUT_FILE')THEN
                 inpline_occured = .false.
                 BACKSPACE(11) !This sets the record pos to the previous line
                 EXIT
              END IF
           END IF
        END DO
     ELSE
        DO
           READ(11,1000,iostat=lstln) line
           IF((trim(adjustl(line(1:10)))=='INPUT_FILE') .AND. (.not.inpline_occured))THEN
              !WRITE(20,1000) trim(adjustl(line))
              inpline_occured = .true.
              cycle
           END IF
           IF(inpline_occured)THEN
              IF(trim(adjustl(line(1:10))) /= 'INPUT_FILE')THEN
                 IF(lstln /= 0)THEN
                    EXIT
                 ELSE
                    !WRITE(20,1000) trim(adjustl(line))
                    CYCLE
                 END IF
              END IF
              IF(trim(adjustl(line(1:10))) == 'INPUT_FILE')THEN
                 inpline_occured = .false.
                 BACKSPACE(11) !This sets the record pos to the previous line
                 EXIT
              END IF
           END IF
        END DO
        
        CYCLE
     END IF
  END DO

  PRINT*, 'Number of data sets whose CC or I/Sig greater than ',cutoff, ' is ',cnt
  IF(ALLOCATED(i_set)) DEALLOCATE(i_set)
  IF(ALLOCATED(cc_i))  DEALLOCATE(cc_i)
  IF(ALLOCATED(ivs_i)) DEALLOCATE(ivs_i)
  CLOSE(10);CLOSE(11);CLOSE(20)

  end subroutine cc_set_remover


END MODULE subroutines
  
