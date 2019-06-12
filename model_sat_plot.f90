program read_test
    use netcdf
    implicit none
    character(len=100) :: region_mask_file, input_path, output
    character(len=100) :: sat_file_name,org_file
    integer(kind=4), parameter :: imt = 600, jmt = 640, thmax=8
    integer(kind=4), dimension(imt,jmt) :: region_mask, temp_i4
    real(kind=8), dimension(imt,jmt) :: tarea, iarea, ithick,hi_lvl
    real(kind=8), dimension(imt,jmt,5) :: hi_tmp, aice_tmp
    integer(kind=4) :: ncid, vid, i, j, myunit, k, l, m
    integer(kind=4), dimension(3) :: init_date, date, end_date
    integer(kind=4) :: date2(6)
    character(len=60) :: filen,box_file,bin_file,command,ln
    logical :: file_exists,first,divid
    logical, allocatable :: box_ex(:)
    integer(kind=1) days_in_month(12)
    integer :: year,month,varid(2),x_dimid,y_dimid
    real(kind=8), dimension(imt,jmt) :: tlon,tlat
    real(kind=8) :: area_tot, vol_tot, area_reg(3), vol_reg(3)
    integer(kind=4), dimension(imt,jmt) :: rect_mask,box_mask,wb
    integer(kind=4) :: iostatus, nmax,nmax2,nmax3
    real(kind=8), dimension(:), allocatable :: hi_box,cnt
    integer(kind=4), dimension(:), allocatable :: box_index
    real(kind=8) :: dummy, ice_thick_dist(thmax),treshold,fillvalue
    character(len=10) :: data_type(2)
    real(kind=8), dimension(:,:), allocatable :: thick_bins
    integer(kind=4), dimension(:), allocatable :: tot_ice_cells,tot_cells,tot_cells2
    integer(kind=4), dimension(:,:), allocatable :: ice_cells
    data ice_thick_dist /0.05,2.5,6.5,12.5,20.5,30.5,40.5,50.0/
    integer, allocatable, dimension(:,:) :: time,tmp_time 
    real(kind=8), allocatable, dimension(:,:) :: bins,tmp_bins 
    real(kind=8), allocatable, dimension(:,:) :: lvl_thick_bins, lvl_ice_cells
    real(kind=8), allocatable, dimension(:) :: hi_lvl_box,tot_lvl_ice_cells
    real(kind=8) :: aice,vice
    character(len=20) :: reg_name(4)

    reg_name = (/"total","Bothnian_Bay","Bothnian_Sea","Gulf_of_Finland"/)
    100 format( i1,a,*(f6.2,a) )

    output = 'test_albedo'
    treshold = ice_thick_dist(1)
    divid = .false.

    days_in_month = (/31,28,31,30,31,30,31,31,30,31,30,31/)

!    open(unit=8, file="rect_corners.txt")
!    nmax = 0
!    do
!        READ(8, *, IOSTAT=iostatus) dummy
!        IF (iostatus < 0) EXIT
!        nmax = nmax + 1
!    end do
!    close(8)
    nmax = 13

    open(unit=8, file="box_vrt.txt")
    nmax3 = 0
    do
        READ(8, *, IOSTAT=iostatus) dummy
        IF (iostatus < 0) EXIT
        nmax3 = nmax3 + 1
    end do
    close(8)

    open(8,file = 'cat.txt',status='replace')
    write(8,100) thmax, achar(9), ( ice_thick_dist(m), achar(9), m=1,thmax )
    close(8)
   
 
    allocate( hi_box(nmax), box_index(nmax), tot_cells(nmax), box_ex(nmax), cnt(nmax) )
    allocate( ice_cells(nmax,thmax), tot_ice_cells(nmax), thick_bins(nmax,thmax) )
    allocate( hi_lvl_box(nmax), lvl_thick_bins(nmax,thmax))
    allocate( tot_lvl_ice_cells(nmax), lvl_ice_cells(nmax,thmax),tot_cells2(nmax))

!###############################################################
!read region mask and tarea
    region_mask_file = '/scratch/lustre/plgjjakacki/PLGRID_NG/&
        &archive/PLGNG002/ocn/hist/PLGNG002.pop.h.2018-11-01-03600.nc'

    region_mask = 0
    call check( nf90_open( trim(region_mask_file), NF90_NOWRITE, ncid ), 100 )
    call check( nf90_inq_varid( ncid, 'REGION_MASK', vid ), 200 )
    call check( nf90_get_var( ncid, vid, temp_i4 ), 300 )
    call check( nf90_inq_varid( ncid, 'TAREA', vid ), 400 )
    call check( nf90_get_var( ncid, vid, tarea ), 500 )
    call check( nf90_close( ncid ), 600 )

    tarea = tarea*0.0001 !conversion from cm^2 to m^2
    do i = 1, imt
        do j = 1, jmt
            if ( temp_i4(i,j) .eq. 1 ) region_mask(i,j) = 1!Bothnian Bay
            if ( temp_i4(i,j) .eq. 2 ) region_mask(i,j) = 2!Bothnian Sea
            if ( temp_i4(i,j) .eq. 3 ) region_mask(i,j) = 3!Gulf of Finland
            if ( temp_i4(i,j) .ne. 0 .and. region_mask(i,j).eq.0 ) region_mask(i,j) = 4!Baltic
        end do
    end do
!###############################################################

    CALL check( nf90_open( 'grid.nc', NF90_NOWRITE, ncid ), 1 )
    CALL check( nf90_inq_varid( ncid, 'TLON', varid(1) ), 2 )
    CALL check( nf90_inq_varid( ncid, 'TLAT', varid(2) ), 3 )
    CALL check( nf90_get_var( ncid, varid(1), tlon ), 5 )
    CALL check( nf90_get_var( ncid, varid(2), tlat ), 6 )
    CALL check( nf90_close( ncid ), 8 )

    box_mask = 0
    where(region_mask.eq.0) box_mask = -99
    call create_box(nmax,nmax3,tlon,tlat,box_mask,box_index,region_mask,tot_cells)
    call create_box(nmax,nmax3,tlon,tlat,rect_mask,box_index,region_mask,tot_cells2)

!###############################################################
CALL check( nf90_create( 'boxes.nc',NF90_CLOBBER,ncid ), 23 )
CALL check( nf90_def_dim( ncid,    'ni', 600, x_dimid ), 24 )
CALL check( nf90_def_dim( ncid,    'nj', 640, y_dimid ), 25 )
CALL check( nf90_def_var( ncid,  'BOXES', NF90_INT, (/x_dimid,y_dimid/), varid(1) ), 27 )
CALL check( nf90_put_att( ncid, varid(1), 'long_name','box number' ), 31 )
CALL check( nf90_put_att( ncid, varid(1), '_FillValue',-99 ), 31 )
CALL check( nf90_enddef(ncid), 44 )
CALL check( nf90_put_var( ncid, varid(1), box_mask ), 45 )
CALL check( nf90_close(ncid), 49 )

wb = 0
where(region_mask.eq.0) wb = -99
CALL check( nf90_create( 'without_boxes.nc',NF90_CLOBBER,ncid ), 23 )
CALL check( nf90_def_dim( ncid,    'ni', 600, x_dimid ), 24 )
CALL check( nf90_def_dim( ncid,    'nj', 640, y_dimid ), 25 )
CALL check( nf90_def_var( ncid,  'wb', NF90_INT, (/x_dimid,y_dimid/),varid(1) ), 27 )
CALL check( nf90_put_att( ncid, varid(1), 'long_name','baltic mask' ), 31 )
CALL check( nf90_put_att( ncid, varid(1), '_FillValue',-99 ), 31 )
CALL check( nf90_enddef(ncid), 44 )
CALL check( nf90_put_var( ncid, varid(1), wb ), 45 )
CALL check( nf90_close(ncid), 49 )

!###############################################################

!STOP
!###############################################################
!year, month, day of beggining and end of date
    init_date(1) = 2017 !year of the beginning
    init_date(2) = 12 !month of the beginning
    init_date(3) = 1 !day of the beginning
    end_date(1) = 2018 !year of the end
    end_date(2) = 5 !month of the end
    end_date(3) = 31 !day of the end
!###############################################################


!###############################################################
!mask for boxes
!    CALL create_rectangles(rect_mask,imt,jmt,nmax,box_index,region_mask,tot_cells)
!    print*, box_index
!###############################################################


!###############################################################
!read ice area and volume
    1000 format( 'ls ',i4.4,'/',i2.2,'/600x640_ice_thickness_baltic*.nc>sat_file' )
    1100 format( 'iceh.',i4.4,'-',i2.2,'-',i2.2,'.nc' )
    1200 format(a,'/model_box',i2.2,'.txt')
    1250 format(a,'/model_bins_box',i2.2,'.txt')
    1203 format(a,'/model_lvl_box',i2.2,'.txt')
    1253 format(a,'/model_lvl_bins_box',i2.2,'.txt')
    1201 format(a,'/sat_box',i2.2,'.txt')
    1251 format(a,'/sat_bins_box',i2.2,'.txt')
    1202 format(a,'/s-1_sat_box',i2.2,'.txt')
    1252 format(a,'/s-1_sat_bins_box',i2.2,'.txt')
    1300 format(i4.4,a,i2.2,a,i2.2,a,f13.6,a,f13.6)
    1400 format(i4.4,a,i2.2,a,i2.2,*(a,f13.6))
    1301 format(i4.4,a,5(i2.2,a),f13.6,a,f13.6)
    1401 format(i4.4,5(a,i2.2),*(a,f13.6))
    1500 format(i4.4,a,i2.2,a,i2.2,a,f13.6,a,f13.6)
    1501 format(i4.4,a,5(i2.2,a),f13.6,a,f13.6)


    input_path = '/scratch/lustre/plgmacmuz/CICE_WKDIR/'//trim(output)//'/history/'
   
    do i = 1, nmax
        write(box_file,1200) trim(output),box_index(i)
        write(bin_file,1250) trim(output),box_index(i)
        myunit=100
        open(myunit+i,file = trim(box_file),status='replace')
        myunit=200
        open(myunit+i,file = trim(bin_file),status='replace')
    end do

    do i = 1, nmax
        write(box_file,1201) trim(output),box_index(i)
        write(bin_file,1251) trim(output),box_index(i)
        myunit=300
        open(myunit+i,file = trim(box_file),status='replace')
        myunit=400
        open(myunit+i,file = trim(bin_file),status='replace')
    end do

    if (divid) then
        do i = 1, nmax
            write(box_file,1202) trim(output),box_index(i)
            write(bin_file,1252) trim(output),box_index(i)
            myunit=500
            open(myunit+i,file = trim(box_file),status='replace')
            myunit=600
            open(myunit+i,file = trim(bin_file),status='replace')
        end do
    end if


    do i = 1, nmax
        write(box_file,1203) trim(output),box_index(i)
        write(bin_file,1253) trim(output),box_index(i)
        myunit=700
        open(myunit+i,file = trim(box_file),status='replace')
        myunit=800
        open(myunit+i,file = trim(bin_file),status='replace')
    end do

    date = init_date
    month = 0
    do 
        if (date(2).ne.month) then
            month = date(2)
!            write(*,*) month

            write(command,1000) date(1), month
!            write(*,*) command
            call execute_command_line( trim(command) )

            open(unit=8, file="sat_file")
            nmax2 = 0
            do
                READ(8, *, IOSTAT=iostatus) dummy
                IF (iostatus < 0) EXIT
                nmax2 = nmax2 + 1
            end do
            close(8)

            open(9, file="sat_file", status='old')
            do l = 1, nmax2
                read(9,'(A)') sat_file_name
                write(*,*) trim(sat_file_name)
!                write(*,*) sat_file_name(38:51)
                org_file = sat_file_name(1:8)//sat_file_name(17:)
!                write(*,*) trim(org_file)
                read(sat_file_name(38:41),*) date2(1)
                read(sat_file_name(42:43),*) date2(2)
                read(sat_file_name(44:45),*) date2(3)
                read(sat_file_name(46:47),*) date2(4)
                read(sat_file_name(48:49),*) date2(5)
                read(sat_file_name(50:51),*) date2(6)
!                do k = 1, 6
!                    write(*,*) date2(k)
!                end do
!                date2(1) = sat_file_name(38:41)
!                date2(2) = sat_file_name(42:43)
!                date2(3) = sat_file_name(44:45)
!                date2(4) = sat_file_name(46:47)
!                date2(5) = sat_file_name(48:49)
!                date2(6) = sat_file_name(50:51)
!                write(*,*) trim(sat_file_name)
                call check( nf90_open( trim(sat_file_name), NF90_NOWRITE, ncid ), 1 )
                call check( nf90_inq_varid( ncid, 'HI' , vid ), 6 )
                call check( nf90_get_var( ncid, vid, ithick ), 12 )
                call check( nf90_get_att( ncid, vid, '_FillValue', fillvalue ), 21 )
!                call check( nf90_get_att( ncid, vid, 'long_name', ln ), 21 )
                call check( nf90_close( ncid ), 23 )

                call check( nf90_open( trim(org_file), NF90_NOWRITE, ncid ), 1 )
                call check( nf90_inq_varid( ncid, 'sea_ice_thickness' , vid ), 6 )
                call check( nf90_get_att( ncid, vid, 'long_name', ln ), 21 )
                call check( nf90_close( ncid ), 23 )
!                write(*,*) trim(sat_file_name),index(ln,'S-1')
!                write(*,*) ln

                hi_box = 0
                tot_ice_cells = 0
                ice_cells = 0
                tot_cells = 0
                box_ex = .true.
                cnt = 0
    
                do i = 1, imt
                    do j = 1, jmt
                        do k = 1, nmax
                            if ( rect_mask(i,j).eq.k .and. ithick(i,j).eq.fillvalue ) then
                                box_ex(k) = .false.
                                cnt(k) = cnt(k)+1
                            end if                
                        end do!boxes
                    end do!j
                end do!i

                do k = 1, nmax
                    if (.not.box_ex(k).and.(cnt(k)/real(tot_cells2(k))).lt.0.09) then
!                        write(*,*) sat_file_name, box_index(k),cnt(k)
                        call fill_empty(imt,jmt,k,rect_mask,ithick,fillvalue)
                        box_ex(k) = .true.
                        if (k.eq.1.) print*, cnt(k)/real(tot_cells2(k))
                    end if
                end do

                do i = 1, imt
                    do j = 1, jmt
                        do k = 1, nmax
                            if ( rect_mask(i,j).eq.k .and. box_ex(k) ) then
                                tot_cells(k) = tot_cells(k) + 1
                                hi_box(k) = hi_box(k)+ithick(i,j)
                                if ( ithick(i,j)*100.gt.treshold ) then
                                    tot_ice_cells(k) = tot_ice_cells(k)+1
                                end if
                                do m = 1, thmax-1
                                    if ( ithick(i,j)*100.gt.ice_thick_dist(m) .and.&
                                        ithick(i,j)*100.le.ice_thick_dist(m+1) ) then
                                        ice_cells(k,m) = ice_cells(k,m)+1         
                                    end if
                                end do
                                if ( ithick(i,j)*100.gt.ice_thick_dist(thmax) ) then
                                    ice_cells(k,thmax) = ice_cells(k,thmax)+1         
                                end if

                            end if
                        end do!boxes
                    end do!j
                end do!i
                 
                do k = 1, nmax
                    hi_box(k) = hi_box(k)/real(tot_cells(k))
                    if (tot_ice_cells(k).gt.0.0) then
                        do m = 1, thmax
                            thick_bins(k,m) = real(ice_cells(k,m))/real(tot_ice_cells(k))
                        end do
                    else
                        hi_box(k) = 0
                        thick_bins(k,:) = 0
                    end if
                end do
    
                do k = 1, nmax
                    if ( box_ex(k) ) then

                        if ( index(ln,'S-1').eq.0 ) then
                           myunit=300  
                        else
                            if (.not.divid) then
                                myunit=300  
                            else
                                myunit=500  
                            end if
                        endif

                        write(myunit+k,1301) date2(1), achar(9), date2(2), achar(9), date2(3),&
                            achar(9), date2(4), achar(9), date2(5), achar(9), date2(6),&
                            achar(9), hi_box(k)*100, achar(9), real(tot_ice_cells(k))/real(tot_cells(k))

                        if ( index(ln,'S-1').eq.0 ) then
                           myunit=400  
                        else
                            if (.not.divid) then
                                myunit=400  
                            else
                                myunit=600  
                            end if
                        endif
                        
                        write(myunit+k,1401) date2(1), achar(9), date2(2), achar(9), date2(3),&
                            achar(9), date2(4), achar(9), date2(5), achar(9), date2(6),&
                            ( achar(9), thick_bins(k,m), m=1,thmax )
!                    else
!                        write(*,*) trim(sat_file_name), box_index(k)
!                        write(*,*) cnt(k)/real(tot_cells(k))
                    end if
                end do

            end do !l
            close(9)
        end if!date(2).ne.month
 
        write(filen,1100) date(1), date(2), date(3)
        INQUIRE(FILE=trim(input_path)//trim(filen),EXIST=file_exists)
        if (file_exists) then
            call check( nf90_open( trim(input_path)//trim(filen), NF90_NOWRITE, ncid),1)
            call check(nf90_inq_varid(ncid, 'vicen' , vid),6)
            call check(nf90_get_var(ncid,vid,hi_tmp),12)
            call check(nf90_inq_varid(ncid, 'vlvl' , vid),6)
            call check(nf90_get_var(ncid,vid,hi_lvl),12)
            call check(nf90_close(ncid),23)
            ithick(:,:) = hi_tmp(:,:,1)+hi_tmp(:,:,2)+hi_tmp(:,:,3)+&
                hi_tmp(:,:,4)+hi_tmp(:,:,5)

            tot_cells = 0

            hi_box = 0
            tot_ice_cells = 0
            ice_cells = 0

            hi_lvl_box = 0
            tot_lvl_ice_cells = 0 
            lvl_ice_cells = 0

            thick_bins = 0
            lvl_thick_bins = 0

            do i = 1, imt
                do j = 1, jmt
                    do k = 1, nmax

                        if ( rect_mask(i,j).eq.k .and. region_mask(i,j) .ne. 0 ) then
                            tot_cells(k) = tot_cells(k) + 1
                            hi_box(k) = hi_box(k)+ithick(i,j)
                            hi_lvl_box(k) = hi_lvl_box(k)+hi_lvl(i,j)

                            if ( ithick(i,j)*100.gt.treshold ) then
                                tot_ice_cells(k) = tot_ice_cells(k)+1
                            end if
                            do m = 1, thmax-1
                                if ( ithick(i,j)*100.gt.ice_thick_dist(m) .and.&
                                    ithick(i,j)*100.le.ice_thick_dist(m+1) ) then
                                    ice_cells(k,m) = ice_cells(k,m)+1         
                                end if
                            end do
                            if ( ithick(i,j)*100.gt.ice_thick_dist(thmax) ) then
                                ice_cells(k,thmax) = ice_cells(k,thmax)+1         
                            end if

                            if ( hi_lvl(i,j)*100.gt.treshold ) then
                                tot_lvl_ice_cells(k) = tot_lvl_ice_cells(k)+1
                            end if
                            do m = 1, thmax-1
                                if ( hi_lvl(i,j)*100.gt.ice_thick_dist(m) .and.&
                                    hi_lvl(i,j)*100.le.ice_thick_dist(m+1) ) then
                                    lvl_ice_cells(k,m) = lvl_ice_cells(k,m)+1         
                                end if
                            end do
                            if ( ithick(i,j)*100.gt.ice_thick_dist(thmax) ) then
                                lvl_ice_cells(k,thmax) = lvl_ice_cells(k,thmax)+1         
                            end if

                        end if

                    end do!boxes
                end do!j
            end do!i
            
            do k = 1, nmax
                hi_box(k) = hi_box(k)/real(tot_cells(k))
                hi_lvl_box(k) = hi_lvl_box(k)/real(tot_cells(k))
                if (tot_ice_cells(k).gt.0) then
                    do m = 1, thmax
                        thick_bins(k,m) = real(ice_cells(k,m))/real(tot_ice_cells(k))
                        lvl_thick_bins(k,m) = real(lvl_ice_cells(k,m))/real(tot_lvl_ice_cells(k))
                    end do
                end if
            end do

            do k = 1, nmax
                myunit=100
                write(myunit+k,1300) date(1), achar(9), date(2), achar(9), date(3),&
                    achar(9), hi_box(k)*100, achar(9), real(tot_ice_cells(k))/real(tot_cells(k))
                myunit=200
                write(myunit+k,1400) date(1), achar(9), date(2), achar(9), date(3),&
                    ( achar(9), thick_bins(k,m), m=1,thmax )
                myunit=700
                write(myunit+k,1300) date(1), achar(9), date(2), achar(9), date(3),&
                    achar(9), hi_lvl_box(k)*100, achar(9),&
                    real(tot_lvl_ice_cells(k))/real(tot_cells(k))
                myunit=800
                write(myunit+k,1400) date(1), achar(9), date(2), achar(9), date(3),&
                    ( achar(9), lvl_thick_bins(k,m), m=1,thmax )
            end do
        else
            write(*,*) 'Not exists:',trim(input_path)//trim(filen)
        end if!file exists

!exit condition
        if ( date(3) .ge. end_date(3) .and.&
            date(2) .ge. end_date(2) .and.&
            date(1) .ge. end_date(1) ) exit
!date +1 day
        date(3) = date(3) + 1
        if ( date(3).gt.days_in_month(date(2)) ) then
            date(3) = 1
            date(2) = date(2) +1
            if ( date(2).gt.12 ) then
                date(2) = 1
                date(1) = date(1) + 1
            end if
        end if


    end do

    do i = 1, nmax
        myunit=100
        close(myunit+i)
        myunit=200
        close(myunit+i)
        myunit=300
        close(myunit+i)
        myunit=400
        close(myunit+i)
        if (divid) then
            myunit=500
            close(myunit+i)
            myunit=600
            close(myunit+i)
        end if
        myunit=700
        close(myunit+i)
        myunit=800
        close(myunit+i)
    end do!units
!###############################################################


    deallocate( hi_box, box_index, tot_cells, box_ex, cnt )
    deallocate( ice_cells, tot_ice_cells, thick_bins )
    deallocate( hi_lvl_box, lvl_thick_bins)
    deallocate( tot_lvl_ice_cells, lvl_ice_cells, tot_cells2)

contains
  subroutine check(status2,ann)
    integer, intent ( in) :: status2,ann

    if(status2 /= nf90_noerr) then
        write(*,*) "error at", ann
      print *, trim(nf90_strerror(status2))
      stop "Stopped"
    end if
  end subroutine check
  
  subroutine create_rectangles(rect_mask,nx,ny,nmax,box_index,region_mask,tot_cells)

    integer*4, intent (in) :: nx,ny,nmax
    integer*4, intent (in) :: region_mask(nx,ny)
    integer*4, intent(out) :: rect_mask(nx,ny),tot_cells(nmax)
    integer*4, intent(out) :: box_index(nmax)
    ! local variables
    integer*4 :: rect(5),ii,p1,p2
    integer*4 :: ti4(nx,ny),i,j,k
    
    rect_mask(:,:) = 0
    tot_cells = 0
    
    open(10,file='rect_corners.txt',status='old')
    do ii=1,nmax
       read(10,*) rect(:)
       box_index(ii) = rect(1)
       p1=rect(2)
       p2=rect(3)
       if (p1 > p2) then 
          rect(2) = p2
          rect(3) = p1
       endif
       p1=rect(4)
       p2=rect(5)
       if (p1 > p2) then 
          rect(4) = p2
          rect(5) = p1
       endif
       rect_mask(rect(2):rect(3),rect(4):rect(5)) = ii
       ti4 = 0
       where (rect_mask == ii)
         ti4 = 1
       endwhere
!       write(*,*) ii,rect(:),maxval(rect_mask),minval(rect_mask),sum(real(ti4))/nx/ny,sum(ti4)
    enddo
    close(10)

    do i = 1, nx
        do j = 1, ny
            do k = 1, nmax
                if ( rect_mask(i,j).eq.k .and. region_mask(i,j) .ne. 0 ) then
                    tot_cells(k) = tot_cells(k)+1
                elseif ( rect_mask(i,j).eq.k .and. region_mask(i,j) .eq. 0 ) then 
                    rect_mask(i,j) = 0
                end if   
            end do
        end do
    end do
  end subroutine create_rectangles


  SUBROUTINE fill_empty(imt,jmt,k,rect_mask,ithick,fillvalue)
    implicit none
    integer(kind=4), intent(in) :: imt,jmt,k
    integer(kind=4), dimension(imt,jmt), intent(in) :: rect_mask
    real(kind=8), intent(in) :: fillvalue
    real(kind=8), intent(inout) :: ithick(imt,jmt)

    integer :: i,j,cnt,mx,done,lbias(8),rbias(8),n
    integer, allocatable :: empty(:,:)
    real(kind=8) :: new

    lbias = (/1,1,0,-1,-1,-1,0,1/)
    rbias = (/0,-1,-1,-1,0,1,1,1/)

    cnt = 0
    do i = 1, imt
        do j = 1, jmt
            if ( rect_mask(i,j).eq.k .and. ithick(i,j).eq.fillvalue ) cnt=cnt+1
        end do
    end do
    allocate(empty(cnt,2))

    cnt = 0
    do i = 1, imt
        do j = 1, jmt
            if ( rect_mask(i,j).eq.k .and. ithick(i,j).eq.fillvalue ) then
                cnt=cnt+1
                empty(cnt,1) = i    
                empty(cnt,2) = j    
            end if
        end do
    end do

    mx = cnt
    cnt = 0
!    write(*,*) mx
    do j = 1, 100
        cnt = cnt+1
        if (cnt.gt.mx) cnt=1
        if (cnt.eq.1) done=0

        if (empty(cnt,1).eq.0) then
            done = done+1
        else
            n = 0
            new = 0
            do i = 1, 8
                if ( ithick(empty(cnt,1)+lbias(i),empty(cnt,2)+rbias(i)).ne.fillvalue ) then
                    n = n+1
                    new = new+ithick(empty(cnt,1)+lbias(i),empty(cnt,2)+rbias(i))
                end if
            end do
            if ( n.ge.1 ) then
!                write(*,*) new/real(n),empty(cnt,1),empty(cnt,2)
                ithick(empty(cnt,1),empty(cnt,2)) = new/real(n)
                empty(cnt,1) = 0
                empty(cnt,2) = 0
!            else
!                write(*,*) empty(cnt,1),empty(cnt,2)
            end if
        end if

        if (done.eq.mx) EXIT 
    end do
!    write(*,*) done 
    deallocate(empty)
  END SUBROUTINE fill_empty

  SUBROUTINE in_convex_polygon(n,x,corners,inside)
    implicit none
    integer, intent(in) :: n
    real(kind=8), intent(in) :: x(2),corners(n,2)
    logical, intent(out) :: inside

    real(kind=8) :: d
    logical :: positive
    integer :: i,j,iter

    d = (x(1) - corners(1,1))*(corners(2,2) - corners(1,2))-&
        (x(2) - corners(1,2))*(corners(2,1) - corners(1,1))
    if (d.ge.0) then
        positive = .true.
    else
        positive = .false.
    end if
    inside = .true.

    do iter = 1, n-1
        i = iter+1
        j = modulo(i,n)+1
        d = (x(1) - corners(i,1))*(corners(j,2) - corners(i,2))-&
            (x(2) - corners(i,2))*(corners(j,1) - corners(i,1))
        if (positive .ne. (d.ge.0)) then
            inside = .false.
            exit
        end if
    end do

  END SUBROUTINE in_convex_polygon

  SUBROUTINE create_box(nob,line_num,lon,lat,box_mask,box_index,region_mask,tot_cells)
    implicit none
    integer, intent(in) :: nob,line_num,region_mask(600,640)
    integer, intent(inout) :: box_mask(600,640)
    integer, intent(out) :: box_index(nob),tot_cells(nob)
    real(kind=8), dimension(600,640), intent(in) :: lon,lat

    integer :: no(line_num),cnt,i,j,k,nxt_cycle
    real(kind=8) :: vrt(line_num,2),x(2),three(3,2)
    real(kind=8), allocatable :: corners(:,:)
    logical :: now, inside, convex

    tot_cells = 0
 
    102 format( i2,',(',f5.2,',',f5.2,')')

    open(10,file='box_vrt.txt',status='old')
        do i = 1, line_num
            read(10,102) no(i),vrt(i,1),vrt(i,2)
        end do
    close(10)
  
    now = .false.
    cnt = 0
    do i = 1, line_num
        cnt = cnt+1
        if ( i.eq.line_num ) then
            now = .true.
        elseif ( no(i+1).ne.no(i) ) then
            now = .true.
        end if
        if (now) then
            allocate( corners(cnt,2) )
            do j = 1, cnt
                corners(j,:) = vrt(i+j-cnt,:)
            end do
            do j = 1, 600
                do k = 1, 640
                    x(1) = lon(j,k)
                    x(2) = lat(j,k)
                    call in_convex_polygon(cnt,x,corners,inside)
                    if (inside .and. region_mask(j,k) .ne. 0 ) box_mask(j,k) = no(i)
                    if (box_mask(j,k).eq.20) box_mask(j,k)=10
                end do
            end do
            deallocate ( corners )
            cnt = 0
            now = .false.
        end if 
    end do

    do k = 1, nob
        box_index(k) = k
    end do

    do i = 1, 600
        do j = 1, 640
            do k = 1, nob
                if ( box_mask(i,j).eq.k .and. region_mask(i,j) .ne. 0 ) then
                    tot_cells(k) = tot_cells(k)+1
                end if
            end do
        end do
    end do
  END SUBROUTINE create_box
         
end program

 function nxt_cycle(start,step,period) result(j)
    integer, intent(in) :: start,step,period ! input
    integer             :: j ! output
    j = modulo(start-1+step,period)+1
 end function nxt_cycle
