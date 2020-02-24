
module ncread_write_mod

use  ncd_define_mod
use      ncfile_mod, only: ncfile_type, switch_define_mode
use    ncvarmap_mod, only: get_varmap, get_varmap_tavg
use      ncaxis_mod, only: get_ncaxis_data, get_ncaxis_no

implicit none
private

public  write_time_coord
public  write_variable_tave
public  write_variable
public  read_time_coord
public  read_variable
public  read_variable_tave
public  axis_values
public  time_axis_values
public  variable_tave_info

!-----------------------------------------------------------------------

interface write_time_coord
   module procedure  write_time_coord_float_0d,   &
                     write_time_coord_float_1d,   &
                     write_time_coord_double_0d,  &
                     write_time_coord_double_1d
end interface

!-----------------------------------------------------------------------

interface write_variable_tave
   module procedure  write_variable_tave_float_0d,   &
                     write_variable_tave_float_1d,   &
                     write_variable_tave_double_0d,  &
                     write_variable_tave_double_1d
end interface

!-----------------------------------------------------------------------

interface write_variable
   module procedure  write_variable_1d_count,               &
                     write_variable_1d, write_variable_2d,  &
                     write_variable_3d, write_variable_4d
end interface

!-----------------------------------------------------------------------

interface read_time_coord
   module procedure  read_time_coord_float_0d,   &
                     read_time_coord_float_1d,   &
                     read_time_coord_double_0d,  &
                     read_time_coord_double_1d
end interface

!-----------------------------------------------------------------------

interface read_variable
   module procedure  read_variable_1d_count,              &
                     read_variable_1d, read_variable_2d,  &
                     read_variable_3d, read_variable_4d
end interface

!-----------------------------------------------------------------------

interface read_variable_tave
   module procedure  read_variable_tave_float_0d,   &
                     read_variable_tave_float_1d,   &
                     read_variable_tave_double_0d,  &
                     read_variable_tave_double_1d
end interface

!-----------------------------------------------------------------------

interface axis_values
   module procedure  axis_values_num, axis_values_str
end interface

!-----------------------------------------------------------------------

interface time_axis_values
   module procedure  time_axis_values_float, time_axis_values_double
end interface

!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------

  subroutine write_time_coord_float_0d (file, start, data)
    type(ncfile_type), intent(inout) :: file
    integer          , intent(in)    :: start
    real(float_type) , intent(in)    :: data

    real(double_type) :: ddata
    integer    ierr

    if (file%define_mode) call switch_define_mode (file)

    ddata = data
    call ncvpt1 (file%id, file%recdim, start, ddata, ierr)

  end subroutine write_time_coord_float_0d

!-----------------------------------------------------------------------

  subroutine write_time_coord_double_0d (file, start, data)
    type(ncfile_type), intent(inout) :: file
    integer          , intent(in)    :: start
    real(double_type), intent(in)    :: data

    integer    ierr

    if (file%define_mode) call switch_define_mode (file)

    call ncvpt1 (file%id, file%recdim, start, data, ierr)

  end subroutine write_time_coord_double_0d

!-----------------------------------------------------------------------

  subroutine write_time_coord_float_1d (file, start, data)
    type(ncfile_type), intent(inout) :: file
    integer          , intent(in)    :: start
    real(float_type) , intent(in)    :: data(:)

    real(double_type), dimension(size(data)) :: ddata
    integer    ierr

    if (file%define_mode) call switch_define_mode (file)

    ddata = data
    call ncvpt (file%id, file%recdim, start, size(data), ddata, ierr)

  end subroutine write_time_coord_float_1d

!-----------------------------------------------------------------------

  subroutine write_time_coord_double_1d (file, start, data)
    type(ncfile_type), intent(inout) :: file
    integer          , intent(in)    :: start
    real(double_type), intent(in)    :: data(:)

    integer    ierr

    if (file%define_mode) call switch_define_mode (file)

    call ncvpt (file%id, file%recdim, start, size(data), data, ierr)

  end subroutine write_time_coord_double_1d

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------

 subroutine write_variable_1d_count (file, var_name, start, count, data)

  type(ncfile_type), intent(inout) :: file
  character(len=*) , intent(in)    :: var_name
  integer          , intent(in)    :: start(:), count(:)
  real             , intent(in)    :: data(:)

  integer   varid
  logical   pack
  real      scale, offset

     call variable_info (file, var_name, varid, pack, scale, offset)

!   ----- add offset, multiply by scale factor (if necessary) -----

     if (file%define_mode) call switch_define_mode (file)

     if (pack) then
        call varput_short_1d (file%id, varid, start, count,  &
                          scale, offset, data)
     else
        call varput_float_1d (file%id, varid, start, count, data)
     endif

 end subroutine write_variable_1d_count

!-----------------------------------------------------------------------

  subroutine write_variable_1d (file, var_name, start, data)

  type(ncfile_type), intent(inout) :: file
  character(len=*) , intent(in)    :: var_name
  integer          , intent(in)    :: start(:)
  real             , intent(in)    :: data(:)

  integer :: varid, kount(size(start))
  logical :: pack
  real    :: scale, offset

     call variable_info (file, var_name, varid, pack, scale, offset)

!   ----- add offset, multiply by scale factor (if necessary) -----

     if (file%define_mode) call switch_define_mode (file)

     kount = 1
     kount(1) = size(data,1)
     
     if (pack) then
        call varput_short_1d (file%id, varid, start, kount,  &
                          scale, offset, data)
     else
        call varput_float_1d (file%id, varid, start, kount, data)
     endif

  end subroutine write_variable_1d

!-----------------------------------------------------------------------

  subroutine write_variable_2d (file, var_name, start, data)

  type(ncfile_type), intent(inout) :: file
  character(len=*) , intent(in)    :: var_name
  integer          , intent(in)    :: start(:)
  real             , intent(in)    :: data(:,:)

  integer :: varid, len(1), kount(size(start))
  logical :: pack
  real    :: scale, offset

     if (size(start) < 2)  call error_handler ('write_variable_2d', &
                                         'invalid size for start', 2)

     call variable_info (file, var_name, varid, pack, scale, offset)

!   ----- add offset, multiply by scale factor (if necessary) -----

     if (file%define_mode) call switch_define_mode (file)

     kount = 1
     kount(1) = size(data,1)
     kount(2) = size(data,2)
     len = (/ size(data,1)*size(data,2) /)

     if (pack) then
        call varput_short_1d (file%id, varid, start, kount,    &
                              scale, offset, reshape(data,len) )
     else
        call varput_float_1d (file%id, varid, start, kount, &
                              reshape(data,len) )
     endif

  end subroutine write_variable_2d

!-----------------------------------------------------------------------

  subroutine write_variable_3d (file, var_name, start, data)

  type(ncfile_type), intent(inout) :: file
  character(len=*) , intent(in)    :: var_name
  integer          , intent(in)    :: start(:)
  real             , intent(out)   :: data(:,:,:)

  integer :: varid, len(1), kount(size(start))
  logical :: pack
  real    :: scale, offset

     if (size(start) < 3)  call error_handler ('write_variable_3d', &
                                         'invalid size for start', 2)

     call variable_info (file, var_name, varid, pack, scale, offset)

!   ----- add offset, multiply by scale factor (if necessary) -----

     if (file%define_mode) call switch_define_mode (file)

     kount = 1
     kount(1) = size(data,1)
     kount(2) = size(data,2)
     kount(3) = size(data,3)
     len = (/ size(data,1)*size(data,2)*size(data,3) /)

     if (pack) then
        call varput_short_1d (file%id, varid, start, kount,  &
                              scale, offset, reshape(data,len) )
     else
        call varput_float_1d (file%id, varid, start, kount,  &
                              reshape(data,len) )
     endif

  end subroutine write_variable_3d

!-----------------------------------------------------------------------

  subroutine write_variable_4d (file, var_name, start, data)

  type(ncfile_type), intent(inout) :: file
  character(len=*) , intent(in)    :: var_name
  integer          , intent(in)    :: start(:)
  real             , intent(out)   :: data(:,:,:,:)

  integer :: varid, len(1), kount(size(start))
  logical :: pack
  real    :: scale, offset

     if (size(start) < 4)  call error_handler ('write_variable_4d', &
                                         'invalid size for start', 2)

     call variable_info (file, var_name, varid, pack, scale, offset)

!   ----- add offset, multiply by scale factor (if necessary) -----

     if (file%define_mode) call switch_define_mode (file)

     kount = 1
     kount(1) = size(data,1)
     kount(2) = size(data,2)
     kount(3) = size(data,3)
     kount(4) = size(data,4)
     len = (/ size(data,1)*size(data,2)*size(data,3)*size(data,4) /)

     if (pack) then
        call varput_short_1d (file%id, varid, start, kount,  &
                              scale, offset, reshape(data,len) )
     else
        call varput_float_1d (file%id, varid, start, kount,  &
                              reshape(data,len) )
     endif

  end subroutine write_variable_4d

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------

subroutine write_variable_tave_float_0d (file, var_name, start,  &
                                         t1, t2, dt)

  type(ncfile_type), intent(inout) :: file
  character(len=*) , intent(in)    :: var_name
  integer          , intent(in)    :: start
  real(float_type) , intent(in)    :: t1, t2, dt
! integer          , intent(in)    :: nitems
  real(double_type)  :: time
  integer            :: ierr, tid, tvarno(3)

      call get_varmap_tavg (file%varmap, file%nvarmap, var_name, tvarno)

      if (file%define_mode) call switch_define_mode (file)

      tid = file%var(tvarno(1))%id
      time = t1
      call ncvpt1 (file%id, tid, start, time, ierr)

      tid = file%var(tvarno(2))%id
      time = t2
      call ncvpt1 (file%id, tid, start, time, ierr)

      tid = file%var(tvarno(3))%id
      time = dt
      call ncvpt1 (file%id, tid, start, time, ierr)

end subroutine write_variable_tave_float_0d

!-----------------------------------------------------------------------

subroutine write_variable_tave_double_0d (file, var_name, start,   &
                                          t1, t2, dt)

  type(ncfile_type), intent(inout) :: file
  character(len=*) , intent(in)    :: var_name
  integer          , intent(in)    :: start
  real(double_type), intent(in)    :: t1, t2, dt
! integer          , intent(in)    :: nitems
  integer            :: ierr, tid, tvarno(3)

      call get_varmap_tavg (file%varmap, file%nvarmap, var_name, tvarno)

      if (file%define_mode) call switch_define_mode (file)

      tid = file%var(tvarno(1))%id
      call ncvpt1 (file%id, tid, start, t1, ierr)

      tid = file%var(tvarno(2))%id
      call ncvpt1 (file%id, tid, start, t2, ierr)

      tid = file%var(tvarno(3))%id
      call ncvpt1 (file%id, tid, start, dt, ierr)

end subroutine write_variable_tave_double_0d

!-----------------------------------------------------------------------

  subroutine write_variable_tave_float_1d (file, var_name, start,  &
                                           t1, t2, dt)

  type(ncfile_type), intent(inout) :: file
  character(len=*) , intent(in)    :: var_name
  integer          , intent(in)    :: start
  real(float_type) , intent(in)    :: t1(:), t2(:), dt(:)
! integer          , intent(in)    :: nitems(:)
  real(double_type), dimension(size(t1))  :: time
  integer            :: ierr, tid, tvarno(3)

      call get_varmap_tavg (file%varmap, file%nvarmap, var_name, tvarno)

      if (file%define_mode) call switch_define_mode (file)

      tid = file%var(tvarno(1))%id
      time = t1
      call ncvpt (file%id, tid, start, size(t1), time, ierr)

      tid = file%var(tvarno(2))%id
      time = t2
      call ncvpt (file%id, tid, start, size(t2), time, ierr)

      tid = file%var(tvarno(3))%id
      time = dt
      call ncvpt (file%id, tid, start, size(dt), time, ierr)

  end subroutine write_variable_tave_float_1d

!-----------------------------------------------------------------------

  subroutine write_variable_tave_double_1d (file, var_name, start,  &
                                            t1, t2, dt)

  type(ncfile_type), intent(inout) :: file
  character(len=*) , intent(in)    :: var_name
  integer          , intent(in)    :: start
  real(double_type), intent(in)    :: t1(:), t2(:), dt(:)
! integer          , intent(in)    :: nitems(:)
  integer            :: ierr, tid, tvarno(3)

      call get_varmap_tavg (file%varmap, file%nvarmap, var_name, tvarno)

      if (file%define_mode) call switch_define_mode (file)

      tid = file%var(tvarno(1))%id
      call ncvpt (file%id, tid, start, size(t1), t1, ierr)

      tid = file%var(tvarno(2))%id
      call ncvpt (file%id, tid, start, size(t2), t2, ierr)

      tid = file%var(tvarno(3))%id
      call ncvpt (file%id, tid, start, size(dt), dt, ierr)

  end subroutine write_variable_tave_double_1d

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------

  subroutine varput_short_1d (fid, vid, start, len, scale, offset, data)
  integer, intent(in) :: fid, vid, start(:), len(:)
  real   , intent(in) :: scale, offset
  real   , intent(in) :: data(:)
  integer(short_type), dimension(size(data,1)) :: sdata
  real               , dimension(size(data,1)) :: pdata
  integer  ierr

    if (scale == 1.0 .and. offset == 0.0) then
        sdata = int(data + sign(0.5,data))
    else
        pdata = (data-offset)/scale
        sdata = int(pdata + sign(0.5,pdata))
    endif

    call ncvpt (fid, vid, start, len, sdata, ierr)

  end subroutine varput_short_1d

!-----------------------------------------------------------------------

  subroutine varput_float_1d (fid, vid, start, len, data)
  integer, intent(in) :: fid, vid, start(:), len(:)
  real   , intent(in) :: data(:)
  real(float_type), dimension(size(data,1)) :: fdata
  integer  ierr

    fdata = data
    call ncvpt (fid, vid, start, len, fdata, ierr)

  end subroutine varput_float_1d

!-----------------------------------------------------------------------
!#################### READ ROUTINES ####################################
!-----------------------------------------------------------------------

  subroutine read_time_coord_float_0d (file, start, data)
    type(ncfile_type), intent(in)  :: file
    integer          , intent(in)  :: start
    real(float_type) , intent(out) :: data
    real(double_type)                 ddata
    integer    ierr

    call ncvgt1 (file%id, file%recdim, start, ddata, ierr)
    data = ddata

  end subroutine read_time_coord_float_0d

!-----------------------------------------------------------------------

  subroutine read_time_coord_double_0d (file, start, data)
    type(ncfile_type), intent(in)  :: file
    integer          , intent(in)  :: start
    real(double_type), intent(out) :: data
    integer    ierr

    call ncvgt1 (file%id, file%recdim, start, data, ierr)

  end subroutine read_time_coord_double_0d

!-----------------------------------------------------------------------

  subroutine read_time_coord_float_1d (file, start, data)
    type(ncfile_type), intent(in)  :: file
    integer          , intent(in)  :: start
    real(float_type) , intent(out) :: data(:)
    real(double_type), dimension(size(data)) :: ddata
    integer    ierr

    call ncvgt (file%id, file%recdim, start, size(data), ddata, ierr)
    data = ddata

  end subroutine read_time_coord_float_1d

!-----------------------------------------------------------------------

  subroutine read_time_coord_double_1d (file, start, data)
    type(ncfile_type), intent(in)  :: file
    integer          , intent(in)  :: start
    real(double_type), intent(out) :: data(:)
    integer    ierr

    call ncvgt (file%id, file%recdim, start, size(data), data, ierr)

  end subroutine read_time_coord_double_1d

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------

  subroutine read_variable_1d_count (file, var_name, start, count, data)

  type(ncfile_type), intent(in)  :: file
  character(len=*) , intent(in)  :: var_name
  integer          , intent(in)  :: start(:), count(:)
  real             , intent(out) :: data(:)

  integer :: varid
  logical :: pack
  real    :: scale, offset

     call variable_info (file, var_name, varid, pack, scale, offset)

!   ----- add offset, multiply by scale factor (if necessary) -----

     if (pack) then
        call varget_short_1d (file%id, varid, start, count,  &
                              scale, offset, data )
     else
        call varget_float_1d (file%id, varid, start, count, data )
     endif

  end subroutine read_variable_1d_count

!-----------------------------------------------------------------------

  subroutine read_variable_1d (file, var_name, start, data)

  type(ncfile_type), intent(in)  :: file
  character(len=*) , intent(in)  :: var_name
  integer          , intent(in)  :: start(:)
  real             , intent(out) :: data(:)

  integer :: varid, kount(size(start))
  logical :: pack
  real    :: scale, offset

     call variable_info (file, var_name, varid, pack, scale, offset)

!   ----- add offset, multiply by scale factor (if necessary) -----

     kount = 1
     kount(1) = size(data,1)

     if (pack) then
        call varget_short_1d (file%id, varid, start, kount,  &
                              scale, offset, data )
     else
        call varget_float_1d (file%id, varid, start, kount, data )
     endif

  end subroutine read_variable_1d

!-----------------------------------------------------------------------

  subroutine read_variable_2d (file, var_name, start, data)

  type(ncfile_type), intent(in)  :: file
  character(len=*) , intent(in)  :: var_name
  integer          , intent(in)  :: start(:)
  real             , intent(out) :: data(:,:)

  integer :: varid, len(2), kount(size(start))
  logical :: pack
  real    :: scale, offset
  real, dimension(size(data,1)*size(data,2)) :: rdata

     if (size(start) < 2)  call error_handler ('read_variable_2d', &
                                         'invalid size for start', 2)

     call variable_info (file, var_name, varid, pack, scale, offset)

!   ----- add offset, multiply by scale factor (if necessary) -----

     kount = 1
     kount(1) = size(data,1)
     kount(2) = size(data,2)
     len = (/ size(data,1),size(data,2) /)

     if (pack) then
        call varget_short_1d (file%id, varid, start, kount,  &
                              scale, offset, rdata)
     else
        call varget_float_1d (file%id, varid, start, kount, rdata)
     endif

     data = reshape (rdata,len)

  end subroutine read_variable_2d

!-----------------------------------------------------------------------

  subroutine read_variable_3d (file, var_name, start, data)

  type(ncfile_type), intent(in)  :: file
  character(len=*) , intent(in)  :: var_name
  integer          , intent(in)  :: start(:)
  real             , intent(out) :: data(:,:,:)

  integer :: varid, len(3), kount(size(start))
  logical :: pack
  real    :: scale, offset
  real, dimension(size(data,1)*size(data,2)*size(data,3)) :: rdata

     if (size(start) < 3)  call error_handler ('read_variable_3d', &
                                         'invalid size for start', 2)

     call variable_info (file, var_name, varid, pack, scale, offset)

!   ----- add offset, multiply by scale factor (if necessary) -----

     kount = 1
     kount(1) = size(data,1)
     kount(2) = size(data,2)
     kount(3) = size(data,3)
     len = (/ size(data,1),size(data,2),size(data,3) /)

     if (pack) then
        call varget_short_1d (file%id, varid, start, kount,  &
                              scale, offset, rdata)
     else
        call varget_float_1d (file%id, varid, start, kount, rdata)
     endif

     data = reshape (rdata,len)

  end subroutine read_variable_3d

!-----------------------------------------------------------------------

  subroutine read_variable_4d (file, var_name, start, data)

  type(ncfile_type), intent(in)  :: file
  character(len=*) , intent(in)  :: var_name
  integer          , intent(in)  :: start(:)
  real             , intent(out) :: data(:,:,:,:)

  integer :: varid, len(4), kount(size(start))
  logical :: pack
  real    :: scale, offset
  real, dimension(size(data,1)*size(data,2)*  &
                  size(data,3)*size(data,4)) :: rdata

     if (size(start) < 4)  call error_handler ('read_variable_4d', &
                                         'invalid size for start', 2)

     call variable_info (file, var_name, varid, pack, scale, offset)

!   ----- add offset, multiply by scale factor (if necessary) -----

     kount = 1
     kount(1) = size(data,1)
     kount(2) = size(data,2)
     kount(3) = size(data,3)
     kount(4) = size(data,4)
     len = (/ size(data,1),size(data,2),size(data,3),size(data,4) /)

     if (pack) then
        call varget_short_1d (file%id, varid, start, kount,  &
                              scale, offset, rdata)
     else
        call varget_float_1d (file%id, varid, start, kount, rdata)
     endif

     data = reshape (rdata,len)

  end subroutine read_variable_4d

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------

  subroutine varget_short_1d (fid, vid, start, len, scale, offset, data)
  integer, intent(in)  :: fid, vid, start(:), len(:)
  real   , intent(in)  :: scale, offset
  real   , intent(out) :: data(:)
  integer(short_type), dimension(size(data,1)) :: sdata
  integer  ierr

    call ncvgt (fid, vid, start, len, sdata, ierr)
    if (scale == 1.0 .and. offset == 0.0) then
       data = real(sdata)
    else
       data = real(sdata)*scale + offset
    endif

  end subroutine varget_short_1d

!-----------------------------------------------------------------------

  subroutine varget_float_1d (fid, vid, start, len, data)
  integer, intent(in)  :: fid, vid, start(:), len(:)
  real   , intent(out) :: data(:)
  real(float_type), dimension(size(data,1)) :: fdata
  integer  ierr

    call ncvgt (fid, vid, start, len, fdata, ierr)
    data = fdata

  end subroutine varget_float_1d

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------

  function read_variable_tave_float_0d (file, var_name, start,  &
                                        t1, t2, dt, nitems)     &
                                result (flag)

  type(ncfile_type), intent(in)  :: file
  character(len=*) , intent(in)  :: var_name
  integer          , intent(in)  :: start
  real(float_type) , intent(out) :: t1, t2, dt
  integer, optional, intent(out) :: nitems
  logical                       :: flag
  real(double_type)             :: time
  integer            :: ierr, tid, tvarno(3)

     flag = variable_tave_info (file, var_name, tvarno)

     if (.not.flag) then
         t1=0.; t2=0.; dt=0.
         if (present(nitems)) nitems=0
         return
     endif

     tid = file%var(tvarno(1))%id
     call ncvgt1 (file%id, tid, start, time, ierr)
     t1 = time

     tid = file%var(tvarno(2))%id
     call ncvgt1 (file%id, tid, start, time, ierr)
     t2 = time

     tid = file%var(tvarno(3))%id
     if (present(nitems)) then
       call ncvgt1 (file%id, tid, start, nitems, ierr)
       dt=0.
     else
       call ncvgt1 (file%id, tid, start, time, ierr)
       dt = time
     endif

  end function read_variable_tave_float_0d

!-----------------------------------------------------------------------

  function read_variable_tave_double_0d (file, var_name, start,  &
                                         t1, t2, dt, nitems)     &
                                 result (flag)

  type(ncfile_type), intent(in)  :: file
  character(len=*) , intent(in)  :: var_name
  integer          , intent(in)  :: start
  real(double_type), intent(out) :: t1, t2, dt
  integer, optional, intent(out) :: nitems
  logical                        :: flag
  integer            :: ierr, tid, tvarno(3)

     flag = variable_tave_info (file, var_name, tvarno)

     if (.not.flag) then
         t1=0.; t2=0.; dt=0.
         if (present(nitems)) nitems=0
         return
     endif

     tid = file%var(tvarno(1))%id
     call ncvgt1 (file%id, tid, start, t1, ierr)

     tid = file%var(tvarno(2))%id
     call ncvgt1 (file%id, tid, start, t2, ierr)

     tid = file%var(tvarno(3))%id
     if (present(nitems)) then
       call ncvgt1 (file%id, tid, start, nitems, ierr)
       dt=0.
     else
       call ncvgt1 (file%id, tid, start, dt, ierr)
     endif

  end function read_variable_tave_double_0d

!-----------------------------------------------------------------------

  function read_variable_tave_float_1d (file, var_name, start,  &
                                        t1, t2, dt, nitems)     &
                                result (flag)

  type(ncfile_type), intent(in)  :: file
  character(len=*) , intent(in)  :: var_name
  integer          , intent(in)  :: start
  real(float_type) , intent(out) :: t1(:), t2(:), dt(:)
  integer, optional, intent(out) :: nitems(:)
  logical                        :: flag
  real(double_type), dimension(size(t1))  :: time
  integer            :: ierr, tid, tvarno(3)

     flag = variable_tave_info (file, var_name, tvarno)

     if (.not.flag) then
         t1=0.; t2=0.; dt=0.
         if (present(nitems)) nitems=0
         return
     endif

     tid = file%var(tvarno(1))%id
     call ncvgt (file%id, tid, start, size(t1), time, ierr)
     t1 = time

     tid = file%var(tvarno(2))%id
     call ncvgt (file%id, tid, start, size(t1), time, ierr)
     t2 = time

     tid = file%var(tvarno(3))%id
     if (present(nitems)) then
       call ncvgt (file%id, tid, start, size(t1), nitems, ierr)
       dt=0.
     else
       call ncvgt (file%id, tid, start, size(t1), dt, ierr)
     endif

  end function read_variable_tave_float_1d

!-----------------------------------------------------------------------

 function read_variable_tave_double_1d (file, var_name, start,  &
                                        t1, t2, dt, nitems)     &
                                result (flag)

  type(ncfile_type), intent(in)  :: file
  character(len=*),  intent(in)  :: var_name
  integer         ,  intent(in)  :: start
  real(double_type), intent(out) :: t1(:), t2(:), dt(:)
  integer, optional, intent(out) :: nitems(:)
  logical                        :: flag
  integer            :: ierr, tid, tvarno(3)

     flag = variable_tave_info (file, var_name, tvarno)

     if (.not.flag) then
         t1=0.; t2=0.; dt=0.
         if (present(nitems)) nitems=0
         return
     endif

     tid = file%var(tvarno(1))%id
     call ncvgt (file%id, tid, start, size(t1), t1, ierr)

     tid = file%var(tvarno(2))%id
     call ncvgt (file%id, tid, start, size(t1), t2, ierr)

     tid = file%var(tvarno(3))%id
     if (present(nitems)) then
       call ncvgt (file%id, tid, start, size(t1), nitems, ierr)
       dt=0.
     else
       call ncvgt (file%id, tid, start, size(t1), dt, ierr)
     endif

 end function read_variable_tave_double_1d

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------

  subroutine variable_info (file, var_name, varid, pack, scale, offset)

     type(ncfile_type), intent(in)  :: file
     character(len=*) , intent(in)  :: var_name
     integer          , intent(out) :: varid
     logical          , intent(out) :: pack
     real             , intent(out) :: scale, offset
     integer   varno

     call get_varmap (file%varmap, file%nvarmap, var_name,  &
                      varno, pack, scale, offset)
     varid = file%var(varno)%id

   end subroutine variable_info

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------

  function variable_tave_info (file, var_name, tvars) result (flag)

     type(ncfile_type), intent(in)  :: file
     character(len=*) , intent(in)  :: var_name
     integer, optional, intent(out) :: tvars(3)
     integer                        :: tvarno(3)
     logical                        :: flag

      call get_varmap_tavg (file%varmap, file%nvarmap, var_name, tvarno)

      if (tvarno(1) == 0 .or. tvarno(2) == 0 .or. tvarno(3) == 0) then
         flag = .false.
      else
         flag = .true.
      endif

      if (present(tvars)) tvars = tvarno

   end function variable_tave_info

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------

  subroutine axis_values_num (file, num, data, start)

     type(ncfile_type), intent(in)  :: file
     integer          , intent(in)  :: num
     real             , intent(out) :: data(:)
     integer, optional, intent(in)  :: start
!    real(float_type),  dimension(size(data)) :: fdata
     real(double_type), dimension(size(data)) :: fdata
     integer :: len, begin

       begin = 1
       if (present(start)) begin = start

       len = size(data)
       if ( len+begin-1 > file%axis(num)%len ) call error_handler &
            ('axis_values_num', 'bad data size or start', 2)

       call get_ncaxis_data (file%axis(num), fdata, begin, len)

       data = fdata

  end subroutine axis_values_num

!-----------------------------------------------------------------------

  subroutine axis_values_str (file, axis_name, data, start)

     type(ncfile_type), intent(in)  :: file
     character(len=*) , intent(in)  :: axis_name
     real             , intent(out) :: data(:)
     integer, optional, intent(in)  :: start
!    real(float_type),  dimension(size(data)) :: fdata
     real(double_type), dimension(size(data)) :: fdata
     integer :: num

       num = get_ncaxis_no (file%axis(1:file%ndim), axis_name)

       call axis_values_num (file, num, data, start)

  end subroutine axis_values_str

!-----------------------------------------------------------------------

  subroutine time_axis_values_float (file, data, start)

     type(ncfile_type), intent(in)  :: file
     real(float_type) , intent(out) :: data(:)
     integer, optional, intent(in)  :: start
     real(double_type), dimension(size(data)) :: ddata
     integer  num, len, begin

       begin = 1
       if (present(start)) begin = start

       num = file%recdim
       len = size(data)
       if ( len+begin-1 > file%axis(num)%len ) call error_handler &
            ('time_axis_values_float', 'bad data size or start', 2)

       call read_time_coord (file, begin, data)

  end subroutine time_axis_values_float

!-----------------------------------------------------------------------

  subroutine time_axis_values_double (file, data, start)

     type(ncfile_type), intent(in)  :: file
     real(double_type), intent(out) :: data(:)
     integer, optional, intent(in)  :: start
     integer  num, len, begin

       begin = 1
       if (present(start)) begin = start

       num = file%recdim
       len = size(data)
       if ( len+begin-1 > file%axis(num)%len ) call error_handler &
            ('time_axis_values_double', 'bad data size or start', 2)

       call read_time_coord (file, begin, data)

  end subroutine time_axis_values_double

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------

  subroutine error_handler (routine, message, level)

    character(len=*), intent(in) :: routine, message
    integer         , intent(in) :: level

       select case (level)
          case (0)
              print *, 'NOTE from ' // trim(routine) // ':'
              print *,                 trim(message)
          case (1)
              print *, 'WARNING from ' // trim(routine) // ':'
              print *,                    trim(message)
          case default
              print *, 'FATAL ERROR from ' // trim(routine) // ':'
              print *,                        trim(message)
       end select

  end subroutine error_handler

!-----------------------------------------------------------------------
!#######################################################################

end module ncread_write_mod

