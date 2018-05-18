
module ncaxis_mod

use ncd_define_mod
use      ncatt_mod, only:  ncatt_type, read_ncatt, write_ncatt, &
                           copy_ncatt
use     ncdata_mod, only: ncdata_type, read_ncdata, write_ncdata,  &
                          get_ncdata,  put_ncdata,  remove_ncdata, &
                          copy_ncdata

implicit none
private

public :: read_ncaxis, write_ncaxis, close_ncaxis, remove_ncaxis_data
public :: get_ncaxis, put_ncaxis, get_ncaxis_data, get_ncaxis_no
public :: copy_ncaxis, ncaxis_type

   type ncaxis_type
       character(len=maxn)     :: name
       integer                 :: datatype, len, natt, ndim,  &
                                  vdims(max_dims)
       type( ncatt_type), pointer ::  Att(:)
       type(ncdata_type)          :: Data
   end type ncaxis_type

!-----------------------------------------------------------------------

   interface put_ncaxis
      module procedure put_ncaxis_float, put_ncaxis_double
   end interface

   interface get_ncaxis_data
      module procedure get_axis_data_float, get_axis_data_double
   end interface

CONTAINS

!-----------------------------------------------------------------------

   function read_ncaxis (fid, varid, rdim) result (Axis)

   integer,        intent(in) :: fid, varid
   logical,        intent(in) :: rdim
   type(ncaxis_type)          :: Axis
   integer  i, ierr
   character(len=maxn) :: aname

     call ncdinq (fid, varid, Axis%name, Axis%len, ierr)
     call ncvinq (fid, varid, aname, Axis%datatype, Axis%ndim,  &
                  Axis%vdims, Axis%natt, ierr)

!    ---- check Axis%name and aname ??? ----

     if (Axis%natt > 0) then
         allocate (Axis%Att(Axis%natt))
         do i=1,Axis%natt
            Axis%Att(i) = read_ncatt (fid, varid, i)
         enddo
     endif

!   ---- do not read axis data for record dimension (pass len=0) ----
     if (Axis%len > 0) then
        if (rdim) then
          Axis%Data = read_ncdata (fid, varid, Axis%datatype, 0)
        else
          Axis%Data = read_ncdata (fid, varid, Axis%datatype, Axis%len)
           endif
     endif

   end function read_ncaxis

!-----------------------------------------------------------------------

   subroutine write_ncaxis (fid, Axis)

   integer,           intent(in) :: fid
   type(ncaxis_type), intent(in) :: Axis
   integer  i, varid, ierr
   integer  ncdid

      call ncddef (fid, trim(Axis%name), Axis%len, ierr)
      call ncvdef (fid, trim(Axis%name), Axis%datatype, Axis%ndim,  &
                   Axis%vdims, ierr)

      varid = ncdid (fid, trim(Axis%name), ierr)

      do i=1,Axis%natt
          call write_ncatt (fid, varid, Axis%Att(i))
      enddo

      if (Axis%len > 0) then
         call ncendf (fid, ierr)
             call write_ncdata (fid, varid, Axis%Data)
         call ncredf (fid, ierr)
      endif

   end subroutine write_ncaxis

!-----------------------------------------------------------------------

   subroutine remove_ncaxis_data (Axis)
   type(ncaxis_type), intent(inout) :: Axis

      Axis%len = 0
      call remove_ncdata (Axis%Data)

   end subroutine remove_ncaxis_data

!-----------------------------------------------------------------------

   subroutine close_ncaxis (fid, varid, Rec_axis)

   integer,           intent(in)    :: fid, varid
   type(ncaxis_type), intent(inout) :: Rec_axis
   integer   ierr
   character(len=maxn) :: name

!  ------ update length of the record dimension ------

     call ncdinq (fid, varid, name, Rec_Axis%len, ierr)

   end subroutine close_ncaxis

!-----------------------------------------------------------------------

   function get_ncaxis_no (Axis, name) result (no)

   type(ncaxis_type), intent(in), dimension(:)  :: Axis
   character(len=*),  intent(in)                :: name
   integer                                      :: no
   integer  i

     no = 0

     do i = 1, size(Axis)
       if ( trim(Axis(i)%name) == trim(name) ) then
           no = i
           exit
       endif
     enddo

   end function get_ncaxis_no

!-----------------------------------------------------------------------

   subroutine get_ncaxis (Axis, naxis, name, datatype, len, natt, id)

   type(ncaxis_type), intent(in), dimension(:) :: Axis
   integer,           intent(in)               :: naxis
   character(len=*),  intent(in)               :: name
   integer         ,  intent(out)              :: datatype, len, natt
   integer,           intent(out)              :: id
   integer  k, n, na

      na = min(size(Axis), naxis)

      id = get_ncaxis_no (Axis(1:na), name)

      if (id > 0) then
               len = Axis(id)%len
              natt = Axis(id)%natt
          datatype = Axis(id)%datatype
      endif

!del  call get_ncdata (Axis(id)%Data, data)

   end subroutine get_ncaxis

!-----------------------------------------------------------------------

   subroutine put_ncaxis_float (Axis, naxis, name, data, id)

   type(ncaxis_type), intent(inout), dimension(:) :: Axis
   integer,           intent(inout)               :: naxis
   character(len=*),  intent(in)                  :: name
   real(float_type),  intent(in)   , dimension(:) :: data
   integer,           intent(out)                 :: id
   integer    i, na, nc

      na = min(size(Axis), naxis)

      id = get_ncaxis_no (Axis(1:na), name)

      if (id == 0) then
          if (naxis+1 <= max_ndim) then
              if (naxis+1 <= size(Axis)) then
                  naxis = naxis + 1
                  id = naxis
                  allocate (Axis(id)%Att(max_natt))
              else
                  print *, 'ERROR: axis too small.'
                  stop 111
              endif
          else
              print *, 'ERROR: max_ndim too small.'
              stop 111
          endif
      endif

!------------- add/replace axis -------------

         do i=1,maxn
           Axis(id)%name(i:i) = ' '
         enddo

         nc = len_trim(name)
         Axis(id)%name(1:nc) = name(1:nc)
         Axis(id)%len        = size(data)
         Axis(id)%datatype   = NCFLOAT
         Axis(id)%natt       = 0
         Axis(id)%ndim       = 1
         Axis(id)%vdims      = (/ id, 0, 0, 0 /)
         Axis(id)%Data       = put_ncdata(data)


   end subroutine put_ncaxis_float

!-----------------------------------------------------------------------

   subroutine put_ncaxis_double (Axis, naxis, name, data, id)

   type(ncaxis_type), intent(inout), dimension(:) :: Axis
   integer,           intent(inout)               :: naxis
   character(len=*),  intent(in)                  :: name
   real(double_type), intent(in)   , dimension(:) :: data
   integer,           intent(out)                 :: id
   integer   na, nc

      na = min(size(Axis), naxis)

      id = get_ncaxis_no (Axis(1:na), name)

      if (id == 0) then
          if (naxis+1 <= max_ndim) then
              if (naxis+1 <= size(Axis)) then
                  naxis = naxis + 1
                  id = naxis
                  allocate (Axis(id)%Att(max_natt))
              else
                  print *, 'ERROR: axis too small.'
                  stop 111
              endif
          else
              print *, 'ERROR: max_ndim too small.'
              stop 111
          endif
      endif

!------------- add/replace axis -------------

!!       do i=1,maxn
!!         Axis(id)%name(i:i) = ' '
!!       enddo

         nc = len_trim(name)
!del     Axis(id)%name(1:nc) = name(1:nc)
         Axis(id)%name       = name
         Axis(id)%len        = size(data)
         Axis(id)%datatype   = NCDOUBLE
         Axis(id)%natt       = 0
         Axis(id)%ndim       = 1
         Axis(id)%vdims      = (/ id, 0, 0, 0 /)
         Axis(id)%Data       = put_ncdata(data)


   end subroutine put_ncaxis_double

!-----------------------------------------------------------------------

   subroutine get_axis_data_float (Axis, data, start, count)

   type(ncaxis_type), intent(in)                :: Axis
   real(float_type),  intent(out), dimension(:) :: data
   integer         ,  intent(in),  optional     :: start, count
   integer   start_local, count_local

      if (Axis%datatype /= NCFLOAT) then
          print *, 'ERROR in get_axis_data_float: &
                   &axis has wrong datatype.'
          stop 111
      endif

      start_local = 1;        if (present(start)) start_local = start
      count_local = Axis%len; if (present(count)) count_local = count

      call get_ncdata (Axis%Data, data, start_local, count_local)

   end subroutine get_axis_data_float

!-----------------------------------------------------------------------

   subroutine get_axis_data_double (Axis, data, start, count)

   type(ncaxis_type), intent(in)                :: Axis
   real(double_type), intent(out), dimension(:) :: data
   integer          , intent(in),  optional     :: start, count
   integer   start_local, count_local

      if (Axis%datatype /= NCDOUBLE) then
          print *, 'ERROR in get_axis_data_float: &
                   &axis has wrong datatype.'
          stop 111
      endif

      start_local = 1;        if (present(start)) start_local = start
      count_local = Axis%len; if (present(count)) count_local = count

      call get_ncdata (Axis%Data, data, start_local, count_local)

   end subroutine get_axis_data_double

!-----------------------------------------------------------------------

   subroutine copy_ncaxis (Axis_in, Axis_out)

   type(ncaxis_type), intent(in)  :: Axis_in
   type(ncaxis_type), intent(out) :: Axis_out
   integer :: i

      Axis_out%name     = Axis_in%name
      Axis_out%datatype = Axis_in%datatype
      Axis_out%len      = Axis_in%len
      Axis_out%ndim     = Axis_in%ndim
      Axis_out%vdims    = Axis_in%vdims
      Axis_out%natt     = Axis_in%natt

      allocate ( Axis_out%Att( Axis_out%natt ) )
      do i = 1, Axis_out%natt
           call copy_ncatt ( Axis_in%Att(i), Axis_out%Att(i) )
      enddo

           call copy_ncdata ( Axis_in%Data, Axis_out%Data )

   end subroutine copy_ncaxis

!-----------------------------------------------------------------------


end module ncaxis_mod

