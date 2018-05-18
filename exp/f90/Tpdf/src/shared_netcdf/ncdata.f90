
module ncdata_mod

use ncd_define_mod

implicit none
private

public :: read_ncdata, write_ncdata, get_ncdata, put_ncdata,  &
          remove_ncdata, copy_ncdata, ncdata_type

   type ncdata_type
      integer             :: datatype, len
      character(len=1)   , pointer :: bdata(:)
      character(len=maxc), pointer :: cdata(:)
      integer(short_type), pointer :: sdata(:)
      integer( long_type), pointer :: ldata(:)
      real   (float_type), pointer :: fdata(:)
      real  (double_type), pointer :: ddata(:)
   end type ncdata_type

!-----------------------------------------------------------------------

   interface get_ncdata
      module procedure get_data_float, get_data_double
   end interface

   interface put_ncdata
      module procedure put_data_float, put_data_double
   end interface

   integer, parameter :: FATAL = 2

CONTAINS

!-----------------------------------------------------------------------

 subroutine alloc_ncdata (type, len, Data)

 integer,           intent(in)  :: type, len
 type(ncdata_type), intent(out) :: Data

!  allocates space in appropriate data type
!  assumes that Data has not been allocated

   select case (type)
      case (NCBYTE)
         allocate (Data%bdata(len))
      case (NCCHAR)
         allocate (Data%cdata(len))
      case (NCSHORT)
         allocate (Data%sdata(len))
      case (NCLONG)
         allocate (Data%ldata(len))
      case (NCFLOAT)
         allocate (Data%fdata(len))
      case (NCDOUBLE)
         allocate (Data%ddata(len))
      case default
         call error_handler ('alloc_ncdata', 'invalid data type', FATAL)
   end select

   Data%datatype = type
   Data%len      = len

 end subroutine alloc_ncdata

!-----------------------------------------------------------------------

   function read_ncdata (id, varid, type, len) result (Data)

   integer,        intent(in)  :: id, varid, type, len
   type(ncdata_type)           :: Data
   integer  n, ierr
   integer :: start = 1

     call alloc_ncdata (type, len, Data)

     if (Data%len == 0) return

     select case (Data%datatype)
        case (NCBYTE)
           call ncvgtc (id, varid, start, Data%len, Data%bdata,  &
                        Data%len, ierr)
        case (NCCHAR)
           n = Data%len * maxc
           call ncvgtc (id, varid, start, Data%len, Data%cdata, n,  &
                        ierr)
        case (NCSHORT)
           call ncvgt  (id, varid, start, Data%len, Data%sdata, ierr)
        case (NCLONG)
           call ncvgt  (id, varid, start, Data%len, Data%ldata, ierr)
        case (NCFLOAT)
           call ncvgt  (id, varid, start, Data%len, Data%fdata, ierr)
        case (NCDOUBLE)
           call ncvgt  (id, varid, start, Data%len, Data%ddata, ierr)
     end select

   end function read_ncdata

!-----------------------------------------------------------------------

   subroutine write_ncdata (id, varid, data)

   integer,           intent(in) :: id, varid
   type(ncdata_type), intent(in) :: data
   integer  n, ierr
   integer :: start = 1

     select case (Data%datatype)
        case (NCBYTE)
           call ncvptc (id, varid, start, Data%len, Data%bdata,  &
                        Data%len, ierr)
        case (NCCHAR)
           n = Data%len * maxc
           call ncvptc (id, varid, start, Data%len, Data%cdata, n,  &
                        ierr)
        case (NCSHORT)
           call ncvpt  (id, varid, start, Data%len, Data%sdata, ierr)
        case (NCLONG)
           call ncvpt  (id, varid, start, Data%len, Data%ldata, ierr)
        case (NCFLOAT)
           call ncvpt  (id, varid, start, Data%len, Data%fdata, ierr)
        case (NCDOUBLE)
           call ncvpt  (id, varid, start, Data%len, Data%ddata, ierr)
        case default
           print *, 'invalid data type'
           stop 111
     end select

     if (ierr /= 0) then
        call error_handler ('write_ncdata',  &
                         'error returned from ncvpt or ncvptc', FATAL)
     endif

   end subroutine write_ncdata

!-----------------------------------------------------------------------

   subroutine get_data_float (Data, data_out, start, count)

   type(ncdata_type), intent(in)                :: Data
   real(float_type),  intent(out), dimension(:) :: data_out
   integer         ,  intent(in) , optional     :: start, count
   integer   is,ie,knt

     if (Data%datatype /= NCFLOAT) then
         call error_handler ('get_data_float',   &
           'wrong data type in get_data_float', FATAL)
     endif

     is = 1; ie = Data%len
     if (present(start)) is = start
     if (present(count)) ie = is + count - 1
     knt = ie-is+1

     if (knt > size(data_out)) then
         print *, 'ERROR size of data_out too small.'
         stop 111
     endif

     data_out(1:knt) = Data%fdata(is:ie)

   end subroutine get_data_float

!-----------------------------------------------------------------------

   subroutine get_data_double (Data, data_out, start, count)

   type(ncdata_type), intent(in)                :: Data
   real(double_type), intent(out), dimension(:) :: data_out
   integer          , intent(in) , optional     :: start, count
   integer   is,ie,knt

     if (Data%datatype /= NCDOUBLE) then
         call error_handler ('get_data_double',   &
           'wrong data type in get_data_float', FATAL)
     endif

     is = 1; ie = Data%len
     if (present(start)) is = start
     if (present(count)) ie = is + count - 1
     knt = ie-is+1

     if (knt > size(data_out)) then
         print *, 'ERROR size of data_out too small.'
         stop 111
     endif

     data_out(1:knt) = Data%ddata(is:ie)

   end subroutine get_data_double

!-----------------------------------------------------------------------

   function put_data_float (data_in) result (data)

   real(float_type), intent(in), dimension(:) :: data_in
   type(ncdata_type)                          :: data

     call alloc_ncdata ( NCFLOAT, size(data_in), Data )
     Data%fdata(1:Data%len) = data_in(1:Data%len)

   end function put_data_float

!-----------------------------------------------------------------------

   function put_data_double (data_in) result (Data)

   real(double_type), intent(in), dimension(:) :: data_in
   type(ncdata_type)                           :: Data

     call alloc_ncdata ( NCDOUBLE, size(data_in), Data )
     Data%ddata(1:Data%len) = data_in(1:Data%len)

   end function put_data_double

!-----------------------------------------------------------------------

   subroutine copy_ncdata (Data_in, Data_out)
   type(ncdata_type), intent(in)  :: Data_in
   type(ncdata_type), intent(out) :: Data_out

     call alloc_ncdata (Data_in%datatype, Data_in%len, Data_out)

     select case (Data_out%datatype)
        case (NCBYTE)
           Data_out%bdata = Data_in%bdata
        case (NCCHAR)
           Data_out%cdata = Data_in%cdata
        case (NCSHORT)
           Data_out%sdata = Data_in%sdata
        case (NCLONG)
           Data_out%ldata = Data_in%ldata
        case (NCFLOAT)
           Data_out%fdata = Data_in%fdata
        case (NCDOUBLE)
           Data_out%ddata = Data_in%ddata
     end select

   end subroutine copy_ncdata

!-----------------------------------------------------------------------

   subroutine remove_ncdata (Data)
   type(ncdata_type), intent(inout) :: Data

     select case (Data%datatype)
        case (NCBYTE)
           deallocate (Data%bdata)
        case (NCCHAR)
           deallocate (Data%cdata)
        case (NCSHORT)
           deallocate (Data%sdata)
        case (NCLONG)
           deallocate (Data%ldata)
        case (NCFLOAT)
           deallocate (Data%fdata)
        case (NCDOUBLE)
           deallocate (Data%ddata)
        case default
           call error_handler ('remove_ncdata',   &
                               'unknown data type', FATAL)
     end select

     Data%datatype = 0
     Data%len      = 0

   end subroutine remove_ncdata

!-----------------------------------------------------------------------

   subroutine error_handler (routine, message, level)

   character(len=*), intent(in) :: routine, message
   integer         , intent(in) :: level

      print *, 'ERROR in ' // trim(routine)
      print *, trim(message)

      call abort

   end subroutine error_handler

!-----------------------------------------------------------------------

end module ncdata_mod

