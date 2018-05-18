
module ncfile_access_mod

use  ncd_define_mod
use      ncfile_mod, only: ncfile_type
use      ncaxis_mod, only: get_ncaxis_no
use    ncvarmap_mod, only: get_varmap_no

implicit none
private

public  file_name, num_axis, num_var, axis_name, var_name
public  time_axis_index, time_axis_length
public  axis_length, num_var_axis, var_axis, var_axis_len
public  axis_index, var_index

!-----------------------------------------------------------------------

interface axis_length
   module procedure  axis_len_num, axis_len_str
end interface

interface num_var_axis
   module procedure  num_var_axis_num, num_var_axis_str
end interface

interface  var_axis
   module procedure var_axis_num, var_axis_str
end interface

interface  var_axis_len
   module procedure var_axis_len_num, var_axis_len_str
end interface


integer, parameter :: FATAL = 2


contains

!-----------------------------------------------------------------------
!-------------- return the file name -----------------------------------

   function file_name (File) result (fname)

   type(ncfile_type) , intent(in)     :: File
   character(len=len_trim(File%name)) :: fname
   integer    nc

       nc = len_trim(File%name)
       fname(1:nc) = File%name(1:nc)

   end function file_name

!-----------------------------------------------------------------------
!---------- return the number of defined axes (including time) ---------

   function num_axis (File) result (num)

   type(ncfile_type), intent(in) :: File
   integer                          num

       num = File%ndim

   end function num_axis

!-----------------------------------------------------------------------
!----------- return the number of defined variables --------------------
!          ------- (including time avg vars) ---------

   function num_var (File) result (num)

   type(ncfile_type), intent(in) :: File
   integer                          num

       num = File%nvarmap

   end function num_var

!-----------------------------------------------------------------------
!-------------- return the axis name -----------------------------------

   function axis_name (File, axis_index) result (name)

   type(ncfile_type) , intent(in)  :: File
   integer           , intent(in)  :: axis_index
   character(len=len_trim(File%Axis(axis_index)%name)) :: name
   integer    nc

   if (axis_index < 1 .or. axis_index > File%ndim) call error_handler  &
                ('axis_name', 'number of dimensions exceeded.', FATAL)

       nc = len_trim(File%Axis(axis_index)%name)
       name(1:nc) = File%Axis(axis_index)%name(1:nc)

   end function axis_name

!-----------------------------------------------------------------------
!-------------- return the length of the requested axis ----------------

   function axis_len_num (File, axis_index) result (length)

   type(ncfile_type), intent(in)  :: File
   integer          , intent(in)  :: axis_index
   integer                         length

       length = File%Axis(axis_index)%len

   end function axis_len_num

!-----------------------------------------------------------------------

   function axis_len_str (File, axis_name) result (length)

   type(ncfile_type), intent(in)  :: File
   character(len=*) , intent(in)  :: axis_name
   integer                          length, num

       num = get_ncaxis_no (File%Axis(1:File%ndim), axis_name)
       length = axis_len_num(File, num)

   end function axis_len_str

!-----------------------------------------------------------------------

   function time_axis_index (File) result(num)

   type(ncfile_type), intent(in)  :: File
   integer                        :: num

       num = File%recdim

   end function time_axis_index

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

   function time_axis_length (File) result(length)

   type(ncfile_type), intent(in)  :: File
   integer                    length, num

       num = File%recdim
       length = File%Axis(num)%len

   end function time_axis_length

!-----------------------------------------------------------------------
!-------------- return the variable name -------------------------------

   function var_name (File, var_index) result (name)

   type(ncfile_type), intent(in)  :: File
   integer          , intent(in)  :: var_index
   character(len=len_trim(File%Varmap(var_index)%name)) :: name
   integer    nc

   if (var_index < 1 .or. var_index > File%nvarmap) call error_handler &
                ('var_name', 'number of variables exceeded.', FATAL)

       nc = len_trim(File%Varmap(var_index)%name)
       name(1:nc) = File%Varmap(var_index)%name(1:nc)

   end function var_name

!-----------------------------------------------------------------------
!-------------- return the variable dimensions -------------------------

   function num_var_axis_num (File, var_index) result (n)

   type(ncfile_type), intent(in) :: File
   integer          , intent(in) :: var_index
   integer                        n, no

       no = File%Varmap(var_index)%varno
       n = File%Var(no)%ndim

   end function num_var_axis_num

!-----------------------------------------------------------------------

   function num_var_axis_str (File, var_name) result (n)

   type(ncfile_type), intent(in) :: File
   character(len=*) , intent(in) :: var_name
   integer                        n, no

       no = get_varmap_no (File%Varmap(1:File%nvarmap), var_name)
       n = num_var_axis_num (File, no)

   end function num_var_axis_str

!-----------------------------------------------------------------------
!-------------- return the variable dimension sizes --------------------

   subroutine var_axis_len_num (File, var_index, len)

   type(ncfile_type), intent(in)  :: File
   integer          , intent(in)  :: var_index
   integer          , intent(out) :: len(:)

   integer :: i, m, nd, no, dims(4)

      m = size(len)
      no = File%Varmap(var_index)%varno
      nd = File%Var(no)%ndim

      len = 0
      dims(1:nd) = File%Var(no)%vdims(1:nd)
      do i=1,min(m,nd)
         if (dims(i) > 0) len(i) = axis_len_num(File, dims(i))
      enddo

   end subroutine var_axis_len_num

!-----------------------------------------------------------------------

   subroutine var_axis_len_str (File, var_name, len)

   type(ncfile_type), intent(in)  :: File
   character(len=*) , intent(in)  :: var_name
   integer          , intent(out) :: len(:)
   integer     no

      no = get_varmap_no (File%Varmap(1:File%nvarmap), var_name)
      if (no > 0) then
        call var_axis_len_num (File, no, len)
      endif

   end subroutine var_axis_len_str

!-----------------------------------------------------------------------
!-------------- return the variable dimension id numbers ---------------

   subroutine var_axis_num (File, var_index, axes)

   type(ncfile_type), intent(in)  :: File
   integer          , intent(in)  :: var_index
   integer          , intent(out) :: axes(:)

   integer :: i, naxes, ndims, no

      naxes = size(axes)
      no = File%Varmap(var_index)%varno
      ndims = File%Var(no)%ndim

      axes = 0
      do i = 1, min(naxes,ndims)
         axes(i) = File%Var(no)%vdims(i)
      enddo

   end subroutine var_axis_num

!-----------------------------------------------------------------------

   subroutine var_axis_str (File, var_name, axes)

   type(ncfile_type),intent(in)  :: File
   character(len=*), intent(in)  :: var_name
   integer        ,  intent(out) :: axes(:)

   integer :: no

      no = get_varmap_no (File%Varmap(1:File%nvarmap), var_name)
      if (no > 0) then
          call var_axis_num (File, no, axes)
      endif

   end subroutine var_axis_str

!-----------------------------------------------------------------------

   function axis_index (File, axis_name) result (index)

   type(ncfile_type),intent(in)  :: File
   character(len=*), intent(in)  :: axis_name
   integer                       :: index

      index = get_ncaxis_no (File%Axis(1:File%ndim), axis_name)

   end function axis_index

!-----------------------------------------------------------------------

   function var_index (File, var_name) result (index)

   type(ncfile_type),intent(in)  :: File
   character(len=*), intent(in)  :: var_name
   integer                       :: index

      index = get_varmap_no (File%Varmap(1:File%nvarmap), var_name)

   end function var_index

!-----------------------------------------------------------------------

   subroutine error_handler (routine, message, level)
   character(len=*), intent(in) :: routine, message
   integer         , intent(in) :: level

      if ( level == FATAL ) then
          print *, 'ERROR in ',trim(routine),' in ncfile_access_mod'
          print *, trim(message)
          call abort
      else
          print *, 'MESSAGE from ',trim(routine),' in ncfile_access_mod'
          print *, trim(message)
      endif

   end subroutine error_handler
!-----------------------------------------------------------------------
!#######################################################################

end module ncfile_access_mod

