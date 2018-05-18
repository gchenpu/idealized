
module netcdf_file_mod

use  ncd_define_mod
use   ncfile_mod, only: ncfile_type, write_ncfile, read_ncfile, &
                        copy_ncfile
use   ncaxis_mod, only: put_ncaxis, get_ncaxis_no, remove_ncaxis_data
use    ncvar_mod, only: put_ncvar, get_ncvar_no
use    ncatt_mod, only: put_ncatt, get_ncatt
use ncvarmap_mod, only: put_varmap, get_varmap, get_varmap_no

implicit none
private

public  define_file, define_axis, define_time_axis,         &
        define_variable, define_variable_tave, define_tave, &
        read_file_header, write_file_header, copy_file_header
public  get_axis_units, get_axis_long_name, get_axis_cart_name, &
        get_axis_edges_name, get_axis_direction, get_calendar_type,  &
        get_var_units, get_var_long_name,       &
        get_missing_value, get_valid_range


!---- overloaded get routines ----

interface get_axis_units
  module procedure get_axis_units_id, get_axis_units_ch
end interface

interface get_axis_long_name
  module procedure get_axis_long_name_id, get_axis_long_name_ch
end interface

interface get_axis_edges_name
  module procedure get_axis_edges_name_id, get_axis_edges_name_ch
end interface

interface get_axis_cart_name
  module procedure get_axis_cart_name_id, get_axis_cart_name_ch
end interface

interface get_axis_direction
  module procedure get_axis_direction_id, get_axis_direction_ch
end interface

interface get_calendar_type
  module procedure get_calendar_type_id, get_calendar_type_ch,  &
                   get_calendar_type_nn
end interface

interface get_var_units
  module procedure get_var_units_id, get_var_units_ch
end interface

interface get_var_long_name
  module procedure get_var_long_name_id, get_var_long_name_ch
end interface

interface get_missing_value
  module procedure get_missing_value_id, get_missing_value_ch
end interface

interface get_valid_range
  module procedure get_valid_range_id, get_valid_range_ch
end interface


!---- time averaging character strings -----

 character(len=32) :: t1_long_name = 'Start time for average period'
 character(len=32) :: t2_long_name = 'End time for average period'
 character(len=32) :: t3_long_name = 'Length of average period'
!character(len=32) :: t3_long_name = 'Number of items in avg period'


integer, parameter :: FATAL = 2

CONTAINS

!-----------------------------------------------------------------------

     subroutine define_axis (File, axis_name, values, units,  &
                             cart_axis, long_name, direction)

     type(ncfile_type), intent(inout) :: File
     character(len=*) , intent(in)    :: axis_name
     real             , intent(in), dimension(:) :: values
     character(len=*) , intent(in), optional     :: units, cart_axis, &
                                                    long_name, direction

     real(float_type), dimension(size(values)) :: coords
     integer  id, iatt

     coords = values
     call put_ncaxis (File%Axis, File%ndim, axis_name, coords, id)

!    -------- generate axis attributes --------

     if (present(long_name)) then
            call put_ncatt (File%Axis(id)%Att, File%Axis(id)%natt,  &
                           'long_name', long_name, iatt)
     endif

     if (present(cart_axis)) then
            call put_ncatt (File%Axis(id)%Att, File%Axis(id)%natt,  &
                           'cartesian_axis', cart_axis, iatt)
     endif

     if (present(units)) then
            call put_ncatt (File%Axis(id)%Att, File%Axis(id)%natt,  &
                           'units', units, iatt)
     endif

     if (present(direction)) then
            call put_ncatt (File%Axis(id)%Att, File%Axis(id)%natt,  &
                           'positive', direction, iatt)
     endif

     end subroutine define_axis

!-----------------------------------------------------------------------

     subroutine define_time_axis (File, axis_name, units, long_name,  &
                                  calendar, date)

     type(ncfile_type), intent(inout) :: File
     character(len=*) , intent(in)    :: axis_name
     character(len=*) , intent(in), optional :: units, long_name, &
                                                calendar
     integer          , intent(in), optional :: date(:)

     real(double_type), dimension(0) :: values
     character(len=maxc)             :: lunits
     integer  i, nc, id, iatt

     call put_ncaxis (File%Axis, File%ndim, axis_name, values, id)

     File%recdim = id

!    -------- generate axis attributes --------

     if (present(long_name)) then
            call put_ncatt (File%Axis(id)%Att, File%Axis(id)%natt,  &
                           'long_name', long_name, iatt)
     endif

            call put_ncatt (File%Axis(id)%Att, File%Axis(id)%natt,  &
                           'cartesian_axis', 'T', iatt)

     if (present(units)) then
            do i=1,maxc; lunits(i:i)=' '; enddo
            nc = len_trim(units)
            if (present(date)) then
                lunits = units // ' since ' // encode_date (date)
                nc = nc + 28
            else
                lunits = units
            endif
            call put_ncatt (File%Axis(id)%Att, File%Axis(id)%natt,  &
                           'units', lunits(1:nc), iatt)
     endif

     if (present(calendar)) then
            call put_ncatt (File%Axis(id)%Att, File%Axis(id)%natt,  &
                           'calendar_type', calendar, iatt)
     endif

     end subroutine define_time_axis

!-----------------------------------------------------------------------

     function encode_date (date) result (label)

     integer, intent(in) :: date(:)
     character(len=21)   :: label
     integer                        :: i, nd, idate(6)
     character(len=2), dimension(6) :: adate

       nd = size(date)
       if (nd < 3) then
           print *, 'ERROR in encode_date: date too small'
           stop 111
       endif
       nd = min (nd, 6)
       idate = 0; idate(1:nd) = date(1:nd)

       do i=2,6
         if (idate(i) < 10) then
            write (adate(i), 10) idate(i)
         else
            write (adate(i), 20) idate(i)
         endif
       enddo

       write (label, 30) idate(1),(adate(i),i=2,6)

 10    format ('0', i1)
 20    format (i2)
 30    format (i4,'-',a2,'-',a2,'-',a2,':',a2,':',a2,'.0')

     end function encode_date

!-----------------------------------------------------------------------

   subroutine define_variable (File, var_name, axes, units, long_name, &
                               miss_data, valid_range, pack )

     type(ncfile_type),  intent(inout)        :: File
     character(len=*), intent(in)           :: var_name, axes(:)
     character(len=*), intent(in), optional :: units, long_name
     real            , intent(in), optional :: valid_range(2), miss_data
     logical         , intent(in), optional :: pack

     real(float_type)    :: vrange(2), mdata, scalef, aoffset
     real                :: sdat, srng(2), scale, offset
     integer(short_type) :: sdata, srange(2)
     integer, dimension(max_dims) :: vdims
     integer  i, ndim, id, iatt, tid(3)
     logical   do_short

     ndim = size(axes)
     
        vdims = 0
     do i = 1, ndim
        vdims(i) = get_ncaxis_no (File%Axis(1:File%ndim), axes(i))
     enddo

!   ------ check for packed 16-bit integers -------
     do_short = .false.
     if (present(pack)) then
         if (present(valid_range)) then
             scale = abs(valid_range(2)-valid_range(1))/65500.
             if (scale <= 1.0) then
                 do_short = pack
                 offset   = 0.5*(valid_range(1)+valid_range(2))
             endif
         else
             scale = 1.0; offset = 0.0
             do_short = pack
         endif
     endif

!   ------- add new variable --------
    if (do_short) then
      call put_ncvar (File%Var, File%nvar, var_name, NCSHORT,  &
                      ndim, vdims, id)
    else
      call put_ncvar (File%Var, File%nvar, var_name, NCFLOAT,  &
                      ndim, vdims, id)
    endif

!    ---------- generate variable attributes ---------

     if (present(long_name)) then
            call put_ncatt (File%Var(id)%Att, File%Var(id)%natt,  &
                           'long_name', long_name, iatt)
     endif

     if (present(units)) then
            call put_ncatt (File%Var(id)%Att, File%Var(id)%natt,  &
                           'units', units, iatt)
     endif

     if (present(valid_range)) then
            if (do_short) then
               if (scale == 1.0 .and. offset == 0.0) then
                   srange = int(valid_range+sign(0.5,valid_range))
               else
                   srng   = (valid_range-offset)/scale
                   srange = int(srng+sign(0.5,srng))
               endif
               call put_ncatt (File%Var(id)%Att, File%Var(id)%natt,  &
                              'valid_range', srange, iatt)
            else
               vrange = valid_range
               call put_ncatt (File%Var(id)%Att, File%Var(id)%natt,  &
                              'valid_range', vrange, iatt)
            endif
     endif

     if (present(miss_data)) then
            if (do_short) then
               if (scale == 1.0 .and. offset == 0.0) then
                   sdata = int(miss_data+sign(0.5,miss_data))
               else
                   sdat  = (miss_data-offset)/scale
                   sdata = int(sdat+sign(0.5,sdat))
               endif
               call put_ncatt (File%Var(id)%Att, File%Var(id)%natt,  &
                              'missing_value', sdata, iatt)
            else
               mdata = miss_data
               call put_ncatt (File%Var(id)%Att, File%Var(id)%natt,  &
                              'missing_value', mdata, iatt)
            endif
     endif

     if (do_short .and. present(valid_range)) then
            scalef = scale; aoffset = offset
            call put_ncatt (File%Var(id)%Att, File%Var(id)%natt,  &
                           'scale_factor', scalef, iatt)
            call put_ncatt (File%Var(id)%Att, File%Var(id)%natt,  &
                           'add_offset',  aoffset, iatt)
     endif


!    ----------- save new variable mapping information ------

     tid = 0
     if (do_short) then
         call put_varmap (File%Varmap, File%nvarmap,   &
                          var_name, id, tid, do_short, scale, offset)
     else
         call put_varmap (File%Varmap, File%nvarmap, var_name, id, tid)
     endif


     end subroutine define_variable

!-----------------------------------------------------------------------

     subroutine define_variable_tave (File, var_name, avg_name )

     type(ncfile_type), intent(inout)        :: File
     character(len=*) , intent(in)           :: var_name
     character(len=*) , intent(in), optional :: avg_name
     integer, dimension(max_dims) :: vdims
     integer  i, n, nc, idim, ndim, id, iatt, iatt_units, tid(3)
     character(len=maxn)      :: aname, avg_name_local
     character(len=maxc)      :: units
     character(len=3*maxn+15) :: time_string
     real     scale, offset
     logical  pack
     integer, dimension(max_dims) :: axes


!    ---------- get variable number ---------

      call get_varmap (File%Varmap, File%nvarmap, var_name,  &
                       id, pack, scale, offset)

      if (id == 0) then
          print *, ' ERROR in define_variable_tave: &
                   & requested variable does not exist.'
          stop 111
      endif

!    -------- make sure this variable has a record dimension ------
!    --------    and then get the record dimension units     ------

      if (File%recdim == 0) then
          print *, ' ERROR in define_variable_tave: &
                   & no record dimension for this file.'
          stop 111
      endif

      idim = 0; iatt_units = 0
      do i = 1, File%Var(id)%ndim
         if (File%Var(id)%vdims(i) == File%recdim) then
            idim = File%recdim
            do n=1,maxc; units(n:n)=' '; enddo
            call get_ncatt (File%Axis(idim)%Att, File%Axis(idim)%natt, &
                            'units', units, iatt_units)
            exit
         endif
      enddo

      if (idim == 0) then
          print *, ' ERROR in define_variable_tave: &
                & requested variable does not have a record dimension.'
          stop 111
      endif

      avg_name_local = var_name
      if (present(avg_name)) avg_name_local = avg_name

!    ----------- time averaging information ----------

      tid = 0
      n = len_trim(avg_name_local)
      nc = 3*n + 11
!del  nc = 3*n + 15
      time_string(1:nc) = avg_name_local(1:n) // '_T1,' //  &
                          avg_name_local(1:n) // '_T2,' //  &
                          avg_name_local(1:n) // '_DT'
!del                      avg_name_local(1:n) // '_NITEMS'
      nc = min(nc, maxc)

!     ----- add/replace time average attribute ----

      call put_ncatt (File%Var(id)%Att, File%Var(id)%natt,  &
                     'time_avg_info', time_string(1:nc), iatt)


      call define_tave (File, avg_name_local, tid)


!    ---------- add variable mapping information ----------

      call put_varmap (File%Varmap, File%nvarmap,   &
                       var_name, id, tid, pack, scale, offset)


     end subroutine define_variable_tave

!-----------------------------------------------------------------------

  subroutine define_tave (File, avg_name, tid)

  type(ncfile_type), intent(inout) :: File
  character(len=*) , intent(in)    :: avg_name
  integer, optional, intent(out)   :: tid(3)

  character(len=maxn)      :: aname
  character(len=maxc)      :: units
  integer :: n, id(3), idim, iatt, iatt_units, vdims(max_dims)

     n = len_trim(avg_name)
     if (n == 0) return

!---- make sure there is a record dimension and then get the units -----

   if (File%recdim == 0) then
       print *, ' ERROR in define_tave: &
                & no record dimension for this file.'
       stop 111
   endif

      idim = File%recdim
      call get_ncatt (File%Axis(idim)%Att, File%Axis(idim)%natt, &
                         'units', units, iatt_units)

     vdims = 0; vdims(1) = File%recdim

     aname(1:n+3) = avg_name(1:n) // '_T1'
     id(1) = get_ncvar_no (File%Var(1:File%nvar), aname(1:n+3))
     if (id(1) == 0) then
            call put_ncvar (File%Var, File%nvar, aname(1:n+3),  &
                            NCDOUBLE, 1, vdims, id(1))
            if (iatt_units > 0) then
               call put_ncatt (File%Var(id(1))%Att,   &
                               File%Var(id(1))%natt,  &
                               'units', trim(units), iatt)
            endif
            call put_ncatt (File%Var(id(1))%Att, File%Var(id(1))%natt, &
                            'long_name', trim(t1_long_name), iatt)
     endif

     aname(1:n+3) = avg_name(1:n) // '_T2'
     id(2) = get_ncvar_no (File%Var(1:File%nvar), aname(1:n+3))
     if (id(2) == 0) then
            call put_ncvar (File%Var, File%nvar, aname(1:n+3),  &
                            NCDOUBLE, 1, vdims, id(2))
            if (iatt_units > 0) then
               call put_ncatt (File%Var(id(2))%Att,   &
                               File%Var(id(2))%natt,  &
                               'units', trim(units), iatt)
            endif
            call put_ncatt (File%Var(id(2))%Att, File%Var(id(2))%natt, &
                            'long_name', trim(t2_long_name), iatt)
     endif

     aname(1:n+3) = avg_name(1:n) // '_DT'
!del aname(1:n+7) = avg_name(1:n) // '_NITEMS'
     id(3) = get_ncvar_no (File%Var(1:File%nvar), aname(1:n+3))
     if (id(3) == 0) then
            call put_ncvar (File%Var, File%nvar, aname(1:n+3),  &
                            NCDOUBLE, 1, vdims, id(3))
            if (iatt_units > 0) then
               call put_ncatt (File%Var(id(3))%Att,   &
                               File%Var(id(3))%natt,  &
                               'units', trim(units), iatt)
!del                           'units', 'none', iatt)
            endif
            call put_ncatt (File%Var(id(3))%Att, File%Var(id(3))%natt, &
                            'long_name', trim(t3_long_name), iatt)
     endif

     if (present(tid)) tid = id

  end subroutine define_tave

!-----------------------------------------------------------------------

     subroutine define_file (File, file_name, title, history, coards)

     type(ncfile_type), intent(inout)        :: File
     character(len=*),  intent(in)           :: file_name
     character(len=*),  intent(in), optional :: title, history
     logical,           intent(in), optional :: coards
     integer   iatt

        allocate (File%Axis  (max_ndim))
        allocate (File%Var   (max_nvar))
        allocate (File%Varmap(max_nvar))
        allocate (File%Att   (max_natt))

        File%id      = 0
        File%name    = file_name
        File%ndim    = 0
        File%nvar    = 0
        File%nvarmap = 0
        File%natt    = 0
        File%recdim  = 0

        File%define_mode = .false.

        if (present(title)) then
                 call put_ncatt (File%Att, File%natt,  &
                                'title', trim(title), iatt)
        endif

        if (present(history)) then
                 call put_ncatt (File%Att, File%natt,  &
                                'history', trim(history), iatt)
        endif

        if (present(coards)) then
            if (coards)  call put_ncatt (File%Att, File%natt,  &
                                        'Conventions', 'COARDS', iatt)
        endif

     end subroutine define_file

!-----------------------------------------------------------------------

 function read_file_header (file_name) result (File)

 character(len=*), intent(in)  :: file_name
 type(ncfile_type)             :: File

 integer  i, no, nc, n1, tno, sno, ano, tavg_no(3)
 real(float_type) :: scalef, aoffset
 real             :: scale ,  offset
 character(len=maxc) :: string, tstring
 logical             :: pack

!---- read header information ----

      File = read_ncfile(file_name)


!----------------------------------------
!---- setup varmap_type information -----
!----------------------------------------
     allocate (File%Varmap(max_nvar))
     File%nvarmap = 0

     do i = 1, File%nvar
       nc = len_trim(File%Var(i)%name)
!      --- skip time averaging fields ---
       if (nc >= 3) then
         n1 = nc-2
         if ( File%Var(i)%name(nc-2:nc) == '_T1' ) cycle
         if ( File%Var(i)%name(nc-2:nc) == '_T2' ) cycle
         if ( File%Var(i)%name(nc-2:nc) == '_DT' ) cycle
       endif
       if (nc >= 7) then
         n1 = nc-6
         if ( File%Var(i)%name(nc-2:nc) == '_NITEMS' ) cycle
       endif

!      call get_ncatt (File%Var(i)%Att, File%Var(i)%natt, 'long_name', &
!                     string, no)
!      if (no > 0) then

            call get_ncatt (File%Var(i)%Att, File%Var(i)%natt,&
                            'time_avg_info', tstring, tno)
            tavg_no = 0
            if (tno > 0) call set_tavg_no (File, tstring, tavg_no)

            call get_ncatt (File%Var(i)%Att, File%Var(i)%natt, &
                            'scale_factor', scalef, sno)

            call get_ncatt (File%Var(i)%Att, File%Var(i)%natt, &
                            'add_offset', aoffset, ano)

            pack = .false.
            if (File%Var(i)%datatype == NCSHORT) pack = .true.

            if (sno > 0 .and. ano > 0) then
               scale = scalef; offset = aoffset; pack = .true.
               call put_varmap (File%Varmap, File%nvarmap,  &
                                File%Var(i)%name, i,        &
                                tavg_no, pack, scale, offset)
            else
               call put_varmap (File%Varmap, File%nvarmap,  &
                                File%Var(i)%name, i,  &
                                tavg_no, pack)
            endif
!      endif
                                       
     enddo

 end function read_file_header

!-----------------------------------------------------------------------

 subroutine write_file_header (File)
 type(ncfile_type), intent(inout) :: File

       call write_ncfile (File)

 end subroutine write_file_header

!-----------------------------------------------------------------------

 function copy_file_header (File_in, file_name_out, title, coards)  &
                    result (File_out)
 type(ncfile_type), intent(in)           :: File_in
character(len=*), intent(in), optional :: file_name_out, title
logical,          intent(in), optional :: coards
 type(ncfile_type)                       :: File_out
  integer  iatt, rdim

!!!!!   File_out = File_in
        call copy_ncfile (File_in, File_out)

!     ------ remove record dimension data ------
        call remove_ncaxis_data (File_out%Axis(File_out%recdim))

        if (present(file_name_out)) File_out%name = file_name_out

        if (present(title)) then
                 call put_ncatt (File_out%Att, File_out%natt,  &
                                'title', trim(title), iatt)
        endif

        if (present(coards)) then
            if (coards)  call put_ncatt (File_out%Att, File_out%natt,  &
                                        'Conventions', 'COARDS', iatt)
        endif


 end function copy_file_header

!-----------------------------------------------------------------------

 subroutine set_tavg_no (File, string, tavg_no)

 type(ncfile_type) , intent(in)  :: File
 character(len=*), intent(in)  :: string
 integer         , intent(out) :: tavg_no(3)

 integer   is(3), ie(3), nc, i, k

      k = 1; is(k) = 1
      do i = 1, len_trim(string)
          if (string(i:i) /= ',') cycle
          ie(k) = i-1
          k = k+1
          if (k == 4) exit
          is(k) = i+1
      enddo
          ie(3) = len_trim(string)

      do k = 1, 3
          tavg_no(k) =  &
              get_ncvar_no (File%Var(1:File%nvar), string(is(k):ie(k)))
      enddo

 end subroutine set_tavg_no

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!---------- LOGICAL FUNCTIONS TO RETRIEVE OPTIONAL ARGUMENTS -----------
!-----------------------------------------------------------------------
!---------------------- axis units -------------------------------------

 function get_axis_units_id (File, axis_num, units) result (flag)

 type(ncfile_type) , intent(in)  :: File
 integer         , intent(in)  :: axis_num
 character(len=*), intent(out) :: units
 logical                       :: flag
 integer  i, iatt, nc, ncs
 character(len=maxc) :: string

   if (axis_num < 1 .or. axis_num > File%ndim) then
        call error_handler ('get_axis_units_id',    &
                            'invalid axis index.', FATAL)
   endif

   call get_ncatt (File%Axis(axis_num)%Att, File%Axis(axis_num)%natt, &
                  'units', string, iatt)

   nc = len(units)
   do i=1,nc; units(i:i)=' '; enddo

   if (iatt > 0) then
       flag = .true.
       ncs = min(nc, len_trim(string))
       units(1:ncs) = string(1:ncs)
   else
       flag = .false.
   endif

 end function get_axis_units_id

!---------------------- axis long_name ---------------------------------

 function get_axis_long_name_id (File, axis_num, long_name)  &
                         result (flag)

 type(ncfile_type), intent(in)  :: File
 integer          , intent(in)  :: axis_num
 character(len=*) , intent(out) :: long_name
 logical                        :: flag
 integer  i, iatt, nc, ncs
 character(len=maxc) :: string

   if (axis_num < 1 .or. axis_num > File%ndim) then
        call error_handler ('get_axis_long_name_id',    &
                            'invalid axis index.', FATAL)
   endif

   call get_ncatt (File%Axis(axis_num)%Att, File%Axis(axis_num)%natt, &
                  'long_name', string, iatt)

   nc = len(long_name)
   do i=1,nc; long_name(i:i)=' '; enddo

   if (iatt > 0) then
       flag = .true.
       ncs = min(nc, len_trim(string))
       long_name(1:ncs) = string(1:ncs)
   else
       flag = .false.
   endif

 end function get_axis_long_name_id

!---------------------- axis long_name ---------------------------------

 function get_axis_edges_name_id (File, axis_num, edges_name)  &
                          result (flag)

 type(ncfile_type), intent(in)  :: File
 integer          , intent(in)  :: axis_num
 character(len=*) , intent(out) :: edges_name
 logical                        :: flag
 integer  i, iatt, nc, ncs
 character(len=maxc) :: string

   if (axis_num < 1 .or. axis_num > File%ndim) then
        call error_handler ('get_axis_edges_name_id',    &
                            'invalid axis index.', FATAL)
   endif

   call get_ncatt (File%Axis(axis_num)%Att, File%Axis(axis_num)%natt, &
                  'edges', string, iatt)

   nc = len(edges_name)
   do i=1,nc; edges_name(i:i)=' '; enddo

   if (iatt > 0) then
       flag = .true.
       ncs = min(nc, len_trim(string))
       edges_name(1:ncs) = string(1:ncs)
   else
       flag = .false.
   endif

 end function get_axis_edges_name_id

!---------------------- axis cart_name ---------------------------------

 function get_axis_cart_name_id (File, axis_num, cart_name)  &
                         result (flag)

 type(ncfile_type), intent(in)  :: File
 integer          , intent(in)  :: axis_num
 character(len=*) , intent(out) :: cart_name
 logical                        :: flag
 integer  i, iatt, nc, ncs
 character(len=maxc) :: string

   if (axis_num < 1 .or. axis_num > File%ndim) then
        call error_handler ('get_axis_cart_name_id',    &
                            'invalid axis index.', FATAL)
   endif

   call get_ncatt (File%Axis(axis_num)%Att, File%Axis(axis_num)%natt, &
                   'cartesian_axis', string, iatt)

   nc = len(cart_name)
   do i=1,nc; cart_name(i:i)=' '; enddo

   if (iatt > 0) then
       flag = .true.
       ncs = min(nc, len_trim(string))
       cart_name(1:ncs) = string(1:ncs)
   else
       flag = .false.
   endif

 end function get_axis_cart_name_id

!---------------------- axis (positive) direction ----------------------

 function get_axis_direction_id (File, axis_num, direction)  &
                         result (flag)

 type(ncfile_type), intent(in)  :: File
 integer          , intent(in)  :: axis_num
 character(len=*) , intent(out) :: direction
 logical                        :: flag
 integer  i, iatt, nc, ncs
 character(len=maxc) :: string

   if (axis_num < 1 .or. axis_num > File%ndim) then
        call error_handler ('get_axis_direction_id',    &
                            'invalid axis index.', FATAL)
   endif

   call get_ncatt (File%Axis(axis_num)%Att, File%Axis(axis_num)%natt, &
                  'positive', string, iatt)

   nc = len(direction)
   do i=1,nc; direction(i:i)=' '; enddo

   if (iatt > 0) then
       flag = .true.
       ncs = min(nc, len_trim(string))
       direction(1:ncs) = string(1:ncs)
   else
       flag = .false.
   endif

 end function get_axis_direction_id

!---------------------- time axis calendar type ------------------------

 function get_calendar_type_id (File, axis_num, calendar)  &
                        result (flag)

 type(ncfile_type), intent(in)  :: File
 integer          , intent(in)  :: axis_num
 character(len=*) , intent(out) :: calendar
 logical                        :: flag
 integer  i, iatt, nc, ncs
 character(len=maxc) :: string

   if (axis_num < 1 .or. axis_num > File%ndim) then
        call error_handler ('get_calendar_type_id',     &
                            'invalid axis index.', FATAL)
   endif

   call get_ncatt (File%Axis(axis_num)%Att, File%Axis(axis_num)%natt, &
                  'calendar_type', string, iatt)

   nc = len(calendar)
   do i=1,nc; calendar(i:i)=' '; enddo

   if (iatt > 0) then
       flag = .true.
       ncs = min(nc, len_trim(string))
       calendar(1:ncs) = string(1:ncs)
   else
       flag = .false.
   endif

 end function get_calendar_type_id

!-----------------------------------------------------------------------
!---------------------- variable units ---------------------------------

 function get_var_units_id (File, var_num, units) result (flag)

 type(ncfile_type) , intent(in)  :: File
 integer         , intent(in)  :: var_num
 character(len=*), intent(out) :: units
 logical                       :: flag
 integer  i, iatt, no, nc, ncs
 character(len=maxc) :: string

   if (var_num < 1 .or. var_num > File%nvarmap) then
        call error_handler ('get_var_units_id',     &
                            'invalid variable index.', FATAL)
   endif

   no = File%Varmap(var_num)%varno
   call get_ncatt (File%Var(no)%Att, File%Var(no)%natt, &
                  'units', string, iatt)

   nc = len(units)
   do i=1,nc; units(i:i)=' '; enddo

   if (iatt > 0) then
       flag = .true.
       ncs = min(nc, len_trim(string))
       units(1:ncs) = string(1:ncs)
   else
       flag = .false.
   endif

 end function get_var_units_id

!---------------------- variable long_name -----------------------------

 function get_var_long_name_id (File, var_num, long_name)  &
                        result (flag)

 type(ncfile_type), intent(in)  :: File
 integer          , intent(in)  :: var_num
 character(len=*) , intent(out) :: long_name
 logical                        :: flag
 integer  i, iatt, no, nc, ncs
 character(len=maxc) :: string

   if (var_num < 1 .or. var_num > File%nvarmap) then
        call error_handler ('get_var_long_name_id',  &
                            'invalid variable index.', FATAL)
   endif

   no = File%Varmap(var_num)%varno
   call get_ncatt (File%Var(no)%Att, File%Var(no)%natt, &
                  'long_name', string, iatt)

   nc = len(long_name)
   do i=1,nc; long_name(i:i)=' '; enddo

   if (iatt > 0) then
       flag = .true.
       ncs = min(nc, len_trim(string))
       long_name(1:ncs) = string(1:ncs)
   else
       flag = .false.
   endif

 end function get_var_long_name_id

!---------------- special value feature --------------------------------

 function get_missing_value_id (File, var_num, value) result (flag)

 type(ncfile_type), intent(in)  :: File
 integer          , intent(in)  :: var_num
 real             , intent(out) :: value
 logical                        :: flag
 integer   varno, iatt
 logical   pack
 real      scale, offset
 real(float_type)    :: fval
 integer(short_type) :: sval

    if (var_num < 1 .or. var_num > File%nvarmap) then
        call error_handler ('get_missing_value_id',  &
                            'invalid variable index.', FATAL)
    endif

    call get_varmap (File%Varmap, File%nvarmap,  &
                     File%Varmap(var_num)%name,  &
                     varno, pack, scale, offset)

    if (pack) then
       call get_ncatt (File%Var(varno)%Att, File%Var(varno)%natt,  &
                       'missing_value', sval, iatt)
       if (iatt > 0) then
          if (scale == 1.0 .and. offset == 0.0) then
              value = real(sval)
          else
              value = real(sval)*scale+offset
          endif
          flag = .true.
       else
          flag = .false.
          value = 0.0
       endif
    else
       call get_ncatt (File%Var(varno)%Att, File%Var(varno)%natt,  &
                       'missing_value', fval, iatt)
       if (iatt > 0) then
          value = fval
          flag = .true.
       else
          flag = .false.
          value = 0.0
       endif
    endif

 end function get_missing_value_id

!---------------- valid range ------------------------------------------

 function get_valid_range_id (File, var_num, vrange) result (flag)

 type(ncfile_type), intent(in)  :: File
 integer          , intent(in)  :: var_num
 real             , intent(out) :: vrange(2)
 logical                        :: flag
 integer   varno, iatt
 logical   pack
 real      scale, offset
 real(float_type)    :: fval(2)
 integer(short_type) :: sval(2)

   if (var_num < 1 .or. var_num > File%nvarmap) then
       call error_handler ('get_valid_range_id',  &
                           'invalid variable index', FATAL)
   endif

!     ---- get scale facor and add offset ----

    call get_varmap (File%Varmap, File%nvarmap,  &
                     File%Varmap(var_num)%name,  &
                     varno, pack, scale, offset)

!     ---- get valid range attribute ----
!     pack=.false.    ! for old version

    if (pack) then
       call get_ncatt (File%Var(varno)%Att, File%Var(varno)%natt,  &
                       'valid_range', sval, iatt)
       if (iatt > 0) then
          if (scale == 1.0 .and. offset == 0.0) then
              vrange = real(sval)
          else
              vrange = real(sval)*scale+offset
          endif
          flag = .true.
       else
          flag = .false.
          vrange = 0.0
       endif
    else
       call get_ncatt (File%Var(varno)%Att, File%Var(varno)%natt,  &
                       'valid_range', fval, iatt)
       if (iatt > 0) then
          vrange = fval
          flag = .true.
       else
          flag = .false.
          vrange = 0.0
       endif
    endif

 end function get_valid_range_id

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!----- OVERLOADED (using char string name) LOGICAL FUNCTIONS -----
!-----------------------------------------------------------------------
!---------------------- axis units -------------------------------------

 function get_axis_units_ch (File, axis_name, units) result (flag)

 type(ncfile_type), intent(in)  :: File
 character(len=*),  intent(in)  :: axis_name
 character(len=*),  intent(out) :: units
 logical :: flag
 integer :: index

   index = get_ncaxis_no (File%Axis(1:File%ndim), axis_name)

   if (index == 0) call error_handler ('get_axis_units_ch',       &
                                       'invalid axis name.', FATAL)

   flag = get_axis_units_id (File, index, units)

 end function get_axis_units_ch

!---------------------- axis long_name ---------------------------------

 function get_axis_long_name_ch (File, axis_name, long_name)  &
                         result (flag)

 type(ncfile_type), intent(in)  :: File
 character(len=*) , intent(in)  :: axis_name
 character(len=*) , intent(out) :: long_name
 logical :: flag
 integer :: index

   index = get_ncaxis_no (File%Axis(1:File%ndim), axis_name)

   if (index == 0) call error_handler ('get_axis_long_name_ch',   &
                                       'invalid axis name.', FATAL)

   flag = get_axis_long_name_id (File, index, long_name)

 end function get_axis_long_name_ch

!---------------------- axis long_name ---------------------------------

 function get_axis_edges_name_ch (File, axis_name, edges_name)  &
                          result (flag)

 type(ncfile_type), intent(in)  :: File
 character(len=*) , intent(in)  :: axis_name
 character(len=*) , intent(out) :: edges_name
 logical :: flag
 integer :: index

   index = get_ncaxis_no (File%Axis(1:File%ndim), axis_name)

   if (index == 0) call error_handler ('get_axis_edges_name_ch',   &
                                       'invalid axis name.', FATAL)

   flag = get_axis_edges_name_id (File, index, edges_name)

 end function get_axis_edges_name_ch

!---------------------- axis cart_name ---------------------------------

 function get_axis_cart_name_ch (File, axis_name, cart_name)  &
                         result (flag)

 type(ncfile_type), intent(in)  :: File
 character(len=*) , intent(in)  :: axis_name
 character(len=*) , intent(out) :: cart_name
 logical :: flag
 integer :: index

   index = get_ncaxis_no (File%Axis(1:File%ndim), axis_name)

   if (index == 0) call error_handler ('get_axis_cart_name_ch',   &
                                       'invalid axis name.', FATAL)

   flag = get_axis_cart_name_id (File, index, cart_name)

 end function get_axis_cart_name_ch

!---------------------- axis (positive) direction ----------------------

 function get_axis_direction_ch (File, axis_name, direction)  &
                         result (flag)

 type(ncfile_type), intent(in)  :: File
 character(len=*) , intent(in)  :: axis_name
 character(len=*) , intent(out) :: direction
 logical :: flag
 integer :: index

   index = get_ncaxis_no (File%Axis(1:File%ndim), axis_name)

   if (index == 0) call error_handler ('get_axis_direction_ch',   &
                                       'invalid axis name.', FATAL)

   flag = get_axis_direction_id (File, index, direction)

 end function get_axis_direction_ch

!---------------------- time axis calendar type ------------------------

 function get_calendar_type_ch (File, axis_name, calendar)  &
                        result (flag)

 type(ncfile_type), intent(in)  :: File
 character(len=*) , intent(in)  :: axis_name
 character(len=*) , intent(out) :: calendar
 logical :: flag
 integer :: index

   index = get_ncaxis_no (File%Axis(1:File%ndim), axis_name)

   if (index == 0) call error_handler ('get_calendar_type_ch',    &
                                       'invalid axis name.', FATAL)

   flag = get_calendar_type_id (File, index, calendar)

 end function get_calendar_type_ch

!---------------------- time axis calendar type ------------------------

 function get_calendar_type_nn (File, calendar) result (flag)

 type(ncfile_type), intent(in)  :: File
 character(len=*) , intent(out) :: calendar
 logical :: flag
 integer :: index

   index = File%recdim

   if (index > 0) then
       flag = get_calendar_type_id (File, index, calendar)
   else
       flag = .false.
   endif

 end function get_calendar_type_nn

!-----------------------------------------------------------------------
!---------------------- variable units ---------------------------------

 function get_var_units_ch (File, var_name, units) result (flag)

 type(ncfile_type), intent(in)  :: File
 character(len=*) , intent(in)  :: var_name
 character(len=*) , intent(out) :: units
 logical :: flag
 integer :: index

   index = get_varmap_no (File%Varmap(1:File%nvarmap), var_name)

   if (index == 0) call error_handler ('get_var_units_ch',      &
                                 'invalid variable name.', FATAL)

   flag = get_var_units_id (File, index, units)

 end function get_var_units_ch

!---------------------- variable long_name -----------------------------

 function get_var_long_name_ch (File, var_name, long_name)  &
                        result (flag)

 type(ncfile_type), intent(in)  :: File
 character(len=*) , intent(in)  :: var_name
 character(len=*) , intent(out) :: long_name
 logical :: flag
 integer :: index

   index = get_varmap_no (File%Varmap(1:File%nvarmap), var_name)

   if (index == 0) call error_handler ('get_var_long_name_ch',  &
                                 'invalid variable name.', FATAL)

   flag = get_var_long_name_id (File, index, long_name)

 end function get_var_long_name_ch

!---------------- special value feature --------------------------------

 function get_missing_value_ch (File, var_name, value) result (flag)

 type(ncfile_type), intent(in)  :: File
 character(len=*) , intent(in)  :: var_name
 real             , intent(out) :: value
 logical :: flag
 integer :: index

   index = get_varmap_no (File%Varmap(1:File%nvarmap), var_name)

   if (index == 0) call error_handler ('get_missing_value_ch',  &
                                 'invalid variable name.', FATAL)

   flag = get_missing_value_id (File, index, value)

 end function get_missing_value_ch

!---------------- valid range ------------------------------------------

 function get_valid_range_ch (File, var_name, vrange) result (flag)

 type(ncfile_type), intent(in)  :: File
 character(len=*) , intent(in)  :: var_name
 real             , intent(out) :: vrange(2)
 logical :: flag
 integer :: index

   index = get_varmap_no (File%Varmap(1:File%nvarmap), var_name)

   if (index == 0) call error_handler ('get_valid_range_ch',    &
                                 'invalid variable name.', FATAL)

   flag = get_valid_range_id (File, index, vrange)

 end function get_valid_range_ch

!-----------------------------------------------------------------------
!#######################################################################

   subroutine error_handler (routine, message, level)
   character(len=*), intent(in) :: routine, message
   integer         , intent(in) :: level

      if ( level == FATAL ) then
          print *, 'ERROR in ',trim(routine),' in netcdf_file_mod'
          print *, trim(message)
          call abort
      else
          print *, 'MESSAGE from ',trim(routine),' in netcdf_file_mod'
          print *, trim(message)
      endif

   end subroutine error_handler

!#######################################################################

end module netcdf_file_mod

