module ncatt_mod

use ncd_define_mod

implicit none
private

public :: get_ncatt, put_ncatt, read_ncatt, write_ncatt
public :: copy_ncatt, remove_ncatt
public :: ncatt_type

   type ncatt_type
      character(len=maxn) :: name
      integer             :: datatype, len
      character(len=1)   , pointer :: batt(:)
      character(len=maxc), pointer :: catt(:)
      integer(short_type), pointer :: satt(:)
      integer( long_type), pointer :: latt(:)
      real   (float_type), pointer :: fatt(:)
      real  (double_type), pointer :: datt(:)
   end type ncatt_type

   integer, parameter :: FATAL = 2

!-----------------------------------------------------------------------

   interface get_ncatt
      module procedure  get_ncatt_byte
      module procedure  get_ncatt_short_0d, get_ncatt_short_1d
      module procedure  get_ncatt_float_0d, get_ncatt_float_1d
   end interface

   interface put_ncatt
      module procedure  put_ncatt_byte
      module procedure  put_ncatt_short_0d, put_ncatt_short_1d
      module procedure  put_ncatt_float_0d, put_ncatt_float_1d
   end interface

contains

!-----------------------------------------------------------------------

   subroutine alloc_ncatt_data (type, len, Att)

   integer,          intent(in)    :: type, len
   type(ncatt_type), intent(inout) :: Att

!  allocates space in appropriate data type
!  assumes that Att has not been allocated

      select case (type)
         case (NCBYTE:NCCHAR)
            allocate (Att%batt(len))
!!       case (NCCHAR)
!!          allocate (Att%catt(len))
         case (NCSHORT)
            allocate (Att%satt(len))
         case (NCLONG)
            allocate (Att%latt(len))
         case (NCFLOAT)
            allocate (Att%fatt(len))
         case (NCDOUBLE)
            allocate (Att%datt(len))
         case default
            call error_handler ('alloc_ncatt_data',                 &
                                'unknown attribute data type', FATAL)
      end select

      Att%datatype = type
      Att%len      = len

   end subroutine alloc_ncatt_data

!-----------------------------------------------------------------------

   subroutine remove_ncatt (Att)

   type(ncatt_type), intent(inout) :: Att

!  deallocates space in Att
!  assumes that Att has been allocated

      select case (Att%datatype)
         case (NCBYTE:NCCHAR)
            deallocate (Att%batt)
!!       case (NCCHAR)
!!          deallocate (Att%catt)
         case (NCSHORT)
            deallocate (Att%satt)
         case (NCLONG)
            deallocate (Att%latt)
         case (NCFLOAT)
            deallocate (Att%fatt)
         case (NCDOUBLE)
            deallocate (Att%datt)
         case default
            call error_handler ('dealloc_ncatt_data',               &
                                'unknown attribute data type', FATAL)
      end select

      Att%name     = ' '
      Att%datatype = 0
      Att%len      = 0

   end subroutine remove_ncatt

!-----------------------------------------------------------------------

   function read_ncatt (fid, varid, attid) result (Att)

      integer,        intent(in)  :: fid, varid, attid
      type(ncatt_type)            :: Att
      integer :: n, ierr, type, len

        call ncanam (fid, varid, attid, Att%name, ierr)
        call ncainq (fid, varid, Att%name, type, len, ierr)

        call alloc_ncatt_data (type, len, Att)

        select case (Att%datatype)
           case (NCBYTE:NCCHAR)
              call ncagtc (fid, varid, Att%name, Att%batt, Att%len,ierr)
!!         case (NCCHAR)
!!            n = Att%len * maxc
!!            call ncagtc (fid, varid, Att%name, Att%catt, n, ierr)
           case (NCSHORT)
              call ncagt  (fid, varid, Att%name, Att%satt, ierr)
           case (NCLONG)
              call ncagt  (fid, varid, Att%name, Att%latt, ierr)
           case (NCFLOAT)
              call ncagt  (fid, varid, Att%name, Att%fatt, ierr)
           case (NCDOUBLE)
              call ncagt  (fid, varid, Att%name, Att%datt, ierr)
!!         case default
!!            print *, 'invalid attribute data type'
!!            stop 111
        end select

   end function read_ncatt

!-----------------------------------------------------------------------

   subroutine write_ncatt (fid, varid, Att)

      integer,          intent(in) :: fid, varid
      type(ncatt_type), intent(in) :: Att
      integer  i, n, nc, ncx, ierr

        nc = len_trim(Att%name)

        select case (Att%datatype)
           case (NCBYTE:NCCHAR)
              call ncaptc (fid, varid, Att%name(1:nc), Att%datatype,  &
                           Att%len, Att%batt, ierr)
!!         case (NCCHAR)
!!            ncx = 0
!!            do i=1,Att%len
!!               ncx = max(ncx, len_trim(Att%catt(i)))
!!            enddo
!!            n = Att%len * ncx
!!            call ncaptc (fid, varid, Att%name(1:nc), Att%datatype,  &
!!                         n, Att%catt(:)(1:ncx), ierr)
           case (NCSHORT)
              call ncapt  (fid, varid, Att%name(1:nc), Att%datatype,  &
                           Att%len, Att%satt, ierr)
           case (NCLONG)
              call ncapt  (fid, varid, Att%name(1:nc), Att%datatype,  &
                           Att%len, Att%latt, ierr)
           case (NCFLOAT)
              call ncapt  (fid, varid, Att%name(1:nc), Att%datatype,  &
                           Att%len, Att%fatt, ierr)
           case (NCDOUBLE)
              call ncapt  (fid, varid, Att%name(1:nc), Att%datatype,  &
                           Att%len, Att%datt, ierr)
           case default
              print *, 'invalid attribute data type'
              stop 111
        end select

   end subroutine write_ncatt

!-----------------------------------------------------------------------

   function get_ncatt_no (Att, name) result (no)

      type(ncatt_type),  intent(in), dimension(:)  :: Att
      character(len=*),  intent(in)                :: name
      integer                                      :: no
      integer  i, nc

        nc = len_trim(name)
        no = 0

        do i = 1, size(Att)
          if (Att(i)%name(1:len_trim(Att(i)%name)) == name(1:nc)) then
              no = i
              exit
          endif
        enddo

   end function get_ncatt_no

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------

   subroutine get_ncatt_byte (Att, natt, name, string, id)

      type(ncatt_type),  intent(in), dimension(:) :: Att
      integer,           intent(in)               :: natt
      character(len=*),  intent(in)               :: name
      character(len=*),  intent(out)              :: string
      integer,           intent(out)              :: id
      integer  k, n, na

         na = min(size(Att), natt)

         id = get_ncatt_no (Att(1:na), name)

         if (id > 0) then
             if (Att(id)%datatype == NCBYTE .or.  &
                 Att(id)%datatype == NCCHAR) then
                      n = min(Att(id)%len, len(string))
                      do k=1,n
                         string(k:k) = Att(id)%batt(k)
                      enddo
                      do k=n+1,len(string)
                         string(k:k) = ' '
                      enddo
             endif
         endif

   end subroutine get_ncatt_byte

!-----------------------------------------------------------------------

   subroutine get_ncatt_short_1d (Att, natt, name, data, id)

      type(ncatt_type),    intent(in), dimension(:)  :: Att
      integer,             intent(in)                :: natt
      character(len=*),    intent(in)                :: name
      integer(short_type), intent(out), dimension(:) :: data
      integer,             intent(out)               :: id
      integer  k, n, na

         na = min(size(Att), natt)

         id = get_ncatt_no (Att(1:na), name)

         if (id > 0) then
             if (Att(id)%datatype == NCSHORT) then
                     n = min(Att(id)%len, size(data))
                     data(1:n) = Att(id)%satt(1:n)
             endif
         endif

   end subroutine get_ncatt_short_1d

!-----------------------------------------------------------------------

   subroutine get_ncatt_float_1d (Att, natt, name, data, id)

      type(ncatt_type),  intent(in), dimension(:)  :: Att
      integer,           intent(in)                :: natt
      character(len=*),  intent(in)                :: name
      real (float_type), intent(out), dimension(:) :: data
      integer,           intent(out)               :: id
      integer  k, n, na

         na = min(size(Att), natt)

         id = get_ncatt_no (Att(1:na), name)

         if (id > 0) then
             if (Att(id)%datatype == NCFLOAT) then
                     n = min(Att(id)%len, size(data))
                     data(1:n) = Att(id)%fatt(1:n)
             endif
         endif

   end subroutine get_ncatt_float_1d

!-----------------------------------------------------------------------

   subroutine get_ncatt_short_0d (Att, natt, name, data, id)

      type(ncatt_type),    intent(in), dimension(:) :: Att
      integer,             intent(in)               :: natt
      character(len=*),    intent(in)               :: name
      integer(short_type), intent(out)              :: data
      integer,             intent(out)              :: id
      integer(short_type)            , dimension(1) :: data_1d

        call get_ncatt_short_1d (Att, natt, name, data_1d, id)
        if (id > 0) data = data_1d(1)

   end subroutine get_ncatt_short_0d

!-----------------------------------------------------------------------

   subroutine get_ncatt_float_0d (Att, natt, name, data, id)

      type(ncatt_type),  intent(in), dimension(:) :: Att
      integer,           intent(in)               :: natt
      character(len=*),  intent(in)               :: name
      real (float_type), intent(out)              :: data
      integer,           intent(out)              :: id
      real (float_type)            , dimension(1) :: data_1d

        call get_ncatt_float_1d (Att, natt, name, data_1d, id)
        if (id > 0) data = data_1d(1)

   end subroutine get_ncatt_float_0d

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------

   subroutine put_ncatt_byte (Att, natt, name, string, id)

      type(ncatt_type),  intent(inout), dimension(:) :: Att
      integer,           intent(inout)               :: natt
      character(len=*),  intent(in)                  :: name
      character(len=*),  intent(in)                  :: string
      integer,           intent(out)                 :: id
      integer  k, n, na

         na = min(size(Att), natt)

         id = get_ncatt_no (Att(1:na), name)

!        ------ replace string -----

         if (id > 0) then
             if (Att(id)%datatype == NCBYTE .or.  &
                 Att(id)%datatype == NCCHAR) then
                     n = len(string)
                     if (Att(id)%len /= n) then
                        deallocate (Att(id)%batt)
                          allocate (Att(id)%batt(n))
                          Att(id)%len = n
                     endif
                     do k=1,n
                       Att(id)%batt(k) = string(k:k)
                     enddo
                     return
             endif
         endif

!        ------ append string -----

         if (natt+1 > max_natt) return
         natt = natt + 1
         id   = natt
         n = len(string)
         Att(id)%name     = name
         Att(id)%datatype = NCCHAR
         Att(id)%len      = n
         allocate (Att(id)%batt(n))
         do k=1,n
            Att(id)%batt(k) = string(k:k)
         enddo

   end subroutine put_ncatt_byte

!-----------------------------------------------------------------------

   subroutine put_ncatt_short_1d (Att, natt, name, data, id)

      type(ncatt_type),    intent(inout), dimension(:) :: Att
      integer,             intent(inout)               :: natt
      character(len=*),    intent(in)                  :: name
      integer(short_type), intent(in)   , dimension(:) :: data
      integer,             intent(out)                 :: id
      integer  k, n, na

         na = min(size(Att), natt)

         id = get_ncatt_no (Att(1:na), name)

!        ------ replace data -----

         if (id > 0) then
             if (Att(id)%datatype == NCSHORT) then
                     n = size(data)
                     if (Att(id)%len /= n) then
                        deallocate (Att(id)%satt)
                          allocate (Att(id)%satt(n))
                          Att(id)%len = n
                     endif
                     Att(id)%satt(1:n) = data(1:n)
                     return
             endif
         endif

!        ------ append string -----

         if (natt+1 > max_natt) return
         natt = natt + 1
         id   = natt
         n = size(data)
         Att(id)%name     = name
         Att(id)%datatype = NCSHORT
         Att(id)%len      = n
         allocate (Att(id)%satt(n))
         Att(id)%satt(1:n) = data(1:n)

   end subroutine put_ncatt_short_1d

!-----------------------------------------------------------------------

   subroutine put_ncatt_float_1d (Att, natt, name, data, id)

      type(ncatt_type),  intent(inout), dimension(:) :: Att
      integer,           intent(inout)               :: natt
      character(len=*),  intent(in)                  :: name
      real (float_type), intent(in)   , dimension(:) :: data
      integer,           intent(out)                 :: id
      integer  k, n, na

         na = min(size(Att), natt)

         id = get_ncatt_no (Att(1:na), name)

!        ------ replace data -----

         if (id > 0) then
             if (Att(id)%datatype == NCFLOAT) then
                     n = size(data)
                     if (Att(id)%len /= n) then
                        deallocate (Att(id)%fatt)
                          allocate (Att(id)%fatt(n))
                          Att(id)%len = n
                     endif
                     Att(id)%fatt(1:n) = data(1:n)
                     return
             endif
         endif

!        ------ append string -----

         if (natt+1 > max_natt) return
         natt = natt + 1
         id   = natt
         n = size(data)
         Att(id)%name     = name
         Att(id)%datatype = NCFLOAT
         Att(id)%len      = n
         allocate (Att(id)%fatt(n))
         Att(id)%fatt(1:n) = data(1:n)

   end subroutine put_ncatt_float_1d

!-----------------------------------------------------------------------

   subroutine put_ncatt_short_0d (Att, natt, name, data, id)

      type(ncatt_type),    intent(inout), dimension(:) :: Att
      integer,             intent(inout)               :: natt
      character(len=*),    intent(in)                  :: name
      integer(short_type), intent(in)                  :: data
      integer,             intent(out)                 :: id
      integer(short_type),                dimension(1) :: data_1d

         data_1d(1) = data
         call put_ncatt_short_1d (Att, natt, name, data_1d, id)

   end subroutine put_ncatt_short_0d

!-----------------------------------------------------------------------

   subroutine put_ncatt_float_0d (Att, natt, name, data, id)

      type(ncatt_type),  intent(inout), dimension(:) :: Att
      integer,           intent(inout)               :: natt
      character(len=*),  intent(in)                  :: name
      real (float_type), intent(in)                  :: data
      integer,           intent(out)                 :: id
      real (float_type),                dimension(1) :: data_1d

         data_1d(1) = data
         call put_ncatt_float_1d (Att, natt, name, data_1d, id)

   end subroutine put_ncatt_float_0d

!-----------------------------------------------------------------------

   subroutine copy_ncatt (Att_in, Att_out)

      type(ncatt_type), intent(in)  :: Att_in
      type(ncatt_type), intent(out) :: Att_out

      Att_out%name = Att_in%name

      call alloc_ncatt_data (Att_in%datatype, Att_in%len, Att_out)

      select case (Att_out%datatype)
         case (NCBYTE:NCCHAR)
            Att_out%batt = Att_in%batt
!!       case (NCCHAR)
!!          n = Att%len * maxc
!!          Att_out%catt = Att_in%catt
         case (NCSHORT)
            Att_out%satt = Att_in%satt
         case (NCLONG)
            Att_out%latt = Att_in%latt
         case (NCFLOAT)
            Att_out%fatt = Att_in%fatt
         case (NCDOUBLE)
            Att_out%datt = Att_in%datt
      end select

   end subroutine copy_ncatt

!-----------------------------------------------------------------------

   subroutine error_handler (routine, message, level)

   character(len=*), intent(in) :: routine, message
   integer         , intent(in) :: level

      print *, 'ERROR in ' // trim(routine)
      print *, trim(message)

      call abort

   end subroutine error_handler

!-----------------------------------------------------------------------

end module ncatt_mod

