
module ncvar_mod

use ncd_define_mod
use      ncatt_mod, only: ncatt_type, read_ncatt, write_ncatt,  &
                          copy_ncatt

implicit none
private

public :: read_ncvar, write_ncvar, get_ncvar, put_ncvar, get_ncvar_no
public :: copy_ncvar, ncvar_type


   type ncvar_type
       character(len=maxn)     :: name
       integer                 :: id, datatype, natt, ndim,  &
                                  vdims(max_dims)
       type(ncatt_type), pointer :: Att(:)
   end type ncvar_type

CONTAINS

!-----------------------------------------------------------------------

      function read_ncvar (fid, varid) result (Var)

      integer,        intent(in) :: fid, varid
      type(ncvar_type)           :: Var
      integer  i, ierr

        Var%id = varid

        call ncvinq (fid, Var%id, Var%name, Var%datatype, Var%ndim,  &
                     Var%vdims, Var%natt, ierr)

        if (Var%natt > 0) then
            allocate (Var%Att(Var%natt))
            do i=1,Var%natt
               Var%Att(i) = read_ncatt (fid, Var%id, i)
            enddo
        endif

      end function read_ncvar

!-----------------------------------------------------------------------

      subroutine write_ncvar (fid, Var)

      integer,          intent(in)    :: fid
      type(ncvar_type), intent(inout) :: Var
      integer  i, ierr
      integer  ncvid

        call ncvdef (fid, trim(Var%name), Var%datatype, Var%ndim,  &
                     Var%vdims, ierr)

        Var%id = ncvid (fid, trim(Var%name), ierr)

         do i=1,Var%natt
             call write_ncatt (fid, Var%id, Var%Att(i))
         enddo

      end subroutine write_ncvar

!-----------------------------------------------------------------------

      function get_ncvar_no (Var, name) result (no)

      type(ncvar_type), intent(in), dimension(:)  :: Var
      character(len=*), intent(in)                :: name
      integer                                     :: no
      integer  i

        no = 0

        do i = 1, size(Var)
          if ( trim(Var(i)%name) == trim(name) ) then
              no = i
              exit
          endif
        enddo

      end function get_ncvar_no

!-----------------------------------------------------------------------

      function get_var_id (Var, name) result (id)

      type(ncvar_type),  intent(in), dimension(:)  :: Var
      character(len=*),  intent(in)                :: name
      integer                                      :: id
      integer  i

        id = 0

        do i = 1, size(Var)
          if ( trim(Var(i)%name) == trim(name) ) then
               id = Var(i)%id
               exit
          endif
        enddo

      end function get_var_id

!-----------------------------------------------------------------------

      subroutine get_ncvar (Var, nvar, name, datatype, ndim, vdims,  &
                               natt, no)

      type(ncvar_type),  intent(in), dimension(:) :: Var
      integer,           intent(in)               :: nvar
      character(len=*),  intent(in)               :: name
      integer         ,  intent(out)              :: ndim, datatype,  &
                                                     vdims(max_dims), &
                                                     natt, no
      integer  nv

         nv = min(size(Var), nvar)

         no = get_ncvar_no (Var(1:nv), name)

         if (no > 0) then
                 ndim = Var(no)%ndim
                vdims = Var(no)%vdims
                 natt = Var(no)%natt
             datatype = Var(no)%datatype
         endif

      end subroutine get_ncvar

!-----------------------------------------------------------------------

      subroutine put_ncvar (Var, nvar, name, datatype, ndim, vdims, no)

      type(ncvar_type), intent(inout), dimension(:) :: Var
      integer,          intent(inout)               :: nvar
      character(len=*), intent(in)               :: name
      integer         , intent(in)               :: datatype, ndim,  &
                                                    vdims(max_dims)
      integer         , intent(out)              :: no
      integer   nv

         nv = min(size(Var), nvar)

         no = get_ncvar_no (Var(1:nv), name)

!--------- add a new variable ----------

         if (no == 0) then
             if (nvar+1 <= max_nvar) then
                 if (nvar+1 <= size(Var)) then
                     nvar = nvar + 1
                     no = nvar
                     allocate (Var(no)%Att(max_natt))
                 else
                     print *, 'ERROR: var too small.'
                     stop 111
                 endif
             else
                 print *, 'ERROR: max_nvar too small.'
                 stop 111
             endif
         endif

!--------- add/replace variable ----------

             Var(no)%name     = name
             Var(no)%id       = 0
             Var(no)%ndim     = ndim
             Var(no)%vdims    = vdims
             Var(no)%natt     = 0
             Var(no)%datatype = datatype


      end subroutine put_ncvar

!-----------------------------------------------------------------------

   subroutine copy_ncvar (Var_in, Var_out)

   type(ncvar_type), intent(in)  :: Var_in
   type(ncvar_type), intent(out) :: Var_out
   integer :: i

     Var_out%name     = Var_in%name
     Var_out%id       = Var_in%id
     Var_out%datatype = Var_in%datatype
     Var_out%ndim     = Var_in%ndim
     Var_out%vdims    = Var_in%vdims
     Var_out%natt     = Var_in%natt

     allocate ( Var_out%Att (Var_out%natt) )
     do i = 1, Var_out%natt
         call copy_ncatt ( Var_in%Att(i), Var_out%Att(i) )
     enddo

   end subroutine copy_ncvar

!-----------------------------------------------------------------------

end module ncvar_mod

