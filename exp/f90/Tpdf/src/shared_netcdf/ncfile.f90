
module ncfile_mod

use ncd_define_mod
use     ncaxis_mod, only: ncaxis_type, read_ncaxis, write_ncaxis, &
                                      close_ncaxis, copy_ncaxis
use      ncvar_mod, only:  ncvar_type, read_ncvar , write_ncvar, &
                           copy_ncvar
use      ncatt_mod, only:  ncatt_type, read_ncatt , write_ncatt, &
                            copy_ncatt
use   ncvarmap_mod, only: ncvarmap_type, copy_varmap

implicit none

!public :: read_ncfile, write_ncfile, close_ncfile
 public :: copy_ncfile

      type ncfile_type
           integer             :: id, ndim, nvar, natt, recdim, nvarmap
           character(len=maxn) :: name
           logical             :: define_mode
           type(  ncaxis_type), pointer ::   Axis(:)
           type(   ncvar_type), pointer ::    Var(:)
           type(   ncatt_type), pointer ::    Att(:)
           type(ncvarmap_type), pointer :: Varmap(:)
      end type ncfile_type

CONTAINS

!-----------------------------------------------------------------------

      function read_ncfile (name) result (File)

      character(len=*), intent(in) :: name
      type(ncfile_type)            :: File
      integer   ::  i, ierr, nv, nvar_total
      integer  ncopn

         File%name = name
         File%id = ncopn (name, NCNOWRIT, ierr)
         call ncinq (File%id, File%ndim, nvar_total, File%natt,  &
                     File%recdim, ierr)

         if (File%ndim > 0) then
!del         allocate (File%Axis(max_ndim))
             allocate (File%Axis(File%ndim))
             do i=1,File%ndim
                if (i == File%recdim) then
                   File%Axis(i) = read_ncaxis (File%id, i, .true.)
                else
                   File%Axis(i) = read_ncaxis (File%id, i, .false.)
                endif
             enddo
         endif

         File%nvar = nvar_total - File%ndim
         if (File%nvar > 0) then
!del         allocate (File%Var(max_nvar))
             allocate (File%Var(File%nvar))
             do i=1,File%nvar
                File%Var(i) = read_ncvar (File%id, i+File%ndim)
             enddo
         endif

         if (File%natt > 0) then
             allocate (File%Att(max_natt))
             do i=1,File%natt
                File%Att(i) = read_ncatt (File%id, NCGLOBAL, i)
             enddo
         endif

      end function read_ncfile

!-----------------------------------------------------------------------

      subroutine write_ncfile (File)

      type(ncfile_type) , intent(inout) :: File
      integer  i, nc, id, ierr, idum
      integer  nccre, ncsfil

         nc = len_trim(File%name)
!!!      File%id = nccre (File%name(1:nc), NCCLOB, ierr)
         id = nccre (File%name(1:nc), NCCLOB, ierr)
         File%id = id
         File%define_mode = .true.

         idum = ncsfil (File%id, NCNOFILL, ierr)

         do i=1,File%ndim
             call write_ncaxis (File%id, File%Axis(i))
         enddo

         do i=1,File%nvar
             call write_ncvar (File%id, File%Var(i))
         enddo

         do i=1,File%natt
             call write_ncatt (File%id, NCGLOBAL, File%Att(i))
         enddo

      end subroutine write_ncfile

!-----------------------------------------------------------------------

      function create_file (name) result (File)

      character(len=*), intent(in) :: name
      type(ncfile_type)            :: File
      integer  ierr, idum
      integer  nccre, ncsfil

           File%id = nccre (name, NCCLOB, ierr)
           idum = ncsfil (File%id, NCNOFILL, ierr)

           File%name    = name
           File%ndim    = 0
           File%nvar    = 0
           File%natt    = 0
           File%recdim  = 0
           File%nvarmap = 0

           File%define_mode = .true.

      end function create_file

!-----------------------------------------------------------------------

  subroutine close_ncfile (File)

   type(ncfile_type), intent(inout) :: File
   integer  ierr

!  ----- update length of record dimension  ------
      if (File%recdim > 0) then
        call close_ncaxis (File%id, File%recdim, File%Axis(File%recdim))
      endif

      call ncclos (File%id, ierr)

  end subroutine close_ncfile

!-----------------------------------------------------------------------

  subroutine empty_ncfile (File)

   type(ncfile_type), intent(inout) :: File

       if (associated(File%Axis  )) deallocate (File%Axis  )
       if (associated(File%Var   )) deallocate (File%Var   )
       if (associated(File%Att   )) deallocate (File%Att   )
       if (associated(File%Varmap)) deallocate (File%Varmap)

  end subroutine empty_ncfile

!-----------------------------------------------------------------------

      subroutine switch_define_mode (File)

      type(ncfile_type), intent(inout) :: File
      integer  ierr

        if (File%define_mode) then
            call ncendf (File%id, ierr)
            File%define_mode = .false.
        else
            call ncredf (File%id, ierr)
            File%define_mode = .true.
        endif

      end subroutine switch_define_mode

!-----------------------------------------------------------------------

   subroutine copy_ncfile (File_in, File_out)

     type(ncfile_type), intent(in)  :: File_in
     type(ncfile_type), intent(out) :: File_out

     integer :: i

!     ---- copy axes ----
      File_out%ndim = File_in%ndim
      allocate (File_out%Axis(File_out%ndim))
      do i = 1, File_out%ndim
         call copy_ncaxis ( File_in%Axis(i), File_out%Axis(i) )
      enddo

!     ---- copy variables ----
      File_out%nvar = File_in%nvar
      allocate (File_out%Var(File_out%nvar))
      do i = 1, File_out%nvar
         call copy_ncvar ( File_in%Var(i), File_out%Var(i) )
      enddo

!     ---- copy global attributes ----
      File_out%natt = File_in%natt
      allocate (File_out%Att(File_out%natt))
      do i = 1, File_out%natt
         call copy_ncatt ( File_in%Att(i), File_out%Att(i) )
      enddo

!     ---- copy variable mapping ----
      File_out%nvarmap = File_in%nvarmap
      allocate (File_out%Varmap(File_out%nvarmap))
      do i = 1, File_out%nvarmap
         call copy_varmap ( File_in%Varmap(i), File_out%Varmap(i) )
      enddo

!     ---- miscellaneous ----
      File_out%recdim = File_in%recdim
      File_out%define_mode = .true.

!     ---- un-copied data ----
!     File_out%id   = 
!     File_out%name = 

   end subroutine copy_ncfile

!-----------------------------------------------------------------------

end module ncfile_mod

