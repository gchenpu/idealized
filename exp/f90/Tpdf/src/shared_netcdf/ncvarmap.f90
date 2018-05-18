
module ncvarmap_mod

use ncd_define_mod

implicit none
private

public :: get_varmap, put_varmap, copy_varmap
public :: get_varmap_no, get_varmap_tavg
public :: ncvarmap_type

type ncvarmap_type
    character(len=maxn) :: name
    integer             :: varno, tavg_varno(3)
    real                :: scale, offset
    logical             :: pack
end type ncvarmap_type

CONTAINS

!-----------------------------------------------------------------------

function get_varmap_no (Varmap, name) result (no)

   type(ncvarmap_type), intent(in), dimension(:) :: Varmap
   character(len=*),    intent(in)               :: name
   integer                                     :: no
   integer  i

   no = 0

   do i = 1, size(Varmap)
      if (trim(Varmap(i)%name) == trim(name)) then
          no = i
          exit
      endif
   enddo

end function get_varmap_no

!-----------------------------------------------------------------------

subroutine get_varmap (Varmap, nvarmap, name, varno,   &
                       pack, scale, offset)

   type(ncvarmap_type), intent(in), dimension(:) :: Varmap
   integer          ,   intent(in)               :: nvarmap
   character(len=*),    intent(in)               :: name
   integer         ,    intent(out)              :: varno
   logical         ,    intent(out)              :: pack
   real            ,    intent(out)              :: scale, offset
   integer   no

   no = get_varmap_no (Varmap(1:nvarmap), name)

   if (no > 0) then
       varno = Varmap(no)%varno
       pack  = Varmap(no)%pack
       scale = Varmap(no)%scale
       offset = Varmap(no)%offset
   else
       varno = 0
       pack  = .false.
       scale = 1.0; offset = 0.0
   endif

end subroutine get_varmap

!-----------------------------------------------------------------------

subroutine get_varmap_tavg (Varmap, nvarmap, name, tavg_varno)

   type(ncvarmap_type), intent(in), dimension(:) :: Varmap
   integer            , intent(in)               :: nvarmap
   character(len=*),    intent(in)               :: name
   integer         ,    intent(out)              :: tavg_varno(3)
   integer   no

   tavg_varno = 0
   no = get_varmap_no (Varmap(1:nvarmap), name)

   if (no > 0) tavg_varno = Varmap(no)%tavg_varno

end subroutine get_varmap_tavg

!-----------------------------------------------------------------------

subroutine put_varmap (Varmap, nvarmap, name, varno, tavg_varno,  &
                       pack, scale, offset)

   type(ncvarmap_type), intent(inout), dimension(:) :: Varmap
   integer            , intent(inout)               :: nvarmap
   character(len=*),  intent(in)                 :: name
   integer         ,  intent(in)                 :: varno, tavg_varno(3)
   logical         ,  intent(in), optional       :: pack
   real            ,  intent(in), optional       :: scale, offset
   integer   no

   no = get_varmap_no (Varmap(1:nvarmap), name)

!  ---- replace varmap value ? ----
   if (no > 0) then
       Varmap(no)%name       = name
       Varmap(no)%varno      = varno
       Varmap(no)%tavg_varno = tavg_varno
       if (present(pack)) then
           Varmap(no)%pack = pack
       else
           Varmap(no)%pack = .false.
       endif
       if (present(scale) .and. present(offset)) then
           Varmap(no)%scale  = scale
           Varmap(no)%offset = offset
       else
           Varmap(no)%scale  = 1.0
           Varmap(no)%offset = 0.0
       endif
       return
   endif

!  ---- add new varmap value -----

   if (nvarmap+1 <= max_nvar) then
       nvarmap = nvarmap+1
       no = nvarmap
       Varmap(no)%name       = name
       Varmap(no)%varno      = varno
       Varmap(no)%tavg_varno = tavg_varno
       if (present(pack)) then
           Varmap(no)%pack = pack
       else
           Varmap(no)%pack = .false.
       endif
       if (present(scale) .and. present(offset)) then
           Varmap(no)%scale  = scale
           Varmap(no)%offset = offset
       else
           Varmap(no)%scale  = 1.0
           Varmap(no)%offset = 0.0
       endif
   else
       print *, 'ERROR nvarmap exceeds max_nvar'
       stop 111
   endif

end subroutine put_varmap

!-----------------------------------------------------------------------

subroutine copy_varmap (Varmap_in, Varmap_out)

   type(ncvarmap_type), intent(in)  :: Varmap_in
   type(ncvarmap_type), intent(out) :: Varmap_out

   Varmap_out%name       = Varmap_in%name
   Varmap_out%varno      = Varmap_in%varno
   Varmap_out%tavg_varno = Varmap_in%tavg_varno
   Varmap_out%scale      = Varmap_in%scale
   Varmap_out%offset     = Varmap_in%offset
   Varmap_out%pack       = Varmap_in%pack

end subroutine copy_varmap

!-----------------------------------------------------------------------

end module ncvarmap_mod

