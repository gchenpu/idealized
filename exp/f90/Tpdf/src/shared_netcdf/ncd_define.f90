
module ncd_define_mod

!------------------ flags for ncopen and nccreate ----------------------

   integer, parameter :: NCNOWRIT =  0, NCWRITE =  1,  &
                         NCNOCLOB =  4, NCCLOB  =  0
!!!vers2                 NCNOCLOB = 15, NCCLOB  = 11

   integer, parameter :: NCNOFILL = 256, NCFILL = 0

   integer, parameter :: NCGLOBAL = 0

   integer, parameter :: MAXVDIMS = 100

!------------------- valid data types ----------------------------------

   integer, parameter :: NCBYTE = 1, NCCHAR  = 2, NCSHORT  = 3,  &
                         NCLONG = 4, NCFLOAT = 5, NCDOUBLE = 6

   integer, parameter ::   short_type = selected_int_kind(4)
   integer, parameter ::    long_type = selected_int_kind(8)

   integer, parameter ::  float_type = selected_real_kind( 6, 37)
   integer, parameter :: double_type = selected_real_kind(15,307)

!------------------- maximum character string lengths ------------------

   integer, parameter :: maxn = 128
   integer, parameter :: maxc = 256

!------------- maximum number of dimensions per variable ---------------

   integer, parameter :: max_dims = 4

!--------- maximum number of dimensions per file -----------------------
!--------- maximum number of variables per file ------------------------
!--------- maximum number of attributes per variable -------------------

              integer :: max_ndim = 12
              integer :: max_nvar = 200
              integer :: max_natt = 10

!----------------- default fill values ---------------------------------

   integer, parameter :: fill_byte   = 129, fill_char   = 0
   integer, parameter :: fill_short  = -32767
   integer, parameter :: fill_long   = -2147483647
   real ,   parameter :: fill_float  = 9.9692099683868690e+36
   real ,   parameter :: fill_double = 9.9692099683868690e+36

!-----------------------------------------------------------------------

end module ncd_define_mod

