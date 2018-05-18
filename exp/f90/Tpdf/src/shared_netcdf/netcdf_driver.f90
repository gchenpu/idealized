
   module netcdf_driver_mod

!---------------------- public interfaces ------------------------------

  use  netcdf_file_mod, only:  define_axis,  define_time_axis,  &
                         define_variable, define_variable_tave,  &
                         define_tave,  define_file,             &
                         read_file_header,  write_file_header,  &
                         copy_file_header, get_axis_units,   &
                         get_axis_long_name, get_axis_cart_name,  &
                         get_axis_edges_name,               &
                         get_var_units, get_var_long_name,  &
                         get_missing_value, get_axis_direction,  &
                         get_calendar_type, get_valid_range

  use  ncread_write_mod, only:  write_time_coord, write_variable_tave, &
                                write_variable, read_time_coord,  &
                                read_variable, read_variable_tave, &
                                axis_values, time_axis_values,  &
                                variable_tave_info

  use  ncfile_access_mod, only:  file_name, num_axis, num_var,       &
                                 axis_name, axis_index, axis_length, &
                                 time_axis_index, time_axis_length,  &
                                 var_name, var_index, num_var_axis,  &
                                 var_axis_len, var_axis

  use       ncfile_mod, only:  ncfile_type, close_ncfile, empty_ncfile

  use   ncd_define_mod, only:  float_type, double_type

!-----------------------------------------------------------------------

  end module netcdf_driver_mod

