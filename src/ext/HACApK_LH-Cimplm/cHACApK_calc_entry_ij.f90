!***cHACApK_entry_ij
function cHACApK_entry_ij(i,j,i_bemv) bind(c,name='cHACApK_entry_ij')
  use iso_c_binding
  use m_HACApK_calc_entry_record
  use m_HACApK_calc_entry_ij
  real(c_double) :: cHACApK_entry_ij
  integer(c_int),value :: i,j,i_bemv
  type(st_HACApK_calc_entry),pointer :: st_bemv
  call HACApK_calc_entry_record_get(i_bemv,st_bemv)
  cHACApK_entry_ij=HACApK_entry_ij(i,j,st_bemv)
end function cHACApK_entry_ij
