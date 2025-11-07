module m_HACApK_calc_entry_record
  use m_HACApK_calc_entry_ij

  private

  public:: HACApK_calc_entry_record_set
  public:: HACApK_calc_entry_record_get
  public:: HACApK_calc_entry_record_clear

!*** type :: st_HACApK_calc_entry_record
  type :: st_HACApK_calc_entry_record
    logical :: used = .false.
    type(st_HACApK_calc_entry), pointer :: st_bemv
  end type st_HACApK_calc_entry_record

  integer, parameter :: MAX_BEMV = 8

  type(st_HACApK_calc_entry_record), save :: st_bemvrecord(MAX_BEMV)

contains

!***HACApK_calc_entry_record_set
  subroutine HACApK_calc_entry_record_set(st_bemv, id)
    implicit none
    type(st_HACApK_calc_entry), intent(in), target :: st_bemv
    integer, intent(out) :: id
    integer :: i
    id = 0
    do i = 1, MAX_BEMV
      if (.not. st_bemvrecord(i)%used) then
        id = i
        exit
      endif
    end do
    if (id == 0) then
      print*,'Error: HACApK_calc_entry_record_set: too many st_bemv set'
      stop
    endif
    st_bemvrecord(id)%st_bemv => st_bemv
    st_bemvrecord(id)%used = .true.
  end subroutine HACApK_calc_entry_record_set

!***HACApK_calc_entry_record_get
  subroutine HACApK_calc_entry_record_get(id, st_bemv)
    implicit none
    integer, intent(in) :: id
    type(st_HACApK_calc_entry), pointer :: st_bemv
    if (id <= 0 .or. MAX_BEMV < id) then
      print*,'Error: HACApK_calc_entry_record_get: id out of range'
      stop
    endif
    if (.not. st_bemvrecord(id)%used) then
      print*,'Error: HACApK_calc_entry_record_get: invalid id'
      stop
    endif
    st_bemv => st_bemvrecord(id)%st_bemv
  end subroutine HACApK_calc_entry_record_get

!***HACApK_calc_entry_record_clear
  subroutine HACApK_calc_entry_record_clear(id)
    implicit none
    integer, intent(in) :: id
    if (id <= 0 .or. MAX_BEMV < id) then
      print*,'Error: HACApK_calc_entry_record_clear: id out of range'
      stop
    endif
    if (.not. st_bemvrecord(id)%used) then
      print*,'Error: HACApK_calc_entry_record_clear: invalid id'
      stop
    endif
    st_bemvrecord(id)%st_bemv => null()
    st_bemvrecord(id)%used = .false.
  end subroutine HACApK_calc_entry_record_clear

end module m_HACApK_calc_entry_record
