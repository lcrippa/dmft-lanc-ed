  if(allocated(first_state_up))deallocate(first_state_up)
  if(allocated(last_state_up))deallocate(last_state_up)
  allocate(first_state_up(DimDw),last_state_up(DimDw))
  first_state_dw = huge(1)
  last_state_dw  = 0
  !
  first_state_up = huge(1)
  last_state_up  = 0
  !
  do i=first_state,last_state
     idw=idw_index(i,DimUp)
     if(idw < first_state_dw) first_state_dw = idw
     if(idw > last_state_dw ) last_state_dw  = idw
     !
     iup=iup_index(i,DimUp)
     if(iup < first_state_up(idw)) first_state_up(idw) = iup
     if(iup > last_state_up(idw) ) last_state_up(idw)  = iup
     !
  enddo

