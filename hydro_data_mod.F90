    module hydro_para_mod
      use global_para_mod
      public
      save
      integer(i4),allocatable::nflow(:)                      ! total number of flow-intervals in a sub-basin    nflow(nc)
      integer(i4),allocatable::ngrid(:,:)                    ! the number of grid in each flow interval                    ----ngrid(nsub,nflow)
      integer(i4),allocatable::grid_row(:,:,:)               ! the row number of each grid                                 ----grid_row(nsub,nflow,ngrid)
      integer(i4),allocatable::grid_col(:,:,:)               ! the col number of each grid                                 ----grid_col(nsub,nflow,ngrid)
      character*6,allocatable::subbasin(:)                  ! name of each sub catchment                                         ----subbasin(nc)
      integer(i4),allocatable::psubbasin(:)                 ! the ids of sub catchments,e.g., 1001,2003                          ----psubbasin
      integer(i4),allocatable::nbasinup(:,:)                ! the rank of each up-basin of each sub catchment, e.g., 1,2,3,4     ----nbasinup(nc,8)
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
    end module hydro_para_mod
