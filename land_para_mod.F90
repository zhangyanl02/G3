    module land_para_mod
      use global_para_mod
      public
      save
      integer(i4)::nland                        ! number of landuse types in whole area
      integer(i4)::nv                           ! number of landuse types
      real(r8),allocatable::    NDVI(:,:,:,:)    ! monthly NDVI  NDVI(nrow,ncol,12,3)
      real(r8),allocatable::    yndvi(:,:,:,:)   ! spot NDVI    yndvi(nrow,ncol,12,3)
      real(r8),allocatable::    LAI(:,:,:,:)     ! monthly LAI    LAI(nrow,ncol,12,4)
      real(r8),allocatable::    Kcanopy(:)           ! vegetation coverage of a landuse type   Kcanopy(nv)
      real(r8),allocatable::    LAImax(:)            ! maximum LAI in a year    LAImax(nv)
      real(r8),allocatable::    root(:)              ! root depth (m)   root(nv)
      integer(i4),allocatable:: landtyp(:)              ! landuse type of the whole area   landtyp(nv)
      real(r8),allocatable::    land_ratio(:,:,:)    ! Area fraction of each landuse type  land_ratio(nrow,ncol,nv)
      integer(i4),allocatable:: land_use(:,:,:)    !land_use(nrow*10,ncol*10,1)
      real(r8),allocatable::    Kcrop(:)             !evaporation coefficient of crop   Kcrop(nv)
    end module land_para_mod
