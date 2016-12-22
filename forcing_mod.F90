    module forcing_mod
    public
    save
      real(r8),allocatable::    raind(:,:,:)       ! daily rainfall (mm)   nrow,ncol,366
      real(r8),allocatable:: rain_daily(:,:,:)  ! daily precipitaiton (mm)     rain_daily(nrow,ncol,31)
      real(r8),allocatable:: pre_hour(:,:,:,:)                  !pre_hour(nrow,ncol,31,24)
      real(r8),allocatable:: tmin_daily(:,:,:)  ! daily minimum air temperature (degree)         tmin_daily(nrow,ncol,31)
      real(r8),allocatable:: tmax_daily(:,:,:)  ! daily minimum air temperature (degree)    tmax_daily(nrow,ncol,31)
      real(r8),allocatable:: evap_daily(:,:,:,:) !daily value of potential evaporation(mm)    evap_daily(nrow,ncol,31,nv)
      real(r8),allocatable:: snow(:,:,:)        ! snow depth (mm water) of the initial year      snow(nrow,ncol,nv)
      
      
    end module forcing_mod
