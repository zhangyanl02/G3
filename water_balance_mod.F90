    module water_balance_mod
      use global_para_mod,only:r8,i4
      public
      save
      real(r8):: annual_rain        ! annual precipitation 
      real(r8):: annual_runoff      ! annual runoff
      real(r8):: annual_Eact        ! annual actual evaporation
      real(r8):: annual_Ecy
      real(r8):: annual_Ecp
      real(r8):: annual_Ese
      real(r8):: annual_Cst         ! annual canopy interception
      real(r8):: annual_Ssts        ! annual surface storage(snow)
      real(r8):: annual_Sstr        ! annual surface storage(rain)
      real(r8):: annual_SBst        ! annual subsurface storage
      real(r8):: annual_Gst         ! annual groundwater storage
      real(r8):: annual_Cst0        ! initial value of annual canopy storage
      real(r8):: annual_Ssts0       ! initial value of annual surface storage(snow)
      real(r8):: annual_Sstr0       ! initial value of annual surface storage(rain)
      real(r8):: annual_SBst0       ! initial value of annual subsurface storage
      real(r8):: annual_Gst0        ! initial value of groundwater storage
      
      real(r8)::     total_rain_d(366)            ! basin mean precipitation
      real(r8)::     total_Ep_d(366)                  ! basin mean potential evaporation
      real(r8)::     total_Eact_d(366)            ! basin mean actual evaporation
      real(r8)::     total_runoff_d(366)            ! basin mean runoff
      
      real(r8),allocatable:: evap_month(:,:)         ! monthly mean evaporation    nrow,ncol
      real(r8),allocatable:: soil_month(:,:)         !nrow,ncol
      real(r8),allocatable:: runoff_month(:,:)       !nrow,ncol
      real(r8),allocatable:: grunoff_month(:,:)      !nrow,ncol
      real(r8),allocatable:: drunoff_month(:,:)      !nrow,ncol
      real(r8),allocatable:: srunoff_month(:,:)      !nrow,ncol
      real(r8),allocatable:: ecy_month(:,:)      !nrow,ncol
      real(r8),allocatable:: ecp_month(:,:)      !nrow,ncol
      real(r8),allocatable:: ese_month(:,:)      !nrow,ncol
      
      real(r8),allocatable::ecy(:,:,:)           !Evaporation from canopy   nrow,ncol,366
      real(r8),allocatable::ecp(:,:,:)           !Evaporation from crop !nrow,ncol,366
      real(r8),allocatable::ese(:,:,:)           !Evaporation from soil surface !nrow,ncol,366
      real(r8),allocatable::eactd(:,:,:)         !daily actual evapotranspiration (mm) actual evaopration nrow,ncol,366
      real(r8),allocatable::epd(:,:,:)         ! daily potential evaporation (mm)   !nrow,ncol,366
    end module water_balance_mod
