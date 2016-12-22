    module water_balance_mod
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
      
      
    end module water_balance_mod
