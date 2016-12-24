    module hydro_data_mod
      use global_para_mod
      public
      save
      integer(i4)::nsub             ! number of sub-catchments  (nc in dims2.inc)
      integer(i4)::nx               ! number of flow interval in a sub-catchment
      integer(i4)::np               ! number of grids in a flow-interval
      integer(i4)::nlayer           ! number of UZ layers
      integer(i4)::inicon           ! the way to specify initial condition
                                    !   inicon=1, input initial condition from files in directory ../simulation
                                    !   inicon=0, arbitrarily given in the program
      integer(i4)::MNgauge          ! Number of Gauge          old
      integer(i4)::Maxgauge         ! Max_Number of Gauge      old +new
      integer(i4)::NNstation        ! rain station select; 1 for old, 2 for new, 3 for old+new

      real(r8)::total_rech,total_qsub,total_ssub
      
      integer(i4),allocatable::nflow(:)                      ! total number of flow-intervals in a sub-basin    nflow(nc)
      integer(i4),allocatable::ngrid(:,:)                    ! the number of grid in each flow interval                    ----ngrid(nsub,nflow)
      integer(i4),allocatable::grid_row(:,:,:)               ! the row number of each grid                                 ----grid_row(nsub,nflow,ngrid)
      integer(i4),allocatable::grid_col(:,:,:)               ! the col number of each grid                                 ----grid_col(nsub,nflow,ngrid)
      character*6,allocatable::subbasin(:)                  ! name of each sub catchment                                         ----subbasin(nc)
      integer(i4),allocatable::psubbasin(:)                 ! the ids of sub catchments,e.g., 1001,2003                          ----psubbasin
      integer(i4),allocatable::nbasinup(:,:)                ! the rank of each up-basin of each sub catchment, e.g., 1,2,3,4     ----nbasinup(nc,8)
      integer(i4),allocatable::inbasin(:,:)                ! the rank of each up-basin of each sub catchment, e.g., 1,2,3,4     ----nbasinup(nc,8)
      
      real(r8),allocatable::   area(:,:)                    ! area of the local grid (m2)------               area(nrow,ncol)
      real(r8),allocatable::   slp(:,:)               ! average slope of the local grid (ND)            slp(nrow,ncol)
      real(r8),allocatable::   length(:,:)   ! average hillslope length (m)        length(nrow,ncol)
      real(r8),allocatable::   Ds(:,:)       ! average depth of the topsoil (m)    Ds(nrow,ncol)
      real(r8),allocatable::   Dg(:,:)       ! average depth of the uncinfined acwuifer (m)     Dg(nrow,ncol)
      
      real(r8),allocatable::   dx(:,:)                          ! (rdx) length of flow intervals (m)                          ----dx(nsub,nflow)
      real(r8),allocatable::   dr(:,:)                          ! (drr) river depth of the flow interval (m)                  ----dr(nsub,nflow)
      real(r8),allocatable::   s0(:,:)                          ! river-bed slope of the flow interval                        ----s0(nsub,nflow)
      real(r8),allocatable::   b(:,:)                           ! (riverb) river width of the flow interval (m)               ----b(nsub,nflow)
      real(r8),allocatable::   roughness(:,:)                   ! river roughness (manning's coefficient) of the flow interval----(nsub,nflow)
      
      real(r8),allocatable::   subarea(:)    ! total basin area     nc
      real(r8)::basinarea,fice, area_frac
      
      real(r8),allocatable::    Sstmax(:)            ! maximum surface water detension (mm)     nv
      real(r8),allocatable::    surfn(:)             ! surface roughness (Manning's coefficient)     nv
      real(r8),allocatable::    kcrop(:)             ! evaporation coefficient of crop     nv
      
      real(r8),allocatable::    Qh(:,:)     ! hourly mean discharge (m3/s)     nc,8800
      integer(i4),allocatable:: layer(:,:)   ! number of UZ layer     nrow,ncol
      real(r8),allocatable::    D(:,:,:)    ! depth of each UZ layer(m)    nrow,ncol,nlayer
      real(r8),allocatable::    k0(:,:,:)   ! saturated hydraulic conductivity (mm/hr)   nrow,ncol,nlayer
      real(r8),allocatable::    w(:,:,:,:) ! soil moisture     nrow,ncol,nv,nlayer
      real(r8),allocatable::    qin(:)             ! lateral inflow into river    nflow
      real(r8),allocatable::    kg(:,:)          ! hydraulic conductivity of groundwater (m/sec)    nrow,ncol
      real(r8),allocatable::    GWcs(:,:)            ! groundwater storage coefficient       nrow,ncol
      real(r8),allocatable::    Cst(:,:,:)           ! canopy storage     nrow,ncol,nv
      real(r8),allocatable::    Sst(:,:,:)           ! surface storage   nrow,ncol,nv
      real(r8),allocatable::    Dgl(:,:,:)           ! depth to groundwater level in the grid    nrow,ncol,nv
      real(r8),allocatable::    GWst(:,:,:)          !   nrow,ncol,nv
      real(r8),allocatable::    Drw(:,:)          ! river water depth   nsub,nflow
      real(r8),allocatable::    discharge(:,:)    ! river flow discharge   nsub,nflow      



      

      
      integer(i4),allocatable::countt2(:)   !nsub
      
    end module hydro_data_mod
