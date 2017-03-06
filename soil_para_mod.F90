    module soil_para_mod
      use global_para_mod
      public
      save
      integer(i4)::ns                              ! number of soil types
      integer(i4)::nlayer
      integer(i4):: nsoil                          ! number of soil types in whole area
      integer(i4),allocatable::soiltyp(:)          ! Soil type of the whole area       soiltyp(ns)
      integer(i4),allocatable:: soil(:, :)         ! soil code of each grid  soil(nrow, ncol)
      real(r8),allocatable::    wsat(:,:)          ! saturated soil moisture   nrow,ncol
      real(r8),allocatable::    wrsd(:,:)          ! residual soil moisture    nrow,ncol
      real(r8),allocatable::    wfld(:,:)    ! soil moisture at field capacity    nrow,ncol
      real(r8),allocatable::    alpha(:,:)   ! soil water parameter     nrow,ncol
      real(r8),allocatable::    watern(:,:)  ! soil water parameter    nrow,ncol
      real(r8),allocatable::    ksat1(:,:)        ! saturated hydraulic conductivity at surface (mm/h)     nrow,ncol
      real(r8),allocatable::    ksat2(:,:)     ! saturated hydraulic conductivity at bottom (mm/h)     nrow,ncol
      real(r8),allocatable::    anik(:)              ! soil anisotropy ratio anik(nv)
    end module soil_para_mod
