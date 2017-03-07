    module global_para_mod
      public
      save
      integer,parameter::i4=4
      integer,parameter::r8=4
      integer(i4)::nrow             ! number of rows of the whole area
      integer(i4)::ncol             ! number of columns of the whole area
      integer(i4)::year,month2,day,month
      integer(i4)::hour             ! hour in a day   (1-24)
      integer(i4)::iflai            ! 1 for using lai 

      integer(i4):: idc            ! contineous day in a year (1-366)
      integer(i4):: ihc            ! contineous hour in a year (1-366*24)
      integer(i4):: hydroyear      ! year for hydro simulation
      integer(i4):: startyear        ! simulation start
      integer(i4):: startmonth     ! start month in the hydro-year
      integer(i4):: startday       ! start day in the first month of the hydro-year
      integer(i4):: endyear          ! simulation end
      integer(i4):: endmonth       ! end   month in the hydro-year
      integer(i4):: endday         ! end   day in the last  month of the hydro-year
      integer(i4):: start, finish  ! 0 - faulse,    1 - true
      integer(i4):: start_sub, end_sub  ! for partial simulation of a basin
      real(r8)::    dt             ! time step (second)
      integer(i4):: dayinmonth(12) ! days of a month
      data dayinmonth /31,28,31,30,31,30,31,31,30,31,30,31/
      data startmonth /1/
      data startday   /1/
      data endmonth   /12/
      data endday     /31/
      character*200::para_dir    ! directory for parameters
      character*200::result2_dir ! directory for storing simulation result
      character*200::simul_dir   ! directory for storing temporal variables
      
      character*200::area_map     ! area of the grid (m2)
      character*200::ele_map       ! mean elevation of the grid (m)
      character*200::slp_map       ! mean slope of hillslope in the grid (ND)
      character*200::slplen_map       ! mean length of hillslope in the grid (m)
      character*200::soil_map      ! soil type of the grid (ND)
      !character*200 landuse_map   ! land use type of the grid (ND)
      character*200 ::ds_map   ! depth of the topsoil (m)
      
    end module global_para_mod
