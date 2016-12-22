    module global_para_mod
      public
      save
      integer,parameter::i4=4
      integer,parameter::r8=4
      integer(i4)::nrow             ! number of rows of the whole area
      integer(i4)::ncol             ! number of columns of the whole area
      integer(i4)::year,month2,day,hour,month
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
      
      character*200::para_dir    ! directory for parameters
      character*200::result2_dir ! directory for storing simulation result
      character*200::simul_dir   ! directory for storing temporal variables
    end module global_para_mod
