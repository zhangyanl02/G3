    program GBHM3
      use hydro_data_mod,only:subbasin,nsub,psubbasin,pbasinup,nbasinup,area,ele,slp,length,Ds,Dg,&
                         nflow,dx,Dr,s0,b,roughness,ngrid,grid_row,grid_col,Sstmax,surfn,&
                         kg,GWcs,nx,np,nc,inicon,MNgauge,Maxgauge,NNstation,nhrpdt
      use forcing_mod,only:rain_daily,pre_hour,tmin_daily,tmax_daily,evap_daily,snow
      use soil_para_mod,only:soiltyp,soil,nsoil,anik,wsat,wfld,wrsd,alpha,watern,ksat1,ksat2,ns,nlayer
      use land_para_mod,only:landtyp,land_ratio,land_use,nland,NDVI,yndvi,LAI,Kcanopy,LAImax,Kcrop,&
                              root,nv
      use global_para_mod,only:i4,r8,dt,year,month,day,hour,dayinmonth,hydroyear,startmonth,endmonth,startday,endday,&
                              idc,ihc,start,finish,start_sub,end_sub,nrow,ncol,startyear,endyear,iflai
      use GBHM3Lib_mod
      implicit none
      integer(i4)::isub,iflow,isoil,iland,year_tmp
      character*200::landuse            ! name of land use
      real(r8)::suction,suction2   ! soil suction
      
!     Parameter file & Input file
      character*100::grid_area     ! area of the grid (m2)
      character*100::top_ele       ! mean elevation of the grid (m)
      character*100::top_slp       ! mean slope of hillslope in the grid (ND)
      character*100::top_len       ! mean length of hillslope in the grid (m)
      character*100::soil_map      ! soil type of the grid (ND)
      !character*200 landuse_map   ! land use type of the grid (ND)
      character*100 ::geology_map   ! depth of the topsoil (m)
      integer(i4):: im             !,iy, id, ih ! for year, month, day and hour loop


      !common /location  /  para_dir,data_dir,result1_dir,	dem_dir,result2_dir,simul_dir
      character*200::  para_dir    ! directory for parameters
      character*200::  data_dir    ! directory for input data
      character*200::  dem_dir    ! directory for input data
      character*200::  result1_dir ! directory for storing simulation result
      character*200::  result2_dir ! directory for storing simulation result
      character*200::  simul_dir   ! directory for storing temporal variables
      character*200::  ndvi_dir,lai_dir ! added by gaobing 20131204
      character*200::  soil_dir  !! added by mqh 
      character*200::  infile
      character*5:: ncols,nrows
      character*9:: xllcorner,yllcorner
      character*8:: cellsize
      character*12:: nodata
      integer(i4):: nnr, nnc
      real(r8):: x0, y0
      real(r8):: gridsize
      real(r8):: znodata
      !common / asc / ncols,nrows,xllcorner,yllcorner,cellsize,nodata, nnr, nnc, x0, y0, gridsize,znodata
      
      character*2::    ch2, a2    ! a 2-character variable
      character*4::    ch4        ! a 4-character variable
      character*200:: atmp
      real::tmp, tmp1, tmp2   ! , value, mean  ! temporary variables
      integer:: i, j, k           ! temporary variables
      integer:: ir, ic, idd            !, iyy, imm
      integer:: ia1, ia2, ib1, ib2, ic1, ic2,ir2
      integer:: itmp ,tmp_num
      real:: tmp_rain,tmp_tmin,tmp_tmax,tmp_evap
      real:: ran0
      integer:: ii1,ig
      

      para_dir   ="/home/zhangyanlin/model/gbhm/model/model_test_data/hydropara/"
      data_dir   ="../data/"
      dem_dir   ="../dem/"
      result1_dir="../result_river/"
      result2_dir="../result_spatial/"
      simul_dir ="../simulation/"
      ndvi_dir ="../read_ndvi/"
      soil_dir ="../data/soil/"  !mqh
      lai_dir = "../lai/"
      dt=3600.0   ! time step of hydrological simulation, unit: second
      grid_area   = "cell_area.asc"
      top_ele     = "elevation.asc"
      top_len     = "slope_length.asc"
      top_slp     = "slope.asc"
      soil_map    = "soil_unit.asc"
      geology_map = "soil_depth.asc"
      
      call initial_subcatchment(para_dir)       !allocate memory for the variables before used
      call read_hydro_para(para_dir,nhrpdt)
      call read_veg_para(para_dir)
      
      start_sub = 1
      end_sub   = 153

!      parameter (nrow = 106, ncol = 78)
!      parameter (nc = 153, nx = 300, np = 500)
!  parameter (nv = 9,  ns = 3,  nlayer = 10)
!  parameter (startyear =2000, endyear =2000)
!  parameter (inicon = 0)
!  parameter (MNgauge = 10)
!  parameter (NNstation = 3)
!    parameter (Maxgauge = 19)
!  parameter (iflai = 0)     

    end program GBHM3

