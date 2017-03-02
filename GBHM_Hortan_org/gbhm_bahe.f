c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######          GBHM2 (The GBHM model of Version 2)         ######
c     ######                       (GBHM2)                        ######
c     ######                                                      ######
c     ######             Developed by Dr. Dawen YANG              ######
c     ######           Department of Civil Engineering            ######
c     ######                University of Tokyo                   ######
c     ######                                                      ######
c     ######                  修改于 2006.5.17                    ######
c     ##################################################################
c     ##################################################################
c
      PROGRAM GBHM2
c
c#######################################################################
c
c   PURPOSE: Hydrological simulation of a large river basin
c
c        (1) Using Horton Basin Numbering scheme             
c        (2) Subgrid parameterization: Flow-interval & Hillslope-river Concept
c        (3) Physically-based hydrological simulation on hillslope       
c        (4) Kinematic wave for river network routing           
c        (5) Crop iirigation: dry crop and paddy  
c
c#######################################################################
c
c     AUTHOR: Dr. Dawen YANG
c     Version 1: 08/2000 
c                                                                               
c     Version 2: 04/2002
c
c#######################################################################
c
c
c#######################################################################
c
c     Variable Declarations.
c
c#######################################################################
c
      implicit none
      include 'dims2.inc'
c
c#######################################################################
c
c     Sub-basin parameters
c
c#######################################################################
c	character*200 basin             ! name of simulated basin
	character*6   subbasin(nc)      ! name of sub-basin ! changed by gaobing 2013
	integer       nsub              ! total number of sub-basins
	integer       isub              ! sub-basin number
	integer       psubbasin(nc)     ! subbasin p number ! added by gaobing 2013
      integer       pbasinup(nc,4)    ! upbasin of p number ! added by gaobing 2013
	integer       nbasinup(nc,4)    ! upbasin of n number ! added by gaobing 2013
c 
c#######################################################################
c
c     Topographic parameters
c
c#######################################################################
c
	real area(nrow,ncol)     ! area of the local grid (m2)
	real ele(nrow,ncol)      ! elevation of the local grid (m)
	real slp(nrow,ncol)      ! average slope of the local grid (ND)
	real length(nrow,ncol)   ! average hillslope length (m)
	real Ds(nrow,ncol)       ! average depth of the topsoil (m)
	real Dg(nrow,ncol)       ! average depth of the unconfined acquifer (m)

	integer nflow(nc)        ! total number of flow-intervals in a sub-basin
	integer iflow            ! flow-interval number
	real    dx(nc,nx)        ! flow interval length (m)
      real    Dr(nc,nx)        ! river depth in the flow interval (m)
	real    s0(nc,nx)        ! slope of river bed (ND)
	real    b(nc,nx)         ! width of river (m)
	real    roughness(nc,nx) ! Manning's roughness

	integer ngrid(nc,nx)       ! number of grids in this flow-interval
	integer grid_row(nc,nx,np) ! row of grids in this flow-interval
	integer grid_col(nc,nx,np) ! column of grids in this flow-interval
c
c#######################################################################
c
c     Atmospheric forcing data
c
c#######################################################################
c
	real rain_daily(nrow,ncol,31)  ! daily precipitation (mm)
	real pre_hour(nrow,ncol,31,24)  ! hourly precipitation (mm)
	real tmin_daily(nrow,ncol,31)  ! daily minimum air temperature (degree)
	real tmax_daily(nrow,ncol,31)  ! daily minimum air temperature (degree)
	real evap_daily(nrow,ncol,31,nv)! daily value of potential evaporation(mm)
	real snow(nrow,ncol,nv)    ! snow depth(mm water) of the initial year
c
c#######################################################################
c
c     Soil-landuse distribution
c
c#######################################################################
c
c	character*200 soilname			 ! name of soil
	character*200 landuse            ! name of land use
      integer soiltyp(ns)              ! Soil type of the whole area
      integer landtyp(nv)              ! landuse type of the whole area
	integer soil(nrow, ncol)         ! soil code of each grid
      real land_ratio(nrow,ncol,nv)    ! Area fraction of each landuse type
	integer land_use(nrow*10,ncol*10,1)   !   added by xujijun
	integer nsoil                    ! number of soil types in whole area
	integer isoil                    ! soil number
	integer nland                    ! number of landuse types in whole area
	integer iland                    ! land-use number

	real    NDVI(nrow,ncol,12,3)    ! monthly NDVI
	real    yndvi(nrow,ncol,12,3)   ! spot NDVI ! 
	real    LAI(nrow,ncol,12,4)      ! fan_LAI added by gaobing 20131204
	real    Kcanopy(nv)           ! vegetation coverage of a landuse type
	real    LAImax(nv)            ! maximum LAI in a year
	real 	Kcrop(nv)             !evaporation coefficient of crop
      real    root(nv)              ! root depth (m)
	real    anik(nv)              ! soil anisotropy ratio 
      real    Sstmax(nv)	          ! maximum surface water detension (mm)
	real    surfn(nv)             ! surface roughness (Manning's coefficient)

c	real    wsat(ns)	        ! saturated soil moisture
c	real    wrsd(ns)		    ! residual soil moisture 
c	real    wfld(ns)            ! soil moisture at field capacity
c	real    alpha(ns)	        ! soil water parameter
c	real    watern(ns)		    ! soil water parameter
c	real    ksat1(ns)	  ! saturated hydraulic conductivity at surface (mm/h)
c      real    ksat2(ns)     ! saturated hydraulic conductivity at bottom (mm/h)
c	real    kg(ns)              ! hydraulic conductivity of groundwater (mm/h)
c      real    GWcs(ns)            ! groundwater storage coefficient
      real    wsat(nrow,ncol)     ! added by gaobing 20131125
      real    wfld(nrow,ncol)     ! added by gaobing 20131125
      real    wrsd(nrow,ncol)     ! added by gaobing 20131125
      real    alpha(nrow,ncol)    ! added by gaobing 20131125
	real    watern(nrow,ncol)   ! added by gaobing 20131125
	real    ksat1(nrow,ncol)    ! added by gaobing 20131125
	real    ksat2(nrow,ncol)    ! added by gaobing 20131125
      real    kg(nrow,ncol)       ! added by gaobing 20131125
	real    GWcs(nrow,ncol)     ! added by gaobing 20131125  
      real    suction,suction2             ! soil suction
c	real    w                   ! soil moisture
c
c#######################################################################
c
c     Parameter file & Input file
c
c#######################################################################
c
	character*100 grid_area     ! area of the grid (m2)
	character*100 top_ele       ! mean elevation of the grid (m)
	character*100 top_slp       ! mean slope of hillslope in the grid (ND)
	character*100 top_len       ! mean length of hillslope in the grid (m)
	character*100 soil_map      ! soil type of the grid (ND)
c	character*200 landuse_map   ! land use type of the grid (ND)
	character*100 geology_map   ! depth of the topsoil (m)
c
c#######################################################################
c
c     Date and time variables
c
c#######################################################################
c
      real    dt             ! time step (second)
	integer year, year_tmp          ! calendar year   (19XX)
	integer month          ! month in a year (1-12)
	integer day            ! day in a month  (1-30)
	integer hour           ! hour in a day   (1-24)
	integer dayinmonth(12) ! days of a month
	integer hydroyear      ! year for hydro simulation
	integer startmonth     ! start month in the hydro-year
	integer endmonth	   ! end   month in the hydro-year
	integer startday       ! start day in the first month of the hydro-year
	integer endday         ! end   day in the last  month of the hydro-year
	integer im             !,iy, id, ih ! for year, month, day and hour loop
	integer idc            ! contineous day in a year (1-366)
	integer ihc            ! contineous hour in a year (1-366*24)

	integer start, finish       ! 0 - faulse,    1 - true
	integer start_sub, end_sub  ! for partial simulation of a basin
c
c#######################################################################
c
c      Pfafstetter numbering system
c
c#######################################################################
c
	integer kfs1(10)         ! key of forward subdivision of level 1
	integer kfs2(10,10)      ! key of forward subdivision of level 2
	integer kfs3(10,10,10)   ! key of forward subdivision of level 3
	integer nlevel          ! total level of subdivision in Pfafstetter scheme
	integer l1, l2, l3
      integer level
	integer p1(9), p2(9), p3(9) ! Pfaftetter basin number of level 1, 2 and 3
c
c#######################################################################
c
c     Other variables
c
c#######################################################################
c
      character*200  para_dir    ! directory for parameters
	character*200  data_dir    ! directory for input data
	character*200  dem_dir    ! directory for input data
	character*200  result1_dir ! directory for storing simulation result
	character*200  result2_dir ! directory for storing simulation result
	character*200  simul_dir   ! directory for storing temporal variables
	character*200  ndvi_dir,lai_dir ! added by gaobing 20131204
	character*200  soil_dir  !! added by mqh 
	character*200  infile

      character*5 ncols,nrows
      character*9 xllcorner,yllcorner
      character*8 cellsize
      character*12 nodata
	integer nnr, nnc
	real x0, y0
	real gridsize
	real znodata
   	character*2    ch2, a2    ! a 2-character variable
	character*4    ch4        ! a 4-character variable

      character*200 atmp
      real    tmp, tmp1, tmp2   ! , value, mean  ! temporary variables
	integer i, j, k           ! temporary variables
c	integer ii, jj
	integer ir, ic, idd            !, iyy, imm
c	real    basinarea		  ! total basin area
c	real    areacheck		  ! for checking basin area
	integer ia1, ia2, ib1, ib2, ic1, ic2,ir2
	integer itmp ,tmp_num
	real tmp_rain,tmp_tmin,tmp_tmax,tmp_evap
	real ran0

cc	real irr_zone(9),irr_area(nc,nx)
cc	integer ig
cc	real re_capacity1(nc,nx) !capacity of reservoir in location(isub,iflow),m3
cc	real re_capacity0(nc,nx)  ! initial capacity of reservoir 
c	
c#######################################################################
c
c     Common blocks.
c
c#######################################################################
c
     	common /date1  / year, month, day, hour, idc,ihc
     	common /date2  / hydroyear,startmonth,startday,endmonth,endday
     	common /date3  / dayinmonth
	common /simulation/  start, finish, dt, start_sub, end_sub
	common /river1    /  nsub, subbasin
	common /river2    /  nflow, dx, Dr
	common /river3    /  s0, b, roughness
	common /river4    /  psubbasin,nbasinup ! added by gaobing 2013

cc	common /reservoir/ re_capacity1,re_capacity0  

	common / weather / rain_daily,tmin_daily, 
     :         	tmax_daily,evap_daily,snow,pre_hour
	common /LAID   / NDVI,yndvi,LAI

	common /grid_attrib/ ngrid,grid_row,grid_col
      common /slope     /  area, length, slp, Ds, Dg
	common /soil_land /  nsoil,soiltyp,nland,
     :                	landtyp,land_ratio,soil,land_use
	common /veg_para / Kcanopy,LAImax,root,anik
     :                  ,Sstmax,surfn,Kcrop
	common /soil_para /  wsat,wfld,wrsd,alpha,
     :                    watern,ksat1,ksat2,kg,GWcs
	common /location  /  para_dir,data_dir,result1_dir,	dem_dir,
     :                     result2_dir,simul_dir
      ! 标记

	common /pfaf/p1,p2,p3          !Pfaftetter basin number


	common / asc / ncols,nrows,xllcorner,yllcorner,cellsize,
     :           nodata, nnr, nnc, x0, y0, gridsize,znodata

      real basinarea2,subarea2(nc)  !added by mqh
      integer ii1,ig

c
c
c#######################################################################
c
c	To define partial simulation of this basin
c
c#######################################################################
c
cc	 修改于2013.05.09  gaobing								! 书签1
	start_sub = 1
	end_sub   = 153  ! total is 124  !mqh
	
c
c#######################################################################
c
c     Initialize the model parameters. Most of them are global parameters 
c     passed through common blocks.
c
c#######################################################################
c
	data dayinmonth /31,28,31,30,31,30,31,31,30,31,30,31/
	data startmonth /1 /
	data startday   /1 /
	data endmonth   /12 /
	data endday     /31/ 
c
	data kfs1/10*0/
	data kfs2/100*0/
	data kfs3/1000*0/

	
	
	
c
cc	para_dir   ="../parameter-2/"

	
	
	
c
cc	para_dir   ="../parameter-2/"

	para_dir   ="../parameter/"
	data_dir   ="../data/"
      dem_dir   ="../dem/"
	result1_dir="../result_river/"
	result2_dir="../result_spatial/"
	simul_dir ="../simulation/"
	ndvi_dir ="../read_ndvi/"
	soil_dir ="../data/soil/"  !mqh
	lai_dir = "../lai/"
	!标记
c
	dt=3600.0   ! time step of hydrological simulation, unit: second
c
	grid_area   = "cell_area.asc"
	top_ele     = "elevation.asc"
	top_len     = "slope_length.asc"
	top_slp     = "slope.asc"
	soil_map    = "soil_unit.asc"
c	landuse_map
	geology_map = "soil_depth.asc"
c
c#######################################################################
c
c           generate subcatchment name from kfs.dat
c                       Pfafstetter Scheme
c
c#######################################################################
c
	call strlen(para_dir,ib1,ib2)
	open(1,file = para_dir(ib1:ib2)//'subbasin.dat', status='old')
c	open(1,file = para_dir(ib1:ib2)//'kfs.dat', status='old')
c	read(1,*) atmp, nlevel
c	read(1,*)
c	do level = 1, nlevel-1
c	  if(level .eq.1) then  ! Level 1
c	    read(1,*)
c	    read(1,*) atmp, (kfs1(l1), l1=1,9)
c	  elseif(level .eq. 2) then ! Level 2
c	    read(1,*)
c	    do l1=1,9
c	      if(kfs1(l1).eq.1) read(1,*) atmp, (kfs2(l1,l2), l2=1,9)
cc	修改于2006.4.14
c	      if(kfs1(l1).ne.1) read(1,*) 
c	    end do
c	  elseif(level .eq. 3) then ! Level 3
c		read(1,*)
c	    do l1=1,9
c	      if(kfs1(l1).eq.1) then
c                do l2=1,9
c                  if(kfs2(l1,l2) .eq. 1)  read(1,*) atmp, 
c     :                 (kfs3(l1,l2,l3), l3=1,9)
c                end do
c	      endif
c	    end do
c	  endif
c	end do
c	close(1)
c      
         do i = 1,nc
        read(1,*)nsub,psubbasin(i),(pbasinup(i,j),j=1,3)

	  enddo
	   close(1)
c	nsub=0
c	do l1 = 9, 1, -1        ! Level 1
c	  if(kfs1(l1) .eq. 0) then
c	    nsub=nsub+1
c	    subbasin(nsub)='ws'//char(48+l1)//
c     :           char(48+0)//char(48+0)
c           print *,nsub,'   ',subbasin(nsub)
c	  else
c	    do l2 = 9, 1, -1    ! Level 2
c              if(kfs2(l1,l2) .eq. 0) then
c	        nsub=nsub+1
c	        subbasin(nsub)='ws'//char(48+l1)//char(48+l2)//
c     :               char(48+0)
c              print *,nsub,'   ',subbasin(nsub)
c	      else
c                do l3 = 9, 1, -1 ! Level 3
c                  nsub=nsub+1
c	            subbasin(nsub)='ws'//char(48+l1)//
c     :                 char(48+l2)//char(48+l3)
c                 print *,nsub,'   ',subbasin(nsub)
c               
c                end do
c              endif
c	    end do
c	  endif
c	end do

           nbasinup=0
        do i =1,nsub
	      write(ch4,'(i4.4)')psubbasin(i)
          subbasin(i)='ws'//ch4
	 	  enddo
          do i = 1,nsub
               do k =1,4
                  if(pbasinup(i,k) >0) then
                      do j =1,nsub
                         if(psubbasin(j)== pbasinup(i,k)) then
                            nbasinup(i,k)=j                    
	                      
	                   endif
	                enddo
	            endif
	         enddo
	    enddo


	print *,'number of sub-catchment:',nsub  ! added by gaobing 2013
      print * ,'ok1, finish of subbasin.dat', para_dir
c	pause 

c
c
c#######################################################################
c
c     Read land vegetation parameters
c
c#######################################################################
c
 	para_dir   ="../parameter/"

      open(1, file = para_dir(ib1:ib2)//'vege_para.dat',
     :        status='old')
	read(1,*)
	read(1,*)
	i = 1
15    read(1,*,end=25) iland, Kcanopy(i), LAImax(i), kcrop(i), 
     :                  root(i),anik(i), Sstmax(i), surfn(i), landuse
c      print *, iland, Kcanopy(i), LAImax(i), kcrop(i), 
c     +                  root(i),anik(i), Sstmax(i), surfn(i), landuse 
c      Sstmax(i)=Sstmax(i)
	i = i + 1
	goto 15
25	close(1)
	nland = i-1
	print *,'landuse types:',nland
      print *,'ok2, finish of vegetation.dat', para_dir
c	pause 

c
c
c#######################################################################
c
c     Read soil parameters
c
c#######################################################################
c
c	open(1, file = para_dir(ib1:ib2) // 'soil_water_para.dat',
c     :       status = 'old')
c	read(1,*)
c	read(1,*)
c	i = 1
c35    read(1,*, end=45) isoil,wsat(i), wrsd(i), alpha(i), watern(i),  
c     :                  ksat1(i), ksat2(i), kg(i), GWcs(i)

C      print *, i,isoil, wsat(i), wrsd(i), alpha(i), watern(i),  
C    +                  ksat1(i), ksat2(i), kg(i), GWcs(i)
C     pause 

	
c
c#######################################################################
c
c	Read base maps
c
c#######################################################################
c
  	call strlen(dem_dir,ib1,ib2)

	do i = 1, 6
	  if(i.eq.1) infile = grid_area
	  if(i.eq.2) infile = top_ele
	  if(i.eq.3) infile = top_len
	  if(i.eq.4) infile = top_slp
	  if(i.eq.5) infile = soil_map
	  if(i.eq.6) infile = geology_map
  	  call strlen(infile,ia1,ia2)
	  open(11,file= dem_dir(ib1:ib2)//infile(ia1:ia2), status='old')
        
c	  print * ,'open  dem_dir=, infile',  dem_dir, infile
c       pause 
	  
	  read(11,'(A, i16)') ncols,nnc
c	  write(*,'(A, i16)') ncols,nnc
	 	 
	  read(11,'(A, i16)') nrows,nnr
	  read(11,'(A, f16.3)') xllcorner,x0
	  read(11,'(A, f16.3)') yllcorner,y0
	  read(11,'(A, f16.0)') cellsize,gridsize
	  read(11,'(A, f8.0)') nodata,znodata

	  do ir=1,nnr
	    if(i.eq.1) read(11,*) (area(ir,ic),   ic=1,nnc)
	    if(i.eq.2) read(11,*) (ele(ir,ic),    ic=1,nnc)
	    if(i.eq.3) read(11,*) (length(ir,ic), ic=1,nnc)
	    if(i.eq.4) read(11,*) (slp(ir,ic),    ic=1,nnc)
		if(i.eq.5) read(11,*) (soil(ir,ic),   ic=1,nnc)
		if(i.eq.6) read(11,*) (Ds(ir,ic),     ic=1,nnc)
	  end do
	  close(11)
	end do
c
	do ir=1,nnr
	  do ic=1,nnc
         if (soil(ir,ic).eq.0) then
	    print *, 'ir,ic,soil(ir,ic)=', ir,ic,soil(ir,ic)
         endif
        if (Ds(ir,ic).eq.1) then
	    Ds(ir,ic)=Ds(ir,ic)*2.5          !added by mqh
         endif
         if(Ds(ir,ic).lt.1.and.Ds(ir,ic).ne.-9999) then
	     print *,'wrong of Ds(ir,ic) ', Ds(ir,ic)
	     pause
	    endif

	     
      if(area(ir,ic).eq.-9999) then
       ele(ir,ic)=-9999
	 length(ir,ic)=-9999
       slp(ir,ic)=-9999
	 soil(ir,ic)=-9999
       Ds(ir,ic)=-9999
      endif

	if(soil(ir,ic).eq.-9999 .and. Ds(ir,ic).ne.-9999) 
     : print *, 'soil1', soil(ir,ic), Ds(ir,ic), ir, ic
	if(soil(ir,ic).ne.-9999 .and. Ds(ir,ic).eq.-9999) 
     : print *, 'soil2', soil(ir,ic), Ds(ir,ic), ir, ic

cc		area(ir,ic)=area(ir,ic)*1000.0*1000.0

          if(area(ir,ic).ne.-9999) then
	       area(ir,ic)=area(ir,ic)*1000.0*1000.0
          endif
	 
	    if(Ds(ir,ic) .ne. -9999) then
c       changed by gaobing
c           Ds(ir,ic)=Ds(ir,ic)
c            if(Ds(ir,ic) .lt.2)  Ds(ir,ic)=Ds(ir,ic)+1.0 !
              if(Ds(ir,ic) .le. 0.0) print *, ir,ic,Ds(ir,ic)
		endif
c end of changed
	    if(slp(ir,ic) .ne. -9999) then
c		  slp(ir,ic)=tan(slp(ir,ic)*3.1415926/180.0)
		  slp(ir,ic) = slp(ir,ic)/100.0
		endif
c	    if (slp(ir,ic).lt.1.0e-6) slp(ir,ic)= 1.0e-6

	  end do
	end do
	
	
	

	do ir=1,nnr
	  do ic=1,nnc
		Dg(ir,ic)=-9999.0
	    if(area(ir,ic).ne.-9999.and.Ds(ir,ic).ne.-9999) then
		  Dg(ir,ic)=10.0*Ds(ir,ic)
		  if(ele(ir,ic) .ge. 3500.0) Dg(ir,ic)=amin1(2.0, Ds(ir,ic))
	      if (Dg(ir,ic) .lt. 0)  print *, "wrong in calculate Dg" 
c	       print *,	area(ir,ic),Ds(ir,ic),Dg(ir,ic),slp(ir,ic)
		endif
	  end do
	end do
c
      print * ,'ok4, finish of cell_area.asc'
      print * ,'ok5, finish of elevation.asc'
      print * ,'ok6, finish of slope_length.asc'
      print * ,'ok7, finish of slope.asc'
      print * ,'ok8, finish of soil_unit.asc'
      print * ,'ok9, finish of soil_depth.asc'
c	pause 

  	call strlen(para_dir,ib1,ib2)
	open(3,file=para_dir(ib1:ib2)//'soil_code.txt', status='old')
	read(3,*)
	i = 1
200	read(3,*,end=250) itmp, soiltyp(i), tmp
	i = i+1
	goto 200
250	close(3)
	nsoil = i-1
	print *,'soil types:',nsoil
c	pause
c
	do ir=1,nnr
	  do ic=1,nnc
	    if(area(ir,ic).ne.-9999.and.soil(ir,ic) .ne. -9999) then
	      k=0
		  do isoil = 1, nsoil
		    if(soil(ir,ic) .eq. soiltyp(isoil)) itmp=isoil !code --> record
	        if(soil(ir,ic) .eq. soiltyp(isoil)) k=k+1
		  end do

	      soil(ir,ic) = itmp

		 if(k .eq. 0) print *,'no this type of soil:',
     :		 ir,ic,soil(ir,ic)
	     if(k.gt.1) print *,'more one code for this type of soil:',
     :	    k, soil(ir,ic)
		endif
	  end do
	end do

      print * ,'ok10, finish of soil_code.txt'
c	pause

c############### ################################
         
c#      read soil parameter

c################################################
      call strlen(soil_dir,ib1,ib2)
c     changed by mqh
	do i = 1, 6
	  if(i.eq.1) infile = 'thsgm1k.asc'
	  if(i.eq.2) infile = 'thr1k.asc'
	  if(i.eq.3) infile = 'k01k.asc'
	  if(i.eq.4) infile = 'k11k.asc'
c	  if(i.eq.5) infile = 'ksg_heihe.asc'
	  if(i.eq.5) infile = 'alpha1k.asc'
        if(i.eq.6) infile = 'n1k.asc'
  	  call strlen(infile,ia1,ia2)
	open(11,file=soil_dir(ib1:ib2)//infile(ia1:ia2),status='old')
       

	  
	  read(11,'(A, i16)') ncols,nnc
 
	  read(11,'(A, i16)') nrows,nnr
	  read(11,'(A, f16.3)') xllcorner,x0
	  read(11,'(A, f16.3)') yllcorner,y0
	  read(11,'(A, f16.0)') cellsize,gridsize
	  read(11,'(A, f8.0)') nodata,znodata

	  do ir=1,nnr
	    if(i.eq.1) read(11,*) (wsat(ir,ic),     ic=1,nnc)
	    if(i.eq.2) read(11,*) (wrsd(ir,ic),     ic=1,nnc)
	    if(i.eq.3) read(11,*) (ksat1(ir,ic),    ic=1,nnc)
	    if(i.eq.4) read(11,*) (ksat2(ir,ic),    ic=1,nnc)
c		if(i.eq.5) read(11,*) (kg(ir,ic),       ic=1,nnc)
		if(i.eq.5) read(11,*) (alpha(ir,ic),    ic=1,nnc)
		if(i.eq.6) read(11,*) (watern(ir,ic),   ic=1,nnc)
	  end do
	  close(11)
	end do ! added by gaobing 20131125


       do ir = 1,nnr
	     do ic = 1,nnc

	     if(area(ir,ic) > 0 .and. ksat1(ir,ic) < 0) then
	 print *, 'error in ksat1' ,ir,ic
	     else
c      ksat1(ir,ic) = ksat1(ir,ic)*0.001/3600.0         !mm/hr--->m/s
      ksat1(ir,ic) = ksat1(ir,ic)*0.001/3600.0/2.4         !cm/day--->mm/hr--->m/s, mqh
	     endif
c	  if(area(ir,ic) > 0 .and. ksat2(ir,ic) < 0) then
c	 print *, 'error in ksat2' ,ir,ic
c	     else
c	ksat2(ir,ic) = ksat2(ir,ic)*0.001/3600.0         !mm/hr--->m/s 
c	     endif
c        if(area(ir,ic) > 0 .and. kg(ir,ic) < 0) then
c	 print *, 'error in kg' ,ir,ic
c	     else
c	kg(ir,ic)    = kg(ir,ic)*0.001/3600.0		     !mm/hr--->m/s
c	  endif

	 if(area(ir,ic) > 0 .and. wsat(ir,ic) < 0) then
	 print *, 'error in wsat' ,ir,ic
	 endif

	if(area(ir,ic) > 0 .and. wrsd(ir,ic) < 0) then
	 print *, 'error in wrsd' ,ir,ic
	 endif

	if(area(ir,ic) > 0 .and. alpha(ir,ic) < 0) then
	 print *, 'error in alpha' ,ir,ic
	 endif
	if(area(ir,ic) > 0 .and. watern(ir,ic) < 0) then
	 print *, 'error in watern' ,ir,ic
	 endif
	     enddo
	 enddo

  
         do ir = 1,nnr
	      do ic = 1,nnc
	     if(area(ir,ic) > 0)  then
	          GWcs(ir,ic)=0.15       
!	          ksat1(ir,ic) = 0.9*ksat1(ir,ic)    
	       ksat2(ir,ic) = 0.15*ksat1(ir,ic)  !from 0.1 to 0.15 ,mqh
!	          ksat2(ir,ic) = 0.09*ksat1(ir,ic)  !from 0.1 to 0.15 ,mqh
	           kg(ir,ic)    = 0.03*ksat1(ir,ic)

	      else
	          GWcs(ir,ic)=-9999           
	          ksat2(ir,ic) = -9999 
	           kg(ir,ic)    = -9999 
	      end if
	      enddo
	   enddo
c	print *,"ksat=",ksat1(i),ksat2(i),kg(i)

	suction = -1.02  ! suction, meter water-head   !*1.2,added by mqh
       suction2 = -1.02 
        do ir = 1,nnr
           do ic = 1,nnc
	     if(area(ir,ic) > 0 )  then
c	if(wsat(ir,ic).lt. 0.48 ) wsat(ir,ic)=0.48   ! test by gaobing
    
      call MoistureFromSuction_V(tmp, wsat(ir,ic), wrsd(ir,ic), 
     :                           watern(ir,ic), alpha(ir,ic), suction2)
      wfld(ir,ic) = tmp
	if(wfld(ir,ic).gt.wsat(ir,ic) .or. wfld(ir,ic).lt.wrsd(ir,ic))
     :        print *,"wfld out range wrsd ~ wsat"
             ksat2(ir,ic)= ksat2(ir,ic)*1.3
	       kg(ir,ic)=kg(ir,ic)*2.0
           endif
        enddo
	  enddo	 
c      if(wfld(ir,ic).lt. 0.60*wsat(ir,ic))wfld(ir,ic) = 0.60*wsat(ir,ic)
c	if(wfld(ir,ic).gt. 0.85*wsat(ir,ic))wfld(ir,ic) = 0.85*wsat(ir,ic)
c	print *, i, isoil, wsat(i), wfld(i)
c      i = i + 1
c	goto 35
c45	close(1)
      print * ,'ok3, finish of soil_water_para', para_dir
c	pause 


c      print *, i-1



c
c
c#######################################################################
c
c	Read land use map
c
c#######################################################################
c
      print *, 'nv=, nland=', nv, nland
  	call strlen(dem_dir,ib1,ib2)
	  infile = 'lu.asc'
	  iland = 1 
	  call strlen(infile,ia1,ia2)
	  open(1,file=dem_dir(ib1:ib2)//infile(ia1:ia2),
     $	status='old')
	  read(1, *) ncols,nnc
	  read(1, *) nrows,nnr
	  read(1, *) xllcorner,x0
	  read(1, *) yllcorner,y0
	  read(1, *) cellsize,gridsize
	  read(1, *) nodata,znodata
	  do ir=1,nnr
c	      print *,ic,ir
	    read(1,*) (land_use(ir,ic,iland), ic=1,nnc)

       
	  end do
	  close(1)

	  print *, 'ncols, nrows=', nnc, nnr
	  print *, 'xmin, ymin, cellsize=', x0, y0, gridsize
  
       print * ,'ok11, finish of read land use map'
c	pause 

c
c
c#######################################################################
c
c	Read land use ratio  distribution
c
c#######################################################################
c
      print *, 'nland=', nland
  	call strlen(para_dir,ib1,ib2)
	do i=1,nland
	  write(ch2, '(i2.2)') i
	  infile = 'lu_class'//ch2//'.asc'
c	  print *,nland,infile
	  call strlen(infile,ia1,ia2)
	  open(1,file=para_dir(ib1:ib2)//infile(ia1:ia2),
     $	status='old')
	  read(1, *) ncols,nnc
	  read(1, *) nrows,nnr
	  read(1, *) xllcorner,x0
	  read(1, *) yllcorner,y0
	  read(1, *) cellsize,gridsize
	  read(1, *) nodata,znodata
	  do ir=1,nnr
	    read(1,*) (land_ratio(ir,ic,i), ic=1,nnc)
	  end do
	  close(1)
	end do
	
      do ir=1,nnr
	  do ic=1,nnc
	    tmp=0.0
	    k=0
	    do iland=1,nland
	      if(land_ratio(ir,ic,iland) .gt.0) then
		    tmp=tmp+land_ratio(ir,ic,iland)
	        k=k+1
	      end if
	    end do
	    if(k.gt.0 .and. abs(tmp-1).gt.0.1E-3) print *,
     :                   "wrong in land_ratio",ir,ic,tmp 
	  end do
	end do
  
c      print * ,'ok11, finish of Lu_class01-10'
c	pause 

c
c

c#######################################################################
c
c	Read subbasin parameters
c
c#######################################################################c
	print *,'Reading catchment parameters ...'
      do isub = 1, nsub
        infile = subbasin(isub) // '_river'
c	  print *,isub,subbasin(isub),infile
c	  pause
        call strlen(infile, ia1, ia2)
        open(3,file=para_dir(ib1:ib2)//infile(ia1:ia2),status='old')
        read(3,*) nflow(isub)
c	  print *, isub, nflow(isub)
c	  pause
        do iflow = 1, nflow(isub)

          read(3,*) ngrid(isub,iflow),
     :              dx(isub,iflow),s0(isub,iflow),b(isub,iflow),
     :              roughness(isub,iflow),Dr(isub,iflow)
      roughness(isub,iflow)=roughness(isub,iflow)*2.0
cc		Dr(isub,iflow) = 1.5 * Dr(isub,iflow)					! 书签6
c	    print *, ngrid(isub,iflow),
c     :              dx(isub,iflow),s0(isub,iflow),b(isub,iflow),
c     :              roughness(isub,iflow),Dr(isub,iflow)
		if( dx(isub,iflow).lt.0.1) dx(isub,iflow) = 3000.0
               dx(isub,iflow) = dx(isub,iflow)*1.50
       if(s0(isub,iflow).eq.0) then   ! added by xujijun 20051104
	     s0(isub,iflow)=0.00001
		 print *, 'wrong in s0',s0(isub,iflow)  
       endif  
	 if(s0(isub,iflow).eq.-9999.0) then ! added by xujijun 20051104
	     s0(isub,iflow)=0.00001 
 		 print *, 'wrong in s0',s0(isub,iflow)  
	 endif


c	    b(isub,iflow)=2.0*b(isub,iflow)
c	    Dr(isub,iflow)=1.5*Dr(isub,iflow)

	    if(ngrid(isub,iflow) .gt. np) then
          print *,'more than np grids:',ngrid(isub,iflow),isub,iflow,np
          endif

	    read(3,*) (grid_row(isub,iflow,j),grid_col(isub,iflow,j),
     :               j = 1, ngrid(isub,iflow))
          if(s0(isub,iflow).le.0.1E-5) s0(isub,iflow) = 0.1E-5
        end do
        close(3)
	end do
      print * ,'ok12, finish of ws000_river'
c	pause 

c 
c
c#######################################################################
c
c	             read initial snow depth
c
c#######################################################################
c	  print *,'reading snow depth ..'
c	  call strlen(dem_dir,ic1,ic2)
c       write(ch4, '(i4.4)') 1993
c        open(3,file= dem_dir(ic1:ic2)//'snow_depth'//'.asc',
c     :         status='old')
c	  read(3, '(A, i16)') ncols,nnc
c	  read(3, '(A, i16)') nrows,nnr
c	  read(3, '(A, f16.3)') xllcorner,x0
c	  read(3, '(A, f16.3)') yllcorner,y0
c	  read(3, '(A, f16.0)') cellsize,gridsize
c	  read(3, '(A, f16.0)') nodata,znodata
c	  do ir=1,nnr
c	    read(3,*) (snow(ir,ic,1), ic=1,nnc)
c	    do ic=1,nnc
c	      if(snow(ir,ic,1).lt.0.0) snow(ir,ic,1)=0.0
c		  do i = 2, nland
c			snow(ir,ic,i) = snow(ir,ic,1)
c		  end do
c	    end do
c	  end do
c	  close(3)
c      print * ,'ok13, finish of snow'
c	pause

c  
c 
c
c#######################################################################
c
c	   read initial meteorological input interpolation
c
c#######################################################################
c

        call Initial_input

c
c
c#######################################################################
c
c     Model initialization
c
c#######################################################################
c
      print * ,'start of simulation'

cc     changed by xujijun 
         


      start  = 1
      finish = 0
      do isub = start_sub, end_sub
	  call hillslope_model(isub)
c        print * ,'fish of hillslope_model, isub=', isub 
c	  pause
	  call river_routing(isub) 
      end do

c      print * ,'fish of river_routing'
c	pause

c
c#######################################################################
c
c                         start simulation
c
c#######################################################################
c
      do itmp=1, 1        
	do 888 hydroyear = startyear, endyear
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc		
c		  !added by mqh  输出子流域面积
c      do isub=1,nsub
c        if(isub.eq.start_sub) basinarea2 = 0
c        subarea2(isub) = 0 
c        do iflow = 1, nflow(isub)
c		do ig = 1, ngrid(isub,iflow)
c		  ir = grid_row(isub,iflow,ig)
c		  ic = grid_col(isub,iflow,ig)
c		  if(area(ir,ic) .eq. -9999.0) print*,'wrong in grid-area'
c		  basinarea2 = basinarea2 + area(ir,ic)
c		  subarea2(isub) = subarea2(isub) + area(ir,ic)
c		end do
c	  end do

c      do ii1 =1,4
c         if(nbasinup(isub,ii1) > 0) then
c      subarea2(isub) = subarea2(isub)
c     $+subarea2(nbasinup(isub,ii1))    
c                endif
c        enddo
c	enddo  !子流域循环
c         call strlen(result1_dir,ib1,ib2)         
c       open(123,file = result1_dir(ib1:ib2)//'basinarea.txt', 
c     $status='unknown')
c      do isub=1,nsub
c      write(123,*) isub,subarea2(isub) 
c       enddo
c        close(123)
        !end of add
ccccccccccccccccccccccccccccccccccccccccccc       
        start  = 0
	  finish = 0
        year = hydroyear
        write(ch4, '(i4.4)') year
c        print *, 'year=', year
c       pause 
c      
c#######################################################################
c	
c	             read monthly LAI
c
c#######################################################################
	  print*,'reading NDVI ...'

	  call strlen(ndvi_dir,ic1,ic2)
	  year_tmp = year
	  if(year .le. 1982)  year_tmp = 1982
         if(year .ge. 2010)  year_tmp = 2010
	   if(year .gt. 1998.and.year .lt.2008) then
	  write(ch4, '(i4.4)') year_tmp
c	  print *, 'start read ndvi'
c	  print *, 'dem_dir', dem_dir, year
c	  pause 
	 	  do i = startmonth, endmonth
        if(i.lt.10) then 
	  open(3,file= ndvi_dir(ic1:ic2)//'ndvi'//ch4//'_'//
     &      char(48+i)//'_0.asc', status='old')
	  else 
   	  open(3,file= ndvi_dir(ic1:ic2)//'ndvi'//ch4//'_1'//
     &	  char(48+i-10)//'_0.asc',  status='old')
	  endif
	   
	    read(3, *) ncols,nnc
	    read(3, *) nrows,nnr
	    read(3, *) xllcorner,x0
	    read(3, *) yllcorner,y0
	    read(3, *) cellsize,gridsize
	    read(3, *) nodata,znodata
c          print *, 'x0=, y0=', x0, y0

	    do ir=1,nnr
	      read(3,*) (yndvi(ir,ic,i,1), ic=1,nnc)
	    end do
	  end do

	  do ir=1,nnr
	    do ic=1,nnc
             do i=startmonth, endmonth  ! added by xujijun
		      if(area(ir,ic).eq.-9999) then 
	             yndvi(ir,ic,i,1) = -9999
	          endif
	       enddo

		end do
	  end do
	  close(3)

	  do i = startmonth, endmonth
        if(i.lt.10) then 
	  open(3,file= ndvi_dir(ic1:ic2)//'ndvi'//ch4//'_'//
     &      char(48+i)//'_1.asc', status='old')
	  else
          open(3,file= ndvi_dir(ic1:ic2)//'ndvi'//ch4//'_1'//
     &	  char(48+i-10)//'_1.asc',  status='old')
	  endif
	   
	    read(3, *) ncols,nnc
	    read(3, *) nrows,nnr
	    read(3, *) xllcorner,x0
	    read(3, *) yllcorner,y0
	    read(3, *) cellsize,gridsize
	    read(3, *) nodata,znodata
c          print *, 'x0=, y0=', x0, y0

	    do ir=1,nnr
	      read(3,*) (yndvi(ir,ic,i,2), ic=1,nnc)
	    end do
	  end do

	  do ir=1,nnr
	    do ic=1,nnc
             do i=startmonth, endmonth  ! added by xujijun
		      if(area(ir,ic).eq.-9999) then 
	             yndvi(ir,ic,i,2) = -9999
	          endif
	       enddo

		end do
	  end do
	  close(3)
	   
	  
	  
	  do i = startmonth, endmonth
        if(i.lt.10) then 
	  open(3,file= ndvi_dir(ic1:ic2)//'ndvi'//ch4//'_'//
     &      char(48+i)//'_2.asc', status='old')
	  else
          open(3,file= ndvi_dir(ic1:ic2)//'ndvi'//ch4//'_1'//
     &	  char(48+i-10)//'_2.asc',  status='old')
	  endif
	   
	    read(3, *) ncols,nnc
	    read(3, *) nrows,nnr
	    read(3, *) xllcorner,x0
	    read(3, *) yllcorner,y0
	    read(3, *) cellsize,gridsize
	    read(3, *) nodata,znodata
c          print *, 'x0=, y0=', x0, y0

	    do ir=1,nnr
	      read(3,*) (yndvi(ir,ic,i,3), ic=1,nnc)
	    end do
	  end do

	  do ir=1,nnr
	    do ic=1,nnc
             do i=startmonth, endmonth  ! added by xujijun
		      if(area(ir,ic).eq.-9999) then 
	             yndvi(ir,ic,i,3) = -9999
	        endif
	       enddo
		  end do
	  end do
	  close(3)
	  
	  else   !标记 ，前面说的>1998的，1982~1998的是剩下的else
	 write(ch4, '(i4.4)') year_tmp
         	  do i = startmonth, endmonth
        if(i.lt.10) then 
	  open(3,file= ndvi_dir(ic1:ic2)//'ndvi'//ch4//
     &   '_'//char(48+i)//'_1.asc', status='old')
	  else 
   	  open(3,file= ndvi_dir(ic1:ic2)//'ndvi'//ch4//'_1'//
     &	  char(48+i-10)//'_1.asc',  status='old')
	  endif
	   
	    read(3, *) ncols,nnc
	    read(3, *) nrows,nnr
	    read(3, *) xllcorner,x0
	    read(3, *) yllcorner,y0
	    read(3, *) cellsize,gridsize
	    read(3, *) nodata,znodata
c          print *, 'x0=, y0=', x0, y0

	    do ir=1,nnr
	      read(3,*) (yndvi(ir,ic,i,1), ic=1,nnc)
	    end do
	  end do

	  do ir=1,nnr
	    do ic=1,nnc
             do i=startmonth, endmonth  ! added by xujijun
		      if(area(ir,ic).eq.-9999) then 
	             yndvi(ir,ic,i,1) = -9999
	          endif
	       enddo

		end do
	  end do
	  close(3)

        	  do i = startmonth, endmonth
        if(i.lt.10) then 
	  open(3,file= ndvi_dir(ic1:ic2)//'ndvi'//ch4//
     &      '_'//char(48+i)//'_2.asc', status='old')
	  else 
   	  open(3,file= ndvi_dir(ic1:ic2)//'ndvi'//ch4//'_1'//
     &	  char(48+i-10)//'_2.asc',  status='old')
	  endif
	   
	    read(3, *) ncols,nnc
	    read(3, *) nrows,nnr
	    read(3, *) xllcorner,x0
	    read(3, *) yllcorner,y0
	    read(3, *) cellsize,gridsize
	    read(3, *) nodata,znodata
c          print *, 'x0=, y0=', x0, y0

	    do ir=1,nnr
	      read(3,*) (yndvi(ir,ic,i,2), ic=1,nnc)
	    end do
	  end do

	  do ir=1,nnr
	    do ic=1,nnc
             do i=startmonth, endmonth  ! added by xujijun
		      if(area(ir,ic).eq.-9999) then 
	             yndvi(ir,ic,i,2) = -9999
	          endif
	       enddo

		end do
	  end do
	  close(3)


	  endif

	  print *, 'finish of read ndvi'
	  print *, 'year=' , year 
c	  pause 


      ! 标记
      if(iflai.eq.1)  then   !added by mqh,我没有用到lai
        print *, 'reading lai'

	call strlen(lai_dir,ic1,ic2)
	  year_tmp = year
	  if(year .le. 2000)  year_tmp = 2000
         if(year .ge. 2012)  year_tmp = 2012
	  
	  write(ch4, '(i4.4)') year_tmp
c	  print *, 'start read ndvi'
c	  print *, 'dem_dir', dem_dir, year
c	  pause 
	 	  do i = startmonth, endmonth
        if(i.lt.10) then 
	  open(3,file= lai_dir(ic1:ic2)//'lai_'//ch4//'0'//
     &      char(48+i)//'1.asc', status='old')
	  else 
   	  open(3,file=lai_dir(ic1:ic2)//'lai_'//ch4//'1'//
     &	  char(48+i-10)//'1.asc',  status='old')
	  endif
	   
	    read(3, *) ncols,nnc
	    read(3, *) nrows,nnr
	    read(3, *) xllcorner,x0
	    read(3, *) yllcorner,y0
	    read(3, *) cellsize,gridsize
	    read(3, *) nodata,znodata
c          print *, 'x0=, y0=', x0, y0

	    do ir=1,nnr
	      read(3,*) (lai(ir,ic,i,1), ic=1,nnc)
	    end do
	     close(3)
	  end do

	  do ir=1,nnr
	    do ic=1,nnc
             do i=startmonth, endmonth  ! added by xujijun
		      if(area(ir,ic) > 0 .and.  lai(ir,ic,i,1) < 0) then 
	             lai(ir,ic,i,1) =0
	          endif
	          
	       enddo

		end do
	  end do
	  

	  do i = startmonth, endmonth
        if(i.lt.10) then 
	  open(3,file= lai_dir(ic1:ic2)//'lai_'//ch4//'0'//
     &      char(48+i)//'2.asc', status='old')
	  else
          open(3,file= lai_dir(ic1:ic2)//'lai_'//ch4//'1'//
     &	  char(48+i-10)//'2.asc',  status='old')
	  endif
	   
	    read(3, *) ncols,nnc
	    read(3, *) nrows,nnr
	    read(3, *) xllcorner,x0
	    read(3, *) yllcorner,y0
	    read(3, *) cellsize,gridsize
	    read(3, *) nodata,znodata
c          print *, 'x0=, y0=', x0, y0

	    do ir=1,nnr
	      read(3,*) (lai(ir,ic,i,2), ic=1,nnc)
	    end do
	         close(3)
	  end do

	  do ir=1,nnr
	    do ic=1,nnc
             do i=startmonth, endmonth  ! added by xujijun
		      if(area(ir,ic) > 0 .and.  lai(ir,ic,i,2) < 0) then 
	             lai(ir,ic,i,2) = 0
	          endif
	       enddo

		end do
	  end do
	  
	   
	  
	  
	  do i = startmonth, endmonth
        if(i.lt.10) then 
	  open(3,file= lai_dir(ic1:ic2)//'lai_'//ch4//'0'//
     &      char(48+i)//'3.asc', status='old')
	  else
          open(3,file= lai_dir(ic1:ic2)//'lai_'//ch4//'1'//
     &	  char(48+i-10)//'3.asc',  status='old')
	  endif
	   
	    read(3, *) ncols,nnc
	    read(3, *) nrows,nnr
	    read(3, *) xllcorner,x0
	    read(3, *) yllcorner,y0
	    read(3, *) cellsize,gridsize
	    read(3, *) nodata,znodata
c          print *, 'x0=, y0=', x0, y0

	    do ir=1,nnr
	      read(3,*) (lai(ir,ic,i,3), ic=1,nnc)
	    end do
	          close(3)
	  end do

	  do ir=1,nnr
	    do ic=1,nnc
             do i=startmonth, endmonth  ! added by xujijun
		      if(area(ir,ic) > 0 .and. lai(ir,ic,i,3) < 0) then 
	             lai(ir,ic,i,3) = 0
	        endif
	       enddo
		  end do
	  end do


	  do i = startmonth, endmonth
        if(i.lt.10) then 
	  open(3,file= lai_dir(ic1:ic2)//'lai_'//ch4//'0'//
     &      char(48+i)//'4.asc', status='old')
	  else
          open(3,file= lai_dir(ic1:ic2)//'lai_'//ch4//'1'//
     &	  char(48+i-10)//'4.asc',  status='old')
	  endif
	   
	    read(3, *) ncols,nnc
	    read(3, *) nrows,nnr
	    read(3, *) xllcorner,x0
	    read(3, *) yllcorner,y0
	    read(3, *) cellsize,gridsize
	    read(3, *) nodata,znodata
c          print *, 'x0=, y0=', x0, y0

	    do ir=1,nnr
	      read(3,*) (lai(ir,ic,i,4), ic=1,nnc)
	    end do
	          close(3)
	  end do

	  do ir=1,nnr
	    do ic=1,nnc
             do i=startmonth, endmonth  ! added by xujijun
		      if(area(ir,ic) > 0 .and. lai(ir,ic,i,4) < 0) then 
	             lai(ir,ic,i,4) = 0
	        endif
	       enddo
		  end do
	  end do
      end if  !added by mqh,因为不用lai，用if略过这段


	  

c
c
c#######################################################################
c
c	    read meteorological data
c
c#######################################################################

cc	2222	删除于2006.4.14	 ! skip input meteorological data

c
c	   print *,"After read data",basinarea
c        print *,"hehe"
        idc=0
	  ihc=0
	  do 777 im = 1, 12    !标记，书签
	    month = im
c#######################################################################
c
c	    read meteorological data
c
c#######################################################################
c          call Initial_input

c          print *, 'start to read meteorological data'
c		print *, 'year=,im=', year, im
c		pause 
          
		call T_P_E_cal(year,im)
c       fill the missing data using nearest available value 

c	  print *, 'finish of T_P_E'
c	  pause 
        
c        print *,'running at1: ',year,month
c	  go to 777
	    
        do idd=1,dayinmonth(im)
	    do ir=1,nrow
	      do ic=1,ncol
	        if(soil(ir,ic).ne.-9999 ) then

c	      fill minimum temperature
     		      if (tmin_daily(ir,ic,idd).le.-99.0) then
	          k=0
	          tmp_tmin=0.0
	          tmp_num=0.0
411             k=k+1
                do ir2=ir-k,ir+k
			    do ic2=ic-k,ic+k
	              if(tmin_daily(ir2,ic2,idd).ne.-9999.0) then
				     tmp_tmin=tmp_tmin+tmin_daily(ir2,ic2,idd)
	                 tmp_num=tmp_num+1
				  end if 
				end do
			  end do
			  if(tmp_num.eq.0) goto 411 
			  tmin_daily(ir,ic,idd)= tmp_tmin/tmp_num
	          end if

c	        fill maximum temperature
     		      if (tmax_daily(ir,ic,idd).le.-99.0) then
	          k=0
	          tmp_tmax=0.0
	          tmp_num=0.0
421             k=k+1
                do ir2=ir-k,ir+k
			    do ic2=ic-k,ic+k
	              if(tmin_daily(ir2,ic2,idd).ne.-9999.0) then
				     tmp_tmax=tmp_tmax+tmax_daily(ir2,ic2,idd)
	                 tmp_num=tmp_num+1
				  end if 
				end do
			  end do
			  if(tmp_num.eq.0) goto 421 
			  tmax_daily(ir,ic,idd)= tmp_tmax/tmp_num
	          end if

c		 fill rainfall
     		    if (rain_daily(ir,ic,idd).lt.0.0) then
	          k=0
	          tmp_rain=0.0
	          tmp_num=0.0
401             k=k+1
                do ir2=ir-k,ir+k
			    do ic2=ic-k,ic+k
	              if(rain_daily(ir2,ic2,idd).ge.0.0) then
				     tmp_rain=tmp_rain+rain_daily(ir2,ic2,idd)
	                 tmp_num=tmp_num+1
				  end if 
				end do
			  end do
			  if(tmp_num.eq.0) goto 401 
			  rain_daily(ir,ic,idd)= tmp_rain/tmp_num
	          end if
c
c	adjust the snowfall according to the elevation
c       tmp1 = 0.5*(tmax_daily(ir,ic,idd)+tmin_daily(ir,ic,idd))
c       tmp2 = 1.0 + ele(ir,ic) / 7500.0
c       if(tmp1.le.-0.5 .and. ele(ir,ic).gt.2000.0)
c     :                rain_daily(ir,ic,idd) = tmp2*rain_daily(ir,ic,idd)
c
c			fill Ep
     		    do i=1,nland  
			  if (evap_daily(ir,ic,idd,i).lt.0.0) then
	          k=0
	          tmp_evap=0.0
	          tmp_num=0.0
431             k=k+1
                do ir2=ir-k,ir+k
			    do ic2=ic-k,ic+k
	              if(evap_daily(ir2,ic2,idd,i).ge.0.0) then
				     tmp_evap=tmp_evap+evap_daily(ir2,ic2,idd,i)
	                 tmp_num=tmp_num+1
				  end if 
				end do
			  end do
			  if(tmp_num.eq.0) goto 431 
			  evap_daily(ir,ic,idd,i)= tmp_evap/tmp_num
	          end if
              end do
			end if       
	      end do   
		end do
	  end do
c
c	    month = im + 3
c	    if(month .gt. 12) month = month - 12
c		if(iy.gt.1 .and. month.eq.1) year = year + 1
          print *,'running at: ',year, month
c         day loop
	    if(mod(year,4).eq.0 .and. month.eq.2) dayinmonth(month)=29
	    if(mod(year,4).ne.0 .and. month.eq.2) dayinmonth(month)=28
	    do 666 day = 1, dayinmonth(month)
            idc = idc + 1
c           print *,'running at: ', year, month,day

c           hour loop
            do 555 hour = 1, 24
              ihc = ihc + 1
c	        print *,'running at: ', year, month,day,hour
           
c#######################################################################
c	
c	        sub-catchment loop	  
c
c#######################################################################
c              isub = 0
c              do l1 = 9, 1, -1                                 ! Level 1
c                if(kfs1(l1) .eq. 0) then
c                  isub = isub+1
c                  subbasin(isub)='ws'//char(48+l1)//
c     :                            char(48+0)//char(48+0)
c	            if (isub.ge.start_sub .and. isub.le. end_sub) then
c				   call hillslope_model(isub)
c                     call river_routing(isub, 1, l1, l2, l3) 
c	            endif
c                else
c                  do l2 = 9, 1, -1                             ! Level 2
c                    if(kfs2(l1,l2) .eq. 0) then
c                      isub = isub+1
c                      subbasin(isub)='ws'//char(48+l1)//char(48+l2)//
c     :                               char(48+0)
c   	                if (isub.ge.start_sub .and. isub.le. end_sub) then
c  				      print *, 'ok1'
c                          print *, isub,subbasin(isub)
c					  call hillslope_model(isub)
c 				      print *, 'ok2'

                    do isub = 1,nsub
                      call hillslope_model(isub)
				      call river_routing(isub)
				    enddo

c  				      print *, 'ok3'

c	                endif
c                    else
c                      do l3 = 9, 1, -1                        ! Level 3
c                        isub = isub+1
c                        subbasin(isub)='ws'//char(48+l1)//
c     :                                  char(48+l2)//char(48+l3)
c	                if (isub.ge.start_sub .and. isub.le. end_sub) then
c					  call hillslope_model(isub)
c					  call river_routing(isub,3,l1,l2,l3)
c	                 endif
c                        p3(l3)=isub
c                      end do
c                    endif
c                    p2(l2)=isub
c                  end do
c                endif
c                p1(l1)=isub
c              end do
555         continue   
666	    continue
777	  continue
888   continue
      end do
999   stop
	end
      
c
c     
      subroutine strlen(str, l1,l2)
      character str*200
      integer i,l1,l2,k
      k=0
      do i = 1, 200
        if(k.eq.0 .and. str(i:i).NE.' ') then
          l1=i
          k=1
        elseif(k.eq.1 .and. str(i:i).EQ.' ') then
          l2 = i-1
          return
        endif
      end do
      l2 = i
      return
      end
 
c
c
c**********************************************************************
c                                                                     *
c                 ***** River Routing Model *****                     *
c                     (Kinematic Wave method)                         *
c                                                                     *
c**********************************************************************
c
	subroutine river_routing(isub)
c
c#######################################################################
c
c   PURPOSE:
c    Solving the kinematic wave model using nonlinear numerical scheme 
c
c#######################################################################
c      
c     AUTHOR: Dr. Dawen YANG
c     08/1998
c
c     MODIFICATION HISTORY: 07/2002
c
c#######################################################################
c
c
c**********************************************************************
c  Variable definition
c
      implicit none
      include 'dims2.inc'
c#######################################################################
c
c     Sub-basin parameters
c
c#######################################################################
c
c      character*100 basin             ! name of simulated basin
	character*6   subbasin(nc)      ! name of sub-basin changed by gaobing 2013
	integer       psubbasin(nc)                   ! added by gaobing  2013
	integer       nbasinup(nc,4)
	integer       nsub              ! total number of sub-basins
	integer       isub              ! sub-basin number
c
c#######################################################################
c
c     Topographic parameters
c
c#######################################################################
c
	integer nflow(nc)    ! total number of flow-intervals in a sub-basin
	integer iflow        ! flow-interval number
	real    dx(nc,nx)    ! flow interval length (m)
      real    Dr(nc,nx)    ! river depth in the flow interval (m)
      real    Drw(nc,nx)   ! water depth in the river of the flow interval (m)
c
	real    s0(nc,nx)        ! slope of river bed (ND)
	real    b(nc,nx)         ! width of river (m)
	real    roughness(nc,nx) ! Manning's roughness
c
c#######################################################################
c
c     Date and time variables
c
c#######################################################################
c
      real    dt             ! time step (hour)
	integer year           ! calendar year   (19XX)
	integer month          ! month in a year (1-12)
	integer day            ! day in a month  (1-30)
	integer hour           ! hour in a day   (1-24)
	integer dayinmonth(12) ! days of a month
	integer hydroyear      ! year for hydro simulation
	integer startmonth     ! start month in the hydro-year
	integer endmonth	   ! end   month in the hydro-year
	integer startday       ! start day in the first month of the hydro-year
	integer endday         ! end   day in the last  month of the hydro-year
	integer iy, im, id ,ih ! for year, month, day and hour loop
      integer idc            ! contineous day in a year (1-366)
	integer ihc            ! contineous hour in a year (1-366*24)

	integer start, finish  ! 0 - faulse,    1 - true
	integer start_sub, end_sub  ! for partial simulation of a basin
c
c#######################################################################
c
c     Discharge variables
c
c#######################################################################
c
      real qin(nx)         ! lateral inflow of one flow interval (m3/sec/m)
	                     !, from 'hydro.inc'
c
c      real beta            ! Manning's equation parameter
      real q1(nc,nx)       ! discharge of last time step (m^3/s)
      real q2(nc,nx)       ! discharge of current time step (m^3/s)
      real qr1(nc,nx)     ! discharge of reservoir flowout last time step(m^3/s)
      real qr2(nc,nx)  ! discharge of reservoir flowout current time step(m^3/s)
      real qlin1(nc,nx)    ! lateral inflow of last time step (m^3/m/s)
      real qlin2(nc,nx)    ! lateral inflow of current time step (m^3/m/s)
	real Qh(nc,8800)     ! hourly mean discharge (m3/s)
      common /qh/ qh   !mqh
      real Qd(nc,366)	     ! daily average discharge (m^3/s)
	real Qm(nc,12)       ! monthly mean discharge (m3/s)
      real Qobs(nc,366)    ! added by gaobing
cc	real re_capacity(nc,nx) !capacity of reservoir in location(isub,iflow), m3
cc	real re_capacity0(nc,nx)  ! initial capacity of reservoir 

      character*200 gauge_name    ! name of discharge gauges
      character*200 gauge_name2    ! name of discharge gauges 
c
c      save q1       ! discharge of last time step (m^3/s)
      save q2       ! discharge of current time step (m^3/s)
      save qr1      ! discharge of reservoir flowout last time step(m^3/s)
c      save qr2      ! discharge of reservoir flowout current time step(m^3/s)
      save qlin1    ! lateral inflow of last time step (m^3/m/s)
      save qlin2    ! lateral inflow of current time step (m^3/m/s)
!	save Qh       ! hourly mean discharge (m3/s)
      save Qd	      ! daily average discharge (m^3/s)
	save Qm       ! monthly mean discharge (m3/s)
	save Qobs     ! added by gaobing 
c      
c#######################################################################
c
c      Pfafstetter numbering system
c
c#######################################################################
c
c	integer kfs1(10)         ! key of forward subdivision of level 1
c	integer kfs2(10,10)      ! key of forward subdivision of level 2
c	integer kfs3(10,10,10)   ! key of forward subdivision of level 3
c	integer nlevel           ! total level of subdivision in Pfafstetter scheme
	integer l1, l2, l3
      integer level
	integer p1(9), p2(9), p3(9)
c
c#######################################################################
c
c     Other variables
c
c#######################################################################
c
      character*200  para_dir    ! directory for parameters
      character*200  dem_dir    ! directory for dem
	character*200  data_dir    ! directory for input data
	character*200  result1_dir ! directory for storing simulation result
	character*200  result2_dir ! directory for storing simulation result
	character*200  simul_dir   ! directory for storing temporal result

	character*4    ch4         ! a 4-character variable
      character*3    ch3
      real    tmp, value         ! temporary variables
	integer i, j, k, m ,l        !, n      ! temporary variables
c	integer jj, is, il         ! temporary varaibles
	integer ia1, ia2, ib1, ib2, ic1, ic2,id1,id2
	real    criterion, h1, h2, f, df
	real    q1_upper, q2_upper,rs    !, i_sub,  lo_dis
	integer ii1, ii2           ! , ii3
	real Qtotal,basinarea
	real	b_tmp              ! adjust riverbed-width during heavy-flood
	real	rs_tmp             ! adjust rs during heavy-flood
							   ! added in 2006.4.30
c
c#######################################################################
c
c     Common blocks.
c
c#######################################################################
c
	common /date1  / year, month, day, hour, idc,ihc
     	common /date2  / hydroyear,startmonth,startday,endmonth,endday
	common /date3  / dayinmonth
	
	common /simulation/  start, finish, dt, start_sub, end_sub
	common /river1    /  nsub, subbasin 
	common /river2    /  nflow, dx, Dr
	common /river3    /  s0, b, roughness
      common /river4    /  psubbasin,nbasinup ! added by gaobing 2013
	common /location  /  para_dir,data_dir,result1_dir,dem_dir,
     :                	result2_dir,simul_dir
	common /TEST   /basinarea

	common /pfaf/p1,p2,p3          !Pfaftetter basin number

      common /lateralflow/  qin
	common /river_status/ Drw, qr2  ! depth of river water 

cc	common /reservoir/ re_capacity ,re_capacity0  
      integer im2,idd2,hour2,jj
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c	data beta/0.6/
c
c
c      print *, isub,level,l1,l2,l3
c	print *, 'start of river_routing'
	 
	if(month.eq.startmonth .and. 
     :            	day.eq.startday .and. hour.eq.1) then
	  do i = 1, nsub
	    do j= 1, 366
	      Qd(i,j)=0.0
	    end do
	  end do
c        added by gaobing
        do i = 1, nsub
           do j= 1, 12
               Qm(i,j)=0.0
	     end do
	  end do 
	 do i = 1, nsub
           do j= 1, 366
               Qobs(i,j)=0.0
	     end do
	  end do
c        end of added by gaobing
	endif
c     
c
c#######################################################################
c
c     Initialize the model parameters.
c
c#######################################################################
c
      if(start .eq. 1) then
	  do i = 1, nsub
	    do j = 1, nflow(isub)
	      qlin1(i,j) = 0.0
	    end do
	  end do
c
c#######################################################################
c
c     Specify the initial condition
c
c#######################################################################
 
        if(inicon.eq.1) then
          call strlen(simul_dir, ia1,ia2)

          open(3,file = simul_dir(ia1:ia2)//subbasin(isub)//'I_flow2',
     $                  status='old')
          read (3,*)(i,q1(isub,j), j=1, nflow(isub))
          close(3)
        else
          q1(isub,1)=0.2    !mqh,changed from 0.5 to 0.2
          do iflow = 2, nflow(isub)
            q1(isub,iflow) = q1(isub,iflow-1)+0.1  !mqh,changed from 0.4 to 0.1
          end do
        endif
c       calculate initial river water depth
	  do iflow = 1, nflow(isub)
	    criterion = 0.01
	    k = 1
	    h1  = q1(isub,iflow) / b(isub,iflow)
5	    tmp = roughness(isub,iflow) * q1(isub,iflow)/
     :          sqrt(s0(isub,iflow))
	    f   = b(isub,iflow)*h1 - 
     :          (tmp**0.6)*((b(isub,iflow)+2.0*h1)**0.4)
	    if(k.gt.1 .and. abs(f) .lt. criterion) goto 8
	    df = b(isub,iflow)-0.8*((tmp/(b(isub,iflow)+2.0*h1))**0.6)
	    h2 = h1-f/df
	    h1 = h2
	    if(k.ge.10) goto 8
	    k = k+1
	    goto 5
 8	    Drw(isub,iflow) = h2		! print *,"initial h=",h2
                         
	  end do
c         initial lateral inflow
c	  do iflow = 1, nflow(isub)
c	    qlin1(isub,iflow) = 0
c	  end do
c         initial qr1 inflow
	  do iflow = 1, nflow(isub)
	    qr1(isub,iflow) = q1(isub,iflow)
	  end do
c	  return
	endif 
c      
c
c#######################################################################
c
c     Initialize the reservoir parameters.
c
c#######################################################################
c





c
	if(start .eq. 1) return
c
c#######################################################################
c
c     lateral inflow
c
c#######################################################################
c
c	print *, 'ok1 river_routing'

	if(month.eq.startmonth.and.day.eq.startday.and.ihc.eq.1) then
        do iflow = 1, nflow(isub)
           qlin1(isub,iflow)= qin(iflow)/dx(isub,iflow)
        end do
 	end if

      do iflow = 1, nflow(isub)
	  qlin2(isub,iflow) = qin(iflow)/dx(isub,iflow)
c        if(month.eq.startmonth .and. day.eq.startday .and. ihc.eq.1)
c     :      print *,qlin1(isub,iflow),qlin2(isub,iflow),isub,iflow
      end do

c	print *, 'ok2 river_routing'

c     
c#######################################################################
c
c     River routing
c
c#######################################################################
c
c      m=int(ddt/dt)
      m = 1
      do 555 j = 1, m                  !time loop
        do 333 iflow = 1, nflow(isub)  !river segment loop
c
c#######################################################################
c
c	junction boundary condition: define the river network
c
c#######################################################################
c
	  q1_upper=0.
	  q2_upper=0.
	  if(iflow .eq. 1) then  ! changed by gaobing 2013
	     if(psubbasin(isub) > 1000 .and. psubbasin(isub) < 2000) then
              q1_upper=0.
	        q2_upper=0.
c	        print *, isub
	     else
              do ii1 =1,4
                if(nbasinup(isub,ii1) > 0) then
	          ii2=nbasinup(isub,ii1)
                q1_upper=q1_upper+qr1(ii2,nflow(ii2))
	          q2_upper=q2_upper+qr2(ii2,nflow(ii2))
c	       print *, isub, ii2,qr1(ii2,nflow(ii2)),qr2(ii2,nflow(ii2))
c	           pause
c	          print *, isub, ii2
	          endif
	        enddo
	     endif
c	    if(level .eq. 1) then
c	      if(l1.eq.9 .or. mod(l1,2).eq.0) then
c	        q1_upper=0.
c	        q2_upper=0.
c	      else
c	        ii1=p1(l1+1)
c	        ii2=p1(l1+2)
c	        q1_upper=qr1(ii1,nflow(ii1))+qr1(ii2,nflow(ii2))
c	        q2_upper=qr2(ii1,nflow(ii1))+qr2(ii2,nflow(ii2))
c	      endif
c	    elseif(level .eq. 2) then
c	      if(l2.eq.9 .or. mod(l2,2).eq.0) then
c	        q1_upper=0.
c	        q2_upper=0.
c	        if(l2.eq.9 .and. (l1.ne.9 .and. mod(l1,2).ne.0)) then 
c     			!taken from upper level
c		      ii1=p1(l1+1)
c		      ii2=p1(l1+2)
c	          q1_upper=qr1(ii1,nflow(ii1))+qr1(ii2,nflow(ii2))
c	          q2_upper=qr2(ii1,nflow(ii1))+qr2(ii2,nflow(ii2))
c	        endif
c	      else
c	        ii1=p2(l2+1)
c	        ii2=p2(l2+2)
c	        q1_upper=qr1(ii1,nflow(ii1))+qr1(ii2,nflow(ii2))
c	        q2_upper=qr2(ii1,nflow(ii1))+qr2(ii2,nflow(ii2))
c	      endif
c	    elseif(level .eq. 3) then
c	      if(l3 .eq. 9 .or. mod(l3,2) .eq. 0) then
c	        q1_upper=0.
c	        q2_upper=0.
c	        if(l3.eq.9 .and. (l2.ne.9 .and. mod(l2,2).ne.0)) then 
c			  !taken from upper level
c		      ii1=p2(l2+1)
c		      ii2=p2(l2+2)
c	          q1_upper=qr1(ii1,nflow(ii1))+qr1(ii2,nflow(ii2))
c	          q2_upper=qr2(ii1,nflow(ii1))+qr2(ii2,nflow(ii2))
c	        endif
c	        if(l3.eq.9 .and. l2.eq.9 .and. 
c    :                    (l1.ne.9 .and. mod(l1,2).ne.0)) then	 
c			  !taken from upper level
c		      ii1=p1(l1+1)
c		      ii2=p1(l1+2)
c	          q1_upper=qr1(ii1,nflow(ii1))+qr1(ii2,nflow(ii2))
c	          q2_upper=qr2(ii1,nflow(ii1))+qr2(ii2,nflow(ii2))
c	        endif
c	        if(l3.eq.9 .and. l2.eq.9 .and. l1.eq.9) then 
c			  !free upper boundary
c	          q1_upper=0.
c	          q2_upper=0.
c	        endif
c	      else
c	        ii1=p3(l3+1)
c	        ii2=p3(l3+2)
c	        q1_upper=qr1(ii1,nflow(ii1))+qr1(ii2,nflow(ii2))
c	        q2_upper=qr2(ii1,nflow(ii1))+qr2(ii2,nflow(ii2))
c	      endif
c	    endif
	  else
	    q1_upper=qr1(isub,iflow-1)
	    q2_upper=qr2(isub,iflow-1)
	  endif	  	               !end of river network definition
c	
c#######################################################################
c
c      routing to next point
c
c#######################################################################
c
c	print *, 'ok21 river_routing'

          rs=roughness(isub,iflow)
c	    print *,"in nkws", dt,dx(isub,iflow),b(isub,iflow),
c     :    s0(isub,iflow),rs,qlin1(isub,iflow),qlin2(isub,iflow),q1_upper
c     :    ,q1(isub,iflow),q2_upper,q2(isub,iflow),Drw(isub,iflow),iflow
   
         if(s0(isub,iflow).eq.0) then   ! added by xujijun 20051104
	     s0(isub,iflow)=0.00001
		 print *, 'wrong in s0',s0(isub,iflow)  
c	     pause
          endif  
	 if(s0(isub,iflow).eq.-9999.0) then ! added by xujijun 20051104
	     s0(isub,iflow)=0.00001 
 		 print *, 'wrong in s0',s0(isub,iflow)  
c         pause
	 endif

	b_tmp = b(isub,iflow)
	rs_tmp = rs
	if(Drw(isub,iflow) .gt. 0.5*Dr(isub,iflow)) then		! 书签2
		b_tmp=1.1*b(isub,iflow)
		rs_tmp = 1.5*rs
	endif
	if(Drw(isub,iflow) .gt. 1.0*Dr(isub,iflow)) then
		b_tmp=1.5*b(isub,iflow)
		rs_tmp = 3.0*rs
	endif
          call nkws(dt,dx(isub,iflow),b_tmp
     :		,s0(isub,iflow),rs_tmp,	qlin1(isub,iflow),
     :        qlin2(isub,iflow),q1_upper,q1(isub,iflow),
     :        q2_upper,q2(isub,iflow),Drw(isub,iflow))

c	print *, 'ok22 river_routing'

          qr2(isub,iflow)=q2(isub,iflow)
		
cc		call reservoir_op(isub,iflow,qr2(isub,iflow)    
cc     :         ,re_capacity(isub,iflow),0)  
 
c 	print *, 'ok23 river_routing'
    
          !1-with; 0-withou reservoir operation model
	if(isub.eq.1 .and. iflow.eq.1 .and. q2_upper.gt.0.0) then     
	    print *,"out nkws", dt,dx(isub,iflow),b(isub,iflow),
     :    s0(isub,iflow),rs,qlin1(isub,iflow),qlin2(isub,iflow),q1_upper
     :   ,q1(isub,iflow),q2_upper,q2(isub,iflow),Drw(isub,iflow),idc,ihc
	    pause
	endif
 333      continue   
          do iflow = 1, nflow(isub)
            qlin1(isub,iflow) = qlin2(isub,iflow)
            q1(isub,iflow)    = q2(isub,iflow)
            qr1(isub,iflow)   = qr2(isub,iflow)
          end do
c
 555  continue

c	print *, 'ok3 river_routing'


c#######################################################################
c
c	discharge output
c
c#######################################################################
c
c      print *,"discharge output"
      Qh(isub,ihc)  = qr2(isub,nflow(isub))
	Qd(isub,idc)  = Qd(isub,idc) + qr2(isub,nflow(isub))/24.0
	!标记2   qr2
      Qm(isub,month)= Qm(isub,month) + qr2(isub,nflow(isub))/
     :                (float(dayinmonth(month))*24.0)
c   
c added by gaobing

c        if(subbasin(isub) .eq. 'ws1001') then
c	  Qobs(isub,idc)=Qobs(isub,idc)+qr2(isub,nflow(isub))/24.0
c	  else if(subbasin(isub) .eq. 'ws535') then
c       Qobs(isub,idc)=Qobs(isub,idc)+qr2(isub,nflow(isub))/24.0
c	  elseif(subbasin(isub) .eq. 'ws760')  then
c       Qobs(isub,idc)=Qobs(isub,idc)+qr2(isub,nflow(isub))/24.0
c	  elseif(subbasin(isub) .eq. 'ws277')  then
c        Qobs(isub,idc)=Qobs(isub,idc)+qr2(isub,nflow(isub))/24.0
c	  elseif(subbasin(isub) .eq. 'ws271')  then
c        Qobs(isub,idc)=Qobs(isub,idc)+qr2(isub,nflow(isub))/24.0
c	  elseif(subbasin(isub) .eq. 'ws610')  then
c        Qobs(isub,idc)=Qobs(isub,idc)+qr2(isub,nflow(isub))/24.0
c        else
        Qobs(isub,idc)=Qobs(isub,idc)+qr2(isub,nflow(isub))/24.0
c	  endif 
c end of added
c     Save the status at end of the year
	Qtotal=0.0 

c      if(day .eq. dayinmonth(month) .and. hour .eq. 24) then
	if(month.eq.endmonth .and. day.eq.endday .and. hour.eq.24) then
        
c	  if(subbasin(isub) .eq. 'ws100') then
	  if(subbasin(isub) .eq. 'ws4017') then ! changed by xujijun   !mqh
	    do i=1,idc
		  Qtotal=Qtotal+Qd(isub,i)*3600*24*1000.0/basinarea
c		  print *,basinarea
c		  pause
		  !标记
		end do  
c	    print *,"Qtotal=",Qtotal
	  end if

	  call strlen(simul_dir, ia1,ia2)
	  open(9,file = simul_dir(ia1:ia2)//subbasin(isub)//'I_flow2',
     $	  status='replace')
c	   print *,	simul_dir(ia1:ia2)//subbasin(isub)//'I_flow'
         do iflow = 1, nflow(isub)
            write(9,110) iflow, q2(isub,iflow)
 110        format(2x,i4,f20.3)
         end do
         close(9)

c	print *, 'ok4 river_routing'


c	save river discharge at selected gauges
	  gauge_name = 'no_gauge'
	  if(subbasin(isub) .eq. 'ws4017') gauge_name = 'maduwang'  !mqh
!	  if(subbasin(isub) .eq. 'ws3025') gauge_name = 'luolicun'
c	  if(subbasin(isub) .eq. 'ws3038') gauge_name = 'qilian'
c	  if(subbasin(isub) .eq. 'ws277') gauge_name = 'dage'
c	  if(subbasin(isub) .eq. 'ws271') gauge_name = 'daiying'
c	  if(subbasin(isub) .eq. 'ws610') gauge_name = 'sandaoying'
c	  if(subbasin(isub) .eq. 'ws650') gauge_name = 'pangduo'
c	  if(subbasin(isub) .eq. 'ws590') gauge_name = 'yangcun'
c	  if(subbasin(isub) .eq. 'ws610') gauge_name = 'baijiachuan'
c	  if(subbasin(isub) .eq. 'ws510') gauge_name = 'longmen'
c	  if(subbasin(isub) .eq. 'ws241') gauge_name = 'zhangjiashan'
c	  if(subbasin(isub) .eq. 'ws250') gauge_name = 'xianyang'
c	  if(subbasin(isub) .eq. 'ws140') gauge_name = 'baimashi'
c	  if(subbasin(isub) .eq. 'ws110') gauge_name = 'yichang'
c


	  call strlen(result1_dir,ib1,ib2)
	  if(year.eq.startyear ) then
	    open(9,file = result1_dir(ib1:ib2)//
     :             subbasin(isub)//'.daily', status='unknown')
	  else
	    open(9,file = result1_dir(ib1:ib2)//
     :    subbasin(isub)//'.daily', access='append', status='old')
	  endif
	  j = 0
        do im = startmonth, endmonth
		do id = 1, dayinmonth(im)
		  j = j+1
            write(9, '(3i6,f16.3)') year,im,id, Qd(isub,j)  !标记1，Qd(isub,j)
          end do
	  end do
        close(9)
     
c     
        
c 
c
        
	  call strlen(gauge_name,ic1,ic2)
	  if(gauge_name(ic1:ic2) .eq. 'no_gauge') goto 999
c
	  call strlen(data_dir,ia1,ia2)
	  call strlen(result1_dir,ib1,ib2)
c        write(ch4, '(i4)') hydroyear
	  if(year.eq.startyear ) then
          open(8,file = result1_dir(ib1:ib2)//
     :		          gauge_name(ic1:ic2)//'.daily', status='unknown')
        else
          open(8,file = result1_dir(ib1:ib2)//          
     :	gauge_name(ic1:ic2)//'.daily', access='append', status='old')
	  endif

	  open(3,file=data_dir(ia1:ia2)//
     $		gauge_name(ic1:ic2)//'.txt', status='old')
	  k=0
 250    read(3,*,end=350) iy,im, id,value
	  if(k.eq.0 .and. iy.gt.hydroyear) goto 350
	  if(iy.eq.hydroyear .and. im.eq.startmonth 
     :	                     .and. id .eq.startday) then
	    k=1
	    j=0	   
	  endif
	  if(k.eq.1) then
	    j=j+1
		if(gauge_name(ic1:ic2) .eq. 'maduwang') then
		write(8,*) iy,im,id, Qd(isub,j)*160./163., value    !mqh
		!163
		  !adjust the area difference 调整流域面积差别
		else
	      write(8,*) iy,im, id,Qobs(isub,j)
		endif
		if(iy. eq. year .and. 
     $       im .eq. endmonth .and. id.eq.endday) goto 100
	  endif
	  goto 250
 350	  if(k.eq.0) write(8,450)
 450	  format(1x,'No discharge data recorded for this year')
 100	  close(8)
	  close(3)  
	endif
c	print *, 'ok5 river_routing'

      	 !added by mqh,to export hourly discharge
	  gauge_name2 = 'no_gauge'
	  if(subbasin(isub) .eq. 'ws4017') gauge_name2 = 'maduwang'
	  if(subbasin(isub) .eq. 'ws3025') gauge_name2 = 'luolicun'


	      call strlen(gauge_name2,id1,id2)
	      if(gauge_name2(id1:id2) .eq. 'no_gauge') goto 999
	      call strlen(result1_dir,ib1,ib2)
          
	if(year.eq.startyear.and.ihc.eq.1 ) then
              open(9,file = result1_dir(ib1:ib2)//
     :		          gauge_name2(id1:id2)//'.hourly', status='unknown')
            else 
               open(9,file = result1_dir(ib1:ib2)//
     :   gauge_name2(id1:id2)//'.hourly', access='append', status='old')
           endif
        if(gauge_name2(id1:id2) .eq. 'maduwang') then
                  write(9,*) year,month,day,hour,Qh(isub,ihc)*160./163.
        elseif  (gauge_name2(id1:id2) .eq. 'luolicun') then
                     write(9,*) year,month,day,hour,Qh(isub,ihc)   
        end if
        close(9)
        !标记  
        !end of added
        
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
       !added by mqh输出子流域的小时流量  ，mqh、2014.5.3
c               goto 19922
! 999    if(year.ge.2000.and.(isub.eq.136.or.isub.eq.153)) then
999    if(year.ge.2000) then
       call strlen(result1_dir,ib1,ib2)
         write(ch3, '(i3.3)') isub 
         if(year.eq.2000.and.ihc.eq.1 ) then
	    open(19,file = result1_dir(ib1:ib2)//'q_'//
     :             ch3//'.hourly', status='unknown')
	  else
	    open(19,file = result1_dir(ib1:ib2)//'q_'//
     :    ch3//'.hourly', access='append', status='old')
         endif
       write(19, '(4i6,f16.6)') year,month,day,hour,Qh(isub,ihc) 

        close(19)
       !end of add
       endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


19922	return

      	 
	end     
c
c
c
c
c
c**********************************************************************
c                                                                     *
c              ***** Reservoir control Model *****                    *
c                                                                     *
c**********************************************************************
c
c	Definition of Variales:
c     Q2       - known, discharge of flow in reservoir in current time step
c	re_capacity  - known, reservoir capacity last time step
c	qout     - unknown, outflow of reservoir
c**********************************************************************

cc	删除于2006.4.16
c
c
c
c
c**********************************************************************
c                                                                     *
c                 ***** Kinematic Wave Model *****                    *
c             (Nonlinear Scheme using Newton's method)                *
c                                                                     *
c**********************************************************************
c
c	Definition of Variales:
c
c            (time)
c              ^ 
c              |
c	          Q1(known)       Q2(unknown)
c
c               Q01(known)      Q02(known)    --> (distance)
c
c	Q01, Q02 - known, discharge of last time step
c	Q1       - known, discharge of current time step
c     Q2       - unknown, discharge of current time step
c	qlin0    - known, lateral inflow of last time step
c	qlin     - known, lateral inflow of current time step
c
c	y  -- water depth
c	p  -- wetted perimeter
c	b  -- river width
c	s0 -- river slope
c
c**********************************************************************
c
	subroutine nkws(dt, dx, b, s0, roughness,
     :			    qlin0, qlin, Q01,Q02,Q1, Q2, y)
	real Q01, Q02, Q1, Q2
	real qlin0, qlin
	real dt
	real dx, b, s0, roughness
	real y, p, beta, criterion
c
c     print *, 'start of nkws'
c	print *, 'dt,dx,b,s0,roughness'
c	print *, dt,dx,b,s0,roughness
c	print *, 'qlin0,qlin,Q01,Q02,Q1,q2,y'
c	print *, qlin0,qlin,Q01,Q02,Q1,q2,y

	beta = 0.6
	criterion = 0.0001
	k   = 1
	h1  = Q02/b
 15	tmp = roughness*Q02/sqrt(s0)
	f   = b*h1 - (tmp**0.6)*((b+2.0*h1)**0.4)
	if(k.gt.1 .and. abs(f) .lt. criterion) goto 18
	df = b - 0.8*((tmp/(b+2.0*h1))**0.6)
	h2 = h1 - f/df
	h1 = h2
	if(k.ge.30) goto 18
	k  = k+1
	goto 15
 18	y = h2
	p = b+2.0*y  ! for rectangular channel
c
c     the initial discharge estimated using linear scheme
	alfa = (roughness*p**(2.0/3.0)/sqrt(s0))**0.6
	if((Q02+Q1).le.0.0) then
	  cc   =0
       else
	  cc   = (0.5*(Q02+Q1))**(beta-1.0)
      endif
	aa   = dt*Q1/dx + alfa*beta*Q02*cc + 0.5*(qlin0+qlin)*dt
	bb   = dt/dx + alfa*beta*cc
	qq1  = aa/bb

	if(qq1 .le. 0.1e-5)  qq1=0.1e-5 
c     
c     Using Newton's method to calculate discharge
	ctmp = dt*Q1/dx + alfa*Q02**beta + 0.5*dt*(qlin+qlin0)
	k= 1
 20	f= dt*qq1/dx + alfa*qq1**beta-ctmp
	if((k.gt.1).and.(abs(f).le.criterion)) goto 30
	df  = dt/dx + alfa*beta*qq1**(beta-1.0)
	qq2 = qq1 - f/df     
c	print *,k,ctmp,f,df,alfa ,qq1,qq2
      
	if(qq2 .le. 0.1e-5) then
	  qq2=0.0
	  goto 30
	endif
	
      qq1 = qq2
	if(k .ge. 30) goto 30
	k = k+1
	goto 20
c
 30	Q2 = qq2
c
	k   = 1
	h1  = Q2/b
 45	tmp = roughness*Q2/sqrt(s0)
	f   = b*h1-(tmp**0.6)*((b+2.0*h1)**0.4)
	if(k.gt.1 .and. abs(f) .lt. criterion) goto 50
	df = b-0.8*((tmp/(b+2.0*h1))**0.6)
	h2 = h1-f/df
	h1 = h2
	if(k.ge.30) goto 50
	k = k+1
	goto 45
 50	y = h2
c	print *,"nkws q=" ,q2,"  drw=",y
	return
	end
c       
c
c
c**********************************************************************
c                                                                     *
c            ***** Hydrological Response on Hillslope  *****          *
c                                                                     *
c**********************************************************************	
c
	subroutine hillslope_model(isub)
c
c#######################################################################
c
c     PURPOSE:
c             A physically-based hillslope hydrological model.   
c
c#######################################################################
c      
c     AUTHOR: Dawen Yang
c     08/1998 
c                                                                               
c     MODIFICATION HISTORY: 09/2000, 11/2001,7/2002,11/2005 11/2013 by gaobing
c
c#######################################################################
c
c
c#######################################################################
c
c     Variable Declarations.
c
c#######################################################################
c
      implicit none
      include 'dims2.inc'
c
c#######################################################################
c
	character*6   subbasin(nc)      ! name of sub-basin ! changed by gaobing 2013
	integer       nsub              ! total number of sub-basins
	integer       psubbasin(nc)
	integer       nbasinup(nc,4)     !  added by gaobing 2013
	integer       isub              ! sub-basin number
	integer       i_sub             ! sub-basin number 
c
c#######################################################################
c
c     Topographic parameters
c
c#######################################################################
c
	real area(nrow,ncol)     ! area of the local grid (m2)
c	real ele(nrow,ncol)      ! elevation of the local grid (m)
	real slp(nrow,ncol)      ! average slope of the local grid (ND)
	real length(nrow,ncol)   ! average hillslope length (m)
	real Ds(nrow,ncol)       ! average depth of the topsoil (m)
	real Dg(nrow,ncol)       ! average depth of the uncinfined acwuifer (m)

	integer nflow(nc)        ! total number of flow-intervals in a sub-basin
	integer iflow            ! flow-interval number
	real    dx(nc,nx)        ! flow interval length (m)
      real    Dr(nc,nx)        ! river depth in the flow interval (m)
	real    s0(nc,nx)        ! slope of river bed (ND)
	real    b(nc,nx)         ! width of river (m)
	real    roughness(nc,nx) ! Manning's roughness

	integer ngrid(nc,nx)       ! number of grids in this flow-interval
	integer grid_row(nc,nx,np) ! row of grids in this flow-interval
	integer grid_col(nc,nx,np) ! column of grids in this flow-interval
c
c#######################################################################
c
c     Atmospheric forcing data
c
c#######################################################################
c
	real rain_daily(nrow,ncol,31)  ! daily precipitaiton (mm)
	real pre_hour(nrow,ncol,31,24)
	real tmin_daily(nrow,ncol,31)  ! daily minimum air temperature (degree)
	real tmax_daily(nrow,ncol,31)  ! daily minimum air temperature (degree)
	real evap_daily(nrow,ncol,31,nv) !daily value of potential evaporation(mm)
	real snow(nrow,ncol,nv)        ! snow depth (mm water) of the initial year
	real smf                        ! snow melting factor
cc	real snow1(nrow,ncol)
c
c#######################################################################
c
c     Soil-landuse distribution
c
c#######################################################################
c
c     character*200 soilname			 ! name of soil
c	character*200 landuse            ! name of land use
      integer soiltyp(ns)              ! Soil type of the whole area
      integer landtyp(nv)              ! landuse type of the whole area
	integer soil(nrow, ncol)         ! soil code of each grid
      real land_ratio(nrow,ncol,nv)    ! Area fraction of each landuse type
	integer land_use(nrow*10,ncol*10,1)
	integer nsoil                    ! number of soil types in whole area
	integer isoil                    ! soil number
	integer nland                    ! number of landuse types in whole area
	integer iland                    ! land-use number

      real    NDVI(nrow,ncol,12,3)    ! monthly NDVI
	real    yndvi(nrow,ncol,12,3) ! spot NDVI
	real    LAI(nrow,ncol,12,4)     ! monthly LAI
	real    Kcanopy(nv)           ! vegetation coverage of a landuse type
	real    LAImax(nv)            ! maximum LAI in a year
	real    root(nv)              ! root depth (m)
	real    anik(nv)              ! soil anisotropy ratio 
      real    Sstmax(nv)            ! maximum surface water detension (mm)
	real    surfn(nv)             ! surface roughness (Manning's coefficient)

	real    kcrop(nv)             ! evaporation coefficient of crop
	real    kmoist                ! evaporation coefficient of soil mositure  

	real    wsat(nrow,ncol)	        ! saturated soil moisture
	real    wrsd(nrow,ncol)		    ! residual soil moisture 
	real    wfld(nrow,ncol)            ! soil moisture at field capacity
	real    alpha(nrow,ncol)	        ! soil water parameter
	real    watern(nrow,ncol)		    ! soil water parameter
	real    ksat1(nrow,ncol)	  ! saturated hydraulic conductivity at surface (mm/h)
      real    ksat2(nrow,ncol)     ! saturated hydraulic conductivity at bottom (mm/h)
c      real    suction             ! soil suction
	real    soil_con_f       ! soil water conservation factor
	real    detension        ! surface detension (mm)
      real    detension2  !MQH
      	real Qh(nc,8800)     ! hourly mean discharge (m3/s)
      common /qh/ qh   !mqh
c
c	save    LAI
c
c#######################################################################
c
c     Date and time variables
c
c#######################################################################
c
      real    dt             ! time step (second)
	integer year           ! calendar year   (19XX)
	integer month          ! month in a year (1-12)
	integer day            ! day in a month  (1-30)
	integer hour           ! hour in a day   (1-24)
	integer dayinmonth(12) ! days of a month
	integer hydroyear      ! year for hydro simulation
	integer startmonth     ! start month in the hydro-year
	integer endmonth	   ! end   month in the hydro-year
	integer startday       ! start day in the first month of the hydro-year
	integer endday         ! end   day in the last  month of the hydro-year
c	integer iy, im, id, ih ! for year, month, day and hour loop
	integer idc            ! contineous day in a year (1-366)
	integer ihc            ! contineous hour in a year (1-366*24)

	integer start, finish  ! 0 - faulse,    1 - true
	integer start_sub, end_sub  ! for partial simulation of a basin
c
c#######################################################################
c
c     UZ parameters
c
c#######################################################################
c
	integer	layer(nrow,ncol)			! number of UZ layer
	real	D(nrow,ncol,nlayer)			! depth of each UZ layer(m)
	real	k0(nrow,ncol,nlayer)	! saturated hydraulic conductivity (mm/hr)
	real	w(nrow,ncol,nv,nlayer)		! soil moisture

	save	layer			            ! number of UZ layer
	save	D			                ! depth of each UZ layer(m)
	save	k0	                        ! saturated hydraulic conductivity (mm/hr)
	save	w                    		! soil moisture
c
c#######################################################################
c
c     runoff  
c
c#######################################################################
c
c	real    exc_inf             ! infiltration excess
c	real	q_hr(np,nv)			! 
	real    q_hillslope         ! surface runoff of one simulation 
                                  ! unit (m3/s/m, flow into river)
	real    qin(nx)             ! lateral inflow into river
c
c#######################################################################
c
c     groundwater parameters
c
c#######################################################################
c
	real    kg(nrow,ncol)          ! hydraulic conductivity of groundwater (m/sec) 
	real    GWcs(nrow,ncol)		! groundwater storage coefficient changed by gaobing 20131125
c
c#######################################################################
c
c     status parameters
c
c#######################################################################
c
	real    Cstmax(nrow,ncol,nv)	    ! maximal canopy storage
	real	Cst(nrow,ncol,nv)			! canopy storage
	real	Sst(nrow,ncol,nv)			! surface storage
c	real	Gst(nrow,ncol)				! groundwaterstorage
	real	Dgl(nrow,ncol,nv)	   	! depth to groundwater level in the grid
c	integer	Isat(nrow,ncol)				! saturation status
	real    GWst(nrow,ncol,nv)	        ! 

	real	prec                ! rainfall (mm/hr)
	real	temper              ! temperature(degree,hourly)
	real    snowmelt            ! snowmelt equvalent water (mm)
      real    Pnet                ! net precipitation (mm)
	real	Ep		            ! potential evaporation(mm/hr)
	real    Etr                 ! reference transpiration (mm)
	real    Es                  ! reference evaporation (mm)

	real    Eact                ! actual evaporation
	real	EfromCanopy			! actual evaporation from canopy storage
	real	EfromCrop			! actual transpiration from crop
c	real	EfromRoot(nlayer)	! actual transpiration calculated for each root layer
	real	EfromSurface		! actual evaporation from soil surface
	real    DeficitCst          ! deficit of canopy storage (mm)

	save	Cst     			! canopy storage
	save	Sst     			! surface storage
	save	Dgl         	   	! depth to groundwater level in the grid
	save    GWst    	        ! 
c
c#######################################################################
c
c     irrigation parameters
c
c#######################################################################
c



c#######################################################################
c
c     irrigation P_M parameters
c
c#######################################################################
c
	integer m_grow(6,7)
	integer d_grow(6,7)  
c
c#######################################################################
c
c     water depth & discharge of river 
c
c#######################################################################
c
	real	Drw(nc,nx)          ! river water depth
	real	discharge(nc,nx)    ! river flow discharge

cc	real    re_capacity(nc,nx)  ! capacity of reservoir
cc	real    re_capacity0(nc,nx)  ! capacity of reservoir
cc	real    re_cap(nc,nx,366)   ! daily capacity of reservoir 
c
c#######################################################################
c
c     water balance
c
c#######################################################################
c
      real raind(nrow,ncol,366)       ! daily rainfall (mm)
	real epd(nrow,ncol,366)         ! daily potential evaporation (mm)
	real eactd(nrow,ncol,366)       ! daily actual evapotranspiration (mm)
	real ecy(nrow,ncol,366)
	real ecp(nrow,ncol,366)
	real ese(nrow,ncol,366)
	real runoffd(nrow,ncol,366)     ! daily runoff (mm)
      real srunoff(nrow,ncol,366)
	real groundoff(nrow,ncol,366)
	real soiloff(nrow,ncol,366)

      real srunoff2(nc,8800)   !mqh ，不同产流
	real groundoff2(nc,8800)
	real soiloff2(nc,8800)
      real soiloff_total
      real srunoff_total
      real groundoff_total

	real annual_rain        ! annual precipitation 
cc      real annual_irrigation  ! annual irrigation water use(mm) 
	real annual_runoff      ! annual runoff
	real annual_Eact        ! annual actual evaporation
	real annual_Ecy
	real annual_Ecp
	real annual_Ese
      real annual_Cst         ! annual canopy interception
	real annual_Ssts        ! annual surface storage(snow)
	real annual_Sstr        ! annual surface storage(rain)
	real annual_SBst        ! annual subsurface storage
	real annual_Gst         ! annual groundwater storage
	real annual_Cst0        ! initial value of annual canopy storage
	real annual_Ssts0       ! initial value of annual surface storage(snow)
	real annual_Sstr0       ! initial value of annual surface storage(rain)
	real annual_SBst0       ! initial value of annual subsurface storage
	real annual_Gst0        ! initial value of groundwater storage

cc	real annual_re          ! annual reservoir capacity
cc	real annual_re0         ! initial value of reservoir capacity

	real unbalance          ! simulated unbalance error
c
	save eactd          	! daily actual evapotranspiration (mm)
	save ecy
	save ecp
	save ese
	save raind
	save runoffd            ! daily runoff (mm)
	save groundoff          ! added by gaobing
	save srunoff
	save soiloff
    	save groundoff2          ! added by 
	save srunoff2
	save soiloff2
	save annual_rain        ! annual precipitation 
cc      save annual_irrigation  !
	save annual_runoff      ! annual runoff
	save annual_Eact        ! annual actual evaporation
	save annual_Ecy
	save annual_Ecp
	save annual_Ese
      save annual_Cst         ! annual canopy interception
cc	save annual_Sst         ! annual surface storage
	save annual_Ssts        ! annual surface storage(snow)
	save annual_Sstr        ! annual surface storage(rain)
	save annual_SBst        ! annual subsurface storage
	save annual_Gst         ! annual groundwater storage
	save annual_Cst0        ! initial value of annual canopy storage
c	save annual_Sst0        ! initial value of annual surface storage
	save annual_Ssts0       ! initial value of annual surface storage(snow)
	save annual_Sstr0       ! initial value of annual surface storage(rain)
	save annual_SBst0       ! initial value of annual subsurface storage
	save annual_Gst0        ! initial value of groundwater storage

cc	save annual_re          !
cc	save annual_re0         !
c
c#######################################################################
c
c     daily catchment-total rainfall and evaporation
c
c#######################################################################
c
	real	total_rain_d(366)		! basin mean precipitation
	real	total_Ep_d(366)			! basin mean potential evaporation
	real	total_Eact_d(366)		! basin mean actual evaporation
	real	total_runoff_d(366)		! basin mean runoff
c
	save	total_rain_d
	save	total_Ep_d
	save	total_Eact_d
	save	total_runoff_d
c
c#######################################################################
c
c     daily subcatchment-total hillslope runoff(for each subcatchment)
c
c#######################################################################
c
c	real	rain_day(nc,366)		! subbasin mean daily precipitation
c	real	Ep_day(nc,366)			! subbasin mean potential evaporation
c	real	srunoff(nc,366)			! subbasin mean surface runoff
c	real	groundoff(nc,366)		! subbasin mean groundwater runoff
c	real	Eact_day(nc,366)		! subbasin mean actual evaporation
c
c	save	rain_day
c	save	Ep_day
c	save	srunoff
c	save	groundoff
c	save	Eact_day
c
c#######################################################################
c
c     one-D distribution of hillslope runoff for each subcatchment(hourly)
c
c#######################################################################
c
c	real rain(nc,nx)		! precipitation
c	real runoff1(nc,nx)		! surface runoff
c	real runoff2(nc,nx)		! subsurface runoff
c
c	save rain		
c	save runoff1
c	save runoff2
c
c#######################################################################
c
c     spatial distribution monthly output
c
c#######################################################################
c
	real evap_month(nrow,ncol)      ! monthly mean evaporation
	real soil_month(nrow,ncol)
	real runoff_month(nrow,ncol)
	real grunoff_month(nrow,ncol)
	real drunoff_month(nrow,ncol)
	real srunoff_month(nrow,ncol)
	real ecy_month(nrow,ncol)
	real ecp_month(nrow,ncol)
	real ese_month(nrow,ncol)
     
	                                ! monthly mean soil moisture of a selected layer
c	real gwtable_month(nrow,ncol)	! monthly mean groundwater table
c
	save evap_month
	save ecy_month
	save ecp_month
	save ese_month
	save soil_month
	save runoff_month
	save grunoff_month
	save drunoff_month
	save srunoff_month
c	save gwtable_month
c
c#######################################################################
c
c     used in subroutine runoff, please refer this subroutine
c
c#######################################################################
c
      integer d1_nuz
	real    d1_land_ratio, d1_slope, d1_length, d1_Ds, d1_Dg,  
     :        d1_Dr, d1_anik, d1_wsat, d1_wrsd, d1_wfld, 
     :        d1_alpha, d1_watern, d1_kg, d1_GWcs, d1_GWst, d1_Dgl,
     :        d1_Drw, d1_Sst, d1_qsub,d1_dt, d1_snow
	real    d1_ssub  ! added by gaobing
      real    d1_deltz(nlayer), d1_k0(nlayer), d1_w(nlayer)
c	
c#######################################################################
c
c     Other variables
c
c#######################################################################
c
      character*200  para_dir    ! directory for parameters
	character*200  data_dir    ! directory for input data
	character*200  dem_dir    ! directory for input dem
	character*200  result1_dir ! directory for storing simulation result
	character*200  result2_dir ! directory for storing simulation result
	character*200  simul_dir   ! directory for storing temporal variables
c	character*200  infile

      character*5 ncols,nrows
      character*9 xllcorner,yllcorner
      character*8 cellsize
      character*12 nodata
	integer nnr, nnc
	real x0, y0
	real gridsize
	real znodata

	character*2    ch2        ! a 2-character variable
	character*3    ch3        ! a 3-character variable
	character*4    ch4        ! a 4-character variable

c      character*200 atmp
      real    tmp,temp          !,mean , value   ! temporary variables
	integer i, j              !, k,n,m temporary variables
	integer ig,nuz			  !,ii, jj
	integer ir, ic            !, iyy, imm, idd
	real    basinarea, subarea(nc)    ! total basin area
c	real    areacheck		          ! for checking basin area
	integer ia1, ia2  , ib1, ib2 !, ic1, ic2
	integer ii2,j1,j2,jj,imm,iidd,nday
	real    tmp1,tmp2, tmp3, D0,Ddt
c	integer itmp
	real    f
c	real    r1,T1,Ep1 
	real  para_r, para_a, y, tmpEp	!tmpinf_max,, tmpinf2,tmpinf1
c	integer i_root
	real  qin_tmp0, qin_tmp
	real water_depth, surface_n, waterhead, power, fice, area_frac
	real tmppnet
c
c
      real  rain_tmp(300,1000), Ep_tmp(300,1000), T_tmp(300,1000)

	real c1,c2,c3,c4,c5,dLAI
      integer kn
	save  subarea, fice, area_frac
      integer isub2  !added by mqh
      integer luolicun2(nc)  !存储哪些是罗李村的上游站   added by mqh

      data luolicun2 /1,1,0,0,0,1,1,1,0,0,1,1,1,1,1,1,0,1,0,
     :1,1,1,1,0,0,1,1,0,1,0,1,1,
     : 1,1,1,1,0,1,0,1,1,0,0,1,0,1,1,1,1,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,
     :0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,
     :0,1,0,0,0,0,1,1,1,1,1,1,1,1,0,1,1,0,0,1,0,0,1,0,0,1,0,1,1,1,1,1,1,
     :1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
c
c#######################################################################
c
c     NDVI-->LAI look-up table
c
c#######################################################################
c
      integer i_ndvi,ii_ndvi
	real T_NDVI(20),Shurb_LAI(20)
	real Forest_LAI(20),IRR_LAI(20),N_IRR_LAI(20)

	data T_NDVI/0.025,0.075,0.125,0.175,0.225,0.275,0.325,
     $             0.375,0.425,0.475,0.525,0.575,0.625,0.675,
     $             0.725,0.775,0.825,0.875,0.925,0.975/
      data Shurb_LAI/0,0,0.2663,0.3456,0.4357,0.5213,0.6057,
     $               0.6951,0.8028,0.9313,1.102,1.31,1.598,
     $               1.932,2.466,3.426,4.638,6.328,6.328,6.328/
      data Forest_LAI/0,0,0.1516,0.1973,0.2686,0.3732,0.5034,
     $                0.6475,0.7641,0.9166,1.091,1.305,1.683,
     $                2.636,3.557,4.761,5.52,6.091,6.091,6.091/
      data IRR_LAI/0,0,0.3199,0.431,0.5437,0.6574,0.7827,0.931,
     $             1.084,1.229,1.43,1.825,2.692,4.229,5.362,
     $             5.903,6.606,6.606,6.606,6.606/
      data N_IRR_LAI/0,0,0.2452,0.3432,0.4451,0.5463,0.6621,
     $             0.7813,0.8868,0.9978,1.124,1.268,1.474,
     $             1.739,2.738,5.349,6.062,6.543,6.543,6.543/
c
c#######################################################################
c
c     irrigation look-up table
c
c#######################################################################
c
cc      data m_grow /3,7,4,9,4,10,
cc     $             3,8,4,9,4,10,
cc     $             3,8,4,9,4,10,
cc     $             3,7,4,9,4,10,
cc     $             3,8,9999,9999,3,9,
cc     $             1,6,6,10,4,11,
cc     $             1,6,6,9,4,11/
cc      data d_grow /16,20,18,10,25,20,
cc     $             22,4,24,10,8,5,
cc     $             25,4,26,10,18,1,
cc     $             18,25,20,10,18,15,
cc     $             22,4,9999,9999,28,21,
cc     $             9,10,13,1,21,11,
cc     $             4,5,11,13,18,15/
c#######################################################################
c
c     Common blocks.
c
c#######################################################################
c
	common /date1  / year, month, day, hour, idc,ihc
     	common /date2  / hydroyear,startmonth,startday,endmonth,endday
	common /date3  / dayinmonth
	common /simulation/  start, finish, dt, start_sub, end_sub
	common /river1    /  nsub, subbasin
	common /river2    /  nflow, dx, Dr
	common /river3    /  s0, b, roughness
      common /river4    /  psubbasin,nbasinup ! added by gaobing
	common / weather / rain_daily,tmin_daily, 
     :         	tmax_daily,evap_daily,snow,pre_hour
     	common /LAID   / NDVI,yndvi,LAI
	common /TEST   /basinarea
c
	common /grid_attrib/ ngrid,grid_row,grid_col
	common /slope     /  area, length, slp, Ds, Dg
	common /soil_land /  nsoil,soiltyp,nland,
     :                	landtyp,land_ratio,soil,land_use
	common /veg_para / Kcanopy,LAImax,root,anik
     :                    ,Sstmax,surfn,Kcrop
	common /soil_para /  wsat,wfld,wrsd,alpha,
     :                    watern,ksat1,ksat2,kg,GWcs
c
	common /location  /  para_dir,data_dir,result1_dir,dem_dir,
     :                     result2_dir,simul_dir
      common /lateralflow/ qin
	common / river_status / Drw,discharge   ! depth of river water        

cc	common /reservoir/ re_capacity, re_capacity0    
	
	common / asc / ncols,nrows,xllcorner,yllcorner,cellsize,
     :           nodata, nnr, nnc, x0, y0, gridsize,znodata
      real snow1(nrow,ncol)
             
c###########added by mqh ,to export soil water############################################################  
       integer countt2(nc),ii1   !计算网格数
      real soil_water1(nc),soil_water2(nc),soil_water3(nc),
     $soil_water4(nc),soil_water5(nc),soil_water6(nc)
      save soil_water1,soil_water2,soil_water3,soil_water4,soil_water5,
     :soil_water6
      save countt2
      real tmp19921,tmp19922,tmp19923,tmp19924,tmp19925,tmp19926
      real result_soilwater1,result_soilwater2,result_soilwater3,
     :result_soilwater4,result_soilwater5,result_soilwater6  !土壤湿度值.



      integer luolicun(nc)  !存储哪些是罗李村的上游站
      
c###########end of add############################################################	
    
      
      
c      save soil_water1,soil_water2,soil_water3,
c     :soil_water4,soil_water5,soil_water6
c      save soil_water12,soil_water22,soil_water32,
c     :soil_water42,soil_water52,soil_water62

      data luolicun /1,1,0,0,0,1,1,1,0,0,1,1,1,1,1,1,0,1,0,
     :1,1,1,1,0,0,1,1,0,1,0,1,1,
     : 1,1,1,1,0,1,0,1,1,0,0,1,0,1,1,1,1,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,
     :0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,
     :0,1,0,0,0,0,1,1,1,1,1,1,1,1,0,1,1,0,0,1,0,0,1,0,0,1,0,1,1,1,1,1,1,
     :1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
     
     

c       soil_water12=0,soil_water22=0,soil_water32=0
c      soil_water42=0,soil_water52=0,soil_water62=0
      
      
!        real,allocatable :: snow1(:,:)
!        allocate (snow1(nrow,ncol))
       
       soil_water1(isub)=0
       soil_water2(isub)=0
       soil_water3(isub)=0
      soil_water4(isub)=0
      soil_water5(isub)=0
      soil_water6(isub)=0
     
     
 
c
c***********************************************************     
c                   Initiate general parameters
c***********************************************************
c
      smf  = 0.1     ! snowmelting factor  !mqh,from 0.15 to 0.1
c	GWcs = 0.15     ! groundwater storage 
c	print *,"hillslope"	,smf,GWcs
c
c
c***********************************************************     
c                   vertical discretization
c***********************************************************
c
c       print * ,'ok0 of hillslope'

       if(start .eq. 1 .and. isub.eq.start_sub) then

       do ir = 1, nrow
         do 10 ic = 1, ncol

           isoil=soil(ir,ic)

c
c--------------- equal interval-----------------
c
c           layer(ir,ic) = 8
c           do j = 1, 8
c              D(ir,ic,j) = Ds(ir,ic)/8.0
c           end do
c
c--------------- variable interval--------------
c                0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0
c
c          print * ,'isoil=,ncol=,soil()', isoil,ir, ic,soil(ir,ic)

c           D(ir,ic,1)   = 0.05
            D(ir,ic,1)   = 0.05  !changed by mqh
           layer(ir,ic) = 1
c           print *, 'layer(ir,ic)=', layer(ir,ic)

           tmp          = D(ir,ic,1)
           do j = 2,10
             if(j.le.4)              D(ir,ic,j) = D(ir,ic,j-1) + 0.05
             if(j.gt.4 .and. j.le.5) D(ir,ic,j) = D(ir,ic,j-1) + 0.10
             if(j.gt.5 .and. j.le.7) D(ir,ic,j) = D(ir,ic,j-1) + 0.20
             if(j.gt.7 .and. j.le.8) D(ir,ic,j) = D(ir,ic,j-1) + 0.30
             if(j.gt.8)              D(ir,ic,j) = D(ir,ic,j-1) + 0.50
             tmp          = tmp + D(ir,ic,j)
             layer(ir,ic) = j
             if(tmp .ge. Ds(ir,ic))    goto 20
           end do
 20        D(ir,ic,j) = D(ir,ic,j) + Ds(ir,ic)-tmp


           if(D(ir,ic,j) .lt. D(ir,ic,j-1)-0.05) then
             D(ir,ic,j-1) = 0.5*(D(ir,ic,j)+D(ir,ic,j-1))
             D(ir,ic,j)   = D(ir,ic,j-1)
           endif
c

cc           if(area(ir,ic).eq.-9999.or.isoil.eq.-9999.
cc     $	                           or.ksat1(isoil).eq.0) goto 10 

           if(area(ir,ic).eq.-9999.or.isoil.eq.-9999) goto 10

		 ! in case of mistake in the soil map
 
c            print * ,'ok004 of hillslope'


           if(layer(ir,ic) .lt. 4)  print *,
     :     "Ds = ",Ds(ir,ic),tmp,"layer = ",layer(ir,ic),isoil, ir, ic

c
c--------------- compute the k0 for each layer --------------
c
           tmp = 0.0
		 do j = 1, layer(ir,ic)
             tmp = tmp + D(ir,ic,j)

c	       f   = -alog(ksat2(isoil)/ksat1(isoil))/5.0
             f   = -alog(ksat2(ir,ic)/ksat1(ir,ic))/Ds(ir,ic)
             k0(ir,ic,j) = ksat1(ir,ic)*exp(-f*tmp)
c		   print *, j, k0(ir,ic,j), D(ir,ic,j), tmp, Ds(ir,ic)
           end do
           if(abs(tmp-Ds(ir,ic)) .gt. 0.01)     print *,  
     :           'wrong in the vertical discretization', tmp, Ds(ir,ic),
     :            layer(ir,ic), D(ir,ic,layer(ir,ic))

 10      continue
 
       end do
      end if
c
c      print * ,'ok1 of hillslope'
c        pause

c

c***********************************************************     
c        read soil and groundwater initial conditions
c***********************************************************
c
      if(start.eq.1) then
c--------------- read from the files --------------
c
	  if(inicon.eq.1) then

	   	call strlen(simul_dir,ia1,ia2)
		open(1,file=simul_dir(ia1:ia2)//subbasin(isub)//
     :		         'I_soil2',status='old')
          do iflow = 1, nflow(isub)
	    do ig = 1, ngrid(isub,iflow)
	      ir    = grid_row(isub,iflow,ig)
	      ic    = grid_col(isub,iflow,ig)
	      isoil = soil(ir,ic)
            read(1,*) ia1, ia2, ib1, ib2
            do iland = 1, ib1
              if(land_ratio(ir,ic,iland) .gt. 0.0) 
     :        read (1,*) (w(ir,ic,iland,j),j=1,ib2), 
     :                    Dgl(ir,ic,iland), Gwst(ir,ic,iland)
              if(Gwst(ir,ic,iland) .ge. Dg(ir,ic)*GWcs(ir,ic)-0.1E-2) 
     :                    Dgl(ir,ic,iland) = Ds(ir,ic)
	        if(Gwst(ir,ic,iland) .le. 0.0) Gwst(ir,ic,iland) = 0.0
	        if(Dgl(ir,ic,iland)  .ge. Ds(ir,ic)+Dg(ir,ic)) 
     :                        Dgl(ir,ic,iland) = Ds(ir,ic)+Dg(ir,ic)
            end do
          end do
          end do
		close(1)
c
c--------------- specify by the following way --------------
	   else
          if(isub.eq.start_sub) then
		  do ir=1,nrow
		    do ic=1,ncol
		      isoil=soil(ir,ic)
                if(isoil .eq. -9999) goto 28
			  do iland=1,nland
			    Dgl(ir,ic,iland)  = Ds(ir,ic)
 	            GWst(ir,ic,iland) = (Ds(ir,ic)+Dg(ir,ic)
     :		                       -Dgl(ir,ic,iland))*Gwcs(ir,ic)
			    tmp = 0.0
                  D0  = Ds(ir,ic)*wrsd(ir,ic)/(wsat(ir,ic)-wrsd(ir,ic))
			    do j= 1, layer(ir,ic)
			      tmp = tmp+D(ir,ic,j)
			      w(ir,ic,iland,j) = 
     :                      wsat(ir,ic)*(D0+tmp)/(D0+Ds(ir,ic))
	              if(isoil.eq.-9999) w(ir,ic,iland,j) = wfld(ir,ic)

	              if(w(ir,ic,iland,j)-wsat(ir,ic) .gt. 0.01) 
     :              print *, wsat(ir,ic),w(ir,ic,iland,j), D0, tmp, 
     :              Ds(ir,ic),'33333333'

	              if(w(ir,ic,iland,j) .lt. wrsd(ir,ic)+ 0.001) 
     :              print *, wrsd(ir,ic),
     :              w(ir,ic,iland,j), D0, tmp, Ds(ir,ic), '9999'

			    end do
		  	  end do
28              continue
		    end do
		  end do
          end if
	  end if
      end if
c
c       print * ,'ok2 of hillslope'
c        pause

c***********************************************************     
c            sub-basin area 面积
c***********************************************************
c
      if(start.eq.1) then
	  if(isub.eq.start_sub) basinarea = 0
        subarea(isub) = 0 
        do iflow = 1, nflow(isub)
		do ig = 1, ngrid(isub,iflow)
		  ir = grid_row(isub,iflow,ig)
		  ic = grid_col(isub,iflow,ig)
		  if(area(ir,ic) .eq. -9999.0) print*,'wrong in grid-area'
		  basinarea = basinarea + area(ir,ic)
		  subarea(isub) = subarea(isub) + area(ir,ic)
		end do
	  end do
	end if

!        print *,basinarea
!        pause

        !end of add
	
c
c
c        print * ,'ok3 of hillslope'
c***********************************************************     
c                     initialize status variables
c***********************************************************
c
	if(isub.eq.start_sub .and. start.eq.1 ) then
        do ir = 1, nrow
         do  ic = 1, ncol
           evap_month(ir,ic) = 0.0
           soil_month(ir,ic) = 0.0
           runoff_month(ir,ic) = 0.0
           grunoff_month(ir,ic) = 0.0
           srunoff_month(ir,ic) = 0.0
	     drunoff_month(ir,ic) = 0.0
          	ecy_month(ir,ic)  = 0.0
	        ecp_month(ir,ic)  = 0.0
              ese_month(ir,ic)  = 0.0

           do iland = 1,nland
             Cst(ir,ic,iland) = 0.0
             Sst(ir,ic,iland) = 0.0
          end do
         end do
        end do
      end if

	if(isub.eq.start_sub .and. (start.eq.1 .or. 
     :           (idc.eq.1 .and.ihc.eq.1))) then
        do ir = 1, nrow
         do 50 ic = 1, ncol
           do i=1,366
             raind(ir,ic,i)     = 0.0
cc	       epd(ir,ic,i)       = 0.0
	       eactd(ir,ic,i)     = 0.0
	       ecy(ir,ic,i)       = 0.0
	       ecp(ir,ic,i)       = 0.0
	       ese(ir,ic,i)       = 0.0
	       runoffd(ir,ic,i)   = 0.0
             srunoff(ir,ic,i)   = 0.0
	       groundoff(ir,ic,i) = 0.0
	       soiloff(ir,ic,i)   = 0.0
		 end do
 50      continue
        end do

	  annual_Cst0  = 0.0
	  annual_Ssts0 = 0.0
	  annual_Sstr0 = 0.0
	  annual_SBst0 = 0.0
	  annual_Gst0  = 0.0
cc	  annual_re0   = 0.0
	end if
c  
	if( start .eq. 1 .or. (idc.eq.1 .and. ihc.eq.1)) then 
        do iflow = 1, nflow(isub)
cc          annual_re0 = annual_re0 + re_capacity(isub,iflow)
		do ig = 1, ngrid(isub,iflow)
	      ir = grid_row(isub,iflow,ig)
	      ic = grid_col(isub,iflow,ig)
c	      isoil = soil(ir,ic)
            do  iland = 1, nland
              if(land_ratio(ir,ic,iland) .gt. 0.0) then
              annual_Cst0 = annual_Cst0 + Cst(ir,ic,iland)*
     :                      0.001*land_ratio(ir,ic,iland)*area(ir,ic) !in m3
              annual_Sstr0= annual_Sstr0 + Sst(ir,ic,iland)*
     :                      0.001*land_ratio(ir,ic,iland)*area(ir,ic) !in m3
              annual_Ssts0= annual_Ssts0 + snow(ir,ic,iland)*
     :                      0.001*land_ratio(ir,ic,iland)*area(ir,ic) !in m3
              annual_Gst0 = annual_Gst0 + GWst(ir,ic,iland)*
     :                      area(ir,ic)*land_ratio(ir,ic,iland)       !in m3
              do j = 1 ,layer(ir,ic)
                annual_SBst0 = annual_SBst0 + w(ir,ic,iland,j) *
     :            D(ir,ic,j)*land_ratio(ir,ic,iland)*area(ir,ic)      !in m3
              end do
	        end if
            end do 
          end do
	  end do
	end if
c	print *,basinarea
c
c        print * ,'ok4 of hillslope'
c***********************************************************     
c                     initialize irrigation variables
c***********************************************************

      if(start .eq. 1) return
c
c****************************************************************
c          start simulation   
c****************************************************************
c
c
c        print * ,'ok5 of hillslope'
c****************************************************************
c          interception, evapotranspiration   
c****************************************************************
c
      do 999 iflow=1,nflow(isub)
c	  print *,iflow,"ngrid=", ngrid(isub,iflow) 
        do ig=1,ngrid(isub,iflow)
	    ir=grid_row(isub,iflow,ig)
	    ic=grid_col(isub,iflow,ig)
	    isoil=soil(ir,ic)
c           changged by gaobing


       
c	      if(isub.ge.10 .and. isub.le.18)
c              fice = (temper+10.0)/10.0
c            fice = amax1(fice, 0.1)
c            fice = amin1(fice, 1.0)
c
c****************************************************************
c variable using for water balance calculation
c****************************************************************
c
c            raind(ir,ic,idc) = raind(ir,ic,idc) + prec
c	      epd(ir,ic,idc)   = epd(ir,ic,idc)   + Ep
c
c****************************************************************

          do 888 iland = 1, nland
            call rain_model(isub, ir, ic, hour, month,
     :                          rain_daily(ir,ic,day),         
     :                          tmin_daily(ir,ic,day),
     :                          tmax_daily(ir,ic,day), 
     :                          evap_daily(ir,ic,day,iland), !
     :                          prec, temper, Ep, iland)
c            print * , month,day,hour,temper
c****************************************************************
c variable using for water balance calculation
c****************************************************************
c
c      changed by xujijun
c	if(isub.ge.10 .and. isub.le.18 )          prec = 3*prec
cc      
            prec = pre_hour(ir,ic,day,hour)
c	      
c         if(year == 2000)  prec = prec*1.4  !标记，什么意思？原来是有的
cc	print *, prec, ir, ic, day, hour
		 
		  if(iland.eq.1) then
              raind(ir,ic,idc) = raind(ir,ic,idc) + prec
		  end if
cc	      if(land_ratio(ir,ic,iland).ge.0) then
cc		    epd(ir,ic,idc)=epd(ir,ic,idc)+Ep*land_ratio(ir,ic,iland) 
cc	      end if

	      pnet           = prec
	      Eact           = 0.0
	      if(land_ratio(ir,ic,iland).le.0.0) goto 888


ccccccccccccccccccccccccc   added by gaobing ccccccccccccccccc
c         print *,j
         !标记，这里本来都是被注释掉的
c        if(iland.eq.1 .or. iland.eq.10) then  !waterbody or frost snow
c	         do kn = 1,4
c	         LAI(ir,ic,j,kn) = 0.0
c	         enddo
c		  end if 
 
c          if(iland.eq.2) then
c             c4 = 0.06
c	    elseif( iland .eq. 3) then
c	       c4=0.2
c          elseif( iland .eq. 4) then
c 	       c4=0.27
c             c4=10.00
c	    elseif( iland .eq. 6) then
c	       c4=0.4
c	    elseif( iland .eq. 7) then
c	       c4=0.3
c	    elseif( iland .eq. 8) then
cc	       c4=0.08
c            c4=0.40
c	    elseif( iland .eq. 9) then
c	       c4=0.3
 
c	    endif
	 

	        if(iflai .eq. 0 ) then
c
c****************************************************************
c          NDVI -> LAI  
c****************************************************************
c      
        
               
		do j=1,12
	        if(year >=1998) then
	     do kn = 1,3	     
	        NDVI(ir,ic,j,kn)=yndvi(ir,ic,j,kn)
	    enddo
	       else
              do kn = 1,2	     
	        NDVI(ir,ic,j,kn)=yndvi(ir,ic,j,kn)
	       enddo
	        NDVI(ir,ic,j,3)=yndvi(ir,ic,j,2)
		   endif

           do kn =1,3
		  if(NDVI(ir,ic,j,kn) .le. 0.0)   NDVI(ir,ic,j,kn) = 0.025
		  if(NDVI(ir,ic,j,kn) .gt. 0.975) NDVI(ir,ic,j,kn) = 0.975
           

cc       changed by xujijun start

            do i_ndvi=1, 19
	        if(NDVI(ir,ic,j,kn) .le. T_NDVI(1)) ii_ndvi = 1
	        if(NDVI(ir,ic,j,kn).gt.T_NDVI(i_ndvi) .and. 
     $           NDVI(ir,ic,j,kn).le.T_NDVI(i_ndvi+1)) then
                 ii_ndvi = i_ndvi
			end if    
	        if(NDVI(ir,ic,j,kn) .gt. T_NDVI(20)) ii_ndvi = 20
		  end do

	      c4 = 0.4  !changed from 0.4 to 0.3
c            c4=0.6             ! changed by gaobing
cc		  if(iland.eq.1 .or. iland.eq.19) then  !waterbody or frost snow
		  if(iland.eq.1 .or. iland.eq.10) then  !waterbody or frost snow
	         LAI(ir,ic,j,kn) = 0.0
		  end if 

	      if(iland.eq.2) then  !urban-area
               if(ii_ndvi.eq.0 .or. ii_ndvi.eq.20) then
	           LAI(ir,ic,j,kn)=Shurb_LAI(ii_ndvi)           
			  else
	           LAI(ir,ic,j,kn) = Shurb_LAI(ii_ndvi) + 
     $           (NDVI(ir,ic,j,kn)-T_NDVI(ii_ndvi))*
     $           (Shurb_LAI(ii_ndvi+1)-
     $           Shurb_LAI(ii_ndvi))/(T_NDVI(ii_ndvi+1)-T_NDVI(ii_ndvi))           
			 end if
	         c4 = 0.06
		  end if 

	      if(iland.eq.3) then  !baresoil
c	         LAI(ir,ic,j) = 1.71*NDVI(ir,ic,j)+0.48
	         LAI(ir,ic,j,kn) = 0.2*NDVI(ir,ic,j,kn)+0.2
	         c4=0.4 ! changed by gaobing
		  end if 

cc	      if(iland.eq.4 .or. iland.eq.5) then  !forest
	      if(iland.eq.4) then  !forest
                 if(ii_ndvi.eq.0 .or. ii_ndvi.eq.20) then
	           LAI(ir,ic,j,kn) = Forest_LAI(ii_ndvi)           
			  else
	           LAI(ir,ic,j,kn)=Forest_LAI(ii_ndvi)+(NDVI(ir,ic,j,kn)
     $           -T_NDVI(ii_ndvi))*(Forest_LAI(ii_ndvi+1)-Forest_LAI
     $           (ii_ndvi))/(T_NDVI(ii_ndvi+1)-T_NDVI(ii_ndvi))           
			 end if
	         if(iland.eq.4) c4 = 0.2
cc	         if(iland.eq.5) c4 = 0.08
		  end if 

cc	      if(iland.ge.6 .and. iland.le.12) then  !irrigated-cropland

	      if(iland.ge.5) then  ! sparse vegetion
c               if(ii_ndvi.eq.0 .or. ii_ndvi.eq.20) then
c	           LAI(ir,ic,j,kn)=IRR_LAI(ii_ndvi)           
c			  else
c	           LAI(ir,ic,j,kn)=IRR_LAI(ii_ndvi)+(NDVI(ir,ic,j,kn)-
c     $           T_NDVI(ii_ndvi))*(IRR_LAI(ii_ndvi+1)-
c    $           IRR_LAI(ii_ndvi))/(T_NDVI(ii_ndvi+1)-T_NDVI(ii_ndvi))           
c			 end if
c	         c4 = 0.45
           LAI(ir,ic,j,kn)=1.71*NDVI(ir,ic,j,kn)+0.48
           c4 = 0.3
		  end if 

cc	      if(iland.eq.13) then  !upland
	      if(iland.eq.6) then  !upland
               if(ii_ndvi.eq.0 .or. ii_ndvi.eq.20) then
	           LAI(ir,ic,j,kn)=N_IRR_LAI(ii_ndvi)           
			  else
	           LAI(ir,ic,j,kn)=N_IRR_LAI(ii_ndvi)+(NDVI(ir,ic,j,kn)-
     $           T_NDVI(ii_ndvi))*(N_IRR_LAI(ii_ndvi+1)-
     $           N_IRR_LAI(ii_ndvi))/(T_NDVI(ii_ndvi+1)-T_NDVI(ii_ndvi))           
			 end if
		  end if 

cc	      if(iland.ge.14 .and. iland.le.16) then  !grassland

	      if(iland.eq.7) then  !grassland
	         LAI(ir,ic,j,kn)=1.71*NDVI(ir,ic,j,kn)+0.48
c			 if(isub .le. 23) then  !Taking special value of LAI for the high plateau
c	           LAI(ir,ic,j)=-alog((1.0-NDVI(ir,ic,j)/0.915)
c     $                            /0.83)/0.96
c	           if(LAI(ir,ic,j) .lt. 0.0)  LAI(ir,ic,j) = 0.0
c	           if(LAI(ir,ic,j) .gt. 2.18) LAI(ir,ic,j) = 2.18
c			 end if 
cc	         if(iland.eq.14) c4 = 0.6
cc	         if(iland.eq.15) c4 = 0.3
cc	         if(iland.eq.16) c4 = 0.1
	         if(iland.eq.7) c4 = 0.3
c              if(iland.eq.7) c4 = 0.6
		  end if 

cc	      if(iland.eq.17) then  !shrub

	      if(iland.eq.8) then  !shrub
               if(ii_ndvi.eq.0 .or. ii_ndvi.eq.20) then
	           LAI(ir,ic,j,kn)=Shurb_LAI(ii_ndvi)           
			  else
	           LAI(ir,ic,j,kn)=Shurb_LAI(ii_ndvi)+(NDVI(ir,ic,j,kn)-
     $           T_NDVI(ii_ndvi))*(Shurb_LAI(ii_ndvi+1)-
     $           Shurb_LAI(ii_ndvi))/(T_NDVI(ii_ndvi+1)-T_NDVI(ii_ndvi))           
			 end if
	         c4 = 0.08
		  end if 

cc	      if(iland.eq.18) then  !wetland

	      if(iland.eq.9) then  !wetland
	         LAI(ir,ic,j,kn) = 1.71*NDVI(ir,ic,j,kn)+0.48
		  end if 
c
            if(LAI(ir,ic,j,kn) .lt. 0.0) LAI(ir,ic,j,kn) = 0.0
		  if(LAI(ir,ic,j,kn) .gt. LAImax(iland)) LAImax(iland) = 
     :                                              LAI(ir,ic,j,kn)
	enddo ! kn
c	     c4=0.2           ! changed by gaobing
          end do


c
c	print *, iland, NDVI(ir,ic,month),LAI(ir,ic,month), month, ir, ic
c
           if(year >=1998) then
           	        
              if(day.lt. 10) then
                dLAI=LAI(ir,ic,month,1)
	         elseif (day. gt. 20) then
                dLAI=LAI(ir,ic,month,3)
	         else
	         dLAI=LAI(ir,ic,month,2)
              end if 
             
	      else

	        if(day.lt. 15) then
                dLAI=LAI(ir,ic,month,1)
	         
	         else
	         dLAI=LAI(ir,ic,month,2)
              end if 

	    endif

	    else ! IFLAI == 1
c         print *, iflai, 'iflai==1'

            if(iland .eq. 7) then
	        do kn = 1,4
              if (LAI(ir,ic,month,kn) > LAImax(iland) ) then
                 LAI(ir,ic,month,kn)=LAImax(iland)-0.2
	        endif
	        enddo
	     endif
            if(iland .eq. 9) then
	          do kn =1,4
              if (LAI(ir,ic,month,kn) > LAImax(iland) ) then
	             LAI(ir,ic,month,kn)=LAImax(iland)-0.2
	        endif
	          enddo
	      endif
	     if(iland .eq. 8) then
	        do kn = 1,4
              if (LAI(ir,ic,month,kn) > LAImax(iland) ) then
                 LAI(ir,ic,month,kn)=LAImax(iland)-0.2
	        endif
	        enddo
	     endif

	    if(iland .eq. 3) then
	        do kn = 1,4
              if (LAI(ir,ic,month,kn) > LAImax(iland) ) then
                 LAI(ir,ic,month,kn)=LAImax(iland)-0.2
	        endif
	        enddo
	     endif

	 if(iland .eq. 1 .or. iland .eq . 6) then
	        do kn = 1,4
              if (LAI(ir,ic,month,kn) > LAImax(iland) ) then
                 LAI(ir,ic,month,kn)=LAImax(iland)
	        endif
	        enddo
	     endif

	 if(iland .eq. 10) then
	        do kn = 1,4
              if (LAI(ir,ic,month,kn) > LAImax(iland) ) then
	        
                 LAI(ir,ic,month,kn)=LAImax(iland)
	      
	        endif
	        enddo
	     endif

	    if(iland .eq. 2) then
	        do kn = 1,4
              if (LAI(ir,ic,month,kn) > LAImax(iland) ) then
                 LAI(ir,ic,month,kn)=LAImax(iland)
	        endif
	        enddo
	     endif

            if(day.lt. 8) then
                dLAI=LAI(ir,ic,month,1)	
			        
	      elseiF(  day .lt. 16 .and. day .ge. 8) then
                dLAI=LAI(ir,ic,month,2)
	      elseif(  day .lt. 24 .and. day .ge. 16) then
                dLAI=LAI(ir,ic,month,3)
	      else
	         dLAI=LAI(ir,ic,month,4)
            end if 
           

	    endif
c
c****************************************************************
c          adjust irrigation LAI based on plant growth stage  
c****************************************************************
c   
c           goto 111  d_grow(1,iland-5)

cc	111	删除于2006.4.14	 !skip d_grow

c
c****************************************************************
c          interception  
c****************************************************************
c   
c		  Cstmax(ir,ic,iland) = 0.20*dLAI   !LAI(ir,ic,month)
             Cstmax(ir,ic,iland) = 0.10*dLAI
cc	      if(iland.eq.4 .or. iland.eq.5 .or.iland.eq.17) then !forest or shurb
	      if(iland.eq.4 .or. iland.eq.8) then !forest or shurb
c		    Cstmax(ir,ic,iland) = 0.25*dLAI !LAI(ir,ic,month)
             Cstmax(ir,ic,iland) = 1.00*dLAI 
		  end if

c		  Cstmax(ir,ic,iland) = 0.2*LAI(ir,ic,month)

		  if(pnet.le. 0.0) goto 250

		  DeficitCst = Cstmax(ir,ic,iland) - Cst(ir,ic,iland)
		  DeficitCst = amax1(DeficitCst,0.0)
		  if(pnet .gt. DeficitCst) then
	        Cst(ir,ic,iland) = Cst(ir,ic,iland) + DeficitCst
			pnet = pnet - DeficitCst 
	       else
	        Cst(ir,ic,iland) = Cst(ir,ic,iland)+pnet
			pnet = 0.0 
		  end if 
c
c****************************************************************
c          calculate snow depth on the ground(mm)
c
c		 calculate snow melt   
c****************************************************************
c		 
c            if(snow(ir,ic,iland).gt.0.0) then
c		    print *, 'ir, ic, iland, snow(ir,ic,iland), temper, pnet'
c              print *, ir, ic, iland, snow(ir,ic,iland), temper, pnet
c             pause
c		  endif
		     
		  if(temper .le. 1.0) then
			snow(ir,ic,iland) = snow(ir,ic,iland) + pnet
			pnet=0.0 

c          added by xujijun	- start	 
c	  if(snow(ir,ic,iland).ge.100) snow(ir,ic,iland) = 100
c          added by xujijun	- end	 
            end if
		  snowmelt=0.0
c		  if(month .ge. 2 .and. 

c		  print *, 'isub, month, ir, ic, pnet'
c		  print *, isub, month, ir, ic, pnet
c           changed by gaobing
		  if(temper .gt. 0.0 
     :                      .and. snow(ir,ic,iland).gt.0.0) then
c          added by xujijun	- gaobing
		  if(month.ge.3.and.month.le.10) smf= 0.15
            if(month.ge.5.and.month.le.8)  smf = 0.15
           
	      

c          added by xujijun	- end	 
	        snowmelt = (smf+pnet/20.0) * (temper-1.5)
c               snowmelt = (smf+pnet/10.0) * temper
	        snowmelt = amin1(snowmelt,snow(ir,ic,iland))
             	   
			snow(ir,ic,iland)=snow(ir,ic,iland)-snowmelt	       
		
			pnet=pnet+snowmelt 
	         
            
c		  if(snow(ir,ic,iland).gt.0.0) then
c		   print *, 'smf, snowmelt, pnet, snow' 
c		   print *, smf, snowmelt, pnet, snow(ir,ic,iland)
c            endif
           
		  end if
            
c	      if(snow(ir,ic,iland) > 60) then  ! added vy gaobing 
c               	snow(ir,ic,iland)=60
c	            snowmelt=snowmelt+snow(ir,ic,iland)-60
c	            pnet=pnet+snow(ir,ic,iland)-60
c	         endif

c        added by gaobing	- end
            Sst(ir,ic,iland) = Sst(ir,ic,iland) + pnet
c	      if(pnet.gt.0) print *, "Sst=,p",Sst(ir,ic,iland),pnet,ir,ic
c 
c****************************************************************
c          (1) evaporation from canopy storage  
c****************************************************************
c    

250         continue
            EfromCanopy = 0.0
            EfromCrop   = 0.0
	      EfromSurface= 0.0


!
!           c2=0.05+0.1*dLAI/LAImax(iland)

           c2=0.05+0.1*dLAI/LAImax(iland)
	      c1=0.31*dLAI

            !c2=0.05+0.10*dLAI/LAImax(iland)
	      !c1=0.31*dLAI

c           changed by xujijun

	      c3 =  0.23+0.1*(dLAI/LAImax(iland))**0.5		! 
		  c5 =	0.23
		if (luolicun2(isub).eq.1) then
		c1=c1*0.5
		c2=c2*0.5
		c3=c3*0.5
		c5=c5*0.4
		endif
		  
        		  


c		if(temper .le. 0.0) Ep=0.0
            Etr = Ep*amin1(c2+c1,1.0)
            Es  = (Ep*(c2+(1-c2)*(1-amin1(c2+c1,1.0))) 
     :             -Etr*(1-amin1(c2+c1,1.0)))
     :             *(1-exp(-c4*dLAI))
     :             +Ep*exp(-c4*dLAI)*(c5+(1.0-c5)*dLAI/LAImax(iland))
     :              **(1.0+c3)
c            Es  = (Ep*(c2+(1-c2)*(1-amin1(c2+c1,1.0))) 
c     :             -Etr*(1-amin1(c2+c1,1.0)))
c     :             *Kcanopy(iland)
c     :             +Ep*(1-Kcanopy(iland))*(dLAI/LAImax(iland))
c     :              **(1.0+c3)
	      
            if(iland.eq.3) then  !baresoil
              Es  = Ep*exp(-c4*dLAI)
c              Es  = (Ep*(c2+(1-c2)*(1-amin1(c2+c1,1.0))) 
c     :             -Etr*(1-amin1(c2+c1,1.0)))
c     :             *Kcanopy(iland)
c     :             +Ep*(1-Kcanopy(iland))
		  end if
		     
            if(iland.eq.1) then  !waterbody
              Es  = Ep           !*(1-Kcanopy(iland))
		  end if   
              Etr = Etr*(1-exp(-c4*dLAI))      !Kcanopy(iland)

c            Etr = Ep*kcrop(iland)*Kcanopy(iland)*
c     :                            LAI(ir,ic,month)/LAImax(iland)
c            Es  = Ep*kcrop(iland) - Etr

c		  tmp = 1.0
c		  if(isub.le.23) then
c		    tmp = fice*0.8
c		    if(year.gt.1960) tmp = fice*0.82
c		    if(year.gt.1970) tmp = fice*0.84
c		    if(year.gt.1980) tmp = fice*0.86
c		    if(year.gt.1990) tmp = fice*0.88
c	      endif
c		  if(isub.ge.10.and. isub.le.69
c     :          .and.iland.ge.6.and.iland.le.12) then
c               tmp =1.05 + 0.05*sin((month-3)*3.1415926/6.0
c     :                - min(2,max(month-6,0))*3.1415926/6.0)
c	      endif
c		  if(isub.gt.69.and.iland.ge.6.and.iland.le.12) then
c               tmp =1.08 + 0.08*sin((month-3)*3.1415926/6.0
c     :                - min(2,max(month-6,0))*3.1415926/6.0)
c	      endif
c
c		  Etr = tmp * Etr
c		  Es  = tmp * Es

         do kn=1,3
	if(kcrop(iland).gt.1.or.(LAI(ir,ic,month,kn)-LAImax(iland)).gt.1e-6)
     :  print *,"kcrop out of range", kcrop(iland), LAI(ir,ic,month,kn), 
     :                                LAImax(iland), iland
         enddo
!		  if(prec. gt. 0.01) goto 290   !mqh

            if(Etr .lt. 0.1E-9) goto 270

		  EfromCanopy = amin1(Etr, Cst(ir,ic,iland)) 
            Cst(ir,ic,iland) = 
     :                  Cst(ir,ic,iland) - EfromCanopy
c	      if(prec .gt. 0.0001) then
c	     print *,"Canopy",EfromCanopy,Etr,Es,Cst(ir,ic,iland),idc,ihc
c	      end if
            if(EfromCanopy.lt.0.0) 
     $      print *,'wrong in E from canopy',ihc,iflow,EfromCanopy
c
c****************************************************************
c          (2) transpiration from vegetation
c****************************************************************
c
            if(root(iland).eq.0.0) goto 270
            i      = 0          ! layer of root	     
            para_r = 0.1        ! root density ratio of deepest layer (<=1.0)
            tmp    = 0.0
            do j = 1, layer(ir,ic)
              tmp = tmp + D(ir,ic,j)
              i   = i + 1
              if(tmp .ge. root(iland)) goto 266
            end do
 266        if(layer(ir,ic) .le. i) i = layer(ir,ic)
            para_a = 1.0/(float(i) - 0.5*(1.0-para_r)*float(i-1))
		      
            do j=1,i
              y = (1.0-(1.0-para_r)*float(j-1)/
     :             float(i))*para_a

              call ETsoil(w(ir,ic,iland,j), 
     :                    wfld(ir,ic), wrsd(ir,ic), kmoist)
              tmp = y * Etr * kmoist
			    
            if(tmp .lt. 0.0) print *,isub,iflow,isoil,iland,j,tmp,'4545'
			temp = amax1(0.0, 
     :                        1000.0*(w(ir,ic,iland,j)-
     :                        wrsd(ir,ic)-1.0E-6)*D(ir,ic,j))
	        tmp=amin1(tmp,temp) 
              w(ir,ic,iland,j) = w(ir,ic,iland,j) - 
     :                           tmp/(1000.0*D(ir,ic,j))
c            由mqh注释   标记1
c	        if(w(ir,ic,iland,j).lt.wrsd(ir,ic)+0.1E-6) then
c	          print *,"crop evap from w",w(ir,ic,iland,j),j,'11111'
c	          pause
c			end if
              EfromCrop = EfromCrop + tmp
c	      if(prec .gt. 0.0001) then
c	     print *,"Crop,w",EfromCrop,Etr,Es,w(ir,ic,iland,j),idc
c	      end if
            end do
c
c****************************************************************     
c          (3) evaporation from soil surface
c****************************************************************
c

 270        if(Es .lt. 0.1E-9) goto 290
		  if(Sst(ir,ic,iland) .le. 0.0) goto 275
		  EfromSurface     = amin1(Es, Sst(ir,ic,iland))
		  Sst(ir,ic,iland) = Sst(ir,ic,iland) - EfromSurface
		  if(EfromSurface .lt. -0.0001)
     :	        print *,'wrong in EfromSurface',
     :            Es, EfromSurface,Sst(ir,ic,iland)
       
275         if(Sst(ir,ic,iland) .lt. 0.0) Sst(ir,ic,iland) = 0.0

            tmpEp = Es-EfromSurface 
		  if(tmpEp .le. 0.1E-9) goto 290
            call ETsoil( w(ir,ic,iland,1), 
     :                    wfld(ir,ic), wrsd(ir,ic), kmoist )
            tmp = tmpEp * kmoist
            if(tmp .lt. 0.0) print *, 'wrong in EfromSoil', tmp
            temp = amax1(0.0,  
     :                  1000.0*(w(ir,ic,iland,1)-wrsd(ir,ic)-1.0E-6)
     :                  *D(ir,ic,1))
            tmp = amin1(tmp,temp)  
            w(ir,ic,iland,1) = w(ir,ic,iland,1) - 
     :                       tmp/(1000.0*D(ir,ic,1))
c            print *,w(ir,ic,iland,1)
c           由mqh注释  标记2
c	      if(w(ir,ic,iland,1).lt.wrsd(ir,ic)+0.1E-6) then
c	        print *,"wrong in calculate surface evap",ir,ic,iland,'1111'
c	        pause
c		  end if
            EfromSurface = EfromSurface + tmp
		      
 290        continue
c


c		  print *, 'es, sst', Es, Sst(ir,ic,iland),ir,ic
c	      pause

c****************************************************************
c variable using for water balance calculation
c****************************************************************
c
	      eactd(ir,ic,idc) = eactd(ir,ic,idc) + 
     :                         (EfromCanopy+EfromCrop+EfromSurface) *
     :                          land_ratio(ir,ic,iland)            ! mm
      ecy(ir,ic,idc)   =  ecy(ir,ic,idc) + EfromCanopy 
     :	*land_ratio(ir,ic,iland) 
      ecp(ir,ic,idc)   =  ecp(ir,ic,idc) + EfromCrop 
     :	*land_ratio(ir,ic,iland) 
      ese(ir,ic,idc)   =  ese(ir,ic,idc) + EfromSurface
     :	*land_ratio(ir,ic,iland)                       ! added by gaobing
            eactd(ir,ic,idc) = eactd(ir,ic,idc) + Eact * 
     :                              land_ratio(ir,ic,iland)        ! mm
      !标记
c****************************************************************
c
c	      print *,"eactd=",eactd(isub,iflow,idc),EfromCanopy
c     :	      ,EfromCrop,EfromSurface,idc,ihc,ir,ic,iland
c	      print *,"rain=",prec,"  Ep=",EfromSurface,idc,ihc
c
 888		continue
	  end do
 999  end do
c
c        print * ,'ok6 of hillslope'

c****************************************************************
c          UZ, GW and surface routing
c****************************************************************
c
      srunoff2(isub,ihc)=0
      groundoff2(isub,ihc)=0
      soiloff2(isub,ihc)=0  !mqh
      do 1999 iflow=1, nflow(isub)
	  qin_tmp0 = 0.0

	  do ig=1,ngrid(isub,iflow)
	    ir=grid_row(isub,iflow,ig)
	    ic=grid_col(isub,iflow,ig)
	    isoil=soil(ir,ic)
		
c            call rain_model(hour, rain_daily(ir,ic,idc),         
c     :                          tmin_daily(ir,ic,idc),
c     :                          tmax_daily(ir,ic,idc),
c     :                          evap_daily(ir,ic,idc),
c     :                          prec, temper, Ep)
            fice = 1.00
	   


               fice=(temper+7.5)/7.5
	         fice=amax1(fice, 0.5)
	         fice=amin1(fice , 1.0)
c		  endif

	    qin_tmp=0.0

	    d1_slope  = slp(ir,ic)           ! m/m
	    d1_length = length(ir,ic)       ! meter           *1.3，added by mqh ,2014.4.17
	    d1_Ds     = Ds(ir,ic)            ! meter
	    d1_Dg     = Dg(ir,ic)		     ! meter
	    d1_Dr     = Dr(isub, iflow)		 ! meter
	    d1_Drw    = Drw(isub, iflow)	 ! meter
	    !mqh  added by mqh
c	    if (luolicun(isub).eq.0) kg(ir,ic)=kg(ir,ic)*2
	    d1_kg     = kg(ir,ic)*fice       ! m/sec
          d1_GWcs   = GWcs(ir,ic)          ! non
          d1_dt     = dt	 

	    d1_nuz    = layer(ir,ic)
          do j = 1, d1_nuz
	      d1_deltz(j) = D(ir,ic,j)
	      d1_k0(j)    = k0(ir,ic,j)*fice      ! m/sec
		end do

	    d1_wsat   = wsat(ir,ic)
	    d1_wrsd   = wrsd(ir,ic)
	    d1_wfld   = wfld(ir,ic)
	    d1_alpha  = alpha(ir,ic)
	    d1_watern = watern(ir,ic)

c         print * ,'ok611 of hillslope'

	    do 1800 iland = 1, nland
            if(land_ratio(ir,ic,iland) .le. 0.0) goto 1800
            d1_land_ratio = land_ratio(ir,ic,iland)
	      d1_anik = anik(iland)
		  if( d1_anik .lt. 0.0)  print *, "wrong anik..." 
		  if( d1_anik .lt. 1.0)  print *, "wrong anik",anik(iland) 
	        
 		  do j = 1, d1_nuz
	        d1_w(j) = w(ir,ic,iland,j)

c	        if(d1_wrsd-d1_w(j).gt.0.005) 
c     :		print *,"d1_w(j) < wrsd", d1_w(j), d1_wrsd
c     :                ,ir,ic,j,idc,ihc

c	        if(d1_w(j)-d1_wsat .gt. 0.005) 
c     :          print *,"d1_w(j)>wsat", d1_w(j),d1_wsat,j

	      end do
		
		  d1_Sst = Sst(ir,ic,iland)/1000.0
            d1_Dgl = Dgl(ir,ic,iland)           
            d1_GWst= GWst(ir,ic,iland)
		  d1_snow= snow(ir,ic,iland)
            
cc	      if(month.gt.3) then
c	      if(temper.gt.1.5.and. d1_snow.gt.0.and.snowmelt.gt.0) then
c            print *,' temper4, sst, d1_Dgl, d1_GWst, d1_snow'
c		  print *,temper, Sst(ir,ic,iland), d1_Dgl, d1_GWst, d1_snow
c            print *, pnet, snowmelt
c            endif
cc            endif

c		  print *,"input data", d1_slope,d1_length,d1_Ds,d1_Dg,d1_Dr,
c     :						  d1_Drw,d1_kg,d1_sst,d1_Dgl,d1_GWst
c         print * ,'ok6111 of hillslope'
c     		print *, 'layer(ir,ic), d1_nuz=',layer(ir,ic), d1_nuz  
!            print *,d1_Drw,d1_Dr,d1_Dgl
!            pause
	      call runoff(isub, month, d1_snow, temper,
     :                  d1_land_ratio,d1_slope, d1_length, 
     :                  d1_Ds, d1_Dg ,d1_Dr,d1_anik,d1_wsat,
     :                  d1_wrsd, d1_wfld,d1_alpha ,d1_watern, 
     :                  d1_nuz,d1_deltz, d1_k0,d1_w,d1_kg,d1_GWcs,
     :              d1_GWst,d1_Dgl,d1_Drw,d1_Sst,d1_dt,d1_qsub,d1_ssub)

c          print * ,'ok6112 of hillslope'

		  do j = 1, d1_nuz
	        w(ir,ic,iland,j) = d1_w(j)
	      end do
c            if(1000.0*d1_Sst-Sst(ir,ic,iland) .gt. 0.01) print *,"Sst"
c     :           ,1000.0*d1_Sst,Sst(ir,ic,iland),ir,ic,iland,idc,ihc

		  Sst(ir,ic,iland) = d1_Sst*1000.0
c	      if(temper.gt.1.5.and. d1_snow.gt.0.and.snowmelt.gt.0) then
c            print *,' temper5, sst, d1_Dgl, d1_GWst, d1_snow'
c		  print *,temper, Sst(ir,ic,iland), d1_Dgl, d1_GWst, d1_snow
c            print *, pnet, snowmelt
c            endif  
            Dgl(ir,ic,iland) = d1_Dgl
	      GWst(ir,ic,iland)= d1_GWst
c 		  print *,"output data",d1_slope,d1_length,d1_Ds,d1_Dg,d1_Dr
c     :     ,d1_Drw,d1_kg,d1_sst,d1_Dgl,d1_GWst,land_ratio(ir,ic,iland)

	      qin_tmp = qin_tmp + d1_qsub*land_ratio(ir,ic,iland)
c
c          print * ,'ok612 of hillslope'

c****************************************************************
c variable using for water balance calculation
c****************************************************************
c
           groundoff2(isub,ihc)=groundoff2(isub,ihc)+ area(ir,ic)* 
     :      (d1_qsub-d1_ssub)*land_ratio(ir,ic,iland)/length(ir,ic)         ! m3/s                    !mqh,计算不同成分产流
           soiloff2(isub,ihc)=soiloff2(isub,ihc)+ area(ir,ic)* 
     :      d1_ssub*land_ratio(ir,ic,iland)/length(ir,ic)         ! m3/s                    !mqh,计算不同成分产流
	      groundoff(ir,ic,idc) = groundoff(ir,ic,idc) + 1000.0* 
     :                             dt*d1_qsub*land_ratio(ir,ic,iland)/
     :                             length(ir,ic)         ! mm
	      soiloff(ir,ic,idc)   = soiloff(ir,ic,idc)+1000.0*d1_ssub* 
     :		                     land_ratio(ir,ic,iland)
c            runoffd(ir,ic,idc)   = runoffd(ir,ic,idc) + 
c     :                             groundoff(ir,ic,idc)  ! mm
c
c****************************************************************
c
c
c****************************************************************
c          surface routing: steady constant sheet flow
c****************************************************************
c
	      soil_con_f = 1.0         ! soil conservation factor
	      
c
		  detension = Sstmax(iland)*soil_con_f
		  detension = detension * 2.0
		   detension = detension * ((dLAI / LAImax(iland))**0.5)
		
cc	修改于2006.4.16
cc     :                  * LAI(ir,ic,month)/LAImax(iland)
c     :                  * sqrt(LAI(ir,ic,month)/LAImax(iland))
c     :                  * sqrt(NDVI(ir,ic,month))
c            detension = amax1(3.0, detension)
c           detension = amax1(3.0, detension).


            detension = amax1(3.0, detension) 

		  water_depth = amax1(0.0, (Sst(ir,ic,iland)-detension) )   ! mm
c            water_depth = amax1(0.0, (Sst(ir,ic,iland)+tmppnet) )		  
c
		  q_hillslope = 0.0
c
            if( water_depth .le. 0.01 !) goto 1800					! 书签4
     :         .or. Drw(isub,iflow).ge. 1.0*Dr(isub,iflow) ) goto 1800

            Sst(ir,ic,iland) = Sst(ir,ic,iland)-water_depth
	      water_depth = 0.001 * water_depth          !in meter, surface runoff

            surface_n = surfn(iland) ! * sqrt(water_depth*1000.)
c     :                  * sqrt(LAI(ir,ic,month)/LAImax(iland))
cc
            waterhead = slp(ir,ic) ! + water_depth/length(ir,ic)
c	      if(waterhead .lt. 0.1E-8) waterhead = 0.1E-8
            power       = 1.6667

	      q_hillslope = dt * sqrt(waterhead)*
     :                    water_depth**power/surface_n     ! m3/m, one hillslope
		  if(q_hillslope .le. 0.1E-20) q_hillslope = 0.0

		  if(iland.eq.6)    ! for agricultural fields
     :         q_hillslope = q_hillslope * 
     :          amax1(0.3, (1.0-LAI(ir,ic,month,2)/LAImax(iland)))

            q_hillslope = amin1(q_hillslope, 
     :                          water_depth*length(ir,ic))

	      water_depth = water_depth - q_hillslope/length(ir,ic)	
			   
            qin_tmp= qin_tmp + q_hillslope/dt * 
     :                    land_ratio(ir,ic,iland)         !m3/s/m, one hillslope
c


        
c	      update surface storage
            Sst(ir,ic,iland) = Sst(ir,ic,iland)+1000.0*water_depth
            water_depth = 0.0
c
c          print * ,'ok613 of hillslope'

c****************************************************************
c variable using for water balance calculation
c****************************************************************
c
	      srunoff2(isub,ihc) = srunoff2(isub,ihc) + 
     :                           q_hillslope*land_ratio(ir,ic,iland)/
     :                           length(ir,ic)*area(ir,ic)/dt         ! m3/s,mqh
!	      srunoff(ir,ic,idc) = srunoff(ir,ic,idc) + 1000.0*
!     :                           q_hillslope*land_ratio(ir,ic,iland)/
!     :                           length(ir,ic)         ! mm
c            runoffd(ir,ic,idc)   = runoffd(ir,ic,idc) + 
c     :                             srunoff(ir,ic,idc)  ! mm
c
c****************************************************************
c
1800	    continue
   
		
		qin_tmp  = qin_tmp *area(ir,ic)/length(ir,ic)! m3/s, one flow-interval
	    qin_tmp0 = qin_tmp0 + qin_tmp
c
          runoffd(ir,ic,idc) = runoffd(ir,ic,idc) + 
     :                       1000.0*qin_tmp*dt/area(ir,ic)

c          print * ,'ok614 of hillslope'

	  end do
c

c          print * ,'ok61 of hillslope'

c	----------------------------------------------------------------------   
c       lateral inflow to the river
c
	  qin(iflow) = qin_tmp0 ! m3/s, total lateral inflow of one flow-interval

cc	修改于2006.1.15
cc	print *, qin(iflow)
c
c
cc	5555	删除于2006.4.14	 !skip irrigation loop

c
c****************************************************************
c    monthly average soil moisture of the root depth
c****************************************************************
c
        do 5444 ig = 1, ngrid(isub,iflow)
            ir    = grid_row(isub,iflow,ig)
            ic    = grid_col(isub,iflow,ig)
            isoil = soil(ir,ic)
          do 4888 iland = 1, nland
            if(land_ratio(ir,ic,iland) .le. 0.0) goto 4888
c  
            tmp = 0.0
            i   = 0
            do j = 1, layer(ir,ic)
              tmp = tmp + D(ir,ic,j)
              i = i+1
              if(tmp .ge. root(iland))     goto 4267  !原来是4266

            end do

cccccccccccccccccccccccccccccccccccccccccc            
            !added by mqh
4267                tmp19921 = 0.0
                tmp19922 = 0.0
                tmp19923 = 0.0
                tmp19924 = 0.0
                tmp19925 = 0.0
                tmp19926 = 0.0
       
            do j = 1, 1    !试试选基层土壤比较合适
              tmp19921 = tmp19921 + D(ir,ic,j)     
             
            end do
            do j = 1, 2    !试试选基层土壤比较合适
              tmp19922 = tmp19922 + D(ir,ic,j)     
             
            end do
            do j = 1, 3    !试试选基层土壤比较合适
              tmp19923 = tmp19923 + D(ir,ic,j)     
             
            end do
            do j = 1, 4    !试试选基层土壤比较合适
              tmp19924 = tmp19924 + D(ir,ic,j)     
             
            end do
            do j = 1, 5    !试试选基层土壤比较合适
              tmp19925 = tmp19925 + D(ir,ic,j)     
             
            end do
            do j = 1, 6    !试试选基层土壤比较合适
              tmp19926 = tmp19926 + D(ir,ic,j)     
             
            end do
            !end of add
cccccccccccccccccccccccccccccccccccccccccc            
4266        do j = 1, i
              soil_month(ir,ic) = soil_month(ir,ic) + 
     $                            land_ratio(ir,ic,iland)*
     $                            (w(ir,ic,iland,j)-wrsd(ir,ic))/
     $                            ((wsat(ir,ic)-wrsd(ir,ic))*24.0)*
     $                             D(ir,ic,j)/tmp
                  
            end do   
            
4888      continue
5444    continue  
c
c          print * ,'ok62 of hillslope'

c****************************************************************
c          
1999	continue


cc        print * ,'ok7 of hillslope'


c     save monthly spatial distribution
       
      if(isub.eq.end_sub .and. 
     :   day.eq.dayinmonth(month) .and. hour.eq.24) then
        write(ch2, '(i2.2)') month
        write(ch4, '(i4.4)') hydroyear
        do ii2 = 1, nsub
          do iflow = 1, nflow(ii2)
            do ig = 1, ngrid(ii2,iflow)
              ir    = grid_row(ii2,iflow,ig)
              ic    = grid_col(ii2,iflow,ig)
              isoil = soil(ir,ic)
		    if(area(ir,ic).lt.0.0) soil_month(ir,ic) = -9999.0
              if(area(ir,ic).lt.0.0) evap_month(ir,ic) = -9999.0
              if(area(ir,ic).lt.0.0) ecy_month(ir,ic)  = -9999.0
	        if(area(ir,ic).lt.0.0) ecp_month(ir,ic)  = -9999.0
	        if(area(ir,ic).lt.0.0) ese_month(ir,ic)  = -9999.0
              if(area(ir,ic).lt.0.0) runoff_month(ir,ic) = -9999.0
              if(area(ir,ic).lt.0.0) grunoff_month(ir,ic) = -9999.0
	        if(area(ir,ic).lt.0.0) drunoff_month(ir,ic) = -9999.0
	        if(area(ir,ic).lt.0.0) srunoff_month(ir,ic) = -9999.0
              soil_month(ir,ic) = soil_month(ir,ic)/
     $                            float( dayinmonth(month) )   ! monthly mean
              do i = idc-dayinmonth(month)+1, idc
                evap_month(ir,ic) = evap_month(ir,ic) + eactd(ir,ic,i)/
     $                              float( dayinmonth(month) ) ! monthly mean
                ecy_month(ir,ic) = ecy_month(ir,ic) + ecy(ir,ic,i)/
     $                              float( dayinmonth(month) ) ! monthly mean
	          ecp_month(ir,ic) = ecp_month(ir,ic) + ecp(ir,ic,i)/
     $                              float( dayinmonth(month) ) ! monthly mean
	          ese_month(ir,ic) = ese_month(ir,ic) + ese(ir,ic,i)/
     $                              float( dayinmonth(month) ) ! monthly mean
                runoff_month(ir,ic)= runoff_month(ir,ic) + 
     $               runoffd(ir,ic,i)/float(dayinmonth(month)) ! monthly mean
	          grunoff_month(ir,ic)= grunoff_month(ir,ic) + 
     $               groundoff(ir,ic,i)/float(dayinmonth(month)) ! monthly mean
	          srunoff_month(ir,ic)= srunoff_month(ir,ic) + 
     $               srunoff(ir,ic,i)/float(dayinmonth(month)) ! monthly mean
	          drunoff_month(ir,ic)= drunoff_month(ir,ic) + 
     $               soiloff(ir,ic,i)/float(dayinmonth(month)) ! monthly mean
              end do
            end do
          end do
        end do
c       print *,'evap.'//month//year
        call strlen(result2_dir, ia1,ia2)

 
        open(33,file = result2_dir(ia1:ia2)//  ! changed by xujijun
     $                 'evap_'//ch4//ch2//'.asc')        
	    write(33,'(A, i16)') ncols,nnc
	    write(33,'(A, i16)') nrows,nnr
	    write(33,'(A, f16.3)') xllcorner,x0
	    write(33,'(A, f16.3)') yllcorner,y0
	    write(33,'(A, f16.0)') cellsize,gridsize
	    write(33,'(A, f8.0)') nodata,znodata
        do ir=1 , nrow
	    do ic=1 , ncol
            if(area(ir,ic).lt.0.0) evap_month(ir,ic) = -9999.0
	    enddo
        enddo
	  do ir=1,nrow
		 write(33,'(295f13.4)') (evap_month(ir,ic), ic=1,ncol)
        end do
        close(33)

!cccccc	 open(33,file = result2_dir(ia1:ia2)//  ! changed by xujijun
!     $                 'ecy_'//ch4//ch2//'.asc')        
!	    write(33,'(A, i16)') ncols,nnc
!	    write(33,'(A, i16)') nrows,nnr
!	    write(33,'(A, f16.3)') xllcorner,x0
!	    write(33,'(A, f16.3)') yllcorner,y0
!	    write(33,'(A, f16.0)') cellsize,gridsize
!	    write(33,'(A, f8.0)') nodata,znodata
!        do ir=1 , nrow
!	    do ic=1 , ncol
!            if(area(ir,ic).lt.0.0) ecy_month(ir,ic) = -9999.0
!	    enddo
!        enddo
!	  do ir=1,nrow
!		 write(33,'(295f13.4)') (ecy_month(ir,ic), ic=1,ncol)
!        end do
!        close(33)

ccccc	open(33,file = result2_dir(ia1:ia2)//  ! changed by xujijun  !mqh7.18
c     $                 'ecp_'//ch4//ch2//'.asc')        
c	    write(33,'(A, i16)') ncols,nnc
c	    write(33,'(A, i16)') nrows,nnr
c	    write(33,'(A, f16.3)') xllcorner,x0
c	    write(33,'(A, f16.3)') yllcorner,y0
c	    write(33,'(A, f16.0)') cellsize,gridsize
c	    write(33,'(A, f8.0)') nodata,znodata
c       do ir=1 , nrow
c	    do ic=1 , ncol
c            if(area(ir,ic).lt.0.0) ecp_month(ir,ic) = -9999.0
c	    enddo
c       enddo
c	  do ir=1,nrow
c		 write(33,'(295f13.4)') (ecp_month(ir,ic), ic=1,ncol)
c       end do
c      close(33)

cccc	open(33,file = result2_dir(ia1:ia2)//  ! changed by xujijun    !mqh7.18
c     $                 'ese_'//ch4//ch2//'.asc')        
c	    write(33,'(A, i16)') ncols,nnc
c	    write(33,'(A, i16)') nrows,nnr
c	    write(33,'(A, f16.3)') xllcorner,x0
c	    write(33,'(A, f16.3)') yllcorner,y0
c	    write(33,'(A, f16.0)') cellsize,gridsize
c	    write(33,'(A, f8.0)') nodata,znodata
c       do ir=1 , nrow
c	    do ic=1 , ncol
c            if(area(ir,ic).lt.0.0) ese_month(ir,ic) = -9999.0
c	    enddo
c       enddo
c	  do ir=1,nrow
c		 write(33,'(295f13.4)') (ese_month(ir,ic), ic=1,ncol)
c       end do
c        close(33)

ccccc      open(33,file = result2_dir(ia1:ia2)//  ! changed by xujijun       !mqh7.18
c     $                 'snow_'//ch4//ch2//'.asc')        
c	    write(33,'(A, i16)') ncols,nnc
c	    write(33,'(A, i16)') nrows,nnr
c	    write(33,'(A, f16.3)') xllcorner,x0
c	    write(33,'(A, f16.3)') yllcorner,y0
c	    write(33,'(A, f16.0)') cellsize,gridsize
c	    write(33,'(A, f8.0)') nodata,znodata
c       do ir=1 , nrow
c	    do ic=1 , ncol
c	      snow1(ir,ic) = 0.0
c		  do i = 1, nland
c	     snow1(ir,ic) = snow1(ir,ic)+snow(ir,ic,i)*land_ratio(ir,ic,i)
c	      enddo
c            if(area(ir,ic).lt.0.0) snow1(ir,ic) = -9999.0
c	    enddo
c        enddo
c	  do ir=1,nrow
c		 write(33,'(295f13.4)') (snow1(ir,ic), ic=1,ncol)
c        end do
cccccc        close(33)

c       print *,'soil.'//month//year
c        open(33,file=result2_dir(ia1:ia2)//
c     $       'soil.'//ch2//ch4,status='replace')

cccccc        open(33,file = result2_dir(ia1:ia2)//  ! changed by xujijun   !mqh7.18
c     $                 'runoff_'//ch4//ch2//'.asc')        
c	    write(33,'(A, i16)') ncols,nnc
c	    write(33,'(A, i16)') nrows,nnr
c	    write(33,'(A, f16.3)') xllcorner,x0
c	    write(33,'(A, f16.3)') yllcorner,y0
c	    write(33,'(A, f16.0)') cellsize,gridsize
c	    write(33,'(A, f8.0)') nodata,znodata
c       do ir=1 , nrow
c	    do ic=1 , ncol
c            if(area(ir,ic).lt.0.0) runoff_month(ir,ic) = -9999.0
c	    enddo
c        enddo
c	  do ir=1,nrow
c		 write(33,'(10f13.4)') (runoff_month(ir,ic), ic=1,ncol)
c        end do
c        close(33)

cccc	open(33,file = result2_dir(ia1:ia2)//  ! changed by xujijun   !mqh7.18
c     $                 'grunoff_'//ch4//ch2//'.asc')        
c	    write(33,'(A, i16)') ncols,nnc
c	    write(33,'(A, i16)') nrows,nnr
c	    write(33,'(A, f16.3)') xllcorner,x0
c	    write(33,'(A, f16.3)') yllcorner,y0
c	    write(33,'(A, f16.0)') cellsize,gridsize
c	    write(33,'(A, f8.0)') nodata,znodata
c       do ir=1 , nrow
c	    do ic=1 , ncol
c            if(area(ir,ic).lt.0.0) grunoff_month(ir,ic) = -9999.0
c	    enddo
c        enddo
c	  do ir=1,nrow
c		 write(33,'(10f13.4)') (grunoff_month(ir,ic), ic=1,ncol)
c        end do
c        close(33)

cccc	open(33,file = result2_dir(ia1:ia2)//  ! changed by xujijun   !mqh7.18
c     $                 'srunoff_'//ch4//ch2//'.asc')        
c	    write(33,'(A, i16)') ncols,nnc
c	    write(33,'(A, i16)') nrows,nnr
c	    write(33,'(A, f16.3)') xllcorner,x0
c	    write(33,'(A, f16.3)') yllcorner,y0
c	    write(33,'(A, f16.0)') cellsize,gridsize
c	    write(33,'(A, f8.0)') nodata,znodata
c       do ir=1 , nrow
c	    do ic=1 , ncol
c           if(area(ir,ic).lt.0.0) srunoff_month(ir,ic) = -9999.0
c	    enddo
c       enddo
c	  do ir=1,nrow
c		 write(33,'(10f13.4)') (srunoff_month(ir,ic), ic=1,ncol)
c        end do
c        close(33)

ccc	open(33,file = result2_dir(ia1:ia2)//  ! changed by xujijun   !mqh
c     $                 'drunoff_'//ch4//ch2//'.asc')        
c	    write(33,'(A, i16)') ncols,nnc
c	    write(33,'(A, i16)') nrows,nnr
c	    write(33,'(A, f16.3)') xllcorner,x0
c	    write(33,'(A, f16.3)') yllcorner,y0
c	    write(33,'(A, f16.0)') cellsize,gridsize
c	    write(33,'(A, f8.0)') nodata,znodata
c       do ir=1 , nrow
c	    do ic=1 , ncol
c            if(area(ir,ic).lt.0.0) drunoff_month(ir,ic) = -9999.0
c	    enddo
c        enddo
c	  do ir=1,nrow
c		 write(33,'(10f13.4)') (drunoff_month(ir,ic), ic=1,ncol)
c       end do
c        close(33)



!        open(33,file=result2_dir(ia1:ia2)//    ! changed by xujijun
!     $       'soil_'//ch4//ch2//'.asc')
!	    write(33,'(A, i16)') ncols,nnc
!	    write(33,'(A, i16)') nrows,nnr
!	    write(33,'(A, f16.3)') xllcorner,x0
!	    write(33,'(A, f16.3)') yllcorner,y0
!	    write(33,'(A, f16.0)') cellsize,gridsize
!	    write(33,'(A, f8.0)') nodata,znodata
!        do ir=1 , nrow
!	    do ic=1 , ncol
!            if(area(ir,ic).lt.0.0) soil_month(ir,ic) = -9999.0
!	    enddo
!        enddo
!  
!        do ir=1,nrow
!          write(33,'(295f13.4)') (soil_month(ir,ic), ic=1,ncol)
!        end do
!        close(33)
!c     
!        do ir = 1, nrow
!          do ic = 1, ncol
!            evap_month(ir,ic) = 0.0
!            soil_month(ir,ic) = 0.0
!            runoff_month(ir,ic) = 0.0
!	      ecy_month(ir,ic) = 0.0
!            ecp_month(ir,ic) = 0.0
!	      ese_month(ir,ic) = 0.0
!	      drunoff_month(ir,ic) = 0.0
!	      grunoff_month(ir,ic) = 0.0
!	      srunoff_month(ir,ic) = 0.0
!          end do
!        end do
      end if
      

  
c        endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!输出不同成分产流机制的大小，单位：mm
!      if (isub.eq.end_sub) then
!        srunoff_total=0
!        soiloff_total=0
!        groundoff_total=0
!        do jj=1,end_sub
!            srunoff_total=srunoff2(jj,ihc)+srunoff_total
!            soiloff_total= soiloff_total+ soiloff2(jj,ihc)
!            groundoff_total=groundoff_total+groundoff2(jj,ihc)
!        enddo
!
!!        print *,srunoff_total,soiloff_total,groundoff_total
!!        pause
!          !!groundoff          
!
!      if (hydroyear.eq.2000.and.ihc.eq.1) then
!       open(69,file=result2_dir(ia1:ia2)//'su_so_g_runoff.dat',
!     $status='unknown')
!		elseif((hydroyear.gt.2000).or.(hydroyear.eq.2000.and.ihc.gt.1)) then
!	open(69,file=result2_dir(ia1:ia2)//'su_so_g_runoff.dat',
!     $access='append',status='old')
!        endif
!       write(69, '(4i6,f16.3)') hydroyear,month,day,hour,srunoff_total
!     :soiloff_total,groundoff_total
!!       print *,'groundoff2'
!       close(69)         
!      endif


      
c 1010 continue
c
c      output monthly irrigation water intake for each sub-catchment
c      if(idc.eq.dayinmonth(month) .and. hour.eq.24) then
c        write(ch2, '(i4.4)') month
c        write(ch4, '(i4.4)') hydroyear
c        call strlen(result1_dir, ia1,ia2)
c        open(10,file=result1_dir(ia1:ia2)//
c     $       subbasin(isub)//'_irr.'//ch4, status='replace')
c        j2=0
c        do imm=1,12
c            tmp1=0.0
c            tmp2=0.0
c            do j1=1,dayinmonth(imm)
c              j2=j2+1
c              tmp1=tmp1+total_intake(isub,j2)
c              tmp2=tmp2+total_req(isub,j2)
c            end do
c            write(10,*) imm, tmp1,tmp2
c        end do
c        close(10)
c      endif
c     
c     output hydrological characteristics at end of simulation year
            
19921    if(month.eq.endmonth .and. day.eq.endday .and. hour.eq.24) then
c 
    
c   (1) sub-basin mean: rain, runoff, Eact
        write(ch4, '(i4.4)') hydroyear     
        call strlen(simul_dir, ia1,ia2)
c        open(10, file = simul_dir(ia1:ia2)//
c     $                  subbasin(isub)//'_runoff.'//ch4, 
c     $                  status='replace')

        open(10, file = simul_dir(ia1:ia2)//         !  changged by xujijun
     $                  subbasin(isub)//'_runoff'//ch4)

        do jj = 1, idc
          tmp1 = 0.0
          tmp2 = 0.0
          tmp3 = 0.0
          do iflow = 1, nflow(isub)
	      do ig = 1, ngrid(isub,iflow)
	        ir = grid_row(isub,iflow,ig)
	        ic = grid_col(isub,iflow,ig)
	        if(area(ir,ic) .gt. 0.0) then
                tmp1 = tmp1 + raind(ir,ic,jj)*area(ir,ic)
                tmp2 = tmp2 + runoffd(ir,ic,jj)*area(ir,ic)
                tmp3 = tmp3 + eactd(ir,ic,jj)*area(ir,ic)
              endif
	      end do
	    end do
          tmp1 = tmp1 / subarea(isub)
          tmp2 = tmp2 / subarea(isub)
          tmp3 = tmp3 / subarea(isub)
          write(10,70,err=75) jj, tmp1, tmp2,tmp3
 70       format(i6, 3f16.5)
        end do
75      close(10)
          
c     
c       (2) whole basin mean: rain, runoff, Eact
c        if(isub.eq.end_sub) then
c          open(12, file = simul_dir(ia1:ia2)//
c     $                    'rain-evap-runoff.'//ch4,status='replace')

        if(isub.eq.end_sub) then
          open(12, file = simul_dir(ia1:ia2)//       ! changed  by xujijun  
     $                    'rain-evap-runoff'//ch4)

        do jj = 1, idc
          tmp1 = 0.0
          tmp2 = 0.0
          tmp3 = 0.0
          do ii2 = 1, nsub
            do iflow = 1, nflow(ii2)
	        do ig = 1, ngrid(ii2,iflow)
	          ir = grid_row(ii2,iflow,ig)
	          ic = grid_col(ii2,iflow,ig)
	          if(area(ir,ic) .gt. 0.0) then
                  tmp1 = tmp1 + raind(ir,ic,jj)*area(ir,ic)
                  tmp2 = tmp2 + runoffd(ir,ic,jj)*area(ir,ic)
                  tmp3 = tmp3 + eactd(ir,ic,jj)*area(ir,ic)
                endif
	        end do
	      end do
          end do
          tmp1 = tmp1 / basinarea
          tmp2 = tmp2 / basinarea
          tmp3 = tmp3 / basinarea
          write(12,130) i,tmp1, tmp2, tmp3
 130      format(2x,i6, 3f16.3)
        end do
        close(12)
c         added by gaobing
c              open(14, file = simul_dir(ia1:ia2)//       ! changed  by xujijun  
c     $                    'rain-yingluoxia'//ch4//'.txt')

c        do jj = 1, idc
c          tmp1 = 0.0
c          tmp2 = 0.0
c          do ii2 = 1, 153
c            do iflow = 1, nflow(ii2)
c	        do ig = 1, ngrid(ii2,iflow)
c	          ir = grid_row(ii2,iflow,ig)
c	          ic = grid_col(ii2,iflow,ig)
c	          if(area(ir,ic) .gt. 0.0) then
c                  tmp1 = tmp1 + raind(ir,ic,jj)*area(ir,ic)                  
c               endif
c	        end do
c	      end do
c	            tmp2 = tmp2 + subarea(ii2)
c          end do
c          tmp1 = tmp1 / tmp2
c        
c          write(14,131) i,tmp1, tmp2
c131      format(2x,i6, 2f16.3) 
c	    end do
c           close(14)


c	open(14, file = simul_dir(ia1:ia2)//       ! changed  by xujijun  
c     $                    'rain-zhamashike'//ch4//'.txt')

c        do jj = 1, idc
c          tmp1 = 0.0
c          tmp2 = 0.0
c          do ii2 = 1, 90
c            do iflow = 1, nflow(ii2)
c	        do ig = 1, ngrid(ii2,iflow)
c	          ir = grid_row(ii2,iflow,ig)
c	          ic = grid_col(ii2,iflow,ig)
c	          if(area(ir,ic) .gt. 0.0) then
c                  tmp1 = tmp1 + raind(ir,ic,jj)*area(ir,ic)                  
c                endif
c	        end do
c	      end do
c	            tmp2 = tmp2 + subarea(ii2)
c          end do
c          tmp1 = tmp1 / tmp2
        
c          write(14,131) i,tmp1, tmp2

c	    end do
c           close(14)
c         end of added

        endif
      end if
c
c       if(isub==124 .and. day == 31 ) print *, 'ok hillslopeend'
c	if(day .eq. dayinmonth(month) .and. hour.eq.24) then
c	if(month.eq.endmonth .and. day.eq.endday .and. hour.eq.24) then
	  if(isub.eq.start_sub) then
	    annual_rain  = 0.0
cc	    annual_irrigation=0.0
		annual_runoff= 0.0
	    annual_Eact  = 0.0
	    annual_Cst   = 0.0
	    annual_Ssts  = 0.0
	    annual_Sstr  = 0.0
	    annual_SBst  = 0.0
	    annual_Gst   = 0.0
	    annual_Ecy   = 0.0
	    annual_Ecp   = 0.0
	    annual_Ese   = 0.0
cc	    annual_re    = 0.0
        endif
  	 
	  do iflow = 1, nflow(isub)
c	     annual_re = annual_re + re_capacity(isub,iflow)
		do ig=1,ngrid(isub,iflow)
	      ir=grid_row(isub,iflow,ig)
	      ic=grid_col(isub,iflow,ig)
	      isoil=soil(ir,ic)
            do i = 1, idc         !366
              annual_runoff = annual_runoff + runoffd(ir,ic,i)*
     :                        area(ir,ic)                           ! mm
              annual_Eact   = annual_Eact + eactd(ir,ic,i)*
     :                        area(ir,ic)     		               ! mm
cc              annual_Ecy   = annual_Ecy + ecy(ir,ic,i)*
cc     :                        area(ir,ic)     		               ! mm
cc              annual_Ecp   = annual_Ecp + ecp(ir,ic,i)*
cc     :                        area(ir,ic)     		               ! mm
cc              annual_Ese   = annual_Ese + ese(ir,ic,i)*
cc     :                        area(ir,ic)     		               ! mm
             annual_rain   = annual_rain + raind(ir,ic,i)*           !rain_daily
     :                        area(ir,ic)              		       ! mm
c              annual_irrigation= annual_irrigation +
c     :                        total_irr(ir,ic,i)*area(ir,ic)        ! mm
            end do
            
            do 2240 iland = 1, nland
              if(land_ratio(ir,ic,iland) .le. 0.0) goto 2240
              annual_Cst = annual_Cst + Cst(ir,ic,iland)*
     :                    0.001*land_ratio(ir,ic,iland)*area(ir,ic)    !in m3
              annual_Sstr = annual_Sstr + Sst(ir,ic,iland)*
     :                    0.001*land_ratio(ir,ic,iland)*area(ir,ic)    !in m3
              annual_Ssts = annual_Ssts + snow(ir,ic,iland)*
     :                    0.001*land_ratio(ir,ic,iland)*area(ir,ic)    !in m3
              annual_Gst = annual_Gst+GWst(ir,ic,iland)*area(ir,ic)
     :                           	*land_ratio(ir,ic,iland)             !in m3
              do j = 1 ,layer(ir,ic)
			  annual_SBst = annual_SBst + w(ir,ic,iland,j) *
     :              D(ir,ic,j)*land_ratio(ir,ic,iland)*area(ir,ic)     !in m3
              end do
2240        continue
          end do
	  end do

	  if(isub .eq. end_sub) then
c	    print *,"annual_SBst=", annual_SBst		  
	    annual_rain  = annual_rain / basinarea
cc	    annual_irrigation=annual_irrigation/basinarea
	    annual_runoff= annual_runoff / basinarea
	    annual_Eact  = annual_Eact / basinarea
	    annual_Ecy  = annual_Ecy / basinarea
	    annual_Ecp  = annual_Ecp / basinarea
	    annual_Ese  = annual_Ese / basinarea
	    annual_Cst  = 1000.0*(annual_Cst-annual_Cst0)/basinarea
	    annual_Ssts = 1000.0*(annual_Ssts-annual_Ssts0)/ basinarea
	    annual_Sstr = 1000.0*(annual_Sstr-annual_Sstr0)/ basinarea
	    annual_SBst = 1000.0*(annual_SBst-annual_SBst0)/ basinarea
	    annual_Gst  = 1000.0*(annual_Gst-annual_Gst0)/basinarea
cc	    annual_re   = 1000.0*(annual_re-annual_re0)/basinarea
	    unbalance   = annual_rain - annual_runoff - annual_Eact - 
     :                   annual_Cst - annual_Ssts- annual_Sstr- 
     :                   annual_SBst -annual_Gst !+ annual_irrigation

	    call strlen(result2_dir,ia1,ia2)
          write(ch4, '(i4.4)') hydroyear
          open(35,file = result2_dir(ia1:ia2)//'waterbalance'//ch4, 
     :                                            status='replace')
          write(35,1110) hydroyear
1110      format(//25x,'Water Balance ',i4, ' (unit: mm)'/)
          write(35,*) '   Rainfall  S-storage(snow)    S-storage(r)',
     :     '    Sub-storage    G-storage     Evaporation     Runoff',
     :     '     error        efromconpy    efromcorp   efromsurface',
     :     '     reservoir     irrigation'
c       由mqh注释掉
c          print *, annual_rain,annual_Ssts,annual_Sstr,annual_SBst,
c    :          annual_Gst,annual_Eact, annual_runoff, unbalance

          write(35,1120) annual_rain, annual_Ssts,annual_Sstr,
     :            annual_SBst,annual_Gst,annual_Eact, annual_runoff
     :            , unbalance,annual_Ecy,annual_Ecp,annual_Ese
c     :            ,annual_re,annual_irrigation
1120      format(4x,f16.3,11(5x,f16.3))
	    close(35)
 
          write(ch2, '(i2.2)') month   ! changed by xujijun
 
	    call strlen(result2_dir,ia1,ia2)
          open(36,file = result2_dir(ia1:ia2)//'waterbalance.'//ch2, 
     :                                              status='replace')
          write(36,1110) month
          write(36,*) '   Rainfall  S-storage(snow)    S-storage(r)',
     :     '    Sub-storage    G-storage     Evaporation     Runoff',
     :     '     error        efromconpy    efromcorp   efromsurface',
     :     '     reservoir    irrigation'
          write(36,1120) annual_rain, annual_Ssts,annual_Sstr,
     :            annual_SBst,annual_Gst,annual_Eact, annual_runoff
     :            , unbalance,annual_Ecy,annual_Ecp,annual_Ese
     :           !  ,annual_re,annual_irrigation
	    close(36)
	  endif
c      endif
c
c       (4a) save the status at end of the simulation
c

	if(month.eq.endmonth.and.day.eq.endday.and.hour.eq.24) then
        call strlen(simul_dir,ia1,ia2)
c	  open(9,file = simul_dir(ia1:ia2)//
c     :                subbasin(isub)//'I_soil',status='replace')

	  open(9,file = simul_dir(ia1:ia2)//      ! changed by xujijun
     :                subbasin(isub)//'I_soil2')

c	  print *, simul_dir(ia1:ia2)//subbasin(isub)//'I_soil'
        do iflow = 1, nflow(isub)
	    do ig=1,ngrid(isub,iflow)
	      ir=grid_row(isub,iflow,ig)
	      ic=grid_col(isub,iflow,ig)
	      isoil=soil(ir,ic)
            write(9,*) iflow, soil(ir,ic), nland, layer(ir,ic)
            do iland = 1, nland
              if(land_ratio(ir,ic,iland) .gt. 0.0) 
     :        write(9,*) (w(ir,ic,iland,j),j=1,layer(ir,ic)), 
     :                    Dgl(ir,ic,iland), Gwst(ir,ic,iland)
            end do
          end do
        end do
        close(9)
	end if
c
c       (4b) save the status at start of the simulation in every year
c

	if(month.eq.startmonth.and.day.eq.startday.and.hour.eq.1) then
        write(ch4,'(i4.4)') year
	  call strlen(simul_dir,ia1,ia2)
	  open(9,file = simul_dir(ia1:ia2)//
     :                subbasin(isub)//'I_soil2'//ch4,status='unknown')
c	  print *, simul_dir(ia1:ia2)//subbasin(isub)//'I_soil'
        do iflow = 1, nflow(isub)
	    do ig=1,ngrid(isub,iflow)
	      ir=grid_row(isub,iflow,ig)
	      ic=grid_col(isub,iflow,ig)
	      isoil=soil(ir,ic)
            write(9,*) iflow, soil(ir,ic), nland, layer(ir,ic)
            do iland = 1, nland
              if(land_ratio(ir,ic,iland) .gt. 0.0) 
     :        write(9,*) (w(ir,ic,iland,j),j=1,layer(ir,ic)), 
     :                    Dgl(ir,ic,iland), Gwst(ir,ic,iland)
            end do
          end do
        end do
        close(9)
	end if

cc     
c       (5) save the irrigation data
c



c       (6) save the reservoir data
c     

c     
c       (7) save the irrigation soil moisture data

!        deallocate (snow1)

      return
	end
c
c
c****************************************************************
c
c	Sub-program -- downscaling forcing data 
c                        daily -> hourly
c
c****************************************************************
c
      subroutine rain_model(isub, ir, ic, idh, month,           
     :           rain_daily, tmin_daily, tmax_daily, evap_daily,
     :           r, T, Ep, iland)

	include "dims2.inc"

	integer isub, ir, ic, iland, month
      integer idh, tmd
      integer rt0(nrow,ncol)
      real    rain_daily, tmin_daily, tmax_daily, evap_daily
	real    r, T, Ep         
c
c------------------------ create the start time ------------
c
	if(idh .eq. 1 .and. iland.eq.1) then
	  tmp = ran0(idum)
	  rt0(ir,ic) = 4 + int((21-4)*tmp+0.45)
	  if(rt0(ir,ic) .gt. 21) rt0(ir,ic) = 21
c	  print *,"rt0(ir,ic)=", rt0(ir,ic), tmp
	endif
c
c------------------------ downscaling Ep -------------------
c
	  if(evap_daily.lt.0.0) evap_daily = 0.0
	  Ep = 0.0
	!------------ 引用自然日入流数据 ----------------!
	  if(idh.eq.6.or.idh.eq.17) Ep = evap_daily*0.03
	  if(idh.eq.7.or.idh.eq.16) Ep = evap_daily*0.06
	  if(idh.eq.8.or.idh.eq.15) Ep = evap_daily*0.08
	  if(idh.eq.9.or.idh.eq.14) Ep = evap_daily*0.10
	  if(idh.eq.10.or.idh.eq.13) Ep = evap_daily*0.11
	  if(idh.eq.11.or.idh.eq.12) Ep = evap_daily*0.12
c
c
c------------------------ downscaling Prec --------------------
c
cc      un_understand 

cc      changed by xujijun

cc	if(isub.ge.22 .and. month .le. 5) then
cc	  r = 0.0
cc	  if(rain_daily.ge.0.0 .and. rain_daily.le. 8.0) then
cc          if(idh.eq.rt0(ir,ic))     r = rain_daily
cc        elseif(rain_daily.gt.8.0 .and. rain_daily.le.15.0) then
cc          if(idh.eq.(rt0(ir,ic)-1)) r = 0.10*rain_daily
cc          if(idh.eq.rt0(ir,ic))     r = 0.80*rain_daily
cc          if(idh.eq.(rt0(ir,ic)+1)) r = 0.10*rain_daily
cc        elseif(rain_daily.gt.15.0 .and. rain_daily.le.25.0) then
cc          if(idh.eq.(rt0(ir,ic)-2)) r = 0.10*rain_daily
cc          if(idh.eq.(rt0(ir,ic)-1)) r = 0.20*rain_daily
cc          if(idh.eq.rt0(ir,ic))     r = 0.40*rain_daily
cc          if(idh.eq.(rt0(ir,ic)+1)) r = 0.20*rain_daily
cc          if(idh.eq.(rt0(ir,ic)+2)) r = 0.10*rain_daily
cc        elseif(rain_daily.gt.25.0 .and. rain_daily.le.50.0) then
cc          if(idh.eq.(rt0(ir,ic)-4)) r = 0.05*rain_daily
cc          if(idh.eq.(rt0(ir,ic)-3)) r = 0.10*rain_daily
cc          if(idh.eq.(rt0(ir,ic)-2)) r = 0.10*rain_daily
cc          if(idh.eq.(rt0(ir,ic)-1)) r = 0.15*rain_daily
cc          if(idh.eq.rt0(ir,ic))     r = 0.20*rain_daily
cc          if(idh.eq.(rt0(ir,ic)+1)) r = 0.15*rain_daily
cc          if(idh.eq.(rt0(ir,ic)+2)) r = 0.10*rain_daily
cc          if(idh.eq.(rt0(ir,ic)+3)) r = 0.10*rain_daily
cc	    if(idh.eq.(rt0(ir,ic)+4)) r = 0.05*rain_daily
cc       elseif(rain_daily.gt.50.0) then
cc		r = rain_daily/24.0
cc	  end if
cc	else
	  r = 0.0
	  if(rain_daily.ge.0.0 .and. rain_daily.le. 5.0) then
          if(idh.eq.rt0(ir,ic))     r = rain_daily
        elseif(rain_daily.gt.5.0 .and. rain_daily.le.10.0) then
          if(idh.eq.(rt0(ir,ic)-1)) r = 0.30*rain_daily
          if(idh.eq.rt0(ir,ic))     r = 0.40*rain_daily
          if(idh.eq.(rt0(ir,ic)+1)) r = 0.30*rain_daily
        elseif(rain_daily.gt.10.0 .and. rain_daily.le.30.0) then
          if(idh.eq.(rt0(ir,ic)-3)) r = 0.10*rain_daily
          if(idh.eq.(rt0(ir,ic)-2)) r = 0.15*rain_daily
          if(idh.eq.(rt0(ir,ic)-1)) r = 0.15*rain_daily
          if(idh.eq.rt0(ir,ic))     r = 0.20*rain_daily
          if(idh.eq.(rt0(ir,ic)+1)) r = 0.15*rain_daily
          if(idh.eq.(rt0(ir,ic)+2)) r = 0.15*rain_daily
          if(idh.eq.(rt0(ir,ic)+3)) r = 0.10*rain_daily
        elseif(rain_daily.gt.30.0) then
		r = rain_daily/24.0
c          if(idh.le.3  .or. idh.ge.22) r = 0.5 * rain_daily/24.0
c          if(idh.ge.11 .and. idh.le.13) r = 2.0 * rain_daily/24.0
	  end if
cc      endif
c
c
c------------------------ downscaling Temp --------------------
c
	  if(tmin_daily.eq.-9999) then
	  	tmin_daily = 0.0
	  	print *,"min temperature is  -9999"
	  end if
	  if(tmax_daily.eq.-9999) then
	  	tmax_daily = 0.0
	  	print *,"max temperature is  -9999"
	  end if

c	  if(idh.ge.1 .and. idh.lt.12) T=tmin_daily+
c     $	         float(idh)*(tmax_daily-tmin_daily)/12.0
c	  if(idh.ge.12 .and. idh.lt.24) T=tmax_daily+
c     $	         float(idh-12)*(tmin_daily-tmax_daily)/12.0
	  if(idh.ge.4 .and. idh.lt.14) then
	     T = tmin_daily+ (tmax_daily-tmin_daily)*
     $	    (1-cos((idh-4)*3.1415926/10.0))/2.0
	   else
	     tmd = idh
		 if(idh.lt.4) tmd = 24+idh
		 T = tmin_daily + (tmax_daily-tmin_daily)*
     $	    (1+cos((tmd-14)*3.1415926/14.0))/2.0
	   end if
      return
	end
c
c    
      FUNCTION ran0(idum)
      INTEGER idum,IA,IM,IQ,IR,MASK
      REAL ran0,AM
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *MASK=123459876)
      INTEGER k
      idum=ieor(idum,MASK)
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      ran0=AM*idum
      idum=ieor(idum,MASK)
      return
      END
c
c
c
c****************************************************************
      subroutine ETsoil(cm,cmf,cmw,Kcm)
	real cm,cmf,cmw,Kcm
	real cmc1,cmc2

c	if(cm.lt.cmw) cm = cmw
c      Kcm = (cm-cmw)/(cmf-cmw)
c      if(cm .ge. cmf)  Kcm=1.0

	c=0.5
	r=0.5
	Kcm=0.0
	cmc2=r*(cmf+cmw)
	cmc1=(cmc2-c*cmf)/(1.0-c)
	if(cm .ge. cmf)  Kcm=1.0
	if(cm.lt.cmf .and. cm.ge.cmc2) Kcm=(cm-cmc1)/(cmf-cmc1)
	if(cm.lt.cmc2 .and. cm.gt.cmw) Kcm=c*(cm-cmw)/(cmc2-cmw)

c	r=0.55
c	Kcm=0.0
c	cmc2=r*(cmf+cmw)
c	cmc1=(cm-cmc2)/(cmf-cmc2)
c
c	if(cm .ge. cmf)  Kcm=1.0
c	if(cm.lt.cmf .and. cm.gt.cmw) then
c	  Kcm = 1.0 /(1.0 + exp(-5.0*cmc1) )
c      end if

	if(Kcm.lt.0.0 .or. Kcm.gt.1.0) 
     $	print *,'wrong in calculating Kcm in ETsoil',cm,cmw,cmf,Kcm
	end
c
c
c
c****************************************************************
c	subroutine rexcess(re,r,k0)
c	real k0,re,r
c	if((k0/r).le.5.0) then
c	  re=r*exp(-1.0*k0/r)
c	  re=amin1(re,r)
c	 else
c	  re=0.0
c	end if
c
c	return
c	end
c
c
c
c**********************************************************************
c                                                                     *
c                   ***** SUBROUTINE RUNOFF  *****                    *
c                                                                     *
c**********************************************************************	
c
      SUBROUTINE runoff(isub, month, snow, temper, 
     :                   land_ratio,slope,length, Ds, Dg, Dr, anik, 
     :                   wsat, wrsd,wfld, alpha, watern,nuz, deltz, 
     :               k0,w, kg, GWcs, GWst, Dgl, Drw, Sst,dt, qsub,ssub)
c
c**********************************************************************	
c
c     PURPOSE:
c    calculation of soil moisture, infiltration excess, recharge to                  
c    groundwater and subsurface flow                
c                                                                               
C#######################################################################
c      
c     AUTHOR: Dawen Yang
c     8/1998
c
c     MODIFICATION HISTORY: 10/2001
c
c
c#######################################################################
c
c     Variable Declarations.
c
c#######################################################################
c
      implicit none
	include 'dims2.inc'
c
c#######################################################################
c
c     hillslope, groundwater, UZ,  variables
c
c#######################################################################
c
	integer isub
	integer month          ! month in a year
	real    temper         ! air temperature (degree)
	real    snow           ! snow depth (mm-water)
      real    land_ratio     ! area ratio of the land-soil type on a hillslope
      real    slope	       ! slope gradient of hillslope (ND)
	real    length		   ! length of hillslope (m)
	real    Ds			   ! depth of topsoil (m)
	real    Dg			   ! depth of unconfined acquifer (m)
	real    Dr			   ! depth of river (m)

	real    wsat           ! saturated soil moisture
	real    wrsd		   ! residual soil moisture
	real    wfld		   ! soil moisture at field capacity
	real    alpha		   ! soil-water parameter
	real    watern		   ! soil-water parameter
	real    anik           ! soil anitropic ratio (>=1.0)

	integer nuz            ! number of UZ layers
	real    deltz(nlayer)  ! depth of each layer (m)
	real    k0(nlayer)	! saturated hydraulic conductivity of each layer (m/s)
	real    w(nlayer)	   ! soil moisture of each layer 
	real    suct(nlayer)   ! suction of each layer (meter H2O)
	real    ksoil(nlayer)  ! hydraulic conductivity of each layer (m/s)
	real    avk(nlayer)	   ! average hydraulic conductivity
	real    avz(nlayer)    ! average depth 
	real    dpdw(nlayer)   ! 

	real    GWcs           ! GW storage coefficient
	real    kg			   ! GW hydraulic conductivity (m/s)
	real    avkg		   ! average GW hydraulic conductivity (m/s)

	real    GWst	       ! Groundwater storage (m) 
	real    Dgl            ! depth to GW level (m)
	real    Drw			   ! water depth in river (m)
	real    qsub		   ! runoff fromsubsurface (m3/s/m)
	real    Sst			   ! surface storage (m)
      real    ssub           ! runoff from soil layer ! added by gaobing
	real    dt			   ! time step (second)
	real    ss_f           ! slope shape (concave or convex) factor
c
c#######################################################################
c
c     temporary variables
c
c#######################################################################
c
      integer i, j, isat, GWkey
	real  tmp_D, tmpqsub,ps1, ps2, rech        !qsub1,qsub2,dtmp, , tw
      real  conductivity_V,excess,deficit,qmin,qmax,cmtmp, tmp,tmp2	 !
	integer MARK
	real    temA(nuz), temB(nuz), temC(nuz), temD(nuz), temq(nuz)
	real  inf, f    ! inf for infiltration? ,mqh
	  
	real total_rech ,total_qsub,total_ssub    !added by mqh ,2014.4.24
	save total_rech,total_qsub,total_ssub
	integer year,month2, day, hour, idc,ihc
	common /date1  / year, month2, day, hour, idc,ihc
  

cc      if(isub.eq.1.and.ihc.eq.1)   total_qsub=0  !mqh
cc      if(isub.eq.1.and.ihc.eq.1)   total_rech=0
cc       if(isub.eq.1.and.ihc.eq.1)   total_ssub=0
c

c      print *, 'nuz=', nuz 
      do i=1,nuz
        temq(i) = 0.0
      end do
c
c------------------------ concave hillslope -------------------------

cc     changed by gaobing

	  ss_f = 0.02    ! linear or convex slope   2013/4/13 23:05
	  
c
c-----------------------Soil property for Loess----------------------
c
c	if(isub.ge.36 .and. isub.le.137) anik = 1.0
c
c--------------------------------------------------------------------
c      GWcs = 0.15	  ! GW storage coefficient
c
C#######################################################################
c
c    calculation of recharge rate to groundwater: rech                                                    
c 
C#######################################################################
c
	inf  = 0.0
      rech = 0.0
	tmp  = 0.0

c      print *,'dgl,ds,dg',dgl,ds,dg
c      pause
      if(Dgl .ge. Ds) then
        call SuctionFromMoisture_V(w(nuz),wsat,wrsd,watern,alpha,ps1)  ! m
        !采用Van Genuchten公式计算土壤吸力
        tmp = wfld*(1.0- (Dgl-Ds)/Dg)
        tmp = amax1(tmp, wrsd)

        call SuctionFromMoisture_V(tmp, wsat,wrsd,watern,alpha,ps2)    ! m
        avz(nuz)=0.5*(deltz(nuz)+(Dgl-Ds))                             ! m
        avk(nuz) = conductivity_V(k0(nuz),wsat,wrsd,watern,w(nuz))
        !采用Van Genuchten公式计算非饱和土壤导水率
        tmp      = 0.5*(deltz(nuz) + (Dgl-Ds))
	  rech     = -avk(nuz)*(ps2-ps1)/tmp + avk(nuz)                  ! m/s
        !Richards计算入渗
        rech = rech*dt
c
c      print *,isub
c      pause

      
        if(rech .ge. 0.0) then
          tmp = (w(nuz)-wfld)*deltz(nuz)
          tmp = amax1(tmp, 0.0)
          tmp2= Dg*GWcs - GWst
	    rech= amin1(rech, tmp, tmp2)
c	    if(rech.lt.-0.1E-6)  then
c		  print *,"rech<0",rech,tmp,tmp2,GWst, GWcs
c	    end if
        !由mqh注释
        else
		rech   = amax1(rech, -1.0*GWst)
          tmp    = amax1( 0.0, (wsat-w(nuz))*deltz(nuz) )
          rech   = amax1(rech, -1.0*tmp)
	    if(rech.gt.1.0E-10) print *,"rech>0" ,rech,tmp,GWst, GWcs
        endif
cc        total_rech=total_rech+rech  !mqh
        if(abs(rech) .lt. 0.1E-20) goto 30
c
c        if(rech.lt.0.0) print *, rech, w(nuz), (ps2-ps1),avz(nuz),'1212'
c
        w(nuz) = w(nuz) - rech/deltz(nuz)
	  GWst   = GWst + rech
	  Dgl    = Dgl  - rech / GWcs
	  if(Dgl.gt.(Ds+Dg+1.0E-8)) then
		print *,"Dgl>Ds+Dg",Dgl,Ds,Dg,rech,ps2,ps1,avk(nuz)
	  end if
	  if(Gwst.lt.0) then
		print *,"Gwst<0",Gwst
	  end if
	  if(Dgl.gt.(Ds+Dg+1.0E-8) .or. Gwst.lt.0 ) then
          w(nuz) = w(nuz) + rech/deltz(nuz)
	    GWst   = GWst - rech
	    Dgl    = Dgl  + rech / GWcs
	  end if
	
        	     
	
c	  if((GWst- Dg*GWcs) .gt. 0.1E-5)  
c     :      print *, "wrong in recharge:", rech, GWst,Dg,Dgl
        !由mqh 注释
 
30        rech = 0.0
	endif
	 
cc	 if(isub.eq.153.and.idc.eq.365.and.hour.eq.24)  then
cc      print *,'total_rech=',isub,total_rech !mqh
c      pause
cc      endif

c       print *, 'ok1 of runoff'
c#######################################################################
c
c     calculate infiltration
c
c#######################################################################
c
	inf    = 0.0
	excess = 0.0
      if(Sst.gt.0.0) then

cc     un understand

cc	if(isub.ge.44 .and. isub.le.94) then  !infiltration excess
cc
cc	  tmp    = 0.0
cc	  tmp2   = (dt*k0(1)/Sst)  * sqrt(w(1)/wsat)
cc        tmp    = Sst*exp(-tmp2)
cc	  excess = amin1(tmp, Sst)
cc	  Sst    = Sst - excess
c	if(excess .gt. 0.1) print *, excess, month, day
cc	endif
  
	  inf  = amax1(k0(1)*dt, 0.0)
	  inf  = amin1(inf,      Sst)
        Sst  = Sst - inf
        tmp  = w(1) + inf/deltz(1)
        w(1) = amin1(wsat, tmp)
        Sst  = Sst + (tmp-w(1))*deltz(1)  

	  Sst    = Sst + excess
	  excess = 0.0
	  
c	  print *,inf
c        pause
	  if(Sst .lt. 0.1E-9) Sst=0.0
	  inf = 0.0  !问题1
       
	endif

c       print *, 'ok2 of runoff'

c#######################################################################
c
c     calculate ksoil
c
c#######################################################################
c
	do i=1,nuz
	  ksoil(i)=conductivity_V(k0(i),wsat,wrsd,watern,w(i))
	  call SuctionFromMoisture_V(w(i),
     :       wsat,wrsd,watern,alpha,suct(i))
	end do
c
C#######################################################################
c
c    calculation of inter-layer exchanges of water due to gravitation           
c    and hydraulic gradient.                       
c
C#######################################################################
c
	do i=1,nuz-1
	  avz(i)=0.5*(deltz(i)+deltz(i+1))
c
        tmp    = (w(i)*w(i+1)*(deltz(i)+deltz(i+1)))/
     :                        (w(i)*deltz(i)+w(i+1)*deltz(i+1))    !看不懂
        tmp2   = (k0(i)*deltz(i)+k0(i+1)*deltz(i+1))/
     :                        (deltz(i)+deltz(i+1))
        avk(i) = conductivity_V(tmp2,wsat,wrsd,watern,tmp)

c        if(w(i) .lt. w(i+1)) avk(i) = 0.1*avk(i)   ! used in SiB2

	  if(abs(w(i+1)-w(i)) .lt. 0.1E-5) then
		dpdw(i) = 0.0
	  else
	    dpdw(i)=(suct(i+1)-suct(i))/(w(i+1)-w(i))
	  endif
      end do

c       print *, 'ok3 of runoff'

C#######################################################################
c
c     backward implicit calculation of flows between soil layers.      !看不太懂         
c
C#######################################################################
c       print *, 'ok300 of runoff'

	 temA(1)=0.0

c       print *, 'ok301 of runoff'

	temB(1)=1.0+avk(1)*dpdw(1)/avz(1) * 
     :        (dt/deltz(1)+dt/deltz(2))
 
c       print *, 'ok301 of runoff'

	temC(1)=-avk(1)*dpdw(1)/avz(1) * dt/deltz(2)
c       print *, 'ok302 of runoff'

	temD(1)=avk(1)*dpdw(1)/avz(1) * 
     :        (w(1)-w(2)+
     :         inf*dt/deltz(1))+avk(1)
 
c       print *, 'ok303 of runoff'

	do i=2,nuz-2, 1
        temA(i)=-avk(i)*dpdw(i)/avz(i) * dt/deltz(i)
        temB(i)=1.0+avk(i)*dpdw(i)/avz(i) * 
     :          (dt/deltz(i)+dt/deltz(i+1))
        temC(i)=-avk(i)*dpdw(i)/avz(i) * dt/deltz(i+1)
        temD(i)=avk(i)*dpdw(i)/avz(i) * 
     :          (w(i)-w(i+1)) + avk(i)
      end do
c       print *, 'ok31 of runoff'


      temA(nuz-1)=-avk(nuz-1)*dpdw(nuz-1)/avz(nuz-1) *
     :             dt/deltz(nuz-1)
      temB(nuz-1)=1.0+avk(nuz-1)*dpdw(nuz-1)/avz(nuz-1)*
     :             (dt/deltz(nuz-1)+dt/deltz(nuz))
	temC(nuz-1)=0.0
	temD(nuz-1)=avk(nuz-1)*dpdw(nuz-1)/avz(nuz-1) *
     :            (w(nuz-1)-w(nuz)+ 
     :            rech*dt/deltz(nuz)) + avk(nuz-1)
c
c       print *, 'ok32 of runoff'

      if(nuz.eq.3) then
       temq(1)=(temD(2)-temD(1)*temB(2)/temC(1))/
     :         (temA(2)-temB(1)*temB(2)/temC(1)) 
       temq(2)=(temD(2)-temD(1)*temA(2)/temB(1))/
     :         (temB(2)-temA(2)*temC(1)/temB(1))
	else
	  call TRDIG (nuz-1, temA, temB, temC, temD, temq, MARK)
	  if(MARK .ne. 1) print *, 'wrong in solving  linear equations...'
	endif
c       print *, 'ok33 of runoff'

c
C#######################################################################
c
c     update moisture of each soil moisture layer due to layer interflow.         
c
C#######################################################################

c       print *, 'ok4 of runoff'

	DO i = 1, nuz-1 
        qmax = amax1((w(i)-wrsd-1.0E-6) * deltz(i)/dt,0.0)   
        qmin = amin1(-(w(i+1)-wrsd-1.0E-6) * deltz(i+1)/dt,0.0)

	  if(temq(i).lt.0.0 .and. temq(i+1).gt.0.0) then
          qmin = 0.4*qmin
        endif
         if(i.ge.2) then
	  if(temq(i-1).lt.0.0 .and. temq(i).gt.0.0) then
          qmax = 0.5*qmax
          end if
        endif

        temq(i) = amin1(temq(i), qmax)                                              
        temq(i) = amax1(temq(i), qmin)
                                   
        w(i)   =  w(i)   - temq(i) * dt / deltz(i) 
        w(i+1) =  w(i+1) + temq(i) * dt / deltz(i+1)     
      END DO
      w(nuz) = w(nuz) - rech*dt/deltz(nuz)
c
c     saturation excess
      GWkey  = 0
      excess = 0.0
      excess = (w(nuz) - wsat) * deltz(nuz)
      excess = amax1(excess, 0.0)
      w(nuz) = w(nuz) - excess / deltz(nuz)
      if(abs(Dgl-Ds).lt.0.1E-5 .and. abs(w(nuz)-wsat).lt.0.1E-5) then 
        Dgl   = Ds - deltz(nuz)
        GWkey = 1
      else 
        GWkey = 0
      endif
      do i=nuz-1, 1, -1
	  w(i)   = w(i) + excess / deltz(i)
        excess = (w(i)-wsat) * deltz(i)
        excess = amax1(excess, 0.0)
        w(i)   = w(i) - excess /deltz(i)
        if(GWkey.eq.1 .and. abs(w(i)-wsat) .lt. 0.1E-5) then 
          Dgl = Dgl-deltz(i)
        else
          GWkey = 0
        endif
	end do
c
	Sst = Sst - inf*dt + excess
c     由mqh注释  标记3  2014/4/14
c	do i = 1, nuz 
c	  if(w(i) .lt. wrsd-0.0001 .or. w(i) .gt. wsat+0.00001)  then
c      print *, i, w(i), wrsd, wsat, temq(i),'22222'
c      endif
c      end do

C#######################################################################
c
c     prevent negative values of www(i)                                         
c
C#######################################################################
c
      DO i = 1,nuz-1                                                           
        deficit  = amax1 (0.0, (1.0e-10 - w(i)))                                  
        w(i)   = w(i) + deficit                                              
        w(i+1) = w(i+1) - deficit * deltz(i)/deltz(i+1)                 
      END DO
      w(nuz) = amax1(w(nuz), 1.0e-10)                                         
c
c       print *, 'ok5 of runoff'

c#######################################################################
c
c    calculation of exchange rate between groundwater and river: qsub                            
c                                                                               
c    Yang D. et al., Hydrological Processes, 14, 403-416, 2000                                          
c
c####################################################################### 
c
	isat    = nuz+1
      qsub    = 0.0
      tmpqsub = 0.0
c
      if( (Dgl-Ds) .ge. 0.0 ) then ! groundwater table is bolow the toplayer
	  avkg = kg
        call gwriv(Dgl,length,slope,Ds,Dg,Dr,Drw,avkg,tmpqsub) ! m3/s/m
        qsub = tmpqsub * dt / length			                 ! m
        if(abs(qsub) .lt. 0.1E-20) qsub = 0.0
c
      else                         ! groundwater table reaches the toplayer
        if(Drw.ge.Dr .and. Dgl.le.0.0) goto 1010
        tmp_D = 0.0
        avkg  = 0.0
        do i = nuz, 1, -1
          tmp_D = tmp_D + deltz(i)
		if( (Ds-tmp_D) .le. 0.5) then
            avkg  = avkg + anik*k0(i)*deltz(i)
		else
            avkg  = avkg + k0(i)*deltz(i)
		endif
          if(tmp_D .ge. (Ds-Dgl)) goto 50
        end do
50      if(Dr .gt. Ds) then
          tmp_D = tmp_D + Dr-Ds
          avkg  = avkg + kg*(Dr-Ds)
        endif
        avkg  = avkg / tmp_D
        call gwriv(Dgl,length,slope,Ds,Dg,Dr,Drw,avkg,tmpqsub) ! m3/s/m
        qsub = tmpqsub * dt / length			                 ! m
        if(abs(qsub) .lt. 0.1E-20) qsub = 0.0
      endif
c
c     Renewing the groundwater table and toplayer soil moisture
      if(qsub .lt. 0.0) then      ! River infiltrates into aquifer
        GWst  = GWst - qsub
c        print *,qsub
        if(GWst .le. Dg*GWcs) then
          Dgl   = Dgl + qsub/GWcs
	    if(Dgl .lt. 0.0) print *, 
     :       'wrong in renewing groundwater table ... 1', Dgl
        else
          tmp    = GWst - Dg*GWcs
          GWst   = Dg*GWcs
          Dgl    = Ds
          w(nuz) = w(nuz) + tmp / deltz(nuz)
          excess = (w(nuz) - wsat) * deltz(nuz)
          excess = amax1(excess, 0.0)
          w(nuz) = w(nuz) - excess / deltz(nuz)
          GWkey  = 0
          if(abs(w(nuz)-wsat) .lt. 0.1E-5) then 
            Dgl = Ds - deltz(nuz)
            GWkey = 1
          else 
            GWkey = 0
          endif
          do i=nuz-1, 1, -1
	      w(i)   = w(i) + excess / deltz(i)
            excess = (w(i)-wsat) * deltz(i)
            excess = amax1(excess, 0.0)
            w(i)   = w(i) - excess /deltz(i)
            if(GWkey.eq.1 .and. abs(w(i)-wsat) .lt. 0.1E-5) then 
              Dgl = Dgl-deltz(i)
            else
              GWkey = 0
            endif
	    end do
        endif
	elseif(qsub .gt. 0.0) then                ! Aquifer flows out
       if( Dgl .ge. Ds-0.01 ) then
          GWst  = GWst - qsub
	    Dgl   = Dgl + qsub/GWcs
	    if(Dgl.gt.(Ds+Dg) .or. GWst.lt.0.0) then
	      qsub= qsub+GWst
		  Dgl = Ds+Dg
            GWst= 0.0
	    end if
	  else
          tmp_D = 0.0
          isat  = nuz+1
          do i = nuz, 1, -1
            tmp_D = tmp_D + deltz(i)
            if( abs(tmp_D - (Ds-Dgl)) .le. 0.1E-3) goto 100
          end do
          goto 110
100       if(abs(w(i)-wsat) .le. 0.1E-3)   isat = i
          goto 120
110       if(abs(w(i+1)-wsat) .le. 0.1E-3) then
            isat = i+1
	    end if   
c	    print *,"wrong in Dgl & qsub>0",isat,Dgl
c		
120       tmp = 0.0
          do i = isat, nuz
            if(w(i)-wsat .gt. 0.1E-3) 
     :      print *, "wrong in GW flow", isat, i,Dgl,w(i),wsat
            tmp  = tmp + GWcs*deltz(i)
            w(i) = w(i) - GWcs
            Dgl  = Dgl + deltz(i)
            if(tmp .ge. qsub) goto 150
          end do
150       if(i.le.nuz) w(i) = w(i) + (tmp-qsub)/deltz(i)
          if(i.gt.nuz) then
            GWst  = GWst - (qsub-tmp)
	      Dgl   = Dgl + (qsub-tmp)/GWcs
	    end if
	    	 if(i.gt.1) then
          if(  w(i-1)-wsat .gt. 0.1E-3) 
     :       print *,"w(i)>wsat", i-1, w(i-1),wsat,isat,nuz
          if(Dgl.gt.Ds .and. abs(GWst-Dg*GWcs).lt.0.1E-3) Dgl=Ds
          end if
        endif
	endif
	
cc	total_qsub=qsub+total_qsub   !added by  mqh,2014.4.24
cc      if(isub.eq.153.and.idc.eq.365.and.hour.eq.24) then
cc        print *,'total_qsub=',total_qsub
c        pause
cc        endif
c       print *, 'ok6 of runoff'

c      goto 1010
c
c#######################################################################
c
c   Subsurface flow from top saturated zone above the groundwater level
c
c#######################################################################
c
      j = min(isat, nuz)
	ssub=0.0 ! added by gaobing
c	j = nuz
	tmp = 0.0
	do 200 i = 1, j
	  tmp = tmp + deltz(i)
        if((w(i)-wfld) .gt. 0.1E-3) then

          ksoil(i) = conductivity_V(k0(i),wsat,wrsd,watern,w(i))
		tmpqsub  = ksoil(i)*slope*dt
		if(tmp .le. 0.5) tmpqsub = anik*tmpqsub

		f = 100.0
		if(tmp.le.1.5) f = -alog(0.1 /1.0)

cc		if(isub.ge.8 .and. isub.le.22) then
cc		  f = -alog(0.5 /1.0) /Ds
cc		endif

cc		if(isub.ge.44 .and. isub.le.128) then   !Weihe 
cc		  if(isub.le.108 .and. tmp.le.1.5) f = (-alog(0.5 /1.0)) /1.5
cc		  if(isub.ge.109 .and. tmp.le.0.8) f = (-alog(0.2/1.0)) /1.0
cc		endif

		tmpqsub = tmpqsub * 
     :              (deltz(i)+exp(-f*tmp)*ss_f*length)/length
c		if(tmpqsub .le. 0.1E-30) tmpqsub = 0.0
c
          tmp     = (w(i)-wfld)*deltz(i)
          tmpqsub = amin1(tmpqsub, tmp)
          tmpqsub = amax1(tmpqsub, 0.0)
          w(i)    = w(i) - tmpqsub/deltz(i)
	    qsub    = qsub + tmpqsub
	    ssub    = ssub+tmpqsub  ! added by gaobing     !ssub只包括了Subsurface flow from top saturated zone above the groundwater level，而qsub还包括了之前的由地下水流向河道的水
cc	    total_ssub=  total_ssub+ssub
cc	  if(isub.eq.153.and.idc.eq.365.and.hour.eq.24) then
cc       print *,'total_ssub=',total_ssub
c        pause
cc        endif  !mqh
	  endif
200   continue
c             print *, 'ok7 of runoff'

1010    ssub = ssub*length / dt              ! m --> m3/sec/m
        qsub = qsub*length / dt              ! m --> m3/sec/m
c


      RETURN                                                                    
      END 
c
c
c
C#######################################################################
c
c     calculate suction from soil moisture by Van Genuchten's equation                                         
c
C#######################################################################
c
	subroutine SuctionFromMoisture_V(w,wsat,wrsd,n,alpha,ps)
      implicit none
	real w,wsat,wrsd,alpha,m,n,se,ps
	real tmpe, tmpps,tmpw
	m    = 1.0-1.0/n
      tmpw = w
	if(tmpw.ge.wsat) then
	  ps=0.0
	elseif(tmpw.lt.wsat) then
	  if(tmpw.lt.wrsd+0.001) tmpw=wrsd+0.001
	  se=(tmpw-wrsd)/(wsat-wrsd)
	  tmpe=se**(1.0/m)
	  tmpps=-(1.0/tmpe-1.0)**(1.0/n)/alpha
        ps = tmpps/100.0 ! cm->m
c     print *,ps
	end if
	return
	end
c
C#######################################################################
c
c     calculate soil moisture from suction by Van Genuchten's equation                                         
c
C#######################################################################
c
	subroutine MoistureFromSuction_V(w,wsat,wrsd,n,alpha,ps)
      implicit none
	real w,wsat,wrsd,alpha,m,n,ps
	real se,tmpe,tmpps
	tmpps=100.0*ps ! m->cm
	m=1.0-1.0/n
	tmpe=1.0+(alpha*abs(tmpps))**n
	se=(1.0/tmpe)**m
	w=se*(wsat-wrsd)+wrsd
!      print *,w
	 if(w .gt. wsat) w = wsat
	 if(w .lt. wrsd) w = wrsd
	 return
	 end
c
C#######################################################################
c
c      calculate soil hydraulic conductivity by Van Genuchten's equation                                 
c
C#######################################################################
c
	function conductivity_V(k0,wsat,wrsd,n,w)
	real k0,wsat,wrsd,n,m,w
      real tmpw
	m=1.0-1.0/n
      tmpw = w
	if(tmpw.gt.wsat) tmpw=wsat
	
      if(tmpw.le.wrsd) then
	  tmpw=wrsd+1.0E-10
	elseif(tmpw.lt.wrsd+0.0001)then
	  tmpw=wrsd+0.0001
	end if

      se=(tmpw-wrsd)/(wsat-wrsd)
	conductivity_V=k0*sqrt(se)*(1.0-(1.0-se**(1.0/m))**m)**2
	if(conductivity_V.lt.0.0 .or. conductivity_V.gt.k0) 
     $  print *,'wrong in calculating conductivity',w,conductivity_V
	return
	end
c
C#######################################################################
c
c    calculation of exchange rate between groundwater and river: qsub
c                                                                               
c    Yang D. et al., Hydrological Processes, 14, 403-416, 2000                                          
c
C#######################################################################	 
	subroutine gwriv(Dtg,length,slope,Ds,Dg,Dr,Drw,kg,Q)
c
c     The datum is sitted at the bottom of the unconfined aquifer
c     Varibles:
c     Dtg:      depth to groundwater (m)
c     length:   length of hillslope (m)
c     slope:    slope of hillslope (m)
c     Ds:       depth of top soil (m)
c     Dg:       depth of unconfined groundwater acquifer below topsoil (m)
c     Dr:       depth of river (m)
c     Drw:      depth of river water (m)
c     kg:       hydraulic conductivity (m/sec)
c     conlen:   contact length between the river and aquifer (m)
c     grad:     gradient of water head
c     Q:        discharge exchanged between aquifer and river (m^3/sec)
c
      implicit none
      real Dtg,length,slope,Ds,Dg,Dr,Drw,kg,Q
	real H1,H2,hs1,hs2,grad, conlen, Hrd
c	print *,'length,slope,Dtg,Ds,Dg,Dr,Drw',length,slope,Dtg,Ds,Dg,Dr,Drw
c	if(drw.ge.0.1) 	then
c	print *,'length,slope,Dtg,Ds,Dg,Dr,Drw',length,slope,Dtg,Ds,Dg,Dr,Drw
c      pause
c      endif
	Drw = amax1(Drw, 0.0)
	Drw = amin1(Dr, Drw)

	Hrd = Ds+Dg-Dr                          ! distance from datum to riverbed (m)
	if(Dtg .lt. Dr-Drw) then      !chengyashui
        H1  = 0.5*length*slope + Ds+Dg-Dtg
      else         !feichengyashui
        H1  = sqrt(0.5*length*slope) + Ds+Dg-Dtg   ! waterhead of groundwater (m)
	endif
      hs1=H1-Hrd                               ! saturated acquifer depth (m)

	if(Hrd.ge.0.0) then
        H2 =Hrd+Drw                      ! waterhead of river (m)
        hs2=Drw                          ! water depth in river (m)
	else
	  hs2=amax1(0.0, Drw)
        H2 =amax1(Hrd+Drw, 0.0)        
	endif
	grad   =(H1-H2)/(0.5*length)     ! gradient of waterhead
	conlen =0.5*(abs(hs1)+hs2)       ! contact length between the river
                                         ! and aquifer	(m)
	Q=kg*grad*conlen			  ! discharge per unit width of a hillslope
                                    ! (m3/sec/m)
c      if (H2.ge.H1.and.q.le.-10e-8) then
c      if(H2.ge.H1) then
c      print *,"H1=",H1,H2,slope,"Drw=",Drw,grad,"q=",Q,conlen
c      PAUSE
c      endif
	if(Drw .le .5E-3 .and. Q.lt.0.0) Q=0.0
	return
	end
c
c
C*****************************************************************
C                                                                *
C     This routine solves a linear system of equations           *
C                  A * X = RS                                    *
C     for a tridiagonal, strongly nonsingular matrix A.          *
C     The matrix is given by the three vectors DL, DM and DU     *
C     which designate the lower co-diagonal, the diagonal and    *
C     the upper co-diagonal elements of A, respectively.         *
C     The system of equations has the form :                     *
C                                                                *
C     DM(1) * X(1) + DU(1) * X(2)                      = RS(1)   *
C                                                                *
C     DL(I) * X(I-1) + DM(I) * X(I) + DU(I) * X(I+1)   = RS(I)   *
C            for I = 2, ..., N-1, and                            *
C                                                                *
C     DL(N) * X(N-1) + DM(N) * X(N)                    = RS(N)   *
C                                                                *
C                                                                *
C     INPUT PARAMETERS:                                          *
C     =================                                          *
C     N    : number of equations, N > 2                          *
C     DL   : N-vector DL(1:N); lower co-diagonal of A            *
C            DL(2), DL(3), ... ,DL(N)                            *
C     DM   : N-vector DM(1:N); the diagonal of A                 *
C            DM(1), DM(2), ... , DM(N)                           *
C     DU   : N-vector DU(1:N); upper co-diagonal of A            *
C            DU(1), DU(2), ... , DU(N-1)                         *
C     RS   : N-vector RS(1:N); the right hand side               *
C                                                                *
C                                                                *
C     OUTPUT PARAMETERS:                                         *
C     ==================                                         *
C     DL   :)                                                    *
C     DM   :)                                                    *
C     DU   :) these are overwritten with auxiliary vectors       *
C     RS   :)                                                    *
C     X    : N-vector X(1:N), containing the solution of the     *
C            system of equations                                 *
C     MARK : error parameter                                     *
C            MARK= 1 : everything is o.k.                        *
C            MARK= 0 : the matrix A is not strongly nonsingular  *
C            MARK=-1 : error on N: N <= 2                        *
C                                                                *
C     NOTE: if MARK = 1, the determinant of A can be calculated  *
C           after this subroutine has run as follows:            *
C              DET A = DM(1) * DM(2) * ... * DM(N)               *
C                                                                *
C----------------------------------------------------------------*
      SUBROUTINE TRDIG (N,DL,DM,DU,RS,X,MARK)
c	parameter (N=3)
      REAL DL(1:N),DM(1:N),DU(1:N),RS(1:N),X(1:N)
      MARK = -1
      IF (N .LT. 3) goto 999
C
C   checking for strong nonsingularity with N=1
      MARK = 0
      ROW = ABS(DM(1)) + ABS(DU(1))
      IF (ROW .EQ. 0.0E0) goto 999
      D = 1.0E0/ROW
      IF (ABS(DM(1))*D .LE. 1.0E-20) goto 999
C
C   factoring A while checking for strong nonsingularity
      DL(1) = 0.0E0
      DU(N) = 0.0E0
      DU(1) = DU(1)/DM(1)
      DO I=2,N,1
         ROW = ABS(DL(I)) + ABS(DM(I)) + ABS(DU(I))
         IF (ROW .EQ. 0.0E0) goto 999
         D = 1.0E0/ROW
         DM(I) = DM(I) - DL(I) * DU(I-1)
         IF (ABS(DM(I))*D .LE. 1.0E-20) goto 999
         IF (I .LT. N) THEN
            DU(I) = DU(I)/DM(I)
         ENDIF
      END DO
      MARK=1
c
999   continue
C
C  If MARK = 1, update the right hand side and solve via backsubstitution
      IF (MARK .EQ. 1) THEN
        RS(1) = RS(1)/DM(1)
        DO I=2,N,1
          RS(I) = (RS(I) - DL(I) * RS(I-1)) / DM(I)
        END DO
C       backsubstitution
        X(N) = RS(N)
        DO I=N-1,1,-1
          X(I) = RS(I) - DU(I) * X(I+1)
        END DO
      ENDIF
c	print *,MARK, X
      RETURN
      END
