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
c     ######                  �޸��� 2006.5.17                    ######
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
      ! ���

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
cc	 �޸���2013.05.09  gaobing								! ��ǩ1
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
	!���
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
cc	�޸���2006.4.14
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
cc		Dr(isub,iflow) = 1.5 * Dr(isub,iflow)					! ��ǩ6
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
c		  !added by mqh  ������������
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
c	enddo  !������ѭ��
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
	  
	  else   !��� ��ǰ��˵��>1998�ģ�1982~1998����ʣ�µ�else
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


      ! ���
      if(iflai.eq.1)  then   !added by mqh,��û���õ�lai
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
      end if  !added by mqh,��Ϊ����lai����if�Թ����


	  

c
c
c#######################################################################
c
c	    read meteorological data
c
c#######################################################################

cc	2222	ɾ����2006.4.14	 ! skip input meteorological data

c
c	   print *,"After read data",basinarea
c        print *,"hehe"
        idc=0
	  ihc=0
	  do 777 im = 1, 12    !��ǣ���ǩ
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
	if(Drw(isub,iflow) .gt. 0.5*Dr(isub,iflow)) then		! ��ǩ2
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
	!���2   qr2
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
		  !���
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
            write(9, '(3i6,f16.3)') year,im,id, Qd(isub,j)  !���1��Qd(isub,j)
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
		  !adjust the area difference ��������������
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
        !���  
        !end of added
        
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
       !added by mqh����������Сʱ����  ��mqh��2014.5.3
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

c
c
c

c


c

c
c

