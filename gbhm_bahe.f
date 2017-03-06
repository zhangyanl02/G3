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

