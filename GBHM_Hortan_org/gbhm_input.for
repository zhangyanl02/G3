c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######      calculete potential evapotranspiration          ######
c     ######                                                      ######
c     ######             Developed by Dr. Chong Li                ######
c     ######           Department of Civil Engineering            ######
c     ######                University of Tokyo                   ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################


c*******************************************************************************
c	To prepare and initial data for interpolation from gauge to grid
c*******************************************************************************

      subroutine Initial_input()
      implicit none
      include 'dims2.inc'
c
c#######################################################################
c     Variable Declarations.
c#######################################################################   
c      integer nr, nc
c      parameter (nr = nrow, nc = ncol)     ! This data is in matrix 
      integer ntypes	
      parameter (ntypes = nv)	               ! number of land-use types

      real      G_elev(nrow,ncol)              ! matrix for hight of site [m]
	real      zone(nrow,ncol)              ! matrix for special zone
      real      lon(nrow,ncol)               ! matrix for longtitude
      real      lat(nrow,ncol)               ! matrix for latitude

      real pi
      parameter (pi=3.1415926)

      integer ir,ic,it,i
  	character*2    ch2        ! a 2-character variable
      character*5 ncols,nrows
      character*9 xllcorner,yllcorner
      character*8 cellsize
      character*12 nodata
	integer nnr, nnc
	real xx0, yy0
	real gridsize
	real znodata
	real znodata1


	character*200, infile, in_dir,dem_dir,file_dir2,infile2 !file_dir2,infile2,mqh
	
c
	real x0(Maxgauge),y0(Maxgauge),elev(Maxgauge)
	real x(nrow,ncol),y(nrow,ncol)

	real land_use(nrow,ncol,nv)
c	real Hbeta(5,12,4),Tmbeta(12,4)

	integer G_id(Maxgauge),G_zone(Maxgauge)
	real xmin,xmax,ymin,ymax
	integer Y_start(Maxgauge), M_start(Maxgauge)
      integer Y_end(Maxgauge),M_end(Maxgauge),D_end,GG_id,T_zone
	integer Y_start1,M_start1,Y_end1,M_end1, D_start
	integer N_gauge,ia1,ia2,ic1,ic2
	real  G_eps(Maxgauge),p_zone     ! added by gaobing 2013 
c
      common / E_con /lon,lat,G_elev,zone,x0,y0,G_id,x,y,
     :          G_zone,Y_Start,M_start,Y_end,M_end,elev,N_gauge,G_eps

	common / asc / ncols,nrows,xllcorner,yllcorner,cellsize,
     :           nodata, nnr, nnc, xx0, yy0, gridsize,znodata

cc      added by xujijun
	real x1(Maxgauge),y1(Maxgauge),z1(Maxgauge)
	integer Y1_start(Maxgauge), M1_start(Maxgauge)
      integer Y1_end(Maxgauge),M1_end(Maxgauge),GG_id1,N_gauge1

	integer G_id1(Maxgauge),TF1(Maxgauge),N_gauge2
      common / E_con1/x1,y1,z1,G_id1,Y1_Start,M1_start,
     $      Y1_end,M1_end,N_gauge1

      !added by mqh
      integer GG_id2  !added by mqh
      integer l  !mqh,year loop
      integer startofmonth(maxgauge,20),startofday(maxgauge,20),
     $ endofmonth(maxgauge,20),endofday(maxgauge,20)  !mqh
      integer GG_id3,year2,startofmonth2,startofday2,
     $      endofmonth2,endofday2,iy   !mqh  
     
      common/start_end/startofmonth,startofday,endofmonth,endofday

cc    end of added

	  nnr = nrow
        nnc = ncol
c
	  in_dir = '../parameter/'
	  dem_dir = '../dem/'

        print *, 'Start to read meteorological data'
c	  pause

	  infile='longlat.txt'
	  call strlen(infile, ia1,ia2)
	  call strlen(dem_dir, ic1,ic2)
	  open(11, file=dem_dir(ic1:ic2)//infile(ia1:ia2), status='old')
        do ir = 1, nnr
          do ic = 1, nnc
            read(11,*) lon(ir,ic),lat(ir,ic)
          end do
        end do
        close(11)
       print *, ' ok1, finish of longlat.txt'

	  xmax = -99999999.0
	  xmin =  99999999.0
	  ymax = -99999999.0
	  ymin =  99999999.0
	  infile='lambert.txt'
	  call strlen(infile, ia1,ia2)
	  call strlen(dem_dir, ic1,ic2)
	  open(16, file=dem_dir(ic1:ic2)//infile(ia1:ia2), status='old')
	  do ir = 1, nrow
	    do ic = 1, ncol
	      read(16,*) x(ir,ic), y(ir,ic)
		  if(x(ir,ic) .ge. xmax) xmax = x(ir,ic)
		  if(x(ir,ic) .le. xmin) xmin = x(ir,ic)
		  if(y(ir,ic) .ge. ymax) ymax = y(ir,ic)
		  if(y(ir,ic) .le. ymin) ymin = y(ir,ic)
	    end do
	  end do
	  close(16)
	  print *, 'watershed extent: ', xmin,xmax,ymin,ymax
        print *, ' ok2, finish of lambert.txt'
c
	  infile='elevation.asc'
	  call strlen(infile, ia1,ia2)
	  call strlen(dem_dir, ic1,ic2)
	  open(12, file=dem_dir(ic1:ic2)//infile(ia1:ia2), status='old')
	  read(12,*) ncols,nnc
	  read(12,*) nrows,nnr
	  read(12,*) xllcorner,xx0
	  read(12,*) yllcorner,yy0
	  read(12,*) cellsize,gridsize
	  read(12,*) nodata,znodata
	  do ir = 1 ,nnr
          read(12,*) (G_elev(ir, ic),ic=1,nnc)
	  end do
	  close (12)

	  print *, 'ncols, nrows=', nnc, nnr
	  print *, 'xmin, ymin, cellsize=', xx0, yy0, gridsize
        print *, ' ok3, finish of elevation.asc'
c
	  infile='zone.asc'
	  call strlen(infile, ia1,ia2)
	  call strlen(dem_dir, ic1,ic2)
	  open(13, file=dem_dir(ic1:ic2)//infile(ia1:ia2), status='old')
	  read(13,*) ncols,nnc
	  read(13,*) nrows,nnr
	  read(13,*) xllcorner,xx0
	  read(13,*) yllcorner,yy0
	  read(13,*) cellsize,gridsize
	  read(13,*) nodata,znodata
	  do ir = 1 ,nnr
          read(13,*) (zone(ir, ic),ic=1,nnc)
	  end do
	  Close (13)
	  print *, 'ncols, nrows=', nnc, nnr
	  print *, 'xmin, ymin, cellsize=', xx0, yy0, gridsize
        print *, ' ok4, finish of zone.asc'

cc	修改于2006.4.14
       
cc        print *, ' ok5, finish of lu_class.asc'
c
	  xmax = -99999999.0
	  xmin =  99999999.0
	  ymax = -99999999.0
	  ymin =  99999999.0

	  N_gauge=0
	  infile='gauge_lambert.dat'
	  call strlen(infile, ia1,ia2)
	  call strlen(in_dir, ic1,ic2)
	  open(15, file=in_dir(ic1:ic2)//infile(ia1:ia2), status='old')
	  do ir = 1, Maxgauge
		  read(15,*,end=166) x0(ir), y0(ir),G_id(ir)
		  !标记
c		  print *, x0(ir), y0(ir),G_id(ir)
c      pause
	      N_gauge=N_gauge+1
		  if(x0(ir) .ge. xmax) xmax = x0(ir)
		  if(x0(ir) .le. xmin) xmin = x0(ir)
		  if(y0(ir) .ge. ymax) ymax = y0(ir)
		  if(y0(ir) .le. ymin) ymin = y0(ir)
	  end do
166	  close(15)
        print *, 'all of old gauge', N_gauge
	  print *, 'Gauge extent:   ',  xmin,xmax,ymin,ymax
        print *, ' ok6, finish of gauge_lambert.txt'

	  infile='gauge_zone.dat'
	  call strlen(infile, ia1,ia2)
	  call strlen(in_dir, ic1,ic2)
	  open(17, file=in_dir(ic1:ic2)//infile(ia1:ia2), status='old')
	  do ir = 1 ,N_gauge
          read(17,*,end=333) GG_id,T_zone
	      do ic=1,N_gauge
	        if(GG_id .eq. G_id(ic)) then
				G_zone(ic)=T_zone
			end if
		  end do
	  end do
333	  Close (17)
	  print *, GG_id, G_zone(N_gauge)
        print *, ' ok7, finish of gauge_zone.txt'

	  infile='gauge_period.dat'
	  call strlen(infile, ia1,ia2)
	  call strlen(in_dir, ic1,ic2)
	  open(18, file=in_dir(ic1:ic2)//infile(ia1:ia2), status='old')
	  do ir = 1, N_gauge
	      read(18,*) Y_start1, M_start1,D_start,Y_end1,
     $                M_end1,D_end,GG_id
	      do ic=1,N_gauge
	        if(GG_id .eq. G_id(ic)) then
				Y_Start(ic)=Y_start1
	            M_start(ic)=M_start1
	            Y_end(ic)=Y_end1
	            M_end(ic)=M_end1
			end if
		  end do
	  end do
	  close(18)
c	  print *, Y_start(N_gauge), M_start(N_gauge)
        print *, ' ok8, finish of gauge_period.txt'
	  infile='gauge_elevation.dat'
	  call strlen(infile, ia1,ia2)
	  call strlen(in_dir, ic1,ic2)
	  open(19, file=in_dir(ic1:ic2)//infile(ia1:ia2), status='old')
	  do ir = 1, N_gauge
	      read(19,*) GG_id,elev(ir)
	  end do
	  close(19)
        print *, ' ok9, finish of gauge_elevation.txt'

c	  print *, 'finish of initial input'
c	  pause


cc     added by xujijun 2006.2.25   ! read added new gauge 
       if(NNstation.ge.2) then
	  print *, 'read added new Gauge '

	  N_gauge1=0
	  infile='gauge_lambert_new.txt'
	  call strlen(infile, ia1,ia2)
	  call strlen(in_dir, ic1,ic2)
	  open(18, file=in_dir(ic1:ic2)//infile(ia1:ia2), status='old')
	  !added by mqh
	  xmax = -99999999.0
	  xmin =  99999999.0
	  ymax = -99999999.0
	  ymin =  99999999.0
	  do ir = 1, Maxgauge
		  read(18,*,end=167) x1(ir), y1(ir),G_id1(ir)
	      N_gauge1=N_gauge1+1
		  if(x1(ir) .ge. xmax) xmax = x1(ir)
		  if(x1(ir) .le. xmin) xmin = x1(ir)
		  if(y1(ir) .ge. ymax) ymax = y1(ir)
		  if(y1(ir) .le. ymin) ymin = y1(ir)
	  end do
167	  close(18)
        print *, 'all of new gauge', N_gauge1
	  print *, 'new Gauge extent:   ', xmin,xmax,ymin,ymax
c
	  xmax = -99999999.0
	  xmin =  99999999.0
	  ymax = -99999999.0
	  ymin =  99999999.0

	  infile='gauge_period_new.txt'
	  call strlen(infile, ia1,ia2)
	  call strlen(in_dir, ic1,ic2)
	  open(15, file=in_dir(ic1:ic2)//infile(ia1:ia2), status='old')
	  do ir = 1, N_gauge1
	      read(15,*) Y_start1, M_start1,D_start,Y_end1,
     $                M_end1,D_end,GG_id1
	      do ic=1,N_gauge1
	        if(GG_id1 .eq. G_id1(ic)) then
				Y1_Start(ic)=Y_start1
	            M1_start(ic)=M_start1
	            Y1_end(ic)=Y_end1
	            M1_end(ic)=M_end1
			end if
		  end do
	  end do
	  close(15)
	  print *, 'Read new gauge period'

       endif
! 标记
ccccccccccccccccccccccccccccccc added by gaobing 2013
c	  infile='gauge_eps.txt'
c	  call strlen(infile, ia1,ia2)
c	  call strlen(in_dir, ic1,ic2)
c	  open(17, file=in_dir(ic1:ic2)//infile(ia1:ia2), status='old')
c	  do ir = 1 ,N_gauge
c          read(17,*,end=367) GG_id,p_zone
c	      do ic=1,N_gauge
c	        if(GG_id .eq. G_id(ic)) then
c				G_eps(ic)=p_zone
c			end if
c		  end do
c	  end do
c 367	  Close (17)
cc	  print *, GG_id, G_eps(N_gauge)
c       print *, ' ok, finish of gauge_eps.txt'

cc     end of added

ccccccccccccccccccccccccccccccccccccccccccccccccadded by mqh
      file_dir2 = '../data/gauge_new/'  !mqh
	infile2='start_end.txt'
	
	  call strlen(infile2, ia1,ia2)
	  call strlen(file_dir2, ic1,ic2)
	  open(15, file=file_dir2(ic1:ic2)//infile2(ia1:ia2), status='old')
	  
	  do ir = 1, N_gauge1
	  do l=1991,1999
	      iy=l-1990
	      startofmonth(ir,iy)=0
	          startofday(ir,iy)=0
	           endofmonth(ir,iy)=0
	           endofday(ir,iy)=0
	           enddo
	  do l=2000,2010
	      iy=l-1990
      read(15,*) GG_id3,year2,startofmonth2,startofday2,
     $      endofmonth2,endofday2
	      do ic=1,N_gauge1

c	      print *, G_id1(ic)
c	      pause
	        if(GG_id3 .eq. G_id1(ic).and.year2.eq.l) then
				startofmonth(ic,iy)=startofmonth2
	          startofday(ic,iy)=startofday2
	           endofmonth(ic,iy)=endofmonth2
	           endofday(ic,iy)=endofday2
c           print *,year2,GG_id3,startofmonth(ic,iy),startofday(ic,iy),
c     :   endofmonth(ic,iy),endofday(ic,iy),iy
c      pause    !怎么只读了一个站点的呢？
			end if
		enddo  
	  enddo
c	  		      pause
	  enddo
	  close(15)


cccccccccccccccccccccccccccccccccccccccccccccccccccccc end of added

	end 
c
c
c*******************************************************************************
c	To interpolate meteorological gauge data to 10km grid
c*******************************************************************************
c
      subroutine T_P_E_cal(iyear,im)
      implicit none
      include 'dims2.inc'
c
c#######################################################################
c     Variable Declarations.
c#######################################################################   
c      integer nr, nc
c      parameter (nr = nrow, nc = ncol)     ! This data is in matrix 
      integer ntypes	
      integer land_type !land type
      parameter (ntypes = nv)	               ! number of land-use types
      real      pre(nrow,ncol,31)         ! matrix for precipitation [mm]
      real      tmp(nrow,ncol,31)            ! matrix for mean tempreture	[degC]
      real      t_max(nrow,ncol,31)        ! matrix for maximum tempreture[degC]
      real      t_min(nrow,ncol,31)        ! matrix for minimum tempreture[degC]
      real      cld(nrow,ncol,31)            ! matrix for cloud cover	[oktas*10]
      real      wind(nrow,ncol,31)           ! matrix for wind speed at 10_m
      real      sunshine(nrow,ncol,31)       ! matrix for sunshine time [hr]
      real      RH(nrow,ncol,31)             ! matrix for relative humidity
      real      G_elev(nrow,ncol)            ! matrix for hight of site [m]
	real      zone(nrow,ncol)              ! matrix for special zone
      real      land_use(nrow,ncol,ntypes)   ! matrix for land-use
      real      lon(nrow,ncol)               ! matrix for longtitude
      real      lat(nrow,ncol)               ! matrix for latitude
	real      evap(nrow,ncol,31,nv)        ! potential evapotranspiration 
	real      snow0(nrow,ncol,nv)   ! snow depth(mm water) of the initial year
	real      SH(nrow,ncol,31)             ! added by gaobing

      real pi
      parameter (pi=3.1415926)

      integer Jd,ir,ic,it
	integer iyear,im,iday,nday,imm
      integer mmonth(12)
      real Ep, R_n
	real latitude

	character*5    ch5        ! a 5-character variable
	character*4    ch4        ! a 4-character variable
	character*2    ch2        ! a 4-character variable
      CHARACTER*9    ch9
      CHARACTER*15   ch15

	
      integer nnr, nnc

      data mmonth/31,28,31,30,31,30,31,31,30,31,30,31/
c
      character  file_dir*200   ! directory for parameters
	real x0(Maxgauge),y0(Maxgauge),z0(Maxgauge),elev(Maxgauge)
	real x(nrow,ncol),y(nrow,ncol)
	real Tm(Maxgauge,31),Tu(Maxgauge,31),Td(Maxgauge,31)
	real Pr(Maxgauge,31),Rhh(Maxgauge,31),Cloud(Maxgauge,31)
	real Sun(Maxgauge,31),Wd(Maxgauge,31)
	real Hbeta(5,12,4),Tmbeta(12,4)
	integer G_id(Maxgauge),G_zone(Maxgauge),TF(Maxgauge)
	real    G_eps(Maxgauge)
	integer Y_start(Maxgauge), M_start(Maxgauge)
      integer Y_end(Maxgauge),M_end(Maxgauge),GG_id
	integer Y_start1,M_start1, D_start      !H_start added by mqh, for torrential flood
 	integer Y_start2,M_start2, D_start2 ,H_start !added by mqh
      character*5 ncols,nrows
      character*9 xllcorner,yllcorner 
      character*8 cellsize
      character*12 nodata
	real xx0, yy0
	real gridsize
	real znodata

	integer mm,ib1,ib2
	integer N_gauge,TFF,ir1
c
	real dd0,dv
	real Tmean,Tmax,Tmin,RH1,Windm,prec,Tsun,Snow,prec2  !prec2 added by mqh, for torrential flood
 
	
      common / E_con /lon,lat,G_elev,zone,x0,y0,G_id,x,y,
     :          G_zone,Y_Start,M_start,Y_end,M_end,elev,N_gauge,G_eps
	common / weather / pre,t_min,t_max,evap,snow0, prehour

	common / asc / ncols,nrows,xllcorner,yllcorner,cellsize,
     :           nodata, nnr, nnc, xx0, yy0, gridsize,znodata

cc      added by xujijun
      character  file_dir1*200   
	character*10    ch10        ! a 10-character variable
      character*8     ch8 
      CHARACTER*30    ch30
	
	real x1(Maxgauge),y1(Maxgauge),z1(Maxgauge)
      real xtemp(Maxgauge),ytemp(Maxgauge),ztemp(Maxgauge)
	integer Y1_start(Maxgauge), M1_start(Maxgauge), N_gaugex,newgauge
      integer Y1_end(Maxgauge),M1_end(Maxgauge),GG_id1,N_gauge1,N_gauge2
      integer GG_id2  !mqh
	integer G_id1(Maxgauge),TF1(Maxgauge),idh
	real   Pr1(Maxgauge,31),prhour(Maxgauge,31,24)  
      integer   np_day(Maxgauge),nh(Maxgauge,31)
	real   prehour(nrow,ncol,31,24), pretemp
	real    result_m, result_month(nrow,ncol), znodata1 
      common / E_con1/x1,y1,z1,G_id1,Y1_Start,M1_start,
     $      Y1_end,M1_end,N_gauge1
    !added by mqh.输出面平均雨量
      integer countt(nc),ih
      real result_hour(31,24,nc)
      save result_hour
c      save countt
      integer       psubbasin(nc)     ! subbasin p number ! added by gaobing 2013
      integer       nbasinup(nc,4)    ! upbasin of n number ! added by gaobing 2013
      common /river4    /  psubbasin,nbasinup ! added by gaobing 2013
      integer ii1,id1,id2
      !end of added
cc      end of added
cccccccccccccccccccccccccccccc输出罗李村平均雨量
      integer  isub,nsub,iflow,nflow(nc),ig
	integer ngrid(nc,nx)       ! number of grids in this flow-interval
	integer grid_row(nc,nx,np) ! row of grids in this flow-interval
	integer grid_col(nc,nx,np) ! column of grids in this flow-interval
      common /river2    /  nflow 
      common /grid_attrib/ ngrid,grid_row,grid_col
      integer start_sub, end_sub
      integer start, finish, dt
      common /simulation/  start, finish, dt, start_sub, end_sub
      character*200  result2_dir ! directory for storing simulation result
      character*3 ch3
      integer year2000
cccccccccccccccccccccccccccccccccccccccccccccccccccc    
      real area(nrow,ncol)
      common /slope     /  area
      real evap_sum(nrow,ncol)   !added by mqh 2014/4/16
      real evap_total
      real countb 
       integer ia1,ia2
     

	  file_dir = '../data/' 
	  file_dir1 = '../data/gauge_new/'



	 do it=1, ntypes
          write(ch2, '(i2.2)') it
	    
          open(14, file='../parameter/'//'lu_class'//ch2//'.asc',
     &                    status='old')
	    read(14,*) ncols,nnc
	    read(14,*) nrows,nnr
	    read(14,*) xllcorner,xx0
	    read(14,*) yllcorner,yy0
	    read(14,*) cellsize,gridsize
	    read(14,*) nodata,znodata
c	    print *,  'nodata,znodata=' , znodata
c	    pause
          do ir=1,nnr
	    read(14,*) (land_use(ir,ic,it), ic=1,nnc)
	    end do
	  close(14)
          close(14)
        end do
          
            
         do ir=1,Maxgauge
	     do iday=1,31
              Tm(ir,iday)=-9999
              Tu(ir,iday)=-9999
              Td(ir,iday)=-9999
              Rhh(ir,iday)=-9999
              Cloud(ir,iday)=-9999
              Wd(ir,iday)=-9999
              Pr(ir,iday)=-9999
              Sun(ir,iday)=-9999
              pr1(ir,iday)=-9999
              do idh=1,24
               Prhour(ir,iday,idh)=0
              end do                        ! added by mqh, for torrential flood
		end do 
	  end do
c            added by gaobing 
            
c            end of       added 


	 call strlen(file_dir,ib1,ib2)
        print *, 'the old gauge number1', N_gauge
	  do ir=1, N_gauge
c	             print *, 'ir=,G_id=', ir,G_id(ir)
c            pause
		if(G_id(ir) .gt.100) then
           write(ch5, '(i5.5)') G_id(ir)
c	     ch9 = trim('day_'//ch5)
	     open(8,file=file_dir(ib1:ib2)//'day_'//ch5//'.dat'
     $          ,status='old')
c		 read(8,*) ch15
  	     do ic=1,50000
      read(8,*,end=177) GG_id,Y_start1,M_start1,D_start,
     $                         Tmean,Tmax,Tmin,RH1,prec,Windm,Tsun	  ! it's changed,mqh 
  

		    if(Y_start1.gt.iyear) goto 177
			if(Y_start1.eq.iyear .and. M_start1.eq.im) then 
			call alter_item(M_start1,prec,Tmean,Tmax,Tmin,
     $                    Windm,RH1,Tsun) ! it's changed,mqh
              Tm(ir,D_start)=Tmean  !M_start1,Y_start1-1950,
              Tu(ir,D_start)=Tmax
              Td(ir,D_start)=Tmin
              Rhh(ir,D_start)=RH1
              Wd(ir,D_start)=Windm
              Pr(ir,D_start)=prec
              Sun(ir,D_start)=Tsun 

      if(NNstation.eq.3) Pr1(ir,D_start)=prec

             
			end if	   
		  end do
177	      close(8)
          end if
	  end do

c       added by xujijun   2006 2 25
c       print *, 'NNstation=', NNstation         
       if(NNstation.eq.1)  goto 175
c         changged of gaobing
c         end of added
	 if(NNstation.eq.2) then
  	   N_gaugex=0
 	   do ir=1, N_gauge1 
           xtemp(ir) = x1(ir)
           ytemp(ir) = y1(ir)
           ztemp(ir) = z1(ir)
         enddo
       endif
	 if(NNstation.eq.3)  then 
	  N_gaugex=N_gauge
c	  N_gauge=10
	  
 	  do ir=1, N_gaugex 
           xtemp(ir) = x0(ir)
           ytemp(ir) = y0(ir)
           ztemp(ir) = z0(ir)
        enddo
        

	  do ir=N_gaugex+1, N_gaugex+N_gauge1 
           xtemp(ir) = x1(ir-N_gaugex)
           ytemp(ir) = y1(ir-N_gaugex)
           ztemp(ir) = z1(ir-N_gaugex)
         enddo
       endif
	 
        call strlen(file_dir1,ib1,ib2)
	  do ir=N_gaugex+1, N_gaugex+N_gauge1
		if(G_id1(ir-N_gaugex) .gt.100) then
          write(ch8, '(i8.8)') G_id1(ir-N_gaugex)
c          write(ch5, '(i5.5)') G_id1(ir-N_gaugex)
c		ch9 = trim('day_'//ch5)
	    open(18,file=file_dir1(ib1:ib2)//'day_'//ch8//'.txt'
     $         ,status='old')
c       	read(18,*) ch30
	    do ic=1,50000
c	  	read(18,*,end=178) GG_id1,Y_start1,M_start1,D_start,prec, ch15
		read(18,*,end=178) GG_id1,Y_start1,M_start1,D_start,prec
c	  	write(*,*) GG_id1,Y_start1,M_start1,D_start,prec, ch15
c            if (GG_id1.eq.41136000.and.prec.ne.0) then
c            print *,GG_id1,Y_start1,M_start1,D_start,prec
c            print *,GG_id1,Y_start1,M_start1,D_start,H_start,prec_hour
c           pause
c            endif

              if (prec > 500) prec=prec*0.1
             if(G_id1(ir-N_gaugex).ne.GG_id1) then
		     print *, 'wrong in read new rain data'  
	         pause
             endif

		   if(Y_start1.gt.iyear) goto 178
		    if(Y_start1.eq.iyear .and. M_start1.eq.im) then 
                if(prec.ge.999)  prec=0 
                 Pr1(ir,D_start)=prec
c                 Prhour(ir,D_start2,H_start)=prec2   !changed by mqh, for torrential flood
cc	flag
			end if	    
		  end do
178	      close(18)
          end if
	  end do
  
	  !added by mqh,重复了一遍之前的，分开读取小时雨量
	    do ir=N_gaugex+1, N_gaugex+N_gauge1
		if(G_id1(ir-N_gaugex) .gt.100) then
         write(ch8, '(i8.8)') G_id1(ir-N_gaugex)
          open(19,file=file_dir1(ib1:ib2)//'h'//ch8//'.dat'
     $         ,status='old')     !added by mqh,for torrential flood

	    do ic=1,100000

		read(19,*,end=1782) GG_id2,Y_start2,M_start2,D_start2,H_start,prec2 !changed by mqh,for torrential flood
            if(G_id1(ir-N_gaugex).ne.GG_id2) then
		     print *, 'wrong in read new rain data'  
	         pause
             endif

		   if(Y_start2.gt.iyear) goto 1782
			if(Y_start2.eq.iyear .and. M_start2.eq.im) then 
                   Prhour(ir,D_start2,H_start)=prec2   !changed by mqh, for torrential flood
			end if	   
c			if( M_start2.eq.10 .and. D_start2.eq.10) then
c		 print *, GG_id2,Y_start2,M_start2,D_start2,H_start,prec2,'111'
c          pause
c          endif

		  end do

c		  pause
1782        close(19)
          end if
	  end do

c      end of added

175    continue



        nnr = nrow
        nnc = ncol
	  TFF=1
        Jd = 0
	  nday = 365
	  if(mod(iyear,4).eq.0) then
	    nday=366
	    mmonth(2)=29
	   else
	    nday=365
	    mmonth(2)=28
	  end if


		
          do ir1=1,N_gauge
	       TF(ir1)=0
		   if (iyear.gt.Y_start(ir1) .or. 
     $          (iyear.eq.Y_start(ir1) .and. im.ge.M_start(ir1))) then
	        TF(ir1)=1
		    if (iyear.eq.Y_start(ir1) .and. im.eq.M_start(ir1)) then
                TFF=1
		    end if
	      end if
		  if (iyear.gt.Y_end(ir1) .or. 
     $          (iyear.eq.Y_end(ir1) .and. im.gt.M_end(ir1))) then
	        TF(ir1)=0
		  end if
	      if (iyear.eq.Y_end(ir1) .and. im.eq.M_end(ir1)) then
c	        TF(ir1)=0
	        TFF=1
		  end if
		end do
          print * ,'ok'
            
c


          TFF=1
		call in_out_T(TFF,TF,x0,y0,z0,x,y,G_elev,elev,zone,Tm,
     $        mmonth(im),mm,dd0,iyear,im,Hbeta,Tmp,N_gauge)
	    print * , 'ok'
		call in_out_T(TFF,TF,x0,y0,z0,x,y,G_elev,elev,zone,Tu,
     $        mmonth(im),mm,dd0,iyear,im,Hbeta,t_max,N_gauge)
		call in_out_T(TFF,TF,x0,y0,z0,x,y,G_elev,elev,zone,Td,
     $        mmonth(im),mm,dd0,iyear,im,Hbeta,t_min,N_gauge)
		call in_out_else(TFF,TF,x0,y0,z0,x,y,G_elev,elev,Rhh,
     $        mmonth(im),mm,dd0,iyear,im,Tmbeta,RH,N_gauge)
		call in_out_else(TFF,TF,x0,y0,z0,x,y,G_elev,elev,Cloud,
     $        mmonth(im),mm,dd0,iyear,im,Tmbeta,cld,N_gauge)
		call in_out_else(TFF,TF,x0,y0,z0,x,y,G_elev,elev,Wd,
     $        mmonth(im),mm,dd0,iyear,im,Tmbeta,wind,N_gauge)
		call in_out_else(TFF,TF,x0,y0,z0,x,y,G_elev,elev,Sun,
     $        mmonth(im),mm,dd0,iyear,im,Tmbeta,sunshine,N_gauge)

cc		call in_out_P(TFF,TF,x0,y0,z0,x,y,G_elev,elev,Pr,
cc     $        mmonth(im),mm,dd0,iyear,im,Tmbeta,pre,N_gauge)

cc       added by xujijun
          if(NNstation.eq.1) then         
            call rain_downscale2(N_gauge, N_gauge1,mmonth,iyear,    
     :                   im, Pr, Prhour)  
     
c        added by gaobing 2013	  

  	     call in_out_P(TFF,TF,x0,y0,z0,x,y,G_elev,elev,Prhour,
     $     mmonth(im),mm,dd0,iyear,im,Tmbeta,prehour,N_gauge,G_eps)
 	    endif 

	    
	    
cc       added by gaobing
c         if(im .le. 4 .or. im .ge. 10) then
c           call rain_downscale2(N_gauge, N_gauge1,mmonth,iyear,             
c     :                   im, Pr, Prhour)
c 	     call in_out_P(TFF,TF,x0,y0,z0,x,y,G_elev,elev,Prhour,
c     $        mmonth(im),mm,dd0,iyear,im,Tmbeta,prehour,N_gauge)
c	    endif 
c         end of added
cc      added by xujijun
   	    TFF=1
          Jd = 0
          Newgauge=0
         if(NNstation.eq.1)  goto 185
c
c         if(im .le. 4 .or. im .ge. 10) goto 185
c
	   if(NNstation.eq.2)  N_gaugex=0
	   if(NNstation.eq.3)  then 
	      N_gaugex=N_gauge
 	      do ir1=1, N_gaugex 
              TF1(ir1) = TF(ir1)
	        if(TF1(ir1).ge.1) Newgauge=Newgauge+1
            enddo
	   print *, 'old gauge number=', Newgauge
         endif

         do ir1=N_gaugex+1,N_gaugex+N_gauge1
            TF1(ir1) = 0
		  if (iyear.gt.Y1_start(ir1-N_gaugex) .or. 
     $           (iyear.eq.Y1_start(ir1-N_gaugex) .and. 
     $           im.ge.M1_start(ir1-N_gaugex))) then
	        TF1(ir1)=1
		    if (iyear.eq.Y1_start(ir1-N_gaugex) .and.
     $            im.eq.M1_start(ir1-N_gaugex)) then
                TFF=1
		    end if
	      end if
		  if (iyear.gt.Y1_end(ir1-N_gaugex) .or. 
     $          (iyear.eq.Y1_end(ir1-N_gaugex) .and. 
     $           im.gt.M1_end(ir1-N_gaugex))) then
	        TF1(ir1)=0
		  end if
	      if (iyear.eq.Y1_end(ir1-N_gaugex) .and. 
     $           im.eq.M1_end(ir1-N_gaugex)) then
	        TFF=1
		  end if
            if(TF1(ir1).ge.1) Newgauge=Newgauge+1
		end do
 	    print *, 'old+new gauge number=', Newgauge


		TFF=1
          if(NNstation.eq.1) goto 185
	    if(NNstation.eq.2) N_gauge2=N_gauge1
	    if(NNstation.eq.3) N_gauge2=N_gauge1+N_gaugex
                call rain_downscale2(N_gauge, N_gauge1, mmonth, iyear,  
     :                   im, Pr1, Prhour)
	   call in_out_P(TFF,TF1,xtemp,ytemp,ztemp,x,y,G_elev,elev,Prhour,
     $        mmonth(im),mm,dd0,iyear,im,Tmbeta,prehour,N_gauge2,G_eps)
	   print * ,'ok'
185      continue 

cc     end of added
c        print *, 'N_gauge, N_gaugex, N_gauge1, N_gauge2, newgauge'  
c        print *, N_gauge, N_gaugex, N_gauge1, N_gauge2, newgauge  

cc		TFF=0
		Jd=0
		if(im.gt.1) then
		  do imm= 1, im-1
	        jd = jd + mmonth(imm)
	      end do
          end if
		do 777 iday=1,mmonth(im)
	      Jd=Jd+1
            do ir = 1, nnr
              do ic = 1,nnc

			  if(t_max(ir,ic,iday).le.-100 .or. t_min(ir,ic,iday)
     $              .le.-100 .or. tmp(ir,ic,iday).le.-100) then
			     do  it=1,ntypes
	               evap(ir,ic,iday,it)=-9999
                   end do
			   else
	          SH(ir,ic,iday) =0.0
			   do  it=1,ntypes
c                latitude = lat(ir,ic)*pi/180. !degree->radian
                  evap(ir,ic,iday,it) =0.0
	            R_n=-9999
                  land_type=it

c	            if(zone(ir,ic).eq.1.and.(im.le.5 .or. im.ge.11)) then

cc	            dv=0.0
cc				if(zone(ir,ic).eq.1.and. snow0(ir,ic,it).ge.1.0) then
cc	              dv= 0.05 + 0.25 * amin1(1.0,snow0(ir,ic,it)/10.0)
cc				  if(snow0(ir,ic,it).gt.10.0) then
cc	                dv=0.0 
cc					land_type=20
cc				  end if   
cc				end if
cc				if(cld(ir,ic,iday).lt.0 .and. 
cc     :               sunshine(ir,ic,iday).lt.0) sunshine(ir,ic,iday)=0
	
     		 
				if(RH(ir,ic,iday).lt.0) RH(ir,ic,iday)=0
				if(wind(ir,ic,iday).lt.0) wind(ir,ic,iday)=0
cc      added by xujijun 2006 3 9                  
	     pretemp=0
           do idh=1,24 
             pretemp=pretemp+Prehour(ir,ic,iday,idh)
		 enddo 
           pre(ir,ic,iday)=pretemp
cc	print *, 'pre(ir,ic,iday)=', ir, ic, iday, pretemp	!flag 2006.4.16
cc      end of  added
				call cal_evp(Jd, land_type, tmp(ir,ic,iday),
     :		    t_max(ir,ic,iday),t_min(ir,ic,iday),cld(ir,ic,iday),  
     :		    sunshine(ir,ic,iday),RH(ir,ic,iday),wind(ir,ic,iday),
     :		    G_elev(ir,ic),lat(ir,ic),Ep, pre(ir,ic,iday),R_n,dv)

                 SH(ir,ic,iday)= SH(ir,ic,iday)+R_n*land_use(ir,ic,it)
c                 if(land_use(ir,ic,it).lt.0.0) land_use(ir,ic,it)=0.0                  
                  evap(ir,ic,iday,it) = Ep             !daily Ep
!           evap_sum(ir,ic)=evap_sum(ir,ic)+evap(ir,ic,iday,it)
!     :*land_use(ir,ic,it)  !added by mqh ,2014/4/16
			   end do
		      end if
			end do   
		  end do 
777     continue 




!      !shuchu evap_sum
!                znodata1 = -9999
!
!      countb=0
!      print *,'nnr,nnc=',nnr,nnc
!	     do ir=1,nnr
!		 do  ic=1,nnc
!         if(area(ir,ic).le.0) then
!         evap_sum(ir,ic)=-9999
!         else
!         countb=countb+1
!         evap_total=evap_total+evap_sum(ir,ic)/mmonth(im)
!      evap_sum(ir,ic)=evap_sum(ir,ic)/mmonth(im)
!      endif
!      enddo
!      enddo 
!      evap_total=evap_total/countb
!      
!      
!        file_dir = '../result/' 
!       call strlen(file_dir,ia1,ia2)
!      if (im.eq.1) then
!      open(012,file=file_dir(ia1:ia2)//'evap_total.dat',
!     $status='unknown')
!      else
!      open(012,file=file_dir(ia1:ia2)//'evap_total.dat',
!     $access='append',status='old')
!      endif
!      write(012,'(2i16,f8.1)') im,countb,evap_total
!      close(012) 
!      
!      
!        if(iyear.eq.2000.and.im.eq.7) then
!           file_dir = '../result/' 
!       call strlen(file_dir,ia1,ia2)
!       do it =1,1
!        	  write(ch2, '(i2.2)') it
!        open(33,file=file_dir(ia1:ia2)//'1234ep200007'//ch2//'.asc',
!     $status='unknown')
!           write(33,'(A, i16)') ncols,nnc
!	    write(33,'(A, i16)') nrows,nnr
!	    write(33,'(A, f16.3)') xllcorner,xx0
!	    write(33,'(A, f16.3)') yllcorner,yy0
!	    write(33,'(A, f16.0)') cellsize,gridsize
!	    write(33,'(A, f8.1)') nodata,znodata1
!	    write(33, '(10f11.4)') ((evap_sum(ir,ic),ic=1,nnc),ir=1,nnr)
!!	     do ir=1,nrow
!!		 write(33,'(295f13.4)') (evap_sum(ir,ic,it), ic=1,ncol)
!!        end do
!        close(33)
!        	    enddo
!	    endif
      
      !--------------------------------mqh
ccc      file_dir = '../result/' 
c       call strlen(file_dir,ia1,ia2)
c        if (iyear.eq.startyear.and.im.eq.1) then
c        open(69,file=file_dir(ia1:ia2)//'evap_sum.dat',
c     $status='unknown')
c		else
c		open(69,file=file_dir(ia1:ia2)//'evap_sum.dat',
c     $access='append',status='old')
c        endif        

c         write(69, '(i6,f16.3)') im,evap_sum
c         close(69)         
     
      




cc        added by xu , in order to test	 --start     2006 2 25
	    file_dir   = '../result/'
          nnc = ncol
	    nnr = nrow 
cc      output pre of grid
		call strlen(file_dir,ib1,ib2)
          write(ch4, '(i4.4)') iyear
          write(ch2, '(i2.2)') im
c		ch10 = trim('Pre_'//ch4//ch2)
ccc	    open(21,file=file_dir(ib1:ib2)//'pre_'//ch4//ch2//'.asc')   !mqh
c 	    do ir=1, nnr
c		  do ic=1, nnc
c			if (G_elev(ir,ic).gt.0.0) then
c		      result_m = 0.0
c			  do iday=1,mmonth(im)
c		        result_m = result_m + pre(ir,ic,iday)
c			  enddo
c			  result_month(ir,ic)= result_m/mmonth(im)
c	        else
c	          result_month(ir,ic)=-999.9
c	        endif
c		  enddo
cccc		enddo  
		


!ccccccccccc输出各子流域的面平均雨量		
c        goto 19923	
		 result2_dir="../result_spatial/"
		 do isub=start_sub,end_sub
            write(ch3, '(i3.3)') isub
            call strlen(result2_dir,id1,id2)
            year2000=0
          if (iyear.eq.2000.and.im.eq.1) then
       open(012,file=result2_dir(id1:id2)//'p'//'_'//ch3//'.dat',
     $status='unknown')
      year2000=1
        elseif (year2000.eq.0.and.iyear.ge.2000) then
	open(012,file=result2_dir(id1:id2)//'p'//'_'//ch3//'.dat',
     $access='append',status='old')
        endif
        	do iday=1,mmonth(im)	
       	      do ih=1,24	

                countt(isub)=0
                result_hour(iday,ih,isub)=0
                do iflow = 1, nflow(isub)
                     do ig = 1, ngrid(isub,iflow)           
              ir    = grid_row(isub,iflow,ig)
              ic    = grid_col(isub,iflow,ig)  
              countt(isub) = countt(isub)+1       
          result_hour(iday,ih,isub) = result_hour(iday,ih,isub) 
     $+ prehour(ir,ic,iday,ih)
                    enddo
                enddo
                
!                if(psubbasin(isub) > 2000)  then
!                do ii1 =1,4
!                   if(nbasinup(isub,ii1) > 0) then
!      result_hour(iday,ih,isub) = result_hour(iday,ih,isub)
!     $+result_hour(iday,ih,nbasinup(isub,ii1))*
!     $countt(nbasinup(isub,ii1))  
!      countt(isub) =countt(isub) +countt(nbasinup(isub,ii1))           
!                endif
!                enddo
!               endif
!               
             
      result_hour(iday,ih,isub) =
     $ result_hour(iday,ih,isub) /countt(isub)    
       write(012, '(4i6,f16.3)') iyear,im,iday,ih, 
     $result_hour(iday,ih,isub)
          enddo  !do ih
            enddo    !do iday
          		close(012) 
		enddo       ! do isub

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc		
c		print *,'end of ouputing the area_average rain'
		!end of added
		
		
		
		
          znodata1 = -999.9
	    write(21,'(a, i16)') ncols,nnc
	    write(21,'(a, i16)') nrows,nnr
	    write(21,'(a, f16.3)') xllcorner,xx0
	    write(21,'(a, f16.3)') yllcorner,yy0
	    write(21,'(a, f16.0)') cellsize,gridsize
	    write(21,'(a, f8.1)') nodata,znodata1
	write(21, '(10f11.4)') ((result_month(ir,ic),ic=1,nnc),ir=1,nnr)
	    close(21)

cc       output tmean of grid
c		ch10 = trim('Tmp_'//ch4//ch2)
cccc       open(21,file=file_dir(ib1:ib2)//'SH_'//ch4//ch2//'.asc')
cc 	    do ir=1, nnr
c		  do ic=1, nnc
c			if (G_elev(ir,ic).gt.0.0) then
c		      result_m = 0.0
c			  do iday=1,mmonth(im)
c		        result_m = result_m + SH(ir,ic,iday)
c			  enddo
c			  result_month(ir,ic)= result_m/mmonth(im)
c	        else
c	          result_month(ir,ic)=-999.9
c	        endif
c		  enddo
c		enddo  
c          znodata1 = -999.9
c	    write(21,'(a, i16)') ncols,nnc
c	    write(21,'(a, i16)') nrows,nnr
c	    write(21,'(a, f16.3)') xllcorner,xx0
c	    write(21,'(a, f16.3)') yllcorner,yy0
c	    write(21,'(a, f16.0)') cellsize,gridsize
c	    write(21,'(a, f8.1)') nodata,znodata1
c	write(21, '(10f11.4)') ((result_month(ir,ic),ic=1,nnc),ir=1,nnr)
c	    close(21)

        goto 2119
!	    open(21,file=file_dir(ib1:ib2)//'tmp_'//ch4//ch2//'.asc')
! 	    do ir=1, nnr
!		  do ic=1, nnc
!			if (G_elev(ir,ic).gt.0.0) then
!		      result_m = 0.0
!			  do iday=1,mmonth(im)
!		        result_m = result_m + Tmp(ir,ic,iday)
!			  enddo
!			  result_month(ir,ic)= result_m/mmonth(im)
!	        else
!	          result_month(ir,ic)=-999.9
!	        endif
!		  enddo
!		enddo  
!	    write(21,'(a, i16)') ncols,nnc
!	    write(21,'(a, i16)') nrows,nnr
!	    write(21,'(a, f16.3)') xllcorner,xx0
!	    write(21,'(a, f16.3)') yllcorner,yy0
!	    write(21,'(a, f16.0)') cellsize,gridsize
!	    write(21,'(a, f8.1)') nodata,znodata1
!	write(21, '(10f11.4)') ((result_month(ir,ic),ic=1,nnc),ir=1,nnr)
!	    close(21)

c		ch10 = trim('Tmn_'//ch4//ch2)
!	    open(21,file=file_dir(ib1:ib2)//'Tmn_'//ch4//ch2//'.asc')
! 	    do ir=1, nnr
!		  do ic=1, nnc
!			if (G_elev(ir,ic).gt.0.0) then
!		      result_m = 0.0
!			  do iday=1,mmonth(im)
!		        result_m = result_m + T_min(ir,ic,iday)
!			  enddo
!			  result_month(ir,ic)= result_m/mmonth(im)
!	        else
!	          result_month(ir,ic)=-999.9
!	        endif
!		  enddo
!		enddo  
!	    write(21,'(a, i16)') ncols,nnc
!	    write(21,'(a, i16)') nrows,nnr
!	    write(21,'(a, f16.3)') xllcorner,xx0
!	    write(21,'(a, f16.3)') yllcorner,yy0
!	    write(21,'(a, f16.0)') cellsize,gridsize
!	    write(21,'(a, f8.1)') nodata,znodata1
! 	write(21, '(10f11.4)') ((result_month(ir,ic),ic=1,nnc),ir=1,nnr)
!	    close(21)

c		ch10 = trim('Tmx_'//ch4//ch2)
!	    open(21,file=file_dir(ib1:ib2)//'Tmx_'//ch4//ch2//'.asc')
! 	    do ir=1, nnr
!		  do ic=1, nnc
!			if (G_elev(ir,ic).gt.0.0) then
!		      result_m = 0.0
!			  do iday=1,mmonth(im)
!		        result_m = result_m + T_max(ir,ic,iday)
!			  enddo
!			  result_month(ir,ic)= result_m/mmonth(im)
!	        else
!	          result_month(ir,ic)=-999.9
!	        endif
!		  enddo
!		enddo  
!	    write(21,'(a, i16)') ncols,nnc
!	    write(21,'(a, i16)') nrows,nnr
!	    write(21,'(a, f16.3)') xllcorner,xx0
!	    write(21,'(a, f16.3)') yllcorner,yy0
!	    write(21,'(a, f16.0)') cellsize,gridsize
!	    write(21,'(a, f8.1)') nodata,znodata1
! 	write(21, '(10f11.4)') ((result_month(ir,ic),ic=1,nnc),ir=1,nnr)
!	    close(21)

cc       output  Wind of grid 
c		ch10 = trim('Wnd_'//ch4//ch2)
!	    open(21,file=file_dir(ib1:ib2)//'Wnd_'//ch4//ch2//'.asc')
! 	    do ir=1, nnr
!		  do ic=1, nnc
!			if (G_elev(ir,ic).gt.0.0) then
!		      result_m = 0.0
!			  do iday=1,mmonth(im)
!		        result_m = result_m + Wind(ir,ic,iday)
!			  enddo
!			  result_month(ir,ic)= result_m/mmonth(im)
!	        else
!	          result_month(ir,ic)=-999.9
!	        endif
!		  enddo
!		enddo  
!	    write(21,'(a, i16)') ncols,nnc
!	    write(21,'(a, i16)') nrows,nnr
!	    write(21,'(a, f16.3)') xllcorner,xx0
!	    write(21,'(a, f16.3)') yllcorner,yy0
!	    write(21,'(a, f16.0)') cellsize,gridsize
!	    write(21,'(a, f8.1)') nodata,znodata1
!	write(21, '(10f11.4)') ((result_month(ir,ic),ic=1,nnc),ir=1,nnr)
!	    close(21)

cc       output  RH of grid 
c		ch10 = trim('RH_'//ch4//ch2)
!	    open(21,file=file_dir(ib1:ib2)//'RH_'//ch4//ch2//'.asc')
! 	    do ir=1, nnr
!		  do ic=1, nnc
!			if (G_elev(ir,ic).gt.0.0) then
!		      result_m = 0.0
!			  do iday=1,mmonth(im)
!		        result_m = result_m + RH(ir,ic,iday)
!			  enddo
!			  result_month(ir,ic)= result_m/mmonth(im)
!	        else
!	          result_month(ir,ic)=-999.9
!	        endif
!		  enddo
!		enddo  
!	    write(21,'(a, i16)') ncols,nnc
!	    write(21,'(a, i16)') nrows,nnr
!	    write(21,'(a, f16.3)') xllcorner,xx0
!	    write(21,'(a, f16.3)') yllcorner,yy0
!	    write(21,'(a, f16.0)') cellsize,gridsize
!	    write(21,'(a, f8.1)') nodata,znodata1
!	write(21, '(10f11.4)') ((result_month(ir,ic),ic=1,nnc),ir=1,nnr)
!	    close(21)
!
!cc       output  cld of grid 
!c		ch10 = trim('Cld_'//ch4//ch2)
!	    open(21,file=file_dir(ib1:ib2)//'Cld_'//ch4//ch2//'.asc')
! 	    do ir=1, nnr
!		  do ic=1, nnc
!			if (G_elev(ir,ic).gt.0.0) then
!		      result_m = 0.0
!			  do iday=1,mmonth(im)
!		        result_m = result_m + cld(ir,ic,iday)
!			  enddo
!			  result_month(ir,ic)= result_m/mmonth(im)
!	        else
!	          result_month(ir,ic)=-999.9
!	        endif
!		  enddo
!		enddo  
!	    write(21,'(a, i16)') ncols,nnc
!	    write(21,'(a, i16)') nrows,nnr
!	    write(21,'(a, f16.3)') xllcorner,xx0
!	    write(21,'(a, f16.3)') yllcorner,yy0
!	    write(21,'(a, f16.0)') cellsize,gridsize
!	    write(21,'(a, f8.1)') nodata,znodata1
!	write(21, '(10f11.4)') ((result_month(ir,ic),ic=1,nnc),ir=1,nnr)
!	    close(21)
!
!cc       output  sunshine of grid 
!c		ch10 = trim('Sun_'//ch4//ch2)
!	    open(21,file=file_dir(ib1:ib2)//'Sun_'//ch4//ch2//'.asc')
! 	    do ir=1, nnr
!		  do ic=1, nnc
!			if (G_elev(ir,ic).gt.0.0) then
!		      result_m = 0.0
!			  do iday=1,mmonth(im)
!		        result_m = result_m + sunshine(ir,ic,iday)
!			  enddo
!			  result_month(ir,ic)= result_m/mmonth(im)
!	        else
!	          result_month(ir,ic)=-999.9
!	        endif
!		  enddo
!		enddo  
!	    write(21,'(a, i16)') ncols,nnc
!	    write(21,'(a, i16)') nrows,nnr
!	    write(21,'(a, f16.3)') xllcorner,xx0
!	    write(21,'(a, f16.3)') yllcorner,yy0
!	    write(21,'(a, f16.0)') cellsize,gridsize
!	    write(21,'(a, f8.1)') nodata,znodata1
!	write(21, '(10f11.4)') ((result_month(ir,ic),ic=1,nnc),ir=1,nnr)
!	    close(21)
2119   continue 
cc        added by xu , in order to test ----end
    	end

c*******************************************************************************
c	To interpolate meteorological gauge data to 10km grid
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*
c 

c      interplate from eight nearest points
c
       subroutine in_out_T(TFF,TF,x0,y0,z0,x,y,G_elev,elev,zone,datain,
     $                 mday,mm,dd0,iyear,im,Hbeta,dataout,N_gauge)
       include 'dims2.inc'

        real x0(Maxgauge),y0(Maxgauge),z0(Maxgauge),elev(Maxgauge)
        real x(nrow,ncol),y(nrow,ncol),G_elev(nrow,ncol),zone(nrow,ncol)
        real datain(Maxgauge,31)
        real dataout(nrow,ncol,31)
        real xx(nrow,ncol,8),yy(nrow,ncol,8)
 	    real ww(nrow,ncol,10),dd(nrow,ncol,8),beta(nrow,ncol,4)
        real Hbeta(5,12,4)
        integer zz(nrow,ncol,8)
        integer TF(Maxgauge)
        integer irr1, irr2, irr3, irr4
        integer irr5, irr6, irr7, irr8
	    real d1, d2, d3, d4, d5, d6, d7, d8,dd0
	    real wk(8),ak(8),sta(8) 
c       double precision wk,ak,sta
	    integer day,im,mm
	    integer N_gauge,TFF
	    integer mday,iyear
c       dimension xp(5000000),yp(5000000),a(6),zp(5000000)
c       double precision xp,yp,a,dt1,dt2,dt3,b,zp

	  If(TFF.eq.1) then
	  do ir = 1, nrow	  
	    do ic = 1, ncol
c################################################################
c		To find out the 8 nearest points
c################################################################
	   if(G_elev(ir, ic).gt.0.0) then
		 do ir1=1,N_gauge
		  z0(ir1)=0
cc		  if(x0(ir1) .ge. (x(ir,ic)-1500000.0) .and.
cc     :      	 x0(ir1) .le. (x(ir,ic)+1500000.0) .and.
cc     :		 y0(ir1) .ge. (y(ir,ic)-1500000.0) .and.
cc     :       y0(ir1) .le. (y(ir,ic)+1500000.0) .and. TF(ir1).eq.1) then
cc	         z0(ir1)=1
		  if(x0(ir1) .ge. (x(ir,ic)-800000.0) .and.	! 修改于2007.9.8
     :      	 x0(ir1) .le. (x(ir,ic)+800000.0) .and.
     :		 y0(ir1) .ge. (y(ir,ic)-800000.0) .and.
     :       y0(ir1) .le. (y(ir,ic)+800000.0) .and. TF(ir1).eq.1) then
	         z0(ir1)=1
		  end if
	    end do 
c		  
		d1 = 10.0E30
		d2 = 10.0E31
		d3 = 10.0E32
		d4 = 10.0E33
	    d5 = 10.0E34
	    d6 = 10.0E35
	    d7 = 10.0E36
	    d8 = 10.0E37
c	    print *,N_gauge
c	    pause

	    do ir1 = 1, N_gauge

		  if(z0(ir1) .eq. 1) then
            d = (x(ir,ic)-x0(ir1))**2+
     :          (y(ir,ic)-y0(ir1))**2
	      if(d .lt. d1) then
	         d8=d7
	         irr8=irr7
	         d7=d6
	         irr7=irr6
	         d6=d5
	         irr6=irr5
	         d5=d4
	         irr5=irr4
	         d4=d3
	         irr4=irr3
	         d3=d2
	         irr3=irr2
	         d2=d1
	         irr2=irr1
	         d1=d
	         irr1=ir1
	         goto 199
	       end if 
	       if(d .ge. d1 .and. d .lt. d2) then
	         d8=d7
	         irr8=irr7
	         d7=d6
	         irr7=irr6
	         d6=d5
	         irr6=irr5
	         d5=d4
	         irr5=irr4
	         d4=d3
	         irr4=irr3
	         d3=d2
	         irr3=irr2
	         d2=d
	         irr2=ir1
	          goto 199
		   end if 
	       if(d .ge. d2 .and. d .lt. d3) then
	         d8=d7
	         irr8=irr7
	         d7=d6
	         irr7=irr6
	         d6=d5
	         irr6=irr5
	         d5=d4
	         irr5=irr4
	         d4=d3
	         irr4=irr3
	         d3=d
	         irr3=ir1
	          goto 199
		   end if 
	       if(d .ge. d3 .and. d .lt. d4) then
	         d8=d7
	         irr8=irr7
	         d7=d6
	         irr7=irr6
	         d6=d5
	         irr6=irr5
	         d5=d4
	         irr5=irr4
	         d4=d
	         irr4=ir1
	          goto 199
		   end if 
	       if(d .ge. d4 .and. d .lt. d5) then
	         d8=d7
	         irr8=irr7
	         d7=d6
	         irr7=irr6
	         d6=d5
	         irr6=irr5
	         d5=d
	         irr5=ir1
	          goto 199
		   end if 
	       if(d .ge. d5 .and. d .lt. d6) then
	         d8=d7
	         irr8=irr7
	         d7=d6
	         irr7=irr6
	         d6=d
	         irr6=ir1
	          goto 199
		   end if 
	       if(d .ge. d6 .and. d .lt. d7) then
	         d8=d7
	         irr8=irr7
	         d7=d
	         irr7=ir1
	          goto 199
		   end if 
	       if(d .ge. d7 .and. d .lt. d8) then
	         d8=d
	         irr8=ir1
	          goto 199
		   end if 
		   end if
199	     continue
	    end do
        !标记
	    xx(ir,ic,1)=x0(irr1)
	    xx(ir,ic,2)=x0(irr2)
	    xx(ir,ic,3)=x0(irr3)
	    xx(ir,ic,4)=x0(irr4)
	    xx(ir,ic,5)=x0(irr5)
	    xx(ir,ic,6)=x0(irr6)
	    xx(ir,ic,7)=x0(irr7)
	    xx(ir,ic,8)=x0(irr8)
	    yy(ir,ic,1)=y0(irr1)
	    yy(ir,ic,2)=y0(irr2)
	    yy(ir,ic,3)=y0(irr3)
	    yy(ir,ic,4)=y0(irr4)
	    yy(ir,ic,5)=y0(irr5)
	    yy(ir,ic,6)=y0(irr6)
	    yy(ir,ic,7)=y0(irr7)
	    yy(ir,ic,8)=y0(irr8)
	    zz(ir,ic,1) =irr1   
	    zz(ir,ic,2) =irr2   
	    zz(ir,ic,3) =irr3   
	    zz(ir,ic,4) =irr4   
	    zz(ir,ic,5) =irr5   
	    zz(ir,ic,6) =irr6   
	    zz(ir,ic,7) =irr7   
	    zz(ir,ic,8) =irr8
c		print *, irr1, irr2, irr3, irr4,irr5, irr6, irr7, irr8,ir,ic,x(ir,ic)
c		print *, d1, d2, d3, d4,d5, d6, d7, d8
c		pause
		dd(ir,ic,1)=sqrt(d1)
	    dd(ir,ic,2)=sqrt(d2)
	    dd(ir,ic,3)=sqrt(d3)
	    dd(ir,ic,4)=sqrt(d4)
	    dd(ir,ic,5)=sqrt(d5)
	    dd(ir,ic,6)=sqrt(d6)
	    dd(ir,ic,7)=sqrt(d7)
	    dd(ir,ic,8)=sqrt(d8)
	    
		mm=4
cc	    dd0=1200.0*1000.0
	    dd0=800.0*1000.0	! 修改于2006.4.26
		do ir1=1, 8
		  wk(ir1)=exp(-mm*dd(ir,ic,ir1)/dd0)
c		  Wk_tmp=-exp(-mm*1.0)
c	      Wk_d=dd(ir,ic,ir1)/dd0
c	      wk(ir1)=exp(-1*mm*Wk_d*Wk_d)+Wk_tmp
c            if(dd(ir,ic,ir1).gt.300.0*1000.0 .and. iyear.lt.1959
c     $                             .and.G_elev(ir, ic).gt.2000.0) then
            if(dd(ir,ic,ir1).gt.10000.0*10000.0) then   !mqh
cc            if(dd(ir,ic,ir1).gt.400.0*1000.0) then
              wk(ir1)=0.0
	      end if 
          end do

		do ir1=1,8
	       ak(ir1)=0.0
	       Sum_wl=0.0
		   do ir2=1,8
	         if(ir2 .ne. ir1) then
                 if(y0(zz(ir,ic,ir2)).ge.y(ir,ic)) then
	             sign1=1
				else
	             sign1=-1
			   end if
                 if(y0(zz(ir,ic,ir1)).ge.y(ir,ic)) then
	             sign2=1
				else
	             sign2=-1
			   end if
c			   print *,x0(ir2),x(ir,ic),dd(ir,ic,ir2)
c			   print *,x0(ir1),x(ir,ic),dd(ir,ic,ir1)
			   sta1=(x0(zz(ir,ic,ir2))-x(ir,ic))/dd(ir,ic,ir2)
			   sta2=(x0(zz(ir,ic,ir1))-x(ir,ic))/dd(ir,ic,ir1)	
			   sta(ir2)=sign1*acos(sta1)
     $                   -sign2*acos(sta2)
		       ak(ir1)=ak(ir1)+wk(ir2)*(1-cos(sta(ir2)))
	           Sum_wl=Sum_wl+wk(ir2)
			 end if  
		   end do
	       if(Sum_wl.gt.0) then
		     ak(ir1)=ak(ir1)/Sum_wl
	        else
		     ak(ir1)=0.0
	       end if
          end do
		
	    Sum_wk=0.0
		do ir1=1,8
		   wk(ir1)=wk(ir1)*(1+ak(ir1))
             Sum_wk=Sum_wk+wk(ir1)
             ww(ir,ic,ir1)=wk(ir1)
          end do
		ww(ir,ic,9)=Sum_wk
	 
	    end if
		end do
	   end do
	  end if 
     
	  do ir = 1, nrow
	    do ic = 1, ncol
          if(G_elev(ir, ic).gt.0.0) then
cc            if(zone(ir,ic).eq.1 .and. iyear.lt.1959) then
c              beta(ir,ic,1)=Hbeta(zone(ir,ic),im,1)
c              beta(ir,ic,2)=Hbeta(zone(ir,ic),im,2)
cc              beta(ir,ic,1)=0
cc              beta(ir,ic,2)=-1.0/1000.0
cc     $                    -2.5/1000.0*sin(3.14159263*(im-1)/12.0)
cc            else 
cc              beta(ir,ic,1)=0
cc              beta(ir,ic,2)=-1.5/1000.0-int(zone(ir,ic)/4)/1000.0
cc     $                    -4.0/1000.0*sin(3.14159263*(im-1)/12.0)

              beta(ir,ic,1)=0
c              beta(ir,ic,2)=-5.0/1000.0-int(zone(ir,ic)/4)/1000.0
	        beta(ir,ic,2)=-5.0/1000.0-int(zone(ir,ic)/4)/1000.0
c               beta(ir,ic,2)=0
cc     $                    -2.0/1000.0*sin(3.14159263*(im-1)/12.0)



cc           end if
       
	      beta(ir,ic,3)=Hbeta(zone(ir,ic),im,3)
            beta(ir,ic,4)=Hbeta(zone(ir,ic),im,4)
		end if
		end do
        end do        
c
	  do day=1,mday
c	   print *, day
	  do ir = 1, nrow
	    do ic = 1, ncol
		dataout(ir,ic,day)=-9999
	    if(G_elev(ir, ic).gt.0.0) then
c
c################################################################
c	To interplate the value using the inverse-distance weight
c################################################################
c         Sum_pr=0.0  
c		Sum_wk=ww(ir,ic,9)
c		do ir1=1,8
c		   if(datain(zz(ir,ic,ir1),iyear-1950,im,day).gt.0) then
c	         Sum_pr=Sum_pr+ww(ir,ic,ir1)*1.0
c			else
c	         if(datain(zz(ir,ic,ir1),iyear-1950,im,day).lt.0) then
c	           Sum_wk=Sum_wk-ww(ir,ic,ir1)
c	          else
c	           Sum_pr=Sum_pr+ww(ir,ic,ir1)*0.0
c	         end if
c		   end if
c          end do
c          if(Sum_wk.gt.0.0001) then
c		  ww(ir,ic,10)=Sum_pr/Sum_wk
c           else
c		  ww(ir,ic,10)=0.0
c	    end if
c
		  dataout(ir,ic,day)=0.0
		  Sum_wk=ww(ir,ic,9)
		  do ir1=1, 8
			if(datain(zz(ir,ic,ir1),day).gt.-80.0) then
                dz=G_elev(ir, ic)-elev(zz(ir,ic,ir1))
c                dy=(y(ir, ic)-y0(zz(ir,ic,ir1)))/10000.0
			  dataout(ir,ic,day)=dataout(ir,ic,day)+ww(ir,ic,ir1)
     $                    *(datain(zz(ir,ic,ir1),day)
     $                             +beta(ir,ic,1)+beta(ir,ic,2)*dz)
	 

c     $                             +beta(ir,ic,3)+beta(ir,ic,4)*dy)
	         else
			  Sum_wk=Sum_wk-ww(ir,ic,ir1)  
	        end if
	      end do
            if(Sum_wk.gt.ww(ir,ic,9)/8.0) then
		    dataout(ir,ic,day)=dataout(ir,ic,day)/Sum_wk
             else
		    dataout(ir,ic,day)=-9999
	        d_temp=0.0
	        n_temp=0
			do ir1=8,1,-1
			 if(datain(zz(ir,ic,ir1),day).gt.-80.0) then 
                 dz=G_elev(ir, ic)-elev(zz(ir,ic,ir1))
                 dy=(y(ir, ic)-y0(zz(ir,ic,ir1)))/10000.0
		       n_temp=n_temp+1
			   d_temp=d_temp+datain(zz(ir,ic,ir1),day)
     $                        +beta(ir,ic,1)+beta(ir,ic,2)*dz
c     $                        +beta(ir,ic,3)+beta(ir,ic,4)*dy
               end if 
	        end do
	        if(n_temp.gt.0) then
			  dataout(ir,ic,day)=d_temp/n_temp
	        end if
		  end if

	    end if
		end do
	  end do

	  end do

	end 
c
c*******************************************************************************
c	To interpolate meteorological gauge data to 10km grid
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*
c	To extract rainfall for the Yellow River basin from the gauge data set
c     interplate from eight nearest points
c
      subroutine in_out_P(TFF,TF,xp0,yp0,zp0,x,y,G_elev,elev,datain1,
     $                mday,mm,dd0,iyear,im,Hbeta,dataout1,NN_gauge,
     :            Gauge_eps)
       include 'dims2.inc'
	real xp0(Maxgauge),yp0(Maxgauge),zp0(Maxgauge),elev(Maxgauge)
	real x(nrow,ncol),y(nrow,ncol),G_elev(nrow,ncol)
	real datain1(Maxgauge,31,24)
	real dataout1(nrow,ncol,31,24)
	real xx(nrow,ncol,8),yy(nrow,ncol,8)
	real ww(nrow,ncol,10),dd(nrow,ncol,8),beta(nrow,ncol,4)
	real Hbeta(12,4),dz
	integer zz(nrow,ncol,8)
	integer TF(Maxgauge)
c
	integer irr1, irr2, irr3, irr4
	integer irr5, irr6, irr7, irr8
	real d1, d2, d3, d4, d5, d6, d7, d8,dd0
	real wk(8),ak(8),sta(8) 
c      double precision wk,ak,sta

	integer day,im,mm,ihour
	integer NN_gauge,TFF
	integer mday,iyear
      real Gauge_eps(Maxgauge)

c     
	  If(TFF.eq.1) then
	  do ir = 1, nrow
	    do ic = 1, ncol
	     if(G_elev(ir, ic).gt.0.0) then    ! add by xujijun 20060324

c################################################################
c		To find out the 8 nearest points
c################################################################
		do ir1=1,NN_gauge
		  zp0(ir1)=0
		  if(xp0(ir1) .ge. (x(ir,ic)-1000000.0) .and.
     :      	 xp0(ir1) .le. (x(ir,ic)+1000000.0) .and.
     :		 yp0(ir1) .ge. (y(ir,ic)-1000000.0) .and.
     :       yp0(ir1) .le. (y(ir,ic)+1000000.0) .and. TF(ir1).eq.1) then
	         zp0(ir1)=1
		  end if
	    end do 
c		  
		d1 = 10.0E30
		d2 = 10.0E31
		d3 = 10.0E32
		d4 = 10.0E33
	    d5 = 10.0E34
	    d6 = 10.0E35
	    d7 = 10.0E36
	    d8 = 10.0E37
	    do ir1 = 1, NN_gauge
		  if(zp0(ir1) .eq. 1) then
            d = (x(ir,ic)-xp0(ir1))**2+
     :          (y(ir,ic)-yp0(ir1))**2
	      if(d .lt. d1) then
	         d8=d7
	         irr8=irr7
	         d7=d6
	         irr7=irr6
	         d6=d5
	         irr6=irr5
	         d5=d4
	         irr5=irr4
	         d4=d3
	         irr4=irr3
	         d3=d2
	         irr3=irr2
	         d2=d1
	         irr2=irr1
	         d1=d
	         irr1=ir1
	         goto 199
	       end if 
	       if(d .ge. d1 .and. d .lt. d2) then
	         d8=d7
	         irr8=irr7
	         d7=d6
	         irr7=irr6
	         d6=d5
	         irr6=irr5
	         d5=d4
	         irr5=irr4
	         d4=d3
	         irr4=irr3
	         d3=d2
	         irr3=irr2
	         d2=d
	         irr2=ir1
	          goto 199
		   end if 
	       if(d .ge. d2 .and. d .lt. d3) then
	         d8=d7
	         irr8=irr7
	         d7=d6
	         irr7=irr6
	         d6=d5
	         irr6=irr5
	         d5=d4
	         irr5=irr4
	         d4=d3
	         irr4=irr3
	         d3=d
	         irr3=ir1
	          goto 199
		   end if 
	       if(d .ge. d3 .and. d .lt. d4) then
	         d8=d7
	         irr8=irr7
	         d7=d6
	         irr7=irr6
	         d6=d5
	         irr6=irr5
	         d5=d4
	         irr5=irr4
	         d4=d
	         irr4=ir1
	          goto 199
		   end if 
	       if(d .ge. d4 .and. d .lt. d5) then
	         d8=d7
	         irr8=irr7
	         d7=d6
	         irr7=irr6
	         d6=d5
	         irr6=irr5
	         d5=d
	         irr5=ir1
	          goto 199
		   end if 
	       if(d .ge. d5 .and. d .lt. d6) then
	         d8=d7
	         irr8=irr7
	         d7=d6
	         irr7=irr6
	         d6=d
	         irr6=ir1
	          goto 199
		   end if 
	       if(d .ge. d6 .and. d .lt. d7) then
	         d8=d7
	         irr8=irr7
	         d7=d
	         irr7=ir1
	          goto 199
		   end if 
	       if(d .ge. d7 .and. d .lt. d8) then
	         d8=d
	         irr8=ir1
	          goto 199
		   end if 
		   end if
199	     continue
	    end do

	    xx(ir,ic,1)=xp0(irr1)
	    xx(ir,ic,2)=xp0(irr2)
	    xx(ir,ic,3)=xp0(irr3)
	    xx(ir,ic,4)=xp0(irr4)
	    xx(ir,ic,5)=xp0(irr5)
	    xx(ir,ic,6)=xp0(irr6)
	    xx(ir,ic,7)=xp0(irr7)
	    xx(ir,ic,8)=xp0(irr8)
	    yy(ir,ic,1)=yp0(irr1)
	    yy(ir,ic,2)=yp0(irr2)
	    yy(ir,ic,3)=yp0(irr3)
	    yy(ir,ic,4)=yp0(irr4)
	    yy(ir,ic,5)=yp0(irr5)
	    yy(ir,ic,6)=yp0(irr6)
	    yy(ir,ic,7)=yp0(irr7)
	    yy(ir,ic,8)=yp0(irr8)
	    zz(ir,ic,1) =irr1   
	    zz(ir,ic,2) =irr2   
	    zz(ir,ic,3) =irr3   
	    zz(ir,ic,4) =irr4   
	    zz(ir,ic,5) =irr5   
	    zz(ir,ic,6) =irr6   
	    zz(ir,ic,7) =irr7   
	    zz(ir,ic,8) =irr8
c		print *, irr1, irr2, irr3, irr4,irr5, irr6, irr7, irr8,ir,ic,x(ir,ic)
c		print *, d1, d2, d3, d4,d5, d6, d7, d8
c        pause

cc	    if(ir.eq.24.and.ic.eq.76) then
cc		 print *, irr1, irr2, irr3, irr4,irr5, irr6, irr7, irr8
cc		 print *, d1, d2, d3, d4,d5, d6, d7, d8
cc         pause
cc	    endif

cc	    if(ir.eq.79.and.ic.eq.86) then
cc		 print *, irr1, irr2, irr3, irr4,irr5, irr6, irr7, irr8
cc		 print *, d1, d2, d3, d4,d5, d6, d7, d8
cc         pause
cc	    endif

		dd(ir,ic,1)=sqrt(d1)
	    dd(ir,ic,2)=sqrt(d2)
	    dd(ir,ic,3)=sqrt(d3)
	    dd(ir,ic,4)=sqrt(d4)
	    dd(ir,ic,5)=sqrt(d5)
	    dd(ir,ic,6)=sqrt(d6)
	    dd(ir,ic,7)=sqrt(d7)
	    dd(ir,ic,8)=sqrt(d8)
	    
		mm=4
          dd0=400.0*1000.0
c	    dd0=40.0*1000.0
		do ir1=1, 8
		   wk(ir1)=exp(-mm*dd(ir,ic,ir1)/dd0)
c		   Wk_tmp=-exp(-mm*1.0)
c		   Wk_d=dd(ir,ic,ir1)/dd0
c		   wk(ir1)=exp(-1*mm*Wk_d*Wk_d)+Wk_tmp
c	       if(dd(ir,ic,ir1).gt.dd0) then
c	         wk(ir1)=0.0
c		   end if 
            if(dd(ir,ic,ir1).gt.60.0*1000.0) then
              wk(ir1)=0.0
	      end if 
          end do
		do ir1=1,8
	       ak(ir1)=0.0
	       Sum_wl=0.0
		   do ir2=1,8
	         if(ir2 .ne. ir1) then
                 if(yp0(zz(ir,ic,ir2)).ge.y(ir,ic)) then
	             sign1=1
				else
	             sign1=-1
			   end if
                 if(yp0(zz(ir,ic,ir1)).ge.y(ir,ic)) then
	             sign2=1
				else
	             sign2=-1
			   end if
c			   print *,xp0(ir2),x(ir,ic),dd(ir,ic,ir2)
c			   print *,xp0(ir1),x(ir,ic),dd(ir,ic,ir1)
			   sta1=(xp0(zz(ir,ic,ir2))-x(ir,ic))/dd(ir,ic,ir2)
			   sta2=(xp0(zz(ir,ic,ir1))-x(ir,ic))/dd(ir,ic,ir1)	
			   sta(ir2)=sign1*acos(sta1)
     $                   -sign2*acos(sta2)
		       ak(ir1)=ak(ir1)+wk(ir2)*(1-cos(sta(ir2)))
	           Sum_wl=Sum_wl+wk(ir2)
			 end if  
		   end do
	       if(Sum_wl.gt.0) then
		     ak(ir1)=ak(ir1)/Sum_wl
	        else
		     ak(ir1)=0.0
	       end if
          end do
		
	    Sum_wk=0.0
		do ir1=1,8
		   wk(ir1)=wk(ir1)*(1+ak(ir1))
             Sum_wk=Sum_wk+wk(ir1)
             ww(ir,ic,ir1)=wk(ir1)
          end do
		ww(ir,ic,9)=Sum_wk
	    
	      endif
	    end do
	  end do
	  end if 
c
c
	  do day=1,mday
         do ihour = 1, 24 
	    do ir = 1, nrow
	     do ic = 1, ncol
	      if(G_elev(ir, ic).gt.0.0) then    ! add by xujijun 20060324

c################################################################
c	To interplate the value using the inverse-distance weight
c################################################################
		  dataout1(ir,ic,day,ihour)=0.0 
		  Sum_wk=ww(ir,ic,9)
		  do ir1=1, 8
	         
			if(datain1(zz(ir,ic,ir1),day,ihour).ge.0.0) then
c	         dz=G_elev(ir, ic)-elev(zz(ir,ic,ir1))	     
               dz=0.0         
c	 if(npd(zz(ir,ic,ir1))> 0 .and. nh(zz(ir,ic,ir1),day) > 0) then
c	       print *, npd(zz(ir,ic,ir1)),nh(zz(ir,ic,ir1),day)
             dataout1(ir,ic,day,ihour)=dataout1(ir,ic,day,ihour)+
     $    ww(ir,ic,ir1)*datain1(zz(ir,ic,ir1),day,ihour)*
     :  (1+Gauge_eps(zz(ir,ic,ir1))*dz/10000.0)
c     :dz*11.9/100.0*p_mm(im) /npd(zz(ir,ic,ir1))/1.0/
c     : nh(zz(ir,ic,ir1),day)/1.0)
c	    else
c              dataout1(ir,ic,day,ihour)=dataout1(ir,ic,day,ihour)+
c     $    ww(ir,ic,ir1)*datain1(zz(ir,ic,ir1),day,ihour)
c	endif
	               
c    +     
c     $  dz*11.9/100.0/24.0/365.0)
	         else
			  Sum_wk=Sum_wk-ww(ir,ic,ir1)  
	        end if
	      end do
            if(Sum_wk.gt.ww(ir,ic,9)/8.0) then
		    dataout1(ir,ic,day,ihour)=dataout1(ir,ic,day,ihour)/Sum_wk

cc	           if(ir.eq.20.and.ic.eq.24) then
cc                 print *, 'ir,ic, pre', ir, ic,dataout1(ir,ic,day,ihour)
cc                 pause
cc	           endif
cc	           if(ir.eq.79.and.ic.eq.86) then
cc                 print *, 'ir,ic, pre', ir, ic,dataout1(ir,ic,day,ihour)
cc                 pause
cc	           endif

             else
		    dataout1(ir,ic,day,ihour)=-9999
	        d_temp=0.0
	        n_temp=0
			do ir1=8,1,-1
			 if(datain1(zz(ir,ic,ir1),day,ihour).ge.0.0) then 
		       n_temp=n_temp+1
c	           dz=G_elev(ir, ic)-elev(zz(ir,ic,ir1))
                 dz=0.0
c      if(npd(zz(ir,ic,ir1))> 0 .and. nh(zz(ir,ic,ir1),day)>0) then
	d_temp=d_temp+datain1(zz(ir,ic,ir1),day,ihour)*
     :  (1+Gauge_eps(zz(ir,ic,ir1))*dz/10000.0)
c     :dz*11.9/100.0*p_mm(im) /npd(zz(ir,ic,ir1))/1.0/
c     : nh(zz(ir,ic,ir1),day)/1.0
	 
c	          else
c         d_temp=d_temp+datain1(zz(ir,ic,ir1),day,ihour)
c	           endif
c	  +
c     :	   dz*11.9/100.0/24.0/365.0
               end if
	        end do
	        if(n_temp.gt.0) then
	          dataout1(ir,ic,day,ihour)=d_temp/n_temp
			end if
		  end if
c                if(G_elev(ir, ic).gt. 3000) then
c          dataout1(ir,ic,day,ihour)=dataout1(ir,ic,day,ihour)
c     :        *(1+(G_elev(ir, ic)-1000)*0.04/1000)
c	           endif
	      endif
	    end do
	   end do
	  enddo
	 end do

	end

c*******************************************************************************
c	To interpolate meteorological gauge data to 10km grid
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*
c	To extract rainfall for the Yellow River basin from the gauge data set
c     interplate from eight nearest points
c
      subroutine in_out_else(TFF,TF,x0,y0,z0,x,y,G_elev,elev,
     $           datain,mday,mm,dd0,iyear,im,Hbeta,dataout,N_gauge)
       include 'dims2.inc'

	real x0(Maxgauge),y0(Maxgauge),z0(Maxgauge),elev(Maxgauge)
	real x(nrow,ncol),y(nrow,ncol),G_elev(nrow,ncol)
	real datain(Maxgauge,31)
	real dataout(nrow,ncol,31)
	real xx(nrow,ncol,8),yy(nrow,ncol,8)
	real ww(nrow,ncol,10),dd(nrow,ncol,8),beta(nrow,ncol,4)
	real Hbeta(12,4)
	integer zz(nrow,ncol,8)
	integer TF(Maxgauge)
c
	integer irr1, irr2, irr3, irr4
	integer irr5, irr6, irr7, irr8
	real d1, d2, d3, d4, d5, d6, d7, d8,dd0
	real wk(8),ak(8),sta(8) 
c      double precision wk,ak,sta

	integer day,im,mm
	integer N_gauge,TFF
	integer mday,iyear
c
c
	  If(TFF.eq.1) then
	  do ir = 1, nrow
	    do ic = 1, ncol
 	     if(G_elev(ir, ic).gt.0.0) then    ! add by xujijun 20060324
             
c################################################################
c		To find out the 8 nearest points
c################################################################
		do ir1=1,N_gauge
		  z0(ir1)=0
		  if(x0(ir1) .ge. (x(ir,ic)-10000000.0) .and.
     :      	 x0(ir1) .le. (x(ir,ic)+10000000.0) .and.
     :		 y0(ir1) .ge. (y(ir,ic)-10000000.0) .and.
     :       y0(ir1) .le. (y(ir,ic)+10000000.0) .and. TF(ir1).eq.1) then
	         z0(ir1)=1
		  end if
	    end do 
		d1 = 10.0E30
		d2 = 10.0E31
		d3 = 10.0E32
		d4 = 10.0E33
	    d5 = 10.0E34
	    d6 = 10.0E35
	    d7 = 10.0E36
	    d8 = 10.0E37
	    do ir1 = 1, N_gauge
		  if(z0(ir1) .eq. 1) then
            d = (x(ir,ic)-x0(ir1))**2+
     :          (y(ir,ic)-y0(ir1))**2
	      if(d .lt. d1) then
	         d8=d7
	         irr8=irr7
	         d7=d6
	         irr7=irr6
	         d6=d5
	         irr6=irr5
	         d5=d4
	         irr5=irr4
	         d4=d3
	         irr4=irr3
	         d3=d2
	         irr3=irr2
	         d2=d1
	         irr2=irr1
	         d1=d
	         irr1=ir1
	         goto 199
	       end if 
	       if(d .ge. d1 .and. d .lt. d2) then
	         d8=d7
	         irr8=irr7
	         d7=d6
	         irr7=irr6
	         d6=d5
	         irr6=irr5
	         d5=d4
	         irr5=irr4
	         d4=d3
	         irr4=irr3
	         d3=d2
	         irr3=irr2
	         d2=d
	         irr2=ir1
	          goto 199
		   end if 
	       if(d .ge. d2 .and. d .lt. d3) then
	         d8=d7
	         irr8=irr7
	         d7=d6
	         irr7=irr6
	         d6=d5
	         irr6=irr5
	         d5=d4
	         irr5=irr4
	         d4=d3
	         irr4=irr3
	         d3=d
	         irr3=ir1
	          goto 199
		   end if 
	       if(d .ge. d3 .and. d .lt. d4) then
	         d8=d7
	         irr8=irr7
	         d7=d6
	         irr7=irr6
	         d6=d5
	         irr6=irr5
	         d5=d4
	         irr5=irr4
	         d4=d
	         irr4=ir1
	          goto 199
		   end if 
	       if(d .ge. d4 .and. d .lt. d5) then
	         d8=d7
	         irr8=irr7
	         d7=d6
	         irr7=irr6
	         d6=d5
	         irr6=irr5
	         d5=d
	         irr5=ir1
	          goto 199
		   end if 
	       if(d .ge. d5 .and. d .lt. d6) then
	         d8=d7
	         irr8=irr7
	         d7=d6
	         irr7=irr6
	         d6=d
	         irr6=ir1
	          goto 199
		   end if 
	       if(d .ge. d6 .and. d .lt. d7) then
	         d8=d7
	         irr8=irr7
	         d7=d
	         irr7=ir1
	          goto 199
		   end if 
	       if(d .ge. d7 .and. d .lt. d8) then
	         d8=d
	         irr8=ir1
	          goto 199
		   end if 
		   end if
199	     continue
	    end do

	    xx(ir,ic,1)=x0(irr1)
	    xx(ir,ic,2)=x0(irr2)
	    xx(ir,ic,3)=x0(irr3)
	    xx(ir,ic,4)=x0(irr4)
	    xx(ir,ic,5)=x0(irr5)
	    xx(ir,ic,6)=x0(irr6)
	    xx(ir,ic,7)=x0(irr7)
	    xx(ir,ic,8)=x0(irr8)
	    yy(ir,ic,1)=y0(irr1)
	    yy(ir,ic,2)=y0(irr2)
	    yy(ir,ic,3)=y0(irr3)
	    yy(ir,ic,4)=y0(irr4)
	    yy(ir,ic,5)=y0(irr5)
	    yy(ir,ic,6)=y0(irr6)
	    yy(ir,ic,7)=y0(irr7)
	    yy(ir,ic,8)=y0(irr8)
	    zz(ir,ic,1) =irr1   
	    zz(ir,ic,2) =irr2   
	    zz(ir,ic,3) =irr3   
	    zz(ir,ic,4) =irr4   
	    zz(ir,ic,5) =irr5   
	    zz(ir,ic,6) =irr6   
	    zz(ir,ic,7) =irr7   
	    zz(ir,ic,8) =irr8
c		print *, irr1, irr2, irr3, irr4,irr5, irr6, irr7, irr8
c		print *, d1, d2, d3, d4,d5, d6, d7, d8
		dd(ir,ic,1)=sqrt(d1)
	    dd(ir,ic,2)=sqrt(d2)
	    dd(ir,ic,3)=sqrt(d3)
	    dd(ir,ic,4)=sqrt(d4)
	    dd(ir,ic,5)=sqrt(d5)
	    dd(ir,ic,6)=sqrt(d6)
	    dd(ir,ic,7)=sqrt(d7)
	    dd(ir,ic,8)=sqrt(d8)
	    
		mm=4
cc	    dd0=400.0*1000.0
	    dd0=300.0*1000.0
		do ir1=1, 8
		   wk(ir1)=exp(-mm*dd(ir,ic,ir1)/dd0)
c		   Wk_tmp=-exp(-mm*1.0)
c		   Wk_d=dd(ir,ic,ir1)/dd0
c		   wk(ir1)=exp(-1*mm*Wk_d*Wk_d)+Wk_tmp
c	       if(dd(ir,ic,ir1).gt.dd0) then
c	         wk(ir1)=0.0
c		   end if 
          end do
		do ir1=1,8
	       ak(ir1)=0.0
	       Sum_wl=0.0
		   do ir2=1,8
	         if(ir2 .ne. ir1) then
                 if(y0(zz(ir,ic,ir2)).ge.y(ir,ic)) then
	             sign1=1
				else
	             sign1=-1
			   end if
                 if(y0(zz(ir,ic,ir1)).ge.y(ir,ic)) then
	             sign2=1
				else
	             sign2=-1
			   end if
c			   print *,x0(ir2),x(ir,ic),dd(ir,ic,ir2)
c			   print *,x0(ir1),x(ir,ic),dd(ir,ic,ir1)
			   sta1=(x0(zz(ir,ic,ir2))-x(ir,ic))/dd(ir,ic,ir2)
			   sta2=(x0(zz(ir,ic,ir1))-x(ir,ic))/dd(ir,ic,ir1)
			   	
	           if(sta1.gt.1.0.or.sta1.lt.-1) then 
			     print *, 'sta1',sta1
                   print *, ir, ic, ir2  
	             print *, x0(zz(ir,ic,ir2)), x(ir,ic)
				 print *, zz(ir,ic,ir2),dd(ir,ic,ir2)
			    endif	   
                 if(sta2.gt.1.0.or.sta2.lt.-1)  then
			      print *, 'sta2',sta2
                    print *, ir, ic, ir1  
	              print *, x0(zz(ir,ic,ir1)), x(ir,ic)
				  print *, zz(ir,ic,ir1),dd(ir,ic,ir1)
                  endif

			   sta(ir2)=sign1*acos(sta1)
     $                   -sign2*acos(sta2)
		       ak(ir1)=ak(ir1)+wk(ir2)*(1-cos(sta(ir2)))
	           Sum_wl=Sum_wl+wk(ir2)
			 end if  
		   end do
	       if(Sum_wl.gt.0) then
		     ak(ir1)=ak(ir1)/Sum_wl
	        else
		     ak(ir1)=0.0
	       end if
          end do
		
	    Sum_wk=0.0
		do ir1=1,8
		   wk(ir1)=wk(ir1)*(1+ak(ir1))
             Sum_wk=Sum_wk+wk(ir1)
             ww(ir,ic,ir1)=wk(ir1)
          end do
		ww(ir,ic,9)=Sum_wk
	     
		 endif
	    end do
	  end do
	  end if 
c
	  do day=1,mday
	  do ir = 1, nrow
	    do ic = 1, ncol
	     if(G_elev(ir, ic).gt.0.0) then    ! add by xujijun 20060324

c################################################################
c	To interplate the value using the inverse-distance weight
c################################################################

		  dataout(ir,ic,day)=0.0 
		  Sum_wk=ww(ir,ic,9)
		  do ir1=1, 8
			if(datain(zz(ir,ic,ir1),day).ge.0.0) then
                dataout(ir,ic,day)=dataout(ir,ic,day)+ww(ir,ic,ir1)
     $                     *datain(zz(ir,ic,ir1),day)
	         else
			  Sum_wk=Sum_wk-ww(ir,ic,ir1)  
	        end if
	      end do
            if(Sum_wk.ge.ww(ir,ic,9)/8.0) then
		    dataout(ir,ic,day)=dataout(ir,ic,day)/Sum_wk
             else
		    dataout(ir,ic,day)=-9999
	        d_temp=0.0
	        n_temp=0
			do ir1=8,1,-1
			 if(datain(zz(ir,ic,ir1),day).ge.0.0) then 
		       n_temp=n_temp+1
			   d_temp=d_temp+datain(zz(ir,ic,ir1),day)
               end if
	        end do
	        if(n_temp.gt.0) then
	          dataout(ir,ic,day)=d_temp/n_temp
			end if
		  end if
	    
		 endif
		end do
	  end do

	  end do
	
	end
c
c*************************************************************************
c
      subroutine T_pcir_coeff(start_year,N_gauge,datain,Hbeta
     $                        ,elev,y0,mmonth,G_zone)
       include 'dims2.inc'
	
	real datain(Maxgauge,12,31),datam(Maxgauge,12)
	real Hbeta(5,12,4)
	real y0(Maxgauge),elev(Maxgauge)
	integer mmonth(12),G_zone(Maxgauge)
      integer start_year,im,iyear,N_gauge,izone
	dimension xp(500000),yp(500000),a(6),zp(500000)
      double precision xp,yp,a,dt1,dt2,dt3,zp
c
c      write(6,*) "Zone    A1      A2       A3      A4"
	do izone=1,5
	  
	  do ir=1,Maxgauge
	      do it=1,12
	        datam(ir,it)=0.0
	      end do
	  end do
	
	do im=1,12
	  N_pcir=0.0
	    if(mod(start_year,4).eq.0) then
	      mmonth(2)=29
	     else
	      mmonth(2)=28
	    end if

	    do iday=1,mmonth(im)    
		  do ir = 1, N_gauge
	        if(datain(ir,im,iday).gt.-60) then
			  datam(ir,im)=datam(ir,im)
     $                +datain(ir,im,iday)/mmonth(im)
	         end if
		  end do
		end do

		do ir = 1, N_gauge
	      do ic = 1, N_gauge
	      if(ic.ne.ir) then
			if(datam(ir,im).gt.-60.0
     $     .and. datam(ic,im).gt.-60.0
c     $     .and. elev(ir).lt.1500.0.and. elev(ic).lt.1500.0) then
c     $     .and. elev(ir).gt.2000.0.and. elev(ic).gt.2000.0) then
     $     .and. G_zone(ir).eq.izone.and. G_zone(ic).eq.izone) then
				N_pcir=N_pcir+1
			    xp(N_pcir)=elev(ir)-elev(ic)
			    yp(N_pcir)=datam(ir,im)
     $                      -datam(ic,im)
	        end if
		  end if
		  end do
		end do

	  n=N_pcir
	  m=2
        call pcir(xp,yp,a,n,m,dt1,dt2,dt3)
	  Hbeta(izone,im,1)=a(1)
	  Hbeta(izone,im,2)=a(2)
c	  print *,"elev  ",a(1),a(2)
	end do
c
      do im=1,12
	  N_pcir=0.0
	    if(mod(start_year,4).eq.0) then
	      mmonth(2)=29
	     else
	      mmonth(2)=28
	    end if
		do ir = 1, N_gauge
	      do ic = 1, N_gauge
	      if(ic.ne.ir) then
			if(datam(ir,im).gt.-60.0
     $     .and. datam(ic,im).gt.-60.0
c     $     .and. elev(ir).lt.1500.0.and. elev(ic).lt.1500.0) then
c     $     .and. elev(ir).gt.2000.0.and. elev(ic).gt.2000.0) then
     $     .and. G_zone(ir).eq.izone.and. G_zone(ic).eq.izone) then
				N_pcir=N_pcir+1
			    zp(N_pcir)=(y0(ir)-y0(ic))/10000.0
                  dz=elev(ir)-elev(ic)
			    yp(N_pcir)=datam(ir,im)
     $                      -datam(ic,im)
     $                      +Hbeta(izone,im,1)+Hbeta(izone,im,2)*dz
 
	        end if
		  end if
		  end do
		end do
	  n=N_pcir
	  m=2
        call pcir(zp,yp,a,n,m,dt1,dt2,dt3)
	  Hbeta(izone,im,3)=a(1)
	  Hbeta(izone,im,4)=a(2)
c	  print *,a(1),a(2)
	end do
	end do

	end
c
c
c*************************************************************************
c
      subroutine alter_item(im,prec,Tmean,Tmax,Tmin,Windm,RH1,Tsun)
      real prec,Tmean,Tmax,Tmin,Windm,RH1,Tsun
      integer im

             If (Tmean .eq. 32744 .or. Tmean .eq. 32766) Then   !changed by jiaoy
c            print *,'im,Tmean=',GG_id,im, Tmean
c            pause
             Tmean = 100
c              Tmean = 10
           endif
             Tmean=0.1*Tmean
c             Tmean=Tmean

             If (Tmax .eq. 32744 .or. Tmax .eq. 32766) Then   !changed by jiaoy
c            print *,'Tmax=', im, Tmax
             Tmax = 200
c              Tmax = 20
           endif
             Tmax=0.1*Tmax
c             Tmax=Tmax

             If (Tmin .eq. 32744 .or. Tmin .eq. 32766) Then   !changed by jiaoy
c            print *,'Tmin=', im, Tmin
             Tmin = 0
            endif
             Tmin=0.1*Tmin
c             Tmin=Tmin

           if(Tmean .lt. -10 .and.
     $       (Tmax .gt. -6 .and. Tmin .gt. -6)) then
                Tmean=(Tmax+Tmin)/2.0
              end if
           if(Tmax .lt. -10 .and.
     $       (Tmean .gt. -6 .and. Tmin .gt. -6)) then
                Tmax=2*Tmean-Tmin
              end if
           if(Tmin .lt. -10 .and.
     $       (Tmean .gt. -6 .and. Tmax .gt. -6)) then
                Tmin=2*Tmean-Tmax
              end if

             If (Rh1 .eq. 32744 .or. Rh1 .eq. 32766 ) Then   !changed by jiaoy
c            print *,'Rh1=', im, Rh1
             Rh1 = 50
c              Rh1 = 0.50
           endif
              RH1= 0.01*RH1
c              RH1=RH1

             If (Windm .eq. 32744 .or. Windm .eq. 32766 ) Then   !changed by jiaoy
c            print *,'Windm=', im, Windm
             Windm = 30
c              Windm = 3.0
           endif
             If (Windm .ge. 1000) Then          !changed by jiaoy
c              print *,'Windm=', im, Windm
               Windm = Windm - 1000
             endif
            Windm=0.1*Windm
c            Windm=Windm

             If (Tsun .ge. 200 ) Then   !changed by xujijun
c            print *,'Tsun=', im, Tsun
               Tsun = 50
c                Tsun = 5
           End If
             Tsun=0.1*Tsun
c             Tsun=Tsun

             If (prec .ge. 32700 ) Then    !changed by jiaoy
c            print *,'prec=', im, prec
               prec = 0
           End If
             If (prec .ge. 32000 ) Then    !changed by xujijun
c             print *,'prec=', im, prec
               prec = prec - 32000
             End If
             If (prec .ge. 31000 ) Then    !changed by xujijun
c             print *,'prec=', im, prec
               prec = prec - 31000
             End if
             If (prec .ge. 30000 ) Then    !changed by xujijun
c             print *,'prec=', im, prec
               prec = prec - 30000
           End If
           prec=0.1*prec
c           prec=prec
      end
c
*	subroutine alter_item(im,prec,Tmean,Tmax,Tmin,Windm,RH1,Tsun)
*	real prec,Tmean,Tmax,Tmin,Windm,RH1,Tsun
*	integer im 

*  		 If (Tmean .le. -8000 ) Then   !changed by xujijun
c	       print *,'im,Tmean=', im, Tmean
c             Tmean = 100
*              Tmean = 10
*           endif 
c	       Tmean=0.1*Tmean
*             Tmean = Tmean
*  		 If (Tmax .le. -8000 ) Then   !changed by xujijun
c	       print *,'Tmax=', im, Tmax
c             Tmax = 200
*              Tmax=20
*           endif 
c  	       Tmax=0.1*Tmax
*             Tmax = Tmax
*  		 If (Tmin .le. -8000 ) Then   !changed by xujijun
c	       print *,'Tmin=', im, Tmin
*             Tmin = 0
*            endif
c	       Tmin=0.1*Tmin
*             Tmin = Tmin
*	     if(Tmean .lt. -100 .and. 
*     $       (Tmax .gt. -60 .and. Tmin .gt. -60)) then
*                Tmean=(Tmax+Tmin)/2.0
*		  end if
*	     if(Tmax .lt. -100 .and. 
*     $       (Tmean .gt. -60 .and. Tmin .gt. -60)) then
*                Tmax=2*Tmean-Tmin
*		  end if
*	     if(Tmin .lt. -100 .and. 
*     $       (Tmean .gt. -60 .and. Tmax .gt. -60)) then
*                Tmin=2*Tmean-Tmax
*		  end if

*  		 If (RH1 .le. -8000 ) Then   !changed by xujijun
c	       print *,'Rh1=', im, Rh1
*             RH1 = 0.50
*           endif 
c		  RH1= 0.01*RH1 
*            RH1= RH1

*  		 If (Windm .le. -8000 ) Then   !changed by xujijun
c	       print *,'Windm=', im, Windm
c             Windm = 30
*              Windm = 3.0
*           endif 
c	      Windm=0.1*Windm
*            Windm = Windm 
*  		 If (Tsun .le. -200 ) Then   !changed by xujijun
c	       print *,'Tsun=', im, Tsun
c               Tsun = 50
*                Tsun = 5
*           End If
c	       Tsun=0.1*Tsun
*             Tsun=Tsun

*		 If (prec .le. -2000 ) Then    !changed by xujijun
*             prec = 0
*           End If
c	     prec=0.1*prec
*           prec=prec
*      end
c
c*********************************************************************************
c     least square regress  y=a1+a2*(x-Xa)+a3*(x-Xa)^2+a4*(x-Xa)^3+a5*(x-Xa)^4+a6*(x-Xa)^5, where Xa is average of all x
c*********************************************************************************
c
      subroutine pcir(x,y,a,n,m,dt1,dt2,dt3)
	dimension x(n),y(n),a(m),s(500000),t(500000),b(500000)
	double precision x,y,a,s,t,b,dt1,dt2,dt3,z,d1,p,c,d2,g,q,dt
	integer i,j,k
      z=0.0
	do 10 i=1,n
	  z=z+x(i)/n
10    continue
      
	b(1)=1.0
	d1=n
	p=0.0
	c=0.0
	do 20 i=1,n
	  p=p+(x(i)-z)
        c=c+y(i)
20    continue
	
	c=c/d1
	p=p/d1
	a(1)=c*b(1)
	if(m.gt.1) then
	  t(2)=1.0
	  t(1)=-p
	  d2=0.0
	  c=0.0
	  g=0.0
	  do 30 i=1,n
	    q=x(i)-z-p
	    d2=d2+q*q
	    c=y(i)*q+c
	    g=(x(i)-z)*q*q+g
30	  continue
        c=c/d2
	  p=g/d2
	  q=d2/d1
	  d1=d2
	  a(2)=c*t(2)
	  a(1)=c*t(1)+a(1) 
	end if

	do 100 j=3,m
	  s(j)=t(j-1)
	  s(j-1)=-p*t(j-1)+t(j-2)
	  if(j .GE. 4) then
	    do 40 k=j-2,2,-1
	      s(k)=-p*t(k)+t(k-1)-q*b(k)
40        continue
	  end if
	  s(1)=-p*t(1)-q*b(1)
	  d2=0.0
	  c=0.0
	  g=0.0
	  do 70 i=1,n
	    q=s(j)
	    do 60 k=j-1,1,-1
	      q=q*(x(i)-z)+s(k)
60        continue
          d2=d2+q*q
		c=y(i)*q+c
		g=(x(i)-z)*q*q+g  
70      continue
        c=c/d2
	  p=g/d2
	  q=d2/d1
	  d1=d2
	  a(j)=c*s(j)
	  t(j)=s(j)
	  do 80 k=j-1,1,-1
	    a(k)=c*s(k)+a(k)
	    b(k)=t(k)
		t(k)=s(k) 
80      continue	   	   
100   continue

      dt1=0.0
	dt2=0.0
	dt3=0.0
	do 120 i=1,n
	  q=a(m)
	  do 110 k=m-1,1,-1
	    q=q*(x(i)-z)+a(k)
110     continue
        dt=q-y(i)
	  if(abs(dt).gt.dt3) dt3=abs(dt)
	  dt1=dt1+dt*dt
	  dt2=dt2+abs(dt) 
120   continue
	return 
	end


c*******************************************************************************
c	To estimate potential evaporation
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*
c 
      subroutine cal_evp(Jd, land_type, t_a, t_max,t_min,
     :                cloud,sunshine,RH,u,z, lat, Ep,prec0, R_n,dv)
      implicit none

      integer land_type !land type
      real t_a      !Ta  mean air temperature(degree)
      real t_max    !	   maximum temperature(degree)
      real t_min    !	   minimum temperature(degree)
	real t_dew    !    dew temperature(degree)
      real cloud	  !    cloud cover [0.**]
      real prec0	  !    rainfall    [mm]
      real sunshine !    sunshine time [hr]
      real RH       !    relative humidity,units 0.###
      real u		  !u(z)wind speed [m/s]
      real z		  !     height of location
	real u_2      !wind speed at 2-m
      real phi      !latitude of site (positive for Nothern hemisphere) [radian]
c	
      real pi
      parameter (pi=3.1415926)
c	
      real lamta !    latent heat of vaporization of water
      real delta !slope of the curve of satuated vapor pressure at Ta [kPa/degC]
      real e_s	  !es  satuated vapor pressure at Ta [kPa]
c	
      real R_n	  !Rn  net radiation [mm/day]
	real S_n      !
	real L_n      !
      real cloud_f  !n/N cloudiness_function [h/h]
      real sigma	  !    Stefan-Boltzmann constant
      real e_d	  !ed  vapor pressure at dew tempreture	[]
      real s_0	  !S0  extraterrestrial radiation [mm/day]
      real d_r      !dr  relative distance between the earth and the sun
      real delta_s  !  solar declination [radian]
      real omega_s  !s sunset hour angle [radian]
c
      real G		  !G   heat flux into soil or snow layer []
      real gamma	  !    psychrometric constant [kPa/degC]
      real gammas	  !    psychrometric constant [kPa/degC]
      real p0       !    atmospheric pressure [kPa]
      real p_s	  !ps   [kPa]
      real r_d	  !Rd   [J/kg/K]
      real k		  !k
      real lgamma	  ! [ /km]
c	
      real D		  !D   vapor pressure deficit [kPa]
	real Ah       !
	real as       !
	real bs       !
	real ac       !
	real bc       !
	real ae       !
	real be       !
	real S_t      !
	real St0      !
	real f        !
	real epocilon !
	real e_smax   !
	real e_smin   !
	real e_sdew   !
c	
      real alpha    ! albedo for each land-uses
      real Ep       ! daily potential E (mm/day)

	integer Jd    !julian day number =day
	real NN       !total day length(hr)
	real lat
      real dv


!      Jd=151
!       land_type=4
!       t_a=13.14776
!        t_max=19.684
!        t_min=6.9236
!       cloud=-9999
!       sunshine=5.168966
!       RH=0.637
!       u=3.47
!       z=453
!        lat=53.52
!         Ep=-999
!         prec0=2.24

c
c      if(land_type .eq. 1)  alpha = 0.08  ! water body
c      if(land_type .eq. 2)  alpha = 0.23  ! urban area <-not certain
c      if(land_type .eq. 3)  alpha = 0.20  ! bare soil
c      if(land_type .eq. 4)  alpha = 0.13  ! forest
c      if(land_type .eq. 5)  alpha = 0.23  ! sparse vegetation
c      if(land_type .eq. 6)  alpha = 0.23  ! agriculture land 
c      if(land_type .eq. 7)  alpha = 0.23  ! grassland
c      if(land_type .eq. 8)  alpha = 0.23  ! shurub
c      if(land_type .eq. 9)  alpha = 0.23  ! wetland area
c      if(land_type .eq. 10) alpha = 0.6   ! snow
c
        if(land_type .eq. 1)  alpha = 0.05  ! water body
      if(land_type .eq. 2)  alpha = 0.20  ! urban area <-not certain
      if(land_type .eq. 3)  alpha = 0.25  ! bare soil
      if(land_type .eq. 4)  alpha = 0.13  ! forest
CC      if(land_type .eq. 5)  alpha = 0.14  ! forest
      if(land_type .eq. 5)  alpha = 0.13  ! irrigated-cropland
cc      if(land_type .eq. 7)  alpha = 0.23  ! irrigated-cropland
cc      if(land_type .eq. 8)  alpha = 0.23  ! irrigated-cropland
cc      if(land_type .eq. 9)  alpha = 0.23  ! irrigated-cropland
cc      if(land_type .eq. 10) alpha = 0.23  ! irrigated-cropland
cc      if(land_type .eq. 11) alpha = 0.23  ! irrigated-cropland
cc      if(land_type .eq. 12) alpha = 0.23  ! irrigated-cropland
      if(land_type .eq. 6) alpha = 0.13  ! upland
cc      if(land_type .eq. 14) alpha = 0.22  ! grassland
      if(land_type .eq. 7) alpha = 0.16  ! grassland
cc      if(land_type .eq. 16) alpha = 0.24  ! grassland
      if(land_type .eq. 8) alpha = 0.17  ! shrub
      if(land_type .eq. 9) alpha = 0.09  ! wetland
cc      if(land_type .eq. 19) alpha = 0.6   ! frost snow
cc      if(land_type .eq. 20) alpha = 0.6   ! others
c      if(prec0.gt.10) then
c      print *,Jd, land_type, t_a, t_max,t_min,
c     :                cloud,sunshine,RH,u,z, lat, Ep,prec0, R_n,dv
c      pause
c      endif  
        

      alpha = alpha + dv
	sigma = 4.903E-9
     	G = 0
	Ah= 4.19*10**(-3)*prec0*t_a
      lamta = 2.501 - 0.00236*t_a
	
      if(t_a.lt.0.0) then
	  lamta = lamta + 0.334
	end if
	 
	phi = lat *pi/180.0

      e_smax   = 0.6108 * exp(17.27*t_max/(237.3+t_max))
      if(t_max.lt.0.0) then
        e_smax   = 0.6108 * exp(21.88*t_max/(265.5+t_max))
	end if
      e_smin   = 0.6108 * exp(17.27*t_min/(237.3+t_min))
      if(t_min.lt.0.0) then
        e_smin   = 0.6108 * exp(21.88*t_min/(265.5+t_min))
	end if
	e_s      = (e_smax+e_smin)/2.0 
	
	if(RH.gt.0) then
	 if(e_s.gt.0.6108)then		! esat(t=0 degree c)=6.11
	  t_dew=237.3*log10(RH*e_s/0.6108)/(7.5-log10(RH*e_s/0.6108))
cc     ! Reverse form of Tetens    	! formula for water surface
	 else
	  t_dew=265.3*log10(RH*e_s/0.6108)/(9.5-log10(RH*e_s/0.6108)) 
	 ! Reverse form of Tetens    	! formula for ice surface
	 end if
	else
	 t_dew = t_min-2.0  !only for arid zone
	end if
	
	if(R_n .eq. -9999) then
	  delta_s = 0.4093*sin(2.0*pi*Jd/365.0-1.405)
	  omega_s = acos(-tan(phi)*tan(delta_s))
	  NN    = 24.0*omega_s/pi
	  d_r   = 1.0+0.033*cos(2.0*pi*Jd/365.0)
        s_0   = 15.392*d_r*(omega_s*sin(phi)*sin(delta_s) 
     :          + cos(phi)*cos(delta_s)*sin(omega_s))       ! mm/day

c	  if(sunshine.gt.NN) print *,'wrong in calculating NN',NN,sunshine !由mqh标注
	  if(sunshine.gt. NN) sunshine = NN
c
        as=0.25
	  bs=0.5
	  ac=1.35
	  bc=-0.35
	  ae=0.34
	  be=-0.14

	  if(sunshine .eq. -9999) then
	    if(cloud.lt.0) then
		  print *, "wrong in data cloud and sunshine"
		  cloud=0.5
		  end if 
		  cloud_f=1-cloud
	  else
	    cloud_f=sunshine/NN
	  end if 
c        print *, 'cloud_f=',cloud_f 

	  S_t = (as+bs*cloud_f)*s_0
	  S_n = S_t*(1.0-alpha)
c	estimating long-wave radiation Ln
	  St0      = (as+bs)*s_0
	  f        = ac*(S_t/St0)+bc
c     	  e_d      = 0.6108 * exp(17.27*t_min/(237.3+t_min))
     	  e_d      = 0.6108 * exp(17.27*t_dew/(237.3+t_dew))
	  epocilon = ae+be*sqrt(e_d)
	  L_n       = -f*epocilon*sigma*(t_a+273.16)**4
	  L_n       = L_n/lamta              !in mm/day
c	total net radiation
	  R_n = S_n + L_n
	  if(R_n .lt. 0.0) R_n = 0.0
	 else
c       input R_n	  
       end if

       e_s   = 0.6108 * exp(17.27*t_a/(237.3+t_a))
       delta = 4098 * e_s / ((237.3+t_a)**2) 
       if(t_a.lt.0.0) then
        e_s   = 0.6108 * exp(21.88*t_a/(265.5+t_a))
        delta = 5809 * e_s / ((265.5+t_a)**2)
	 end if

 	p0  = 101.325 ! kPa
c      r_d = 287.0
c      lgamma = 6.5  ![ ?? /km]	for z=0-11km
c      k      = 9.8066 / (r_d*lgamma)
c      p_s = p0 * ((t_a+273.15-lgamma*z)/(t_a+273.15))**k
      p_s = p0 * ((293-0.0065*z)/293)**5.256
      gamma = 0.0016286*p_s/lamta

c      e_smax   = 0.6108 * exp(17.27*t_max/(237.3+t_max))
c      e_smin   = 0.6108 * exp(17.27*t_min/(237.3+t_min))
	e_sdew   = 0.6108 * exp(17.27*t_dew/(237.3+t_dew))
      if((e_smax+e_smin)/2.0 .lt. 0.6108) then
        e_sdew   = 0.6108 * exp(21.88*t_dew/(265.5+t_dew))
	end if
      
	if(RH .le. 0.0) then
	  D=(e_smax+e_smin)/2.0-e_sdew
       else
	  D=(e_smax+e_smin)/2.0*(1.0-RH)   
      end if
	
	if(u.gt.0) then
	  u_2=0.749*u
	 else
	  u_2=0.0
      end if

	gammas=gamma*(1+0.33*u_2)

      if(land_type .eq. 1) then
	  Ep=delta/(delta+gamma)*(R_n+Ah)
     $     +gamma/(delta+gamma)*6.43*(1+0.536*u_2)*D/lamta
	 else
	  Ep=delta/(delta+gammas)*(R_n-G)
     $     +gamma/(delta+gammas)*900.0*u_2*D/(t_a+275)
      end if

	if(Ep.lt.0.0) then
	  print *,'wrong in Ep',Ep
      end if
c	print *, im, Jd, evap_daily
!            print *,'Jd, land_type, t_a, t_max,t_min,
!     :                cloud,sunshine,RH,u,z, lat, Ep,prec0, R_n,dv',
!     :Jd, land_type, t_a, t_max,t_min,
!     :                cloud,sunshine,RH,u,z, lat, Ep,prec0, R_n,dv
!      pause

	end
           

      subroutine rain_downscale2(N_gauge,N_gauge1,
     $   monthd, iyear, im,Ppr,Pprhour)

      include 'dims2.inc'
      integer N_gauge, N_gauge1, im, N_g, rt1(Maxgauge,31)
	integer monthd(12),igauge,iday,idh,idum, iyear,idayx
	real Ppr(Maxgauge,31), Pprhour(Maxgauge,31,24)
      real tmp
      integer iy2
      integer startofmonth(maxgauge,20),startofday(maxgauge,20),
     $ endofmonth(maxgauge,20),endofday(maxgauge,20)  !mqh
      common/start_end/startofmonth,startofday,endofmonth,endofday
      
      iy2=iyear-1990
      
c            N_g = N_gauge1+N_gauge
c	       do iday = 1, monthd(im)
c              do igauge=1, N_g  
           
c        if (Ppr(igauge,iday).ne.0.and.Ppr(igauge,iday).ne.-9999) then  
c             print *,Ppr(igauge,iday),igauge,iday
c             endif
      
c           enddo
c            enddo
c             pause

      
c        print *, 'N_gauge,N_gauge1=', N_gauge,N_gauge1
c        print *, 'rain_downscale,im,monthd(im)=', im,monthd(im)
        if(NNstation.eq.1) N_g = N_gauge
        if(NNstation.eq.2) N_g = N_gauge1
        if(NNstation.eq.3) N_g = N_gauge1+N_gauge

c           do igauge=N_gauge+1, N_g   
c           print *,startofmonth(igauge-N_gauge,iy2),
c     :startofday(igauge-N_gauge,iy2),endofmonth(igauge-N_gauge,iy2),
c     :endofday(igauge-N_gauge,iy2),igauge
c            enddo
c            pause
    
    
c      print *, 'N-g=', N_g
        do igauge=1, N_g 
      
c     if(igauge.gt.N_gauge)   print *,startofmonth(igauge-N_gauge,iy2),
c     :startofday(igauge-N_gauge,iy2),endofmonth(igauge-N_gauge,iy2),
c     :endofday(igauge-N_gauge,iy2),igauge           
         do  iday = 1, monthd(im)          
         
          if(igauge.gt.N_gauge) then   !added by mqh   
           if(im.lt.endofmonth(igauge-N_gauge,iy2)) then
                if (im.gt.startofmonth(igauge-N_gauge,iy2)) then
c                print *,'01',im,iday
                  goto 1992
                endif
           endif
           if (im.eq.startofmonth(igauge-N_gauge,iy2)) then
                  if(iday.ge.startofday(igauge-N_gauge,iy2)) then
c                  print *,'011',im,iday
c                     pause
                  goto 1992
                   endif
           endif
          if (im.eq.endofmonth(igauge-N_gauge,iy2)) then
                  if(iday.le.endofday(igauge-N_gauge,iy2)) then
c                  print *,'0111',im,iday
c                     pause
                    goto 1992        
                 endif
           endif
          endif
          !end of added,mqh
            do idh=1, 24 
                Pprhour(igauge,iday,idh)=0.0
	      enddo
c	   IF(igauge.gt.10)   print *,igauge,im,iday,'hkhhhkhjkhj'
1992	   enddo    !1992,added by mqh
	  enddo

	  
c	  print *,idayx

        do igauge=1, N_g   
         do iday = 1, monthd(im)
	      if(Ppr(igauge,iday).lt.0.0) Ppr(igauge,iday)=0.0 
	      if(Ppr(igauge,iday).gt.10000) Ppr(igauge,iday)=0.0 
	   enddo
	  enddo

        do igauge=1, N_g   
         do  iday = 1, monthd(im)
            if(igauge.gt.N_gauge) then   !added by mqh   
            if(im.lt.endofmonth(igauge-N_gauge,iy2)) then
                if (im.gt.startofmonth(igauge-N_gauge,iy2)) then
c                print*,'1',iday,im
                   goto 0404 
                endif
            endif
            if (im.eq.startofmonth(igauge-N_gauge,iy2)) then
                  if(iday.ge.startofday(igauge-N_gauge,iy2)) then
c                   print*,'2',iday,im
                  goto 0404
                   endif
             endif
           if (im.eq.endofmonth(igauge-N_gauge,iy2)) then
                  if(iday.le.endofday(igauge-N_gauge,iy2)) then
c                   print*,'3',iday,im
                    goto 0404        
                  endif
           endif
          endif

          !end of added,mqh
         
         
         
c------------------------ create the start time ------------
            idayx=iday
            tmp = ran0(idum)
            rt1(igauge,iday) = 4 + int((21-4)*tmp+0.45)
            if(rt1(igauge,iday).gt. 21) rt1(igauge,iday) = 21
c           print *,"rt1(igauge,idayx)=", rt1(igauge,idayx), tmp,idum

c------------------------ downscaling Prec --------------------
c      print *,'igauge, iday', igauge, iday, Ppr(igauge,iday)

       if(NNstation.eq.2.or.NNstation.eq.3.and.igauge.gt.N_gauge) then
c        rt1(igauge,iday) = rt1(igauge,iday)+8

!雨量站是8-：00-8：00，气象站是20-20
c  在资料中，雨量站的降雨峰值要比气象站的日期提前，相差12小时而不是8小时。

         rt1(igauge,iday) = rt1(igauge,iday) !已经修改了
       if(rt1(igauge,iday).gt.21) then
            idayx=iday+1
              if(idayx.gt.monthd(im)) idayx = iday
            if(iday.ge.2) then
            if(rt1(igauge,iday-1).le.21) then
                 Ppr(igauge,idayx) = Ppr(igauge,idayx)+Ppr(igauge,iday)
                 Ppr(igauge,iday)  = 0
              else
               Ppr(igauge,idayx) = Ppr(igauge,idayx)+
     &                            Ppr(igauge,iday)-Ppr(igauge,iday-1)
                 Ppr(igauge,iday)  = Ppr(igauge,iday-1)
               rt1(igauge,iday)  = 6
              endif
              endif
            endif
       endif

       do idh=1, 24
       if(im.le. 5) then
c       Pprhour(igauge,iday,idh) = 0.0
       if(Ppr(igauge,iday).ge.0.0.and.Ppr(igauge,iday).le. 8.0) then
          if(idh.eq.rt1(igauge,iday))
     &       Pprhour(igauge,iday,idh) =  Ppr(igauge,iday)
        elseif(Ppr(igauge,iday).gt.8.0.and.
     &                                Ppr(igauge,iday).le.15.0) then
          if(idh.eq.rt1(igauge,iday)-1)
     &         Pprhour(igauge,iday,idh) = 0.1*Ppr(igauge,iday)
          if(idh.eq.rt1(igauge,iday))
     &         Pprhour(igauge,iday,idh) = 0.8*Ppr(igauge,iday)
          if(idh.eq.rt1(igauge,iday)+1)
     &         Pprhour(igauge,iday,idh) = 0.1*Ppr(igauge,iday)
        elseif(Ppr(igauge,iday).gt.15.0.and.
     &                               Ppr(igauge,iday).le.25.0) then
          if(idh.eq.rt1(igauge,iday)-2)
     &             Pprhour(igauge,iday,idh) = 0.1*Ppr(igauge,iday)
          if(idh.eq.rt1(igauge,iday)-1)
     &         Pprhour(igauge,iday,idh) = 0.2*Ppr(igauge,iday)
          if(idh.eq.rt1(igauge,iday))
     &         Pprhour(igauge,iday,idh) = 0.4*Ppr(igauge,iday)
          if(idh.eq.rt1(igauge,iday)+1)
     &         Pprhour(igauge,iday,idh) = 0.2*Ppr(igauge,iday)
          if(idh.eq.rt1(igauge,iday)+2)
     &         Pprhour(igauge,iday,idh) = 0.1*Ppr(igauge,iday)
        elseif(Ppr(igauge,iday).gt.25.0.and.
     &                              Ppr(igauge,iday).le.50.0) then
          if(idh.eq.rt1(igauge,iday)-4)
     &         Pprhour(igauge,iday,idh) = 0.05*Ppr(igauge,iday)
          if(idh.eq.rt1(igauge,iday)-3)
     &         Pprhour(igauge,iday,idh) = 0.1*Ppr(igauge,iday)
          if(idh.eq.rt1(igauge,iday)-2)
     &         Pprhour(igauge,iday,idh) = 0.1*Ppr(igauge,iday)
          if(idh.eq.rt1(igauge,iday)-1)
     &         Pprhour(igauge,iday,idh) = 0.15*Ppr(igauge,iday)
          if(idh.eq.rt1(igauge,iday))
     &         Pprhour(igauge,iday,idh) = 0.2*Ppr(igauge,iday)
          if(idh.eq.rt1(igauge,iday)+1)
     &         Pprhour(igauge,iday,idh) = 0.15*Ppr(igauge,iday)
          if(idh.eq.rt1(igauge,iday)+2)
     &         Pprhour(igauge,iday,idh) = 0.1*Ppr(igauge,iday)
          if(idh.eq.rt1(igauge,iday)+3)
     &         Pprhour(igauge,iday,idh) = 0.1*Ppr(igauge,iday)
          if(idh.eq.rt1(igauge,iday)+4)
     &         Pprhour(igauge,iday,idh) = 0.05*Ppr(igauge,iday)
        elseif(Ppr(igauge,iday).gt.50.0) then
                   Pprhour(igauge,iday,idh) = Ppr(igauge,iday)/24.0
        end if
      else
c        Pprhour(igauge,iday,idh) = 0.0
        if(Ppr(igauge,iday).ge.0.0.and.Ppr(igauge,iday).le. 5.0) then
          if(idh.eq.rt1(igauge,iday))
     &           Pprhour(igauge,iday,idh) = Ppr(igauge,iday)
         elseif(Ppr(igauge,iday).gt.5.0.and.
     &                             Ppr(igauge,iday).le.10.0) then
          if(idh.eq.(rt1(igauge,iday)-1))
     &           Pprhour(igauge,iday,idh) = 0.30*Ppr(igauge,iday)
          if(idh.eq.rt1(igauge,iday))
     &           Pprhour(igauge,iday,idh) = 0.40*Ppr(igauge,iday)
          if(idh.eq.(rt1(igauge,iday)+1))
     &           Pprhour(igauge,iday,idh) = 0.30*Ppr(igauge,iday)
         elseif(Ppr(igauge,iday).gt.10.0.and.
     &                             Ppr(igauge,iday).le.30.0) then
          if(idh.eq.(rt1(igauge,iday)-3))
     &           Pprhour(igauge,iday,idh) = 0.10*Ppr(igauge,iday)
          if(idh.eq.(rt1(igauge,iday)-2))
     &           Pprhour(igauge,iday,idh) = 0.15*Ppr(igauge,iday)
          if(idh.eq.(rt1(igauge,iday)-1))
     &           Pprhour(igauge,iday,idh) = 0.15*Ppr(igauge,iday)
          if(idh.eq.rt1(igauge,iday))
     &           Pprhour(igauge,iday,idh) = 0.20*Ppr(igauge,iday)
          if(idh.eq.(rt1(igauge,iday)+1))
     &           Pprhour(igauge,iday,idh) = 0.15*Ppr(igauge,iday)
          if(idh.eq.(rt1(igauge,iday)+2))
     &           Pprhour(igauge,iday,idh) = 0.15*Ppr(igauge,iday)
          if(idh.eq.(rt1(igauge,iday)+3))
     &           Pprhour(igauge,iday,idh) = 0.10*Ppr(igauge,iday)
        elseif(Ppr(igauge,iday).gt.30.0) then
            Pprhour(igauge,iday,idh) = Ppr(igauge,iday)/24.0
c          if(idh.le.3  .or. idh.ge.22) r = 0.5 * Ppr(igauge,iday)/24.0
c          if(idh.ge.11 .and. idh.le.13) r = 2.0 * Ppr(igauge,iday)/24.0
	  end if
      endif
	 enddo
0404   enddo
	 enddo 
	 
	       
c	        do igauge=1, N_g 
c	       do iday = 1, monthd(im)              
c                do idh=1, 24 
c              if (Pprhour(igauge,iday,idh).gt.1.2.and.igauge.gt.10) then  
c             print *,Pprhour(igauge,iday,idh),igauge,iday,idh
c             endif
c             enddo
c             enddo
c             enddo
c        pause      
	 
      return
	end           






c    
cc     FUNCTION ran0(idum)
cc      INTEGER idum,IA,IM,IQ,IR,MASK
cc      REAL ran0,AM
cc      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
cc     *MASK=123459876)
cc      INTEGER k
cc      idum=ieor(idum,MASK)
cc      k=idum/IQ
cc      idum=IA*(idum-k*IQ)-IR*k
cc      if (idum.lt.0) idum=idum+IM
cc      ran0=AM*idum
cc      idum=ieor(idum,MASK)
cc      return
cc      END