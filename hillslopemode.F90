!            ***** Hydrological Response on Hillslope  *****          *
      subroutine hillslope_model(isub)
        implicit none
!       include 'dims2.inc'
        


        integer       isub              ! sub-basin number
        integer       i_sub             ! sub-basin number 
        integer iflow            ! flow-interval number
        real smf                        ! snow melting factor
        integer nsoil                    ! number of soil types in whole area
        integer isoil                    ! soil number
        integer nland                    ! number of landuse types in whole area
        integer iland                    ! land-use number

        real    kmoist                ! evaporation coefficient of soil mositure  
        real    soil_con_f       ! soil water conservation factor
        real    detension        ! surface detension (mm)
        real    detension2  !MQH


        real    dt             ! time step (second)
        integer year           ! calendar year   (19XX)
        integer month          ! month in a year (1-12)
        integer day            ! day in a month  (1-30)
        integer hour           ! hour in a day   (1-24)
        integer hydroyear      ! year for hydro simulation
        integer startmonth     ! start month in the hydro-year
        integer endmonth       ! end   month in the hydro-year
        integer startday       ! start day in the first month of the hydro-year
        integer endday         ! end   day in the last  month of the hydro-year
        integer idc            ! contineous day in a year (1-366)
        integer ihc            ! contineous hour in a year (1-366*24)
        integer start, finish  ! 0 - faulse,    1 - true
        integer start_sub, end_sub  ! for partial simulation of a basin

        real      D(nrow,ncol,nlayer)    ! depth of each UZ layer(m)
        real      k0(nrow,ncol,nlayer)   ! saturated hydraulic conductivity (mm/hr)
        real      w(nrow,ncol,nv,nlayer) ! soil moisture
        save      layer                              ! number of UZ layer
        save      D                                  ! depth of each UZ layer(m)
        save      k0                              ! saturated hydraulic conductivity (mm/hr)
        save      w                                ! soil moisture
        real    q_hillslope         ! surface runoff of one simulation 
                                  ! unit (m3/s/m, flow into river)
        real    qin(nx)             ! lateral inflow into river
        real    kg(nrow,ncol)          ! hydraulic conductivity of groundwater (m/sec) 
        real    GWcs(nrow,ncol)            ! groundwater storage coefficient changed by gaobing 20131125
        real    Cstmax(nrow,ncol,nv)          ! maximal canopy storage
        real      Cst(nrow,ncol,nv)                  ! canopy storage
        real      Sst(nrow,ncol,nv)                  ! surface storage
        real      Dgl(nrow,ncol,nv)               ! depth to groundwater level in the grid
        real    GWst(nrow,ncol,nv)              ! 
        real      prec                ! rainfall (mm/hr)
        real      temper              ! temperature(degree,hourly)
        real    snowmelt            ! snowmelt equvalent water (mm)
        real    Pnet                ! net precipitation (mm)
        real      Ep                        ! potential evaporation(mm/hr)
        real    Etr                 ! reference transpiration (mm)
        real    Es                  ! reference evaporation (mm)
        real    Eact                ! actual evaporation
        real      EfromCanopy                  ! actual evaporation from canopy storage
        real      EfromCrop                  ! actual transpiration from crop
        real      EfromSurface            ! actual evaporation from soil surface
        real    DeficitCst          ! deficit of canopy storage (mm)
        save      Cst                       ! canopy storage
        save      Sst                       ! surface storage
        save      Dgl                        ! depth to groundwater level in the grid
        save    GWst                  ! 
        integer m_grow(6,7)
        integer d_grow(6,7)  
        real      Drw(nc,nx)          ! river water depth
        real      discharge(nc,nx)    ! river flow discharge
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
        real srunoff2(nc,8800)   !mqh £¬²»Í¬²úÁ÷
        real groundoff2(nc,8800)
        real soiloff2(nc,8800)
        real soiloff_total
        real srunoff_total
        real groundoff_total
        real annual_rain        ! annual precipitation 
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
        real unbalance          ! simulated unbalance error
        save eactd                ! daily actual evapotranspiration (mm)
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
        save annual_runoff      ! annual runoff
        save annual_Eact        ! annual actual evaporation
        save annual_Ecy
        save annual_Ecp
        save annual_Ese
        save annual_Cst         ! annual canopy interception
        save annual_Ssts        ! annual surface storage(snow)
        save annual_Sstr        ! annual surface storage(rain)
        save annual_SBst        ! annual subsurface storage
        save annual_Gst         ! annual groundwater storage
        save annual_Cst0        ! initial value of annual canopy storage
        save annual_Ssts0       ! initial value of annual surface storage(snow)
        save annual_Sstr0       ! initial value of annual surface storage(rain)
        save annual_SBst0       ! initial value of annual subsurface storage
        save annual_Gst0        ! initial value of groundwater storage
        real      total_rain_d(366)            ! basin mean precipitation
        real      total_Ep_d(366)                  ! basin mean potential evaporation
        real      total_Eact_d(366)            ! basin mean actual evaporation
        real      total_runoff_d(366)            ! basin mean runoff
        save      total_rain_d
        save      total_Ep_d
        save      total_Eact_d
        save      total_runoff_d
        real evap_month(nrow,ncol)      ! monthly mean evaporation
        real soil_month(nrow,ncol)
        real runoff_month(nrow,ncol)
        real grunoff_month(nrow,ncol)
        real drunoff_month(nrow,ncol)
        real srunoff_month(nrow,ncol)
        real ecy_month(nrow,ncol)
        real ecp_month(nrow,ncol)
        real ese_month(nrow,ncol)
        save evap_month
        save ecy_month
        save ecp_month
        save ese_month
        save soil_month
        save runoff_month
        save grunoff_month
        save drunoff_month
        save srunoff_month
        integer d1_nuz
        real    d1_land_ratio, d1_slope, d1_length, d1_Ds, d1_Dg
        real    d1_Dr, d1_anik, d1_wsat, d1_wrsd, d1_wfld
        real    d1_alpha, d1_watern, d1_kg, d1_GWcs, d1_GWst, d1_Dgl
        real    d1_Drw, d1_Sst, d1_qsub,d1_dt, d1_snow
        real    d1_ssub  ! added by gaobing
        real    d1_deltz(nlayer), d1_k0(nlayer), d1_w(nlayer)
 
        character*200  para_dir    ! directory for parameters
        character*200  data_dir    ! directory for input data
        character*200  dem_dir    ! directory for input dem
        character*200  result1_dir ! directory for storing simulation result
        character*200  result2_dir ! directory for storing simulation result
        character*200  simul_dir   ! directory for storing temporal variables
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
        real    tmp,temp          !,mean , value   ! temporary variables
        integer i, j              !, k,n,m temporary variables
        integer ig,nuz                    !,ii, jj
        integer ir, ic            !, iyy, imm, idd
        real    basinarea, subarea(nc)    ! total basin area
        integer ia1, ia2  , ib1, ib2 !, ic1, ic2
        integer ii2,j1,j2,jj,imm,iidd,nday
        real    tmp1,tmp2, tmp3, D0,Ddt
        real    f
        real  para_r, para_a, y, tmpEp      !tmpinf_max,, tmpinf2,tmpinf1
        real  qin_tmp0, qin_tmp
        real water_depth, surface_n, waterhead, power, fice, area_frac
        real tmppnet
        real  rain_tmp(300,1000), Ep_tmp(300,1000), T_tmp(300,1000)
        real c1,c2,c3,c4,c5,dLAI
        integer kn
        save  subarea, fice, area_frac
        integer isub2  !added by mqh
        integer luolicun2(nc)  !ŽæŽ¢ÄÄÐ©ÊÇÂÞÀîŽåµÄÉÏÓÎÕŸ   added by mqh

        data luolicun2 /1,1,0,0,0,1,1,1,0,0,1,1,1,1,1,1,0,1,0,
     :1,1,1,1,0,0,1,1,0,1,0,1,1,
     : 1,1,1,1,0,1,0,1,1,0,0,1,0,1,1,1,1,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,
     :0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,
     :0,1,0,0,0,0,1,1,1,1,1,1,1,1,0,1,1,0,0,1,0,0,1,0,0,1,0,1,1,1,1,1,1,
     :1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/

!#######################################################################
!     NDVI-->LAI look-up table
!#######################################################################
        integer i_ndvi,ii_ndvi
        real T_NDVI(20),Shurb_LAI(20)
        real Forest_LAI(20),IRR_LAI(20),N_IRR_LAI(20)

      data T_NDVI/0.025,0.075,0.125,0.175,0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575,0.625,0.675,&
                 0.725,0.775,0.825,0.875,0.925,0.975/
      data Shurb_LAI/0,0,0.2663,0.3456,0.4357,0.5213,0.6057,0.6951,0.8028,0.9313,1.102,1.31,1.598,&
           1.932,2.466,3.426,4.638,6.328,6.328,6.328/
      data Forest_LAI/0,0,0.1516,0.1973,0.2686,0.3732,0.5034,0.6475,0.7641,0.9166,1.091,1.305,1.683,&
           2.636,3.557,4.761,5.52,6.091,6.091,6.091/
      data IRR_LAI/0,0,0.3199,0.431,0.5437,0.6574,0.7827,0.931,1.084,1.229,1.43,1.825,2.692,4.229,5.362,&
           5.903,6.606,6.606,6.606,6.606/
      data N_IRR_LAI/0,0,0.2452,0.3432,0.4451,0.5463,0.6621,0.7813,0.8868,0.9978,1.124,1.268,1.474,&
            1.739,2.738,5.349,6.062,6.543,6.543,6.543/

      common /date1  / year, month, day, hour, idc,ihc
      common /date2  / hydroyear,startmonth,startday,endmonth,endday
      common /date3  / dayinmonth
      common /simulation/  start, finish, dt, start_sub, end_sub
      common /river1    /  nsub, subbasin
      common /river2    /  nflow, dx, Dr
      common /river3    /  s0, b, roughness
      common /river4    /  psubbasin,nbasinup ! added by gaobing
      common / weather / rain_daily,tmin_daily, tmax_daily,evap_daily,snow,pre_hour
      common /LAID   / NDVI,yndvi,LAI
      common /TEST   /basinarea
      common /grid_attrib/ ngrid,grid_row,grid_col
      common /slope     /  area, length, slp, Ds, Dg
      common /soil_land /  nsoil,soiltyp,nland, landtyp,land_ratio,soil,land_use
      common /veg_para / Kcanopy,LAImax,root,anik,Sstmax,surfn,Kcrop
      common /soil_para /  wsat,wfld,wrsd,alpha,watern,ksat1,ksat2,kg,GWcs
      common /location  /  para_dir,data_dir,result1_dir,dem_dir,result2_dir,simul_dir
      common /lateralflow/ qin
      common / river_status / Drw,discharge   ! depth of river water        
      common / asc / ncols,nrows,xllcorner,yllcorner,cellsize, nodata, nnr, nnc, x0, y0, gridsize,znodata
      real snow1(nrow,ncol)
             
      integer countt2(nc),ii1   !ŒÆËãÍøžñÊý
      real soil_water1(nc),soil_water2(nc),soil_water3(nc),soil_water4(nc),soil_water5(nc),soil_water6(nc)
      save soil_water1,soil_water2,soil_water3,soil_water4,soil_water5,soil_water6
      save countt2
      real tmp19921,tmp19922,tmp19923,tmp19924,tmp19925,tmp19926
      real result_soilwater1,result_soilwater2,result_soilwater3,result_soilwater4,result_soilwater5,result_soilwater6
      integer luolicun(nc)  !ŽæŽ¢ÄÄÐ©ÊÇÂÞÀîŽåµÄÉÏÓÎÕŸ
        data luolicun /1,1,0,0,0,1,1,1,0,0,1,1,1,1,1,1,0,1,0,
     :1,1,1,1,0,0,1,1,0,1,0,1,1,
     : 1,1,1,1,0,1,0,1,1,0,0,1,0,1,1,1,1,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,
     :0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,
     :0,1,0,0,0,0,1,1,1,1,1,1,1,1,0,1,1,0,0,1,0,0,1,0,0,1,0,1,1,1,1,1,1,
     :1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/


      soil_water1(isub)=0
      soil_water2(isub)=0
      soil_water3(isub)=0
      soil_water4(isub)=0
      soil_water5(isub)=0
      soil_water6(isub)=0


      smf  = 0.1     ! snowmelting factor  !mqh,from 0.15 to 0.1
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
cc     $                                 or.ksat1(isoil).eq.0) goto 10 

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

c             f   = -alog(ksat2(isoil)/ksat1(isoil))/5.0
             f   = -alog(ksat2(ir,ic)/ksat1(ir,ic))/Ds(ir,ic)
             k0(ir,ic,j) = ksat1(ir,ic)*exp(-f*tmp)
c               print *, j, k0(ir,ic,j), D(ir,ic,j), tmp, Ds(ir,ic)
           end do
           if(abs(tmp-Ds(ir,ic)) .gt. 0.01)     print *,  
     :           'wrong in the vertical discretization', tmp, Ds(ir,ic),
     :            layer(ir,ic), D(ir,ic,layer(ir,ic))

 10      continue
 
       end do
      end if

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
     :                     'I_soil2',status='old')
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
     :                                   -Dgl(ir,ic,iland))*Gwcs(ir,ic)
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
c            sub-basin area Ãæ»ý
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
cc             epd(ir,ic,i)       = 0.0
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
cc        annual_re0   = 0.0
      end if
c  
      if( start .eq. 1 .or. (idc.eq.1 .and. ihc.eq.1)) then 
        do iflow = 1, nflow(isub)
cc          annual_re0 = annual_re0 + re_capacity(isub,iflow)
            do ig = 1, ngrid(isub,iflow)
            ir = grid_row(isub,iflow,ig)
            ic = grid_col(isub,iflow,ig)
c            isoil = soil(ir,ic)
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
c      print *,basinarea
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
c        print *,iflow,"ngrid=", ngrid(isub,iflow) 
        do ig=1,ngrid(isub,iflow)
          ir=grid_row(isub,iflow,ig)
          ic=grid_col(isub,iflow,ig)
          isoil=soil(ir,ic)
c           changged by gaobing


       
c            if(isub.ge.10 .and. isub.le.18)
c              fice = (temper+10.0)/10.0
c            fice = amax1(fice, 0.1)
c            fice = amin1(fice, 1.0)
c
c****************************************************************
c variable using for water balance calculation
c****************************************************************
c
c            raind(ir,ic,idc) = raind(ir,ic,idc) + prec
c            epd(ir,ic,idc)   = epd(ir,ic,idc)   + Ep
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
c      if(isub.ge.10 .and. isub.le.18 )          prec = 3*prec
cc      
            prec = pre_hour(ir,ic,day,hour)
c            
c         if(year == 2000)  prec = prec*1.4  !±êŒÇ£¬Ê²ÃŽÒâËŒ£¿Ô­ÀŽÊÇÓÐµÄ
cc      print *, prec, ir, ic, day, hour
             
              if(iland.eq.1) then
              raind(ir,ic,idc) = raind(ir,ic,idc) + prec
              end if
cc            if(land_ratio(ir,ic,iland).ge.0) then
cc                epd(ir,ic,idc)=epd(ir,ic,idc)+Ep*land_ratio(ir,ic,iland) 
cc            end if

            pnet           = prec
            Eact           = 0.0
            if(land_ratio(ir,ic,iland).le.0.0) goto 888


ccccccccccccccccccccccccc   added by gaobing ccccccccccccccccc
c         print *,j
         !±êŒÇ£¬ÕâÀï±ŸÀŽ¶ŒÊÇ±»×¢ÊÍµôµÄ
c        if(iland.eq.1 .or. iland.eq.10) then  !waterbody or frost snow
c               do kn = 1,4
c               LAI(ir,ic,j,kn) = 0.0
c               enddo
c              end if 
 
c          if(iland.eq.2) then
c             c4 = 0.06
c          elseif( iland .eq. 3) then
c             c4=0.2
c          elseif( iland .eq. 4) then
c              c4=0.27
c             c4=10.00
c          elseif( iland .eq. 6) then
c             c4=0.4
c          elseif( iland .eq. 7) then
c             c4=0.3
c          elseif( iland .eq. 8) then
cc             c4=0.08
c            c4=0.40
c          elseif( iland .eq. 9) then
c             c4=0.3
 
c          endif
       

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
cc              if(iland.eq.1 .or. iland.eq.19) then  !waterbody or frost snow
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
c               LAI(ir,ic,j) = 1.71*NDVI(ir,ic,j)+0.48
               LAI(ir,ic,j,kn) = 0.2*NDVI(ir,ic,j,kn)+0.2
               c4=0.4 ! changed by gaobing
              end if 

cc            if(iland.eq.4 .or. iland.eq.5) then  !forest
            if(iland.eq.4) then  !forest
                 if(ii_ndvi.eq.0 .or. ii_ndvi.eq.20) then
                 LAI(ir,ic,j,kn) = Forest_LAI(ii_ndvi)           
                    else
                 LAI(ir,ic,j,kn)=Forest_LAI(ii_ndvi)+(NDVI(ir,ic,j,kn)
     $           -T_NDVI(ii_ndvi))*(Forest_LAI(ii_ndvi+1)-Forest_LAI
     $           (ii_ndvi))/(T_NDVI(ii_ndvi+1)-T_NDVI(ii_ndvi))           
                   end if
               if(iland.eq.4) c4 = 0.2
cc               if(iland.eq.5) c4 = 0.08
              end if 

cc            if(iland.ge.6 .and. iland.le.12) then  !irrigated-cropland

            if(iland.ge.5) then  ! sparse vegetion
c               if(ii_ndvi.eq.0 .or. ii_ndvi.eq.20) then
c                 LAI(ir,ic,j,kn)=IRR_LAI(ii_ndvi)           
c                    else
c                 LAI(ir,ic,j,kn)=IRR_LAI(ii_ndvi)+(NDVI(ir,ic,j,kn)-
c     $           T_NDVI(ii_ndvi))*(IRR_LAI(ii_ndvi+1)-
c    $           IRR_LAI(ii_ndvi))/(T_NDVI(ii_ndvi+1)-T_NDVI(ii_ndvi))           
c                   end if
c               c4 = 0.45
           LAI(ir,ic,j,kn)=1.71*NDVI(ir,ic,j,kn)+0.48
           c4 = 0.3
              end if 

cc            if(iland.eq.13) then  !upland
            if(iland.eq.6) then  !upland
               if(ii_ndvi.eq.0 .or. ii_ndvi.eq.20) then
                 LAI(ir,ic,j,kn)=N_IRR_LAI(ii_ndvi)           
                    else
                 LAI(ir,ic,j,kn)=N_IRR_LAI(ii_ndvi)+(NDVI(ir,ic,j,kn)-
     $           T_NDVI(ii_ndvi))*(N_IRR_LAI(ii_ndvi+1)-
     $           N_IRR_LAI(ii_ndvi))/(T_NDVI(ii_ndvi+1)-T_NDVI(ii_ndvi))           
                   end if
              end if 

cc            if(iland.ge.14 .and. iland.le.16) then  !grassland

            if(iland.eq.7) then  !grassland
               LAI(ir,ic,j,kn)=1.71*NDVI(ir,ic,j,kn)+0.48
c                   if(isub .le. 23) then  !Taking special value of LAI for the high plateau
c                 LAI(ir,ic,j)=-alog((1.0-NDVI(ir,ic,j)/0.915)
c     $                            /0.83)/0.96
c                 if(LAI(ir,ic,j) .lt. 0.0)  LAI(ir,ic,j) = 0.0
c                 if(LAI(ir,ic,j) .gt. 2.18) LAI(ir,ic,j) = 2.18
c                   end if 
cc               if(iland.eq.14) c4 = 0.6
cc               if(iland.eq.15) c4 = 0.3
cc               if(iland.eq.16) c4 = 0.1
               if(iland.eq.7) c4 = 0.3
c              if(iland.eq.7) c4 = 0.6
              end if 

cc            if(iland.eq.17) then  !shrub

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

cc            if(iland.eq.18) then  !wetland

            if(iland.eq.9) then  !wetland
               LAI(ir,ic,j,kn) = 1.71*NDVI(ir,ic,j,kn)+0.48
              end if 
c
            if(LAI(ir,ic,j,kn) .lt. 0.0) LAI(ir,ic,j,kn) = 0.0
              if(LAI(ir,ic,j,kn) .gt. LAImax(iland)) LAImax(iland) = 
     :                                              LAI(ir,ic,j,kn)
      enddo ! kn
c           c4=0.2           ! changed by gaobing
          end do


c
c      print *, iland, NDVI(ir,ic,month),LAI(ir,ic,month), month, ir, ic
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

cc      111      ÉŸ³ýÓÚ2006.4.14       !skip d_grow

c
c****************************************************************
c          interception  
c****************************************************************
c   
c              Cstmax(ir,ic,iland) = 0.20*dLAI   !LAI(ir,ic,month)
             Cstmax(ir,ic,iland) = 0.10*dLAI
cc            if(iland.eq.4 .or. iland.eq.5 .or.iland.eq.17) then !forest or shurb
            if(iland.eq.4 .or. iland.eq.8) then !forest or shurb
c                Cstmax(ir,ic,iland) = 0.25*dLAI !LAI(ir,ic,month)
             Cstmax(ir,ic,iland) = 1.00*dLAI 
              end if

c              Cstmax(ir,ic,iland) = 0.2*LAI(ir,ic,month)

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
c             calculate snow melt   
c****************************************************************
c             
c            if(snow(ir,ic,iland).gt.0.0) then
c                print *, 'ir, ic, iland, snow(ir,ic,iland), temper, pnet'
c              print *, ir, ic, iland, snow(ir,ic,iland), temper, pnet
c             pause
c              endif
                 
              if(temper .le. 1.0) then
                  snow(ir,ic,iland) = snow(ir,ic,iland) + pnet
                  pnet=0.0 

c          added by xujijun      - start       
c        if(snow(ir,ic,iland).ge.100) snow(ir,ic,iland) = 100
c          added by xujijun      - end       
            end if
              snowmelt=0.0
c              if(month .ge. 2 .and. 

c              print *, 'isub, month, ir, ic, pnet'
c              print *, isub, month, ir, ic, pnet
c           changed by gaobing
              if(temper .gt. 0.0 
     :                      .and. snow(ir,ic,iland).gt.0.0) then
c          added by xujijun      - gaobing
              if(month.ge.3.and.month.le.10) smf= 0.15
            if(month.ge.5.and.month.le.8)  smf = 0.15
           
            

c          added by xujijun      - end       
              snowmelt = (smf+pnet/20.0) * (temper-1.5)
c               snowmelt = (smf+pnet/10.0) * temper
              snowmelt = amin1(snowmelt,snow(ir,ic,iland))
                      
                  snow(ir,ic,iland)=snow(ir,ic,iland)-snowmelt             
            
                  pnet=pnet+snowmelt 
               
            
c              if(snow(ir,ic,iland).gt.0.0) then
c               print *, 'smf, snowmelt, pnet, snow' 
c               print *, smf, snowmelt, pnet, snow(ir,ic,iland)
c            endif
           
              end if
            
c            if(snow(ir,ic,iland) > 60) then  ! added vy gaobing 
c                     snow(ir,ic,iland)=60
c                  snowmelt=snowmelt+snow(ir,ic,iland)-60
c                  pnet=pnet+snow(ir,ic,iland)-60
c               endif

c        added by gaobing      - end
            Sst(ir,ic,iland) = Sst(ir,ic,iland) + pnet
c            if(pnet.gt.0) print *, "Sst=,p",Sst(ir,ic,iland),pnet,ir,ic
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

            c3 =  0.23+0.1*(dLAI/LAImax(iland))**0.5            ! 
              c5 =      0.23
            if (luolicun2(isub).eq.1) then
            c1=c1*0.5
            c2=c2*0.5
            c3=c3*0.5
            c5=c5*0.4
            endif
              
                      


c            if(temper .le. 0.0) Ep=0.0
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

c              tmp = 1.0
c              if(isub.le.23) then
c                tmp = fice*0.8
c                if(year.gt.1960) tmp = fice*0.82
c                if(year.gt.1970) tmp = fice*0.84
c                if(year.gt.1980) tmp = fice*0.86
c                if(year.gt.1990) tmp = fice*0.88
c            endif
c              if(isub.ge.10.and. isub.le.69
c     :          .and.iland.ge.6.and.iland.le.12) then
c               tmp =1.05 + 0.05*sin((month-3)*3.1415926/6.0
c     :                - min(2,max(month-6,0))*3.1415926/6.0)
c            endif
c              if(isub.gt.69.and.iland.ge.6.and.iland.le.12) then
c               tmp =1.08 + 0.08*sin((month-3)*3.1415926/6.0
c     :                - min(2,max(month-6,0))*3.1415926/6.0)
c            endif
c
c              Etr = tmp * Etr
c              Es  = tmp * Es

         do kn=1,3
      if(kcrop(iland).gt.1.or.(LAI(ir,ic,month,kn)-LAImax(iland)).gt.1e-6)
     :  print *,"kcrop out of range", kcrop(iland), LAI(ir,ic,month,kn), 
     :                                LAImax(iland), iland
         enddo
!              if(prec. gt. 0.01) goto 290   !mqh

            if(Etr .lt. 0.1E-9) goto 270

              EfromCanopy = amin1(Etr, Cst(ir,ic,iland)) 
            Cst(ir,ic,iland) = 
     :                  Cst(ir,ic,iland) - EfromCanopy
c            if(prec .gt. 0.0001) then
c           print *,"Canopy",EfromCanopy,Etr,Es,Cst(ir,ic,iland),idc,ihc
c            end if
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
c            ÓÉmqh×¢ÊÍ   ±êŒÇ1
c              if(w(ir,ic,iland,j).lt.wrsd(ir,ic)+0.1E-6) then
c                print *,"crop evap from w",w(ir,ic,iland,j),j,'11111'
c                pause
c                  end if
              EfromCrop = EfromCrop + tmp
c            if(prec .gt. 0.0001) then
c           print *,"Crop,w",EfromCrop,Etr,Es,w(ir,ic,iland,j),idc
c            end if
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
     :              print *,'wrong in EfromSurface',
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
c           ÓÉmqh×¢ÊÍ  ±êŒÇ2
c            if(w(ir,ic,iland,1).lt.wrsd(ir,ic)+0.1E-6) then
c              print *,"wrong in calculate surface evap",ir,ic,iland,'1111'
c              pause
c              end if
            EfromSurface = EfromSurface + tmp
                  
 290        continue
c


c              print *, 'es, sst', Es, Sst(ir,ic,iland),ir,ic
c            pause

c****************************************************************
c variable using for water balance calculation
c****************************************************************
c
            eactd(ir,ic,idc) = eactd(ir,ic,idc) + 
     :                         (EfromCanopy+EfromCrop+EfromSurface) *
     :                          land_ratio(ir,ic,iland)            ! mm
      ecy(ir,ic,idc)   =  ecy(ir,ic,idc) + EfromCanopy 
     :      *land_ratio(ir,ic,iland) 
      ecp(ir,ic,idc)   =  ecp(ir,ic,idc) + EfromCrop 
     :      *land_ratio(ir,ic,iland) 
      ese(ir,ic,idc)   =  ese(ir,ic,idc) + EfromSurface
     :      *land_ratio(ir,ic,iland)                       ! added by gaobing
            eactd(ir,ic,idc) = eactd(ir,ic,idc) + Eact * 
     :                              land_ratio(ir,ic,iland)        ! mm
      !±êŒÇ
c****************************************************************
c
c            print *,"eactd=",eactd(isub,iflow,idc),EfromCanopy
c     :            ,EfromCrop,EfromSurface,idc,ihc,ir,ic,iland
c            print *,"rain=",prec,"  Ep=",EfromSurface,idc,ihc
c
 888            continue
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
c              endif

          qin_tmp=0.0

          d1_slope  = slp(ir,ic)           ! m/m
          d1_length = length(ir,ic)       ! meter           *1.3£¬added by mqh ,2014.4.17
          d1_Ds     = Ds(ir,ic)            ! meter
          d1_Dg     = Dg(ir,ic)                 ! meter
          d1_Dr     = Dr(isub, iflow)             ! meter
          d1_Drw    = Drw(isub, iflow)       ! meter
          !mqh  added by mqh
c          if (luolicun(isub).eq.0) kg(ir,ic)=kg(ir,ic)*2
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

c              if(d1_wrsd-d1_w(j).gt.0.005) 
c     :            print *,"d1_w(j) < wrsd", d1_w(j), d1_wrsd
c     :                ,ir,ic,j,idc,ihc

c              if(d1_w(j)-d1_wsat .gt. 0.005) 
c     :          print *,"d1_w(j)>wsat", d1_w(j),d1_wsat,j

            end do
            
              d1_Sst = Sst(ir,ic,iland)/1000.0
            d1_Dgl = Dgl(ir,ic,iland)           
            d1_GWst= GWst(ir,ic,iland)
              d1_snow= snow(ir,ic,iland)
            
cc            if(month.gt.3) then
c            if(temper.gt.1.5.and. d1_snow.gt.0.and.snowmelt.gt.0) then
c            print *,' temper4, sst, d1_Dgl, d1_GWst, d1_snow'
c              print *,temper, Sst(ir,ic,iland), d1_Dgl, d1_GWst, d1_snow
c            print *, pnet, snowmelt
c            endif
cc            endif

c              print *,"input data", d1_slope,d1_length,d1_Ds,d1_Dg,d1_Dr,
c     :                                      d1_Drw,d1_kg,d1_sst,d1_Dgl,d1_GWst
c         print * ,'ok6111 of hillslope'
c                 print *, 'layer(ir,ic), d1_nuz=',layer(ir,ic), d1_nuz  
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
c            if(temper.gt.1.5.and. d1_snow.gt.0.and.snowmelt.gt.0) then
c            print *,' temper5, sst, d1_Dgl, d1_GWst, d1_snow'
c              print *,temper, Sst(ir,ic,iland), d1_Dgl, d1_GWst, d1_snow
c            print *, pnet, snowmelt
c            endif  
            Dgl(ir,ic,iland) = d1_Dgl
            GWst(ir,ic,iland)= d1_GWst
c               print *,"output data",d1_slope,d1_length,d1_Ds,d1_Dg,d1_Dr
c     :     ,d1_Drw,d1_kg,d1_sst,d1_Dgl,d1_GWst,land_ratio(ir,ic,iland)

            qin_tmp = qin_tmp + d1_qsub*land_ratio(ir,ic,iland)
c
c          print * ,'ok612 of hillslope'

c****************************************************************
c variable using for water balance calculation
c****************************************************************
c
           groundoff2(isub,ihc)=groundoff2(isub,ihc)+ area(ir,ic)* 
     :      (d1_qsub-d1_ssub)*land_ratio(ir,ic,iland)/length(ir,ic)         ! m3/s                    !mqh,ŒÆËã²»Í¬³É·Ö²úÁ÷
           soiloff2(isub,ihc)=soiloff2(isub,ihc)+ area(ir,ic)* 
     :      d1_ssub*land_ratio(ir,ic,iland)/length(ir,ic)         ! m3/s                    !mqh,ŒÆËã²»Í¬³É·Ö²úÁ÷
            groundoff(ir,ic,idc) = groundoff(ir,ic,idc) + 1000.0* 
     :                             dt*d1_qsub*land_ratio(ir,ic,iland)/
     :                             length(ir,ic)         ! mm
            soiloff(ir,ic,idc)   = soiloff(ir,ic,idc)+1000.0*d1_ssub* 
     :                                 land_ratio(ir,ic,iland)
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
            
cc      ÐÞžÄÓÚ2006.4.16
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
            if( water_depth .le. 0.01 !) goto 1800                              ! ÊéÇ©4
     :         .or. Drw(isub,iflow).ge. 1.0*Dr(isub,iflow) ) goto 1800

            Sst(ir,ic,iland) = Sst(ir,ic,iland)-water_depth
            water_depth = 0.001 * water_depth          !in meter, surface runoff

            surface_n = surfn(iland) ! * sqrt(water_depth*1000.)
c     :                  * sqrt(LAI(ir,ic,month)/LAImax(iland))
cc
            waterhead = slp(ir,ic) ! + water_depth/length(ir,ic)
c            if(waterhead .lt. 0.1E-8) waterhead = 0.1E-8
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


        
c            update surface storage
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
!            srunoff(ir,ic,idc) = srunoff(ir,ic,idc) + 1000.0*
!     :                           q_hillslope*land_ratio(ir,ic,iland)/
!     :                           length(ir,ic)         ! mm
c            runoffd(ir,ic,idc)   = runoffd(ir,ic,idc) + 
c     :                             srunoff(ir,ic,idc)  ! mm
c
c****************************************************************
c
1800          continue
   
            
            qin_tmp  = qin_tmp *area(ir,ic)/length(ir,ic)! m3/s, one flow-interval
          qin_tmp0 = qin_tmp0 + qin_tmp
c
          runoffd(ir,ic,idc) = runoffd(ir,ic,idc) + 
     :                       1000.0*qin_tmp*dt/area(ir,ic)

c          print * ,'ok614 of hillslope'

        end do
c

c          print * ,'ok61 of hillslope'

c      ----------------------------------------------------------------------   
c       lateral inflow to the river
c
        qin(iflow) = qin_tmp0 ! m3/s, total lateral inflow of one flow-interval

cc      ÐÞžÄÓÚ2006.1.15
cc      print *, qin(iflow)
c
c
cc      5555      ÉŸ³ýÓÚ2006.4.14       !skip irrigation loop

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
              if(tmp .ge. root(iland))     goto 4267  !Ô­ÀŽÊÇ4266

            end do

cccccccccccccccccccccccccccccccccccccccccc            
            !added by mqh
4267                tmp19921 = 0.0
                tmp19922 = 0.0
                tmp19923 = 0.0
                tmp19924 = 0.0
                tmp19925 = 0.0
                tmp19926 = 0.0
       
            do j = 1, 1    !ÊÔÊÔÑ¡»ù²ãÍÁÈÀ±ÈœÏºÏÊÊ
              tmp19921 = tmp19921 + D(ir,ic,j)     
             
            end do
            do j = 1, 2    !ÊÔÊÔÑ¡»ù²ãÍÁÈÀ±ÈœÏºÏÊÊ
              tmp19922 = tmp19922 + D(ir,ic,j)     
             
            end do
            do j = 1, 3    !ÊÔÊÔÑ¡»ù²ãÍÁÈÀ±ÈœÏºÏÊÊ
              tmp19923 = tmp19923 + D(ir,ic,j)     
             
            end do
            do j = 1, 4    !ÊÔÊÔÑ¡»ù²ãÍÁÈÀ±ÈœÏºÏÊÊ
              tmp19924 = tmp19924 + D(ir,ic,j)     
             
            end do
            do j = 1, 5    !ÊÔÊÔÑ¡»ù²ãÍÁÈÀ±ÈœÏºÏÊÊ
              tmp19925 = tmp19925 + D(ir,ic,j)     
             
            end do
            do j = 1, 6    !ÊÔÊÔÑ¡»ù²ãÍÁÈÀ±ÈœÏºÏÊÊ
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
1999      continue


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

!cccccc       open(33,file = result2_dir(ia1:ia2)//  ! changed by xujijun
!     $                 'ecy_'//ch4//ch2//'.asc')        
!          write(33,'(A, i16)') ncols,nnc
!          write(33,'(A, i16)') nrows,nnr
!          write(33,'(A, f16.3)') xllcorner,x0
!          write(33,'(A, f16.3)') yllcorner,y0
!          write(33,'(A, f16.0)') cellsize,gridsize
!          write(33,'(A, f8.0)') nodata,znodata
!        do ir=1 , nrow
!          do ic=1 , ncol
!            if(area(ir,ic).lt.0.0) ecy_month(ir,ic) = -9999.0
!          enddo
!        enddo
!        do ir=1,nrow
!             write(33,'(295f13.4)') (ecy_month(ir,ic), ic=1,ncol)
!        end do
!        close(33)

ccccc      open(33,file = result2_dir(ia1:ia2)//  ! changed by xujijun  !mqh7.18
c     $                 'ecp_'//ch4//ch2//'.asc')        
c          write(33,'(A, i16)') ncols,nnc
c          write(33,'(A, i16)') nrows,nnr
c          write(33,'(A, f16.3)') xllcorner,x0
c          write(33,'(A, f16.3)') yllcorner,y0
c          write(33,'(A, f16.0)') cellsize,gridsize
c          write(33,'(A, f8.0)') nodata,znodata
c       do ir=1 , nrow
c          do ic=1 , ncol
c            if(area(ir,ic).lt.0.0) ecp_month(ir,ic) = -9999.0
c          enddo
c       enddo
c        do ir=1,nrow
c             write(33,'(295f13.4)') (ecp_month(ir,ic), ic=1,ncol)
c       end do
c      close(33)

cccc      open(33,file = result2_dir(ia1:ia2)//  ! changed by xujijun    !mqh7.18
c     $                 'ese_'//ch4//ch2//'.asc')        
c          write(33,'(A, i16)') ncols,nnc
c          write(33,'(A, i16)') nrows,nnr
c          write(33,'(A, f16.3)') xllcorner,x0
c          write(33,'(A, f16.3)') yllcorner,y0
c          write(33,'(A, f16.0)') cellsize,gridsize
c          write(33,'(A, f8.0)') nodata,znodata
c       do ir=1 , nrow
c          do ic=1 , ncol
c            if(area(ir,ic).lt.0.0) ese_month(ir,ic) = -9999.0
c          enddo
c       enddo
c        do ir=1,nrow
c             write(33,'(295f13.4)') (ese_month(ir,ic), ic=1,ncol)
c       end do
c        close(33)

ccccc      open(33,file = result2_dir(ia1:ia2)//  ! changed by xujijun       !mqh7.18
c     $                 'snow_'//ch4//ch2//'.asc')        
c          write(33,'(A, i16)') ncols,nnc
c          write(33,'(A, i16)') nrows,nnr
c          write(33,'(A, f16.3)') xllcorner,x0
c          write(33,'(A, f16.3)') yllcorner,y0
c          write(33,'(A, f16.0)') cellsize,gridsize
c          write(33,'(A, f8.0)') nodata,znodata
c       do ir=1 , nrow
c          do ic=1 , ncol
c            snow1(ir,ic) = 0.0
c              do i = 1, nland
c           snow1(ir,ic) = snow1(ir,ic)+snow(ir,ic,i)*land_ratio(ir,ic,i)
c            enddo
c            if(area(ir,ic).lt.0.0) snow1(ir,ic) = -9999.0
c          enddo
c        enddo
c        do ir=1,nrow
c             write(33,'(295f13.4)') (snow1(ir,ic), ic=1,ncol)
c        end do
cccccc        close(33)

c       print *,'soil.'//month//year
c        open(33,file=result2_dir(ia1:ia2)//
c     $       'soil.'//ch2//ch4,status='replace')

cccccc        open(33,file = result2_dir(ia1:ia2)//  ! changed by xujijun   !mqh7.18
c     $                 'runoff_'//ch4//ch2//'.asc')        
c          write(33,'(A, i16)') ncols,nnc
c          write(33,'(A, i16)') nrows,nnr
c          write(33,'(A, f16.3)') xllcorner,x0
c          write(33,'(A, f16.3)') yllcorner,y0
c          write(33,'(A, f16.0)') cellsize,gridsize
c          write(33,'(A, f8.0)') nodata,znodata
c       do ir=1 , nrow
c          do ic=1 , ncol
c            if(area(ir,ic).lt.0.0) runoff_month(ir,ic) = -9999.0
c          enddo
c        enddo
c        do ir=1,nrow
c             write(33,'(10f13.4)') (runoff_month(ir,ic), ic=1,ncol)
c        end do
c        close(33)

cccc      open(33,file = result2_dir(ia1:ia2)//  ! changed by xujijun   !mqh7.18
c     $                 'grunoff_'//ch4//ch2//'.asc')        
c          write(33,'(A, i16)') ncols,nnc
c          write(33,'(A, i16)') nrows,nnr
c          write(33,'(A, f16.3)') xllcorner,x0
c          write(33,'(A, f16.3)') yllcorner,y0
c          write(33,'(A, f16.0)') cellsize,gridsize
c          write(33,'(A, f8.0)') nodata,znodata
c       do ir=1 , nrow
c          do ic=1 , ncol
c            if(area(ir,ic).lt.0.0) grunoff_month(ir,ic) = -9999.0
c          enddo
c        enddo
c        do ir=1,nrow
c             write(33,'(10f13.4)') (grunoff_month(ir,ic), ic=1,ncol)
c        end do
c        close(33)

cccc      open(33,file = result2_dir(ia1:ia2)//  ! changed by xujijun   !mqh7.18
c     $                 'srunoff_'//ch4//ch2//'.asc')        
c          write(33,'(A, i16)') ncols,nnc
c          write(33,'(A, i16)') nrows,nnr
c          write(33,'(A, f16.3)') xllcorner,x0
c          write(33,'(A, f16.3)') yllcorner,y0
c          write(33,'(A, f16.0)') cellsize,gridsize
c          write(33,'(A, f8.0)') nodata,znodata
c       do ir=1 , nrow
c          do ic=1 , ncol
c           if(area(ir,ic).lt.0.0) srunoff_month(ir,ic) = -9999.0
c          enddo
c       enddo
c        do ir=1,nrow
c             write(33,'(10f13.4)') (srunoff_month(ir,ic), ic=1,ncol)
c        end do
c        close(33)

ccc      open(33,file = result2_dir(ia1:ia2)//  ! changed by xujijun   !mqh
c     $                 'drunoff_'//ch4//ch2//'.asc')        
c          write(33,'(A, i16)') ncols,nnc
c          write(33,'(A, i16)') nrows,nnr
c          write(33,'(A, f16.3)') xllcorner,x0
c          write(33,'(A, f16.3)') yllcorner,y0
c          write(33,'(A, f16.0)') cellsize,gridsize
c          write(33,'(A, f8.0)') nodata,znodata
c       do ir=1 , nrow
c          do ic=1 , ncol
c            if(area(ir,ic).lt.0.0) drunoff_month(ir,ic) = -9999.0
c          enddo
c        enddo
c        do ir=1,nrow
c             write(33,'(10f13.4)') (drunoff_month(ir,ic), ic=1,ncol)
c       end do
c        close(33)



!        open(33,file=result2_dir(ia1:ia2)//    ! changed by xujijun
!     $       'soil_'//ch4//ch2//'.asc')
!          write(33,'(A, i16)') ncols,nnc
!          write(33,'(A, i16)') nrows,nnr
!          write(33,'(A, f16.3)') xllcorner,x0
!          write(33,'(A, f16.3)') yllcorner,y0
!          write(33,'(A, f16.0)') cellsize,gridsize
!          write(33,'(A, f8.0)') nodata,znodata
!        do ir=1 , nrow
!          do ic=1 , ncol
!            if(area(ir,ic).lt.0.0) soil_month(ir,ic) = -9999.0
!          enddo
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
!            ecy_month(ir,ic) = 0.0
!            ecp_month(ir,ic) = 0.0
!            ese_month(ir,ic) = 0.0
!            drunoff_month(ir,ic) = 0.0
!            grunoff_month(ir,ic) = 0.0
!            srunoff_month(ir,ic) = 0.0
!          end do
!        end do
      end if
      

  
c        endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Êä³ö²»Í¬³É·Ö²úÁ÷»úÖÆµÄŽóÐ¡£¬µ¥Î»£ºmm
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
!            elseif((hydroyear.gt.2000).or.(hydroyear.eq.2000.and.ihc.gt.1)) then
!      open(69,file=result2_dir(ia1:ia2)//'su_so_g_runoff.dat',
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
c              do ig = 1, ngrid(ii2,iflow)
c                ir = grid_row(ii2,iflow,ig)
c                ic = grid_col(ii2,iflow,ig)
c                if(area(ir,ic) .gt. 0.0) then
c                  tmp1 = tmp1 + raind(ir,ic,jj)*area(ir,ic)                  
c               endif
c              end do
c            end do
c                  tmp2 = tmp2 + subarea(ii2)
c          end do
c          tmp1 = tmp1 / tmp2
c        
c          write(14,131) i,tmp1, tmp2
c131      format(2x,i6, 2f16.3) 
c          end do
c           close(14)


c      open(14, file = simul_dir(ia1:ia2)//       ! changed  by xujijun  
c     $                    'rain-zhamashike'//ch4//'.txt')

c        do jj = 1, idc
c          tmp1 = 0.0
c          tmp2 = 0.0
c          do ii2 = 1, 90
c            do iflow = 1, nflow(ii2)
c              do ig = 1, ngrid(ii2,iflow)
c                ir = grid_row(ii2,iflow,ig)
c                ic = grid_col(ii2,iflow,ig)
c                if(area(ir,ic) .gt. 0.0) then
c                  tmp1 = tmp1 + raind(ir,ic,jj)*area(ir,ic)                  
c                endif
c              end do
c            end do
c                  tmp2 = tmp2 + subarea(ii2)
c          end do
c          tmp1 = tmp1 / tmp2
        
c          write(14,131) i,tmp1, tmp2

c          end do
c           close(14)
c         end of added

        endif
      end if
c
c       if(isub==124 .and. day == 31 ) print *, 'ok hillslopeend'
c      if(day .eq. dayinmonth(month) .and. hour.eq.24) then
c      if(month.eq.endmonth .and. day.eq.endday .and. hour.eq.24) then
        if(isub.eq.start_sub) then
          annual_rain  = 0.0
cc          annual_irrigation=0.0
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
cc          annual_re    = 0.0
        endif
         
        do iflow = 1, nflow(isub)
c           annual_re = annual_re + re_capacity(isub,iflow)
            do ig=1,ngrid(isub,iflow)
            ir=grid_row(isub,iflow,ig)
            ic=grid_col(isub,iflow,ig)
            isoil=soil(ir,ic)
            do i = 1, idc         !366
              annual_runoff = annual_runoff + runoffd(ir,ic,i)*
     :                        area(ir,ic)                           ! mm
              annual_Eact   = annual_Eact + eactd(ir,ic,i)*
     :                        area(ir,ic)                                ! mm
cc              annual_Ecy   = annual_Ecy + ecy(ir,ic,i)*
cc     :                        area(ir,ic)                                ! mm
cc              annual_Ecp   = annual_Ecp + ecp(ir,ic,i)*
cc     :                        area(ir,ic)                                ! mm
cc              annual_Ese   = annual_Ese + ese(ir,ic,i)*
cc     :                        area(ir,ic)                                ! mm
             annual_rain   = annual_rain + raind(ir,ic,i)*           !rain_daily
     :                        area(ir,ic)                                 ! mm
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
     :                                 *land_ratio(ir,ic,iland)             !in m3
              do j = 1 ,layer(ir,ic)
                    annual_SBst = annual_SBst + w(ir,ic,iland,j) *
     :              D(ir,ic,j)*land_ratio(ir,ic,iland)*area(ir,ic)     !in m3
              end do
2240        continue
          end do
        end do

        if(isub .eq. end_sub) then
c          print *,"annual_SBst=", annual_SBst              
          annual_rain  = annual_rain / basinarea
cc          annual_irrigation=annual_irrigation/basinarea
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
cc          annual_re   = 1000.0*(annual_re-annual_re0)/basinarea
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
c       ÓÉmqh×¢ÊÍµô
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
c        open(9,file = simul_dir(ia1:ia2)//
c     :                subbasin(isub)//'I_soil',status='replace')

        open(9,file = simul_dir(ia1:ia2)//      ! changed by xujijun
     :                subbasin(isub)//'I_soil2')

c        print *, simul_dir(ia1:ia2)//subbasin(isub)//'I_soil'
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
c        print *, simul_dir(ia1:ia2)//subbasin(isub)//'I_soil'
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
      return
      end subroutine 
