













    subroutine hillslope_model(isub)
      implicit none

        integer(i4)::isub              ! sub-basin number
        integer(i4)::iflow             ! flow-interval number     
        integer(i4)::isoil                    ! soil number
        integer(i4):: ir,ic            
        integer(i4):: ig                   
        integer(i4):: iland   
        
        real(r8)::Cstmax               ! maximal canopy storage           
      
      
      




        integer(i4):: iland                    ! land-use number
        real(r8)::    kmoist                ! evaporation coefficient of soil mositure  
        real(r8)::    soil_con_f       ! soil water conservation factor
        real(r8)::    detension        ! surface detension (mm)
        real(r8)::    q_hillslope         ! surface runoff of one simulation  ! unit (m3/s/m, flow into river)

        real(r8)::      prec                ! rainfall (mm/hr)
        real(r8)::      temper              ! temperature(degree,hourly)
        real(r8)::    snowmelt            ! snowmelt equvalent water (mm)
        real(r8)::    Pnet                ! net precipitation (mm)
        real(r8)::      Ep                        ! potential evaporation(mm/hr)
        real(r8)::    Etr                 ! reference transpiration (mm)
        real(r8)::    Es                  ! reference evaporation (mm)
        real(r8)::    Eact                ! actual evaporation
        real(r8)::      EfromCanopy                  ! actual evaporation from canopy storage
        real(r8)::      EfromCrop                  ! actual transpiration from crop
        real(r8)::      EfromSurface            ! actual evaporation from soil surface
        real(r8)::    DeficitCst          ! deficit of canopy storage (mm)
        real(r8):: unbalance          ! simulated unbalance error


        integer(i4)::d1_nuz
        real(r8)::    d1_land_ratio, d1_slope, d1_length, d1_Ds, d1_Dg
        real(r8)::    d1_Dr, d1_anik, d1_wsat, d1_wrsd, d1_wfld
        real(r8)::    d1_alpha, d1_watern, d1_kg, d1_GWcs, d1_GWst, d1_Dgl
        real(r8)::    d1_Drw, d1_Sst, d1_qsub,d1_dt, d1_snow
        real(r8)::    d1_ssub  ! added by gaobing
        real(r8)::    d1_deltz(nlayer), d1_k0(nlayer), d1_w(nlayer)

        character*2::    ch2        ! a 2-character variable
        character*4::    ch4        ! a 4-character variable
        real(r8):: tmp,temp          !,mean , value   ! temporary variables
        integer(i4):: i, j              !, k,n,m temporary variables



        integer(i4):: ia1, ia2  , ib1, ib2 !, ic1, ic2
        integer(i4):: ii2,j1,j2,jj,imm,iidd,nday
        real(r8)::    tmp1,tmp2, tmp3, D0,Ddt
        real(r8)::    f
        real(r8)::  para_r, para_a, y, tmpEp      !tmpinf_max,, tmpinf2,tmpinf1
        real(r8)::  qin_tmp0, qin_tmp
        real(r8):: water_depth, surface_n, waterhead, power
        real(r8):: tmppnet
        real(r8)::  rain_tmp(300,1000), Ep_tmp(300,1000), T_tmp(300,1000)
        real(r8):: c1,c2,c3,c4,c5,dLAI
        integer(i4):: kn
        integer(i4):: i_ndvi,ii_ndvi
        real(r8):: T_NDVI(20),Shurb_LAI(20)
        real(r8):: Forest_LAI(20),IRR_LAI(20),N_IRR_LAI(20)

        integer(i4):: luolicun2(nc)  !ŽæŽ¢ÄÄÐ©ÊÇÂÞÀîŽåµÄÉÏÓÎÕŸ   added by mqh
        integer(i4)::ii1   !ŒÆËãÍøžñÊý


        real(r8):: tmp19921,tmp19922,tmp19923,tmp19924,tmp19925,tmp19926
        real(r8):: result_soilwater1,result_soilwater2,result_soilwater3,result_soilwater4,result_soilwater5,result_soilwater6
        integer(i4):: luolicun(nc)  !ŽæŽ¢ÄÄÐ©ÊÇÂÞÀîŽåµÄÉÏÓÎÕŸ


!     NDVI-->LAI look-up table
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







!UZ, GW and surface routing
      srunoff2(isub,ihc)=0
      groundoff2(isub,ihc)=0
      soiloff2(isub,ihc)=0  !mqh
      do 1999 iflow=1, nflow(isub)
        qin_tmp0 = 0.0
        do ig=1,ngrid(isub,iflow)
          ir=grid_row(isub,iflow,ig)
          ic=grid_col(isub,iflow,ig)
          isoil=soil(ir,ic)
          fice = 1.00
          fice=(temper+7.5)/7.5
          fice=amax1(fice, 0.5)
          fice=amin1(fice , 1.0)
          qin_tmp=0.0
          d1_slope  = slp(ir,ic)           ! m/m
          d1_length = length(ir,ic)       ! meter           *1.3£¬added by mqh ,2014.4.17
          d1_Ds     = Ds(ir,ic)            ! meter
          d1_Dg     = Dg(ir,ic)                 ! meter
          d1_Dr     = Dr(isub, iflow)             ! meter
          d1_Drw    = Drw(isub, iflow)       ! meter
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
          do 1800 iland = 1, nland
            if(land_ratio(ir,ic,iland) .le. 0.0) goto 1800
            d1_land_ratio = land_ratio(ir,ic,iland)
            d1_anik = anik(iland)
            if( d1_anik .lt. 0.0)  print *, "wrong anik..." 
            if( d1_anik .lt. 1.0)  print *, "wrong anik",anik(iland) 
            do j = 1, d1_nuz
              d1_w(j) = w(ir,ic,iland,j)
            end do
            d1_Sst = Sst(ir,ic,iland)/1000.0
            d1_Dgl = Dgl(ir,ic,iland)           
            d1_GWst= GWst(ir,ic,iland)
            d1_snow= snow(ir,ic,iland)
            call runoff(isub, month, d1_snow, temper,d1_land_ratio,d1_slope, d1_length, d1_Ds,&
                 d1_Dg ,d1_Dr,d1_anik,d1_wsat,d1_wrsd, d1_wfld,d1_alpha ,&
                 d1_watern, d1_nuz,d1_deltz, d1_k0,d1_w,d1_kg,d1_GWcs,d1_GWst,d1_Dgl,d1_Drw,d1_Sst,d1_dt,d1_qsub,d1_ssub)

            do j = 1, d1_nuz
              w(ir,ic,iland,j) = d1_w(j)
            end do
            Sst(ir,ic,iland) = d1_Sst*1000.0
            Dgl(ir,ic,iland) = d1_Dgl
            GWst(ir,ic,iland)= d1_GWst
            qin_tmp = qin_tmp + d1_qsub*land_ratio(ir,ic,iland)

            groundoff2(isub,ihc)=groundoff2(isub,ihc)+ area(ir,ic)* (d1_qsub-d1_ssub)*land_ratio(ir,ic,iland)/length(ir,ic)         ! m3/s    
            soiloff2(isub,ihc)=soiloff2(isub,ihc)+ area(ir,ic)* d1_ssub*land_ratio(ir,ic,iland)/length(ir,ic)         ! m3/s          
            groundoff(ir,ic,idc) = groundoff(ir,ic,idc) + 1000.0* dt*d1_qsub*land_ratio(ir,ic,iland)/length(ir,ic)         ! mm
            soiloff(ir,ic,idc)   = soiloff(ir,ic,idc)+1000.0*d1_ssub* land_ratio(ir,ic,iland)

!          surface routing: steady constant sheet flow
            soil_con_f = 1.0         ! soil conservation factor
            detension = Sstmax(iland)*soil_con_f
            detension = detension * 2.0
            detension = detension * ((dLAI / LAImax(iland))**0.5)
            detension = amax1(3.0, detension) 

            water_depth = amax1(0.0, (Sst(ir,ic,iland)-detension) )   ! mm
              q_hillslope = 0.0
            if( water_depth .le. 0.01 .or. Drw(isub,iflow).ge. 1.0*Dr(isub,iflow) ) goto 1800

            Sst(ir,ic,iland) = Sst(ir,ic,iland)-water_depth
            water_depth = 0.001 * water_depth          !in meter, surface runoff

            surface_n = surfn(iland) ! * sqrt(water_depth*1000.)
            waterhead = slp(ir,ic) ! + water_depth/length(ir,ic)
            power       = 1.6667

            q_hillslope = dt * sqrt(waterhead)*water_depth**power/surface_n     ! m3/m, one hillslope
            if(q_hillslope .le. 0.1E-20) q_hillslope = 0.0

            if(iland.eq.6)    then! for agricultural fields
              q_hillslope = q_hillslope * amax1(0.3, (1.0-LAI(ir,ic,month,2)/LAImax(iland)))
            endif

            q_hillslope = amin1(q_hillslope, water_depth*length(ir,ic))

            water_depth = water_depth - q_hillslope/length(ir,ic)      
            qin_tmp= qin_tmp + q_hillslope/dt * land_ratio(ir,ic,iland)         !m3/s/m, one hillslope
!            update surface storage
            Sst(ir,ic,iland) = Sst(ir,ic,iland)+1000.0*water_depth
            water_depth = 0.0
            srunoff2(isub,ihc) = srunoff2(isub,ihc) + q_hillslope*land_ratio(ir,ic,iland)/length(ir,ic)*area(ir,ic)/dt         ! m3/s,mqh
1800          continue
            qin_tmp  = qin_tmp *area(ir,ic)/length(ir,ic)! m3/s, one flow-interval
            qin_tmp0 = qin_tmp0 + qin_tmp
            runoffd(ir,ic,idc) = runoffd(ir,ic,idc) + 1000.0*qin_tmp*dt/area(ir,ic)
       end do

!       lateral inflow to the river
        qin(iflow) = qin_tmp0 ! m3/s, total lateral inflow of one flow-interval
!    monthly average soil moisture of the root depth
        do 5444 ig = 1, ngrid(isub,iflow)
            ir    = grid_row(isub,iflow,ig)
            ic    = grid_col(isub,iflow,ig)
            isoil = soil(ir,ic)
          do 4888 iland = 1, nland
            if(land_ratio(ir,ic,iland) .le. 0.0) goto 4888
            tmp = 0.0
            i   = 0
            do j = 1, layer(ir,ic)
              tmp = tmp + D(ir,ic,j)
              i = i+1
              if(tmp .ge. root(iland))     goto 4267  !Ô­ÀŽÊÇ4266

            end do
4267            tmp19921 = 0.0
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
4266        do j = 1, i
              soil_month(ir,ic) = soil_month(ir,ic) + 
     $                            land_ratio(ir,ic,iland)*
     $                            (w(ir,ic,iland,j)-wrsd(ir,ic))/
     $                            ((wsat(ir,ic)-wrsd(ir,ic))*24.0)*
     $                             D(ir,ic,j)/tmp
            end do   
4888      continue
5444    continue  
1999      continue
!     save monthly spatial distribution
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
      end if
19921    if(month.eq.endmonth .and. day.eq.endday .and. hour.eq.24) then
!   (1) sub-basin mean: rain, runoff, Eact
        write(ch4, '(i4.4)') hydroyear     
        call strlen(simul_dir, ia1,ia2)

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




        endif
      end if
        if(isub.eq.start_sub) then
          annual_rain  = 0.0
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
        endif
         
        do iflow = 1, nflow(isub)
            do ig=1,ngrid(isub,iflow)
            ir=grid_row(isub,iflow,ig)
            ic=grid_col(isub,iflow,ig)
            isoil=soil(ir,ic)
            do i = 1, idc         !366
              annual_runoff = annual_runoff + runoffd(ir,ic,i)*
     :                        area(ir,ic)                           ! mm
              annual_Eact   = annual_Eact + eactd(ir,ic,i)*
     :                        area(ir,ic)                                ! mm
             annual_rain   = annual_rain + raind(ir,ic,i)*           !rain_daily
     :                        area(ir,ic)                                 ! mm
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
          annual_rain  = annual_rain / basinarea
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

          write(35,1120) annual_rain, annual_Ssts,annual_Sstr,
     :            annual_SBst,annual_Gst,annual_Eact, annual_runoff
     :            , unbalance,annual_Ecy,annual_Ecp,annual_Ese
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


!       (4a) save the status at end of the simulation
      if(month.eq.endmonth.and.day.eq.endday.and.hour.eq.24) then
        call strlen(simul_dir,ia1,ia2)
        open(9,file = simul_dir(ia1:ia2)//      ! changed by xujijun
     :                subbasin(isub)//'I_soil2')
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

!       (4b) save the status at start of the simulation in every year
      if(month.eq.startmonth.and.day.eq.startday.and.hour.eq.1) then
        write(ch4,'(i4.4)') year
        call strlen(simul_dir,ia1,ia2)
        open(9,file = simul_dir(ia1:ia2)//
     :                subbasin(isub)//'I_soil2'//ch4,status='unknown')
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
