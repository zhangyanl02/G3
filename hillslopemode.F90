!            ***** Hydrological Response on Hillslope  *****          *
      subroutine hillslope_model(isub)
        implicit none
!       include 'dims2.inc'
        


        integer(i4)::isub              ! sub-basin number
        integer(i4)::iflow            ! flow-interval number
        real(r8):: smf                        ! snow melting factor
        integer(i4)::isoil                    ! soil number
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
        integer(i4):: ig,nuz                    !,ii, jj
        integer(i4):: ir, ic            !, iyy, imm, idd

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

        data luolicun2 /1,1,0,0,0,1,1,1,0,0,1,1,1,1,1,1,0,1,0,
     :1,1,1,1,0,0,1,1,0,1,0,1,1,
     : 1,1,1,1,0,1,0,1,1,0,0,1,0,1,1,1,1,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,
     :0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,
     :0,1,0,0,0,0,1,1,1,1,1,1,1,1,0,1,1,0,0,1,0,0,1,0,0,1,0,1,1,1,1,1,1,
     :1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/

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


        data luolicun /1,1,0,0,0,1,1,1,0,0,1,1,1,1,1,1,0,1,0,
     :1,1,1,1,0,0,1,1,0,1,0,1,1,
     : 1,1,1,1,0,1,0,1,1,0,0,1,0,1,1,1,1,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,
     :0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,
     :0,1,0,0,0,0,1,1,1,1,1,1,1,1,0,1,1,0,0,1,0,0,1,0,0,1,0,1,1,1,1,1,1,
     :1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/



      smf  = 0.1     ! snowmelting factor  !mqh,from 0.15 to 0.1
      if(start .eq. 1 .and. isub.eq.start_sub) then
        do ir = 1, nrow
          do 10 ic = 1, ncol
            isoil=soil(ir,ic)
            D(ir,ic,1)   = 0.05  !changed by mqh
            layer(ir,ic) = 1
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
 20         D(ir,ic,j) = D(ir,ic,j) + Ds(ir,ic)-tmp
            if(D(ir,ic,j) .lt. D(ir,ic,j-1)-0.05) then
              D(ir,ic,j-1) = 0.5*(D(ir,ic,j)+D(ir,ic,j-1))
              D(ir,ic,j)   = D(ir,ic,j-1)
            endif

           if(area(ir,ic).eq.-9999.or.isoil.eq.-9999) goto 10

!--------------- compute the k0 for each layer --------------
           tmp = 0.0
             do j = 1, layer(ir,ic)
             tmp = tmp + D(ir,ic,j)
             f   = -alog(ksat2(ir,ic)/ksat1(ir,ic))/Ds(ir,ic)
             k0(ir,ic,j) = ksat1(ir,ic)*exp(-f*tmp)
           end do
 10      continue 
       end do
      end if


!        read soil and groundwater initial conditions
      if(start.eq.1) then
!--------------- read from the files --------------
        if(inicon.eq.1) then
          call strlen(simul_dir,ia1,ia2)
          open(1,file=simul_dir(ia1:ia2)//subbasin(isub)//'I_soil2',status='old')
          do iflow = 1, nflow(isub)
            do ig = 1, ngrid(isub,iflow)
              ir    = grid_row(isub,iflow,ig)
              ic    = grid_col(isub,iflow,ig)
              isoil = soil(ir,ic)
              read(1,*) ia1, ia2, ib1, ib2
              do iland = 1, ib1
                if(land_ratio(ir,ic,iland) .gt. 0.0) read (1,*) (w(ir,ic,iland,j),j=1,ib2),Dgl(ir,ic,iland), Gwst(ir,ic,iland)
                if(Gwst(ir,ic,iland) .ge. Dg(ir,ic)*GWcs(ir,ic)-0.1E-2) Dgl(ir,ic,iland) = Ds(ir,ic)
                if(Gwst(ir,ic,iland) .le. 0.0) Gwst(ir,ic,iland) = 0.0
                if(Dgl(ir,ic,iland)  .ge. Ds(ir,ic)+Dg(ir,ic)) Dgl(ir,ic,iland) = Ds(ir,ic)+Dg(ir,ic)
              end do
            end do
          end do
          close(1)
!--------------- specify by the following way --------------
        else
          if(isub.eq.start_sub) then
            do ir=1,nrow
              do ic=1,ncol
                isoil=soil(ir,ic)
                if(isoil .eq. -9999) goto 28
                do iland=1,nland
                  Dgl(ir,ic,iland)  = Ds(ir,ic)
                  GWst(ir,ic,iland) = (Ds(ir,ic)+Dg(ir,ic)-Dgl(ir,ic,iland))*Gwcs(ir,ic)
                      tmp = 0.0
                  D0  = Ds(ir,ic)*wrsd(ir,ic)/(wsat(ir,ic)-wrsd(ir,ic))
                  do j= 1, layer(ir,ic)
                    tmp = tmp+D(ir,ic,j)
                    w(ir,ic,iland,j) = wsat(ir,ic)*(D0+tmp)/(D0+Ds(ir,ic))
                    if(isoil.eq.-9999) w(ir,ic,iland,j) = wfld(ir,ic)
                  end do
                end do
28            continue
              end do
            end do
          end if
        end if
      end if
!            sub-basin area Ãæ»ý
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

!                     initialize status variables
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

      if(isub.eq.start_sub .and. (start.eq.1 .or.(idc.eq.1 .and.ihc.eq.1))) then
        do ir = 1, nrow
          do 50 ic = 1, ncol
            do i=1,366
              raind(ir,ic,i)     = 0.0
              eactd(ir,ic,i)     = 0.0
              ecy(ir,ic,i)       = 0.0
              ecp(ir,ic,i)       = 0.0
              ese(ir,ic,i)       = 0.0
              runoffd(ir,ic,i)   = 0.0
              srunoff(ir,ic,i)   = 0.0
              groundoff(ir,ic,i) = 0.0
              soiloff(ir,ic,i)   = 0.0
            end do
 50       continue
        end do
        annual_Cst0  = 0.0
        annual_Ssts0 = 0.0
        annual_Sstr0 = 0.0
        annual_SBst0 = 0.0
        annual_Gst0  = 0.0
      end if

      if( start .eq. 1 .or. (idc.eq.1 .and. ihc.eq.1)) then 
        do iflow = 1, nflow(isub)
!          annual_re0 = annual_re0 + re_capacity(isub,iflow)
            do ig = 1, ngrid(isub,iflow)
            ir = grid_row(isub,iflow,ig)
            ic = grid_col(isub,iflow,ig)
!            isoil = soil(ir,ic)
            do  iland = 1, nland
              if(land_ratio(ir,ic,iland) .gt. 0.0) then
              annual_Cst0 = annual_Cst0 + Cst(ir,ic,iland)*0.001*land_ratio(ir,ic,iland)*area(ir,ic) !in m3
              annual_Sstr0= annual_Sstr0 + Sst(ir,ic,iland)* 0.001*land_ratio(ir,ic,iland)*area(ir,ic) !in m3
              annual_Ssts0= annual_Ssts0 + snow(ir,ic,iland)*0.001*land_ratio(ir,ic,iland)*area(ir,ic) !in m3
              annual_Gst0 = annual_Gst0 + GWst(ir,ic,iland)*area(ir,ic)*land_ratio(ir,ic,iland)       !in m3
              do j = 1 ,layer(ir,ic)
                annual_SBst0 = annual_SBst0 + w(ir,ic,iland,j) * D(ir,ic,j)*land_ratio(ir,ic,iland)*area(ir,ic)      !in m3
              end do
              end if
            end do 
          end do
        end do
      end if

     
!                     initialize irrigation variables
      if(start .eq. 1) return
!     start simulation
!          interception, evapotranspiration   
      do 999 iflow=1,nflow(isub)
        do ig=1,ngrid(isub,iflow)
          ir=grid_row(isub,iflow,ig)
          ic=grid_col(isub,iflow,ig)
          isoil=soil(ir,ic)
! variable using for water balance calculation
          do 888 iland = 1, nland
            call rain_model(isub, ir, ic, hour, month,rain_daily(ir,ic,day),tmin_daily(ir,ic,day),&
                tmax_daily(ir,ic,day), evap_daily(ir,ic,day,iland),prec, temper, Ep, iland)
            prec = pre_hour(ir,ic,day,hour)
            if(iland.eq.1) then
               raind(ir,ic,idc) = raind(ir,ic,idc) + prec
            end if
            pnet           = prec
            Eact           = 0.0
            if(land_ratio(ir,ic,iland).le.0.0) goto 888
            if(iflai .eq. 0 ) then
!          NDVI -> LAI  
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
                do i_ndvi=1, 19
                  if(NDVI(ir,ic,j,kn) .le. T_NDVI(1)) ii_ndvi = 1
                  if(NDVI(ir,ic,j,kn).gt.T_NDVI(i_ndvi) .and.NDVI(ir,ic,j,kn).le.T_NDVI(i_ndvi+1)) then
                    ii_ndvi = i_ndvi
                  end if    
                  if(NDVI(ir,ic,j,kn) .gt. T_NDVI(20)) ii_ndvi = 20
                end do
                c4 = 0.4  !changed from 0.4 to 0.3
                if(iland.eq.1 .or. iland.eq.10) then  !waterbody or frost snow
                  LAI(ir,ic,j,kn) = 0.0
                end if 

                if(iland.eq.2) then  !urban-area
                   if(ii_ndvi.eq.0 .or. ii_ndvi.eq.20) then
                     LAI(ir,ic,j,kn)=Shurb_LAI(ii_ndvi)           
                   else
                     LAI(ir,ic,j,kn) = Shurb_LAI(ii_ndvi) + (NDVI(ir,ic,j,kn)-T_NDVI(ii_ndvi))*&
                       (Shurb_LAI(ii_ndvi+1)-Shurb_LAI(ii_ndvi))/(T_NDVI(ii_ndvi+1)-T_NDVI(ii_ndvi))           
                   end if
                   c4 = 0.06
                end if 

                if(iland.eq.3) then  !baresoil
                  LAI(ir,ic,j,kn) = 0.2*NDVI(ir,ic,j,kn)+0.2
                  c4=0.4 ! changed by gaobing
                end if 

                if(iland.eq.4) then  !forest
                  if(ii_ndvi.eq.0 .or. ii_ndvi.eq.20) then
                    LAI(ir,ic,j,kn) = Forest_LAI(ii_ndvi)           
                  else
                    LAI(ir,ic,j,kn)=Forest_LAI(ii_ndvi)+(NDVI(ir,ic,j,kn)-T_NDVI(ii_ndvi))*(Forest_LAI(ii_ndvi+1)-&
                        Forest_LAI(ii_ndvi))/(T_NDVI(ii_ndvi+1)-T_NDVI(ii_ndvi))           
                  end if
                  if(iland.eq.4) c4 = 0.2
                end if 

                if(iland.ge.5) then  ! sparse vegetion
                  LAI(ir,ic,j,kn)=1.71*NDVI(ir,ic,j,kn)+0.48
                  c4 = 0.3
                end if 

                if(iland.eq.6) then  !upland
                  if(ii_ndvi.eq.0 .or. ii_ndvi.eq.20) then
                    LAI(ir,ic,j,kn)=N_IRR_LAI(ii_ndvi)           
                  else
                    LAI(ir,ic,j,kn)=N_IRR_LAI(ii_ndvi)+(NDVI(ir,ic,j,kn)-T_NDVI(ii_ndvi))*(N_IRR_LAI(ii_ndvi+1)-&
                      N_IRR_LAI(ii_ndvi))/(T_NDVI(ii_ndvi+1)-T_NDVI(ii_ndvi))           
                  end if
                end if 

                if(iland.eq.7) then  !grassland
                  LAI(ir,ic,j,kn)=1.71*NDVI(ir,ic,j,kn)+0.48
                  if(iland.eq.7) c4 = 0.3
                end if 

               if(iland.eq.8) then  !shrub
                 if(ii_ndvi.eq.0 .or. ii_ndvi.eq.20) then
                   LAI(ir,ic,j,kn)=Shurb_LAI(ii_ndvi)           
                 else
                   LAI(ir,ic,j,kn)=Shurb_LAI(ii_ndvi)+(NDVI(ir,ic,j,kn)-T_NDVI(ii_ndvi))*(Shurb_LAI(ii_ndvi+1)-&
                      Shurb_LAI(ii_ndvi))/(T_NDVI(ii_ndvi+1)-T_NDVI(ii_ndvi))           
                 end if
                 c4 = 0.08
               end if 

               if(iland.eq.9) then  !wetland
                 LAI(ir,ic,j,kn) = 1.71*NDVI(ir,ic,j,kn)+0.48
               end if 

               if(LAI(ir,ic,j,kn) .lt. 0.0) LAI(ir,ic,j,kn) = 0.0
               if(LAI(ir,ic,j,kn) .gt. LAImax(iland)) then
                 LAImax(iland) = LAI(ir,ic,j,kn)
               end if
             enddo ! kn
           end do

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

!          interception  
             Cstmax(ir,ic,iland) = 0.10*dLAI
             if(iland.eq.4 .or. iland.eq.8) then !forest or shurb
               Cstmax(ir,ic,iland) = 1.00*dLAI 
             end if
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
       
             if(temper .le. 1.0) then
               snow(ir,ic,iland) = snow(ir,ic,iland) + pnet
               pnet=0.0 
             end if
             snowmelt=0.0

             if(temper .gt. 0.0 .and. snow(ir,ic,iland).gt.0.0) then
               if(month.ge.3.and.month.le.10) smf= 0.15
               if(month.ge.5.and.month.le.8)  smf = 0.15
               snowmelt = (smf+pnet/20.0) * (temper-1.5)
               snowmelt = amin1(snowmelt,snow(ir,ic,iland))
               snow(ir,ic,iland)=snow(ir,ic,iland)-snowmelt             
               pnet=pnet+snowmelt 
             end if

             Sst(ir,ic,iland) = Sst(ir,ic,iland) + pnet
!            (1) evaporation from canopy storage  
250         continue
            EfromCanopy = 0.0
            EfromCrop   = 0.0
            EfromSurface= 0.0
            c2=0.05+0.1*dLAI/LAImax(iland)
            c1=0.31*dLAI
            c3 =  0.23+0.1*(dLAI/LAImax(iland))**0.5            ! 
            c5 =      0.23
            if (luolicun2(isub).eq.1) then
              c1=c1*0.5
              c2=c2*0.5
              c3=c3*0.5
              c5=c5*0.4
            endif


            Etr = Ep*amin1(c2+c1,1.0)
            Es  = (Ep*(c2+(1-c2)*(1-amin1(c2+c1,1.0))) -Etr*(1-amin1(c2+c1,1.0)))*&
            (1-exp(-c4*dLAI))+Ep*exp(-c4*dLAI)*(c5+(1.0-c5)*dLAI/LAImax(iland))**(1.0+c3)
            if(iland.eq.3) then  !baresoil
              Es  = Ep*exp(-c4*dLAI)
            end if
                 
            if(iland.eq.1) then  !waterbody
              Es  = Ep           !*(1-Kcanopy(iland))
            end if   
            Etr = Etr*(1-exp(-c4*dLAI))      !Kcanopy(iland)

            if(Etr .lt. 0.1E-9) goto 270

            EfromCanopy = amin1(Etr, Cst(ir,ic,iland)) 
            Cst(ir,ic,iland) = Cst(ir,ic,iland) - EfromCanopy
            if(EfromCanopy.lt.0.0) print *,'wrong in E from canopy',ihc,iflow,EfromCanopy

!          (2) transpiration from vegetation
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
              y = (1.0-(1.0-para_r)*float(j-1)/float(i))*para_a
              call ETsoil(w(ir,ic,iland,j), wfld(ir,ic), wrsd(ir,ic), kmoist)
              tmp = y * Etr * kmoist
              tmp=amin1(tmp,temp) 
              w(ir,ic,iland,j) = w(ir,ic,iland,j) -  tmp/(1000.0*D(ir,ic,j))
              EfromCrop = EfromCrop + tmp
            end do

!          (3) evaporation from soil surface

 270        if(Es .lt. 0.1E-9) goto 290
            if(Sst(ir,ic,iland) .le. 0.0) goto 275
            EfromSurface     = amin1(Es, Sst(ir,ic,iland))
275         if(Sst(ir,ic,iland) .lt. 0.0) Sst(ir,ic,iland) = 0.0
            tmpEp = Es-EfromSurface 
            if(tmpEp .le. 0.1E-9) goto 290
            call ETsoil( w(ir,ic,iland,1), wfld(ir,ic), wrsd(ir,ic), kmoist )
            tmp = tmpEp * kmoist
            if(tmp .lt. 0.0) print *, 'wrong in EfromSoil', tmp
            temp = amax1(0.0, 1000.0*(w(ir,ic,iland,1)-wrsd(ir,ic)-1.0E-6)*D(ir,ic,1))
            tmp = amin1(tmp,temp)  
            w(ir,ic,iland,1) = w(ir,ic,iland,1) - tmp/(1000.0*D(ir,ic,1))
            EfromSurface = EfromSurface + tmp
290        continue

            eactd(ir,ic,idc) = eactd(ir,ic,idc) + (EfromCanopy+EfromCrop+EfromSurface) *land_ratio(ir,ic,iland)            ! mm
            ecy(ir,ic,idc)   =  ecy(ir,ic,idc) + EfromCanopy *land_ratio(ir,ic,iland) 
            ecp(ir,ic,idc)   =  ecp(ir,ic,idc) + EfromCrop *land_ratio(ir,ic,iland) 
            ese(ir,ic,idc)   =  ese(ir,ic,idc) + EfromSurface*land_ratio(ir,ic,iland)                       ! added by gaobing
            eactd(ir,ic,idc) = eactd(ir,ic,idc) + Eact * land_ratio(ir,ic,iland)        ! mm
 888      continue
        end do
 999  end do

!          UZ, GW and surface routing
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
            groundoff(ir,ic,idc) = groundoff(ir,ic,idc) + 1000.0* 
     :                             dt*d1_qsub*land_ratio(ir,ic,iland)/
     :                             length(ir,ic)         ! mm
            soiloff(ir,ic,idc)   = soiloff(ir,ic,idc)+1000.0*d1_ssub* 
     :                                 land_ratio(ir,ic,iland)

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
