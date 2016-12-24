













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








!     save monthly spatial distribution
      if(isub.eq.end_sub .and. day.eq.dayinmonth(month) .and. hour.eq.24) then
        write(ch2, '(i2.2)') month
        write(ch4, '(i4.4)') hydroyear


 
        open(33,file = result2_dir(ia1:ia2)//'evap_'//ch4//ch2//'.asc')  ! changed by xujijun
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
