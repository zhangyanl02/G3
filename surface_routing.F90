    subroutine surface_routing(nlayer,isub,year,month,day,hour)
      use global_para_mod,only:i4,r8,dt,idc,ihc
      use water_balance_mod,only:srunoff2,groundoff2,soiloff2,groundoff,soiloff,runoffd
      use hydro_data_mod,only:nflow,ngrid,grid_row,grid_col,slp,length,Ds,Dg,Dr,Drw,kg,GWcs,layer,&
                              w,Sst,Dgl,GWst,Sstmax,surfn,area,k0,D,qin
      use soil_para_mod,only:soil,wsat,wrsd,wfld,alpha,watern,anik
      use land_para_mod,only:land_ratio,LAI,nland,LAImax
      use forcing_mod,only:snow,temp_hour
      
      
      implicit none
      integer(i4),intent(in)::nlayer,year,month,day,hour,isub
      integer(i4)::iflow,ig,ir,ic,isoil,iland
      real(r8)::qin_tmp0
      real(r8)::fice
      real(r8)::d1_slope,d1_length,d1_Ds,d1_Dg,d1_Dr,d1_Drw,d1_kg,d1_GWcs,d1_dt
      real(r8)::d1_wsat,d1_wrsd,d1_wfld,d1_alpha,d1_watern
      integer(i4)::j,d1_nuz
      real(r8)::d1_deltz(nlayer),d1_k0(nlayer),d1_w(nlayer)
      real(r8)::d1_land_ratio,d1_anik
      real(r8)::d1_Sst,d1_Dgl,d1_GWst,d1_snow
      real(r8)::soil_con_f,detension,water_depth,q_hillslope
      real(r8)::surface_n
      real(r8)::waterhead
      real(r8)::power
      real(r8)::dLAI
      real(r8)::d1_qsub,d1_ssub
      real(r8)::temper
      real(r8)::qin_tmp
      srunoff2(isub,ihc)=0
      groundoff2(isub,ihc)=0
      soiloff2(isub,ihc)=0  !mqh
      do iflow=1,nflow(isub)
        qin_tmp0=0.0
        do ig=1,ngrid(isub,iflow)
          ir=grid_row(isub,iflow,ig)
          ic=grid_col(isub,iflow,ig)
          isoil=soil(ir,ic)
          fice = 1.00
          temper=temp_hour(ir,ic,day,hour)
          fice=(temper+7.5)/7.5
          fice=amax1(fice, 0.5)
          fice=amin1(fice , 1.0)
          qin_tmp=0.0
          d1_slope  = slp(ir,ic)         !m/m
          d1_length = length(ir,ic)       ! meter
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
          
          do iland = 1, nland
            if(land_ratio(ir,ic,iland) .gt. 0.0) then
              d1_land_ratio = land_ratio(ir,ic,iland)
              d1_anik = anik(iland)
              dLAI=LAI(ir,ic,iland,day)
              if( d1_anik .lt. 0.0) then
                print *, "wrong anik..."
                stop
              end if
              if( d1_anik .lt. 1.0)  then
                print *, "wrong anik",anik(iland) 
              end if
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
              if( water_depth .gt. 0.01 .and. Drw(isub,iflow).lt. 1.0*Dr(isub,iflow) ) then
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
!               update surface storage
                Sst(ir,ic,iland) = Sst(ir,ic,iland)+1000.0*water_depth
                water_depth = 0.0
                srunoff2(isub,ihc) = srunoff2(isub,ihc) + q_hillslope*land_ratio(ir,ic,iland)/length(ir,ic)*area(ir,ic)/dt         ! m3/s,mqh
              end if
            end if
          end do   ! iland 
            qin_tmp  = qin_tmp *area(ir,ic)/length(ir,ic)! m3/s, one flow-interval
            qin_tmp0 = qin_tmp0 + qin_tmp
            runoffd(ir,ic,idc) = runoffd(ir,ic,idc) + 1000.0*qin_tmp*dt/area(ir,ic)
        end do  !ig
!       lateral inflow to the river
        qin(iflow) = qin_tmp0 ! m3/s, total lateral inflow of one flow-interval
      end do
    end subroutine
