    subroutine Evaportranspiration
      use global,only:r8,i4
      implicit none
      real(r8):: smf                        ! snow melting factor
 
      smf  = 0.1     ! snowmelting factor  !mqh,from 0.15 to 0.1
      do iflow=1,nflow(isub)
        do ig=1,ngrid(isub,iflow)
           ir=grid_row(isub,iflow,ig)
           ic=grid_col(isub,iflow,ig)
           isoil=soil(ir,ic)
           prec  = pre_hour (ir,ic,day,hour)
           temper= temp_hour(ir,ic,day,hour)
           Ep    = Ep_hour  (ir,ic,day,hour)
          do iland = 1, nland
            if(iland.eq.1) then
              raind(ir,ic,idc) = raind(ir,ic,idc) + prec
            end if
            pnet           = prec
            Eact           = 0.0
            if(land_ratio(ir,ic,iland).gt.0.0) then
              dLAI=LAI(ir,ic,iland,day)
!interception
              Cstmax = 0.10*dLAI
              if(iland.eq.4 .or. iland.eq.8) then !forest or shurb
                Cstmax = 1.00*dLAI 
              end if            
              if(pnet.gt. 0.0) then
                DeficitCst = Cstmax- Cst(ir,ic,iland)
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
              end if
            
!(1) evaporation from canopy storage
              EfromCanopy = 0.0
              EfromCrop   = 0.0
              EfromSurface= 0.0
              c2=0.05+0.1*dLAI/LAImax(iland)
              c1=0.31*dLAI
              c3=0.23+0.1*(dLAI/LAImax(iland))**0.5 
              c5=0.23

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

! (2) transpiration from vegetation
              if(Etr .ge. 0.1E-9) then
                EfromCanopy = amin1(Etr, Cst(ir,ic,iland)) 
                Cst(ir,ic,iland) = Cst(ir,ic,iland) - EfromCanopy
                if(EfromCanopy.lt.0.0) print *,'wrong in E from canopy',ihc,iflow,EfromCanopy
                if(root(iland).ne.0.0) then
                  i      = 0          ! layer of root           
                  para_r = 0.1        ! root density ratio of deepest layer (<=1.0)
                  tmp    = 0.0
                  do j = 1, layer(ir,ic)
                    tmp = tmp + D(ir,ic,j)
                    i   = i + 1
                    if(tmp .ge. root(iland)) then
                      exit
                    end if
                  end do
                  if(layer(ir,ic) .le. i) i = layer(ir,ic)
                  para_a = 1.0/(float(i) - 0.5*(1.0-para_r)*float(i-1))
                  
                  do j=1,i
                    y = (1.0-(1.0-para_r)*float(j-1)/float(i))*para_a
                    call ETsoil(w(ir,ic,iland,j), wfld(ir,ic), wrsd(ir,ic), kmoist)
                    tmp = y * Etr * kmoist
                    tmp=amin1(tmp,temp) 
                    w(ir,ic,iland,j) = w(ir,ic,iland,j) -  tmp/(1000.0*D(ir,ic,j))
                    EfromCrop = EfromCrop + tmp
                 end do
                endif
              endif
!(3) evaporation from soil surface
              if(Es .ge. 0.1E-9) then
                if(Sst(ir,ic,iland) .gt. 0.0) then
                  EfromSurface     = amin1(Es, Sst(ir,ic,iland))
                end if
                if(Sst(ir,ic,iland) .lt. 0.0) Sst(ir,ic,iland) = 0.0
                tmpEp = Es-EfromSurface 
                if(tmpEp .gt. 0.1E-9) then
                  call ETsoil( w(ir,ic,iland,1), wfld(ir,ic), wrsd(ir,ic), kmoist )
                  tmp = tmpEp * kmoist
                  if(tmp .lt. 0.0) print *, 'wrong in EfromSoil', tmp
                  temp = amax1(0.0, 1000.0*(w(ir,ic,iland,1)-wrsd(ir,ic)-1.0E-6)*D(ir,ic,1))
                  tmp = amin1(tmp,temp)  
                  w(ir,ic,iland,1) = w(ir,ic,iland,1) - tmp/(1000.0*D(ir,ic,1))
                  EfromSurface = EfromSurface + tmp
                end if
              end if
              eactd(ir,ic,idc) = eactd(ir,ic,idc) + (EfromCanopy+EfromCrop+EfromSurface) *land_ratio(ir,ic,iland)            ! mm
              ecy(ir,ic,idc)   = ecy(ir,ic,idc) + EfromCanopy *land_ratio(ir,ic,iland) 
              ecp(ir,ic,idc)   = ecp(ir,ic,idc) + EfromCrop *land_ratio(ir,ic,iland) 
              ese(ir,ic,idc)   = ese(ir,ic,idc) + EfromSurface*land_ratio(ir,ic,iland)                       ! added by gaobing
              eactd(ir,ic,idc) = eactd(ir,ic,idc) + Eact * land_ratio(ir,ic,iland)        ! mm
            end if
          enddo
        end do
      end do
    end subroutine Evaportranspiration
