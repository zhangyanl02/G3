      subroutine Initianize()
      implicit none
      
      
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
      end subroutine Initianize
