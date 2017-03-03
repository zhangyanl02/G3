    module hydro_mod
      use global_para_mod,only:i4,r8
      public
      save
    contains
    
    
    subroutine read_hydro_para(para_dir,nhrpdt)
      use shr_kind_mod, only:r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8,r4=> shr_kind_r4      
      implicit none
      character*200,intent(in)::para_dir
      integer(i4),intent(in)::nhrpdt
      dthydro=nhrpdt*3600.0
      !call initial_subcatchment(para_dir)
      call read_river_para(para_dir)
      return 
    end subroutine read_hydro_para
    
    
    
    
!##########################
!##########################
    subroutine initial_subcatchment(para_dir)
      use shr_kind_mod, only:r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8,r4=> shr_kind_r4
      implicit none
      character*200,intent(in)::para_dir
      integer::i,j,k,fileunit
      character(200)::atmp
      integer(i4)::l1,l2,subcount,iostatus
      character*4::ch4 
      character*200::subbasinfile
      call strlen(para_dir,l1,l2)
      subbasinfile=trim(para_dir(l1:l2))//'subbasin.dat'
      call getunit(fileunit)
      open(fileunit,file = subbasinfile, status='old')
      subcount=0
      iostatus=1
      do while(iostatus.ge.0)	  
        read(fileunit,*,iostat=iostatus) atmp
        subcount=subcount+1
      end do
      close(fileunit)
      subcount=subcount-1
      nsub=subcount

      call strlen(para_dir, l1, l2)
      open(fileunit,file=trim(para_dir(l1:l2))//'riverpara/'//'maxnp',status='old')
      read(fileunit,*) nflowmax,ngridmax
      close(fileunit)

      allocate(subbasin(nsub))
      allocate(psubbasin(nsub))
      allocate(pbasinup(nsub,8))
      allocate(nbasinup(nsub,8))
      !allocate(nextbasinrow(nsub),nextbasincol(nsub),nextbasin(nsub))
      allocate(nflow(nsub))
      allocate(ngrid(nsub,nflowmax))
      allocate(dx(nsub,nflowmax))
      allocate(s0(nsub,nflowmax))
      allocate(b(nsub,nflowmax))
      allocate(roughness(nsub,nflowmax))
      allocate(Dr(nsub,nflowmax))
      allocate(grid_row(nsub,nflowmax,ngridmax))
      allocate(grid_col(nsub,nflowmax,ngridmax))
      allocate(hrr(nsub,nflowmax))
      allocate(q1(nsub,nflowmax))
      allocate(q2(nsub,nflowmax))	
      allocate(qlin1(nsub,nflowmax))
      allocate(qlin2(nsub,nflowmax))
      allocate(qin(nsub))
      allocate(area_sub(nsub))
  
      open(fileunit,file = subbasinfile, status='old')
      do i = 1,nsub
        read(fileunit,*)nsub,psubbasin(i),(pbasinup(i,j),j=1,8)!,nextbasinrow(i),nextbasincol(i),nextbasin(i)
      enddo
      close(fileunit)


      do i =1,nsub
        write(ch4,'(i4.4)') psubbasin(i)
        subbasin(i)='ws'//ch4
      enddo
      
      nbasinup=0     
      do i = 1,nsub
        do k =1,8
          if(pbasinup(i,k)>0) then
            do j =1,nsub
              if(psubbasin(j)== pbasinup(i,k)) then
                nbasinup(i,k)=j 
              endif
            enddo
          endif
        enddo
      enddo  
      call retunit(fileunit)
    end subroutine




    subroutine read_river_para(para_dir)
      use shr_kind_mod, only:r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8,r4=> shr_kind_r4
      implicit none
      character*200,intent(in)::para_dir
      integer(i4)::isub,l1,l2,ll1,ll2,iflow,j,fileunit
      character*200::infile
      call getunit(fileunit)
      do isub = 1, nsub
          infile = subbasin(isub)//'_river'
        !print *,isub,subbasin(isub),"  ",infile
          call strlen(infile, l1, l2)
          call strlen(para_dir, ll1, ll2)
          open(fileunit,file=trim(para_dir(ll1:ll2))//'riverpara/'//infile(l1:l2),status='old')
          read(fileunit,*) nflow(isub)
          do iflow = 1, nflow(isub)
            read(fileunit,*) ngrid(isub,iflow),dx(isub,iflow),s0(isub,iflow),b(isub,iflow),roughness(isub,iflow),Dr(isub,iflow)
            if(dx(isub,iflow) .lt. 0.1) dx(isub,iflow) = 2000.0
            if(s0(isub,iflow) .le. 0 ) then   
              s0(isub,iflow)=0.00001
              print *, 'wrong in s0',isub,iflow,s0(isub,iflow)  
            endif  
            if(s0(isub,iflow).eq.-9999.0) then 
              s0(isub,iflow)=0.00001 
              print *, 'wrong in s0',isub,iflow,s0(isub,iflow)  
            endif
            if(s0(isub,iflow).le.0.1E-5) s0(isub,iflow) = 0.1E-5

            if(ngrid(isub,iflow) .gt. ngridmax) then
              print *,'more than np grids:',ngrid(isub,iflow),isub,iflow,ngridmax
            endif
            read(fileunit,*) (grid_row(isub,iflow,j),grid_col(isub,iflow,j),j = 1, ngrid(isub,iflow))
          end do
          close(fileunit)
          call retunit(fileunit)
      end do
    end subroutine



!      calculate soil hydraulic conductivity by Van Genuchten's equation                                 
      function conductivity_V(k0,wsat,wrsd,n,w)
        implicit none
        real(r8)::conductivity_V
        real(r8),intent(in)::k0          ! Saturated hydraulic conductivity
        real(r8),intent(in)::wsat        ! Saturated water content
        real(r8),intent(in)::wrsd        ! residual water content
        real(r8),intent(in)::n           ! VG n
        real(r8),intent(in)::w           ! water content
        
        real(r8)::m           ! VG m
        real(r8)::tmpw,se        
        m=1.0-1.0/n
        tmpw = w
        if(tmpw.gt.wsat) tmpw=wsat
        if(tmpw.le.wrsd) then
          tmpw=wrsd+1.0E-10
        elseif(tmpw.lt.wrsd+0.0001)then
          tmpw=wrsd+0.0001
        end if

        se=(tmpw-wrsd)/(wsat-wrsd)
        conductivity_V=k0*sqrt(se)*(1.0-(1.0-se**(1.0/m))**m)**2
        if(conductivity_V.lt.0.0 .or. conductivity_V.gt.k0) then
          print *,'wrong in calculating conductivity',w,conductivity_V
        end if
       return
     end function conductivity_V


    SUBROUTINE runoff(isub,month,snow,temper,land_ratio,slope,length,Ds,Dg,Dr,anik,wsat,wrsd,wfld,&
                      alpha,watern,nuz,deltz,k0,w,kg,GWcs,GWst,Dgl,Drw,Sst,dt,qsub,ssub)
      implicit none
      integer(i4),intent(in)::isub        ! the indicator for a subcatchment
      integer(i4),intent(in)::month       ! month in a year
      integer(i4),intent(in)::nuz         ! number of UZ layers
      real(r8),intent(in)::snow           ! snow depth (mm-water)
      real(r8),intent(in)::temper         ! air temperature (degree)
      real(r8),intent(in)::land_ratio     ! area ratio of the land-soil type on a hillslope
      real(r8),intent(in)::slope          ! slope gradient of hillslope (ND)
      real(r8),intent(in)::length         ! length of a hillslope (m)
      real(r8),intent(in)::Ds             ! depth of topsoil (m)
      real(r8),intent(in)::Dg             ! depth of unconfined acquifer (m)
      real(r8),intent(in)::Dr             ! depth of river (m)
      real(r8),intent(in)::anik           ! soil anitropic ratio (>=1.0)
      real(r8),intent(in)::wsat           ! saturated soil moisture, ?????????????????what is the unit?
      real(r8),intent(in)::wrsd           ! residual soil moisture
      real(r8),intent(in)::wfld           ! soil moisture at field capacity
      real(r8),intent(in)::alpha          ! soil-water parameter in VG
      real(r8),intent(in)::watern         ! soil-water parameter in VG
      real(r8),intent(in)::deltz(nuz)     ! depth of each layer (m)
      real(r8),intent(in)::k0(nuz)        ! saturated hydraulic conductivity of each layer (m/s)
      real(r8),intent(in)::dt             ! time step (second)
      real(r8),intent(in)::kg             ! GW hydraulic conductivity (m/s)
      real(r8),intent(in)::GWcs           ! GW storage coefficient
      real(r8),intent(inout)::GWst        ! Groundwater storage (m) 
      real(r8),intent(inout)::Dgl         ! depth to GW level (m)      
      real(r8),intent(inout)::Drw            ! water depth in river (m)
      real(r8),intent(inout)::Sst         ! surface storage (m)
      real(r8),intent(inout)::qsub        ! runoff from subsurface (m3/s/m)
      real(r8),intent(inout)::ssub        ! runoff from soil layer             !added by gaobing
      real(r8),intent(inout)::w(nuz)      ! soil moisture of each layer 
      
      
      real(r8)::suct(nuz)                 ! suction of each layer (meter H2O)
      real(r8)::ksoil(nuz)                ! hydraulic conductivity of each layer (m/s)
      real(r8)::avk(nuz)                  ! average hydraulic conductivity
      real(r8)::avz(nuz)                  ! average depth 
      real(r8)::dpdw(nuz)                 ! 
      real(r8)::avkg                      ! average GW hydraulic conductivity (m/s)
      real(r8)::ss_f                      ! slope shape (concave or convex) factor
      real(r8)::tmp_D
      real(r8)::tmpqsub
      real(r8)::ps1,ps2
      real(r8)::rech                      !recharge rate to groundwater
      real(r8)::conductivity_V
      real(r8)::excess                    ! Saturation excess
      real(r8)::deficit
      real(r8)::qmin
      real(r8)::qmax
      real(r8)::cmtmp,tmp,tmp2
      integer(i4)::mark
      real(r8)::temA(nuz),temB(nuz),temC(nuz),temD(nuz),temq(nuz)
      real(r8)::inf                       !infiltration
      real(r8)::f                         !mqh
      integer(i4)::i,j,isat,GWkey



      temq(:)=0.0
!------------------------ concave hillslope -------------------------
      ss_f=0.02            !linear or convex slope
!     calculation of recharge rate to groundwater: rech                                                    
      inf=0.0
      rech=0.0
      tmp=0.0
      if( Dgl .ge. Ds) then
        call SuctionFromMoisture_V(w(nuz),wsat,wrsd,watern,alpha,ps1)    ! m
        tmp=wfld*(1.0-(Dgl-Ds)/Dg)
        tmp=amax1(tmp,wrsd)
        call SuctionFromMoisture_V(tmp,wsat,wrsd,watern,alpha,ps2)       ! m
        avz(nuz)=0.5*(deltz(nuz)+(Dgl-Ds))                               ! m
        avk(nuz)=conductivity_V(k0(nuz),wsat,wrsd,watern,w(nuz))
        tmp=0.5*(deltz(nuz)+(Dgl-Ds))
        rech=-avk(nuz)*(ps2-ps1)/tmp + avk(nuz)                          ! m/s
        rech=rech*dt

        if(rech .ge. 0.0) then
          tmp = (w(nuz)-wfld)*deltz(nuz)
          tmp = amax1(tmp, 0.0)
          tmp2= Dg*GWcs - GWst
          rech= amin1(rech, tmp, tmp2)
        else
          rech   = amax1(rech, -1.0*GWst)
          tmp    = amax1( 0.0, (wsat-w(nuz))*deltz(nuz) )
          rech   = amax1(rech, -1.0*tmp)
          !if(rech.gt.1.0E-10) print *,"rech>0" ,rech,tmp,GWst, GWcs
        endif
        !total_rech=total_rech+rech         !mqh
        if(abs(rech) .lt. 0.1E-20) then   ! goto 30
          rech = 0.0
        else        
          w(nuz) = w(nuz) - rech/deltz(nuz)
          GWst   = GWst + rech
          Dgl    = Dgl  - rech / GWcs
          if(Dgl.gt.(Ds+Dg+1.0E-8)) then
            print *,"Dgl>Ds+Dg",Dgl,Ds,Dg,rech,ps2,ps1,avk(nuz)
          end if
          if(Gwst.lt.0) then
            print *,"Gwst<0",Gwst
          end if
          if(Dgl.gt.(Ds+Dg+1.0E-8) .or. Gwst.lt.0 ) then
            w(nuz) = w(nuz) + rech/deltz(nuz)
            GWst   = GWst - rech
            Dgl    = Dgl  + rech / GWcs
          end if
        end if
!30        rech = 0.0           the orignial source code.was wrong. This need to be checked.
      endif


!     calculate infiltration
      inf    = 0.0
      excess = 0.0
      if(Sst.gt.0.0) then
!     ununderstand
!     if(isub.ge.44 .and. isub.le.94) then  !infiltration excess
!       tmp    = 0.0
!       tmp2   = (dt*k0(1)/Sst)  * sqrt(w(1)/wsat)
!       tmp    = Sst*exp(-tmp2)
!       excess = amin1(tmp, Sst)
!       Sst    = Sst - excess
!       if(excess .gt. 0.1) print *, excess, month, day
!     endif
        inf  = amax1(k0(1)*dt, 0.0)
        inf  = amin1(inf,Sst)
        Sst  = Sst - inf
        tmp  = w(1) + inf/deltz(1)         ! the unit of w is volummetric water content
        w(1) = amin1(wsat, tmp)
        Sst  = Sst + (tmp-w(1))*deltz(1)
        Sst  = Sst + excess
        excess = 0.0
        if(Sst .lt. 0.1E-9) Sst=0.0
        inf = 0.0  !ÎÊÌâ1
      endif


!     calculate ksoil
      do i=1,nuz
        ksoil(i)=conductivity_V(k0(i),wsat,wrsd,watern,w(i))
        call SuctionFromMoisture_V(w(i),wsat,wrsd,watern,alpha,suct(i))
      end do

!    calculation of inter-layer exchanges of water due to gravitation and hydraulic gradient
      do i=1,nuz-1
        avz(i)=0.5*(deltz(i)+deltz(i+1))
        tmp   = (w(i)*w(i+1)*(deltz(i)+deltz(i+1)))/(w(i)*deltz(i)+w(i+1)*deltz(i+1))
        tmp2  = (k0(i)*deltz(i)+k0(i+1)*deltz(i+1))/(deltz(i)+deltz(i+1))
        avk(i)= conductivity_V(tmp2,wsat,wrsd,watern,tmp)
!       if(w(i) .lt. w(i+1)) avk(i) = 0.1*avk(i)   ! used in SiB2

        if(abs(w(i+1)-w(i)) .lt. 0.1E-5) then
          dpdw(i) = 0.0
        else
          dpdw(i)=(suct(i+1)-suct(i))/(w(i+1)-w(i))
        endif
      end do

!     backward implicit calculation of flows between soil layers.      !¿Ž²»Ì«¶®         
      temA(1)=0.0
      temB(1)=1.0+avk(1)*dpdw(1)/avz(1)*(dt/deltz(1)+dt/deltz(2))
      temC(1)=-avk(1)*dpdw(1)/avz(1) * dt/deltz(2)
      temD(1)=avk(1)*dpdw(1)/avz(1) * (w(1)-w(2)+ inf*dt/deltz(1))+avk(1)
      do i=2,nuz-2, 1
        temA(i)=-avk(i)*dpdw(i)/avz(i) * dt/deltz(i)
        temB(i)=1.0+avk(i)*dpdw(i)/avz(i) * (dt/deltz(i)+dt/deltz(i+1))
        temC(i)=-avk(i)*dpdw(i)/avz(i) * dt/deltz(i+1)
        temD(i)=avk(i)*dpdw(i)/avz(i) * (w(i)-w(i+1)) + avk(i)
      end do
      temA(nuz-1)=-avk(nuz-1)*dpdw(nuz-1)/avz(nuz-1) * dt/deltz(nuz-1)
      temB(nuz-1)=1.0+avk(nuz-1)*dpdw(nuz-1)/avz(nuz-1)*(dt/deltz(nuz-1)+dt/deltz(nuz))
      temC(nuz-1)=0.0
      temD(nuz-1)=avk(nuz-1)*dpdw(nuz-1)/avz(nuz-1) * (w(nuz-1)-w(nuz)+ rech*dt/deltz(nuz)) + avk(nuz-1)
      if(nuz.eq.3) then
         temq(1)=(temD(2)-temD(1)*temB(2)/temC(1))/(temA(2)-temB(1)*temB(2)/temC(1)) 
         temq(2)=(temD(2)-temD(1)*temA(2)/temB(1))/(temB(2)-temA(2)*temC(1)/temB(1))
      else
        call TRDIG (nuz-1, temA, temB, temC, temD, temq, MARK)
        if(MARK .ne. 1) then
           print *, 'wrong in solving  linear equations...'
           stop
        end if
      endif


!     update moisture of each soil moisture layer due to layer interflow.         
      do i = 1, nuz-1 
        qmax = amax1((w(i)-wrsd-1.0E-6) * deltz(i)/dt,0.0)   
        qmin = amin1(-(w(i+1)-wrsd-1.0E-6) * deltz(i+1)/dt,0.0)
        if(temq(i).lt.0.0 .and. temq(i+1).gt.0.0) then
          qmin = 0.4*qmin
        endif
        if(i.ge.2) then
          if(temq(i-1).lt.0.0 .and. temq(i).gt.0.0) then
            qmax = 0.5*qmax
          end if
        endif
        temq(i) = amin1(temq(i), qmax)  
        temq(i) = amax1(temq(i), qmin)
        w(i)   =  w(i)   - temq(i) * dt / deltz(i) 
        w(i+1) =  w(i+1) + temq(i) * dt / deltz(i+1)     
      enddo
      w(nuz) = w(nuz) - rech*dt/deltz(nuz)
!     saturation excess
      GWkey  = 0
      excess = 0.0
      excess = (w(nuz) - wsat) * deltz(nuz)
      excess = amax1(excess, 0.0)
      w(nuz) = w(nuz) - excess / deltz(nuz)
      if(abs(Dgl-Ds).lt.0.1E-5 .and. abs(w(nuz)-wsat).lt.0.1E-5) then 
        Dgl   = Ds - deltz(nuz)
        GWkey = 1
      else 
        GWkey = 0
      endif
      do i=nuz-1, 1, -1
        w(i)   = w(i) + excess / deltz(i)
        excess = (w(i)-wsat) * deltz(i)
        excess = amax1(excess, 0.0)
        w(i)   = w(i) - excess /deltz(i)
        if(GWkey.eq.1 .and. abs(w(i)-wsat) .lt. 0.1E-5) then 
          Dgl = Dgl-deltz(i)
        else
          GWkey = 0
        endif
      end do

      Sst = Sst - inf*dt + excess

!     prevent negative values of www(i)                                         
      do i = 1,nuz-1                                                           
        deficit  = amax1 (0.0, (1.0e-10 - w(i)))                                  
        w(i)   = w(i) + deficit                                              
        w(i+1) = w(i+1) - deficit * deltz(i)/deltz(i+1)                 
      END DO
      w(nuz) = amax1(w(nuz), 1.0e-10)                                         


!    calculation of exchange rate between groundwater and river: qsub                            
      isat    = nuz+1
      qsub    = 0.0
      tmpqsub = 0.0

      if( (Dgl-Ds) .ge. 0.0 ) then ! groundwater table is bolow the toplayer
        avkg = kg
        call gwriv(Dgl,length,slope,Ds,Dg,Dr,Drw,avkg,tmpqsub) ! m3/s/m
        qsub = tmpqsub * dt / length                                   ! m
        if(abs(qsub) .lt. 0.1E-20) qsub = 0.0
      else                         ! groundwater table reaches the toplayer
        if(Drw.ge.Dr .and. Dgl.le.0.0) goto 1010
        tmp_D = 0.0
        avkg  = 0.0
        do i = nuz, 1, -1
          tmp_D = tmp_D + deltz(i)
          if( (Ds-tmp_D) .le. 0.5) then
            avkg  = avkg + anik*k0(i)*deltz(i)
          else
            avkg  = avkg + k0(i)*deltz(i)
          endif
          if(tmp_D .ge. (Ds-Dgl)) goto 50
        end do
50      if(Dr .gt. Ds) then
          tmp_D = tmp_D + Dr-Ds
          avkg  = avkg + kg*(Dr-Ds)
        endif
        avkg  = avkg / tmp_D
        call gwriv(Dgl,length,slope,Ds,Dg,Dr,Drw,avkg,tmpqsub) ! m3/s/m
        qsub = tmpqsub * dt / length                 ! m
        if(abs(qsub) .lt. 0.1E-20) qsub = 0.0
      endif

!     Renewing the groundwater table and toplayer soil moisture
      if(qsub .lt. 0.0) then      ! River infiltrates into aquifer
        GWst  = GWst - qsub
        if(GWst .le. Dg*GWcs) then
          Dgl   = Dgl + qsub/GWcs
          if(Dgl .lt. 0.0) then
            print *, 'wrong in renewing groundwater table ... 1', Dgl
          end if
        else
          tmp    = GWst - Dg*GWcs
          GWst   = Dg*GWcs
          Dgl    = Ds
          w(nuz) = w(nuz) + tmp / deltz(nuz)
          excess = (w(nuz) - wsat) * deltz(nuz)
          excess = amax1(excess, 0.0)
          w(nuz) = w(nuz) - excess / deltz(nuz)
          GWkey  = 0
          if(abs(w(nuz)-wsat) .lt. 0.1E-5) then 
            Dgl = Ds - deltz(nuz)
            GWkey = 1
          else 
            GWkey = 0
          endif
          do i=nuz-1, 1, -1
            w(i)   = w(i) + excess / deltz(i)
            excess = (w(i)-wsat) * deltz(i)
            excess = amax1(excess, 0.0)
            w(i)   = w(i) - excess /deltz(i)
            if(GWkey.eq.1 .and. abs(w(i)-wsat) .lt. 0.1E-5) then 
              Dgl = Dgl-deltz(i)
            else
              GWkey = 0
            endif
          end do
        endif
      elseif(qsub .gt. 0.0) then                ! Aquifer flows out
        if( Dgl .ge. Ds-0.01 ) then
          GWst  = GWst - qsub
          Dgl   = Dgl + qsub/GWcs
          if(Dgl.gt.(Ds+Dg) .or. GWst.lt.0.0) then
            qsub= qsub+GWst
            Dgl = Ds+Dg
            GWst= 0.0
          end if
        else
          tmp_D = 0.0
          isat  = nuz+1
          do i = nuz, 1, -1
            tmp_D = tmp_D + deltz(i)
            if( abs(tmp_D - (Ds-Dgl)) .le. 0.1E-3) goto 100
          end do
          goto 110
100       if(abs(w(i)-wsat) .le. 0.1E-3)   isat = i
          goto 120
110       if(abs(w(i+1)-wsat) .le. 0.1E-3) then
            isat = i+1
          end if   
120       tmp = 0.0
          do i = isat, nuz
            if(w(i)-wsat .gt. 0.1E-3) then
              print *, "wrong in GW flow", isat, i,Dgl,w(i),wsat
            end if
            tmp  = tmp + GWcs*deltz(i)
            w(i) = w(i) - GWcs
            Dgl  = Dgl + deltz(i)
            if(tmp .ge. qsub) goto 150
          end do
150       if(i.le.nuz) w(i) = w(i) + (tmp-qsub)/deltz(i)
          if(i.gt.nuz) then
            GWst  = GWst - (qsub-tmp)
            Dgl   = Dgl + (qsub-tmp)/GWcs
          end if
          if(i.gt.1) then
            if(  w(i-1)-wsat .gt. 0.1E-3) then
              print *,"w(i)>wsat", i-1, w(i-1),wsat,isat,nuz
            end if
            if(Dgl.gt.Ds .and. abs(GWst-Dg*GWcs).lt.0.1E-3) Dgl=Ds
          end if
        endif
      endif
 

!     Subsurface flow from top saturated zone above the groundwater level
      j = min(isat, nuz)
      ssub=0.0 ! added by gaobing
      tmp = 0.0
      do 200 i = 1, j
        tmp = tmp + deltz(i)
        if((w(i)-wfld) .gt. 0.1E-3) then
          ksoil(i) = conductivity_V(k0(i),wsat,wrsd,watern,w(i))
          tmpqsub  = ksoil(i)*slope*dt
          if(tmp .le. 0.5) tmpqsub = anik*tmpqsub
          f = 100.0
          if(tmp.le.1.5) f = -alog(0.1 /1.0)
          tmpqsub = tmpqsub * (deltz(i)+exp(-f*tmp)*ss_f*length)/length
          tmp     = (w(i)-wfld)*deltz(i)
          tmpqsub = amin1(tmpqsub, tmp)
          tmpqsub = amax1(tmpqsub, 0.0)
          w(i)    = w(i) - tmpqsub/deltz(i)
          qsub    = qsub + tmpqsub
          ssub    = ssub+tmpqsub  ! added by gaobing     !ssubÖ»°üÀšÁËSubsurface flow from top saturated zone above the groundwater level£¬¶øqsub»¹°üÀšÁËÖ®Ç°
        endif
200   continue
1010  ssub = ssub*length / dt              ! m --> m3/sec/m
      qsub = qsub*length / dt              ! m --> m3/sec/m
      return
      end subroutine 

!     calculate suction from soil moisture by Van Genuchten's equation                                         
      subroutine SuctionFromMoisture_V(w,wsat,wrsd,n,alpha,ps)
        implicit none
        real(r8),intent(in)::w
        real(r8),intent(in)::wsat
        real(r8),intent(in)::wrsd
        real(r8),intent(in)::alpha
        real(r8),intent(in)::n
        real(r8),intent(inout)::ps
        
        real(r8)::m
        real(r8)::se
        real(r8)::tmpe,tmpps,tmpw
      
        m    = 1.0-1.0/n
        tmpw = w
        if(tmpw.ge.wsat) then
          ps=0.0
        elseif(tmpw.lt.wsat) then
          if(tmpw.lt.wrsd+0.001) tmpw=wrsd+0.001
            se=(tmpw-wrsd)/(wsat-wrsd)
            tmpe=se**(1.0/m)
            tmpps=-(1.0/tmpe-1.0)**(1.0/n)/alpha
            ps = tmpps/100.0 ! cm->m   !!why??
          end if
        return
      end subroutine SuctionFromMoisture_V
      



!     calculate soil moisture from suction by Van Genuchten's equation                                         
      subroutine MoistureFromSuction_V(w,wsat,wrsd,n,alpha,ps)
        implicit none
        real(r8),intent(inout)::w          ! Volumetric water content
        real(r8),intent(in)::wsat       ! Saturated water content
        real(r8),intent(in)::wrsd       ! Residual water content
        real(r8),intent(in)::alpha      ! VG alpha
        real(r8),intent(in)::n          ! VG n
        real(r8),intent(in)::ps         ! suction,pressure
        real(r8)::se,tmpe,tmpps
        real(r8)::m          ! VG m
        
        tmpps=100.0*ps ! m->cm
        m=1.0-1.0/n
        tmpe=1.0+(alpha*abs(tmpps))**n
        se=(1.0/tmpe)**m
        w=se*(wsat-wrsd)+wrsd
        if(w .gt. wsat) w = wsat
        if(w .lt. wrsd) w = wrsd
        return
      end subroutine MoistureFromSuction_V
      


!    calculation of exchange rate between groundwater and river: qsub
     subroutine gwriv(Dtg,length,slope,Ds,Dg,Dr,Drw,kg,Q)
       implicit none
       real(r8),intent(in)::Dtg                !depth to groundwater (m)
       real(r8),intent(in)::length             !length of hillslope (m)
       real(r8),intent(in)::slope              !slope of hillslope (m)
       real(r8),intent(in)::Ds                 !depth of top soil (m)
       real(r8),intent(in)::Dg                 !depth of unconfined groundwater acquifer below topsoil (m)
       real(r8),intent(in)::Dr                 !depth of river (m)
       real(r8),intent(in)::kg                 !hydraulic conductivity (m/sec)
       real(r8),intent(inout)::Q                  !discharge exchanged between aquifer and river (m^3/sec)! discharge per unit width of a hillslope (m3/sec/m)
       real(r8),intent(inout)::Drw                !depth of river water (m)

       real(r8)::H1                 !waterhead of groundwater (m)
       real(r8)::H2                 !waterhead of river (m)
       real(r8)::hs1                !saturated acquifer depth (m)
       real(r8)::hs2                !water depth in river (m)
       real(r8)::Hrd                !distance from datum to riverbed (m)
       real(r8)::conlen             !contact length between the river and aquifer (m)
       real(r8)::grad               !gradient of water head
       !The datum is sitted at the bottom of the unconfined aquifer
       Drw = amax1(Drw, 0.0)
       Drw = amin1(Dr, Drw)
       Hrd = Ds+Dg-Dr
       if(Dtg .lt. Dr-Drw) then
         H1  = 0.5*length*slope + Ds+Dg-Dtg
       else         !feichengyashui
         H1  = sqrt(0.5*length*slope) + Ds+Dg-Dtg
       endif
       hs1=H1-Hrd
       if(Hrd.ge.0.0) then
         H2 =Hrd+Drw
         hs2=Drw
       else
        hs2=amax1(0.0, Drw)
        H2 =amax1(Hrd+Drw, 0.0)
       endif
       grad   =(H1-H2)/(0.5*length)
       conlen =0.5*(abs(hs1)+hs2)
       Q=kg*grad*conlen
       if(Drw .le.5E-3 .and. Q.lt.0.0) Q=0.0
       return
     end subroutine



!                                                                *
!     This routine solves a linear system of equations           *
!                  A * X = RS                                    *
!     for a tridiagonal, strongly nonsingular matrix A.          *
!     The matrix is given by the three vectors DL, DM and DU     *
!     which designate the lower co-diagonal, the diagonal and    *
!     the upper co-diagonal elements of A, respectively.         *
!     The system of equations has the form :                     *
!                                                                *
!     DM(1) * X(1) + DU(1) * X(2)                      = RS(1)   *
!                                                                *
!     DL(I) * X(I-1) + DM(I) * X(I) + DU(I) * X(I+1)   = RS(I)   *
!            for I = 2, ..., N-1, and                            *
!                                                                *
!     DL(N) * X(N-1) + DM(N) * X(N)                    = RS(N)   *
!                                                                *
!                                                                *
!     INPUT PARAMETERS:                                          *
!     =================                                          *
!     N    : number of equations, N > 2                          *
!     DL   : N-vector DL(1:N); lower co-diagonal of A            *
!            DL(2), DL(3), ... ,DL(N)                            *
!     DM   : N-vector DM(1:N); the diagonal of A                 *
!            DM(1), DM(2), ... , DM(N)                           *
!     DU   : N-vector DU(1:N); upper co-diagonal of A            *
!            DU(1), DU(2), ... , DU(N-1)                         *
!     RS   : N-vector RS(1:N); the right hand side               *
!                                                                *
!                                                                *
!     OUTPUT PARAMETERS:                                         *
!     ==================                                         *
!     DL   :)                                                    *
!     DM   :)                                                    *
!     DU   :) these are overwritten with auxiliary vectors       *
!     RS   :)                                                    *
!     X    : N-vector X(1:N), containing the solution of the     *
!            system of equations                                 *
!     MARK : error parameter                                     *
!            MARK= 1 : everything is o.k.                        *
!            MARK= 0 : the matrix A is not strongly nonsingular  *
!            MARK=-1 : error on N: N <= 2                        *
!                                                                *
!     NOTE: if MARK = 1, the determinant of A can be calculated  *
!           after this subroutine has run as follows:            *
!              DET A = DM(1) * DM(2) * ... * DM(N)               *
!                                                                *
!----------------------------------------------------------------*
      subroutine TRDIG(N,DL,DM,DU,RS,X,mark)
        implicit none
        integer(i4),intent(in)::N
        real(r8),intent(inout)::DL(N)
        real(r8),intent(inout)::DM(N)
        real(r8),intent(inout)::DU(N)
        real(r8),intent(inout)::RS(N)
        real(r8),intent(inout)::X(N)
        integer(i4)::mark,markmark,i
        real(r8)::row,d
        MARK = 0
        markmark=0
        row  = abs(DM(1)) + abs(DU(1))
        if (N.Ge.3 .and. row .ne. 0.0)then
           D=1.0E0/ROW
           if(abs(DM(1))*D .gt. 1.0E-20) then
!          checking for strong nonsingularity with N=1
!          factoring A while checking for strong nonsingularity
             DL(1) = 0.0E0
             DU(N) = 0.0E0
             DU(1) = DU(1)/DM(1)
             do I=2,N,1
               row = abs(DL(I)) + abs(DM(I)) + abs(DU(I))
               if (row .eq. 0.0E0) then
                 exit
                 markmark=1
               end if
               D = 1.0E0/ROW
               DM(I) = DM(I) - DL(I) * DU(I-1)
               if (ABS(DM(I))*D .LE. 1.0E-20) then
                 exit
                 markmark=1
               end if
               if (I .LT. N) then
                 DU(I) = DU(I)/DM(I)
               endif
             enddo
           if (markmark.ne.1) then
             MARK=1
           else
             mark=0
           endif
        end if
      end if
!     If MARK = 1, update the right hand side and solve via backsubstitution
      if (MARK .eq. 1) then
        RS(1) = RS(1)/DM(1)
        do I=2,N,1
          RS(I) = (RS(I) - DL(I) * RS(I-1)) / DM(I)
        enddo
!       backsubstitution
        X(N) = RS(N)
        do I=N-1,1,-1
          X(I) = RS(I) - DU(I) * X(I+1)
        enddo
      end if
      return
    end subroutine
    

      subroutine ETsoil(cm,cmf,cmw,Kcm)
        implicit none
        real(r8),intent(in)::cm,cmf,cmw
        real(r8),intent(inout)::Kcm
        real(r8)::cmc1,cmc2
        real(r8)::r,c
!       if(cm.lt.cmw) cm = cmw
!       Kcm = (cm-cmw)/(cmf-cmw)
!       if(cm .ge. cmf)  Kcm=1.0
        c=0.5
        r=0.5
        Kcm=0.0
        cmc2=r*(cmf+cmw)
        cmc1=(cmc2-c*cmf)/(1.0-c)
        if(cm .ge. cmf)  Kcm=1.0
        if(cm.lt.cmf .and. cm.ge.cmc2) Kcm=(cm-cmc1)/(cmf-cmc1)
        if(cm.lt.cmc2 .and. cm.gt.cmw) Kcm=c*(cm-cmw)/(cmc2-cmw)
        if(Kcm.lt.0.0 .or. Kcm.gt.1.0) then
          print *,'wrong in calculating Kcm in ETsoil',cm,cmw,cmf,Kcm
        end if
      end subroutine
      
     function ran0(idum)
       implicit none
       integer(i4),intent(inout)::idum
       integer(i4)::IA,IM,IQ,IR,MASK
       real(r8)::ran0
       real(r8)::AM
       integer(i4)::k
       IA=16807
       IM=2147483647
       IQ=127773
       IR=2836
       MASK=123459876
       AM=1./IM
       idum=ieor(idum,MASK)
       k=idum/IQ
       idum=IA*(idum-k*IQ)-IR*k
       if (idum.lt.0) idum=idum+IM
       ran0=AM*idum
       idum=ieor(idum,MASK)
       return
     end function
     
     subroutine rexcess(re,r,k0)
       real k0,re,r
       if((k0/r).le.5.0) then
         re=r*exp(-1.0*k0/r)
         re=amin1(re,r)
       else
         re=0.0
       end if
       return
     end subroutine


!    Sub-program -- downscaling forcing data daily -> hourly
      subroutine rain_model(nrow,ncol,isub,ir,ic,iland,idh,rain_daily, tmin_daily, &
           tmax_daily, evap_daily,r,T,Ep)
        implicit none
        integer(i4),intent(in)::isub, ir, ic, iland,nrow,ncol
        integer(i4),intent(in)::idh
        real(r8),intent(in):: rain_daily
        real(r8),intent(in)::tmin_daily, tmax_daily, evap_daily
        real(r8),intent(inout):: r, T, Ep
                
        integer(i4)::rt0(nrow,ncol)
        real(r8)::tmp
        integer(i4)::idum, tmd


        if(evap_daily.lt.0.0) then
          print *,"daily ET is less than 0",isub,ir,ic,evap_daily
          stop
        end if
        if(tmin_daily.eq.-9999) then
          print *,"min temperature is  -9999",isub,ir,ic,tmin_daily
          stop
        end if
        if(tmax_daily.eq.-9999) then
          print *,"max temperature is  -9999",isub,ir,ic,tmax_daily
          stop
        end if
!------------------------ create the start time ------------
        if(idh .eq. 1 .and. iland.eq.1) then
          tmp = ran0(idum)
          rt0(ir,ic) = 4 + int((21-4)*tmp+0.45)
          if(rt0(ir,ic) .gt. 21) rt0(ir,ic) = 21
        endif
               
        Ep = 0.0
        if(idh.eq.6.or.idh.eq.17)  Ep = evap_daily*0.03
        if(idh.eq.7.or.idh.eq.16)  Ep = evap_daily*0.06
        if(idh.eq.8.or.idh.eq.15)  Ep = evap_daily*0.08
        if(idh.eq.9.or.idh.eq.14)  Ep = evap_daily*0.10
        if(idh.eq.10.or.idh.eq.13) Ep = evap_daily*0.11
        if(idh.eq.11.or.idh.eq.12) Ep = evap_daily*0.12
!------------------------ downscaling Prec --------------------
        r = 0.0
        if(rain_daily.ge.0.0 .and. rain_daily.le. 5.0) then
          if(idh.eq.rt0(ir,ic))     r = rain_daily
        elseif(rain_daily.gt.5.0 .and. rain_daily.le.10.0) then
          if(idh.eq.(rt0(ir,ic)-1)) r = 0.30*rain_daily
          if(idh.eq.rt0(ir,ic))     r = 0.40*rain_daily
          if(idh.eq.(rt0(ir,ic)+1)) r = 0.30*rain_daily
        elseif(rain_daily.gt.10.0 .and. rain_daily.le.30.0) then
          if(idh.eq.(rt0(ir,ic)-3)) r = 0.10*rain_daily
          if(idh.eq.(rt0(ir,ic)-2)) r = 0.15*rain_daily
          if(idh.eq.(rt0(ir,ic)-1)) r = 0.15*rain_daily
          if(idh.eq.rt0(ir,ic))     r = 0.20*rain_daily
          if(idh.eq.(rt0(ir,ic)+1)) r = 0.15*rain_daily
          if(idh.eq.(rt0(ir,ic)+2)) r = 0.15*rain_daily
          if(idh.eq.(rt0(ir,ic)+3)) r = 0.10*rain_daily
        elseif(rain_daily.gt.30.0) then
          r = rain_daily/24.0
        end if
!------------------------ downscaling Temp --------------------
        if(idh.ge.4 .and. idh.lt.14) then
          T = tmin_daily+ (tmax_daily-tmin_daily)*(1-cos((idh-4)*3.1415926/10.0))/2.0
        else
          tmd = idh
          if(idh.lt.4) tmd = 24+idh
          T = tmin_daily + (tmax_daily-tmin_daily)*(1+cos((tmd-14)*3.1415926/14.0))/2.0
        end if
        return
      end subroutine rain_model
     
     
      subroutine nkws(dt,dx, b, s0, roughness,qlin0, qlin, Q01,Q02,Q1, Q2, y)
        implicit none
        real(r8),intent(in)::Q01                   !known,discharge of last time step
        real(r8),intent(in)::Q02                   !know,discharge of last time step
        real(r8),intent(in)::Q1                    !known,discharge of current time step
        real(r8),intent(in)::qlin0                 !lateral inflow of last time step
        real(r8),intent(in)::qlin                  !lateral inflow of current time step
        real(r8),intent(in)::dt                    ! time interval
        real(r8),intent(in)::dx                    !width of flow interval
        real(r8),intent(in)::b                     !river width
        real(r8),intent(in)::s0                    !river slope
        real(r8),intent(in)::roughness
        real(r8),intent(inout)::y                     !water depth
        real(r8),intent(inout)::Q2                    !unknown, discharge of current time step
        
        real(r8)::p                     !wetted perimeter
        real(r8)::beta
        real(r8)::criterion
        
        integer(i4)::k
        real(r8)::h1,tmp,f,aa,bb,qq1,alfa,cc,ctmp,df,h2,qq2
        
        beta = 0.6
        criterion = 0.0001
        k   = 1
        h1  = Q02/b
        do k=1,30
          tmp = roughness*Q02/sqrt(s0)
          f   = b*h1 - (tmp**0.6)*((b+2.0*h1)**0.4)
          if(k.gt.1 .and. abs(f) .lt. criterion) then
            exit
          end if
          df = b - 0.8*((tmp/(b+2.0*h1))**0.6)
          h2 = h1 - f/df
          h1 = h2
        end do
        y = h2
        p = b+2.0*y ! for rectangular channel
!       the initial discharge estimated using linear scheme
        alfa = (roughness*p**(2.0/3.0)/sqrt(s0))**0.6
        if((Q02+Q1).le.0.0) then
          cc   =0
        else
          cc   = (0.5*(Q02+Q1))**(beta-1.0)
        endif
        aa   = dt*Q1/dx + alfa*beta*Q02*cc + 0.5*(qlin0+qlin)*dt
        bb   = dt/dx + alfa*beta*cc
        qq1  = aa/bb

        if(qq1 .le. 0.1e-5) then
          qq1=0.1e-5 
        end if
!       Using Newton's method to calculate discharge
        ctmp = dt*Q1/dx + alfa*Q02**beta + 0.5*dt*(qlin+qlin0)
        do k=1,30
          f= dt*qq1/dx + alfa*qq1**beta-ctmp
          if((k.gt.1).and.(abs(f).le.criterion)) then
            exit
          end if
          df  = dt/dx + alfa*beta*qq1**(beta-1.0)
          qq2 = qq1 - f/df
          if(qq2 .le. 0.1e-5) then
            qq2=0.0
            exit
          endif
          qq1 = qq2
        end do
        Q2 = qq2
        
        k   = 1
        h1  = Q2/b
        do k=1,30
          tmp = roughness*Q2/sqrt(s0)
          f   = b*h1-(tmp**0.6)*((b+2.0*h1)**0.4)
          if(k.gt.1 .and. abs(f) .lt. criterion) then
            exit
          end if
          df = b-0.8*((tmp/(b+2.0*h1))**0.6)
          h2 = h1-f/df
          h1 = h2          
        end do
        y = h2
        return
      end subroutine nkws

    end module hydro_mod


