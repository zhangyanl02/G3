      subroutine river_routing(isub)
      implicit none

      character*6   subbasin(nc)      ! name of sub-basin changed by gaobing 2013
      integer       psubbasin(nc)                   ! added by gaobing  2013
      integer       nbasinup(nc,4)
      integer       nsub              ! total number of sub-basins
      integer       isub              ! sub-basin number
      integer nflow(nc)    ! total number of flow-intervals in a sub-basin
      integer iflow        ! flow-interval number
      real    dx(nc,nx)    ! flow interval length (m)
      real    Dr(nc,nx)    ! river depth in the flow interval (m)
      real    Drw(nc,nx)   ! water depth in the river of the flow interval (m)
      real    s0(nc,nx)        ! slope of river bed (ND)
      real    b(nc,nx)         ! width of river (m)
      real    roughness(nc,nx) ! Manning's roughness
      real    dt             ! time step (hour)
      integer year           ! calendar year   (19XX)
      integer month          ! month in a year (1-12)
      integer day            ! day in a month  (1-30)
      integer hour           ! hour in a day   (1-24)
      integer dayinmonth(12) ! days of a month
      integer hydroyear      ! year for hydro simulation
      integer startmonth     ! start month in the hydro-year
      integer endmonth         ! end   month in the hydro-year
      integer startday       ! start day in the first month of the hydro-year
      integer endday         ! end   day in the last  month of the hydro-year
      integer iy, im, id ,ih ! for year, month, day and hour loop
      integer idc            ! contineous day in a year (1-366)
      integer ihc            ! contineous hour in a year (1-366*24)

      integer start, finish  ! 0 - faulse,    1 - true
      integer start_sub, end_sub  ! for partial simulation of a basin
c
c#######################################################################
c
c     Discharge variables
c
c#######################################################################
c
      real qin(nx)         ! lateral inflow of one flow interval (m3/sec/m)
                           !, from 'hydro.inc'
c
c      real beta            ! Manning's equation parameter
      real q1(nc,nx)       ! discharge of last time step (m^3/s)
      real q2(nc,nx)       ! discharge of current time step (m^3/s)
      real qr1(nc,nx)     ! discharge of reservoir flowout last time step(m^3/s)
      real qr2(nc,nx)  ! discharge of reservoir flowout current time step(m^3/s)
      real qlin1(nc,nx)    ! lateral inflow of last time step (m^3/m/s)
      real qlin2(nc,nx)    ! lateral inflow of current time step (m^3/m/s)
      real Qh(nc,8800)     ! hourly mean discharge (m3/s)
      common /qh/ qh   !mqh
      real Qd(nc,366)           ! daily average discharge (m^3/s)
      real Qm(nc,12)       ! monthly mean discharge (m3/s)
      real Qobs(nc,366)    ! added by gaobing
cc      real re_capacity(nc,nx) !capacity of reservoir in location(isub,iflow), m3
cc      real re_capacity0(nc,nx)  ! initial capacity of reservoir 

      character*200 gauge_name    ! name of discharge gauges
      character*200 gauge_name2    ! name of discharge gauges 
c
c      save q1       ! discharge of last time step (m^3/s)
      save q2       ! discharge of current time step (m^3/s)
      save qr1      ! discharge of reservoir flowout last time step(m^3/s)
c      save qr2      ! discharge of reservoir flowout current time step(m^3/s)
      save qlin1    ! lateral inflow of last time step (m^3/m/s)
      save qlin2    ! lateral inflow of current time step (m^3/m/s)
!      save Qh       ! hourly mean discharge (m3/s)
      save Qd            ! daily average discharge (m^3/s)
      save Qm       ! monthly mean discharge (m3/s)
      save Qobs     ! added by gaobing 

      integer l1, l2, l3
      integer level
      integer p1(9), p2(9), p3(9)

      character*200  para_dir    ! directory for parameters
      character*200  dem_dir    ! directory for dem
      character*200  data_dir    ! directory for input data
      character*200  result1_dir ! directory for storing simulation result
      character*200  result2_dir ! directory for storing simulation result
      character*200  simul_dir   ! directory for storing temporal result

      character*4    ch4         ! a 4-character variable
      character*3    ch3
      real    tmp, value         ! temporary variables
      integer i, j, k, m ,l        !, n      ! temporary variables
c      integer jj, is, il         ! temporary varaibles
      integer ia1, ia2, ib1, ib2, ic1, ic2,id1,id2
      real    criterion, h1, h2, f, df
      real    q1_upper, q2_upper,rs    !, i_sub,  lo_dis
      integer ii1, ii2           ! , ii3
      real Qtotal,basinarea
      real      b_tmp              ! adjust riverbed-width during heavy-flood
      real      rs_tmp             ! adjust rs during heavy-flood
                                             ! added in 2006.4.30

      common /date1  / year, month, day, hour, idc,ihc
           common /date2  / hydroyear,startmonth,startday,endmonth,endday
      common /date3  / dayinmonth
      
      common /simulation/  start, finish, dt, start_sub, end_sub
      common /river1    /  nsub, subbasin 
      common /river2    /  nflow, dx, Dr
      common /river3    /  s0, b, roughness
      common /river4    /  psubbasin,nbasinup ! added by gaobing 2013
      common /location  /  para_dir,data_dir,result1_dir,dem_dir,
     :                      result2_dir,simul_dir
      common /TEST   /basinarea

      common /pfaf/p1,p2,p3          !Pfaftetter basin number

      common /lateralflow/  qin
      common /river_status/ Drw, qr2  ! depth of river water 

      integer im2,idd2,hour2,jj

      if(month.eq.startmonth .and. day.eq.startday .and. hour.eq.1) then
        do i = 1, nsub
          do j= 1, 366
            Qd(i,j)=0.0
          end do
        end do
        do i = 1, nsub
           do j= 1, 12
               Qm(i,j)=0.0
           end do
        end do 
        do i = 1, nsub
           do j= 1, 366
               Qobs(i,j)=0.0
           end do
        end do
      endif

      if(start .eq. 1) then
        do i = 1, nsub
          do j = 1, nflow(isub)
            qlin1(i,j) = 0.0
          end do
        end do

        if(inicon.eq.1) then
          call strlen(simul_dir, ia1,ia2)
          open(3,file = simul_dir(ia1:ia2)//subbasin(isub)//'I_flow2',status='old')
          read (3,*)(i,q1(isub,j), j=1, nflow(isub))
          close(3)
        else
          q1(isub,1)=0.2    !mqh,changed from 0.5 to 0.2
          do iflow = 2, nflow(isub)
            q1(isub,iflow) = q1(isub,iflow-1)+0.1  !mqh,changed from 0.4 to 0.1
          end do
        endif
        
        
        do iflow=1,nflow(isub)
          criterion=0.01
          k=1
          h1=q1(isub,iflow)/b(isub,iflow) !b: river width of the flow interval (m)
          do k=1,10
            tmp=roughness(isub,iflow)*q1(isub,iflow)/sqrt(s0(isub,iflow))
            f=b(isub,iflow)*h1-(tmp**0.6)*((b(isub,iflow)+2.0*h1)**0.4)
            if (k .GT. 1 .AND. abs(f) .LT. criterion) then
              exit
            else
              df=b(isub,iflow)-0.8*((tmp/(b(isub,iflow)+2.0*h1))**0.6)
              h2=h1-f/df
              h1=h2
            end if
          end do
          h=h2
          hhr(isub,iflow)=h    ! depth of river water
        end do
        
        do iflow = 1, nflow(isub)
          criterion = 0.01
          k = 1
          h1  = q1(isub,iflow) / b(isub,iflow)
          do k=1,10
            tmp = roughness(isub,iflow) * q1(isub,iflow)/sqrt(s0(isub,iflow))
            f   = b(isub,iflow)*h1 - (tmp**0.6)*((b(isub,iflow)+2.0*h1)**0.4)
            if(k.gt.1 .and. abs(f) .lt. criterion) then
              exit
            else
              df = b(isub,iflow)-0.8*((tmp/(b(isub,iflow)+2.0*h1))**0.6)
              h2 = h1-f/df
              h1 = h2
          end do
          Drw(isub,iflow) = h2            ! print *,"initial h=",h2
        end do
        do iflow = 1, nflow(isub)
          qr1(isub,iflow) = q1(isub,iflow)
        end do
      endif 

      if(start .eq. 1) return
      if(month.eq.startmonth.and.day.eq.startday.and.ihc.eq.1) then
        do iflow = 1, nflow(isub)
           qlin1(isub,iflow)= qin(iflow)/dx(isub,iflow)
        end do
      end if

      do iflow = 1, nflow(isub)
        qlin2(isub,iflow) = qin(iflow)/dx(isub,iflow)
      end do

      m = 1
      do 555 j = 1, m                  !time loop
        do 333 iflow = 1, nflow(isub)  !river segment loop
          q1_upper=0.
          q2_upper=0.
          if(iflow .eq. 1) then  ! changed by gaobing 2013
            if(psubbasin(isub) > 1000 .and. psubbasin(isub) < 2000) then
              q1_upper=0.
              q2_upper=0.
            else
              do ii1 =1,4
                if(nbasinup(isub,ii1) > 0) then
                  ii2=nbasinup(isub,ii1)
                  q1_upper=q1_upper+qr1(ii2,nflow(ii2))
                  q2_upper=q2_upper+qr2(ii2,nflow(ii2))
                endif
              enddo
            endif
          else
            q1_upper=qr1(isub,iflow-1)
            q2_upper=qr2(isub,iflow-1)
          endif               !end of river network definition
          rs=roughness(isub,iflow)
          if(s0(isub,iflow).eq.0) then   ! added by xujijun 20051104
            s0(isub,iflow)=0.00001
            print *, 'wrong in s0',s0(isub,iflow)  
          endif  
          if(s0(isub,iflow).eq.-9999.0) then ! added by xujijun 20051104
            s0(isub,iflow)=0.00001 
            print *, 'wrong in s0',s0(isub,iflow)  
          endif
          b_tmp = b(isub,iflow)
          rs_tmp = rs
          if(Drw(isub,iflow) .gt. 0.5*Dr(isub,iflow)) then            ! ÊéÇ©2
            b_tmp=1.1*b(isub,iflow)
            rs_tmp = 1.5*rs
          endif
          if(Drw(isub,iflow) .gt. 1.0*Dr(isub,iflow)) then
            b_tmp=1.5*b(isub,iflow)
            rs_tmp = 3.0*rs
          endif
          call nkws(dt,dx(isub,iflow),b_tmp,s0(isub,iflow),rs_tmp,qlin1(isub,iflow),qlin2(isub,iflow),q1_upper,q1(isub,iflow),&
             q2_upper,q2(isub,iflow),Drw(isub,iflow))
          qr2(isub,iflow)=q2(isub,iflow)
        end do
        do iflow = 1, nflow(isub)
          qlin1(isub,iflow) = qlin2(isub,iflow)
          q1(isub,iflow)    = q2(isub,iflow)
          qr1(isub,iflow)   = qr2(isub,iflow)
        end do
      end do


      Qh(isub,ihc)  = qr2(isub,nflow(isub))
      Qd(isub,idc)  = Qd(isub,idc) + qr2(isub,nflow(isub))/24.0
      Qm(isub,month)= Qm(isub,month) + qr2(isub,nflow(isub))/(float(dayinmonth(month))*24.0)
      Qobs(isub,idc)=Qobs(isub,idc)+qr2(isub,nflow(isub))/24.0
      Qtotal=0.0 

      if(month.eq.endmonth .and. day.eq.endday .and. hour.eq.24) then
        if(subbasin(isub) .eq. 'ws4017') then ! changed by xujijun   !mqh
          do i=1,idc
            Qtotal=Qtotal+Qd(isub,i)*3600*24*1000.0/basinarea
          end do  
        end if
        call strlen(simul_dir, ia1,ia2)
        open(9,file = simul_dir(ia1:ia2)//subbasin(isub)//'I_flow2',status='replace')
        do iflow = 1, nflow(isub)
          write(9,110) iflow, q2(isub,iflow)
 110      format(2x,i4,f20.3)
        end do
        close(9)
        gauge_name = 'no_gauge'
        if(subbasin(isub) .eq. 'ws4017') gauge_name = 'maduwang'  !mqh

        call strlen(result1_dir,ib1,ib2)
        if(year.eq.startyear ) then
          open(9,file = result1_dir(ib1:ib2)//subbasin(isub)//'.daily', status='unknown')
        else
          open(9,file = result1_dir(ib1:ib2)//subbasin(isub)//'.daily', access='append', status='old')
        endif
        j = 0
        do im = startmonth, endmonth
          do id = 1, dayinmonth(im)
            j = j+1
            write(9, '(3i6,f16.3)') year,im,id, Qd(isub,j)  !±êŒÇ1£¬Qd(isub,j)
          end do
        end do
        close(9)
      end if
        
      if(year.ge.2000) then
        call strlen(result1_dir,ib1,ib2)
        write(ch3, '(i3.3)') isub 
        if(year.eq.2000.and.ihc.eq.1 ) then
          open(19,file = result1_dir(ib1:ib2)//'q_'//ch3//'.hourly', status='unknown')
        else
          open(19,file = result1_dir(ib1:ib2)//'q_'//ch3//'.hourly', access='append', status='old')
        endif
        write(19, '(4i6,f16.6)') year,month,day,hour,Qh(isub,ihc) 
        close(19)
      endif
    end subroutine river_routing
