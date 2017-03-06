    module GBHM3Lib_mod
      use global_para_mod,only:i4,r8
    contains
    
    subroutine strlen(str,l1,l2)
      implicit none
      character(200)::str
      integer::i,l1,l2,k
      k=0
      do i = 1, 200
        if(k.eq.0 .and. str(i:i).NE.' ') then
            l1=i
            k=1
        elseif(k.eq.1 .and. str(i:i).EQ.' ') then
            l2 = i-1
            exit
        endif
      end do
    end subroutine
      

    subroutine read_river_para(para_dir)
      use hydro_data_mod,only:nflow,ngrid,dx,s0,b,roughness,dr,grid_row,grid_col,ngridmax,subbasin,nsub
      implicit none
      character*200,intent(in)::para_dir
      integer(i4)::isub,l1,l2,ll1,ll2,iflow,j
      character*200::infile
      do isub = 1, nsub
          infile = subbasin(isub)//'_river'
          print *,isub,subbasin(isub),"  ",infile
          call strlen(infile, l1, l2)
          call strlen(para_dir, ll1, ll2)
          open(14,file=trim(para_dir(ll1:ll2))//'riverpara/'//infile(l1:l2),status='old')
          read(14,*) nflow(isub)
          do iflow = 1, nflow(isub)
            read(14,*) ngrid(isub,iflow),dx(isub,iflow),s0(isub,iflow),b(isub,iflow),roughness(isub,iflow),Dr(isub,iflow)
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
            read(14,*) (grid_row(isub,iflow,j),grid_col(isub,iflow,j),j = 1, ngrid(isub,iflow))
            print*, (grid_row(isub,iflow,j),grid_col(isub,iflow,j),j = 1, ngrid(isub,iflow))
          end do
          close(14)
      end do
    end subroutine read_river_para

    subroutine read_hydro_para(para_dir,nhrpdt)
      use global_para_mod,only:i4,r8
      use hydro_data_mod,only:dthydro
      implicit none
      character*200,intent(in)::para_dir
      integer(i4),intent(in)::nhrpdt
      dthydro=nhrpdt*3600.0
      !call initial_subcatchment(para_dir)
      call read_river_para(para_dir)
      return 
    end subroutine read_hydro_para


    subroutine read_veg_para(para_dir)
      use global_para_mod,only:i4,r8
      use land_para_mod,only:nland,Kcanopy,LAImax,kcrop,root
      use hydro_data_mod,only:Sstmax,surfn
      use soil_para_mod,only:anik
      implicit none
      character*200,intent(in)::para_dir
      integer(i4)::l1,l2,subcount,iostatus
      character(200)::atmp,landusename
      integer(i4)::i,iland
      call strlen(para_dir, l1, l2)
      open(14, file = para_dir(l1:l2)//'vege_para.dat',status='old')
      read(1,*)
      read(1,*)
      i = 1
      do i=1,nland
        read(1,*) iland, Kcanopy(i), LAImax(i), kcrop(i),root(i),anik(i), Sstmax(i), surfn(i), landusename
      end do
      close(14)
      print *,'landuse types:',nland
      print *,'ok2, finish of vegetation.dat', para_dir
    end subroutine read_veg_para


    subroutine initial_subcatchment(para_dir)
      use global_para_mod,only:i4,r8
      use hydro_data_mod,only:nsub,ngridmax,nflowmax,nrow,ncol,subbasin,psubbasin,pbasinup,nbasinup,&
            nflow,ngrid,dx,s0,b,roughness,Dr,grid_row,grid_col,Drw,q1,q2,qlin1,qlin2,qin,area_sub,&
            nhrpdt,area,ele,slp,length,ds,dg,layer,d,k0,w,kg,gwcs,cst,sst,dgl,gwst,nlayer
      use land_para_mod,only:nv,nland,NDVI,yndvi,LAI,Kcanopy,root,landtyp,land_ratio,land_use,Kcrop
      use soil_para_mod,only:soiltyp,soil,wsat,wrsd,wfld,alpha,watern,ksat1,ksat2,anik,ns
      implicit none
      character*200,intent(in)::para_dir
      integer::i,j,k,fileunit
      character(200)::atmp
      integer(i4)::l1,l2,subcount,iostatus
      character*4::ch4 
      character*200::subbasinfile
      character*(10)::atemp
      call strlen(para_dir,l1,l2)
      subbasinfile=trim(para_dir(l1:l2))//'subbasin.dat'
      open(14,file = subbasinfile, status='old')
      subcount=0
      iostatus=1
      do while(iostatus.ge.0)
        read(14,*,iostat=iostatus) atmp
        subcount=subcount+1
      end do
      close(14)
      subcount=subcount-1
      nsub=subcount
      call strlen(para_dir, l1, l2)
      open(14,file=trim(para_dir(l1:l2))//'riverpara/'//'maxnp',status='old')
      read(14,*) nflowmax,ngridmax
      close(14)
      ! get the nrow and ncol from the ws.asc file
      call strlen(para_dir, l1, l2)
      open(14,file=trim(para_dir(l1:l2))//'ws.asc',status='old')
      read(14,*) atemp,ncol
      read(14,*) atemp,nrow
      close(14)
      open(14, file = para_dir(l1:l2)//'vege_para.dat',status='old')
      read(14,*)
      read(14,*)
      subcount=0
      iostatus=1
      do while(iostatus.ge.0)
        read(14,*,iostat=iostatus) atmp
        subcount=subcount+1
      end do
      close(14)
      subcount=subcount-1
      nland=subcount
      
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
      allocate(Drw(nsub,nflowmax))
      allocate(q1(nsub,nflowmax))
      allocate(q2(nsub,nflowmax))	
      allocate(qlin1(nsub,nflowmax))
      allocate(qlin2(nsub,nflowmax))
      allocate(qin(nsub))
      allocate(area_sub(nsub))
      if(.not. allocated(area)) allocate(area(nrow,ncol))
      if(.not. allocated(ele)) allocate(ele(nrow,ncol))
      if(.not. allocated(slp)) allocate(slp(nrow,ncol))
      if(.not. allocated(length)) allocate(length(nrow,ncol))
      if(.not. allocated(Ds)) allocate(Ds(nrow,ncol))
      if(.not. allocated(Dg)) allocate(Dg(nrow,ncol))
      if(.not. allocated(layer)) allocate(layer(nrow,ncol))
      if(.not. allocated(D)) allocate(D(nrow,ncol,nlayer))
      if(.not. allocated(k0)) allocate(k0(nrow,ncol,nlayer))
      if(.not. allocated(w)) allocate(w(nrow,ncol,nv,nlayer))
      if(.not. allocated(kg)) allocate(kg(nrow,ncol))
      if(.not. allocated(GWcs)) allocate(GWcs(nrow,ncol))
      if(.not. allocated(Cst)) allocate(Cst(nrow,ncol,nv))
      if(.not. allocated(Sst)) allocate(Sst(nrow,ncol,nv))
      if(.not. allocated(Dgl)) allocate(Dgl(nrow,ncol,nv))
      if(.not. allocated(GWst)) allocate(GWst(nrow,ncol,nv))
      
      ! Allocate the space for the vegetation variables
      if(.not. allocated(NDVI)) allocate(NDVI(nrow,ncol,12,3))
      if(.not. allocated(yndvi)) allocate(yndvi(nrow,ncol,12,3))
      if(.not. allocated(LAI)) allocate(LAI(nrow,ncol,12,4))
      if(.not. allocated(Kcanopy)) allocate(Kcanopy(nland))
      if(.not. allocated(root)) allocate(root(nland))
      if(.not. allocated(landtyp)) allocate(landtyp(nland))
      if(.not. allocated(land_ratio)) allocate(land_ratio(nrow,ncol,nland))
      if(.not. allocated(land_use)) allocate(land_use(nrow*10,ncol*10,1))
      if(.not. allocated(Kcrop)) allocate(Kcrop(nland))

      ! Allocate the space for the soil parameter variables
      if(.not. allocated(soiltyp)) allocate(soiltyp(ns))
      if(.not. allocated(soil)) allocate(soil(nrow,ncol))
      if(.not. allocated(wsat)) allocate(wsat(nrow,ncol))
      if(.not. allocated(wrsd)) allocate(wrsd(nrow,ncol))
      if(.not. allocated(alpha)) allocate(alpha(nrow,ncol))
      if(.not. allocated(wfld)) allocate(wfld(nrow,ncol))
      if(.not. allocated(watern)) allocate(watern(nrow,ncol))
      if(.not. allocated(ksat1)) allocate(ksat1(nrow,ncol))
      if(.not. allocated(ksat2)) allocate(ksat2(nrow,ncol))
      if(.not. allocated(anik)) allocate(anik(nland))

      open(14,file = subbasinfile, status='old')
      do i = 1,nsub
        read(14,*)nsub,psubbasin(i),(pbasinup(i,j),j=1,8)!,nextbasinrow(i),nextbasincol(i),nextbasin(i)
      enddo
      close(14)
      
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
      print *,'number of sub-catchment:',nsub
      print *,"the nrow and ncol are",nrow,ncol
      print * ,'ok1, finish of subbasin.dat', para_dir
    end subroutine
    
    end module GBHM3Lib_mod
