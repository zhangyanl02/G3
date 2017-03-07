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
    
    subroutine read_soil_para()
      use global_para_mod,only:i4,r8
      use soil_para_mod,only:ksat1
      implicit none
      integer(i4)::i,l1,l2
      
      call strlen(soil_para_file,l1,l2)
	do i = 1, 6
	  if(i.eq.1) infile = 'thsgm1k.asc'
	  if(i.eq.2) infile = 'thr1k.asc'
	  if(i.eq.3) infile = 'k01k.asc'
	  if(i.eq.4) infile = 'k11k.asc'
c	  if(i.eq.5) infile = 'ksg_heihe.asc'
	  if(i.eq.5) infile = 'alpha1k.asc'
        if(i.eq.6) infile = 'n1k.asc'
  	  call strlen(infile,ia1,ia2)
	open(11,file=soil_dir(ib1:ib2)//infile(ia1:ia2),status='old')
       

	  
	  read(11,'(A, i16)') ncols,nnc
 
	  read(11,'(A, i16)') nrows,nnr
	  read(11,'(A, f16.3)') xllcorner,x0
	  read(11,'(A, f16.3)') yllcorner,y0
	  read(11,'(A, f16.0)') cellsize,gridsize
	  read(11,'(A, f8.0)') nodata,znodata

	  do ir=1,nnr
	    if(i.eq.1) read(11,*) (wsat(ir,ic),     ic=1,nnc)
	    if(i.eq.2) read(11,*) (wrsd(ir,ic),     ic=1,nnc)
	    if(i.eq.3) read(11,*) (ksat1(ir,ic),    ic=1,nnc)
	    if(i.eq.4) read(11,*) (ksat2(ir,ic),    ic=1,nnc)
c		if(i.eq.5) read(11,*) (kg(ir,ic),       ic=1,nnc)
		if(i.eq.5) read(11,*) (alpha(ir,ic),    ic=1,nnc)
		if(i.eq.6) read(11,*) (watern(ir,ic),   ic=1,nnc)
	  end do
	  close(11)
	end do ! added by gaobing 20131125


       do ir = 1,nnr
	     do ic = 1,nnc

	     if(area(ir,ic) > 0 .and. ksat1(ir,ic) < 0) then
	 print *, 'error in ksat1' ,ir,ic
	     else
c      ksat1(ir,ic) = ksat1(ir,ic)*0.001/3600.0         !mm/hr--->m/s
      ksat1(ir,ic) = ksat1(ir,ic)*0.001/3600.0/2.4         !cm/day--->mm/hr--->m/s, mqh
	     endif
c	  if(area(ir,ic) > 0 .and. ksat2(ir,ic) < 0) then
c	 print *, 'error in ksat2' ,ir,ic
c	     else
c	ksat2(ir,ic) = ksat2(ir,ic)*0.001/3600.0         !mm/hr--->m/s 
c	     endif
c        if(area(ir,ic) > 0 .and. kg(ir,ic) < 0) then
c	 print *, 'error in kg' ,ir,ic
c	     else
c	kg(ir,ic)    = kg(ir,ic)*0.001/3600.0		     !mm/hr--->m/s
c	  endif

	 if(area(ir,ic) > 0 .and. wsat(ir,ic) < 0) then
	 print *, 'error in wsat' ,ir,ic
	 endif

	if(area(ir,ic) > 0 .and. wrsd(ir,ic) < 0) then
	 print *, 'error in wrsd' ,ir,ic
	 endif

	if(area(ir,ic) > 0 .and. alpha(ir,ic) < 0) then
	 print *, 'error in alpha' ,ir,ic
	 endif
	if(area(ir,ic) > 0 .and. watern(ir,ic) < 0) then
	 print *, 'error in watern' ,ir,ic
	 endif
	     enddo
	 enddo

  
         do ir = 1,nnr
	      do ic = 1,nnc
	     if(area(ir,ic) > 0)  then
	          GWcs(ir,ic)=0.15       
!	          ksat1(ir,ic) = 0.9*ksat1(ir,ic)    
	       ksat2(ir,ic) = 0.15*ksat1(ir,ic)  !from 0.1 to 0.15 ,mqh
!	          ksat2(ir,ic) = 0.09*ksat1(ir,ic)  !from 0.1 to 0.15 ,mqh
	           kg(ir,ic)    = 0.03*ksat1(ir,ic)

	      else
	          GWcs(ir,ic)=-9999           
	          ksat2(ir,ic) = -9999 
	           kg(ir,ic)    = -9999
	      end if
	      enddo
	   enddo
c	print *,"ksat=",ksat1(i),ksat2(i),kg(i)

	suction = -1.02  ! suction, meter water-head   !*1.2,added by mqh
       suction2 = -1.02 
        do ir = 1,nnr
           do ic = 1,nnc
	     if(area(ir,ic) > 0 )  then
c	if(wsat(ir,ic).lt. 0.48 ) wsat(ir,ic)=0.48   ! test by gaobing
    
      call MoistureFromSuction_V(tmp, wsat(ir,ic), wrsd(ir,ic), 
     :                           watern(ir,ic), alpha(ir,ic), suction2)
      wfld(ir,ic) = tmp
	if(wfld(ir,ic).gt.wsat(ir,ic) .or. wfld(ir,ic).lt.wrsd(ir,ic))
     :        print *,"wfld out range wrsd ~ wsat"
             ksat2(ir,ic)= ksat2(ir,ic)*1.3
	       kg(ir,ic)=kg(ir,ic)*2.0
           endif
        enddo
	  enddo	 
c      if(wfld(ir,ic).lt. 0.60*wsat(ir,ic))wfld(ir,ic) = 0.60*wsat(ir,ic)
c	if(wfld(ir,ic).gt. 0.85*wsat(ir,ic))wfld(ir,ic) = 0.85*wsat(ir,ic)
c	print *, i, isoil, wsat(i), wfld(i)
c      i = i + 1
c	goto 35
c45	close(1)
      print * ,'ok3, finish of soil_water_para', para_dir
c	pause 


c      print *, i-1
    end subroutine read_soil_para

    subroutine read_soil_code(para_dir)
      use global_para_mod,only:i4,nrow,ncol
      use hydro_data_mod,only:area
      use soil_para_mod,only:soiltyp,nsoil,soil,soiltyp
      implicit none
      character*200,intent(in)::para_dir
      integer(i4)::l1,l2,i,itmp,tmp,ir,ic,isoil,k
      call strlen(para_dir,l1,l2)
      open(14,file=para_dir(l1:l2)//'soil_code.txt', status='old')
      read(14,*)
      do i=1,nsoil
        read(14,*) itmp, soiltyp(i), tmp
      end do
      close(14)
      do ir=1,nrow
        do ic=1,ncol
          if(area(ir,ic).ne.-9999.and.soil(ir,ic) .ne. -9999) then
              k=0
              do isoil = 1, nsoil
                if(soil(ir,ic) .eq. soiltyp(isoil)) itmp=isoil
                if(soil(ir,ic) .eq. soiltyp(isoil)) k=k+1
              end do
              soil(ir,ic) = itmp
              if(k .eq. 0) print *,'no this type of soil:', ir,ic,soil(ir,ic)
              if(k.gt.1) print *,'more one code for this type of soil:',k, soil(ir,ic)
            endif
        end do
      end do
    end subroutine

    subroutine read_map_files()
      use global_para_mod,only:nrow,ncol,i4,area_map,ele_map,slp_map,slplen_map,soil_map,ds_map
      use hydro_data_mod,only:area,ele,length,slp,Ds,dg
      use soil_para_mod,only:soil
      implicit none
      integer(i4)::l1,l2
      character*5:: ncols,nrows
      character*9:: xllcorner,yllcorner
      character*8:: cellsize
      character*12:: nodata
      integer(i4):: nnr, nnc,ir,ic
      real(r8):: x0, y0
      real(r8):: gridsize
      real(r8):: znodata
      
      !read the grid_area mapfile
      call strlen(area_map,l1,l2)
      open(14,file=area_map(l1:l2), status='old')
      read(14,'(A, i16)') ncols,nnc
      read(14,'(A, i16)') nrows,nnr
      read(14,'(A, f16.3)') xllcorner,x0
      read(14,'(A, f16.3)') yllcorner,y0
      read(14,'(A, f16.0)') cellsize,gridsize
      read(14,'(A, f8.0)') nodata,znodata
      if(nnc.ne.ncol .or. nnr.ne.nrow) then
         print*, "The extent of file"//area_map(l1:l2)//"does not match the model extent",nnr,nnc,ncol,nrow
         stop
      else
         do ir=1,nrow
           read(14,*) (area(ir,ic),ic=1,ncol)
         end do
      end if
      close(14)
        
      !read the elevation mapfile
      call strlen(ele_map,l1,l2)
      open(14,file=ele_map(l1:l2), status='old')
      read(14,'(A, i16)') ncols,nnc
      read(14,'(A, i16)') nrows,nnr
      read(14,'(A, f16.3)') xllcorner,x0
      read(14,'(A, f16.3)') yllcorner,y0
      read(14,'(A, f16.0)') cellsize,gridsize
      read(14,'(A, f8.0)') nodata,znodata
      if(nnc.ne.ncol .or. nnr.ne.nrow) then
         print*, "The extent of file"//ele_map(l1:l2)//"does not match the model extent",nnr,nnc,ncol,nrow
         stop
      else
         do ir=1,nrow
           read(14,*) (ele(ir,ic),ic=1,ncol)
         end do
      end if
      close(14)
      
      !read the slope len mapfile
      call strlen(slplen_map,l1,l2)
      open(14,file=slplen_map(l1:l2), status='old')
      read(14,'(A, i16)') ncols,nnc
      read(14,'(A, i16)') nrows,nnr
      read(14,'(A, f16.3)') xllcorner,x0
      read(14,'(A, f16.3)') yllcorner,y0
      read(14,'(A, f16.0)') cellsize,gridsize
      read(14,'(A, f8.0)') nodata,znodata
      if(nnc.ne.ncol .or. nnr.ne.nrow) then
         print*, "The extent of file"//slplen_map(l1:l2)//"does not match the model extent",nnr,nnc,ncol,nrow
         stop
      else
         do ir=1,nrow
           read(14,*) (length(ir,ic),ic=1,ncol)
         end do
      end if
      close(14)
      
      !read the slp mapfile
      call strlen(slp_map,l1,l2)
      open(14,file=slp_map(l1:l2), status='old')
      read(14,'(A, i16)') ncols,nnc
      read(14,'(A, i16)') nrows,nnr
      read(14,'(A, f16.3)') xllcorner,x0
      read(14,'(A, f16.3)') yllcorner,y0
      read(14,'(A, f16.0)') cellsize,gridsize
      read(14,'(A, f8.0)') nodata,znodata
      if(nnc.ne.ncol .or. nnr.ne.nrow) then
         print*, "The extent of file"//slp_map(l1:l2)//"does not match the model extent",nnr,nnc,ncol,nrow
         stop
      else
         do ir=1,nrow
           read(14,*) (slp(ir,ic),ic=1,ncol)
         end do
      end if
      close(14)
      
      !read the soil mapfile
      call strlen(soil_map,l1,l2)
      open(14,file=soil_map(l1:l2), status='old')
      read(14,'(A, i16)') ncols,nnc
      read(14,'(A, i16)') nrows,nnr
      read(14,'(A, f16.3)') xllcorner,x0
      read(14,'(A, f16.3)') yllcorner,y0
      read(14,'(A, f16.0)') cellsize,gridsize
      read(14,'(A, f8.0)') nodata,znodata
      if(nnc.ne.ncol .or. nnr.ne.nrow) then
         print*, "The extent of file"//soil_map(l1:l2)//"does not match the model extent",nnr,nnc,ncol,nrow
         stop
      else
         do ir=1,nrow
           read(14,*) (soil(ir,ic),ic=1,ncol)
         end do
      end if
      close(14)
      !read the soil depth mapfile
      call strlen(ds_map,l1,l2)
      open(14,file=ds_map(l1:l2), status='old')
      read(14,'(A, i16)') ncols,nnc
      read(14,'(A, i16)') nrows,nnr
      read(14,'(A, f16.3)') xllcorner,x0
      read(14,'(A, f16.3)') yllcorner,y0
      read(14,'(A, f16.0)') cellsize,gridsize
      read(14,'(A, f8.0)') nodata,znodata
      if(nnc.ne.ncol .or. nnr.ne.nrow) then
         print*, "The extent of file"//ds_map(l1:l2)//"does not match the model extent",nnr,nnc,ncol,nrow
         stop
      else
         do ir=1,nrow
           read(14,*) (Ds(ir,ic),ic=1,ncol)
         end do
      end if
      close(14)
      
      do ir=1,nnr
        do ic=1,nnc
          if (soil(ir,ic).eq.0) then
            print *, 'ir,ic,soil(ir,ic)=', ir,ic,soil(ir,ic)
          endif
          if(area(ir,ic).eq.-9999) then
            ele(ir,ic)=-9999
            length(ir,ic)=-9999
            slp(ir,ic)=-9999
            soil(ir,ic)=-9999
            Ds(ir,ic)=-9999
          endif
          if(soil(ir,ic).eq.-9999 .and. Ds(ir,ic).ne.-9999) then
            print *, 'soil1', soil(ir,ic), Ds(ir,ic), ir, ic
          end if
          if(soil(ir,ic).ne.-9999 .and. Ds(ir,ic).eq.-9999) then
            print *, 'soil2', soil(ir,ic), Ds(ir,ic), ir, ic
          end if
          if(area(ir,ic).ne.-9999) then
            area(ir,ic)=area(ir,ic)*1000.0*1000.0
          endif
          if(slp(ir,ic) .ne. -9999) then
            slp(ir,ic) = slp(ir,ic)/100.0
          endif
        end do
      end do
      
      do ir=1,nnr
        do ic=1,nnc
            Dg(ir,ic)=-9999.0
            if(area(ir,ic).ne.-9999.and.Ds(ir,ic).ne.-9999) then
              Dg(ir,ic)=10.0*Ds(ir,ic)
              if(ele(ir,ic) .ge. 3500.0) Dg(ir,ic)=amin1(2.0, Ds(ir,ic))
              if(Dg(ir,ic) .lt. 0)  print *, "wrong in calculate Dg" 
            endif
        end do
      end do
      
    end subroutine read_map_files


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
      use soil_para_mod,only:soiltyp,soil,wsat,wrsd,wfld,alpha,watern,ksat1,ksat2,anik,ns,nsoil
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


      open(14,file=para_dir(l1:l2)//'soil_code.txt', status='old')
      read(14,*)
      subcount=0
      iostatus=1
      do while(iostatus.ge.0)
        read(14,*,iostat=iostatus) atmp
        subcount=subcount+1
      end do
      close(14)
      subcount=subcount-1
      nsoil=subcount
      ns=nsoil
      

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
