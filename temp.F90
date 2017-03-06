    subroutine initial_subcatchment(para_dir)
      use global_para_mod,only:i4,r8
      use hydro_data_mod,only:nsub
      implicit none
      character*200,intent(in)::para_dir
      integer::i,j,k,fileunit
      character(200)::atmp
      integer(i4)::l1,l2,subcount,iostatus
      character*4::ch4 
      character*200::subbasinfile
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
    end subroutine
