            if(iflai .eq. 0 ) then
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
