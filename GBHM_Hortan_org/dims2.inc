c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######                     DIMS2.INC                        ######
c     ######              Developed by Dr. Dawen YANG             ######
c     ######          Department of Civil Engineering             ######
c     ######                University of Tokyo                   ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
c
c
c#######################################################################
c
c     PURPOSE: Using for GBHM2
c
c     Include file 'dims2.inc', pass global paramter 
c
C#######################################################################
c      
c     AUTHOR: Dr. D. Yang
c                                                                               
c     MODIFICATION HISTORY:
c
c#######################################################################
c
c#######################################################################
c
c     Define dimension size
c
c#######################################################################
c
      integer nrow          ! number of rows of the whole area
      integer ncol          ! number of columns of the whole area
c
      integer nc            ! number of sub-catchments 
      integer nx            ! number of flow interval in a sub-catchment
      integer np            ! number of grids in a flow-interval
c
      integer nv            ! number of landuse types
      integer ns            ! number of soil types
      integer nlayer        ! number of UZ layers
c
	  integer startyear		! simulation start
	  integer endyear       ! simulation end
	  integer inicon        ! the way to specify initial condition
	                        ! 1 - by data file in directory ../simulation
							! 0 - arbitrary given in the program
	  integer MNgauge		! Number of Gauge          old     
	  integer Maxgauge      ! Max_Number of Gauge      old +new
	  integer NNstation     ! rainstation select           added by xujijun
	                        ! 1 for old 
	                        ! 2 for new 
	                        ! 3 for old+new 
	  integer iflai         ! 1 for using lai 
          
c
      parameter (nrow = 106, ncol = 78)
      parameter (nc = 153, nx = 300, np = 500)
	  parameter (nv = 9,  ns = 3,  nlayer = 10)
	  parameter (startyear =2000, endyear =2000)
	  parameter (inicon = 0)
	  parameter (MNgauge = 10)
	  parameter (NNstation = 3)
  	  parameter (Maxgauge = 19)
	  parameter (iflai = 0)      






