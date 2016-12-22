    module global_para_mod
      public
      save
      integer,parameter::i4=4
      integer,parameter::r8=4
      integer(i4)::nrow             ! number of rows of the whole area
      integer(i4)::ncol             ! number of columns of the whole area
      integer(i4)::year,month2,day,hour,month
      integer(i4)::iflai            ! 1 for using lai 
      
      
    end module global_para_mod
