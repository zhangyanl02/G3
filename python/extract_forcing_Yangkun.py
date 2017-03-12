from netCDF4 import Dataset
import numpy as np
from gisutil import *
import os

def SpecificHumidity2RelativeHumidity(SH,T,P):
      RH=0.263*P*SH/(np.exp(17.67*(T-273.15)/(T-29.65)))
      return RH



def ExtVarFromYangKun(ncfilename,varname,lat,lon,InterpolationType):
      ncfile = Dataset(ncfilename.strip(), "r")
      nlat=len(ncfile.dimensions["lat"])
      nlon=len(ncfile.dimensions["lon"])
      ntime=len(ncfile.dimensions["time"])
      lat0=ncfile.variables["lat"][0]
      latn=ncfile.variables["lat"][nlat-1]
      lon0=ncfile.variables["lon"][0]
      lonn=ncfile.variables["lon"][nlon-1]
      times=ncfile.variables["time"][:]
      InLatCellsize=(latn-lat0)/(nlat-1)
      InLonCellsize=(lonn-lon0)/(nlon-1)
      dims=np.shape(lon)
      extVar=np.zeros((ntime,dims[0],dims[1]))
      row=np.zeros((dims[0],dims[1]),dtype=np.int)
      col=np.zeros((dims[0],dims[1]),dtype=np.int)
      for i in range(dims[0]):
            for j in range(dims[1]):
                  row[i,j]=int((lat[i,j]-lat0+0.5*InLonCellsize)/InLatCellsize)
                  col[i,j]=int((lon[i,j]-lon0+0.5*InLatCellsize)/InLonCellsize)
      minrow=np.max([np.min(row)-1,0])
      maxrow=np.min([np.max(row)+1,nlat-1])
      mincol=np.max([np.min(col)-1,0])
      maxcol=np.min([np.max(col)+1,nlon-1])      

      tempdata=ncfile.variables[varname][:,minrow:maxrow+1,mincol:maxcol+1]
      latul=ncfile.variables["lat"][minrow]
      lonul=ncfile.variables["lon"][mincol]
      ncfile.close()

      for r in range(dims[0]):
            for c in range(dims[1]):            
                  # No interpolation
                  if(InterpolationType==0):
                        extVar[:,r,c]=tempdata[:,row[r,c]-minrow,col[r,c]-mincol]                  
                  else:
                        IX=lon[r,c]
                        IY=lat[r,c]
                        InCOL=np.int((IX-lonul+0.5*InLonCellsize)/InLonCellsize)
                        InRow=np.int((IY-latul+0.5*InLatCellsize)/InLatCellsize)
                        DeltaX=lonul+InCOL*InLonCellsize-IX
                        DeltaY=latul+InRow*InLatCellsize-IY
                        if(DeltaX<0.0):
                              INCOL2=InCOL+1
                        else:
                              INCOL2=InCOL-1
                        if(DeltaY<0.0):
                              INROW2=InRow+1
                        else:
                              INROW2=InRow-1
                        
                        if(INROW2<0):
                              INROW2=InRow+1
                        if(INCOL2<0):
                              INCOL2=InCOL+1
                        if(INROW2>=maxrow+1-minrow):
                              INROW2=InRow-1
                        if(INCOL2>=maxcol+1-mincol):
                              INCOL2=InCOL-1                        

                        X1=lonul+InCOL*InLonCellsize
                        X2=lonul+INCOL2*InLonCellsize
                        Y1=latul+InRow*InLatCellsize
                        Y2=latul+INROW2*InLatCellsize

                        # Inverse distance interpolation
                        if(InterpolationType==1 and IX-X1!=0.0 and IY-Y1!=0 and IY-Y2!=0 and IX-X2!=0):                  
                              R12=1.0/((IX-X1)**2.0+(IY-Y1)**2.0)
                              R22=1.0/((IX-X1)**2.0+(IY-Y2)**2.0)
                              R32=1.0/((IX-X2)**2.0+(IY-Y1)**2.0)
                              R42=1.0/((IX-X2)**2.0+(IY-Y2)**2.0)
                              Rsigma=R12+R22+R32+R42
                              extVar[:,r,c]=(R12*tempdata[:,InRow,InCOL]+R22*tempdata[:,INROW2,InCOL]+R32*tempdata[:,InRow,INCOL2]+R42*tempdata[:,INROW2,INCOL2])/Rsigma
                        else:
                        # Bilinear interpolation
                              R1=((X2-IX)*tempdata[:,InRow,InCOL]+(IX-X1)*tempdata[:,InRow,INCOL2])/(X2-X1)
                              R2=((X2-IX)*tempdata[:,INROW2,InCOL]+(IX-X1)*tempdata[:,INROW2,INCOL2])/(X2-X1)
                              extVar[:,r,c]=((Y2-IY)*R1+(IY-Y1)*R2)/(Y2-Y1)                  
      return times,extVar


def ConvertThreeHourlyVar2Hourly(varsin):
      ntime=np.shape(varsin)[0]
      nlat=np.shape(varsin)[1]
      nlon=np.shape(varsin)[2]
      HourlyVar=np.zeros((3*ntime,nlat,nlon))
      for t in range(ntime-1):
            HourlyVar[3*t,:,:]=varsin[t,:,:]
            HourlyVar[3*t+1,:,:]=varsin[t,:,:]*2.0/3.0+varsin[t+1,:,:]*1.0/3.0            
            HourlyVar[3*t+2,:,:]=varsin[t,:,:]*1.0/3.0+varsin[t+1,:,:]*2.0/3.0
      for i in range(3):
            HourlyVar[3*(ntime-1)+i,:,:]=varsin[ntime-1,:,:]
      return HourlyVar
def GetHourlyTimes(ThreeHourlyTimes):
      ntime=np.shape(ThreeHourlyTimes)[0]
      timeHourly=np.zeros((3*ntime),dtype=np.int)
      for t in range(ntime-1):
            timeHourly[3*t]=ThreeHourlyTimes[t]
            timeHourly[3*t+1]=ThreeHourlyTimes[t]+1                  
            timeHourly[3*t+2]=ThreeHourlyTimes[t]+2                 
      for i in range(3):
            timeHourly[3*(ntime-1)+i]=ThreeHourlyTimes[ntime-1]+i
      return timeHourly   




def GetTheLast8HoursDataInPreviousYear(SrcDir,VariableName,year,lat,lon,interpolationType,timezone):
      nlat=np.shape(lat)[0]
      nlon=np.shape(lat)[1]
      ncfilename=SrcDir+VariableName+"/"+VariableName+"_ITPCAS-CMFD_V0106_B-01_03hr_010deg_"+str(year*100+12)+".nc"
      timesout,VarData=ExtVarFromYangKun(ncfilename,VariableName,lat,lon,interpolationType)
      VarHourly=ConvertThreeHourlyVar2Hourly(VarData)            
      timeHourly=GetHourlyTimes(timesout)
      timeHourly=timeHourly+timezone
      VarLastMonth=np.zeros((timezone,nlat,nlon))
      timeLastMonth=np.zeros((timezone),dtype=np.int)            
      VarLastMonth[:,:,:]=VarHourly[-timezone:,:,:]
      timeLastMonth[:]=timeHourly[-timezone:]
      return VarLastMonth,timeLastMonth



def GetHourlyVariable(SrcDir,VariableName,datestr,lat,lon,interpolationType,timezone):
            ncfilename=SrcDir+VariableName+"/"+VariableName+"_ITPCAS-CMFD_V0106_B-01_03hr_010deg_"+datestr+".nc"
            timesout,VarData=ExtVarFromYangKun(ncfilename,VariableName,lat,lon,interpolationType)
            VarHourly=ConvertThreeHourlyVar2Hourly(VarData)            
            timeHourly=GetHourlyTimes(timesout)
            timeHourly=timeHourly+timezone
            return VarHourly,timeHourly




def ExtractHourlyDataFromYangKunDataset(SrcDir,ForingOutputDir,startYear,EndYear,lat,lon,interpolationType,timezone,compress):
      nlat=np.shape(lat)[0]
      nlon=np.shape(lat)[1]
      precsLastMonth=np.zeros((timezone,nlat,nlon))
      tempsLastMonth=np.zeros((timezone,nlat,nlon))
      lradsLastMonth=np.zeros((timezone,nlat,nlon))            
      sradsLastMonth=np.zeros((timezone,nlat,nlon))            
      pressLastMonth=np.zeros((timezone,nlat,nlon))            
      shumsLastMonth=np.zeros((timezone,nlat,nlon))            
      windsLastMonth=np.zeros((timezone,nlat,nlon))            
      rhsLastMonth=np.zeros((timezone,nlat,nlon))            
      timeLastMonth=np.zeros((timezone),dtype=np.int)
      precsLastMonth[:,:,:]=-32767
      tempsLastMonth[:,:,:]=-32767
      lradsLastMonth[:,:,:]=-32767      
      sradsLastMonth[:,:,:]=-32767         
      pressLastMonth[:,:,:]=-32767  
      shumsLastMonth[:,:,:]=-32767
      windsLastMonth[:,:,:]=-32767
      rhsLastMonth[:,:,:]=-32767
      timeLastMonth[:]=-32767

      if(startYear>1979):
            precsLastMonth,timeLastMonth=GetTheLast8HoursDataInPreviousYear(SrcDir,"prec",startYear-1,lat,lon,interpolationType,timezone)
            tempsLastMonth,timeLastMonth=GetTheLast8HoursDataInPreviousYear(SrcDir,"temp",startYear-1,lat,lon,interpolationType,timezone)
            lradsLastMonth,timeLastMonth=GetTheLast8HoursDataInPreviousYear(SrcDir,"lrad",startYear-1,lat,lon,interpolationType,timezone)
            sradsLastMonth,timeLastMonth=GetTheLast8HoursDataInPreviousYear(SrcDir,"srad",startYear-1,lat,lon,interpolationType,timezone)    
            pressLastMonth,timeLastMonth=GetTheLast8HoursDataInPreviousYear(SrcDir,"pres",startYear-1,lat,lon,interpolationType,timezone)      
            shumsLastMonth,timeLastMonth=GetTheLast8HoursDataInPreviousYear(SrcDir,"shum",startYear-1,lat,lon,interpolationType,timezone)     
            windsLastMonth,timeLastMonth=GetTheLast8HoursDataInPreviousYear(SrcDir,"wind",startYear-1,lat,lon,interpolationType,timezone)     
            rhsLastMonth=SpecificHumidity2RelativeHumidity(shumsLastMonth,tempsLastMonth,pressLastMonth)
      
      for y in range(startYear,EndYear):
            jd=1
            Dir=ForingOutputDir+str(y)
            if not os.path.exists(Dir):
                  os.makedirs(Dir)
            for m in range(1,13):
                  datestr=str(y*100+m)
                  print datestr
                  precsHourly,timeHourly=GetHourlyVariable(SrcDir,"prec",datestr,lat,lon,interpolationType,timezone)
                  tempsHourly,timeHourly=GetHourlyVariable(SrcDir,"temp",datestr,lat,lon,interpolationType,timezone)
                  lradsHourly,timeHourly=GetHourlyVariable(SrcDir,"lrad",datestr,lat,lon,interpolationType,timezone)
                  sradsHourly,timeHourly=GetHourlyVariable(SrcDir,"srad",datestr,lat,lon,interpolationType,timezone)
                  pressHourly,timeHourly=GetHourlyVariable(SrcDir,"pres",datestr,lat,lon,interpolationType,timezone)
                  shumsHourly,timeHourly=GetHourlyVariable(SrcDir,"shum",datestr,lat,lon,interpolationType,timezone)
                  windsHourly,timeHourly=GetHourlyVariable(SrcDir,"wind",datestr,lat,lon,interpolationType,timezone)
                  rhsHourly=SpecificHumidity2RelativeHumidity(shumsHourly,tempsHourly,pressHourly)

                  ndays=np.int(np.shape(timeHourly)[0]/24)
                  for d in range(ndays):
                        outputncfilename=Dir+"/"+str(y)+"-"+str(y*1000+jd)[-3:]+".nc"
                        nc4file = Dataset(outputncfilename, "w", format="NETCDF4")
                        NY     = nc4file.createDimension("NY", nlat)
                        NX     = nc4file.createDimension("NX", nlon)
                        time    = nc4file.createDimension("time",24)
      
                        lats= nc4file.createVariable("lats","f4",("NY","NX",))
                        lats.units="degrees_north"
                        lats.long_name="Latitude"
                        lats[:,:]=lat[:,:]
            
                        lons= nc4file.createVariable("lons","f4",("NY","NX",))
                        lons.units="degrees_east"
                        lons.long_name="Longitude"
                        lons[:,:]=lon[:,:]

                        times= nc4file.createVariable("time","i4",("time",))
                        times.units="hours since 1900-01-01 00:00:0.0"
                        times.long_name="Time"

                        prec= nc4file.createVariable("prec","f4",("time","NY","NX",),fill_value=-32767,zlib=compress)
                        #prec.scale_factor=0.0025
                        #prec.add_offset=50.0
                        prec.missing_value=-32767
                        prec.units="mm hr-1"
                        prec.long_name="Precipitation rate"

                        temp= nc4file.createVariable("temp","f4",("time","NY","NX",),fill_value=-32767,zlib=compress)
                        #temp.scale_factor=0.01
                        #temp.add_offset=0.0
                        temp.missing_value=-32767
                        temp.units="K"
                        temp.long_name="Near surface air temperature"

                        lrad= nc4file.createVariable("lrad","f4",("time","NY","NX",),fill_value=-32767,zlib=compress)
                        #lrad.scale_factor=0.25
                        #lrad.add_offset=685.0
                        lrad.missing_value=-32767
                        lrad.units="W m-2"
                        lrad.long_name="Surface downward longwave radiation"
                  
                        srad= nc4file.createVariable("srad","f4",("time","NY","NX",),fill_value=-32767,zlib=compress)
                        #srad.scale_factor=0.25
                        #srad.add_offset=685.0
                        srad.missing_value=-32767
                        srad.units="W m-2"
                        srad.long_name="Surface downward shortwave radiation"
                  
                        pres= nc4file.createVariable("pres","f4",("time","NY","NX",),fill_value=-32767,zlib=compress)
                        #pres.scale_factor=2.0
                        #pres.add_offset=63500.
                        pres.missing_value=-32767
                        pres.units="Pa"
                        pres.long_name="Near surface air pressure"                  

                        shum= nc4file.createVariable("shum","f4",("time","NY","NX",),fill_value=-32767,zlib=compress)
                        #shum.scale_factor=0.000001
                        #shum.add_offset=0.025
                        shum.missing_value=-32767
                        shum.units="kg kg-1"
                        shum.long_name="Near surface air specific humidity"
                        
                        wind= nc4file.createVariable("wind","f4",("time","NY","NX",),fill_value=-32767,zlib=compress)
                        #wind.scale_factor=0.002
                        #wind.add_offset=60.0
                        wind.missing_value=-32767
                        wind.units="m s-1"
                        wind.long_name="Near surface wind speed"

                        RH= nc4file.createVariable("rhum","f4",("time","NY","NX",),fill_value=-32767,zlib=compress)
                        #RH.scale_factor=0.01            
                        #RH.missing_value=-32767
                        RH.units="%"
                        RH.long_name="Near surface Relative Humidity"

                        if (d==0):                  
                              times[0:timezone]=timeLastMonth[:]
                              prec[0:timezone,:,:]=precsLastMonth[:,:,:]
                              temp[0:timezone,:,:]=tempsLastMonth[:,:,:]
                              lrad[0:timezone,:,:]=lradsLastMonth[:,:,:]
                              srad[0:timezone,:,:]=sradsLastMonth[:,:,:]
                              pres[0:timezone,:,:]=pressLastMonth[:,:,:]
                              shum[0:timezone,:,:]=shumsLastMonth[:,:,:]
                              wind[0:timezone,:,:]=windsLastMonth[:,:,:]
                              RH[0:timezone,:,:]=rhsLastMonth[:,:,:]                        
                              times[timezone:24]=timeHourly[d*24:(d+1)*24-timezone]
                              prec[timezone:24,:,:]=precsHourly[d*24:(d+1)*24-timezone,:,:]
                              temp[timezone:24,:,:]=tempsHourly[d*24:(d+1)*24-timezone,:,:]
                              lrad[timezone:24,:,:]=lradsHourly[d*24:(d+1)*24-timezone,:,:]
                              srad[timezone:24,:,:]=sradsHourly[d*24:(d+1)*24-timezone,:,:]
                              pres[timezone:24,:,:]=pressHourly[d*24:(d+1)*24-timezone,:,:]
                              shum[timezone:24,:,:]=shumsHourly[d*24:(d+1)*24-timezone,:,:]
                              wind[timezone:24,:,:]=windsHourly[d*24:(d+1)*24-timezone,:,:]
                              RH[timezone:24,:,:]=rhsHourly[d*24:(d+1)*24-timezone,:,:]                        
                        else:
                              times[:]=timeHourly[d*24-timezone:(d+1)*24-timezone]
                              prec[:,:,:]=precsHourly[d*24-timezone:(d+1)*24-timezone,:,:]
                              temp[:,:,:]=tempsHourly[d*24-timezone:(d+1)*24-timezone,:,:]
                              lrad[:,:,:]=lradsHourly[d*24-timezone:(d+1)*24-timezone,:,:]
                              srad[:,:,:]=sradsHourly[d*24-timezone:(d+1)*24-timezone,:,:]
                              pres[:,:,:]=pressHourly[d*24-timezone:(d+1)*24-timezone,:,:]
                              shum[:,:,:]=shumsHourly[d*24-timezone:(d+1)*24-timezone,:,:]
                              wind[:,:,:]=windsHourly[d*24-timezone:(d+1)*24-timezone,:,:]
                              RH[:,:,:]=rhsHourly[d*24-timezone:(d+1)*24-timezone,:,:]                        
                        if(d==ndays-1 and timezone>0):
                              precsLastMonth[:,:,:]=precsHourly[-timezone:,:,:]
                              tempsLastMonth[:,:,:]=tempsHourly[-timezone:,:,:]
                              lradsLastMonth[:,:,:]=lradsHourly[-timezone:,:,:]
                              sradsLastMonth[:,:,:]=sradsHourly[-timezone:,:,:]
                              pressLastMonth[:,:,:]=pressHourly[-timezone:,:,:]
                              shumsLastMonth[:,:,:]=shumsHourly[-timezone:,:,:]
                              windsLastMonth[:,:,:]=windsHourly[-timezone:,:,:]
                              rhsLastMonth[:,:,:]=rhsHourly[-timezone:,:,:]
                              timeLastMonth[:]=timeHourly[-timezone:]                             
                        nc4file.close()
                        jd=jd+1






