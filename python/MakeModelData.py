from extract_soil_mod import *
from gisutil import *
from extract_lai import *
from extract_landuse import *
from extract_forcing_Yangkun import *
from tsavg import *



#  Get the latitude and longitude
model_para_dir="/home/zhangyanlin-t7610/model/Horton_SHAWDHM/BBH/input/"
InputDemFileName="/home/zhangyanlin-t7610/model/GBHM/preprocess/Funy/examples/bbh/BBH1km.tif"
wsfile="/home/zhangyanlin-t7610/model/Horton_SHAWDHM/BBH/input/hydropara/ws.asc"
proj4string=""
latfilename=model_para_dir+"hydropara/lat.asc"
lonfilename=model_para_dir+"hydropara/lon.asc"
nc,nr,xllcorner,yllcorner,cellsize,nodata=readaschead(wsfile)
lat,lon=GetGridLatLon2(InputDemFileName,model_para_dir,proj4string,latfilename,lonfilename)

soil=False
if(soil==True):
      #  Extract Soil Parameters
      sourceSoilDir="/disk2/data/soil/dai/nc4/"
      OutputSoilFile="/home/zhangyanlin-t7610/model/Horton_SHAWDHM/BBH/input/soilpara.nc"
      SoilNodeDepths=[0.025,0.1,0.2,0.3,0.5,0.7,1.0,1.5,2.4,2.5]
      print "extracting soil parameter"
      ExtractSoilDataFromDaiYongJiu(sourceSoilDir,OutputSoilFile,SoilNodeDepths,lat,lon)

landuse=False
if (landuse==True):
      #  Extract land use cover
      print "extracting landuse cover"
      SourceLanduse="/disk2/data/landuse/landuse.nc"
      OutputLanduse="/home/zhangyanlin-t7610/model/Horton_SHAWDHM/BBH/input/landuse.nc"
      landuse=Extract_Landuse(SourceLanduse,lat,lon)
      writeascfile("/home/zhangyanlin-t7610/model/Horton_SHAWDHM/BBH/input/landuse.asc",nr,nc,xllcorner,yllcorner,cellsize,nodata,landuse,"int")
      dims=np.shape(landuse)
      landusefile = Dataset(OutputLanduse,"w",format="NETCDF4")
      nlat   = landusefile.createDimension("NY", dims[0])
      nlon   = landusefile.createDimension("NX", dims[1])
      landuses= landusefile.createVariable("landuse","i1",("NY","NX",))
      landuses.missing_value="-9999"
      landuses.units=""
      landuses.long_name="USGS Land cover"
      landuses[:,:]=landuse[:,:]
      landusefile.close()

lai=False
if (lai==True):
      #  Extract lai
      print "extracting lai"
      StartYear=2010
      EndYear=2015
      for y in range(StartYear,EndYear): 
            print y
            SourceLAIname="/disk2/data/LAI/nc4/global_30s_"+str(y)+".nc"
            LAI=Extract_LAI(SourceLAIname,lat,lon)
            dims=np.shape(LAI)
            laifile = Dataset("/home/zhangyanlin-t7610/model/Horton_SHAWDHM/BBH/input/lai/"+str(y)+".nc","w",format="NETCDF4")
            nlat   = laifile.createDimension("NY", dims[1])
            nlon   = laifile.createDimension("NX", dims[2])
            ntime  = laifile.createDimension("time", dims[0])
            Lais= laifile.createVariable("LAI","f4",("time","NY","NX",))
            Lais.missing_value="-9999"
            Lais.units="-"
            Lais.long_name="Leaf area index"
            Lais[:,:,:]=LAI[:,:,:]
            laifile.close()



forcing=False
if(forcing==True):
      ForingOutputDir="/home/zhangyanlin-t7610/model/Horton_SHAWDHM/BBH/input/forcing/Yangkun/"
      SrcDir="/disk2/data/Forcing/Yangkun/Data_forcing_03hr_010deg_nc4/"
      interpolationType=2
      timezone=8
      compress=False
      ExtractHourlyDataFromYangKunDataset(SrcDir,ForingOutputDir,StartYear,EndYear,lat,lon,interpolationType,timezone,compress)

SourceDir="/home/zhangyanlin-t7610/model/Horton_SHAWDHM/BBH/input/forcing/Yangkun/"
SourceDir2="/home/zhangyanlin-t7610/model/Horton_SHAWDHM/BBH/output/"
tsavg=True
StartYear=2010
EndYear=2014
if(tsavg==True):
      #tsavg=ExtTsavg(SourceDir,StartYear,EndYear,nr,nc)
      tsavg=ExtTsavgFromCalculatedSoil(SourceDir2,StartYear,EndYear,nr,nc)
      writeascfile("/home/zhangyanlin-t7610/model/Horton_SHAWDHM/BBH/input/tsavg5.asc",nr,nc,xllcorner,yllcorner,cellsize,nodata,tsavg,"real")

