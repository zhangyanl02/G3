from netCDF4 import Dataset
from gisutil import *
import numpy as np

def DetermineLayers(Depths,depth_in_org_dataset):
      depth_in_org_dataset[:]=depth_in_org_dataset[:]/100.0
      outputlayer=np.shape(Depths)[0]
      sourcelayer=np.shape(depth_in_org_dataset)[0]
      layer=np.zeros((outputlayer,2*sourcelayer+1),dtype=np.int)
      upper=np.zeros((outputlayer,1),dtype=np.float)
      lower=np.zeros((outputlayer,1),dtype=np.float)
      upper1=np.zeros((sourcelayer,1),dtype=np.float)
      lower1=np.zeros((sourcelayer,1),dtype=np.float)
      dz=np.zeros((outputlayer,sourcelayer),dtype=np.float)
      upper[0]=0.0
      for l in range(1,outputlayer):
            upper[l]=0.5*(Depths[l-1]+Depths[l])
            lower[l-1]=upper[l]
      lower[outputlayer-1]=2*Depths[outputlayer-1]-upper[outputlayer-1]
      upper1[0]=0.0
      lower1[0]=depth_in_org_dataset[0]
      for l in range(1,sourcelayer):
            upper1[l]=depth_in_org_dataset[l-1]
            lower1[l]=depth_in_org_dataset[l]
      for l in range(outputlayer):
            for d in range(sourcelayer):
                  if(upper1[d]>=upper[l] and upper1[d]<=lower[l]):
                        dz[l,d]=np.min([lower1[d],lower[l]])-upper1[d]
                  if(lower1[d]>=upper[l] and lower1[d]<=lower[l]):
                        dz[l,d]=lower1[d]-np.max([upper1[d],upper[l]])
                  if(upper1[d]<=upper[l] and lower1[d]>=lower[l]):
                        dz[l,d]=lower[l]-upper[l]
                  if(upper1[sourcelayer-2]<=upper[l] and lower1[sourcelayer-1]>=upper[l]):
                        if(lower[l]>=lower1[sourcelayer-1]):
                              dz[l,sourcelayer-1]=lower[l]-upper[l]
                  if(upper1[sourcelayer-1]<=upper[l] and upper1[sourcelayer-1]<lower[l]):
                        if(lower[l]>=lower1[sourcelayer-1]):
                              dz[l,sourcelayer-1]=lower[l]-upper[l]                        
            dz[l,:]=dz[l,:]/(lower[l]-upper[l])
            if(np.sum(dz[l,:])==0):
                  dz[l,sourcelayer-1]=1.0
            if(np.sum(dz[l,:])!=1):
                  print "error in determine layer in layer:",l
      return dz
      
def ExtractBDFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon):
      BDfile=sourceDir+"BD.nc"
      BDgroup = Dataset(BDfile, "r")
      nlat=len(BDgroup.dimensions["lat"])
      nlon=len(BDgroup.dimensions["lon"])
      ndepth=len(BDgroup.dimensions["depth"])
      lat0=BDgroup.variables["lat"][0]
      latn=BDgroup.variables["lat"][nlat-1]
      lon0=BDgroup.variables["lon"][0]
      lonn=BDgroup.variables["lon"][nlon-1]
      depth_in_org_dataset=BDgroup.variables["depth"][:]
      layerfactor=DetermineLayers(Depths,depth_in_org_dataset)
      latcellsize=(latn-lat0)/(nlat-1)
      loncellsize=(lonn-lon0)/(nlon-1)
      dims=np.shape(lon)
      extBD=np.zeros((ndepth,dims[0],dims[1]))
      tmpBD=np.zeros((Nlayer,dims[0],dims[1]))
      row=np.zeros((dims[0],dims[1]),dtype=np.int)
      col=np.zeros((dims[0],dims[1]),dtype=np.int)
      for i in range(dims[0]):
            for j in range(dims[1]):
                  row[i,j]=int((lat[i,j]-lat0)/latcellsize)
                  col[i,j]=int((lon[i,j]-lon0)/loncellsize)      
      minrow=np.min(row)
      maxrow=np.max(row)
      mincol=np.min(col)
      maxcol=np.max(col)
      tempdata=BDgroup.variables["BD"][:,minrow:maxrow+1,mincol:maxcol+1]
      BDgroup.close()
      for i in range(dims[0]):
            for j in range(dims[1]):
                  for d in range(ndepth):
                        extBD[d,i,j]=tempdata[d,row[i,j]-minrow,col[i,j]-mincol]
                        if(np.isnan(extBD[d,i,j]) or extBD[d,i,j]==-999):
                              extBD[d,i,j]=-9999
                        if(extBD[d,i,j]!=-9999):
                              extBD[d,i,j]=extBD[d,i,j]*1000.0
                        else:
                              extBD[d,i,j]=extBD[d-1,i,j]
                  for nl in range(Nlayer):
                        tmpBD[nl,i,j]=np.inner(layerfactor[nl,:],extBD[:,i,j])
      return tmpBD
      
def ExtractSAFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon):
      SAfile=sourceDir+"SA.nc"
      SAgroup = Dataset(SAfile, "r")
      nlat=len(SAgroup.dimensions["lat"])
      nlon=len(SAgroup.dimensions["lon"])
      ndepth=len(SAgroup.dimensions["depth"])
      lat0=SAgroup.variables["lat"][0]
      latn=SAgroup.variables["lat"][nlat-1]
      lon0=SAgroup.variables["lon"][0]
      lonn=SAgroup.variables["lon"][nlon-1]
      depth_in_org_dataset=SAgroup.variables["depth"][:]
      layerfactor=DetermineLayers(Depths,depth_in_org_dataset)
      latcellsize=(latn-lat0)/(nlat-1)
      loncellsize=(lonn-lon0)/(nlon-1)
      dims=np.shape(lon)
      extSA=np.zeros((ndepth,dims[0],dims[1]))
      tmpSA=np.zeros((Nlayer,dims[0],dims[1]))
      row=np.zeros((dims[0],dims[1]),dtype=np.int)
      col=np.zeros((dims[0],dims[1]),dtype=np.int)
      for i in range(dims[0]):
            for j in range(dims[1]):
                  row[i,j]=int((lat[i,j]-lat0)/latcellsize)
                  col[i,j]=int((lon[i,j]-lon0)/loncellsize)      
      minrow=np.min(row)
      maxrow=np.max(row)
      mincol=np.min(col)
      maxcol=np.max(col)
      tempdata=SAgroup.variables["SA"][:,minrow:maxrow+1,mincol:maxcol+1]
      SAgroup.close()
      for i in range(dims[0]):
            for j in range(dims[1]):
                  for d in range(ndepth):
                        extSA[d,i,j]=tempdata[d,row[i,j]-minrow,col[i,j]-mincol]
                        if(np.isnan(extSA[d,i,j]) or extSA[d,i,j]==-999):
                              extSA[d,i,j]=-9999
                        if(extSA[d,i,j]==-9999):
                              extSA[d,i,j]=extSA[d-1,i,j]
                  for nl in range(Nlayer):
                        tmpSA[nl,i,j]=np.inner(layerfactor[nl,:],extSA[:,i,j])
      return tmpSA

def ExtractCLFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon):
      CLfile=sourceDir+"CL.nc"
      CLgroup = Dataset(CLfile, "r")
      nlat=len(CLgroup.dimensions["lat"])
      nlon=len(CLgroup.dimensions["lon"])
      ndepth=len(CLgroup.dimensions["depth"])
      lat0=CLgroup.variables["lat"][0]
      latn=CLgroup.variables["lat"][nlat-1]
      lon0=CLgroup.variables["lon"][0]
      lonn=CLgroup.variables["lon"][nlon-1]
      depth_in_org_dataset=CLgroup.variables["depth"][:]
      layerfactor=DetermineLayers(Depths,depth_in_org_dataset)
      latcellsize=(latn-lat0)/(nlat-1)
      loncellsize=(lonn-lon0)/(nlon-1)
      dims=np.shape(lon)
      extCL=np.zeros((ndepth,dims[0],dims[1]))
      tmpCL=np.zeros((Nlayer,dims[0],dims[1]))
      row=np.zeros((dims[0],dims[1]),dtype=np.int)
      col=np.zeros((dims[0],dims[1]),dtype=np.int)
      for i in range(dims[0]):
            for j in range(dims[1]):
                  row[i,j]=int((lat[i,j]-lat0)/latcellsize)
                  col[i,j]=int((lon[i,j]-lon0)/loncellsize)      
      minrow=np.min(row)
      maxrow=np.max(row)
      mincol=np.min(col)
      maxcol=np.max(col)
      tempdata=CLgroup.variables["CL"][:,minrow:maxrow+1,mincol:maxcol+1]
      CLgroup.close()
      for i in range(dims[0]):
            for j in range(dims[1]):
                  for d in range(ndepth):
                        extCL[d,i,j]=tempdata[d,row[i,j]-minrow,col[i,j]-mincol]
                        if(np.isnan(extCL[d,i,j]) or extCL[d,i,j]==-999):
                              extCL[d,i,j]=-9999
                        if(extCL[d,i,j]==-9999):
                              extCL[d,i,j]=extCL[d-1,i,j]
                  for nl in range(Nlayer):
                        tmpCL[nl,i,j]=np.inner(layerfactor[nl,:],extCL[:,i,j])
      return tmpCL

def ExtractK_SCHFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon):
      K_SCHfile=sourceDir+"K_SCH.nc"
      K_SCHgroup = Dataset(K_SCHfile, "r")
      nlat=len(K_SCHgroup.dimensions["lat"])
      nlon=len(K_SCHgroup.dimensions["lon"])
      ndepth=len(K_SCHgroup.dimensions["depth"])
      lat0=K_SCHgroup.variables["lat"][0]
      latn=K_SCHgroup.variables["lat"][nlat-1]
      lon0=K_SCHgroup.variables["lon"][0]
      lonn=K_SCHgroup.variables["lon"][nlon-1]
      depth_in_org_dataset=K_SCHgroup.variables["depth"][:]
      layerfactor=DetermineLayers(Depths,depth_in_org_dataset)
      latcellsize=(latn-lat0)/(nlat-1)
      loncellsize=(lonn-lon0)/(nlon-1)
      dims=np.shape(lon)
      extK_SCH=np.zeros((ndepth,dims[0],dims[1]))
      tmpK_SCH=np.zeros((Nlayer,dims[0],dims[1]))
      row=np.zeros((dims[0],dims[1]),dtype=np.int)
      col=np.zeros((dims[0],dims[1]),dtype=np.int)
      for i in range(dims[0]):
            for j in range(dims[1]):
                  row[i,j]=int((lat[i,j]-lat0)/latcellsize)
                  col[i,j]=int((lon[i,j]-lon0)/loncellsize)      
      minrow=np.min(row)
      maxrow=np.max(row)
      mincol=np.min(col)
      maxcol=np.max(col)
      tempdata=K_SCHgroup.variables["K_SCH"][:,minrow:maxrow+1,mincol:maxcol+1]
      K_SCHgroup.close()
      for i in range(dims[0]):
            for j in range(dims[1]):
                  for d in range(ndepth):
                        extK_SCH[d,i,j]=tempdata[d,row[i,j]-minrow,col[i,j]-mincol]
                        if(np.isnan(extK_SCH[d,i,j]) or extK_SCH[d,i,j]==-999):
                              extK_SCH[d,i,j]=-9999
                        if(extK_SCH[d,i,j]!=-9999):
                              extK_SCH[d,i,j]=extK_SCH[d,i,j]/100.0/24.0/3600.0
                        else:
                              extK_SCH[d,i,j]=extK_SCH[d-1,i,j]
                  for nl in range(Nlayer):
                        tmpK_SCH[nl,i,j]=np.inner(layerfactor[nl,:],extK_SCH[:,i,j])
      return tmpK_SCH



def ExtractLAMBDFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon):
      LAMBDfile=sourceDir+"LAMBDA.nc"
      LAMBDgroup = Dataset(LAMBDfile, "r")
      nlat=len(LAMBDgroup.dimensions["lat"])
      nlon=len(LAMBDgroup.dimensions["lon"])
      ndepth=len(LAMBDgroup.dimensions["depth"])
      lat0=LAMBDgroup.variables["lat"][0]
      latn=LAMBDgroup.variables["lat"][nlat-1]
      lon0=LAMBDgroup.variables["lon"][0]
      lonn=LAMBDgroup.variables["lon"][nlon-1]
      depth_in_org_dataset=LAMBDgroup.variables["depth"][:]
      layerfactor=DetermineLayers(Depths,depth_in_org_dataset)
      latcellsize=(latn-lat0)/(nlat-1)
      loncellsize=(lonn-lon0)/(nlon-1)
      dims=np.shape(lon)
      extLAMBD=np.zeros((ndepth,dims[0],dims[1]))
      tmpLAMBD=np.zeros((Nlayer,dims[0],dims[1]))
      row=np.zeros((dims[0],dims[1]),dtype=np.int)
      col=np.zeros((dims[0],dims[1]),dtype=np.int)
      for i in range(dims[0]):
            for j in range(dims[1]):
                  row[i,j]=int((lat[i,j]-lat0)/latcellsize)
                  col[i,j]=int((lon[i,j]-lon0)/loncellsize)      
      minrow=np.min(row)
      maxrow=np.max(row)
      mincol=np.min(col)
      maxcol=np.max(col)
      tempdata=LAMBDgroup.variables["LAMBDA"][:,minrow:maxrow+1,mincol:maxcol+1]
      LAMBDgroup.close()
      for i in range(dims[0]):
            for j in range(dims[1]):
                  for d in range(ndepth):
                        extLAMBD[d,i,j]=tempdata[d,row[i,j]-minrow,col[i,j]-mincol]
                        if(np.isnan(extLAMBD[d,i,j]) or extLAMBD[d,i,j]==-999):
                              extLAMBD[d,i,j]=-9999
                        if(extLAMBD[d,i,j]==-9999):
                              extLAMBD[d,i,j]=extLAMBD[d-1,i,j]
                  for nl in range(Nlayer):
                        tmpLAMBD[nl,i,j]=np.inner(layerfactor[nl,:],extLAMBD[:,i,j])
      return tmpLAMBD



def ExtractPORFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon):
      PORfile=sourceDir+"POR.nc"
      PORgroup = Dataset(PORfile, "r")
      nlat=len(PORgroup.dimensions["lat"])
      nlon=len(PORgroup.dimensions["lon"])
      ndepth=len(PORgroup.dimensions["depth"])
      lat0=PORgroup.variables["lat"][0]
      latn=PORgroup.variables["lat"][nlat-1]
      lon0=PORgroup.variables["lon"][0]
      lonn=PORgroup.variables["lon"][nlon-1]
      depth_in_org_dataset=PORgroup.variables["depth"][:]
      layerfactor=DetermineLayers(Depths,depth_in_org_dataset)
      latcellsize=(latn-lat0)/(nlat-1)
      loncellsize=(lonn-lon0)/(nlon-1)
      dims=np.shape(lon)
      extPOR=np.zeros((ndepth,dims[0],dims[1]))
      tmpPOR=np.zeros((Nlayer,dims[0],dims[1]))
      row=np.zeros((dims[0],dims[1]),dtype=np.int)
      col=np.zeros((dims[0],dims[1]),dtype=np.int)
      for i in range(dims[0]):
            for j in range(dims[1]):
                  row[i,j]=int((lat[i,j]-lat0)/latcellsize)
                  col[i,j]=int((lon[i,j]-lon0)/loncellsize)      
      minrow=np.min(row)
      maxrow=np.max(row)
      mincol=np.min(col)
      maxcol=np.max(col)
      tempdata=PORgroup.variables["POR"][:,minrow:maxrow+1,mincol:maxcol+1]
      PORgroup.close()
      for i in range(dims[0]):
            for j in range(dims[1]):
                  for d in range(ndepth):
                        extPOR[d,i,j]=tempdata[d,row[i,j]-minrow,col[i,j]-mincol]
                        if(np.isnan(extPOR[d,i,j]) or extPOR[d,i,j]==-999):
                              extPOR[d,i,j]=-9999
                        if(extPOR[d,i,j]==-9999):
                              extPOR[d,i,j]=extPOR[d-1,i,j]
                  for nl in range(Nlayer):
                        tmpPOR[nl,i,j]=np.inner(layerfactor[nl,:],extPOR[:,i,j])
      return tmpPOR


def ExtractPSI_SFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon):
      PSI_Sfile=sourceDir+"PSI_S.nc"
      PSI_Sgroup = Dataset(PSI_Sfile, "r")
      nlat=len(PSI_Sgroup.dimensions["lat"])
      nlon=len(PSI_Sgroup.dimensions["lon"])
      ndepth=len(PSI_Sgroup.dimensions["depth"])
      lat0=PSI_Sgroup.variables["lat"][0]
      latn=PSI_Sgroup.variables["lat"][nlat-1]
      lon0=PSI_Sgroup.variables["lon"][0]
      lonn=PSI_Sgroup.variables["lon"][nlon-1]
      depth_in_org_dataset=PSI_Sgroup.variables["depth"][:]
      layerfactor=DetermineLayers(Depths,depth_in_org_dataset)
      latcellsize=(latn-lat0)/(nlat-1)
      loncellsize=(lonn-lon0)/(nlon-1)
      dims=np.shape(lon)
      extPSI_S=np.zeros((ndepth,dims[0],dims[1]))
      tmpPSI_S=np.zeros((Nlayer,dims[0],dims[1]))
      row=np.zeros((dims[0],dims[1]),dtype=np.int)
      col=np.zeros((dims[0],dims[1]),dtype=np.int)
      for i in range(dims[0]):
            for j in range(dims[1]):
                  row[i,j]=int((lat[i,j]-lat0)/latcellsize)
                  col[i,j]=int((lon[i,j]-lon0)/loncellsize)      
      minrow=np.min(row)
      maxrow=np.max(row)
      mincol=np.min(col)
      maxcol=np.max(col)
      tempdata=PSI_Sgroup.variables["PSI_S"][:,minrow:maxrow+1,mincol:maxcol+1]
      PSI_Sgroup.close()
      for i in range(dims[0]):
            for j in range(dims[1]):
                  for d in range(ndepth):
                        extPSI_S[d,i,j]=tempdata[d,row[i,j]-minrow,col[i,j]-mincol]
                        if(np.isnan(extPSI_S[d,i,j]) or extPSI_S[d,i,j]==-999):
                              extPSI_S[d,i,j]=-9999
                        if(extPSI_S[d,i,j]!=-9999):
                              extPSI_S[d,i,j]=extPSI_S[d,i,j]/100.0
                        else:
                              extPSI_S[d,i,j]=extPSI_S[d-1,i,j]
                  for nl in range(Nlayer):
                        tmpPSI_S[nl,i,j]=np.inner(layerfactor[nl,:],extPSI_S[:,i,j])
      return tmpPSI_S


def ExtractSIFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon):
      SIfile=sourceDir+"SI.nc"
      SIgroup = Dataset(SIfile, "r")
      nlat=len(SIgroup.dimensions["lat"])
      nlon=len(SIgroup.dimensions["lon"])
      ndepth=len(SIgroup.dimensions["depth"])
      lat0=SIgroup.variables["lat"][0]
      latn=SIgroup.variables["lat"][nlat-1]
      lon0=SIgroup.variables["lon"][0]
      lonn=SIgroup.variables["lon"][nlon-1]
      depth_in_org_dataset=SIgroup.variables["depth"][:]
      layerfactor=DetermineLayers(Depths,depth_in_org_dataset)
      latcellsize=(latn-lat0)/(nlat-1)
      loncellsize=(lonn-lon0)/(nlon-1)
      dims=np.shape(lon)
      extSI=np.zeros((ndepth,dims[0],dims[1]))
      tmpSI=np.zeros((Nlayer,dims[0],dims[1]))
      row=np.zeros((dims[0],dims[1]),dtype=np.int)
      col=np.zeros((dims[0],dims[1]),dtype=np.int)
      for i in range(dims[0]):
            for j in range(dims[1]):
                  row[i,j]=int((lat[i,j]-lat0)/latcellsize)
                  col[i,j]=int((lon[i,j]-lon0)/loncellsize)      
      minrow=np.min(row)
      maxrow=np.max(row)
      mincol=np.min(col)
      maxcol=np.max(col)
      tempdata=SIgroup.variables["SI"][:,minrow:maxrow+1,mincol:maxcol+1]
      SIgroup.close()
      for i in range(dims[0]):
            for j in range(dims[1]):
                  for d in range(ndepth):
                        extSI[d,i,j]=tempdata[d,row[i,j]-minrow,col[i,j]-mincol]
                        if(np.isnan(extSI[d,i,j]) or extSI[d,i,j]==-999):
                              extSI[d,i,j]=-9999
                        if(extSI[d,i,j]==-9999):
                              extSI[d,i,j]=extSI[d-1,i,j]
                  for nl in range(Nlayer):
                        tmpSI[nl,i,j]=np.inner(layerfactor[nl,:],extSI[:,i,j])
      return tmpSI

def ExtractSOMFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon):
      SOMfile=sourceDir+"SOM.nc"
      SOMgroup = Dataset(SOMfile, "r")
      nlat=len(SOMgroup.dimensions["lat"])
      nlon=len(SOMgroup.dimensions["lon"])
      ndepth=len(SOMgroup.dimensions["depth"])
      lat0=SOMgroup.variables["lat"][0]
      latn=SOMgroup.variables["lat"][nlat-1]
      lon0=SOMgroup.variables["lon"][0]
      lonn=SOMgroup.variables["lon"][nlon-1]
      depth_in_org_dataset=SOMgroup.variables["depth"][:]
      layerfactor=DetermineLayers(Depths,depth_in_org_dataset)
      latcellsize=(latn-lat0)/(nlat-1)
      loncellsize=(lonn-lon0)/(nlon-1)
      dims=np.shape(lon)
      extSOM=np.zeros((ndepth,dims[0],dims[1]))
      tmpSOM=np.zeros((Nlayer,dims[0],dims[1]))
      row=np.zeros((dims[0],dims[1]),dtype=np.int)
      col=np.zeros((dims[0],dims[1]),dtype=np.int)
      for i in range(dims[0]):
            for j in range(dims[1]):
                  row[i,j]=int((lat[i,j]-lat0)/latcellsize)
                  col[i,j]=int((lon[i,j]-lon0)/loncellsize)      
      minrow=np.min(row)
      maxrow=np.max(row)
      mincol=np.min(col)
      maxcol=np.max(col)
      tempdata=SOMgroup.variables["SOM"][:,minrow:maxrow+1,mincol:maxcol+1]
      SOMgroup.close()
      for i in range(dims[0]):
            for j in range(dims[1]):
                  for d in range(ndepth):
                        extSOM[d,i,j]=tempdata[d,row[i,j]-minrow,col[i,j]-mincol]
                        if(np.isnan(extSOM[d,i,j]) or extSOM[d,i,j]==-999):
                              extSOM[d,i,j]=-9999
                        if(extSOM[d,i,j]==-9999):
                              extSOM[d,i,j]=extSOM[d-1,i,j]
                  for nl in range(Nlayer):
                        tmpSOM[nl,i,j]=np.inner(layerfactor[nl,:],extSOM[:,i,j])
      return tmpSOM


def ExtractTHETASFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon):
      THSCHfile=sourceDir+"THSCH.nc"
      THSCHgroup = Dataset(THSCHfile, "r")
      nlat=len(THSCHgroup.dimensions["lat"])
      nlon=len(THSCHgroup.dimensions["lon"])
      ndepth=len(THSCHgroup.dimensions["depth"])
      lat0=THSCHgroup.variables["lat"][0]
      latn=THSCHgroup.variables["lat"][nlat-1]
      lon0=THSCHgroup.variables["lon"][0]
      lonn=THSCHgroup.variables["lon"][nlon-1]
      depth_in_org_dataset=THSCHgroup.variables["depth"][:]
      layerfactor=DetermineLayers(Depths,depth_in_org_dataset)
      latcellsize=(latn-lat0)/(nlat-1)
      loncellsize=(lonn-lon0)/(nlon-1)
      dims=np.shape(lon)
      extTHSCH=np.zeros((ndepth,dims[0],dims[1]))
      tmpTHSCH=np.zeros((Nlayer,dims[0],dims[1]))
      row=np.zeros((dims[0],dims[1]),dtype=np.int)
      col=np.zeros((dims[0],dims[1]),dtype=np.int)
      for i in range(dims[0]):
            for j in range(dims[1]):
                  row[i,j]=int((lat[i,j]-lat0)/latcellsize)
                  col[i,j]=int((lon[i,j]-lon0)/loncellsize)      
      minrow=np.min(row)
      maxrow=np.max(row)
      mincol=np.min(col)
      maxcol=np.max(col)
      tempdata=THSCHgroup.variables["THSCH"][:,minrow:maxrow+1,mincol:maxcol+1]
      THSCHgroup.close()
      for i in range(dims[0]):
            for j in range(dims[1]):
                  for d in range(ndepth):
                        extTHSCH[d,i,j]=tempdata[d,row[i,j]-minrow,col[i,j]-mincol]
                        if(np.isnan(extTHSCH[d,i,j]) or extTHSCH[d,i,j]==-999):
                              extTHSCH[d,i,j]=-9999
                        if(extTHSCH[d,i,j]==-9999):
                              extTHSCH[d,i,j]=extTHSCH[d-1,i,j]
                  for nl in range(Nlayer):
                        tmpTHSCH[nl,i,j]=np.inner(layerfactor[nl,:],extTHSCH[:,i,j])
      return tmpTHSCH

def ExtractTHETARFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon):
      THRfile=sourceDir+"THR.nc"
      THRgroup = Dataset(THRfile, "r")
      nlat=len(THRgroup.dimensions["lat"])
      nlon=len(THRgroup.dimensions["lon"])
      ndepth=len(THRgroup.dimensions["depth"])
      lat0=THRgroup.variables["lat"][0]
      latn=THRgroup.variables["lat"][nlat-1]
      lon0=THRgroup.variables["lon"][0]
      lonn=THRgroup.variables["lon"][nlon-1]
      depth_in_org_dataset=THRgroup.variables["depth"][:]
      layerfactor=DetermineLayers(Depths,depth_in_org_dataset)
      latcellsize=(latn-lat0)/(nlat-1)
      loncellsize=(lonn-lon0)/(nlon-1)
      dims=np.shape(lon)
      extTHR=np.zeros((ndepth,dims[0],dims[1]))
      tmpTHR=np.zeros((Nlayer,dims[0],dims[1]))
      row=np.zeros((dims[0],dims[1]),dtype=np.int)
      col=np.zeros((dims[0],dims[1]),dtype=np.int)
      for i in range(dims[0]):
            for j in range(dims[1]):
                  row[i,j]=int((lat[i,j]-lat0)/latcellsize)
                  col[i,j]=int((lon[i,j]-lon0)/loncellsize)      
      minrow=np.min(row)
      maxrow=np.max(row)
      mincol=np.min(col)
      maxcol=np.max(col)
      tempdata=THRgroup.variables["THR"][:,minrow:maxrow+1,mincol:maxcol+1]
      THRgroup.close()
      for i in range(dims[0]):
            for j in range(dims[1]):
                  for d in range(ndepth):
                        extTHR[d,i,j]=tempdata[d,row[i,j]-minrow,col[i,j]-mincol]
                        if(np.isnan(extTHR[d,i,j]) or extTHR[d,i,j]==-999):
                              extTHR[d,i,j]=-9999
                        if(extTHR[d,i,j]==-9999):
                              extTHR[d,i,j]=extTHR[d-1,i,j]
                  for nl in range(Nlayer):
                        tmpTHR[nl,i,j]=np.inner(layerfactor[nl,:],extTHR[:,i,j])
      return tmpTHR

def ExtractFCFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon):
      FCfile=sourceDir+"TH33.nc"
      FCgroup = Dataset(FCfile, "r")
      nlat=len(FCgroup.dimensions["lat"])
      nlon=len(FCgroup.dimensions["lon"])
      ndepth=len(FCgroup.dimensions["depth"])
      lat0=FCgroup.variables["lat"][0]
      latn=FCgroup.variables["lat"][nlat-1]
      lon0=FCgroup.variables["lon"][0]
      lonn=FCgroup.variables["lon"][nlon-1]
      depth_in_org_dataset=FCgroup.variables["depth"][:]
      layerfactor=DetermineLayers(Depths,depth_in_org_dataset)
      latcellsize=(latn-lat0)/(nlat-1)
      loncellsize=(lonn-lon0)/(nlon-1)
      dims=np.shape(lon)
      extFC=np.zeros((ndepth,dims[0],dims[1]))
      tmpFC=np.zeros((Nlayer,dims[0],dims[1]))
      row=np.zeros((dims[0],dims[1]),dtype=np.int)
      col=np.zeros((dims[0],dims[1]),dtype=np.int)
      for i in range(dims[0]):
            for j in range(dims[1]):
                  row[i,j]=int((lat[i,j]-lat0)/latcellsize)
                  col[i,j]=int((lon[i,j]-lon0)/loncellsize)      
      minrow=np.min(row)
      maxrow=np.max(row)
      mincol=np.min(col)
      maxcol=np.max(col)
      tempdata=FCgroup.variables["TH33"][:,minrow:maxrow+1,mincol:maxcol+1]
      FCgroup.close()
      for i in range(dims[0]):
            for j in range(dims[1]):
                  for d in range(ndepth):
                        extFC[d,i,j]=tempdata[d,row[i,j]-minrow,col[i,j]-mincol]
                        if(np.isnan(extFC[d,i,j]) or extFC[d,i,j]==-999):
                              extFC[d,i,j]=-9999
                        if(extFC[d,i,j]==-9999):
                              extFC[d,i,j]=extFC[d-1,i,j]
                  for nl in range(Nlayer):
                        tmpFC[nl,i,j]=np.inner(layerfactor[nl,:],extFC[:,i,j])
      return tmpFC


def ExtractSoilDataFromDaiYongJiu(sourceDir,OutputFile,Depths,lat,lon):
      dims=np.shape(lon)
      Nlayer=len(Depths)
      print Nlayer
      BD=ExtractBDFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon)
      SI=ExtractSIFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon)
      SA=ExtractSAFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon)
      CL=ExtractCLFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon)
      SOM=ExtractSOMFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon)
      PSI_S=ExtractPSI_SFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon)
      LAMBD=ExtractLAMBDFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon)
      THETAS=ExtractTHETASFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon)
      KS=ExtractK_SCHFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon) 
      THR=ExtractTHETARFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon) 
      FC=ExtractFCFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon) 
      
      
      soilpara = Dataset(OutputFile,"w",format="NETCDF4")
      lat   = soilpara.createDimension("NY", dims[0])
      lon   = soilpara.createDimension("NX", dims[1])
      depth = soilpara.createDimension("depth", Nlayer)
      
      BDs= soilpara.createVariable("bd","f4",("depth","NY","NX",))
      BDs.missing_value="-9999"
      BDs.units="kg/m3"
      BDs.long_name="Bulk Density"
      BDs[:,:,:]=BD[:,:,:]
      
      SAs= soilpara.createVariable("sa","f4",("depth","NY","NX",))
      SAs.missing_value="-9999"
      SAs.units="0-100 %"
      SAs.long_name="Sand content"
      SAs[:,:,:]=SA[:,:,:]
      
      SIs= soilpara.createVariable("si","f4",("depth","NY","NX",))
      SIs.missing_value="-9999"
      SIs.units="0-100 %"
      SIs.long_name="Silt content"
      SIs[:,:,:]=SI[:,:,:]   
      
      CLs= soilpara.createVariable("cl","f4",("depth","NY","NX",))
      CLs.missing_value="-9999"
      CLs.units="0-100 %"
      CLs.long_name="Clay content"
      CLs[:,:,:]=CL[:,:,:]   
      
      SOMs= soilpara.createVariable("om","f4",("depth","NY","NX",))
      SOMs.missing_value="-9999"
      SOMs.units="0-100 %"
      SOMs.long_name="Organic content"
      SOMs[:,:,:]=SOM[:,:,:]

      KSs= soilpara.createVariable("ks","f4",("depth","NY","NX",))
      KSs.missing_value="-9999"
      KSs.units="m/s"
      KSs.long_name="Saturate hydraulic conductivity"
      KSs[:,:,:]=KS[:,:,:]    
      
      PSI_Ss= soilpara.createVariable("psi_s","f4",("depth","NY","NX",))
      PSI_Ss.missing_value="-9999"
      PSI_Ss.units="m"
      PSI_Ss.long_name="Saturated capillary potential"
      PSI_Ss[:,:,:]=PSI_S[:,:,:]    
      
      LAMBDs= soilpara.createVariable("b","f4",("depth","NY","NX",))
      LAMBDs.missing_value="-9999"
      LAMBDs.units="-"
      LAMBDs.long_name="Pore size distribution index"
      LAMBDs[:,:,:]=1.0/LAMBD[:,:,:]        
          
      THETASs= soilpara.createVariable("theta_s","f4",("depth","NY","NX",))
      THETASs.missing_value="-9999"
      THETASs.units="cm3 cm-3"
      THETASs.long_name="Saturated water content"
      THETASs[:,:,:]=THETAS[:,:,:]
      
      THRs= soilpara.createVariable("thr","f4",("depth","NY","NX",))
      THRs.missing_value="-9999"
      THRs.units="cm3 cm-3"
      THRs.long_name="Residual water content"
      THRs[:,:,:]=THR[:,:,:]
      
      FCs= soilpara.createVariable("fc","f4",("depth","NY","NX",))
      FCs.missing_value="-9999"
      FCs.units="cm3 cm-3"
      FCs.long_name="Field water capacity"
      FCs[:,:,:]=FC[:,:,:]
      
      ZSs= soilpara.createVariable("zs","f4",("depth",))
      ZSs.missing_value="-9999"
      ZSs.units="m/s"
      ZSs.long_name="Saturate hydraulic conductivity"
      ZSs[:]=Depths[:]                          
      soilpara.close()

