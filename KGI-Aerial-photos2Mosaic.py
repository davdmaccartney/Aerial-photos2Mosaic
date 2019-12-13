# import the necessary packages
import os
import gc
from imutils import paths
import argparse
import cv2
import numpy as np
from osgeo import gdal
import time
import easygui
import fnmatch
import time
import sys
import re
import ogr
import osr
import subprocess
import mmap


def striplist(l):
    return([x.strip() for x in l])

def cls():
    os.system('cls' if os.name=='nt' else 'clear')

def compLat_Long(degs, mins, secs, comp_dir):
    return (degs + (mins / 60) + (secs / 3600)) * comp_dir

def mapcount(filename):
    f = open(filename, "r+")
    buf = mmap.mmap(f.fileno(), 0)
    lines = 0
    readline = buf.readline
    while readline():
        lines += 1
    return lines-1


def getSignOf(chifre):
    if chifre >= 0:
        return 1
    else:
        return -1


def hrp2opk(Roll, Pitch, heading):
    Roll = np.deg2rad(Roll)
    Pitch   = np.deg2rad(Pitch)
    heading = np.deg2rad(heading)

    A_SINH = np.sin(heading)
    A_SINR = np.sin(Roll)
    A_SINP = np.sin(Pitch)

    A_COSH = np.cos(heading)
    A_COSR = np.cos(Roll)
    A_COSP = np.cos(Pitch)

    MX = np.zeros((3, 3))
    MX[0][0] =  (A_COSH *A_COSR) + (A_SINH*A_SINP*A_SINR)
    MX[0][1] =  (-A_SINH*A_COSR)+(A_COSH*A_SINP*A_SINR)
    MX[0][2] =   -A_COSP*A_SINR

    MX[1][0] = A_SINH*A_COSP
    MX[1][1] = A_COSH*A_COSP
    MX[1][2] = A_SINP


    MX[2][0] = (A_COSH*A_SINR)-(A_SINH*A_SINP*A_COSR)
    MX[2][1] = (-A_SINH*A_SINR)-(A_COSH*A_SINP*A_COSR)
    MX[2][2] =  A_COSP*A_COSR

    P = np.zeros((3, 3))
    P[0][0] = MX[0][0]
    P[0][1] = MX[1][0]
    P[0][2] = MX[2][0]
    
    P[1][0] = MX[0][1]
    P[1][1] = MX[1][1]
    P[1][2] = MX[2][1]
    
    P[2][0] = MX[2][0]
    P[2][1] = MX[1][2]
    P[2][2] = MX[2][2]

    Omega = 0
    Phi   = 0
    Kappa = 0

    Omega = np.arctan(-P[2][1]/P[2][2])
    Phi = np.arcsin(P[2][2])
    Kappa = np.arctan(-P[1][0]/P[0][0])

    Phi   = abs(np.arcsin(P[2][0]))
    Phi = Phi * getSignOf(P[2][0])
    Omega = abs(np.arccos((P[2][2] / np.cos(Phi))))
    Omega = Omega * (getSignOf(P[2][1] / P[2][2]*-1))
    Kappa = np.arccos(P[0][0] / np.cos(Phi))

    if getSignOf(P[0][0]) == getSignOf((P[1][0] / P[0][0])):
        Kappa = Kappa * -1

    Omega = np.rad2deg(Omega)
    Phi = np.rad2deg(Phi)
    Kappa = np.rad2deg(Kappa)
   
    return(Omega,Phi,Kappa)

def update_progress(progress):
    barLength = 30 # Modify this to change the length of the progress bar
    
    block = int(round(barLength*progress))
    text = "\rPercent: [{0}] {1}% ".format( "#"*block + "-"*(barLength-block), int(progress*100))
    sys.stdout.write(text)
    sys.stdout.flush()
    print('\n')

def transform_wgs84_to_utm(lon, lat, alt):    
    def get_utm_zone(longitude):
        return (int(1+(longitude+180.0)/6.0))

    def is_northern(latitude):
        """
        Determines if given latitude is a northern for UTM
        """
        if (latitude < 0.0):
            return 0
        else:
            return 1

    utm_coordinate_system = osr.SpatialReference()
    utm_coordinate_system.SetWellKnownGeogCS("WGS84") # Set geographic coordinate system to handle lat/lon  
    utm_coordinate_system.SetUTM(get_utm_zone(lon), is_northern(lat))
    wgs84_coordinate_system = utm_coordinate_system.CloneGeogCS() # Clone ONLY the geographic coordinate system 

    # create transform component
    wgs84_to_utm_transform = osr.CoordinateTransformation(wgs84_coordinate_system, utm_coordinate_system) # (<from>, <to>)
    return wgs84_to_utm_transform.TransformPoint(lon, lat, alt) # returns easting, northing, altitude    

def epsgtoepsg(EPSG_in, lon, lat, EPSG_out):


 InSR = osr.SpatialReference()
 InSR.ImportFromEPSG(EPSG_in)      
 OutSR = osr.SpatialReference()
 OutSR.ImportFromEPSG(EPSG_out)     

 Point = ogr.Geometry(ogr.wkbPoint)
 Point.AddPoint(lon,lat) # use your coordinates here
 Point.AssignSpatialReference(InSR)    # tell the point what coordinates it's in
 Point.TransformTo(OutSR)              # project it to the out spatial reference

 return (Point.GetX(),Point.GetY())


def altitude(x, y, DTM):
   
    TL_x, x_res, _, TL_y, _, y_res = DTM.GetGeoTransform()
    x_index = (x - TL_x) / x_res
    y_index = (y - TL_y) / y_res
    array = DTM.ReadAsArray()
    pixel_val = array[int(y_index), int(x_index)]
    del array


    return pixel_val

def rotation(nx, ny, z0, omega, phi, kappa):
      #heading : ψ, roll: Φ, pitch:Θ 

      # apply yaw (around z) / yaw kappa
      xr1= nx * np.cos(kappa) - ny * np.sin(kappa)
      yr1= nx * np.sin(kappa) + ny * np.cos(kappa)
      zr1 = z0

      #apply pitch (around x) / pitch / omega
      xr2 = xr1
      yr2 = yr1 * np.cos(omega) - zr1 * np.sin(omega)
      zr2 = yr1 * np.sin(omega) + zr1 * np.cos(omega)

      # apply roll (around y) / roll / ohi
      xr3 = xr2 * np.cos(phi) - zr2 * np.sin(phi)
      yr3 = yr2
      zr3 = xr2 * np.sin(phi) + zr2 * np.cos(phi)

      return xr3, yr3, zr3

dirname = easygui.diropenbox(msg=None, title="Please select the directory", default=None )
file_Dir = os.path.basename(dirname)
DTM = gdal.Open(dirname+"\DEM_Copernicus\eu_dem_V11__B2.tif")
ci=0
cls()


# definition de la camera
Micron = 4.6000
foc = 51.6611
pRot = 0
xwidth=640
yheight=480
# reduction des photos
oVer=.20


if not os.path.isfile(dirname+"/ExifLog.csv"):
       print('File does not exist')
       sys.exit(0)

total_con = mapcount(dirname+"/ExifLog.csv")


dir = os.path.join(dirname,"mosaic")
if not os.path.exists(dir):
    os.mkdir(dir)

f = open(dir+"\list.txt", "w")

with open(dirname+"/ExifLog.csv") as fp:
   for line in fp:
     data = line.strip().split(',')
     data = list(filter(None, data))
     if data[0] !='Filename':
      nom = data[0]
      nom = nom.replace('.IIQ','')
 
      Longitude = float(data[4])
      Latitude = float(data[5])

      Altitude = data[6]
      Altitude = Altitude.replace(' m','')
      z0 = float(Altitude)
      x0,y0 = epsgtoepsg(4326, Longitude, Latitude, 2154)

      xt,yt = epsgtoepsg(2154, x0, y0, 3035)
      z0 = z0 - altitude(xt,yt,DTM)
      #Pitch,Roll,Yaw
      p = float(data[7])
      r = float(data[8])
      y = float(data[9])

      Omega,Phi,Kappa= hrp2opk(r, p, -y)
      Omega = -Omega
      #Kappa = 90-Kappa
      print(Omega,Phi,Kappa)
      
      Omega = np.deg2rad(Omega)
      Phi = np.deg2rad(Phi)
      Kappa = np.deg2rad(Kappa)
 
      scale_fact = 11608/640
      scale = (z0/float(foc))*1000
      scale = scale*scale_fact

      GSD = (float(scale)*Micron/10000)/100

      iWEmini = int((xwidth*oVer)//2)
      iWNmini = int((yheight*oVer)//2)

      #imgWorldE=xwidth*GSD
      #imgWorldN=yheight*GSD

      imgWorldE=(xwidth-(iWEmini*2))*GSD
      imgWorldN=(yheight-(iWNmini*2))*GSD


      # haut gauche
      x1 = (x0-(imgWorldE/2))
      y1 = (y0+(imgWorldN/2))
      # haut droite
      x2 = (x0+(imgWorldE/2))
      y2 = (y0+(imgWorldN/2))
      # bas droite
      x3 = (x0+(imgWorldE/2))
      y3 = (y0-(imgWorldN/2))
      # bas gauche
      x4 = (x0-(imgWorldE/2))
      y4 = (y0-(imgWorldN/2))


      # all rotations counter clockwise

      # point haut gauche
      # tranlation
      nx = x1-x0
      ny = y1-y0

      # rotations
      xr3, yr3, zr3 = rotation(nx, ny, z0, Omega, Phi, Kappa)

      # tranlation back
      nx1r = xr3 + x0
      ny1r = yr3 + y0

      # point haut droite
      # tranlation
      nx = x2-x0
      ny = y2-y0

       # rotations
      xr3, yr3, zr3 = rotation(nx, ny, z0, Omega, Phi, Kappa)

      # tranlation back
      nx2r = xr3 + x0
      ny2r = yr3 + y0
      

      # point bas droite
      # tranlation
      nx = x3-x0
      ny = y3-y0

      # rotations
      xr3, yr3, zr3 = rotation(nx, ny, z0, Omega, Phi, Kappa)

      # tranlation back
      nx3r = xr3 + x0
      ny3r = yr3 + y0
      

      # point bas gauche
      # tranlation
      nx = x4-x0
      ny = y4-y0

      # rotations rotation(nx, ny, z0, pitch, roll, yaw) 
      xr3, yr3, zr3 = rotation(nx, ny, z0, Omega, Phi, Kappa)

      # tranlation back
      nx4r = xr3 + x0
      ny4r = yr3 + y0


      ci  += 1

     
      tranlatecrop = 'gdal_translate -of GTiff -srcwin ' + str(iWEmini) + ' ' + str(iWNmini) + ' ' + str(xwidth-(iWEmini*2)) + ' ' + str(yheight-(iWNmini*2)) + \
      ' '+ dirname + '\\Miniatures' + '\\' + nom +'.jpg' +' '+ dir + '\\' + nom +'_C.tif'

      tranlate = 'gdal_translate -of GTiff -a_srs EPSG:2154 -gcp 0 0 '+ str(nx1r)+' '+ str(ny1r) + \
      ' -gcp '+str(xwidth-(iWEmini*2))+' 0 '+ str(nx2r) +' '+ str(ny2r) + \
      ' -gcp '+ str(xwidth-(iWEmini*2)) +' ' + str(yheight-(iWNmini*2)) + ' ' + str(nx3r) +' '+ str(ny3r) +\
      ' -gcp 0 ' + str(yheight-(iWNmini*2)) + ' ' + str(nx4r) +' '+ str(ny4r) +\
      ' '+ dir + '\\' + nom +'_C.tif' +' '+ dir + '\\' + nom +'_T.tif'


      warp = 'gdalwarp -co COMPRESS=NONE -tr 2 2 -r bilinear -s_srs EPSG:2154 -dstalpha -multi -of GTiff '+ dir + '\\' + nom +'_T.tif' + ' ' + dir + '\\' + nom +'.tif'

      os.system('set OSGEO4W_ROOT=C:\OSGEO4~1')
      os.system('set path=%OSGEO4W_ROOT%\bin;%WINDIR%\system32;%WINDIR%;%WINDIR%\system32\WBem;C:\LAStools\bin\laszip.exe')
      os.system('SET GDAL_DATA=C:\OSGEO4~1\share\gdal')
      os.system('SET GDAL_DATA=%OSGEO4W_ROOT%\share\gdal')
      os.system('SET GDAL_DRIVER_PATH=%OSGEO4W_ROOT%\bin\gdalplugins')
      os.system('SET GDAL_DRIVER_PATH=C:\OSGEO4~1\bin\gdalplugins')
      os.system('SET GEOTIFF_CSV=C:\OSGEO4~1\share\epsg_csv')
      os.system('set JPEGMEM=100000')
      os.system('SET GDAL_DATA=%OSGEO4W_ROOT%\share\gdal')
      os.system('SET GDAL_DRIVER_PATH=%OSGEO4W_ROOT%\bin\gdalplugins')
      os.system('SET PROJ_LIB=%OSGEO4W_ROOT%\share\proj')

      os.system(tranlatecrop)
      os.system(tranlate)
      os.remove(dir + '\\' + nom +'_C.tif')
      print('\n')
      os.system(warp)
      os.remove(dir + '\\' + nom +'_T.tif')

      f.write(dir + '\\' + nom +'.tif'+'\n')
      f.flush
      cls()
      update_progress(ci/total_con)

cls()
f.close
f = open(dir + '\list.txt')
f.close()
del DTM

merging = 'python C:/OSGeo4W64/bin/gdal_merge.py -of GTiff -co TFW=YES -o ' + dirname + '\merged.tif --optfile '+ dir + '\list.txt'
print('merging\n')
os.system(merging)
delforlder = 'rmdir /S /Q '+ dir
os.system(delforlder)
cls()
print('Done')