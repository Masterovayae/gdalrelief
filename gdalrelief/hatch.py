#!/usr/bin/env python
# a bar plot with errorbars
import numpy as np
import numpy.ma as ma
import numpy.ma as maa
from osgeo import gdal
import matplotlib.pyplot as plt
import gdalrelief.diff as diff
from scipy.signal import argrelextrema
import collections
from osgeo.gdalconst import *
from common import *
from scipy import interpolate
from matplotlib import cm
from mpl_toolkits.mplot3d.axes3d import get_test_data
from itertools import product


TESTRASTER_GOLOUSTNOYE="../data/Goloustnoye/ALTITUDE 1Trim.grd"
TESTRASTER_OLKHON="../data/Olkhon/dem.gtiff"
SAFE_GDAL_FORMAT="GTiff"


class RasterSection(RasterProcessor):
    """
    """

    def __init__(self, raster, hatch=None):
        RasterProcessor.__init__(self, raster)
        self.set_hatch(hatch)

    def set_hatch(self, hatch):
        self.hatch=hatch

    def __call__(self, layer):
        for x,y in self.line():
            sec.append(band[x,y])
        return np.array(sec)


    def scan_sq(self, layer):
        band=RasterProcessor.__call__(self, layer)
        a=self.current_band=band
        h,w=a.shape
        print (h,w)
        yield from self.ss(0,0,w,h)

    DD=((0,1),(0,-1),(1,0),(-1,0),(1,1),(1,-1),(-1,1),(-1,-1))

    def mat_scan(self, layer):
        band=RasterProcessor.__call__(self, layer)
        a=self.current_band=band
        h,w=a.shape
        print (h,w)
        mi1=np.ones((h,w), dtype=bool)
        ma1=np.ones((h,w), dtype=bool)
        mi1=ma.array(mi1,mask=a.mask)
        ma1=ma.array(ma1,mask=a.mask)
        #mi1=ma1
        for dr,dc in self.__class__.DD:
            self.d(a,mi1,True , dr, dc)
            self.d(a,ma1,False, dr, dc)

        return mi1, ma1

    def d(self, a,m,miop, dr,dc):
        b=a[:,:]
        c=a[:,:]
        q=m[:,:]
        if dc>0:
            c=c[:-1,:]
            q=q[:-1,:]
            b=b[1:,:]
        elif dc<0:
            c=c[1:,:]
            q=q[1:,:]
            b=b[:-1,:]

        if dr>0:
            c=c[:,:-1]
            q=q[:,:-1]
            b=b[:,1:]
        else:
            c=c[:,1:]
            q=q[:,1:]
            b=b[:,:-1]
        if miop:
            d=c<b
        else:
            d=c>b
        q &= d


    def check_ext(self, r,c):
        """
        """
        a=self.current_band
        h,w=a.shape
        h=hmin=hmax=a[r,c]
        if np.isnan(h):
            return
        ma=mi=True
        for dc,dr in self.__class__.DD:
            cr=r+dr
            cc=c+dc
            if cc>=w:
                continue
            if cr>=h:
                continue
            if cc<0:
                continue
            if cr<0:
                continue
            print (cr,cc)
            h=a[cr,cc]
            if np.isnan(h):
                continue

            if h>hmax:
                ma=False
            if h<hmin:
                mi=False

        if ma:
            yield (True, r,c,hmax) # True means max
        if mi:
            yield (False, r,c,hmin) # False means min

    def vol_extr(self, rs, ae):
        a=rs.current_band
        h,w=a.shape
        #a=np.array(a)
        #print ("a", a)

        interpol='cubic'
        f = interpolate.interp2d(ae[:,0], ae[:,0], ae[:,0], kind=interpol)

        #ma.masked_values(band, nan_value)
        #mas=self.band
        #m=mas.mask
        #print('mas',m)
        #self.graf_3d(f,h,w)
        prec=1
        xnew = np.arange(0, h, prec)
        ynew = np.arange(0, w, prec)
        znew = f(xnew, ynew)
        #~ for r,c in product(range(h),range(w)):
            #~ z=a[r,c]
            #~ if znew[r,c]==True:
                #~ ad.append((c,r,z))
            #~ if a[r,c]==True:
                #~ a.append((c,r,z))

        znew = znew.transpose()
        #znew=np.array(znew)
        znew=ma.array(znew, mask=a.mask)
        #print ("znew", znew)
        v=a-znew
        #print (a)
        v1=0
        #for i in range(len(v)):
        #    for j in range(len(v[i])):
        #        if v[i,j]>0:
        #
        #            v1=v1+v[i][j]
        v1=np.sum(v)
        #print (v1)
        return v1,v

    def graf_3d(self,f,h,w):
        prec=20
        xnew = np.arange(0, w, prec)
        ynew = np.arange(0, h, prec)
        znew = f(xnew, ynew)
        xnew, ynew = np.meshgrid(xnew, ynew)

        fig = plt.figure(figsize=plt.figaspect(1))

        ax = fig.add_subplot(1, 1, 1, projection='3d')

        surf = ax.plot_surface(xnew, ynew, znew, rstride=1, cstride=1, cmap=cm.coolwarm,
                               linewidth=0, antialiased=True)
        ax.set_zlim3d(400, 900)

        fig.colorbar(surf, shrink=0.5, aspect=10)


        plt.show()
        return

# TEST


def test_2():
    name=TESTRASTER_GOLOUSTNOYE
    rs=RasterSection(raster=name)
    rs.info()

    mi,ma=rs.mat_scan(4)
    a=rs.current_band

    b=maa.array(a, mask=~mi)
    rs.save("mi", b)
    b=maa.array(a, mask=~ma)
    rs.save("ma", b)
    print('mi',mi)
    a=np.array(a)
    h,w=a.shape
    #print(mi,ma)
   # m=rs.bant(4)
    au=[]
    ad=[]

    for r,c in product(range(h),range(w)):
        z=a[r,c]
        if mi[r,c]==True:
            ad.append((c,r,z))
        if ma[r,c]==True:
            au.append((c,r,z))

    au=np.array(au)
    ad=np.array(ad)

    vu,ru=rs.vol_extr(rs,au)
    rs.save("{}-{}-upper".format(name,vu), ru)
    print('vu=',vu)
    vd,rd=rs.vol_extr(rs,ad)
    rs.save("{}-{}-lower".format(name,vd), rd)
    print('vd=',vd)

    vt=vu-vd
    print('vt=',vt)


if __name__=="__main__":
    # register all of the GDAL drivers
    gdal.AllRegister()

    test_2()

    quit()