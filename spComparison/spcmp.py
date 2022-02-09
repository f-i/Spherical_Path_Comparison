#!/usr/bin/python3

'''-------------------------------------------------------------------------###
Created on 5May2016
Modified on 9Feb2022

@__author__	:	Chenjian Fu
@__email__	:	cfu3@kent.edu
@__purpose__	:	To quantitatively compare paleomagnetic APWPs
@__version__	:	0.8.5
@__license__	:	GNU General Public License v3.0

Spherical Path Comparison (spComparison) package is created for quantitatively
measuring similarity of spherical paths, particularly the paleomagnetic apparent
polar wander paths (APWPs) of tectonic plates. It is powered by PmagPy
(https://pmagpy.github.io/) and GMT (http://gmt.soest.hawaii.edu/).
Copyright (C) 2016-2021 @__author__

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <https://www.gnu.org/licenses/>.
-------------------------------------------------------------------------------
Environment:
    *NIX(-like) Bash Shell with "Python3, NumPy, GMT, bc" installed
-------------------------------------------------------------------------------
TODO:
    1. Tidy functions up into classes
###-------------------------------------------------------------------------'''

from os import makedirs, path
from subprocess import Popen, PIPE
from random import random
import re, numpy as np

PLATE_V_MAX_PAST=30  #according to Swanson-Hysell etal.2009, Kulakov etal.2014; today it's about 15.44cm/yr, DeMets etal.2010
#divisor for mean length dif=PLATE_V_MAX_PAST/11.1195051975  #i.e. about 2.7 degree/myr, magnitude of velocity
POL_WAND_DIR_DIF_MAX=180

#Returns azimuth of geodesic ({0},{1}) to ({2},{3}) Source: @__author__, Oct2015
#/home/g/Desktop/git/gpmdb/src/core/azi {0} {1} {2} {3}
AZI="""
gmt mapproject -Af{0}/{1} -fg -o2 <<< '{2} {3}'
"""

#Returns absolute value of (({0}-{1}) % 360)        Source: @__author__, Jan2016
SMALLER_ANGLE_REMAINDER="""
gmt math -Q {0} {1} SUB 360 FMOD ABS =
"""

#Returns one sample within error ellipse of pole (location ({0},{1}), azimuth of
#major axis {4})                                    Source: @__author__, Oct2016
ASSIGN_AZI4ROTATED_ELLIP="""
gmt set PROJ_ELLIPSOID Sphere
angle=`gmt vector -S0/0 -TD -fg <<< "{0} {1}"`
if [ 1 -eq "$(echo "{0} <= 0 || {0} > 180" | bc)" ]; then
p=`gmt fitcircle <<! -L1 -Fs -fg |gmt math STDIN -N3 -C2 $angle ADD --IO_COL_SEPARATOR="/" =
0 0
{0} {1}
!`
else
p=`gmt fitcircle <<! -L1 -Fn -fg |gmt math STDIN -N3 -C2 $angle ADD --IO_COL_SEPARATOR="/" =
0 0
{0} {1}
!`
fi
az1=`gmt mapproject -Af{0}/{1} -fg -o2 <<< $(echo "1 0" |gmt backtracker -E$p -o0,1)`
dir_ex=`echo $az1 - {4} | bc -l`
echo "{2} {3}" |gmt backtracker -E$p -o0,1 |gmt backtracker -E{0}/{1}/$dir_ex -o0,1
"""

#Returns multiple samples within error ellipse of pole (location ({0},{1}),
#azimuth of major axis {4})                         Source: @__author__, Oct2016
ASSIGN_AZI4ROTATED_ELLIPS="""
gmt set PROJ_ELLIPSOID Sphere
agl=`gmt vector -S0/0 -TD -fg <<< "{0} {1}"`
if [ 1 -eq "$(echo "{0} <= 0 || {0} > 180" | bc)" ]; then
p=`gmt fitcircle <<! -L1 -Fs -fg |gmt math STDIN -N3 -C2 $agl ADD --IO_COL_SEPARATOR="/" =
0 0
{0} {1}
!`
else
p=`gmt fitcircle <<! -L1 -Fn -fg |gmt math STDIN -N3 -C2 $agl ADD --IO_COL_SEPARATOR="/" =
0 0
{0} {1}
!`
fi
az1=`gmt mapproject -Af{0}/{1} -fg -o2 <<< $(echo "1 0" |gmt backtracker -E$p -o0,1)`
dir_ex=`echo $az1 - {4} | bc -l`
echo "{2}" |tr -s ' '  '\n' |sed 's/\(\[\|\]\)//g;/^[[:space:]]*$/d' >/tmp/tmp1.d
echo "{3}" |tr -s ' '  '\n' |sed 's/\(\[\|\]\)//g;/^[[:space:]]*$/d' >/tmp/tmp2.d
paste /tmp/tmp1.d /tmp/tmp2.d |gmt backtracker -E$p -o0,1 |gmt backtracker -E{0}/{1}/$dir_ex -o0,1
"""

#1st, sample a whole great circle through both points ({0},{1}) and ({2},{3})
#every $accurac units; 2nd, sample half a great circle through both points
#({6},{7}) and ({4},{5}), starting from point ({6},{7}) pointing to point
#({4},{5}) by length 0.001 to 180 degrees; 3rd, find intersections of the whole
#great circle line and the half great circle line using gmt spatial
#Source: @__author__, Jan-Jun2018
INTERSECTION_BETW2DIRECTIONAL_GEODESICS="""
accurac=1E-2
accura2=1E-3
if [ ! -f /tmp/{8}.d ]; then
gmt project -C{0}/{1} -E{2}/{3} -G$accurac -L-180/180 > /tmp/{8}.d
fi
gmt project -C{6}/{7} -E{4}/{5} -G$accura2 -L$accura2/.009 > /tmp/half.gc
gmt project -C{6}/{7} -E{4}/{5} -G$accurac -L$accurac/180 >> /tmp/half.gc
gmt spatial /tmp/{8}.d /tmp/half.gc -Ie -Fl -o0,1
"""

#Source: @__author__, Jan-Jun2018
RELATIVE_LOC_INTERSECTION2NEXT_GEODESIC="""
az1=`gmt mapproject -Af{0}/{1} -fg -o2 <<< '{2} {3}'`
az2=`gmt mapproject -Af{0}/{1} -fg -o2 <<< '{4} {5}'`
a1=`gmt math -Q $az1 360 FMOD -fx --FORMAT_GEO_OUT=D =`
a2=`gmt math -Q $az2 360 FMOD -fx --FORMAT_GEO_OUT=D =`
tst=`gmt math -Q $a1 $a2 SUB ABS =`
echo $tst
"""

#Source: @__author__, Jan2018
POINT_AHEAD_GEODESIC="""
gcd=`gmt vector -S{0}/{1} -TD -fg <<< "{2} {3}" | gmt math STDIN CEIL 10 ADD =`
gmt project -C{0}/{1} -E{2}/{3} -G1 -L$gcd/`echo 1 + $gcd | bc -l` | head -n1 |gmt math STDIN -o0,1 --IO_COL_SEPARATOR="	" =
"""

def wgs2cart(dis):
    '''converts list/array of vector poles (WGS84 coordinates), in degrees, to
    array of Cartesian coordinates, in -1 <= x,y <= 1
    Source: @__author__, Apr2018'''
    ints=np.ones(len(dis)).transpose()  #get an array of ones to plug into lon,lat pairs
    _d_=np.array(dis)
    rad=np.pi/180.
    if len(_d_.shape)>1:  # array of vectors
        lon,lat=(((_d_[:,0]-180)%360)-180)*rad,_d_[:,1]*rad  #make sure -180,180
        if _d_.shape[1]==3: ints=_d_[:,2]  #take the given lengths assigned in 3rd column in _d_
    else:  # single pole(lon,lat) vector
        lon,lat=(((np.array(float(_d_[0]))-180)%360)-180)*rad,np.array(float(_d_[1]))*rad
        if len(_d_)==3: ints=np.array(_d_[2])
        else: ints=np.array([1.])
    return np.array([ints*lon/np.pi,ints*np.cos(np.pi/2.-lat)]).transpose()

def s2f(_x_):
    """convert x from str to float                  Source: @__author__, 2016"""
    try: return float(_x_)
    except ValueError: return _x_

def s2i(_x_):
    """convert x from str to integer             Source: @__author__, Jan2018"""
    try: return int(_x_)
    except ValueError: return _x_

def rm_unpaired_rows_ar(ar1,ar2,col):
    """remove unpaired row(s) in numpy array ar1 based on their common culumn(s),
    when compared to ar2. col: a number or string, e.g. 1, 'age', etc., or a
    list of numbers, e.g. [0,1,2], ['age','n']   Source: @__author__, Feb2018"""
    return ar1[np.in1d(ar1[:][col],ar2[:][col])]

def ppf(fname,pnh=1):
    """parse path file, from ASCII data to numpy array
    Source: @__author__, Feb2018"""
    pdl="\t"
    puc=(0,1,2,3,4,5,6,7)  #so far not using age range (8,9) yet
    pdt=[('dec','<f8'),('inc','<f8'),('age','<f2'),('dm','<f8'),('dp','<f8'),
         ('dm_azi','<f8'),('k','<f4'),('n','<i2'),('possib_loest_age','<f4'),
         ('possib_hiest_age','<f4')]  #f4:float32,f8:float64,i2:i:int16
    return np.genfromtxt(fname,delimiter=pdl,usecols=puc,dtype=pdt,
                         skip_header=pnh)

def ppf0(fname,pnh=0):
    """parse model-predicted path file, from ASCII data to numpy array
    Source: @__author__, Feb-Sep2018"""
    pdl="\t"
    puc=(0,1,2,3,4,5,2)
    pdt=[('dec','<f8'),('inc','<f8'),('age','<f2'),('dm_azi','<f8'),
         ('dm','<f8'),('dp','<f8'),('n','<i2')]
    vgp=np.genfromtxt(fname,delimiter=pdl,usecols=puc,dtype=pdt,skip_header=pnh)
    vgp['dm']/=2  #model-predicted pole from GMT backtracker with DIAMETER axes, but dm/dp should be radius
    vgp['dp']/=2
    vgp['n']=1
    return vgp

def ppf1(fname):
    """parse raw vgp file without weights, from ASCII data to numpy array
    Source: @__author__, Feb-Sep2018"""
    pdl="\t"
    puc=(0,1,2,4,5)
    pdt="f8,f8,i2,f8,f8"
    vgp=np.genfromtxt(fname,delimiter=pdl,usecols=puc,dtype=pdt,names=True)
    vgp['ED95']=1  #col name 'ED95' is used to store weights
    return vgp

def ppf2(fname):
    """parse raw vgp file with weights, from ASCII data to numpy array
    Source: @__author__, Dec2018"""
    pdl="\t"
    puc=(0,1,11,4,5)
    pdt="f8,f8,f8,f8,f8"
    return np.genfromtxt(fname,delimiter=pdl,usecols=puc,dtype=pdt,names=True)

def cumulative_sum4list(lst):
    """generate a list that stores the cumulative sum of integers, e.g., if my
    set of integers is [0,1,2,3], my generated list would be [0,1,3,6]
    Source: @__author__, Feb2018"""
    return list(np.cumsum(lst))

def structured_array(lst,nms,fmt):
    """lst: layers of tuples (e.g. single tuple (a,b,c)); nms: names of columns;
    fmt: dtypes of columns, e.g. int8,int16,int32,int64 can be replaced by
    equivalent string 'i1','i2','i4','i8'        Source: @__author__, Sep2018"""
    return np.array(lst,np.dtype({'names':nms,'formats':fmt}))

def ellipsenrmdev_1gen(lon,lat,azi,hma,hmi,dros=26,axis_unit=1):
    """originate random points around (0,0) on 2D plane, then rotate (0,0) with
    all these points to a specified location (lon,lat) on Earth for modeling 3D
    spherical surface ellipse; hma/hmi MUST be RADIUS (half of DIAMETER);
    azimuth's unit can be switched by axis_unit (0: have hma,hmi in kilometer;
    1: degree). Note that alpha95 derived from fisher_mean, or dp/dm, is radius,
    not diameter (Chp6, Butler98)                   Source: @__author__, 2016"""
    #sigma: population standard deviation (SD); square of sigma: population variance (e.g. v_1, v_2); length of semi-m[aj/in]or axis = 1.96 SDs
    if axis_unit==0:
        v_1,v_2=((hma/111.195051975)/1.96)**2,((hmi/111.195051975)/1.96)**2
    else: v_1,v_2=(hma/1.96)**2,(hmi/1.96)**2
    #cov: covariance matrix, which is diagonal
    pts=np.random.multivariate_normal(mean=(0,0),cov=[[v_1,0],[0,v_2]],size=dros)  #to rotate the ellipse, multiply a matrix; read fig/gaussians.pdf (http://cs229.stanford.edu/section/gaussians.pdf)
    rnd=PMAGPY3().vector_mean(np.c_[pts,np.ones(dros)])[0][:2]
    rdl=run_sh(ASSIGN_AZI4ROTATED_ELLIP.format(lon,lat,rnd[0],rnd[1],azi))  #see more info from https://pyformat.info/
    rdloc=re.split(r'\t+',rdl.decode().rstrip('\n'))
    try: return float(rdloc[0]),float(rdloc[1])
    except ValueError as err: print("error",err,"on rdloc",rdloc)  #used for debugging

def elips_nrmdev_gen_n(lon,lat,azi,hma,hmi,siz=1000,axis_unit=1):
    """see more details in function ellipsenrmdev_1gen, the difference is here
    it's generating a certain number of random points
    Source: @__author__, Oct2016"""
    if axis_unit==0:  #variance=sigma square
        v_1,v_2=((hma/111.195051975)/1.96)**2,((hmi/111.195051975)/1.96)**2
    else: v_1,v_2=(hma/1.96)**2,(hmi/1.96)**2
    rlo,rla=[],[]
    for _ in range(siz):
        pts=np.random.multivariate_normal(mean=(0,0),cov=[[v_1,0],[0,v_2]],size=26)  #be careful when size is larger than 2738, https://stackoverflow.com/questions/52587499/possible-random-multivariate-normal-bug-when-size-is-too-large
        #if siz>500:
        #    np.savetxt("/tmp/{}.txt".format(str(siz)),pts,fmt='%.9g',delimiter=',')
        #    pts=np.genfromtxt("/tmp/{}.txt".format(str(siz)),delimiter=',')
        rmd=PMAGPY3().vector_mean(np.c_[pts,np.ones(26)])[0][:2]
        rlo.append(rmd[0])
        rla.append(rmd[1])
    r_l=run_sh(ASSIGN_AZI4ROTATED_ELLIPS.format(lon,lat,rlo,rla,azi))  #see more info from https://pyformat.info/
    _r_=np.array([s.strip().split('\t') for s in r_l.decode().splitlines()])
    return list(map(float,_r_[:,0])),list(map(float,_r_[:,1]))

def get_bounds(d_i):
    """Source: Chris Rowan, 2016"""
    bounds=[]  #2sigma bounds
    cart=PMAGPY3().dir2cart(d_i).transpose()  #convert to cartesian coordinates
    mim=int(.025*len(cart[0]))
    mam=int(.975*len(cart[0]))
    for i in range(3):
        comp=cart[i]
        comp.sort()
        bounds.append([comp[mim],comp[mam]])
    return bounds

def get_bounds2d(d_i):
    """For pole vectors (WGS84 coordinates). Note that the function 'get_bounds'
    is for declination/inclination vectors; Better not used for pole vectors.
    Source: @__author__, Apr2018"""
    bounds=[]  #2sigma bounds
    cart=wgs2cart(d_i).transpose()  #convert to cartesian coordinates
    mim=int(.025*len(cart[0]))
    mam=int(.975*len(cart[0]))
    for i in range(2):
        comp=cart[i]
        comp.sort()
        bounds.append([comp[mim],comp[mam]])
    return bounds

def get_fsh(dire):
    """generate Fisher distributed points according to the supplied parameters
    (D,I,N,k) in a numpy void. Source: Chris Rowan and @__author__, 2016-2018"""
    _d_,_i_=[],[]
    for _ in range(int(dire['n'])):
        dec,inc=PMAGPY3().fshdev(dire['k'])
        drot,irot=PMAGPY3().dodirot(dec,inc,dire['dec'],dire['inc'])
        _d_.append(drot)
        _i_.append(irot)
    return np.column_stack((_d_,_i_))

def common_dir_elliptical(po1,po2,boots=1000,fn1='file1',fn2='file2'):
    """po1/2 (pole1/2 in Path1/2): numpy void; da1/2: a nested list of
    directional data [dec,inc] (a di_block). boots=1000 should be sufficient;
    Large number for boots means more computing time. See related discussion here:
    https://stats.stackexchange.com/questions/86040/rule-of-thumb-for-number-of-bootstrap-samples
    For N>25, prepare raw paleopoles beforehand in specified dir, e.g. /tmp/
    Source: @__author__ and Chris Rowan, 2016-Sep2019"""
    if po1['n']>25:
        with open('/tmp/{:s}/{:s}.txt'.format(str(fn1),str(po1['age'])),encoding='UTF-8') as _f_:
            da1=[[s2f(x) for x in line.split()] for line in _f_]
        bdi1=PMAGPY3().di_boot(da1)
    elif po1['n']<=25 and po1['n']>1:
        bdi1=[]
        for _ in range(boots):
            dir1=po1[['dec','inc','k','n']]
            fpars1=PMAGPY3().fisher_mean(get_fsh(dir1))  #need to re-think, why we need a fisher_mean step here? Just for suppressing the dispersion? But this suppressing is too much
            bdi1.append([fpars1['dec'],fpars1['inc']])
    else:
        lons1,lats1=elips_nrmdev_gen_n(po1['dec'],po1['inc'],po1['dm_azi'],
                                       po1['dm'],po1['dp'])
        bdi1=list(zip(lons1,lats1))
    if po2['n']>25:
        with open('/tmp/{:s}/{:s}.txt'.format(str(fn2),str(po2['age'])),encoding='UTF-8') as _f_:
            da2=[[s2f(x) for x in line.split()] for line in _f_]
        bdi2=PMAGPY3().di_boot(da2)
    elif po2['n']<=25 and po2['n']>1:
        bdi2=[]
        for _ in range(boots):
            dir2=po2[['dec','inc','k','n']]
            fpars2=PMAGPY3().fisher_mean(get_fsh(dir2))
            bdi2.append([fpars2['dec'],fpars2['inc']])
    else:
        lons2,lats2=elips_nrmdev_gen_n(po2['dec'],po2['inc'],po2['dm_azi'],
                                       po2['dm'],po2['dp'])
        #print(len(lons2),len(lats2))  #for debugging; should =siz in elips_nrmdev_gen_n function
        bdi2=list(zip(lons2,lats2))
        #np.savetxt("/tmp/ssd.txt",bdi2,delimiter='	',fmt='%.9g')
    ##-------********Save*To*File*For*Debugging********-------------------------
    #print(type(bdi1))
    #print(type(bdi2))
    np.savetxt("/tmp/1cloud.txt",np.array(bdi1),delimiter='	',fmt='%.9g')
    np.savetxt("/tmp/2cloud.txt",np.array(bdi2),delimiter='	',fmt='%.9g')
    ##-------********Save*To*File*For*Debugging********----------------END------
    #now check if pass or fail -pass only if error bounds overlap in x,y, and z
    bounds1,bounds2=get_bounds(bdi1),get_bounds(bdi2)  #type: list
    #**When*2*Poles*Are*Close*To*Geographic*North*Pole--z-axis-cmp-problematic--
    #print(bounds1,bounds2)  #for debugging
    if min(po1['dm'],po1['dp'])>PMAGPY3().angle((po1['dec'],po1['inc']),(0,90)):  #only for Apparent NORTH Polar Wander Paths
        bounds1[2][1]=1
    if min(po2['dm'],po2['dp'])>PMAGPY3().angle((po2['dec'],po2['inc']),(0,90)):  #only for Apparent NORTH Polar Wander Paths
        bounds2[2][1]=1
    #print(bounds1,bounds2)  #for debugging
    #**When*2*Poles*Are*Close*To*Geographic*North*Pole--------------END---------
    out=[]
    for i,j in zip(bounds1,bounds2):
        out.append(1 if i[0]>j[1] or i[1]<j[0] else 0)
    #add a tolerance for when only one axis differentiated but actually really really close to each other
    if sum(out)==1:
        out=[]
        for i,j in zip(bounds1,bounds2):
            out.append(1 if (i[0]>j[1] and abs(i[0]-j[1])>.008) or (i[1]<j[0] and abs(j[0]-i[1])>.008) else 0)  #0.008 is equivalent to about 28 arc minutes
    _o_=1 if sum(out)>=1 else 0  #1 distinguishable/separate, 0 indistinguishable/overlap
    _a_=PMAGPY3().angle((po1['dec'],po1['inc']),(po2['dec'],po2['inc']))
    return _o_,_a_[0]

def common_dir_elliptica_(po1,po2,boots=1000,fn1='file1',fn2='file2'):
    """derived from the function 'common_dir_elliptical'; the difference here is
    this function is using 'get_bounds2d' func to more accurately process the
    bounds of random pole vectors                Source: @__author__, Apr2018"""
    if po1['n']>25:
        with open('/tmp/{:s}/{:s}.txt'.format(str(fn1),str(po1['age'])),encoding='UTF-8') as _f_:
            da1=[[s2f(x) for x in line.split()] for line in _f_]
        bdi1=PMAGPY3().di_boot(da1)
    elif po1['n']<=25 and po1['n']>1:
        bdi1=[]
        for _ in range(boots):
            dir1=po1[['dec','inc','k','n']]
            fpars1=PMAGPY3().fisher_mean(get_fsh(dir1))
            bdi1.append([fpars1['dec'],fpars1['inc']])
    else:
        lons1,lats1=elips_nrmdev_gen_n(po1['dec'],po1['inc'],po1['dm_azi'],
                                       po1['dm'],po1['dp'])
        bdi1=list(zip(lons1,lats1))
    if po2['n']>25:
        with open('/tmp/{:s}/{:s}.txt'.format(str(fn2),str(po2['age'])),encoding='UTF-8') as _f_:
            da2=[[s2f(x) for x in line.split()] for line in _f_]
        bdi2=PMAGPY3().di_boot(da2)
    elif po2['n']<=25 and po2['n']>1:
        bdi2=[]
        for _ in range(boots):
            dir2=po2[['dec','inc','k','n']]
            fpars2=PMAGPY3().fisher_mean(get_fsh(dir2))
            bdi2.append([fpars2['dec'],fpars2['inc']])
    else:
        lons2,lats2=elips_nrmdev_gen_n(po2['dec'],po2['inc'],po2['dm_azi'],
                                       po2['dm'],po2['dp'])
        bdi2=list(zip(lons2,lats2))
    #now check if pass or fail -pass only if error bounds overlap in x,y, and z
    bounds1,bounds2=get_bounds2d(bdi1),get_bounds2d(bdi2)
    out=[]
    for i,j in zip(bounds1,bounds2):
        out.append(1 if i[0]>j[1] or i[1]<j[0] else 0)
    _o_=1 if sum(out)==2 else 0  #1 distinguishable, 0 indistinguishable
    _a_=PMAGPY3().angle((po1['dec'],po1['inc']),(po2['dec'],po2['inc']))
    return _o_,_a_[0]

def common_dir_ellip1gen(point,folder='traj1'):
    """Given fisher parameters of one pole, output a randome point in this
    pole's error ellipse; derived from func 'common_dir_elliptical' for numpy
    void.                                Source: @__author__, Nov2017-Feb2018"""
    if point['n']>25:
        makedirs('/tmp/{0}'.format(folder),exist_ok=True)
        with open('/tmp/{0}/{1}.txt'.format(folder,point['age']),encoding='UTF-8') as _f_:
            d_l=[[s2f(x) for x in line.split()] for line in _f_]
        bdi=PMAGPY3().di_boot(d_l,nob=1)
        lon,lat=bdi[0][0],bdi[0][1]
    elif point['n']<=25 and point['n']>1:
        _d_=point[['dec','inc','k','n']]
        fpars=PMAGPY3().fisher_mean(get_fsh(_d_))
        lon,lat=fpars['dec'],fpars['inc']
    else:
        lon,lat=ellipsenrmdev_1gen(point['dec'],point['inc'],point['dm_azi'],
                                   point['dm'],point['dp'])
    return lon,lat

def ang4_2suc_disp_direc_gdesics(lo1,la1,lo2,la2,lo3,la3):
    """Calculate direction change for 2 successive displacement directional
    geodesics that describe pole wandering (lo1,la1)->(lo2,la2) and then
    (lo2,la2)->(lo3,la3), so it is clear that they intersect at (lo2,la2)
    Source: @__author__, 2017"""
    agl=s2f(run_sh(AZI.format(lo2,la2,lo3,la3)).decode().rstrip('\n'))-\
        s2f(run_sh(AZI.format(lo2,la2,lo1,la1)).decode().rstrip('\n'))
    if agl<-180: agl,sign=360+agl,-1
    elif -180<=agl<0: agl,sign=-agl,1
    elif 0<=agl<180: sign=-1
    else: agl,sign=360-agl,1
    return (180-abs(agl))*sign

def ang_len4_1st_seg(p1x,p1y,p2x,p2y):
    """Calculate direction change and length for the first segment of APWP (a
    directional geodesic, point 1 [p1x,p1y] pointing to point 2 [p2x,p2y])
    Source: @__author__, Jan2018"""
    #segment angle change, compared to the 1st seg, here itself, so always 0
    if p1x==p2x and p1y==p2y: azi,gcd=0,0
    else:
        azi=s2f(run_sh(AZI.format(p1x,p1y,p2x,p2y)).decode().rstrip('\n'))
        gcd=s2f(PMAGPY3().angle((p1x,p1y),(p2x,p2y)))  #segment length
    return azi,gcd

def ang_len4_2nd_seg(p1x,p1y,p2x,p2y,p3x,p3y,apr):
    """Calculate orientation change and length for the 2nd segment of APWP (a
    directional geodesic, point 2 [p2x,p2y] pointing to point 3 [p3x,p3y]),
    compared to the 1st segment (point 1 [p1x,p1y] pointing to point 2
    [p2x,p2y]); also regarded as two connected directional geodesics with their
    intersection located right at point 2 [p2x,p2y]; apr is the angle change of
    its previous segment, here i.e. the 1st seg  Source: @__author__, Jan2018"""
    if p2x==p3x and p2y==p3y: agc,leh=apr,0.
    else:
        agc=ang4_2suc_disp_direc_gdesics(p1x,p1y,p2x,p2y,p3x,p3y)
        leh=s2f(PMAGPY3().angle((p2x,p2y),(p3x,p3y)))
    return agc,leh

def ang4_2sep_direc_gdesics(lo1,la1,lo2,la2,lom,lam,lon,lan,filname):
    """Calculate direction change for 2 SEPARATE directional geodesics, which
    are start(lo1,la1)->end(lo2,la2), and start(lom,lam)->end(lon,lan), so note
    that the 2 geodesics DO NOT intersect at any of these four end points. The
    key here is not only correctly determining one from the 2 intersection
    candicates, but also detecting the relative location of this intersection to
    the next directional geodesic, and further determining the 3rd point in a
    correct direction of the next geodesic.      Source: @__author__, Jan2018"""
    ise=run_sh(INTERSECTION_BETW2DIRECTIONAL_GEODESICS.format(lo1,la1,lo2,la2,
                                                              lom,lam,lon,lan,
                                                              filname))  #see more info from https://pyformat.info/
    lca=re.split(r'\n+',ise.decode("utf-8").rstrip('\n'))
    lcasiz=len(lca)
    isd=np.ones((lcasiz,),dtype=[('on12','f8'),('onmn','f8'),('lon','f8'),
                                 ('lat','f8'),('i2n','f8')])
    isdd=[]
    for j,i in enumerate(lca):
        i_i=re.split(r'\t+',i)
        lcx,lcy=s2f(i_i[0]),s2f(i_i[1])  #intersection longitude, latitude
        if abs(lcy-la2)<1E-3:
            i2n=s2f(run_sh(RELATIVE_LOC_INTERSECTION2NEXT_GEODESIC.format(lo1,la1,lcx,lcy,lo2,
                                                                          la2)).decode().rstrip('\n'))
        elif abs(lcy-la1)<1E-3:
            i2n=s2f(run_sh(RELATIVE_LOC_INTERSECTION2NEXT_GEODESIC.format(lo2,la2,lcx,lcy,lo1,
                                                                          la1)).decode().rstrip('\n'))
        else:
            i2n=s2f(run_sh(RELATIVE_LOC_INTERSECTION2NEXT_GEODESIC.format(lcx,lcy,lo1,la1,lo2,
                                                                          la2)).decode().rstrip('\n'))
        if abs(lcy-lan)<1E-3:
            i2_=s2f(run_sh(RELATIVE_LOC_INTERSECTION2NEXT_GEODESIC.format(lom,lam,lcx,lcy,lon,
                                                                          lan)).decode().rstrip('\n'))
        elif abs(lcy-lam)<1E-3:
            i2_=s2f(run_sh(RELATIVE_LOC_INTERSECTION2NEXT_GEODESIC.format(lon,lan,lcx,lcy,lom,
                                                                          lam)).decode().rstrip('\n'))
        else:
            i2_=s2f(run_sh(RELATIVE_LOC_INTERSECTION2NEXT_GEODESIC.format(lcx,lcy,lom,lam,lon,
                                                                          lan)).decode().rstrip('\n'))
        if (i2n<1 or 179<i2n<181) and (i2_<1 or 179<i2_<181):
            if i2n<1: on12=1-i2n
            else: on12=abs(180-i2n)
            if i2_<1: onmn=1-i2_
            else: onmn=abs(180-i2_)
            isd[j]=(on12,onmn,lcx,lcy,i2n)
        else: isdd.append(j)
    isd=np.delete(isd,isdd,0)
    imi=min(isd['on12'].argsort()[0],isd['onmn'].argsort()[0])  #if isd is empty, IndexError: index 0 is out of bounds for axis 0 with size 0, which normally should not happen
    lcx,lcy=isd['lon'][imi],isd['lat'][imi]  #true intersection longitude, latitude
    #in case the intersection is the same as or extremely close to the arrow point of the 1st geodesic
    if s2f(PMAGPY3().angle((lo2,la2),(lcx,lcy)))<1E-2:
        #2nd point is pole long/lat (the correct one of two intersections)
        agl=ang4_2suc_disp_direc_gdesics(lo1,la1,lcx,lcy,lon,lan)
    #determine the relative location of the intersection to the 1st geodesic
    else:
        if isd['i2n'][imi]<1: agl=ang4_2suc_disp_direc_gdesics(lo1,la1,lcx,lcy,lon,lan)
        elif isd['i2n'][imi]>179 and isd['i2n'][imi]<181:
            hd1=run_sh(POINT_AHEAD_GEODESIC.format(lo1,la1,lcx,lcy))  #must be lo1,la1
            p31=re.split(r'\t+',hd1.decode("utf-8").rstrip('\n'))
            agl=ang4_2suc_disp_direc_gdesics(s2f(p31[0]),s2f(p31[1]),lcx,lcy,
                                             lon,lan)  #(p31[0],p31[1]) is the auxiliary point
    return agl

def ang_len4ge3rd_seg(p1x,p1y,p2x,p2y,pmx,pmy,pnx,pny,apr,filname):
    """Calculate direction change and length for the third and any later segment
    of APWP (Seg No is greater than or equal to 3); apr is the direction change
    of its previous segment                      Source: @__author__, Jan2018"""
    if pmx==pnx and pmy==pny: agc,leh=apr,0.
    else:
        agc=ang4_2sep_direc_gdesics(p1x,p1y,p2x,p2y,pmx,pmy,pnx,pny,filname)
        leh=s2f(PMAGPY3().angle((pmx,pmy),(pnx,pny)))
    return agc,leh

def shape_dif_course(trj1,trj2,fmt1='textfile',fmt2='textfile',pnh1=1,pnh2=0):
    """Directional difference defined using 'course' (accumulative azimuth
    wrt the very beginning, so it could be beyond 180,360) is different from
    azimuth                                         Source: @__author__, 2016"""
    ar1=trj1 if fmt1=='ar' else ppf(trj1,pnh1)  #sep default as tab
    ar2=trj2 if fmt2=='ar' else ppf(trj2,pnh2)
    ar1,ar2=rm_unpaired_rows_ar(ar1,ar2,'age'),rm_unpaired_rows_ar(ar2,ar1,'age')  #remove age-unpaired rows in both arrays
    len1,len2,tt1,tt2,w_s,w_a=0,0,0,0,1/3,1/3  #len1/2 accumulative dif of length for trj 1/2
    accum_seg_a,accum_seg_l,accum_seg_a_dt,ma_seg_l,accum_ma_seg_l=0,0,0,0,0  #max of one section length of trj 1&2; accumulative all ma_seg_l
    n_row=len(ar2)
    lst=[]
    lst_seg_d_a,lst_seg_d_l=[],[]
    for i in range(1,n_row):
        if ar1[i-1]['dec']==ar1[i]['dec'] and ar1[i-1]['inc']==ar1[i]['inc']:
            ds1=0.
            eta1=tt1
        else:
            ds1=s2f(PMAGPY3().angle((ar1[i-1]['dec'],ar1[i-1]['inc']),
                                    (ar1[i]['dec'],ar1[i]['inc'])))
            eta1=s2f(run_sh(AZI.format(ar1[i-1]['dec'],ar1[i-1]['inc'],ar1[i]['dec'],
                                       ar1[i]['inc'])).decode().rstrip('\n'))
        if ar2[i-1]['dec']==ar2[i]['dec'] and ar2[i-1]['inc']==ar2[i]['inc']:
            ds2=0.
            eta2=tt2
        else:
            ds2=s2f(PMAGPY3().angle((ar2[i-1]['dec'],ar2[i-1]['inc']),
                                    (ar2[i]['dec'],ar2[i]['inc'])))
            eta2=s2f(run_sh(AZI.format(ar2[i-1]['dec'],ar2[i-1]['inc'],ar2[i]['dec'],
                                       ar2[i]['inc'])).decode().rstrip('\n'))
        if tt1>eta1 and tt1-eta1>180.: eta1=eta1+360.  #needs to brainstorm for a while, but now it is right
        if eta1>tt1 and eta1-tt1>180. and i>1: eta1=eta1-360.
        if tt2>eta2 and tt2-eta2>180.: eta2=eta2+360.  #i.e. eta1/2 (Course) could be greater than 360
        if eta2>tt2 and eta2-tt2>180. and i>1: eta2=eta2-360.
        #ang=abs(eta2-eta1)  #according to Course (Xie etal2003 described it is accumulation of ZhuanJiao(Rotating Angle),set CW or CCW as positive); Through re-thinking, this is a
        #good solution for closed polygons' comparison, not good for 2 trajectories
        ang=s2f(run_sh(SMALLER_ANGLE_REMAINDER.format(eta2,eta1)).decode().rstrip('\n'))
        if ang>180.: ang=360.-ang
        leh=abs(ds1-ds2)
        seg_a_dt=ang*abs(ds1)
        tt1,tt2=eta1,eta2
        lst_seg_d_a.append(ang)
        lst_seg_d_l.append(leh)
        accum_seg_a=accum_seg_a+ang  #angular difference
        accum_seg_l=accum_seg_l+leh  #length difference
        accum_seg_a_dt=accum_seg_a_dt+seg_a_dt  #Function (9) in Qi16
        len1+=ds1
        len2+=ds2
        ma_seg_l=max(ds1,ds2)
        accum_ma_seg_l+=ma_seg_l
        lst.append((i,ang,eta1,eta2,cumulative_sum4list(lst_seg_d_a)[-1],accum_seg_a_dt,leh,cumulative_sum4list(lst_seg_d_l)[-1]))
    divisor_l=accum_ma_seg_l
    divisor_a=180.*(n_row-1)
    #Qi16 functions might referred to Su15
    print(w_a*accum_seg_a/divisor_a + (1.-w_s-w_a)*accum_seg_l/divisor_l)
    print('Attn: For total dif, weights Ws & Wa are {} and {} repectively'.format(w_s,w_a))
    return structured_array(lst,['00_no','20_ang_seg_dif','22_course_seg1',
                                 '23_course_seg2','24_ang_seg_dif_accum',
                                 '25_ang_seg_dif_dt_accum','30_len_seg_dif',
                                 '34_len_seg_dif_accum'],
                            [np.uint8,np.float64,np.float64,np.float64,
                             np.float64,np.float64,np.float64,np.float64])

def shape_dif(trj1,trj2,fmt1='textfile',fmt2='textfile',whole='n',pnh1=1,pnh2=0):
    """shape difference includes both angular and length difs; Similar to the
    function 'spa_angpre_len_dif'. But here there is no significance tests on
    angular and length difference; Also angular difference value is stored in
    the row of the starting point of each segment (it's stored in the row of the
    ending point in 'spa_angpre_len_dif')        Source: @__author__, Nov2017"""
    ar1=trj1 if fmt1=='ar' else ppf(trj1,pnh1)  #sep default as tab
    ar2=trj2 if fmt2=='ar' else ppf(trj2,pnh2)
    ar1,ar2=rm_unpaired_rows_ar(ar1,ar2,'age'),rm_unpaired_rows_ar(ar2,ar1,'age')  #remove age-unpaired rows in both arrays
    w_a,w_l,tt1,tt2,len1,len2=1/2,1/2,0,0,0,0  #tt1/2 intermedium segment azimuth for trj 1/2
    accum_seg_l,accum_seg_a_dt=0,0
    n_row=min(len(ar1),len(ar2))
    lst=[]
    lst_seg_d_a=[]  #directional diff
    lst_seg_d_l=[]  #segment length diff
    for i in range(1,n_row):  #when n=0, no sense of shape cuz only a pair of poles exist
        if i==n_row-1:
            eta1,eta2,ang=0.,0.,0.
            ds1=0. if ar1[i-1]['dec']==ar1[i]['dec'] and ar1[i-1]['inc']==ar1[i]['inc'] else s2f(PMAGPY3().angle((ar1[i-1]['dec'],ar1[i-1]['inc']),(ar1[i]['dec'],ar1[i]['inc'])))
            ds2=0. if ar2[i-1]['dec']==ar2[i]['dec'] and ar2[i-1]['inc']==ar2[i]['inc'] else s2f(PMAGPY3().angle((ar2[i-1]['dec'],ar2[i-1]['inc']),(ar2[i]['dec'],ar2[i]['inc'])))
            leh=abs(ds1-ds2)
            dt_=abs(ar1[i]['age']-ar1[i-1]['age'])
            seg_a_dt=0.
            tt1,tt2=eta1,eta2
        else:
            if ar1[i-1]['dec']==ar1[i]['dec'] and ar1[i-1]['inc']==ar1[i]['inc']:
                eta1=tt1
                ds1=0.
            else:
                eta1=ang4_2suc_disp_direc_gdesics(ar1[i-1]['dec'],
                                                  ar1[i-1]['inc'],
                                                  ar1[i]['dec'],
                                                  ar1[i]['inc'],
                                                  ar1[i+1]['dec'],
                                                  ar1[i+1]['inc'])
                ds1=s2f(PMAGPY3().angle((ar1[i-1]['dec'],ar1[i-1]['inc']),
                                        (ar1[i]['dec'],ar1[i]['inc'])))
            if ar2[i-1]['dec']==ar2[i]['dec'] and ar2[i-1]['inc']==ar2[i]['inc']:
                eta2=tt2
                ds2=0.
            else:
                eta2=ang4_2suc_disp_direc_gdesics(ar2[i-1]['dec'],
                                                  ar2[i-1]['inc'],
                                                  ar2[i]['dec'],
                                                  ar2[i]['inc'],
                                                  ar2[i+1]['dec'],
                                                  ar2[i+1]['inc'])
                ds2=s2f(PMAGPY3().angle((ar2[i-1]['dec'],ar2[i-1]['inc']),
                                        (ar2[i]['dec'],ar2[i]['inc'])))
            ang=360-abs(eta2-eta1) if abs(eta2-eta1)>180 else abs(eta2-eta1)
            leh=abs(ds1-ds2)
            dt_=abs(ar1[i]['age']-ar1[i-1]['age'])
            seg_a_dt=ang*dt_
            tt1,tt2=eta1,eta2
        len1+=ds1
        len2+=ds2
        lst_seg_d_a.append(format(ang,'.7f').rstrip('0') if ang<.1 else ang)
        accum_seg_a_dt+=seg_a_dt  #similar to function (9) in Qi16
        lst_seg_d_l.append(format(leh,'.7f').rstrip('0') if leh<.1 else leh)
        accum_seg_l+=leh  #length difference
        (mean_seg_a_dt,mean_seg_l)=(0.,0.) if i==0 else (accum_seg_a_dt/abs(ar1[i]['age']-ar1[0]['age']),accum_seg_l/abs(ar1[i]['age']-ar1[0]['age']))
        divisor_l=PLATE_V_MAX_PAST/11.1195051975
        s_a,s_l=mean_seg_a_dt/POL_WAND_DIR_DIF_MAX,mean_seg_l/divisor_l
        d_shp=0. if i==0 else w_a*s_a+w_l*s_l
        lst.append((i,ar1[i]['age'],ang,eta1,eta2,
                    cumulative_sum4list(lst_seg_d_a)[-1],accum_seg_a_dt,
                    mean_seg_a_dt,leh,cumulative_sum4list(lst_seg_d_l)[-1],
                    mean_seg_l,format(d_shp,'.7f').rstrip('0') if d_shp<.1 else d_shp))
    if whole=='y': return d_shp,s_a,s_l
    else:
        print('Attn: For shape dif, weights Ws & Wl are {} and {} repectively'.format(w_a,w_l))
        return structured_array(lst,['00_no','01_tstop','20_ang_seg_dif',
                                     '22_course_seg1','23_course_seg2',
                                     '24_ang_seg_dif_accum',
                                     '25_ang_seg_dif_dt_accum',
                                     '26_ang_seg_dif_dt_mean',
                                     '30_len_seg_dif','34_len_seg_dif_accum',
                                     '35_len_seg_dif_mean','41_shape_dif'],
                                [np.uint8,'<f2','<f8','<f8','<f8','<f8','<f8',
                                 '<f8','<f8','<f8','<f8','<f8'])

def spa_angpre_len_dif(trj1,trj2,fmt1='textfile',fmt2='textfile',pnh1=1,pnh2=0,dfn1='',dfn2=''):
    """Apply sig tests seperately on per-pair-of-coeval-poles' spacial dif
    (distance) and per-pair-of-coeval-segments' angular and length difs;
    Here each segment's directional change is always relative to its previous
    segment (in fact eventually relative to the 1st segment in 2D space)
    trj1 or trj2: an ASCII text file, or a numpy array;
    fmt1: trj1 data format; fmt2: trj2 data format;
    pnh1 or pnh2: header line number, 1 means first line, 0 means no header line;
    dfn1/2: trj1/2 data file name, for creating a folder under the same number
            containing raw VGP data when N>25 during calculation
    Source: @__author__, Jan2018"""
    if fmt1=='textfile': filname1=re.split('/|\.',trj1)[-2]
    if fmt2=='textfile': filname2=re.split('/|\.',trj2)[-2]
    if fmt1=='ar': filname1=dfn1
    if fmt2=='ar': filname2=dfn2
    ar1=trj1 if fmt1=='ar' else ppf(trj1,pnh1)  #sep default as tab
    ar2=trj2 if fmt2=='ar' else ppf(trj2,pnh2)
    ar1=ar1[np.in1d(ar1[:]['age'],ar2[:]['age'])]  #remove age-unpaired rows in both arrays
    ar2=ar2[np.in1d(ar2[:]['age'],ar1[:]['age'])]
    tt1,tt2=0,0  #tt1/2 intermedium segment azimuth for trj 1/2
    n_row=min(len(ar1),len(ar2))
    lst=[]
    #print('00_no\t01_tstop\t10_spa_pol_dif\t11_spa_pol_tes\t20_ang_seg_dif\t21_ang_seg_tes\t30_len_seg_dif\t31_len_seg_tes\t22_course_seg1\t23_course_seg2\t32_len_seg1\t33_len_seg2')  #for ipynb demo
    for i in range(n_row):  # [17]
        #ind,sgd=0,0
        ind,sgd=common_dir_elliptical(ar1[i],ar2[i],fn1=filname1,fn2=filname2)
        #store Nones in the row for the 1st pole, cuz for only the 1st pole, angle change, length and their dif have no meaning except only spacial dif
        if i==0:
            eta1,eta2,ds1,ds2=np.nan,np.nan,np.nan,np.nan
            seg_d_a=np.nan
            seg_d_l=np.nan
            seg0a1=np.nan  #no specific meaning for 1st pole
            seg0l1=np.nan  #no specific meaning for 1st pole
        #store 2 angle changes, 1 ang dif, 2 lengths, 1 len dif of the 1st coeval segments in the row for the 2nd pole
        elif i==1:
            eta1,ds1=ang_len4_1st_seg(ar1[i-1]['dec'],ar1[i-1]['inc'],
                                      ar1[i]['dec'],ar1[i]['inc'])  #ds1/2 segment length for trj 1/2
            eta2,ds2=ang_len4_1st_seg(ar2[i-1]['dec'],ar2[i-1]['inc'],
                                      ar2[i]['dec'],ar2[i]['inc'])
            seg0a1=np.nan  #making the ang dif betw the 1st coeval seg pair always be 0, ie, dif not influenced by rotation models, and 2 paths don't need to be rotated into same frame
            leh=abs(ds1-ds2)  #------------------------i==1-START--------------#
            lst_d_leh_a_ras,lst_d_leh_ras_rbs=[],[]
            for _ in range(1000):  #1000
                a1x,a1y=common_dir_ellip1gen(ar1[i-1],folder=filname1)
                a2x,a2y=common_dir_ellip1gen(ar1[i],folder=filname1)
                b1x,b1y=common_dir_ellip1gen(ar2[i-1],folder=filname2)
                b2x,b2y=common_dir_ellip1gen(ar2[i],folder=filname2)
                _,ds1r=ang_len4_1st_seg(a1x,a1y,a2x,a2y)
                _,ds2r=ang_len4_1st_seg(b1x,b1y,b2x,b2y)
                lst_d_leh_a_ras.append(abs(ds1r-ds1))
                lst_d_leh_ras_rbs.append(abs(ds2r-ds1r))
            _u_=np.percentile(lst_d_leh_a_ras,97.5)
            _l_=np.percentile(lst_d_leh_ras_rbs,2.5)
            #sometimes both 0 and 1 appear for multiple 1000-runs; then e.g. run 100 1000-runs, as long as at least 5 1000-runs get 1, choose 1; because 0 should be strict to achieve.
            seg0l1=0 if _u_>_l_ else 1  #--------------------END-i==1----------#
            seg_d_a=np.nan
            seg_d_l=format(leh,'.7f').rstrip('0') if leh<.1 else leh
            tt1,tt2=eta1,eta2  #if eta1,eta2=0,0, this line is useless; kept here in case we want to measure ang dif betw the 1st coeval seg pair
        else:
            eta1,ds1=ang_len4_2nd_seg(ar1[i-2]['dec'],ar1[i-2]['inc'],
                                      ar1[i-1]['dec'],ar1[i-1]['inc'],
                                      ar1[i]['dec'],ar1[i]['inc'],tt1)
            eta2,ds2=ang_len4_2nd_seg(ar2[i-2]['dec'],ar2[i-2]['inc'],
                                      ar2[i-1]['dec'],ar2[i-1]['inc'],
                                      ar2[i]['dec'],ar2[i]['inc'],tt2)
            ang=360-abs(eta2-eta1) if abs(eta2-eta1)>180 else abs(eta2-eta1)
            leh=abs(ds1-ds2)  #------------------------i>=2-START--------------#
            lst_d_ang_a_ras,lst_d_ang_ras_rbs,lst_d_leh_a_ras,lst_d_leh_ras_rbs=[],[],[],[]
            for _ in range(1000):  #1000
                a0x,a0y=common_dir_ellip1gen(ar1[i-2],folder=filname1)
                a1x,a1y=common_dir_ellip1gen(ar1[i-1],folder=filname1)
                a2x,a2y=common_dir_ellip1gen(ar1[i],folder=filname1)
                b0x,b0y=common_dir_ellip1gen(ar2[i-2],folder=filname2)
                b1x,b1y=common_dir_ellip1gen(ar2[i-1],folder=filname2)
                b2x,b2y=common_dir_ellip1gen(ar2[i],folder=filname2)
                eta1r,ds1r=ang_len4_2nd_seg(a0x,a0y,a1x,a1y,a2x,a2y,tt1)
                eta2r,ds2r=ang_len4_2nd_seg(b0x,b0y,b1x,b1y,b2x,b2y,tt2)
                lst_d_ang_a_ras.append(abs(eta1r-eta1))
                lst_d_ang_ras_rbs.append(abs(eta2r-eta1r))
                lst_d_leh_a_ras.append(abs(ds1r-ds1))
                lst_d_leh_ras_rbs.append(abs(ds2r-ds1r))
            #####################SAVE#TrajI#VS#TrajI#########START########
            #i_i=open('/tmp/{0}i_i.d'.format(i),'w',encoding='UTF-8')
            #for j in lst_d_leh_a_ras: i_i.write("%s\n" % j)
            #i_i.close()  ########SAVE#TrajI#VS#TrajI##########END########
            #####################SAVE#TrajI#VS#TrajII########START########
            #i_ii=open('/tmp/{0}i_ii.d'.format(i),'w',encoding='UTF-8')
            #for j in lst_d_leh_ras_rbs: i_ii.write("%s\n" % j)
            #i_ii.close()  #######SAVE#TrajI#VS#TrajII#########END########
            au_=np.percentile(lst_d_ang_a_ras,97.5)
            al_=np.percentile(lst_d_ang_ras_rbs,2.5)
            #print(au_,al_)
            seg0a1=0 if au_>al_ else 1
            lu_=np.percentile(lst_d_leh_a_ras,97.5)
            ll_=np.percentile(lst_d_leh_ras_rbs,2.5)
            seg0l1=0 if lu_>ll_ else 1  #--------------------END-i>=2----------#
            seg_d_a=format(ang,'.7f').rstrip('0') if ang<.1 else ang
            seg_d_l=format(leh,'.7f').rstrip('0') if leh<.1 else leh
            tt1,tt2=eta1,eta2
        #cuz so far synchronized ages for 2 APWPs are required, so ar2[i]['age'] is also ok
        lst.append((i,ar1[i]['age'],sgd,ind,seg_d_a,seg0a1,seg_d_l,seg0l1,eta1,eta2,ds1,ds2))
        print("{0}\t{1:g}\t{2:.9g}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}".format(i,ar1[i]['age'],sgd,ind,seg_d_a,seg0a1,
                                                                                          seg_d_l,seg0l1,eta1,eta2,ds1,ds2))
    return structured_array(lst,['00_no','01_tstop','10_spa_pol_dif',
                                 '11_spa_pol_tes','20_ang_seg_dif',
                                 '21_ang_seg_tes','30_len_seg_dif',
                                 '31_len_seg_tes','22_course_seg1',
                                 '23_course_seg2',
                                 '32_len_seg1','33_len_seg2'],
                            ['<i2','<f2','<f8','<i1','<f8','<f4','<f8','<f4',
                             '<f8','<f8','<f8','<f8'])

def spa_ang1st_len_dif(trj1,trj2,fmt1='textfile',fmt2='textfile',pnh1=1,pnh2=0):
    """Apply sig tests seperately on per-pair-of-coeval-poles' spacial dif
    (distance) and per-pair-of-coeval-segments' angular and length difs; Here
    each segment's directional change is always relative to the 1st segment,
    which is more complex than "relative to its previous segment"; Requires
    statistical tests on d_a; without tests, could give unexpected d_a, e.g.
    compr of pairC&D, so deeper tweaks are needed in the future.
    Source: @__author__ and Chris Rowan, Nov2017-Feb2018"""
    if fmt1=='textfile': filname1=re.split('/|\.',trj1)[-2]
    if fmt2=='textfile': filname2=re.split('/|\.',trj2)[-2]
    ar1=trj1 if fmt1=='ar' else ppf(trj1,pnh1)  #sep default as tab
    ar2=trj2 if fmt2=='ar' else ppf(trj2,pnh2)
    ar1=ar1[np.in1d(ar1[:]['age'],ar2[:]['age'])]  #remove age-unpaired rows in both arrays
    ar2=ar2[np.in1d(ar2[:]['age'],ar1[:]['age'])]
    tt1,tt2=0,0  #tt1/2 intermedium segment azimuth for trj 1/2
    n_row=min(len(ar1),len(ar2))
    lst=[]
    print('00_no\t01_tstop\t10_spa_pol_dif\t11_spa_pol_tes\t20_ang_seg_dif\t21_ang_seg_tes\t30_len_seg_dif\t31_len_seg_tes\t22_course_seg1\t23_course_seg2\t32_len_seg1\t33_len_seg2')  #for ipynb demo
    for i in range(0,n_row):
        ind,sgd=common_dir_elliptical(ar1[i],ar2[i],fn1=filname1,fn2=filname2)
        #store Nones in the row for the 1st pole, cuz for only the 1st pole, angle change, length and their dif have no meaning except only spacial dif
        if i==0:
            eta1,eta2,ds1,ds2=np.nan,np.nan,np.nan,np.nan
            seg_d_a=np.nan  #coeval segment directional(a)/length(l) diff
            seg_d_l=np.nan
            seg0a1=np.nan  #no specific meaning for 1st pole
            seg0l1=np.nan  #no specific meaning for 1st pole
        #store 2 angle changes, 1 ang dif, 2 lengths, 1 len dif of the 1st coeval segments in the row for the 2nd pole
        elif i==1:
            eta1,ds1=ang_len4_1st_seg(ar1[i-1]['dec'],ar1[i-1]['inc'],
                                      ar1[i]['dec'],ar1[i]['inc'])  #ds1/2 segment length for trj 1/2
            eta2,ds2=ang_len4_1st_seg(ar2[i-1]['dec'],ar2[i-1]['inc'],
                                      ar2[i]['dec'],ar2[i]['inc'])
            #ang option 1: directly assigned 0 so that 1st seg as a reference
            ang=0
            #ang option 2: azi dif
            #ang=360-abs(eta2-eta1) if abs(eta2-eta1)>180 else abs(eta2-eta1)
            #ang option 3: actual dif
            #ang=ang4_2sep_direc_gdesics(ar1[i-1]['dec'],ar1[i-1]['inc'],
            #                            ar1[i]['dec'],ar1[i]['inc'],
            #                            ar2[i-1]['dec'],ar2[i-1]['inc'],
            #                            ar2[i]['dec'],ar2[i]['inc'],filname1)
            seg0a1=0  #making the ang dif betw the 1st coeval seg pair always be 0, ie, dif not influenced by rotation models, and 2 paths don't need to be rotated into same frame
            leh=abs(ds1-ds2)  #------------------------i==1-START--------------#
            lst_d_leh_a_ras,lst_d_leh_ras_rbs=[],[]
            for _ in range(1000):
                a1x,a1y=common_dir_ellip1gen(ar1[i-1],folder=filname1)
                a2x,a2y=common_dir_ellip1gen(ar1[i],folder=filname1)
                b1x,b1y=common_dir_ellip1gen(ar2[i-1],folder=filname2)
                b2x,b2y=common_dir_ellip1gen(ar2[i],folder=filname2)
                _,ds1r=ang_len4_1st_seg(a1x,a1y,a2x,a2y)
                _,ds2r=ang_len4_1st_seg(b1x,b1y,b2x,b2y)
                lst_d_leh_a_ras.append(abs(ds1r-ds1))
                lst_d_leh_ras_rbs.append(abs(ds2r-ds1r))
            _u_=np.percentile(lst_d_leh_a_ras,97.5)
            _l_=np.percentile(lst_d_leh_ras_rbs,2.5)
            seg0l1=0 if _u_>_l_ else 1  #--------------------END-i==1----------#
            seg_d_a=format(ang,'.7f').rstrip('0') if ang<.1 else ang
            seg_d_l=format(leh,'.7f').rstrip('0') if leh<.1 else leh
            tt1,tt2=eta1,eta2  #if eta1,eta2=0,0, this line is useless; kept here in case we want to measure ang dif betw the 1st coeval seg pair; if above used ang=ang4_2sep_direc_gdesics(?), then delete this line
        elif i==2:
            eta1,ds1=ang_len4_2nd_seg(ar1[i-2]['dec'],ar1[i-2]['inc'],
                                      ar1[i-1]['dec'],ar1[i-1]['inc'],
                                      ar1[i]['dec'],ar1[i]['inc'],tt1)  #if above used ang=ang4_2sep_direc_gdesics(?), then replace tt1 with 0
            eta2,ds2=ang_len4_2nd_seg(ar2[i-2]['dec'],ar2[i-2]['inc'],
                                      ar2[i-1]['dec'],ar2[i-1]['inc'],
                                      ar2[i]['dec'],ar2[i]['inc'],tt2)  #if above used ang=ang4_2sep_direc_gdesics(?), then replace tt2 with 0
            ang=360-abs(eta2-eta1) if abs(eta2-eta1)>180 else abs(eta2-eta1)
            leh=abs(ds1-ds2)  #------------------------i==2-START--------------#
            lst_d_ang_a_ras,lst_d_ang_ras_rbs,lst_d_leh_a_ras,lst_d_leh_ras_rbs=[],[],[],[]
            for _ in range(1000):
                a0x,a0y=common_dir_ellip1gen(ar1[i-2],folder=filname1)
                a1x,a1y=common_dir_ellip1gen(ar1[i-1],folder=filname1)
                a2x,a2y=common_dir_ellip1gen(ar1[i],folder=filname1)
                b0x,b0y=common_dir_ellip1gen(ar2[i-2],folder=filname2)
                b1x,b1y=common_dir_ellip1gen(ar2[i-1],folder=filname2)
                b2x,b2y=common_dir_ellip1gen(ar2[i],folder=filname2)
                eta1r,ds1r=ang_len4_2nd_seg(a0x,a0y,a1x,a1y,a2x,a2y,tt1)  #if in i==1 used ang=ang4_2sep_direc_gdesics(?), then replace tt1 with 0
                eta2r,ds2r=ang_len4_2nd_seg(b0x,b0y,b1x,b1y,b2x,b2y,tt2)  #if in i==1 used ang=ang4_2sep_direc_gdesics(?), then replace tt2 with 0
                lst_d_ang_a_ras.append(abs(eta1r-eta1))
                lst_d_ang_ras_rbs.append(abs(eta2r-eta1r))
                lst_d_leh_a_ras.append(abs(ds1r-ds1))
                lst_d_leh_ras_rbs.append(abs(ds2r-ds1r))
            au_=np.percentile(lst_d_ang_a_ras,97.5)
            al_=np.percentile(lst_d_ang_ras_rbs,2.5)
            seg0a1=0 if au_>al_ else 1  #1 distinguishable(different), 0 indistinguishable
            lu_=np.percentile(lst_d_leh_a_ras,97.5)
            ll_=np.percentile(lst_d_leh_ras_rbs,2.5)
            seg0l1=0 if lu_>ll_ else 1  #--------------------END-i==2----------#
            seg_d_a=format(ang,'.7f').rstrip('0') if ang<.1 else ang
            seg_d_l=format(leh,'.7f').rstrip('0') if leh<.1 else leh
            tt1,tt2=eta1,eta2
        else:
            if i==3 and ar1[2]['dec']==ar1[1]['dec'] and ar1[2]['inc']==ar1[1]['inc']:
                eta1,ds1=ang_len4_2nd_seg(ar1[0]['dec'],ar1[0]['inc'],
                                          ar1[1]['dec'],ar1[1]['inc'],
                                          ar1[i]['dec'],ar1[i]['inc'],tt1)  #if 2nd&3rd poles are exactly the same
            else:
                try:
                    eta1,ds1=ang_len4ge3rd_seg(ar1[0]['dec'],ar1[0]['inc'],ar1[1]['dec'],
                                               ar1[1]['inc'],ar1[i-1]['dec'],ar1[i-1]['inc'],
                                               ar1[i]['dec'],ar1[i]['inc'],tt1,filname1)
                except (UnboundLocalError,IndexError):
                    print("No intersection of trj1 {0} 1st seg {1},{2}-{3},{4} and {5}th seg {6},{7}-{8},{9} found".format(filname1,ar1[0]['dec'],ar1[0]['inc'],ar1[1]['dec'],ar1[1]['inc'],
                                                                                                                           i,ar1[i-1]['dec'],ar1[i-1]['inc'],ar1[i]['dec'],ar1[i]['inc']))  #if fail, report
                    break
            if i==3 and ar2[2]['dec']==ar2[1]['dec'] and ar2[2]['inc']==ar2[1]['inc']:
                eta2,ds2=ang_len4_2nd_seg(ar2[0]['dec'],ar2[0]['inc'],
                                          ar2[1]['dec'],ar2[1]['inc'],
                                          ar2[i]['dec'],ar2[i]['inc'],tt2)  #if 2nd&3rd poles are exactly the same
            else:
                try:
                    eta2,ds2=ang_len4ge3rd_seg(ar2[0]['dec'],ar2[0]['inc'],ar2[1]['dec'],
                                               ar2[1]['inc'],ar2[i-1]['dec'],ar2[i-1]['inc'],
                                               ar2[i]['dec'],ar2[i]['inc'],tt2,filname2)
                except (UnboundLocalError,IndexError):
                    print("No intersection of trj2 {0} 1st seg {1},{2}-{3},{4} and {5}th seg {6},{7}-{8},{9} found".format(filname2,ar2[0]['dec'],ar2[0]['inc'],ar2[1]['dec'],ar2[1]['inc'],
                                                                                                                           i,ar2[i-1]['dec'],ar2[i-1]['inc'],ar2[i]['dec'],ar2[i]['inc']))  #if fail, report
                    break
            ang=360-abs(eta2-eta1) if abs(eta2-eta1)>180 else abs(eta2-eta1)
            leh=abs(ds1-ds2)  #------------------------i>=3-START--------------#
            lst_d_ang_a_ras,lst_d_ang_ras_rbs,lst_d_leh_a_ras,lst_d_leh_ras_rbs=[],[],[],[]
            for _ in range(1000):  #for _ in range(2):  #for only wanting to see the d_angular
                a1x,a1y=common_dir_ellip1gen(ar1[i-1],folder=filname1)
                b1x,b1y=common_dir_ellip1gen(ar2[i-1],folder=filname2)
                #if fails, run through the failed iteration of the loop again
                while True:
                    try:
                        a2x,a2y=common_dir_ellip1gen(ar1[i],folder=filname1)  #resample a2x,a2y for saving time; can resample both a1(x,y) and a2
                        eta1r,ds1r=ang_len4ge3rd_seg(ar1[0]['dec'],ar1[0]['inc'],
                                                     ar1[1]['dec'],ar1[1]['inc'],
                                                     a1x,a1y,a2x,a2y,tt1,filname1)
                    except (UnboundLocalError,IndexError): continue  #if fail, reiterate again
                    break
                while True:
                    try:
                        b2x,b2y=common_dir_ellip1gen(ar2[i],folder=filname2)
                        eta2r,ds2r=ang_len4ge3rd_seg(ar2[0]['dec'],ar2[0]['inc'],
                                                     ar2[1]['dec'],ar2[1]['inc'],
                                                     b1x,b1y,b2x,b2y,tt2,filname2)
                    except (UnboundLocalError,IndexError): continue
                    break
                lst_d_ang_a_ras.append(abs(eta1r-eta1))
                lst_d_ang_ras_rbs.append(abs(eta2r-eta1r))
                lst_d_leh_a_ras.append(abs(ds1r-ds1))
                lst_d_leh_ras_rbs.append(abs(ds2r-ds1r))
            au_=np.percentile(lst_d_ang_a_ras,97.5)
            al_=np.percentile(lst_d_ang_ras_rbs,2.5)
            seg0a1=0 if au_>al_ else 1
            lu_=np.percentile(lst_d_leh_a_ras,97.5)
            ll_=np.percentile(lst_d_leh_ras_rbs,2.5)
            seg0l1=0 if lu_>ll_ else 1  #--------------------END-i>=3----------#
            seg_d_a=format(ang,'.7f').rstrip('0') if ang<.1 else ang
            seg_d_l=format(leh,'.7f').rstrip('0') if leh<.1 else leh
            tt1,tt2=eta1,eta2
        lst.append((i,ar1[i]['age'],sgd,ind,seg_d_a,seg0a1,seg_d_l,seg0l1,eta1,eta2,ds1,ds2))
        print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}".format(i,ar1[i]['age'],sgd,ind,seg_d_a,seg0a1,seg_d_l,seg0l1,eta1,eta2,ds1,ds2))
    return structured_array(lst,['00_no','01_tstop','10_spa_pol_dif',
                                 '11_spa_pol_tes','20_ang_seg_dif',
                                 '21_ang_seg_tes','30_len_seg_dif',
                                 '31_len_seg_tes','22_course_seg1',
                                 '23_course_seg2','32_len_seg1','33_len_seg2'],
                            ['<i2','<f2','<f8','<i1','<f8','<f4','<f8','<f4',
                             '<f8','<f8','<f8','<f8'])

def fit_quality(trj1,trj2,rtn_grd=1):
    """Calculate Fit Quality of a pair of polar wander paths
    trj1 or trj2: a numpy array.                 Source: @__author__, Sep2019"""
    ar1=trj1[np.in1d(trj1[:]['age'],trj2[:]['age'])]  #remove age-unpaired rows in both arrays
    ar2=trj2[np.in1d(trj2[:]['age'],trj1[:]['age'])]
    m11=(ar1['dm']+ar1['dp'])/2<=5  #a mask
    sub11=ar1['dm'][m11]
    sub11[:]=1
    m12=((ar1['dm']+ar1['dp'])/2>5)&((ar1['dm']+ar1['dp'])/2<=10)
    sub12=ar1['dm'][m12]
    sub12[:]=2
    m13=((ar1['dm']+ar1['dp'])/2>10)&((ar1['dm']+ar1['dp'])/2<=20)
    sub13=ar1['dm'][m13]
    sub13[:]=3
    m14=(ar1['dm']+ar1['dp'])/2>20
    sub14=ar1['dm'][m14]
    sub14[:]=4
    sc1=(sub11.sum()+sub12.sum()+sub13.sum()+sub14.sum())/len(ar1)
    if sc1<1.5: sy1,sy_1='A',3
    elif 1.5<=sc1<2.5: sy1,sy_1='B',2
    elif 2.5<=sc1<3.5: sy1,sy_1='C',1
    else: sy1,sy_1='D',0
    m21=(ar2['dm']+ar2['dp'])/2<=5
    sub21=ar2['dm'][m21]
    sub21[:]=1
    m22=((ar2['dm']+ar2['dp'])/2>5)&((ar2['dm']+ar2['dp'])/2<=10)
    sub22=ar2['dm'][m22]
    sub22[:]=2
    m23=((ar2['dm']+ar2['dp'])/2>10)&((ar2['dm']+ar2['dp'])/2<=20)
    sub23=ar2['dm'][m23]
    sub23[:]=3
    m24=(ar2['dm']+ar2['dp'])/2>20
    sub24=ar2['dm'][m24]
    sub24[:]=4
    sc2=(sub21.sum()+sub22.sum()+sub23.sum()+sub24.sum())/len(ar2)
    if sc2<1.5: sy2,sy_2='A',3
    elif 1.5<=sc2<2.5: sy2,sy_2='B',2
    elif 2.5<=sc2<3.5: sy2,sy_2='C',1
    else: sy2,sy_2='D',0
    if rtn_grd==0: return sc1,sc2
    elif rtn_grd==1: return sy1,sy2
    else: return sy_1,sy_2

def run_sh(script,stdin=None):
    """Raises error on non-zero return code
    Source: http://stackoverflow.com/questions/2651874/embed-bash-in-python"""
    #Note: by using a list here (['bash', ...]) we avoid quoting issues, as the
    #arguments are passed in exactly this order (spaces, quotes, and newlines
    #won't cause problems):
    with Popen(['bash','-c',script],stdout=PIPE,stderr=PIPE,stdin=PIPE) as proc:
        stdout,stderr=proc.communicate()
        if proc.returncode:
            raise ScriptException(proc.returncode,stdout,stderr,script)
        del stdin
        return stdout

class ScriptException(Exception):
    """Source: http://stackoverflow.com/questions/2651874/embed-bash-in-python
    descriptor '__init__' requires a 'Exception' object, revised by @__author__,
    Dec2017"""
    def __init__(self,returncode,stdout,stderr,script):
        self.returncode,self.stdout,self.stderr=returncode,stdout,stderr
        Exception.__init__(self,returncode,stdout,stderr,script)

def ppf3(fname,pnh=1):
    """Parses prints from function 'spa_ang1st/pre_len_dif' (ASCII data);
    Returns a numpy array                        Source: @__author__, Jun2018"""
    pdl="\t"
    puc=(0,1,2,3,4,5,6,7)
    pdt=[('00_no','<i2'),('01_tstop','<i2'),('10_spa_pol_dif','<f8'),
         ('11_spa_pol_tes','<i1'),('20_ang_seg_dif','<f8'),
         ('21_ang_seg_tes','<i1'),('30_len_seg_dif','<f8'),
         ('31_len_seg_tes','<i1')]
    return np.genfromtxt(fname,delimiter=pdl,usecols=puc,dtype=pdt,skip_header=pnh)

def calc(pair,anystr,t_0=0,t_1=530,oldform=1):
    """Prints difference measurements multiply test result 0 or 1; For t_0 != 0,
    d_a should be based on a new beginning seg
    pair: a numpy array; anystr: user string inputs; t_0: beginning point's age;
    t_1: end point's age; if oldform is 1, d_s is defined as dsicribed in the
    paper, else defined as fraction of n coeval pole pairs that are
    statistically distinguishable from each other
    Source: @__author__, Jun2018"""
    tmp=pair[(pair[:]['01_tstop']<=t_1) & (pair[:]['01_tstop']>=t_0)]
    siz=tmp[:]['10_spa_pol_dif'].size
    d_s=(tmp[:]['10_spa_pol_dif']*tmp[:]['11_spa_pol_tes']).sum()/(50*siz) if oldform==1 else tmp[:]['11_spa_pol_tes'].sum()/siz
    d_a=(tmp[2:]['20_ang_seg_dif']*tmp[2:]['21_ang_seg_tes']).sum()/(180*(siz-2))
    #30cm/yr is about 2.697961777 degree/myr
    d_l=(tmp[1:]['30_len_seg_dif']*tmp[1:]['31_len_seg_tes']).sum()/(2.697961777*(tmp[:]['01_tstop'].max()-tmp[:]['01_tstop'].min()))
    print(f'{d_s:g}\t{d_a:g}\t{d_l:g}\t{anystr:s}\t{t_0:g}\t{t_1:g}\t{(d_s+d_a+d_l)/3.:g}\t{siz:d}',
          end="\n")

def calc_nt(pair,anystr,t_0=0,t_1=530,oldform=1):
    """Derived from the function 'calc'; Prints difference measurements without
    test influencing                             Source: @__author__, Jun2018"""
    tmp=pair[(pair[:]['01_tstop']<=t_1) & (pair[:]['01_tstop']>=t_0)]
    siz=tmp[:]['10_spa_pol_dif'].size
    d_s=(tmp[:]['10_spa_pol_dif']).sum()/(50*siz) if oldform==1 else tmp[:]['11_spa_pol_tes'].sum()/siz
    d_a=tmp[2:]['20_ang_seg_dif'].sum()/(180*(siz-2))
    d_l=(tmp[1:]['30_len_seg_dif']).sum()/(2.697961777*(tmp[:]['01_tstop'].max()-tmp[:]['01_tstop'].min()))
    print(f'{d_s:g}\t{d_a:g}\t{d_l:g}\t{anystr:s}\t{t_0:g}\t{t_1:g}\t{(d_s+d_a+d_l)/3.:g}\t{siz:d}')


class PMAGPY3():
    """PmagPy (https://pmagpy.github.io/) functions below are converted to run
    correctly with Python3.6/7"""

    def __init__(self):
        pass

    @staticmethod
    def dir2cart(dis):
        '''converts list/array of vector directions, in degrees, to array of
        Cartesian coordinates, in x,y,z             Source: pmag.py of PmagPy'''
        ints=np.ones(len(dis)).transpose()  #get an array of ones to plug into dec,inc pairs
        _d_=np.array(dis)
        rad=np.pi/180.
        if len(_d_.shape)>1:  # array of vectors
            decs,incs=_d_[:,0]*rad,_d_[:,1]*rad
            if _d_.shape[1]==3: ints=_d_[:,2]  #take the given lengths assigned in 3rd column in _d_
        else:  #single vector
            decs,incs=np.array(float(_d_[0]))*rad,np.array(float(_d_[1]))*rad
            if len(_d_)==3: ints=np.array(_d_[2])
            else: ints=np.array([1.])
        return np.array([ints*np.cos(decs)*np.cos(incs),
                         ints*np.sin(decs)*np.cos(incs),
                         ints*np.sin(incs)]).transpose()

    @staticmethod
    def cart2dir(cart):
        """converts a direction to Cartesian coordinates, takes an array of
        [x,y,z])                                    Source: pmag.py of PmagPy"""
        cart=np.array(cart)
        rad=np.pi/180.  #constant to convert degrees to radians
        if len(cart.shape)>1: x_s,y_s,z_s=cart[:,0],cart[:,1],cart[:,2]
        else: x_s,y_s,z_s=cart[0],cart[1],cart[2]  #single vector
        r_s=np.sqrt(x_s**2+y_s**2+z_s**2)  #calculate resultant vector length
        decs=(np.arctan2(y_s,x_s)/rad)%360.  #calculate declination taking care of correct quadrants (arctan2) and making modulo 360.
        try: incs=np.arcsin(z_s/r_s)/rad  #calculate inclination (converting to degrees)
        except ZeroDivisionError:
            print("Error: division by 0. z_s is {0}, r_s is {1}".format(z_s,r_s))  # r_s=0
            return np.zeros(3)
        if incs==np.nan:
            print("z_s {0} > r_s {1}, arcsin value does not exist!".format(z_s,r_s))
            return np.zeros(3)
        return np.array([decs,incs,r_s]).transpose()  #return the directions list

    def fisher_mean(self,data):
        """Source: pmag.py of PmagPy
        Used enumerate instead of range+len by @__author__, 2017"""
        #calculates fisher parameters for data
        xbar,fpars=[0,0,0],{}
        _n_=len(data)  #origianl code line, mainly for many directions, i.e. dec,inc
        if _n_<2: return fpars
        _x_=self.dir2cart(data)
        #for rec in data: _x_.append(dir2cart([rec[0],rec[1],1.]))
        for i,_ in enumerate(_x_):
            for j in range(3): xbar[j]+=_x_[i][j]
        _r_=np.sqrt(sum(xbar[c]**2 for c in range(3)))
        for k in range(3): xbar[k]/=_r_
        dire=self.cart2dir(xbar)
        fpars["dec"],fpars["inc"],fpars["n"],fpars["r"]=dire[0],dire[1],_n_,_r_
        if abs(_n_-_r_)>1E-8:
            k=(_n_-1.)/(_n_-_r_)
            fpars["k"]=k
            csd=81./np.sqrt(k)
        else: fpars['k'],csd=100000,0.  #fpars['k']='inf'
        _b_=20.**(1./(_n_-1.))-1
        _a_=1-_b_*(_n_-_r_)/_r_
        #fpars["a"]=_a_  #check _a_
        if _a_<-1: _a_=-1
        a95=np.arccos(_a_)*180./np.pi
        fpars["alpha95"],fpars["csd"]=a95,csd  #estimated angular standard deviation, CSD
        if _a_<0: fpars["alpha95"]=180.  #alpha95 is radius, not diameter
        return fpars

    def angle(self,d_1,d_2):
        """call to angle(d_1,d_2) returns array of angles between lists of 2
        directions d_1,d_2 where d_1 is, for example,
        [[Dec1,Inc1],[Dec2,Inc2],etc.]                 Source: pmag.py of PmagPy
        Modified a bit in style by @__author__, 2017"""
        d_1,d_2=np.array(d_1),np.array(d_2)
        d_1=d_1[:,0:2] if len(d_1.shape)>1 else d_1[:2]  #strip off intensity
        d_2=d_2[:,0:2] if len(d_2.shape)>1 else d_2[:2]  #strip off intensity
        x_1,x_2=self.dir2cart(d_1),self.dir2cart(d_2)  #convert to cartesian from polar
        angles=[]  #set up a list for angles
        for k in range(x_1.shape[0]):  #single vector
            angles.append((np.arccos(np.dot(x_1[k],x_2[k]))*180./np.pi)%360.)  #take the dot product
        return np.array(angles)  #returns in arc degree, not radian

    @staticmethod
    def pseudo(dis):
        """draw a bootstrap sample of DIrectionS    Source: pmag.py of PmagPy"""
        siz=len(dis)
        return np.array(dis)[np.random.randint(siz,size=siz)]

    def di_boot(self,dis,nob=5000):
        """returns bootstrap parameters for Directional data     Source: pmag.py
        of PmagPy                Modified a bit in style by @__author__, 2017"""
        #fpars=self.fisher_mean(dis)  #get average DI for whole dataset
        #now do bootstrap to collect d_i bootstrap means
        d_i=[]  # list of bootstrap directions
        for _ in range(nob):  #repeat nob (number of bootstraps) times
            #if k%50==0: print(k,' out of ',nob)
            pdis=self.pseudo(dis)  #get a pseudosample
            bfpars=self.fisher_mean(pdis)  #get bootstrap mean bootstrap sample
            d_i.append([bfpars['dec'],bfpars['inc']])
        return d_i

    def dogeo(self,dec,inc,azi,plg):
        """rotates dec,inc into geographic coordinates using azi,plg as azimuth
        and plunge of _x_ direction                    Source: pmag.py of PmagPy
        Modified a bit in style by @__author__, 2017"""
        dir_=[dec,inc,1.]  #put dec inc in direction list and set length to unity
        _x_=self.dir2cart(dir_)  #get cartesian coordinates
        # set up rotation matrix
        a_1,a_2,a_3=self.dir2cart([azi,plg,1.]),self.dir2cart([azi+90.,0,1.]),self.dir2cart([azi-180.,90.-plg,1.])
        # do rotation
        xp_=a_1[0]*_x_[0]+a_2[0]*_x_[1]+a_3[0]*_x_[2]
        yp_=a_1[1]*_x_[0]+a_2[1]*_x_[1]+a_3[1]*_x_[2]
        zp_=a_1[2]*_x_[0]+a_2[2]*_x_[1]+a_3[2]*_x_[2]
        # transform back to dec,inc
        dir_geo=self.cart2dir([xp_,yp_,zp_])
        return dir_geo[0],dir_geo[1]    # send back declination and inclination

    def dodirot(self,dec,inc,dbar,ibar):
        """dec=declination,inc=inclination, dbar/ibar are the desired mean
        direction; Returns the rotated Dec/Inc pair
        Source: pmag.py of PmagPy"""
        _d_,irot=self.dogeo(dec,inc,dbar,90.-ibar)
        drot=_d_-180.
        #drot,irot=dogeo(dec,inc,Dbar,Ibar)
        if drot<360.: drot=drot+360.
        if drot>360.: drot=drot-360.
        return drot,irot

    @staticmethod
    def fshdev(kap):
        """kap is kappa, returns a direction from distribution with mean
        declination of 0, inclination of 90 and kappa of kap
        Source: pmag.py of PmagPy"""
        r_1,r_2=random(),random()
        _l_=np.exp(-2*kap)
        fac=np.sqrt((-np.log(r_1*(1-_l_)+_l_))/(2*kap))
        return 2*np.pi*r_2*180./np.pi,90.-2*np.arcsin(fac)*180./np.pi

    def vector_mean(self,data):
        """calculates the vector mean of a given set of vectors;
        Source: pmag.py of PmagPy
        Used enumerate instead of range+len by @__author__, 2017"""
        xbar,_x_=[0,0,0],[]
        for rec in data: _x_.append(self.dir2cart(rec))
        for i,_ in enumerate(_x_):
            for j in range(3): xbar[j]+=_x_[i][j]
        _r_=np.sqrt(sum(xbar[c]**2 for c in range(3)))
        for k in range(3): xbar[k]/=_r_
        dire=self.cart2dir(xbar)
        return dire,_r_

def main1():
    """Run the algorithm on real-world examples of pmag paths (at reduced data
    density) vs. modeled one"""
    plid=101
    tbin=10
    step=5
    #-------------Prepare Model Predicted APWP----------------------------------
    modl_pp=ppf0('/home/i/Desktop/git/public/making_of_reliable_APWPs/data/{}FHS120predictPWP{}{}.d'.format(plid,tbin,step))
    modl_pp[:]['dm']/=111.195051975
    modl_pp[:]['dp']/=111.195051975  #--------------------------------END-------

    #-NA APWPs from Different Algorithms, versus FHS Model Predicted APWP-------
    modl='ay18'
    pid='{}comb'.format(plid)
    wer='/mnt/g/tmp/Fay18_101comb_10_5_1_2/_4/New'
    for fod in range(57,200):
        for mav in [1]:  #[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]:
            for wgt in [2]:  #[0,1,2,3,4,5]:
                pmag_pp=ppf('{0}/{1:03d}/{2}_{3}/{2}_{3}_{4}_{5}_{6}_{7}.txt'.format(wer,fod,modl,pid,tbin,
                                                                                     step,mav,wgt),
                            pnh=1)
                for i in pmag_pp:
                    #print(i['age'])
                    if i['n']>25:
                        makedirs('/tmp/{0}_{1}_{2}_{3}_{4}_{5}'.format(modl,pid,tbin,
                                                                       step,mav,wgt),
                                 exist_ok=True)
                        raw_dir1='{0}/{1:03d}/{2}/{3}/{2}_{3}_{4}_{5}_{6}_{7}/{8}_{9}.d'.format(wer,fod,modl,pid,tbin,step,mav,wgt,
                                                                                                round(i['age']+step,1),round(i['age']-step,1))
                        raw_dir2='{0}/{1:03d}/{2}/{3}/{2}_{3}_{4}_{5}_{6}_{7}/{8}_{9}_WT.d'.format(wer,fod,modl,pid,tbin,step,mav,wgt,
                                                                                                   round(i['age']+step,1),round(i['age']-step,1))
                        raw_pls=ppf1(raw_dir1) if path.isfile(raw_dir1) else ppf2(raw_dir2)
                        #print(raw_pls)
                        np.savetxt("/tmp/{}_{}_{}_{}_{}_{}/{}.txt".format(modl,pid,str(tbin),str(step),str(mav),str(wgt),str(i['age'])),
                                   raw_pls,delimiter='	',fmt='%.9g')
                print('-----------{}----DOUBLE-CHECK----{}----------'.format(mav,wgt))
                makedirs('{0}/{1:03d}/{2}_{3}/{4}_{5}_simil'.format(wer,fod,modl,pid,tbin,step),
                         exist_ok=True)
                #print(pmag_pp)
                #print(modl_pp)
                simil=spa_angpre_len_dif(pmag_pp,modl_pp,'ar','ar',dfn1='{0}_{1}_{2}_{3}_{4}_{5}'.format(modl,pid,tbin,step,mav,wgt))
                np.savetxt('{0}/{1:03d}/{2}_{3}/{4}_{5}_simil/{2}_{3}_{4}_{5}_{6}_{7}.d'.format(wer,fod,modl,pid,str(tbin),str(step),str(mav),str(wgt)),
                           simil,delimiter='	',fmt='%.9g',comments='',
                           header="00_no	01_tstop	10_spa_pol_dif	11_spa_pol_tes	20_ang_seg_dif	21_ang_seg_tes	30_len_seg_dif	31_len_seg_tes	22_course_seg1	23_course_seg2	32_len_seg1	33_len_seg2")
                fq_=fit_quality(pmag_pp,modl_pp)
                print('---FQ:{}---{}----DOUBLE-CHECK----{}---END----'.format(fq_,mav,wgt))

def main():
    """Run the algorithm on real-world examples of pmag paths vs. modeled one"""
    plid=501
    tbin=2
    step=1
    #-------------Prepare Model Predicted APWP----------------------------------
    modl_pp=ppf0('/home/g/Desktop/git/public/making_of_reliable_APWPs/data/{}FHS120predictPWP{}{}.d'.format(plid,tbin,step))
    modl_pp[:]['dm']/=111.195051975
    modl_pp[:]['dp']/=111.195051975  #--------------------------------END-------

    #-Pmag APWPs from Different Algorithms, versus HS Model Predicted APWP------
    modl='ay18'
    pid='{}comb'.format(plid)
    wer='/tmp'
    for mav in range(1,2):  #[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]:
        for wgt in range(4,6):  #[0,1,2,3,4,5]:
            pmag_pp=ppf('{0}/{1}_{2}/{1}_{2}_{3}_{4}_{5}_{6}.txt'.format(wer,modl,pid,tbin,
                                                                         step,mav,wgt),
                        pnh=1)
            for i in pmag_pp:
                #print(i['age'])
                if i['n']>25:
                    makedirs('/tmp/{0}_{1}_{2}_{3}_{4}_{5}'.format(modl,pid,tbin,
                                                                   step,mav,wgt),
                             exist_ok=True)
                    raw_dir1='{0}/{1}/{2}/{1}_{2}_{3}_{4}_{5}_{6}/{7}_{8}.d'.format(wer,modl,pid,tbin,step,mav,wgt,
                                                                                    round(i['age']+step,1),round(i['age']-step,1))
                    raw_dir2='{0}/{1}/{2}/{1}_{2}_{3}_{4}_{5}_{6}/{7}_{8}_WT.d'.format(wer,modl,pid,tbin,step,mav,wgt,
                                                                                       round(i['age']+step,1),round(i['age']-step,1))
                    raw_pls=ppf1(raw_dir1) if path.isfile(raw_dir1) else ppf2(raw_dir2)
                    #print(raw_pls)
                    np.savetxt('/tmp/{}_{}_{}_{}_{}_{}/{}.txt'.format(modl,pid,str(tbin),str(step),str(mav),str(wgt),str(i['age'])),
                               raw_pls,delimiter='	',fmt='%.9g')
            #print('-----------{}----DOUBLE-CHECK----{}----------'.format(mav,wgt))
            makedirs('{0}/{1}_{2}/{3}_{4}_simil'.format(wer,modl,pid,tbin,step),
                     exist_ok=True)
            #print(pmag_pp)
            #print(modl_pp)
            simil=spa_angpre_len_dif(pmag_pp,modl_pp,'ar','ar',dfn1='{0}_{1}_{2}_{3}_{4}_{5}'.format(modl,pid,tbin,step,mav,wgt))
            np.savetxt('{0}/{1}_{2}/{3}_{4}_simil/{1}_{2}_{3}_{4}_{5}_{6}.d'.format(wer,modl,pid,str(tbin),str(step),str(mav),str(wgt)),
                       simil,delimiter='	',fmt='%.9g',comments='',
                       header="00_no	01_tstop	10_spa_pol_dif	11_spa_pol_tes	20_ang_seg_dif	21_ang_seg_tes	30_len_seg_dif	31_len_seg_tes	22_course_seg1	23_course_seg2	32_len_seg1	33_len_seg2")
            #fq_=fit_quality(pmag_pp,modl_pp)
            #print('---FQ:{}---{}----DOUBLE-CHECK----{}---END----'.format(fq_,mav,wgt))

def main_fq():
    """Calculate Fit Quality scores for the 168 pairs of pmag path vs. modeled
    path"""
    plid=501
    tbin=10
    step=5
    #-------------Prepare Model Predicted APWP----------------------------------
    modl_pp=ppf0('/home/g/Desktop/git/public/making_of_reliable_APWPs/data/{}FHS120predictPWP{}{}.d'.format(plid,tbin,step))
    modl_pp[:]['dm']/=111.195051975
    modl_pp[:]['dp']/=111.195051975  #--------------------------------END-------

    #-Pmag APWPs from Different Algorithms, versus HS Model Predicted APWP------
    modl='ay18'
    pid='{}comb'.format(plid)
    wer='/tmp'
    for mav in range(28):  #[17]:
        for wgt in range(6):  #[4]:
            pmag_pp=ppf('{0}/{1}_{2}/{1}_{2}_{3}_{4}_{5}_{6}.txt'.format(wer,modl,pid,tbin,
                                                                         step,mav,wgt),
                        pnh=1)
            fq_=fit_quality(pmag_pp,modl_pp)
            print('{}\t{}\t{}\t{}'.format(mav,wgt,fq_[0],fq_[1]))

def test():
    """Test if the functions work"""
    #Example APWP 1
    pwp1="""> 1. time window length: 10 m.y.; step: 5.0 m.y.
80.14199557582	84.96247253771	0.	3.915449	3.915449	-88.3113304535	28.01	49	0.	7.
80.42135537758	84.95277471805	5.	3.886849	3.886849	-87.67128574119	26.28	53	0.	9.
33.38843294652	86.88873057981	10.	15.183887	15.183887	232.941877527	14.26	8	0.	20.
355.39679874686	82.70569183131	15.	13.310493	13.310493	205.627887959	12.73	11	9.	25.
308.98665645271	77.0051507776	20.	13.980808	13.980808	163.147828636	6.74	19	12.	27.
274.60831938936	84.77268712744	25.	11.846867	11.846867	116.878107816	7.83	22	2.	56.
157.75717380674	83.16102644444	30.	5.154164	5.154164	-10.5390047879	45.95	18	2.	65.
199.00176522396	81.48747812327	35.	7.02344	7.02344	26.745796008	48.27	10	2.	65.
186.54455769797	82.51841548799	40.	5.367151	5.367151	11.627064047	92.98	9	15.	65.
200.72373633575	79.69498729973	45.	12.463375	12.463375	23.294640395	14.38	11	34.	56.
199.64756368863	81.11023939111	50.	12.284287	12.284287	26.704068203	14.78	11	35.	56.
155.77817403083	79.88648221927	55.	5.699119	5.699119	-14.8857382382	72.81	10	40.	72.
141.77335287233	80.25809073133	60.	7.607781	7.607781	-27.8924228568	24.48	16	40.	72.
169.59568379585	78.1602830698	65.	7.646845	7.646845	-2.404197301	25.97	15	45.	85.
232.46368881052	81.25270600104	70.	11.643872	11.643872	68.378030809	16.34	11	56.	84.
232.20995172143	78.99518546502	75.	10.40819	10.40819	68.556969641	15.55	14	50.	100.
196.59319393588	72.15855230297	80.	14.521823	14.521823	20.786364683	10.85	11	65.	100.
165.17925380398	59.98672310777	85.	22.438845	22.438845	-17.9793429138	31.24	3	65.	100.
194.57748128034	71.48924647882	100.	10.272066	10.272066	-2.6994187522	80.97	4	50.	150.
183.38270826627	66.99981750984	105.	7.825817	7.825817	-6.9778981134	74.25	6	60.	150.
188.77797283482	74.84849853192	110.	12.694145	12.694145	13.061898626	23.56	7	65.	146.
179.70207466514	84.32495998057	115.	16.944904	16.944904	0.6564554017	30.37	4	100.	123.
194.5452608166	70.45888170695	120.	3.315562	3.315562	-19.1052823603	213.26	10	100.	146.
194.99486400335	71.47385732389	125.	2.508328	2.508328	-17.080972447	488.66	8	100.	146.
182.79517239359	65.83240478355	130.	15.039353	15.039353	-38.909980903	38.29	4	100.	160.
198.35928957379	66.30096825434	135.	23.931109	23.931109	-28.6682630383	111.01	2	125.	140.
203.1000061	58.	140.	3.8	3.8	-19.2411071267	550.	1	138.	146.
"""
    #Example APWP 2
    pwp2="""-179	90	0	8.99321e-05	8.99321e-05	0	0	0	0	0
-174.427839776	89.6139557807	5	1.59806	0.890964	115.520357868	0	0	0	0
-175.123593158	89.2733285766	10	3.15048	1.6682	121.876635273	0	0	0	0
-154.679734239	88.4257673069	15	2.51776	1.3208	146.615864191	0	0	0	0
-157.327747143	87.4544636774	20	3.09328	1.79648	147.822057651	0	0	0	0
-171.05977315	86.6942827357	25	4.92677	1.86982	143.815966436	0	0	0	0
173.12464727	86.2293292854	30	3.39455	2.56919	127.064961091	0	0	0	0
-170.381724003	85.3636547838	35	4.22046	2.08658	135.388037234	0	0	0	0
165.843502228	85.5098622356	40	8.69388	3.89947	102.669343822	0	0	0	0
-164.018491324	83.3247134078	45	5.96303	2.55664	140.49103191	0	0	0	0
-154.310824855	80.5602801087	50	7.93168	3.5791	143.466274302	0	0	0	0
-158.050597628	79.7302132383	55	6.76195	4.1792	135.949203446	0	0	0	0
-158.204140077	78.228301634	60	9.64978	6.2734	120.059389548	0	0	0	0
-157.4434009	77.468353925	65	6.59471	3.90179	140.681537736	0	0	0	0
-168.404415378	78.8343512249	70	9.16825	5.23025	127.201198244	0	0	0	0
-159.456958728	75.9292342122	75	5.24019	2.77688	143.129239623	0	0	0	0
-176.264371811	78.2350501922	80	8.66873	5.04817	117.665818285	0	0	0	0
-171.615832678	75.187957123	85	4.75738	2.38559	131.493143696	0	0	0	0
-169.755372271	75.1780263693	90	7.32862	2.69115	135.148011206	0	0	0	0
174.166497258	74.5207842853	95	4.62704	3.12618	116.468706884	0	0	0	0
159.462801164	72.5547669906	100	5.86798	5.2042	61.8779217729	0	0	0	0
164.484202345	73.652846477	105	4.641	3.23742	85.2590846134	0	0	0	0
164.346545735	70.8703145084	110	8.07616	4.17059	89.0351101171	0	0	0	0
141.907723703	74.5711373906	115	8.11677	4.33071	53.6664237377	0	0	0	0
138.370516992	73.4335877076	120	14.4025	7.30031	45.9082446966	0	0	0	0
125.186971113	72.3245139848	125	4.6285	1.20312	41.6310240315	0	0	0	0
123.526143186	71.6343221466	130	5.96749	1.1516	38.4366740563	0	0	0	0
122.113254332	71.0169053138	135	4.60901	1.1039	39.6420966533	0	0	0	0
120.580076057	70.3853096532	140	6.48623	1.44101	38.8357093934	0	0	0	0
"""
    #Example raw paleopoles for N>25, since for apwp1's 0 and 5 Ma, their N.25
    pwp1_0ma="""170.3999939	73.30000305	1	2.5
82.5	79.69999695	1	.5
158.	81.	1	1
88.69999695	73.	1	.5
117.8000031	83.19999695	1	1.5
324.7000122	86.80000305	1	1
266.8999939	86.09999847	1	.5
23.89999962	37.20000076	1	1.5
96.5	70.5	1	.5
91.59999847	65.59999847	1	.5
317.7000122	81.59999847	1	.5
82.	80.80000305	1	.5
237.8000031	73.40000153	1	.5
70.30000305	70.90000153	1	.5
351.	67.	1	.5
35.79999924	81.90000153	1	2
72.80000305	74.69999695	1	.5
192.3000031	72.59999847	1	.5
83.5	71.80000305	1	.5
36.	87.19999695	1	.5
83.	83.	1	1
63.79999924	87.80000305	1	1
93.90000153	85.80000305	1	4
166.	87.	1	1
311.	47.	1	1
192.	85.	1	1
170.	84.	1	1
261.	83.	1	1
200.	80.	1	1
39.	83.	1	1
128.	78.	1	.5
93.09999847	71.40000153	1	.5
54.40000153	76.59999847	1	.5
228.	80.	1	.5
133.	63.	1	.5
31.39999962	71.19999695	1	.5
323.	69.40000153	1	1
62.40000153	84.30000305	1	1
98.90000153	83.09999847	1	3.5
123.6999969	77.80000305	1	1.5
27.5	68.30000305	1	1.5
116.1999969	79.19999695	1	4
50.	87.	1	3.5
300.	86.	1	1.5
134.1999969	85.19999695	1	2.5
82.69999695	86.19999695	1	1.5
119.3000031	78.80000305	1	3.5
136.1999969	69.30000305	1	1.5
142.8000031	83.0	1	2
"""

    pwp1_5ma="""170.3999939	73.30000305	1	2.5
354.2000122	75.59999847	1	5.5
4.099999905	81.	1	7.5
4.5	80.5	1	8.5
82.5	79.69999695	1	.5
147.6999969	51.	1	7.5
158.0	81.	1	1
88.69999695	73.	1	.5
117.8000031	83.19999695	1	1.5
324.7000122	86.80000305	1	1
266.8999939	86.09999847	1	.5
23.89999962	37.20000076	1	1.5
96.5	70.5	1	.5
91.59999847	65.59999847	1	.5
317.7000122	81.59999847	1	.5
82.0	80.80000305	1	.5
237.8000031	73.40000153	1	.5
70.30000305	70.90000153	1	.5
351.	67.	1	.5
35.79999924	81.90000153	1	2
72.80000305	74.69999695	1	.5
192.3000031	72.59999847	1	.5
83.5	71.80000305	1	.5
36.	87.19999695	1	.5
83.	83.	1	1
63.79999924	87.80000305	1	1
93.90000153	85.80000305	1	4
166.	87.	1	1
311.	47.	1	1
192.	85.	1	1
170.	84.	1	1
261.	83.	1	1
200.	80.	1	1
39.	83.	1	1
128.	78.	1	0.5
93.09999847	71.40000153	1	.5
54.40000153	76.59999847	1	.5
228.	80.	1	.5
133.	63.	1	.5
31.39999962	71.19999695	1	.5
323.	69.40000153	1	1
62.40000153	84.30000305	1	1
98.90000153	83.09999847	1	3.5
123.6999969	77.80000305	1	1.5
27.5	68.30000305	1	1.5
116.1999969	79.19999695	1	4
50.	87.	1	3.5
300.	86.	1	1.5
134.1999969	85.19999695	1	2.5
82.69999695	86.19999695	1	1.5
119.3000031	78.80000305	1	3.5
136.1999969	69.30000305	1	1.5
142.8000031	83.	1	2
"""

    with open('/tmp/1.d','w',encoding='UTF-8') as apwp1:
        apwp1.write(pwp1)
        apwp1.close()
    with open('/tmp/2.d','w',encoding='UTF-8') as apwp2:
        apwp2.write(pwp2)
        apwp2.close()
    makedirs('/tmp/1',exist_ok=True)
    with open('/tmp/1/0.txt','w',encoding='UTF-8') as apwp1_0:
        apwp1_0.write(pwp1_0ma)
        apwp1_0.close()
    with open('/tmp/1/5.txt','w',encoding='UTF-8') as apwp1_5:
        apwp1_5.write(pwp1_5ma)
        apwp1_5.close()

    diff=spa_angpre_len_dif('/tmp/1.d','/tmp/2.d')
    print(diff)

if __name__=="__main__": main()
