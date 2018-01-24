#!/usr/bin/python3

'''--------------------------------------------------------------------------###
Created on 5May2016
Modified on 23Jan2018

@__author__	:	Chenjian Fu
@__email__		:	cfu3@kent.edu
@__purpose__	:	To quantitatively compare paleomagnetic APWPs
@__version__	:	0.4.4
@__license__	:	GNU General Public License v3.0

Spherical Path Comparison (spComparison) Package is developed for quantitatively
measuring similarity of spherical paths, especially the paleomagnetic apparent
polar wander paths (APWPs) of tectonic plates. It is powered by GMT
(http://gmt.soest.hawaii.edu/) and PmagPy (https://pmagpy.github.io/).
Copyright (C) 2016-2018 @__author__

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.
--------------------------------------------------------------------------------
Environment:
    Python3.6 + NumPy + Pandas + (Numba, only if using a NVIDIA graphics card)
    GMT + *NIX(-like) Shell                     (PmagPy installation not needed)
--------------------------------------------------------------------------------
TODO:
    1. Add functions for per-segment sig tests in function "ang_len_dif";
       edit the corresponding in "apwp_dif_with_ang_len_sig_tests"
    2. Tidy functions up into classes
###--------------------------------------------------------------------------'''

import random,subprocess,os,re,uuid, pandas as pd, numpy as np
from numba import jit  #to accelerate python codes

PLATE_V_MAX_PAST=30  #according to Swanson-Hysell etal.2009, Kulakov etal.2014; today it's about 15.44cm/yr, DeMets etal.2010
POL_WAND_DIR_DIF_MAX=180

#Source: @__author__, Oct2015
AZI="""
gmt mapproject -Af{0}/{1} -fg -o2 <<< '{2} {3}'
"""

#Source: @__author__, Jan2016
SMALLER_ANGLE_REMAINDER="""
gmt math -Q {0} {1} SUB 360 FMOD ABS =
"""

#Source: @__author__, Oct2016
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
dir_ex=`gmt math -Q $az1 {4} SUB =`
echo "{2} {3}" |gmt backtracker -E$p -o0,1 |gmt backtracker -E{0}/{1}/$dir_ex -o0,1
"""

#Source: @__author__, Oct2016
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
dir_ex=`gmt math -Q $az1 {4} SUB =`
echo "{2}" |tr -s ' '  '\n' |sed 's/\(\[\|\]\)//g;/^[[:space:]]*$/d' >/tmp/tmp1.d
echo "{3}" |tr -s ' '  '\n' |sed 's/\(\[\|\]\)//g;/^[[:space:]]*$/d' >/tmp/tmp2.d
paste /tmp/tmp1.d /tmp/tmp2.d |gmt backtracker -E$p -o0,1 |gmt backtracker -E{0}/{1}/$dir_ex -o0,1
"""

#Source: @__author__, Jan2018
INTERSECTION_BETW2DIRECTIONAL_GEODESICS="""
accurac=1E-2
gmt project -C{0}/{1} -E{2}/{3} -G$accurac -L$accurac/180 > /tmp/half_gc.d
gmt project -C{4}/{5} -E{6}/{7} -G$accurac -L-180/`gmt math -Q 180 $accurac SUB =` > /tmp/gc.d
gmt spatial /tmp/half_gc.d /tmp/gc.d -Ie -Fl |head -n1 |gmtmath STDIN -o0,1 --IO_COL_SEPARATOR="	" =
"""

#Source: @__author__, Jan2018
INTERSECTION_BETW2DIRECTIONAL_GEODESI2S="""
gmt spatial /tmp/half_gc.d /tmp/gc.d -Ie -Fl |sed '/^\s*$/d' |sed '2q;d' |gmtmath STDIN -o0,1 --IO_COL_SEPARATOR="	" =
"""

#Source: @__author__, Jan2018
INTERSECTION_BETW2DIRECTIONAL_GEODESI3S="""
gmt spatial /tmp/half_gc.d /tmp/gc.d -Ie -Fl |sed '/^\s*$/d' |sed '3q;d' |gmtmath STDIN -o0,1 --IO_COL_SEPARATOR="	" =
"""

#Source: @__author__, Jan2018
RELATIVE_LOC_INTERSECTION2NEXT_GEODESIC="""
az1=`gmt mapproject -Af{0}/{1} -fg -o2 <<< '{2} {3}'`
az2=`gmt mapproject -Af{0}/{1} -fg -o2 <<< '{4} {5}'`
a1=`gmt math -Q $az1 360 FMOD -fx --FORMAT_GEO_OUT=D =`
a2=`gmt math -Q $az2 360 FMOD -fx --FORMAT_GEO_OUT=D =`
tst=`gmt math -Q $a1 $a2 SUB ABS =`
if [ $(gmt math -Q $tst 1 LT =) -eq 1 ]; then echo 0
elif [ $(gmt math -Q $tst 179 GT =) -eq 1 ]; then echo 1
else echo $tst
fi
"""

#Source: @__author__, Jan2018
POINT_AHEAD_GEODESIC="""
gcd=`gmt vector -S{0}/{1} -TD -fg <<< "{2} {3}" | gmt math STDIN CEIL 10 ADD =`
gmt project -C{0}/{1} -E{2}/{3} -G1 -L$gcd/`gmt math -Q 1 $gcd ADD =` | head -n1 |gmtmath STDIN -o0,1 --IO_COL_SEPARATOR="	" =
"""

def s2f(_x_):
    """convert x from str to float                  Source: @__author__, 2016"""
    try: return float(_x_)
    except ValueError: return _x_

def s2i(_x_):
    """convert x from str to integer             Source: @__author__, Jan2018"""
    try: return int(_x_)
    except ValueError: return _x_

def txt2df_awk(txt,_h_=None,_c_='>'):
    """convert text to pandas dataframe             Source: @__author__, 2017"""
    return pd.read_table(txt,header=_h_,comment=_c_,dtype=str).apply(pd.to_numeric)

def ellipsenrmdev_1gen(lon,lat,azi,maj,mio,dros=26,axis_unit=1):
    """2D plane originning at (0,0); then rotate all points to any specific
    location on the Earth surface for modeling 3D spherical surface ellipse;
    random source; maj/mio must be diameter, not radius, if radius, don't forget
    multiply by 2 beforehand (https://en.wikipedia.org/wiki/68%E2%80%9395%E2%80%9399.7_rule)
    azimuth is in degree(s). When axis_unit=0/1, maj,mio in kilometers/degrees.
    Note that alpha95 from fisher_mean is radius, not diameter (Chp6, Butler98)
    Source: @__author__, 2016"""
    #sigma: standard deviation; variance1/2=square of sigma1/2; length of semi-maj(in)or axis = 1.96 standard deviations
    if axis_unit==0:
        v_1,v_2=((maj/111.195051975)/1.96)**2,((mio/111.195051975)/1.96)**2
    else: v_1,v_2=(maj/1.96)**2,(mio/1.96)**2
    #cov: covariance; the covariance matrix is diagonal
    pts=np.random.multivariate_normal(mean=(0,0),cov=[[v_1,0],[0,v_2]],size=dros) #to rotate the ellipse, multiply a matrix; read fig/gaussians.pdf (http://cs229.stanford.edu/section/gaussians.pdf)
    rnd=PMAGPY36().vector_mean(np.c_[pts,np.ones(dros)])[0][:2]
    rdl=run_sh(ASSIGN_AZI4ROTATED_ELLIP.format(lon,lat,rnd[0],rnd[1],azi))  #see more info from https://pyformat.info/
    rdloc=re.split(r'\t+',rdl.decode().rstrip('\n'))
    try: return float(rdloc[0]),float(rdloc[1])
    except ValueError as err: print("error",err,"on rdloc",rdloc)	#used for debugging

def elips_nrmdev_gen_n(lon,lat,azi,maj,mio,dros=1000,axis_unit=1):
    """when axis_unit=0/1, maj,mio are in kilometers/degrees.
    alpha95 derived from fisher_mean is radius, no diameter (Chp6, Butler98)
    Source: @__author__, Oct2016"""
    if axis_unit==0:	#variance=sigma square
        v_1,v_2=((maj/111.195051975)/1.96)**2,((mio/111.195051975)/1.96)**2
    else: v_1,v_2=(maj/1.96)**2,(mio/1.96)**2
    pts=np.random.multivariate_normal(mean=(0,0),cov=[[v_1,0],[0,v_2]],size=dros)
    r_l=run_sh(ASSIGN_AZI4ROTATED_ELLIPS.format(lon,lat,pts[:,0],pts[:,1],azi))  #see more info from https://pyformat.info/
    _r_=np.array([s.strip().split('\t') for s in r_l.decode().splitlines()])
    return list(map(float,_r_[:,0])),list(map(float,_r_[:,1]))

@jit(nopython=False,parallel=True)
def get_bounds(d_i):
    """Source: Chris Rowan, 2016"""
    bounds=[] #2sigma bounds
    cart=PMAGPY36().dir2cart(d_i).transpose() #convert to cartesian coordinates
    mim=int(.025*len(cart[0]))
    mam=int(.975*len(cart[0]))
    for i in range(3):
        comp=cart[i]
        comp.sort()
        bounds.append([comp[mim],comp[mam]])
    return bounds

def get_fish(dire):
    """generate fisher distributed points according to the supplied parameters
    (includes D,I,N,k) in a pandas object.          Source: Chris Rowan, 2016"""
    _d_,_i_=[],[]
    for _ in range(int(dire.n)):
        dec,inc=PMAGPY36().fshdev(dire.k)
        drot,irot=PMAGPY36().dodirot(dec,inc,dire.dec,dire.inc)
        _d_.append(drot)
        _i_.append(irot)
    return np.column_stack((_d_,_i_))

def common_dir_elliptical(po1,po2,dros=1000,bootsteps=5000,fn1='file1',fn2='file2'):
    """da1/2 : a nested list of directional data [dec,inc] (a di_block)
    Note that bootsteps=1000 could be insufficient for when 1<N<=25, at least
    5000 is needed to make sure the result is >95% significantly robust through
    tests. But at the same time, it needs much more computing time. See further
    discussion here:
    https://stats.stackexchange.com/questions/86040/rule-of-thumb-for-number-of-bootstrap-samples
    Source: @__author__ and Chris Rowan, 2016"""
    #-----pandas DataFrame for Path 1, denoted as po1---------------------
    #for N>25, prepare raw paleopoles beforehand in specified dir, e.g. /tmp/
    if po1.N>25:
        with open('/tmp/{:s}/{:s}.txt'.format(str(fn1),str(int(po1.Age)))) as _f_:
            da1=[[s2f(x) for x in line.split()] for line in _f_]
        bdi1=PMAGPY36().di_boot(da1)
    elif po1.N<=25 and po1.N>1:
        bdi1=[]
        for _ in range(bootsteps):
            dir1=po1[['Lon','Lat','kappa','N']]
            dir1.rename(index={"Lon":"dec","Lat":"inc","kappa":"k","N":"n"},
                        inplace=True)
            fpars1=PMAGPY36().fisher_mean(get_fish(dir1))
            bdi1.append([fpars1['dec'],fpars1['inc']])
    else:
        lons1,lats1=elips_nrmdev_gen_n(po1.Lon,po1.Lat,po1.Az,
                                       po1.Major,po1.Minor,dros)
        bdi1=np.column_stack((lons1,lats1))
    #-----pandas DataFrame for Path 2, denoted as po2---------------------
    if po2.N>25:
        with open('/tmp/{:s}/{:s}.txt'.format(str(fn2),str(int(po2.Age)))) as _f_:
            da2=[[ s2f(x) for x in line.split()] for line in _f_]
        bdi2=PMAGPY36().di_boot(da2)
    elif po2.N<=25 and po2.N>1:
        bdi2=[]
        for _ in range(bootsteps):
            dir2=po2[['Lon','Lat','kappa','N']]
            dir2.rename(index={"Lon":"dec","Lat":"inc","kappa":"k","N":"n"},
                        inplace=True)
            fpars2=PMAGPY36().fisher_mean(get_fish(dir2))
            bdi2.append([fpars2['dec'],fpars2['inc']])
    else:
        lons2,lats2=elips_nrmdev_gen_n(po2.Lon,po2.Lat,po2.Az,
                                       po2.Major,po2.Minor,dros)
        bdi2=np.column_stack((lons2,lats2))
    #now check if pass or fail -pass only if error bounds overlap in x,y, and z
    bounds1,bounds2=get_bounds(bdi1),get_bounds(bdi2)
    out=[]
    for i,j in zip(bounds1,bounds2):
        out.append(0 if i[0]>j[1] or i[1]<j[0] else 1)
    _o_=1. if sum(out)==3 else 0.  #1 'Pass' 0 'Fail'
    _a_=PMAGPY36().angle((po1['Lon'],po1['Lat']),(po2['Lon'],po2['Lat']))
    #return pd.Series([_o_,_a_[0]],index=['Outcome','angle'])
    return _o_,_a_[0]

def common_dir_elliptical_1gen(point,folder='traj1'):
    """derived from func 'common_dir_elliptical' Source: @__author__, Nov2017"""
    if point.N>25:
        os.makedirs(os.path.join('/tmp',folder),exist_ok=True)
        dat=os.path.join('/tmp',folder,'{:s}.txt')
        with open(dat.format(str(int(point.Age)))) as _f_:
            d_l=[[s2f(x) for x in line.split()] for line in _f_]
        bdi=PMAGPY36().di_boot(d_l,nob=1)
        lon,lat=bdi[0][0],bdi[0][1]
    elif point.N<=25 and point.N>1:
        _d_=point[['Lon','Lat','kappa','N']]
        _d_.rename(index={"Lon":"dec","Lat":"inc","kappa":"k","N":"n"},
                   inplace=True)
        fpars=PMAGPY36().fisher_mean(get_fish(_d_))
        lon,lat=fpars['dec'],fpars['inc']
    else:
        lon,lat=ellipsenrmdev_1gen(point.Lon,point.Lat,point.Az,
                                   point.Major,point.Minor)
    return lon,lat

def btr(trj1,trj2,fmt1='textfile',fmt2='textfile',fn1='file1',fn2='file2'):
    """Source: @__author__, 2016-2017"""
#-------Method-Of-Removing-Lines-Start-With-">"
#     f1=open(trj1,'r')
#     lines1=f1.readlines()
#     f1.close()
#     f1=open("/tmp/1.txt","w")
#     for line in lines1:
#         if not line.startswith(">"): f1.write(line)
#     f1.close()
#     f2=open(trj2,'r')
#     lines2=f2.readlines()
#     f2.close()
#     f2=open("/tmp/2.txt","w")
#     for line in lines2:
#         if not line.startswith(">"): f2.write(line)
#     f2.close()
#---END-Method-Of-Removing-Lines-Start-With-">"
    df1=trj1 if fmt1=='df' else txt2df_awk(trj1)  #sep default as tab
    df2=trj2 if fmt2=='df' else txt2df_awk(trj2)
    df1,df2=df1[df1[2].isin(df2[2])],df2[df2[2].isin(df1[2])]  #remove age-unpaired rows in both dataframes
    df1.to_csv('/tmp/1.txt',sep='\t',
               encoding='utf-8',header=False,index=False,float_format='%g')
    df2.to_csv('/tmp/2.txt',sep='\t',
               encoding='utf-8',header=False,index=False,float_format='%g')
    i,lst_no,lst_t,frq_l,ang_l=0,[],[],[],[]
    for li1,li2 in zip(open("/tmp/1.txt",'r'),open("/tmp/2.txt",'r')):
        if not li1.startswith(">") and not li2.startswith(">"):
            lst_no.append(i)
            i=i+1
            #split each line on space to get records
            rc1=li1.split('\t') if '\t' in li1 else li1.split()
            d_1,i_1,t_1,ma1,mi1,az1,k_1,n_1=s2f(rc1[0]),s2f(rc1[1]),s2f(rc1[2]),\
                                            s2f(rc1[3]),s2f(rc1[4]),s2f(rc1[5]),\
                                            s2f(rc1[6]),int(rc1[7])
            rc2=li2.split('\t') if '\t' in li2 else li2.split()
            d_2,i_2,t_2,ma2,mi2,az2,k_2,n_2=s2f(rc2[0]),s2f(rc2[1]),s2f(rc1[2]),\
                                            s2f(rc2[3]),s2f(rc2[4]),s2f(rc2[5]),\
                                            s2f(rc2[6]),int(rc2[7])
            if t_1-t_2==0:
                pl1=pd.Series([d_1,i_1,t_1,ma1,mi1,az1,k_1,n_1,'green'],
                              index=['Lon','Lat','Age','Major','Minor',
                                     'Az','kappa','N','PointColor'])
                pl2=pd.Series([d_2,i_2,t_2,ma2,mi2,az2,k_2,n_2,'orange'],
                              index=['Lon','Lat','Age','Major','Minor',
                                     'Az','kappa','N','PointColor'])
                outcom,angl=common_dir_elliptical(pl1,pl2,fn1=fn1,fn2=fn2)
                frq_l.append(outcom)
                ang_l.append(angl)
                lst_t.append(t_1)
    #return 1.-sum(frq_l)/i
    return pd.DataFrame({'31_spa_pol_dif':frq_l,'32_spa_pol_ang':ang_l,
                         '90_no':lst_no,'94_tstop':lst_t})

@jit(nopython=False,parallel=True)
def btr_(trj,trj2,fmt1='textfile',fmt2='textfile',fna='file'):
    """derived from func 'btr'                   Source: @__author__, Nov2017"""
    df_=trj if fmt1=='df' else txt2df_awk(trj)  #sep default as tab
    df2=trj2 if fmt2=='df' else txt2df_awk(trj2)
    df_=df_[df_[2].isin(df2[2])]  #remove age-unpaired rows in this dataframe
    df_.to_csv('/tmp/1.txt',sep='\t',
               encoding='utf-8',header=False,index=False,float_format='%g')
    o_lo,o_la,_t_=[],[],[]
    for lin in open("/tmp/1.txt",'r'):
        if not lin.startswith(">"):
            #split each line on space to get records
            rec=lin.split('\t') if '\t' in lin else lin.split()
            dec,inc,tim,maj,mio,azi,kap,_n_=s2f(rec[0]),s2f(rec[1]),s2f(rec[2]),\
                                            s2f(rec[3]),s2f(rec[4]),s2f(rec[5]),\
                                            s2f(rec[6]),int(rec[7])
            pl1=pd.Series([dec,inc,tim,maj,mio,azi,kap,_n_,'green'],
                          index=['Lon','Lat','Age','Major','Minor',
                                 'Az','kappa','N','PointColor'])
            lon,lat=common_dir_elliptical_1gen(pl1,folder=fna)
            o_lo.append(lon)
            o_la.append(lat)
            _t_.append(tim)
    return pd.DataFrame({0:o_lo,1:o_la,2:_t_,3:0,4:0,5:0,6:1E8,7:0,8:_t_,9:_t_})  #'dec','inc','age','dm','dp','dm_azi','k','size','possib_loest_age','possib_hiest_age'

def ang4_2suc_disp_direc_gdesics(lo1,la1,lo2,la2,lo3,la3):
    """Angle change calculation for 2 successive displacement directional
    geodesics that describe pole wandering (lo1,la1)->(lo2,la2) and then
    (lo2,la2)->(lo3,la3), so it is clear that they intersect at (lo2,la2)
    Source: @__author__, 2017"""
    agl=s2f(run_sh(AZI.format(lo2,la2,lo3,la3)).decode().rstrip('\n'))-\
        s2f(run_sh(AZI.format(lo2,la2,lo1,la1)).decode().rstrip('\n'))
    if agl<-180: agl,sign=360+agl,-1
    elif agl>=-180 and agl<0: agl,sign=-agl,1
    elif agl>=0 and agl<180: sign=-1
    else: agl,sign=360-agl,1
    return (180-abs(agl))*sign

def ang_len4_1st_seg(p1x,p1y,p2x,p2y):
    """Angle change and length calculation for the first segment of APWP (a
    directional geodesic, point 1 [p1x,p1y] pointing to point 2 [p2x,p2y])
    Source: @__author__, Jan2018"""
    #segment angle change, compared to the 1st seg, here itself, so always 0
    ds1=0. if p1x==p2x and p1y==p2y else s2f(PMAGPY36().angle((p1x,p1y),
                                                              (p2x,p2y)))  #segment length
    return 0.,ds1

def ang_len4_2nd_seg(p1x,p1y,p2x,p2y,p3x,p3y,apr):
    """Angle change and length calculation for the second segment of APWP (a
    directional geodesic, point 2 [p2x,p2y] pointing to point 3 [p3x,p3y]),
    compared to the first segment (point 1 [p1x,p1y] pointing to point 2
    [p2x,p2y]); also regarded as two connected directional geodesics with their
    intersection located right at point 2 [p2x,p2y]; apr is the angle change of
    its previous segment, here i.e. the 1st seg  Source: @__author__, Jan2018"""
    if p2x==p3x and p2y==p3y: agc,leh=apr,0.
    else:
        agc=ang4_2suc_disp_direc_gdesics(p1x,p1y,p2x,p2y,p3x,p3y)
        leh=s2f(PMAGPY36().angle((p2x,p2y),(p3x,p3y)))
    return agc,leh

def ang4_2sep_direc_gdesics(lo1,la1,lo2,la2,lom,lam,lon,lan):
    """Angle change calculation for 2 SEPARATE directional geodesics, which are
    (lo1,la1)->(lo2,la2) and (lom,lam)->(lon,lan), so note that the two
    geodesics DO NOT intersect at any of these four end points.
    The key point here is not only correctly determining the right one from the
    two intersection candicates, but also detecting the relative location of
    this intersection to the next directional geodesic, and further determining
    the 3rd point in the right direction of the next geodesic.
    Source: @__author__, Jan2018"""
    ise=run_sh(INTERSECTION_BETW2DIRECTIONAL_GEODESICS.format(lo1,la1,lo2,la2,
                                                              lom,lam,lon,lan))  #see more info from https://pyformat.info/
    lca=re.split(r'\t+',ise.decode("utf-8").rstrip('\n'))
    lcx,lcy=s2f(lca[0]),s2f(lca[1])
    i2n=s2i(run_sh(RELATIVE_LOC_INTERSECTION2NEXT_GEODESIC.format(lom,lam,lon,
                                                                  lan,lcx,
                                                                  lcy)).decode().rstrip('\n'))
    if i2n not in (0, 1):
        ise=run_sh(INTERSECTION_BETW2DIRECTIONAL_GEODESI2S.format(lo1,la1,lo2,la2,
                                                                  lom,lam,lon,lan))
        lca=re.split(r'\t+',ise.decode("utf-8").rstrip('\n'))
        lcx,lcy=s2f(lca[0]),s2f(lca[1])
        i2n=s2i(run_sh(RELATIVE_LOC_INTERSECTION2NEXT_GEODESIC.format(lom,lam,lon,
                                                                      lan,lcx,
                                                                      lcy)).decode().rstrip('\n'))
    if i2n not in (0, 1):
        ise=run_sh(INTERSECTION_BETW2DIRECTIONAL_GEODESI3S.format(lo1,la1,lo2,la2,
                                                                  lom,lam,lon,lan))
        lca=re.split(r'\t+',ise.decode("utf-8").rstrip('\n'))
        lcx,lcy=s2f(lca[0]),s2f(lca[1])
        i2n=s2i(run_sh(RELATIVE_LOC_INTERSECTION2NEXT_GEODESIC.format(lom,lam,lon,
                                                                      lan,lcx,
                                                                      lcy)).decode().rstrip('\n'))
    if i2n not in (0, 1):
        print('Angle is{0}, Points are {1},{2},{3},{4},{5},{6},{7},{8}'.format(i2n,lo1,la1,lo2,la2,lom,lam,lon,lan))
    #in case the intersection is the same as or extremely close to the starting point of the next geodesic
    if s2f(PMAGPY36().angle((lom,lam),(lcx,lcy)))<1E-2:
        #2nd point is pole long/lat (the right one of two intersections)
        agl=ang4_2suc_disp_direc_gdesics(lo1,la1,lcx,lcy,lon,lan)
    #determine the relative location of the intersection to the next geodesic
    else:
        if i2n==0:
            hd1=run_sh(POINT_AHEAD_GEODESIC.format(lom,lam,lcx,lcy))
            p31=re.split(r'\t+',hd1.decode("utf-8").rstrip('\n'))
            agl=ang4_2suc_disp_direc_gdesics(lo1,la1,lcx,lcy,s2f(p31[0]),
                                             s2f(p31[1]))
        #2nd point is pole long/lat (the right one of two intersections)
        elif i2n==1: agl=ang4_2suc_disp_direc_gdesics(lo1,la1,lcx,lcy,lom,lam)
        else: print("Angle betw CD(next segment)&CI(starting point of next seg pointing to intersection) {0} ShouldBe 0°/180°. Check func RELATIVE_LOC_INTERSECTION2NEXT_GEODESIC.".format(i2n))
    return agl

def ang_len4ge3rd_seg(p1x,p1y,p2x,p2y,pmx,pmy,pnx,pny,apr):
    """Angle change and length calculation for the third and any later segment
    of APWP (Seg No is greater than or equal to 3); apr is the angle change of
    its previous segment                         Source: @__author__, Jan2018"""
    if pmx==pnx and pmy==pny: agc,leh=apr,0.
    else:
        agc=ang4_2sep_direc_gdesics(p1x,p1y,p2x,p2y,pmx,pmy,pnx,pny)
        leh=s2f(PMAGPY36().angle((pmx,pmy),(pnx,pny)))
    return agc,leh

def shape_dif(trj1,trj2,fmt1='textfile',fmt2='textfile',whole='n'):
    """Derived from function 'apwp_dif'; shape dif includes both angular and
    length difs                                  Source: @__author__, Nov2017"""
    df1=trj1 if fmt1=='df' else txt2df_awk(trj1)  #sep default as tab
    df2=trj2 if fmt2=='df' else txt2df_awk(trj2)
    df1,df2=df1[df1[2].isin(df2[2])],df2[df2[2].isin(df1[2])]  #remove age-unpaired rows in both dataframes
    w_a,w_l,tt1,tt2,len1,len2=1/2,1/2,0,0,0,0  #tt1/2 intermedium segment azimuth for trj 1/2; len1/2 total length of trj 1/2
    accum_seg_a,accum_seg_l,accum_seg_a_dt=0,0,0
    n_row=min(len(df1.index),len(df2.index))
    lst_seg_a,lst_accum_seg_a,lst_accum_seg_a_dt,lst_mean_seg_a_dt=[],[],[],[]  #directional diff
    lst_seg_l,lst_accum_seg_l,lst_mean_seg_l=[],[],[]  #segment length diff
    lst_no,lst_t,lst_eta1,lst_eta2,d_shp_l=[],[],[],[],[]
    for i in range(0,n_row):
        if i==0:
            eta1,eta2,ang,ds1,ds2,leh=0.,0.,0.,0.,0.,0.  #ds1/2 intermedium segment GCD for trj 1/2
            dt_,seg_a_dt=abs(df1.iloc[i+1][2]-df1.iloc[i][2]),0.  #what it is doesn't matter
        elif i==n_row-1:
            eta1,eta2,ang=0.,0.,0.
            ds1=0. if df1.iloc[i-1][0]==df1.iloc[i][0] and df1.iloc[i-1][1]==df1.iloc[i][1] else s2f(PMAGPY36().angle((df1.iloc[i-1][0],df1.iloc[i-1][1]),(df1.iloc[i][0],df1.iloc[i][1])))
            ds2=0. if df2.iloc[i-1][0]==df2.iloc[i][0] and df2.iloc[i-1][1]==df2.iloc[i][1] else s2f(PMAGPY36().angle((df2.iloc[i-1][0],df2.iloc[i-1][1]),(df2.iloc[i][0],df2.iloc[i][1])))
            leh=abs(ds1-ds2)
            dt_=abs(df1.iloc[i][2]-df1.iloc[i-1][2])
            seg_a_dt=0.
            tt1,tt2=eta1,eta2
        else:
            if df1.iloc[i-1][0]==df1.iloc[i][0] and df1.iloc[i-1][1]==df1.iloc[i][1]:
                eta1=tt1
                ds1=0.
            else:
                eta1=ang4_2suc_disp_direc_gdesics(df1.iloc[i-1][0],
                                                  df1.iloc[i-1][1],
                                                  df1.iloc[i][0],
                                                  df1.iloc[i][1],
                                                  df1.iloc[i+1][0],
                                                  df1.iloc[i+1][1])
                ds1=s2f(PMAGPY36().angle((df1.iloc[i-1][0],df1.iloc[i-1][1]),
                                         (df1.iloc[i][0],df1.iloc[i][1])))
            if df2.iloc[i-1][0]==df2.iloc[i][0] and df2.iloc[i-1][1]==df2.iloc[i][1]:
                eta2=tt2
                ds2=0.
            else:
                eta2=ang4_2suc_disp_direc_gdesics(df2.iloc[i-1][0],
                                                  df2.iloc[i-1][1],
                                                  df2.iloc[i][0],
                                                  df2.iloc[i][1],
                                                  df2.iloc[i+1][0],
                                                  df2.iloc[i+1][1])
                ds2=s2f(PMAGPY36().angle((df2.iloc[i-1][0],df2.iloc[i-1][1]),
                                         (df2.iloc[i][0],df2.iloc[i][1])))
            ang=360-abs(eta2-eta1) if abs(eta2-eta1)>180 else abs(eta2-eta1)
            leh=abs(ds1-ds2)
            dt_=abs(df1.iloc[i][2]-df1.iloc[i-1][2])
            seg_a_dt=ang*dt_
            tt1,tt2=eta1,eta2
        lst_eta1.append(eta1)
        lst_eta2.append(eta2)
        len1+=ds1
        len2+=ds2
        lst_no.append(i)
        lst_t.append(df1.iloc[i][2])  #because that df1&2 ages are synchronized is required here, so df2.iloc[i][2] is also ok
        lst_seg_a.append(format(ang,'.7f').rstrip('0') if ang<.1 else ang)
        accum_seg_a+=ang  #angular difference
        accum_seg_a_dt+=seg_a_dt  #similar to function (9) in Qi16
        lst_seg_l.append(format(leh,'.7f').rstrip('0') if leh<.1 else leh)
        accum_seg_l+=leh  #length difference
        (mean_seg_a_dt,mean_seg_l)=(0.,0.) if i==0 else (accum_seg_a_dt/abs(df1.iloc[i][2]-df1.iloc[0][2]),accum_seg_l/abs(df1.iloc[i][2]-df1.iloc[0][2]))
        lst_accum_seg_a.append(accum_seg_a)
        lst_accum_seg_a_dt.append(accum_seg_a_dt)
        lst_mean_seg_a_dt.append(mean_seg_a_dt)
        lst_accum_seg_l.append(accum_seg_l)
        lst_mean_seg_l.append(mean_seg_l)
        divisor_l=PLATE_V_MAX_PAST/11.1195051975  #i.e. about 2.7 degree/myr, magnitude of velocity
        s_a,s_l=mean_seg_a_dt/POL_WAND_DIR_DIF_MAX,mean_seg_l/divisor_l
        d_shp=0. if i==0 else w_a*s_a+w_l*s_l
        d_shp_l.append(format(d_shp,'.7f').rstrip('0') if d_shp<.1 else d_shp)
    if whole=='y': return d_shp,s_a,s_l
    else:
        print('Attn: For shape dif, weights Ws & Wl are {} and {} repectively'.format(w_a,w_l))
        return pd.DataFrame({'01_shape_dif':d_shp_l,'11_ang_seg_dif':lst_seg_a,
                             '12_ang_seg_dif_accum':lst_accum_seg_a,
                             '13_ang_seg_dif_dt_accum':lst_accum_seg_a_dt,
                             '14_ang_seg_dif_dt_mean':lst_mean_seg_a_dt,
                             '21_len_seg_dif':lst_seg_l,
                             '22_len_seg_dif_accum':lst_accum_seg_l,
                             '23_len_seg_dif_mean':lst_mean_seg_l,
                             '33_spa_dif':0,'90_no':lst_no,
                             '91_course1':lst_eta1,'92_course2':lst_eta2,
                             '94_tstop':lst_t})

def ang_len_dif(trj1,trj2,fmt1='textfile',fmt2='textfile',whole='n'):
    """Derived from function 'shape_dif'; applying sig tests seperately on
    per-segment's angular and length difs; angular difference was originally
    just azimuth difference, which here is modified to always relative to the
    1st segment                                  Source: @__author__, Jan2018"""
    df1=trj1 if fmt1=='df' else txt2df_awk(trj1)  #sep default as tab
    df2=trj2 if fmt2=='df' else txt2df_awk(trj2)
    df1,df2=df1[df1[2].isin(df2[2])],df2[df2[2].isin(df1[2])]  #remove age-unpaired rows in both dataframes
    tt1,tt2,len1,len2=0,0,0,0  #tt1/2 intermedium segment azimuth for trj 1/2; len1/2 total length of trj 1/2
    accum_seg_a,accum_seg_l,accum_seg_a_dt=0,0,0
    n_row=min(len(df1.index),len(df2.index))
    lst_seg_a,lst_accum_seg_a,lst_accum_seg_a_dt,lst_mean_seg_a_dt=[],[],[],[]  #directional diff
    lst_seg_l,lst_accum_seg_l,lst_mean_seg_l=[],[],[]  #segment length diff
    lst_no,lst_t,lst_eta1,lst_eta2=[],[],[],[]
    for i in range(0,n_row):
        #store 0s in the row for the 1st pole, cuz for the 1st pole, only 1, angle, length and their dif have no meaning except only spacial dif
        if i==0:
            eta1,eta2,ang,ds1,ds2,leh=0.,0.,0.,0.,0.,0.  #ds1/2 intermedium segment GCD for trj 1/2
            dt_,seg_a_dt=abs(df1.iloc[i+1][2]-df1.iloc[i][2]),0.  #what it is doesn't matter
        #store angle change, ang dif, length, len dif of the 1st segment in the row for the 2nd pole
        elif i==1:
            eta1,ds1=ang_len4_1st_seg(df1.iloc[i-1][0],df1.iloc[i-1][1],
                                      df1.iloc[i][0],df1.iloc[i][1])
            eta2,ds2=ang_len4_1st_seg(df2.iloc[i-1][0],df2.iloc[i-1][1],
                                      df2.iloc[i][0],df2.iloc[i][1])
            ang=0 #making the ang dif betw the 1st coeval seg pair always be 0, ie, dif not influenced by rotation models, and 2 paths don't need to be rotated into same frame
            leh=abs(ds1-ds2)
            dt_=abs(df1.iloc[i][2]-df1.iloc[i-1][2])
            seg_a_dt=ang*dt_
            tt1,tt2=eta1,eta2  #if eta1,eta2=0,0, this line is useless; kept here in case we want to measure ang dif betw the 1st coeval seg pair
        elif i==2:
            eta1,ds1=ang_len4_2nd_seg(df1.iloc[i-2][0],df1.iloc[i-2][1],
                                      df1.iloc[i-1][0],df1.iloc[i-1][1],
                                      df1.iloc[i][0],df1.iloc[i][1],tt1)
            eta2,ds2=ang_len4_2nd_seg(df2.iloc[i-2][0],df2.iloc[i-2][1],
                                      df2.iloc[i-1][0],df2.iloc[i-1][1],
                                      df2.iloc[i][0],df2.iloc[i][1],tt2)
            ang=360-abs(eta2-eta1) if abs(eta2-eta1)>180 else abs(eta2-eta1)
            leh=abs(ds1-ds2)
            dt_=abs(df1.iloc[i][2]-df1.iloc[i-1][2])
            seg_a_dt=ang*dt_
            tt1,tt2=eta1,eta2
        else:
            eta1,ds1=ang_len4ge3rd_seg(df1.iloc[0][0],df1.iloc[0][1],
                                       df1.iloc[1][0],df1.iloc[1][1],
                                       df1.iloc[i-1][0],df1.iloc[i-1][1],
                                       df1.iloc[i][0],df1.iloc[i][1],tt1)
            eta2,ds2=ang_len4ge3rd_seg(df2.iloc[0][0],df2.iloc[0][1],
                                       df2.iloc[1][0],df2.iloc[1][1],
                                       df2.iloc[i-1][0],df2.iloc[i-1][1],
                                       df2.iloc[i][0],df2.iloc[i][1],tt2)
            ang=360-abs(eta2-eta1) if abs(eta2-eta1)>180 else abs(eta2-eta1)
            leh=abs(ds1-ds2)
            dt_=abs(df1.iloc[i][2]-df1.iloc[i-1][2])
            seg_a_dt=ang*dt_
            tt1,tt2=eta1,eta2
        lst_eta1.append(eta1)
        lst_eta2.append(eta2)
        len1+=ds1
        len2+=ds2
        lst_no.append(i)
        lst_t.append(df1.iloc[i][2])  #because that df1&2 ages are synchronized is required here, so df2.iloc[i][2] is also ok
        lst_seg_a.append(format(ang,'.7f').rstrip('0') if ang<.1 else ang)
        accum_seg_a+=ang  #angular difference
        accum_seg_a_dt+=seg_a_dt  #similar to function (9) in Qi16
        lst_seg_l.append(format(leh,'.7f').rstrip('0') if leh<.1 else leh)
        accum_seg_l+=leh  #length difference
        (mean_seg_a_dt,mean_seg_l)=(0.,0.) if i==0 else (accum_seg_a_dt/abs(df1.iloc[i][2]-df1.iloc[0][2]),accum_seg_l/abs(df1.iloc[i][2]-df1.iloc[0][2]))
        lst_accum_seg_a.append(accum_seg_a)
        lst_accum_seg_a_dt.append(accum_seg_a_dt)
        lst_mean_seg_a_dt.append(mean_seg_a_dt)
        lst_accum_seg_l.append(accum_seg_l)
        lst_mean_seg_l.append(mean_seg_l)
        divisor_l=PLATE_V_MAX_PAST/11.1195051975  #i.e. about 2.7 degree/myr, magnitude of velocity
        s_a,s_l=mean_seg_a_dt/POL_WAND_DIR_DIF_MAX,mean_seg_l/divisor_l
    if whole=='y': return s_a,s_l
    else:
        return pd.DataFrame({'11_ang_seg_dif':lst_seg_a,
                             '12_ang_seg_dif_accum':lst_accum_seg_a,
                             '13_ang_seg_dif_dt_accum':lst_accum_seg_a_dt,
                             '14_ang_seg_dif_dt_mean':lst_mean_seg_a_dt,
                             '21_len_seg_dif':lst_seg_l,
                             '22_len_seg_dif_accum':lst_accum_seg_l,
                             '23_len_seg_dif_mean':lst_mean_seg_l,
                             '33_spa_dif':0,'90_no':lst_no,
                             '91_course1':lst_eta1,'92_course2':lst_eta2,
                             '94_tstop':lst_t})

@jit(nopython=False,parallel=True)
def sig_test_shape_dif(_a_,_b_,fmt1='textfile',fmt2='textfile',fn1='file1',fn2='file2'):
    """d_a_tested,d_l_tested:
             angular and length differences with sig test on geometry metrics
       d_a,d_l:
             angular and length differences without sig test on geometry metrics
    Source: @__author__ and Chris Rowan, Nov2017"""
    _,d_a,d_l=shape_dif(_a_,_b_,fmt1,fmt2,whole='y')
    lst_d_shp_a_ras,lst_d_shp_ras_rbs=[],[]
    d_a_tested=d_a
    d_l_tested=d_l
    for _ in np.arange(1000):
        ra_=btr_(_a_,_b_,fmt1=fmt1,fmt2=fmt2,fna=fn1)
        rb_=btr_(_b_,_a_,fmt1=fmt2,fmt2=fmt1,fna=fn2)
        d_shp_a_ra,_,_=shape_dif(_a_,ra_,fmt1,fmt2='df',whole='y')
        d_shp_ra_rb,_,_=shape_dif(ra_,rb_,fmt1='df',fmt2='df',whole='y')
        lst_d_shp_a_ras.append(d_shp_a_ra)
        lst_d_shp_ras_rbs.append(d_shp_ra_rb)
    _u_=np.percentile(lst_d_shp_a_ras,97.5)
    _l_=np.percentile(lst_d_shp_ras_rbs,2.5)
    if _u_>_l_: d_a_tested,d_l_tested=0.,0.
    return d_a_tested,d_l_tested,d_a,d_l

def apwp_dif_azi_without_same_fromp(trj1,trj2,fmt1='textfile',fmt2='textfile'):
    """APWP difference based on azimuth without same from p
    ToDo: mapproject -Af depends on the "FROM" point location, so revisition
    is needed                                  Source: @__author__, 2016-2017"""
    df1=trj1 if fmt1=='df' else txt2df_awk(trj1)  #sep default as tab
    df2=trj2 if fmt2=='df' else txt2df_awk(trj2)
    df1,df2=df1[df1[2].isin(df2[2])],df2[df2[2].isin(df1[2])]  #remove age-unpaired rows in both dataframes
    df_space_score=btr(trj1,trj2,fmt1,fmt2)  #significant spacial distance/difference
    w_s,w_a,tt1,tt2,len1,len2=1/3,1/3,0,0,0,0  #tt1/2 intermedium segment azimuth for trj 1/2; len1/2 total length of trj 1/2
    accum_seg_a,accum_seg_l,accum_seg_a_dt=0,0,0
    n_row=min(len(df1.index),len(df2.index))
    lst_no,lst_t,lst_d_s,lst_d,lst_eta1,lst_eta2=[],[],[],[],[],[]
    lst_seg_a,lst_accum_seg_a,lst_accum_seg_a_dt,lst_mean_seg_a_dt=[],[],[],[]  #directional diff
    lst_seg_l,lst_accum_seg_l,lst_mean_seg_l=[],[],[]  #segment length diff
    for i in range(0,n_row):
        if i==0:
            ang,leh,seg_a_dt,ds1,ds2,eta1,eta2=0.,0.,0.,0.,0.,0.,0.  #ds1/2 intermedium segment GCD for trj 1/2
            dt_=abs(df1.iloc[i+1][2]-df1.iloc[i][2])  #what it is doesn't matter
        else:
            if df1.iloc[i-1][0]==df1.iloc[i][0] and df1.iloc[i-1][1]==df1.iloc[i][1]:
                ds1=0.
                eta1=tt1
            else:
                ds1=s2f(PMAGPY36().angle((df1.iloc[i-1][0],df1.iloc[i-1][1]),
                                         (df1.iloc[i][0],df1.iloc[i][1])))
                eta1=s2f(run_sh(AZI.format(df1.iloc[i-1][0],
                                           df1.iloc[i-1][1],df1.iloc[i][0],
                                           df1.iloc[i][1])).decode().rstrip('\n'))
            if df2.iloc[i-1][0]==df2.iloc[i][0] and df2.iloc[i-1][1]==df2.iloc[i][1]:
                ds2=0.
                eta2=tt2
            else:
                ds2=s2f(PMAGPY36().angle((df2.iloc[i-1][0],df2.iloc[i-1][1]),
                                         (df2.iloc[i][0],df2.iloc[i][1])))
                eta2=s2f(run_sh(AZI.format(df2.iloc[i-1][0],
                                           df2.iloc[i-1][1],df2.iloc[i][0],
                                           df2.iloc[i][1])).decode().rstrip('\n'))
            if tt1>eta1 and tt1-eta1>180.: eta1+=360.  #needs to brainstorm for a while, but now it should be correct
            if eta1>tt1 and eta1-tt1>180. and i>1: eta1-=360.
            if tt2>eta2 and tt2-eta2>180.: eta2+=360.  #i.e. eta1/2 (Course) can be more than 360
            if eta2>tt2 and eta2-tt2>180. and i>1: eta2-=360.
            dt_=abs(df1.iloc[i][2]-df1.iloc[i-1][2])
            #ang=abs(eta2-eta1)  #according to Course (谢等2003描述是转角的叠加,设置CW or CCW为正); Through re-thinking, this is a
            #good solution for closed polygons' comparison, not good for 2 trajectories
            ang=s2f(run_sh(SMALLER_ANGLE_REMAINDER.format(eta2,eta1)).decode().rstrip('\n'))
            if ang>180.: ang=360.-ang
            leh=abs(ds1-ds2)
            seg_a_dt=ang*dt_
            tt1,tt2=eta1,eta2
        lst_no.append(i)
        lst_t.append(df1.iloc[i][2])  #because that df1 ages and df2 ages are synchronized is required here, so df2.iloc[i][2] is also ok
        lst_seg_a.append(ang)
        accum_seg_a+=ang  #angular difference
        accum_seg_a_dt+=seg_a_dt  #similar to function (9) in Qi16
        lst_seg_l.append(leh)
        accum_seg_l+=leh  #length difference
        (mean_seg_a_dt,mean_seg_l)=(0.,0.) if i==0 else (accum_seg_a_dt/abs(df1.iloc[i][2]-df1.iloc[0][2]),accum_seg_l/abs(df1.iloc[i][2]-df1.iloc[0][2]))
        lst_accum_seg_a.append(accum_seg_a)
        lst_accum_seg_a_dt.append(accum_seg_a_dt)
        lst_mean_seg_a_dt.append(mean_seg_a_dt)
        lst_accum_seg_l.append(accum_seg_l)
        lst_mean_seg_l.append(mean_seg_l)
        lst_eta1.append(eta1)
        lst_eta2.append(eta2)
        len1+=ds1
        len2+=ds2
        divisor_l=PLATE_V_MAX_PAST/11.1195051975
        d_s=1.-df_space_score.loc[0:i,"31_spa_pol_dif"].sum()/(i+1)
        lst_d_s.append(d_s)
        dif=0. if i==0 else w_a*mean_seg_a_dt/POL_WAND_DIR_DIF_MAX+(1.-w_s-w_a)*mean_seg_l/divisor_l+w_s*d_s
        #Qi16 functions might referred to Su15
        lst_d.append(dif)
    print('Attn: For total dif, weights Ws & Wa are {} and {} repectively'.format(w_s,w_a))
    return pd.DataFrame({'00_dif':lst_d,'11_ang_seg_dif':lst_seg_a,
                         '12_ang_seg_dif_accum':lst_accum_seg_a,
                         '13_ang_seg_dif_dt_accum':lst_accum_seg_a_dt,
                         '14_ang_seg_dif_dt_mean':lst_mean_seg_a_dt,
                         '21_len_seg_dif':lst_seg_l,
                         '22_len_seg_dif_accum':lst_accum_seg_l,
                         '23_len_seg_dif_mean':lst_mean_seg_l,
                         '33_spa_dif':lst_d_s,'90_no':lst_no,
                         '91_course1':lst_eta1,'92_course2':lst_eta2,
                         '94_tstop':lst_t})

def apwp_dif(trj1,trj2,fmt1='textfile',fmt2='textfile'):
    """Initial angular difference was original azimuth difference, which here is
    modified to relative to previous segment (finally actually to the first
    srgment), which is inspired by a talk with Dr Chris Rowan
    Source: @__author__, 2016-2017"""
    df1=trj1 if fmt1=='df' else txt2df_awk(trj1)  #sep default as tab
    df2=trj2 if fmt2=='df' else txt2df_awk(trj2)
    df1,df2=df1[df1[2].isin(df2[2])],df2[df2[2].isin(df1[2])]  #remove age-unpaired rows in both dataframes
    df_space_score=btr(trj1,trj2,fmt1,fmt2)  #significant spacial distance/difference
    w_s,w_a,tt1,tt2,len1,len2=1/3,1/3,0,0,0,0  #tt1/2 intermedium segment azimuth for trj 1/2; len1/2 total length of trj 1/2
    accum_seg_a,accum_seg_l,accum_seg_a_dt=0,0,0
    n_row=min(len(df1.index),len(df2.index))
    lst_no,lst_t,lst_d_s,lst_d,lst_eta1,lst_eta2=[],[],[],[],[],[]
    lst_seg_a,lst_accum_seg_a,lst_accum_seg_a_dt,lst_mean_seg_a_dt=[],[],[],[]  #directional diff
    lst_seg_l,lst_accum_seg_l,lst_mean_seg_l=[],[],[]  #segment length diff
    for i in range(0,n_row):
        if i==0:
            ang,leh,seg_a_dt,ds1,ds2,eta1,eta2=0.,0.,0.,0.,0.,0.,0.  #ds1/2 intermedium segment GCD for trj 1/2
            dt_=abs(df1.iloc[i+1][2]-df1.iloc[i][2])  #what it is doesn't matter
        elif i==n_row-1:
            if df1.iloc[i-1][0]==df1.iloc[i][0] and df1.iloc[i-1][1]==df1.iloc[i][1]:
                ds1=0.
            else:
                ds1=s2f(PMAGPY36().angle((df1.iloc[i-1][0],df1.iloc[i-1][1]),
                                         (df1.iloc[i][0],df1.iloc[i][1])))
            if df2.iloc[i-1][0]==df2.iloc[i][0] and df2.iloc[i-1][1]==df2.iloc[i][1]:
                ds2=0.
            else:
                ds2=s2f(PMAGPY36().angle((df2.iloc[i-1][0],df2.iloc[i-1][1]),
                                         (df2.iloc[i][0],df2.iloc[i][1])))
            dt_=abs(df1.iloc[i][2]-df1.iloc[i-1][2])
            eta1,eta2,ang,seg_a_dt=0.,0.,0.,0.
            leh=abs(ds1-ds2)
            tt1,tt2=eta1,eta2
        else:
            if df1.iloc[i-1][0]==df1.iloc[i][0] and df1.iloc[i-1][1]==df1.iloc[i][1]:
                ds1=0.
                eta1=tt1
            else:
                ds1=s2f(PMAGPY36().angle((df1.iloc[i-1][0],df1.iloc[i-1][1]),
                                         (df1.iloc[i][0],df1.iloc[i][1])))
                eta1=ang4_2suc_disp_direc_gdesics(df1.iloc[i-1][0],
                                                  df1.iloc[i-1][1],
                                                  df1.iloc[i][0],
                                                  df1.iloc[i][1],
                                                  df1.iloc[i+1][0],
                                                  df1.iloc[i+1][1])
            if df2.iloc[i-1][0]==df2.iloc[i][0] and df2.iloc[i-1][1]==df2.iloc[i][1]:
                ds2=0.
                eta2=tt2
            else:
                ds2=s2f(PMAGPY36().angle((df2.iloc[i-1][0],df2.iloc[i-1][1]),
                                         (df2.iloc[i][0],df2.iloc[i][1])))
                eta2=ang4_2suc_disp_direc_gdesics(df2.iloc[i-1][0],
                                                  df2.iloc[i-1][1],
                                                  df2.iloc[i][0],
                                                  df2.iloc[i][1],
                                                  df2.iloc[i+1][0],
                                                  df2.iloc[i+1][1])
            ang=360-abs(eta2-eta1) if abs(eta2-eta1)>180 else abs(eta2-eta1)
            leh=abs(ds1-ds2)
            dt_=abs(df1.iloc[i][2]-df1.iloc[i-1][2])
            seg_a_dt=ang*dt_
            tt1,tt2=eta1,eta2
        lst_no.append(i)
        lst_t.append(df1.iloc[i][2])  #because that df1 ages and df2 ages are synchronized is required here, so df2.iloc[i][2] is also ok
        lst_seg_a.append(format(ang,'.7f').rstrip('0') if ang<.1 else ang)
        accum_seg_a+=ang  #angular difference
        accum_seg_a_dt+=seg_a_dt  #similar to function (9) in Qi16 *ds
        lst_seg_l.append(format(leh,'.7f').rstrip('0') if leh<.1 else leh)
        accum_seg_l+=leh  #length difference
        (mean_seg_a_dt,mean_seg_l)=(0.,0.) if i==0 else (accum_seg_a_dt/abs(df1.iloc[i][2]-df1.iloc[0][2]),accum_seg_l/abs(df1.iloc[i][2]-df1.iloc[0][2]))
        lst_accum_seg_a.append(accum_seg_a)
        lst_accum_seg_a_dt.append(accum_seg_a_dt)
        lst_mean_seg_a_dt.append(mean_seg_a_dt)
        lst_accum_seg_l.append(accum_seg_l)
        lst_mean_seg_l.append(mean_seg_l)
        lst_eta1.append(eta1)
        lst_eta2.append(eta2)
        len1+=ds1
        len2+=ds2
        divisor_l=PLATE_V_MAX_PAST/11.1195051975
        d_s=1.-df_space_score.loc[0:i,"31_spa_pol_dif"].sum()/(i+1)
        lst_d_s.append(d_s)
        dif=0. if i==0 else w_a*mean_seg_a_dt/POL_WAND_DIR_DIF_MAX+(1.-w_s-w_a)*mean_seg_l/divisor_l+w_s*d_s
        lst_d.append(format(dif,'.7f').rstrip('0') if dif<.1 else dif)
    print('Attn: For total dif, weights Ws & Wa are {} and {} repectively'.format(w_s,w_a))
    return pd.DataFrame({'00_dif':lst_d,'11_ang_seg_dif':lst_seg_a,
                         '12_ang_seg_dif_accum':lst_accum_seg_a,
                         '13_ang_seg_dif_dt_accum':lst_accum_seg_a_dt,
                         '14_ang_seg_dif_dt_mean':lst_mean_seg_a_dt,
                         '21_len_seg_dif':lst_seg_l,
                         '22_len_seg_dif_accum':lst_accum_seg_l,
                         '23_len_seg_dif_mean':lst_mean_seg_l,
                         '33_spa_dif':lst_d_s,'90_no':lst_no,
                         '91_course1':lst_eta1,'92_course2':lst_eta2,
                         '94_tstop':lst_t})

def apwp_dif_with_shape_sig_test(trj1,trj2,fmt1='textfile',fmt2='textfile',lag=0,hag=530):
    """For each path, directional change of each segment is from this segment's
    previous segment; shape_sig_test is based on whole path comparison
    Source: @__author__ and Chris Rowan, 2016-2017"""
    filname1=re.split('/|\.',trj1)[-2] if fmt1=='textfile' else str(uuid.uuid4())
    filname2=re.split('/|\.',trj2)[-2] if fmt2=='textfile' else str(uuid.uuid4())
    df1=trj1 if fmt1=='df' else txt2df_awk(trj1)  #sep default as tab
    df2=trj2 if fmt2=='df' else txt2df_awk(trj2)
    df1,df2=df1[(df1[2].isin(df2[2])) & (df1[2]>=lag) & (df1[2]<=hag)],df2[(df2[2].isin(df1[2])) & (df2[2]>=lag) & (df2[2]<=hag)]  #remove age-unpaired rows in both dataframes
    df_space_score=btr(df1,df2,fmt1='df',fmt2='df',fn1=filname1,fn2=filname2)  #significant spacial distance/difference
    w_s,w_a=1/3,1/3
    n_row=min(len(df1.index),len(df2.index))
    lst_no,lst_t0,lst_t,lst_d_s,lst_d_a,lst_d_l,lst_d=[],[],[],[],[],[],[]
    lst_d_a_tested,lst_d_l_tested=[],[]
    print("Dspa\tDang\tDangT\tDlen\tDlenT\tTstart(Ma)\tTstop(Ma)")
    inc=n_row-1  #1(can be invisible) 2 ...; n_row-1 makes the measurement once on whole path
    for i in range(0,n_row,inc):
        lst_no.append(i)
        lst_t0.append(df1.iloc[0][2])
        lst_t.append(df1.iloc[i][2])  #because that df1 ages and df2 ages are synchronized is required here, so df2.iloc[i][2] is also ok
        (d_a_tested,d_l_tested,d_a,d_l)=(0.,0.,0.,0.) if i==0 else sig_test_shape_dif(df1.iloc[0:i+1],df2.iloc[0:i+1],fmt1='df',fmt2='df',fn1=filname1,fn2=filname2)
        d_s=1.-df_space_score.loc[0:i,"31_spa_pol_dif"].sum()/(i+1)  #path space difference
        lst_d_a.append(d_a)
        lst_d_a_tested.append(d_a_tested)
        lst_d_l.append(d_l)
        lst_d_l_tested.append(d_l_tested)
        lst_d_s.append(d_s)
        dif=0. if i==0 else w_s*d_s+w_a*d_a+(1-w_s-w_a)*d_l
        lst_d.append(format(dif,'.7f').rstrip('0') if dif<.1 else dif)
        if i>0: print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(d_s,d_a,d_a_tested,d_l,d_l_tested,df1.iloc[0][2],df1.iloc[i][2]))
    print('Attn: For total dif, weights Ws & Wa are {} and {} repectively'.format(w_s,w_a))
    return pd.DataFrame({'00_dif':lst_d,'33_spa_dif':lst_d_s,
                         '34_ang_dif':lst_d_a,'35_len_dif':lst_d_l,
                         '36_ang_dif_tested':lst_d_a_tested,
                         '37_len_dif_tested':lst_d_l_tested,
                         '90_no':lst_no,'93_tstart':lst_t0,'94_tstop':lst_t})

def apwp_dif_with_ang_len_sig_tests(trj1,trj2,fmt1='textfile',fmt2='textfile',lag=0,hag=530):
    """Derived from function "apwp_dif_with_shape_sig_test"
    Source: @__author__ and Chris Rowan, Jan2018"""
    filname1=re.split('/|\.',trj1)[-2] if fmt1=='textfile' else str(uuid.uuid4())
    filname2=re.split('/|\.',trj2)[-2] if fmt2=='textfile' else str(uuid.uuid4())
    df1=trj1 if fmt1=='df' else txt2df_awk(trj1)  #sep default as tab
    df2=trj2 if fmt2=='df' else txt2df_awk(trj2)
    df1,df2=df1[(df1[2].isin(df2[2])) & (df1[2]>=lag) & (df1[2]<=hag)],df2[(df2[2].isin(df1[2])) & (df2[2]>=lag) & (df2[2]<=hag)]  #remove age-unpaired rows in both dataframes
    df_space_score=btr(df1,df2,fmt1='df',fmt2='df',fn1=filname1,fn2=filname2)  #significant spacial distance/difference
    w_s,w_a=1/3,1/3
    n_row=min(len(df1.index),len(df2.index))
    lst_no,lst_t0,lst_t,lst_d_s,lst_d_a,lst_d_l,lst_d=[],[],[],[],[],[],[]
    lst_d_a_tested,lst_d_l_tested=[],[]
    print("Dspa\tDang\tDangT\tDlen\tDlenT\tTstart(Ma)\tTstop(Ma)")
    inc=n_row-1  #1(can be invisible) 2
    for i in range(0,n_row,inc):
        lst_no.append(i)
        lst_t0.append(df1.iloc[0][2])
        lst_t.append(df1.iloc[i][2])  #because that df1 ages and df2 ages are synchronized is required here, so df2.iloc[i][2] is also ok
        (d_a_tested,d_l_tested,d_a,d_l)=(0.,0.,0.,0.) if i==0 else sig_test_shape_dif(df1.iloc[0:i+1],df2.iloc[0:i+1],fmt1='df',fmt2='df',fn1=filname1,fn2=filname2)
        d_s=1.-df_space_score.loc[0:i,"31_spa_pol_dif"].sum()/(i+1)  #path space difference
        lst_d_a.append(d_a)
        lst_d_a_tested.append(d_a_tested)
        lst_d_l.append(d_l)
        lst_d_l_tested.append(d_l_tested)
        lst_d_s.append(d_s)
        dif=0. if i==0 else w_s*d_s+w_a*d_a+(1-w_s-w_a)*d_l
        lst_d.append(format(dif,'.7f').rstrip('0') if dif<.1 else dif)
        if i>0: print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(d_s,d_a,d_a_tested,d_l,d_l_tested,df1.iloc[0][2],df1.iloc[i][2]))
    print('Attn: For total dif, weights Ws & Wa are {} and {} repectively'.format(w_s,w_a))
    return pd.DataFrame({'00_dif':lst_d,'33_spa_dif':lst_d_s,
                         '34_ang_dif':lst_d_a,'35_len_dif':lst_d_l,
                         '36_ang_dif_tested':lst_d_a_tested,
                         '37_len_dif_tested':lst_d_l_tested,
                         '90_no':lst_no,'93_tstart':lst_t0,'94_tstop':lst_t})

def run_sh(script,stdin=None):
    """Raises error on non-zero return code
    Source: http://stackoverflow.com/questions/2651874/embed-bash-in-python"""
    #Note: by using a list here (['bash', ...]) we avoid quoting issues, as the
    #arguments are passed in exactly this order (spaces, quotes, and newlines
    #won't cause problems):
    proc=subprocess.Popen(['bash','-c',script],stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,stdin=subprocess.PIPE)
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



class PMAGPY36(object):
    """PmagPy (https://pmagpy.github.io/) functions below are converted to run
    correctly with Python3.6.x"""

    def __init__(self):
        pass

    @staticmethod
    @jit(nopython=False,parallel=True)
    def dir2cart(dis):
        '''converts list/array of vector directions, in degrees, to array of
        Cartesian coordinates, in x,y,z             Source: pmag.py of PmagPy'''
        ints=np.ones(len(dis)).transpose() #get an array of ones to plug into dec,inc pairs
        _d_=np.array(dis)
        rad=np.pi/180.
        if len(_d_.shape)>1: # array of vectors
            decs,incs=_d_[:,0]*rad,_d_[:,1]*rad
            if _d_.shape[1]==3: ints=_d_[:,2] #take the given lengths assigned in 3rd column in _d_
        else: # single vector
            decs,incs=np.array(_d_[0])*rad,np.array(_d_[1])*rad
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
        rad=np.pi/180. # constant to convert degrees to radians
        if len(cart.shape)>1: x_s,y_s,z_s=cart[:,0],cart[:,1],cart[:,2]
        else: x_s,y_s,z_s=cart[0],cart[1],cart[2] #single vector
        r_s=np.sqrt(x_s**2+y_s**2+z_s**2) # calculate resultant vector length
        decs=(np.arctan2(y_s,x_s)/rad)%360. # calculate declination taking care of correct quadrants (arctan2) and making modulo 360.
        try: incs=np.arcsin(z_s/r_s)/rad # calculate inclination (converting to degrees)
        except ZeroDivisionError:
            print("Error: division by 0. z_s is {0}, r_s is {1}".format(z_s,r_s))  # r_s=0
            return np.zeros(3)
        if incs==np.nan:
            print("z_s {0} > r_s {1}, arcsin value does not exist!".format(z_s,r_s))
            return np.zeros(3)
        return np.array([decs,incs,r_s]).transpose() # return the directions list

    def fisher_mean(self,data):
        """Source: pmag.py of PmagPy
        Used enumerate instead of range+len by @__author__, 2017"""
        #calculates fisher parameters for data
        xbar,fpars=[0,0,0],{}
        _n_=len(data) #origianl code line, mainly for many directions, i.e. dec,inc
        if _n_<2: return fpars
        _x_=self.dir2cart(data)
        #for rec in data: _x_.append(dir2cart([rec[0],rec[1],1.]))
        for i,_ in enumerate(_x_):
            for j in range(3): xbar[j]+=_x_[i][j]
        _r_=np.sqrt(sum(xbar[c]**2 for c in range(3)))
        for k in range(3): xbar[k]/=_r_
        dire=self.cart2dir(xbar)
        fpars["dec"],fpars["inc"],fpars["n"],fpars["r"]=dire[0],dire[1],_n_,_r_
        if abs(_n_-_r_)>1E-08:
            k=(_n_-1.)/(_n_-_r_)
            fpars["k"]=k
            csd=81./np.sqrt(k)
        else: fpars['k'],csd=100000,0.  #fpars['k']='inf'
        _b_=20.**(1./(_n_-1.))-1
        _a_=1-_b_*(_n_-_r_)/_r_
        if _a_<-1: _a_=-1
        a95=np.arccos(_a_)*180./np.pi
        fpars["alpha95"],fpars["csd"]=a95,csd #estimated angular standard deviation, CSD
        if _a_<0: fpars["alpha95"]=180.
        return fpars

    @jit(nopython=False,parallel=True)
    def angle(self,d_1,d_2):
        """call to angle(d_1,d_2) returns array of angles between lists of 2
        directions d_1,d_2 where d_1 is, for example,
        [[Dec1,Inc1],[Dec2,Inc2],etc.]                 Source: pmag.py of PmagPy
        Modified a bit in style by @__author__, 2017"""
        d_1,d_2=np.array(d_1),np.array(d_2)
        d_1=d_1[:,0:2] if len(d_1.shape)>1 else d_1[:2]	# strip off intensity
        d_2=d_2[:,0:2] if len(d_2.shape)>1 else d_2[:2]	# strip off intensity
        x_1,x_2=self.dir2cart(d_1),self.dir2cart(d_2) # convert to cartesian from polar
        angles=[] # set up a list for angles
        for k in range(x_1.shape[0]): # single vector
            angles.append((np.arccos(np.dot(x_1[k],x_2[k]))*180./np.pi)%360.) # take the dot product
        return np.array(angles)

    @staticmethod
    @jit(nopython=False,parallel=True)
    def pseudo(dis):
        """draw a bootstrap sample of DIrectionS    Source: pmag.py of PmagPy"""
        siz=len(dis)
        return np.array(dis)[np.random.randint(siz,size=siz)]

    @jit(nopython=False,parallel=True)
    def di_boot(self,dis,nob=5000):
        """returns bootstrap parameters for Directional data     Source: pmag.py
        of PmagPy                Modified a bit in style by @__author__, 2017"""
        #fpars=self.fisher_mean(dis)  #get average DI for whole dataset
        #now do bootstrap to collect d_i bootstrap means
        d_i=[]  # list of bootstrap directions
        for _ in range(nob): # repeat nob (number of bootstraps) times
            #if k%50==0: print(k,' out of ',nob)
            pdis=self.pseudo(dis) # get a pseudosample
            bfpars=self.fisher_mean(pdis) # get bootstrap mean bootstrap sample
            d_i.append([bfpars['dec'],bfpars['inc']])
        return d_i

    @jit(nopython=False,parallel=True)
    def dogeo(self,dec,inc,azi,plg):
        """rotates dec,inc into geographic coordinates using azi,plg as azimuth
        and plunge of _x_ direction                    Source: pmag.py of PmagPy
        Modified a bit in style by @__author__, 2017"""
        dir_=[dec,inc,1.] # put dec inc in direction list and set  length to unity
        _x_=self.dir2cart(dir_) # get cartesian coordinates
        # set up rotation matrix
        a_1,a_2,a_3=self.dir2cart([azi,plg,1.]),self.dir2cart([azi+90.,0,1.]),self.dir2cart([azi-180.,90.-plg,1.])
        # do rotation
        xp_=a_1[0]*_x_[0]+a_2[0]*_x_[1]+a_3[0]*_x_[2]
        yp_=a_1[1]*_x_[0]+a_2[1]*_x_[1]+a_3[1]*_x_[2]
        zp_=a_1[2]*_x_[0]+a_2[2]*_x_[1]+a_3[2]*_x_[2]
        # transform back to dec,inc
        dir_geo=self.cart2dir([xp_,yp_,zp_])
        return dir_geo[0],dir_geo[1]    # send back declination and inclination

    @jit(nopython=False,parallel=True)
    def dodirot(self,dec,inc,dbar,ibar):
        """dec=declination,inc=inclination, dbar/ibar are the desired mean direction.
        Returns the rotated Dec/Inc pair            Source: pmag.py of PmagPy"""
        _d_,irot=self.dogeo(dec,inc,dbar,90.-ibar)
        drot=_d_-180.
        #drot,irot=dogeo(dec,inc,Dbar,Ibar)
        if drot<360.: drot=drot+360.
        if drot>360.: drot=drot-360.
        return drot,irot

    @staticmethod
    @jit(nopython=True,parallel=True)
    def fshdev(kap):
        """kap is kappa, returns a direction from distribution with mean
        declination of 0, inclination of 90 and kappa of kap
        Source: pmag.py of PmagPy"""
        r_1,r_2=random.random(),random.random()
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

def main():
    """Run this module"""
    #-------------Prepare Model Predicted APWP----------------------------------
    #-------Following 6 lines can be replaced by just one line of awk code:-----
    #--awk '{print $1,$2,$3,$5/111.195051975,$6/111.195051975,$4,0,0,0,0}' \
    #OFS='\t' 101FHS140predictPWP105.d >/tmp/101FHS140predictPWP105.txt---------
    model_apwp=txt2df_awk('/home/i/Desktop/git/digivisual/tmp/701FHS120predictPWP105.d')
    model_apwp[4]/=111.195051975
    model_apwp[5]/=111.195051975
    model_apwp[6],model_apwp[7],model_apwp[8],model_apwp[9]=0,0,0,0
    model_apwp=model_apwp[[0,1,2,4,5,3,6,7,8,9]]
    model_apwp.rename(columns={4:3,5:4,3:5},inplace=True)  #----------END-------

    #-NA APWPs from Different Algorithms, versus FHS Model Predicted APWP-------
    tbin=10		#18,10,2
    tstep=5	#9,5,1
    platemodel='dm16'
    plateid='701comb'
    root_o_dir='/run/media/i/s'
    for mav in [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]:
        for wgt in [0,1,2,3,4,5]:
            pmag_apwp=txt2df_awk(os.path.join(root_o_dir,platemodel+"_"+plateid,
                                              platemodel+"_"+plateid+"_"+str(tbin)+"_"+str(tstep)+"_"+str(mav)+"_"+str(wgt)+".txt"))
            os.makedirs('/tmp/traj1',exist_ok=True)
            for row in pmag_apwp.itertuples():
                print(row[3])
                if row[8]>25:
                    raw_dir1=os.path.join(root_o_dir,platemodel,plateid,
                                          platemodel+"_"+plateid+"_"+str(tbin)+"_"+str(tstep)+"_"+str(mav)+"_"+str(wgt),
                                          str(round(row[3]+tstep,1))+"_"+str(round(row[3]-tstep,1))+".d")
                    raw_dir2=os.path.join(root_o_dir,platemodel,plateid,
                                          platemodel+"_"+plateid+"_"+str(tbin)+"_"+str(tstep)+"_"+str(mav)+"_"+str(wgt),
                                          str(round(row[3]+tstep,1))+"_"+str(round(row[3]-tstep,1))+"_WT.d")
                    raw_poles=pd.read_table(raw_dir1,
                                            header=0,
                                            comment='>',
                                            dtype=str) if os.path.isfile(raw_dir1) else pd.read_table(raw_dir2,
                                                                                                      header=0,
                                                                                                      comment='>',
                                                                                                      dtype=str)
                    raw_poles[12]=1.0
                    raw_poles[13]=(pd.to_numeric(raw_poles['LOMAGAGE'])+
                                   pd.to_numeric(raw_poles['HIMAGAGE']))/2
                    raw_poles=raw_poles[['PLONG','PLAT',12,13]]
                    print(raw_poles)
                    raw_poles.to_csv('/tmp/traj1/{:s}.txt'.format(str(int(row[3]))),
                                     sep='\t',encoding='utf-8',header=False,
                                     index=False)
            print('-----------{}----DOUBLE-CHECK----{}----------'.format(mav,wgt))
            os.makedirs(os.path.join(root_o_dir,platemodel+"_"+plateid,
                                     str(tbin)+"_"+str(tstep)+"_simil"),exist_ok=True)
            print(pmag_apwp)
            print(model_apwp)
            simil=apwp_dif(pmag_apwp,model_apwp,'df','df')
            simil.to_csv(os.path.join(root_o_dir,platemodel+"_"+plateid,
                                      str(tbin)+"_"+str(tstep)+"_simil",
                                      platemodel+"_"+plateid+"_"+str(tbin)+"_"+str(tstep)+"_"+str(mav)+"_"+str(wgt)+".d"),
                         index=False,header=True,sep='	')
            print('-----------{}----DOUBLE-CHECK----{}---END----'.format(mav,wgt))

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

    apwp1=open('/tmp/1.d','w')
    apwp1.write(pwp1)
    apwp1.close()
    apwp2=open('/tmp/2.d','w')
    apwp2.write(pwp2)
    apwp2.close()
    os.makedirs('/tmp/traj1',exist_ok=True)
    apwp1_0=open('/tmp/traj1/0.txt','w')
    apwp1_0.write(pwp1_0ma)
    apwp1_0.close()
    apwp1_5=open('/tmp/traj1/5.txt','w')
    apwp1_5.write(pwp1_5ma)
    apwp1_5.close()

    diff=apwp_dif('/tmp/1.d','/tmp/2.d')
    print(diff)

if __name__=="__main__": main()
