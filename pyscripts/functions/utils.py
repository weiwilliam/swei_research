__all__ = ['ndate','setup_cmap','cnbestF','latlon_news','lat_ns','lon_we','gen_eqs_by_stats',
           'find_cnlvs']
import numpy as np

def ndate(hinc,cdate):
    from datetime import datetime
    from datetime import timedelta
    yy=int(str(cdate)[:4])
    mm=int(str(cdate)[4:6])
    dd=int(str(cdate)[6:8])
    hh=int(str(cdate)[8:10])
    dstart=datetime(yy,mm,dd,hh)
    dnew=dstart+timedelta(hours=hinc)
    dnewint=int(str('%4.4d' % dnew.year)+str('%2.2d' % dnew.month)+
                str('%2.2d' %dnew.day)+str('%2.2d' % dnew.hour))
    return dnewint
#
# Set colormap through NCL colormap and index
#
def setup_cmap(name,idxlst):
    import os, platform
    os_name=platform.system()
    if (os_name=='Darwin'):
        rootpath='/Users/weiwilliam'
    elif (os_name=='Windows'):
        rootpath='F:\GoogleDrive_NCU\Albany'
    elif (os_name=='Linux'):
        if (os.path.exists('/glade')):
            rootpath='/glade/u/home/swei/research/pyscripts'
        if (os.path.exists('/scratch')):
            rootpath='/home/swei/research/pyscripts'
    import matplotlib.colors as mpcrs
    import numpy as np
    if (os_name!='Linux'):
        nclcmap=rootpath+'/AlbanyWork/Utility/colormaps'
    else:
        nclcmap=rootpath+'/colormaps'
    
    cmapname=name
    f=open(nclcmap+'/'+cmapname+'.rgb','r')
    a=[]
    for line in f.readlines():
        if ('ncolors' in line):
            clnum=int(line.split('=')[1])
        a.append(line)
    f.close()
    b=a[-clnum:]
    c=[]
    selidx=np.array(idxlst,dtype='int')
    if ('MPL' in name):
       for i in selidx[:]:
          if (i==0):
             c.append(tuple(float(y) for y in [1,1,1]))
          elif (i==1):
             c.append(tuple(float(y) for y in [0,0,0]))
          elif (i==-1):
             c.append(tuple(float(y) for y in [0.5,0.5,0.5]))
          else:
             c.append(tuple(float(y) for y in b[i-2].split('#',1)[0].split()))
    else:
       for i in selidx[:]:
          if (i==0):
             c.append(tuple(float(y)/255. for y in [255,255,255]))
          elif (i==1):
             c.append(tuple(float(y)/255. for y in [0,0,0]))
          elif (i==-1):
             c.append(tuple(float(y)/255. for y in [128,128,128]))
          else:
             c.append(tuple(float(y)/255. for y in b[i-2].split('#',1)[0].split()))

    d=mpcrs.LinearSegmentedColormap.from_list(name,c,selidx.size)
    return d

def cnbestF(data):
    import numpy as np
    std=np.nanstd(data)
    mean=np.nanmean(data)
    vmax=np.nanmax(abs(data))
    if (vmax>5*(mean+std*3)):
        cnvmax=mean+std*4
    else:
        cnvmax=vmax
    ccnvmax='%e'%(cnvmax)
    tmp1=ccnvmax.find('-')
    tmp2=ccnvmax.find('+')
    if (tmp1<0):
        tmp=tmp2
    if (tmp2<0):
        tmp=tmp1
    d=int(ccnvmax[tmp:])
    cnmaxF=np.ceil(float(ccnvmax[:tmp-1]))*10**d
    return cnmaxF

def latlon_news(plat,plon):
    deg_sym=u'\u00B0'
    if (plat > 0.):
        ns='N'
    elif (plat < 0.):
        ns='S'
    else:
        ns=''
    if (plon > 0.):
        we='E'
    elif (plon < 0.):
        we='W'
    else:
        we=''
    txlat='%.2f%s %s'%(abs(plat),deg_sym,ns)
    txlon='%.2f%s %s'%(abs(plon),deg_sym,we)
    return txlat,txlon

def lat_ns(plat):
    deg_sym=u'\u00B0'
    if (plat > 0.):
        ns='N'
    elif (plat < 0.):
        ns='S'
    else:
        ns=''
    txlat='%.f%s %s'%(abs(plat),deg_sym,ns)
    return txlat

def lon_we(plon):
    deg_sym=u'\u00B0'
    if (plon > 0.):
       we='E'
    elif (plon < 0.):
       we='W'
    else:
       we=''
    txlon='%.f%s %s'%(abs(plon),deg_sym,we)
    return txlon

def gen_eqs_by_stats(stats_in):
    if (stats_in.intercept<0):
       fiteqs='$y=%.2fx–%.2f$' %(stats_in.slope,abs(stats_in.intercept))
    elif (stats_in.intercept>0):
       fiteqs='$y=%.2fx+%.2f$' %(stats_in.slope,abs(stats_in.intercept))
    else:
       fiteqs='y=%.2f*x' %(stats_in.slope)
    return fiteqs

def find_cnlvs(indata,ntcks=None,eqside=None):
    if not ntcks: ntcks=21
    if not eqside: eqside=0
    tmpmax=np.nanquantile(indata,0.997)
    tmpmin=np.nanquantile(indata,0.003)
    print(tmpmin,tmpmax)
    if ( abs(tmpmax)<1. and tmpmax!=0. ):
       ndecimals=int(abs(np.floor(np.log10(abs(tmpmax)))))
       cnlvmax=round(tmpmax,ndecimals)
    else:
       cnlvmax=np.sign(tmpmax)*(np.ceil(abs(tmpmax)))
    if ( abs(tmpmin)<1. and tmpmin!=0. ):
       ndecimals=int(abs(np.floor(np.log10(abs(tmpmin)))))
       cnlvmin=round(tmpmin,ndecimals)
    else:
       cnlvmin=np.sign(tmpmin)*(np.ceil(abs(tmpmin)))
    print(cnlvmin,cnlvmax)
    if (cnlvmax*cnlvmin<0):
       if (eqside):
          cnlvmax=np.max((abs(cnlvmin),abs(cnlvmax)))
          cnlvs=np.linspace(-cnlvmax,cnlvmax,ntcks)
       else:
          h_ntcks=int(ntcks*0.5)
          if ( np.mod(ntcks,2)==0 ):
             neg_lvs=np.linspace(cnlvmin,0,h_ntcks,endpoint=False)
             pos_int=(abs(cnlvmax)/int(ntcks*0.5))
             pos_lvs=np.arange(0+pos_int,cnlvmax+pos_int,pos_int)
          else:
             neg_lvs=np.linspace(cnlvmin,0,h_ntcks,endpoint=False)
             pos_lvs=np.linspace(0,cnlvmax,h_ntcks+1)
          cnlvs=np.append(neg_lvs,pos_lvs)
    else:
       if (eqside): print('Warning equal side is not applicable because max=%.f, min=%.f' %(cnlvmax,cnlvmin))
       cnlvs=np.linspace(cnlvmin,cnlvmax,ntcks)
    print(cnlvs)
    return cnlvs

