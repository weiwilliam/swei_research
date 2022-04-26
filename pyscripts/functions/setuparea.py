#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 15:10:08 2018

@author: weiwilliam

Note:
    longitude must setup between 0~360
"""

def setarea(areaname):
    if (areaname=='Glb'):
        minlat=-90.;maxlat=90.;minlon=0.;maxlon=360.
    elif (areaname=='NH'):
        minlat=20.;maxlat=80.;minlon=0.;maxlon=360.
    elif (areaname=='SH'):
        minlat=-80.;maxlat=-20.;minlon=0.;maxlon=360.
    elif (areaname == 'NAfr'):
        minlat=5.;maxlat=35.;minlon=342.;maxlon=45.
    elif (areaname == 'NAfrC'):
        minlat=5.;maxlat=35.;minlon=330.;maxlon=350.
    elif (areaname == 'NAtl'):
        minlat=5. ; maxlat=65. ; minlon=265.  ; maxlon=10.
    elif (areaname == 'SAtl'):
        minlat=-20. ; maxlat=5. ; minlon=315.  ; maxlon=10.
    elif (areaname == 'TAfr'):
        minlat=-30. ; maxlat=20. ; minlon=350.  ; maxlon=50.
    elif (areaname == 'Asia'):
        minlat=-20. ; maxlat=70. ; minlon=50. ; maxlon=160.
    elif (areaname == 'EAsia'):
        minlat=10. ; maxlat=45. ; minlon=100. ; maxlon=150.
    elif (areaname == 'SEAsia'):
        minlat=-10. ; maxlat=45. ; minlon=95. ; maxlon=160.
    elif (areaname == 'Indo'):
        minlat=-20. ; maxlat=30. ; minlon=70. ; maxlon=140.
    elif (areaname == 'ArbS'):
        minlat=0.   ; maxlat=25. ; minlon=50  ; maxlon=75.
    elif (areaname == 'SPacAtl'):
        minlat=-90. ; maxlat=10. ; minlon=180. ; maxlon=300.
    elif (areaname == 'SAmer'):
        minlat=-60. ; maxlat=10. ; minlon=275. ; maxlon=335.
    elif (areaname == 'NAmer'):
        minlat=10. ; maxlat=70. ; minlon=220. ; maxlon=310.
    elif (areaname == 'NPO'):
        minlat=50.;maxlat=90.;minlon=0.;maxlon=360.
    elif (areaname == 'NML'):
        minlat=20.;maxlat=50.;minlon=0.;maxlon=360.
    elif (areaname == 'TRO'):
        minlat=-20.;maxlat=20.;minlon=0.;maxlon=360.
    elif (areaname == 'SML'):
        minlat=-50.;maxlat=-20.;minlon=0.;maxlon=360.
    elif (areaname == 'SPO'):
        minlat=-90.;maxlat=-50.;minlon=0.;maxlon=360.
    elif (areaname == 'SO'):
        minlat=-70.;maxlat=-30.;minlon=0.;maxlon=360.
    elif (areaname == 'GLK'):
        minlat=40.;maxlat=50.;minlon=268.;maxlon=285.
    elif (areaname == 'EUO'):
        minlat=30.;maxlat=50.;minlon=0.;maxlon=60.
    elif (areaname == 'NEU'):
        minlat=65.;maxlat=75.;minlon=50.;maxlon=90.
    elif (areaname == 'Tri'):
        minlat=-10.;maxlat=10.;minlon=130.;maxlon=160.
    elif (areaname == 'Pir'):
        minlat=-20.;maxlat=20.;minlon=300.;maxlon=360.
    elif (areaname == 'Indi'):
        minlat=-20.;maxlat=20.;minlon=60.;maxlon=100.
    elif (areaname == 'NMLAtl'):
        minlat=25.;maxlat=35.;minlon=290.;maxlon=310.
    elif (areaname == 'NMLPac'):
        minlat=25.;maxlat=35.;minlon=170.;maxlon=190.
    elif (areaname == 'SMLInd'):
        minlat=-35.;maxlat=-25.;minlon=50.;maxlon=70.
    elif (areaname == 'SMLPac'):
        minlat=-35.;maxlat=-25.;minlon=190.;maxlon=210.
    elif (areaname == 'TIWA'):
        minlat=-8.;maxlat=10.;minlon=190.;maxlon=280.
    elif (areaname == 'NWPac'):
        minlat=48.;maxlat=60.;minlon=140.;maxlon=180.
    elif (areaname == 'Phil'):
        minlat=-15.;maxlat=15.;minlon=90.;maxlon=150.
    elif (areaname == 'HudB'):
        minlat=55.;maxlat=65.;minlon=268.;maxlon=285.
    elif (areaname == 'SF1'):
        minlat=-0.5;maxlat=0.5;minlon=189.5;maxlon=190.5
    elif (areaname == 'SF2'):
        minlat=-2.5;maxlat=-1.5;minlon=219.5;maxlon=220.5
    elif (areaname == 'SF3'):
        minlat=1.5;maxlat=2.5;minlon=204.5;maxlon=205.5
    elif (areaname == 'SF4'):
        minlat=1.5;maxlat=2.5;minlon=234.5;maxlon=235.5
    elif (areaname == 'SF5'):
        minlat=-0.5;maxlat=0.5;minlon=264.5;maxlon=265.5
    elif (areaname == 'ST1'):
        minlat=-2.5;maxlat=-1.5;minlon=155.5;maxlon=156.5
    elif (areaname == 'ST2'):
        minlat=-0.5;maxlat=0.5;minlon=155.5;maxlon=156.5
    elif (areaname == 'ST3'):
        minlat=1.5;maxlat=2.5;minlon=155.5;maxlon=156.5
    elif (areaname == 'ST4'):
        minlat=1.5;maxlat=2.5;minlon=146.5;maxlon=147.5
    elif (areaname == 'r2o1'):
        minlat=-10.;maxlat=40.;minlon=280.;maxlon=10.
    elif (areaname == 'r2o2'):
        minlat=-10.;maxlat=30.;minlon=320.;maxlon=10.
    elif (areaname == 'r2o3'):
        minlat=0.;maxlat=20.;minlon=115.;maxlon=145.
    elif (areaname == 'r2o4'):
        minlat=0.;maxlat=30.;minlon=290.;maxlon=10.
    elif (areaname == 'r2o5'):
        minlat=0.;maxlat=40.;minlon=340.;maxlon=20.
    elif (areaname == 'r2o6'):
        minlat=0.;maxlat=40.;minlon=320.;maxlon=20.
    elif (areaname == 'r2o7'):
        minlat=0.;maxlat=30.;minlon=270.;maxlon=320.
    elif (areaname == 'r2o8'):
        minlat=30.;maxlat=45.;minlon=350.;maxlon=40.
    elif (areaname == 'r2o9'):
        minlat=0.;maxlat=30.;minlon=330.;maxlon=10.
    elif (areaname == 'r2o10'):
        minlat=0.;maxlat=40.;minlon=330.;maxlon=20.
    elif (areaname == 'Y20Smk1'): # 2020091112
        minlat=37.;maxlat=50.;minlon=227.;maxlon=240.
    elif (areaname == 'Y20Smk2'): # 2020091500
        minlat=42.;maxlat=55.;minlon=254.;maxlon=267.
    elif (areaname == 'NYS'): # 2020091500
        minlat=40.;maxlat=50.;minlon=270.;maxlon=290.
    elif (areaname == 'CONUS'): # 2020091500
        minlat=28.; maxlat=52.; minlon=235.; maxlon=285.
    else:
        print('Unknown area! Available lists: '
              'Glb | NH | SH | NAfr | TAfr  | Asia | EAsia|'
              'Indo| SPacAtl | SAmer| NAmer |'
              'NPO | NML| TRO| SML  | SPO   | GLK  | EUO  |'
              'NEU | Tri| Pir| Indi | NMLAtl|NMLPac|SMLInd|'
              'SMLPac   |TIWA|NWPac | Phil  | HudB | SF1~5|'
              'ST1~4| r2o1~10')
        pass
    
    if (areaname!='Glb' and minlon > maxlon):
        crosszero=True
        minlon=minlon-360.
    else:
        crosszero=False

    if (minlon==0. and maxlon==360.):
       cyclic=True
    else:
       cyclic=False
        
    return minlon,maxlon,minlat,maxlat,crosszero,cyclic
