import numpy as np
import matplotlib
import matplotlib.pyplot as plt
# import netCDF4 as nc
import datetime
import time
import os
import sys
import csv
import glob
import seawater as sw
from scipy.signal import find_peaks

def alamo_qc(fpin):
    ### Alamo QC: Apply the prescribed QC test for Alamo profiles.
    # # The criteria follows the standard Argo float
    # Note: look at on profile with the same time at a time
    
    # Input
    # fpin - float dictionary
    
    # Ouput
    # fpp - fpin, but with QC record fields updated
    # Author: Nan-Hsun Chi Aug 2023

    # del fp
    fp = fpin
    # Note: Note: Flags are set according to Ref Table 2 and sec 2.1 of the manual.
    # The flag value is not allowed to be reduced - eg if already set to 4, must
    # not override to 3. This is implemented in the algorithms below.

    # Work through this specific profile
    # Initialise QC variables where needed:
    #  0 = no QC done
    #  1 = good value
    #  6 = air temperature
    #  9 = missing value
    # first, get trap for missing profiles:

    if 'pos_qc' not in fp:
        fp['pos_qc'] = np.zeros( (len(fp['lat']),), dtype=np.int8 )
    if len(fp['pos_qc']) < len(fp['lat']):
        fp['pos_qc'] = np.zeros( (len(fp['lat']),), dtype=np.int8 )
    # if np.any(fp['p_raw']): # np.array returns false if the array is empty else it returns true
    if len(fp['p_raw']) > 0:
        fp['p_qc'] = np.ones( (len(fp['p_raw']),), dtype=np.int16 )
        jj = np.isnan( fp['p_raw'] )
        fp['p_qc'][jj] = 9
    if len(fp['t_raw']) > 0:
        fp['t_qc'] = np.ones( (len(fp['t_raw']),), dtype=np.int16 )
        jj = np.isnan( fp['t_raw'] )
        fp['t_qc'][jj] = 9
    if len(fp['s_raw']) > 0:
        fp['s_qc'] = np.ones( (len(fp['s_raw']),), dtype=np.int16 )
        jj = np.isnan( fp['s_raw'] )
        fp['s_qc'][jj] = 9
    
    #set all tests to zeros before starting
    fp['testsperformed'] = np.zeros( (16,) )
    fp['testsfailed'] = np.zeros( (16,) )

    nlev = len(fp['p_raw'])

    # Test 1:  Platform Identification
    # Because of the way we check our platforms, this will always be OK
    fp['testsperformed'][0] = 1

    # Test 2: Impossible Date Test:
    fp['testsperformed'][1] = 1
    #check date is between
    jref = datetime.datetime(1950,1,1)
    j1 = (datetime.datetime(1997,1,1)-jref).total_seconds()/86400 # days after jref
    j2 = (datetime.datetime.now()-jref).total_seconds()/86400 # days after jref
    # print(j1, j2)

    ibad = (fp['jday']< j1) | (fp['jday'] > j2) | (~np.isnan(fp['jday']))
    if np.any(ibad): # np.array returns false if the array is empty else it returns true
        fp['testsfailed'][1] = 1
        # fp['jday_qc(ibad)'] = 3 ### do not see this later
    
    # position tests, first flag all as ones, then missing locations
    jj = np.isnan(fp['lat'])
    if np.any(jj):
        fp['pos_qc'][jj] = 9

    # Test 3: Impossible Location Test:
    # We have done this test earlier
    fp['testsperformed'][2] = 1
    jj = ( (fp['lat'] < -90) | (fp['lat']>90) | (fp['lon']< 0) | (fp['lon'] > 360) )
    
    if np.any(jj):
        if fp['pos_qc'] != 8: # interpolated <-- This is directly translated from matlab. This line will have problem if any(jj) is true
            fp['testsfailed'][2] = 1
            fp['pos_qc'][jj] = 4
    jj = np.isnan(fp['lat'])
    if np.any(jj):
        fp['pos_qc'][jj] = 9
    
    # Test 4: Position on Land Test:
    # I did not do it since the trajectory of Alamo float is pretty obvious

    # Test 5: Impossible Speed Test: 
    # Do not use for Argo floats, use test 20 instead. 
    # Test speed between profiles. If it looks apparently wrong,
    # try some dvariant tests and maybe remove our present 1st fix if it appears wrong. 
    # Could test more combinations of previous profiles and fix numbers, 
    # but probably best to just eyball any cases where this test fails. 
    # Only look at previous profile, don't go backwards to others.
    fp['testsperformed'][4] = 0 # can't perform if current/ previous positions are not available. 
    # Now set the pos_qc flags that have passed the tests to 1
    if np.any( fp['testsperformed'][1:4]==1 ):
        jj = (fp['pos_qc'] == 0)
        fp['pos_qc'][jj] = 1
    
    # Test 6: Global Range Test:
    fp['testsperformed'][5] = 1

    ip = np.where( fp['p_raw'] < -5 )[0]
    jj = np.where( ((fp['t_raw'] <= -2.5) | (fp['t_raw'] > 40)) )[0]
    kk = np.where( ((fp['s_raw'] < 25.0) | (fp['s_raw'] > 41.)) )[0] ### adjusted smin for Arctic
    if len(ip) > 0:
        newv = np.repeat(4,len(ip));
        fp['p_qc'][ip] = 4#np.max( np.concatenate( (fp['p_qc'][ip], newv) ) )
        fp['testsfailed'][5] = 1
    if len(jj) > 0:
        newv = np.repeat(4,len(jj))
        fp['t_qc'][jj] = 4#np.max( np.concatenate( (fp['t_qc'][jj], newv) ) )
        fp['testsfailed'][5] = 1
    if len(kk) > 0:
        newv = np.repeat(4,len(kk))
        fp['s_qc'][kk] = 4#np.max( np.concatenate( (fp['s_qc'][kk], newv)) );
        fp['testsfailed'][5] = 1
    
    # Test 7: Regional Parameter Test
    # I did this in Test 6.

    # Test 8: Pressure Increasing Test
    fp['testsperformed'][7] = 1
    gg = np.where( ~np.isnan(fp['p_raw']) )[0]

    if np.any( np.diff(fp['p_raw'][gg])==0 ):
        fp['testsfailed'][7] = 1
        jj = ( np.diff(fp['p_raw'][gg])==0 )
        newv = np.repeat(4,jj.sum())
        if len(newv) > 0:
            fp['p_qc'][jj] = 4#np.max( np.concatenate( (fp['p_qc'][jj], newv) ) )
            fp['t_qc'][jj] = 4#np.max( np.concatenate( (fp['t_qc'][jj], newv) ) )
            fp['s_qc'][jj] = 4#np.max( np.concatenate( (fp['s_qc'][jj], newv) ) )

    # new process from here: Might need some updating at some stage, seems a bit clunky.
    bb = []
    kk = np.where( np.diff(fp['p_raw']) <0 )[0]
    if len(kk) > 0:
        for jj in range( len(kk) ):
            for l in range( kk[jj],kk[jj]+1 ):    #max(2,kk(jj)):min(length(fp.p_raw)-2,kk(jj)+1)
                if l >= len(fp['p_raw'])-1:
                    bb.extend( [np.min( [len(fp['p_raw']),l+1] )] )
                elif l == 1:
                    if fp['p_raw'][l] < fp['p_raw'][l+2]:
                        bb.extend([l])
                    else:
                        bb.extend([l+1])
                elif (fp['p_raw'][l] >= fp['p_raw'][l-1]) | (fp['p_raw'][l] <= fp['p_raw'][l+2]):
                    bb.extend([l])
        newv = np.repeat(4,len(bb))
        fp['s_qc'][bb] = np.max( np.concatenate( (fp['s_qc'][bb], newv) ) )
        fp['t_qc'][bb] = np.max( np.concatenate( (fp['t_qc'][bb], newv) ) )
        fp['p_qc'][bb] = np.max( np.concatenate( (fp['p_qc'][bb], newv) ) )

    # Test 9: Spike Test
    # testv is distance of v(n) outside the range of values v(n+1) and v(n-1).
    # If -ve, v(n) is inside the range of those adjacent points.
    fp['testsperformed'][8] = 1
    # PosPeaks,_ = find_peaks(fp['t_raw'],threshold=0.5)
    # NegPeaks,_ = find_peaks(-fp['t_raw'],threshold=0.5)
    # bdt = np.concatenate( ( PosPeaks,NegPeaks ) )
    bdt = findspike(fp['t_raw'],fp['p_raw'],'t')
    if len(bdt) > 0:
        newv = np.repeat(4,len(bdt))
        fp['t_qc'][bdt] = 4#np.max( np.concatenate( (fp['t_qc'][bdt], newv) ) )
        fp['s_qc'][bdt] = 4#np.max( np.concatenate( (fp['s_qc'][bdt], newv) ) )
        fp['testsfailed'][8] = 1
        
    bds = findspike(fp['s_raw'],fp['p_raw'],'s')
    if len(bds) > 0:
        # print(bds)
        newv = np.repeat(4,len(bds))
        fp['s_qc'][bds] = 4#np.max( np.concatenate( (fp['s_qc'][bds], newv) ) )
        fp['testsfailed'][8] = 1

    # Test 10: Top and Bottom Spike Test
    # Argo Quality Control Manual V2.8 (Jan 3, 2013) states
    # that this test is obsolete

    # Test 11: Gradient Test
    if nlev >= 3: #nlev is the number of p_raw values
        fp['testsperformed'][10] = 1
        jj = np.arange(1,nlev-1) #2:(nlev-1);        
        testv = np.abs(fp['t_raw'][jj] - (fp['t_raw'][jj+1]+fp['t_raw'][jj-1])/2);
        kk = np.where( (testv > 9) | ((fp['p_raw'][jj]>500) & (testv>3)))[0]
        if len(kk) > 0:
            newv = np.repeat(4,len(kk))
            fp['t_qc'][kk+1] = 4#np.max( np.concatenate( (fp['t_qc'][kk+1], newv) ) )
            fp['s_qc'][kk+1] = 4#np.max( np.concatenate( (fp['s_qc'][kk+1], newv) ) )
            fp['testsfailed'][10] = 1
        testv = np.abs(fp['s_raw'][jj] - (fp['s_raw'][jj+1]+fp['s_raw'][jj-1])/2)
        kk = np.where( (testv > 1.5) | ((fp['p_raw'][jj]>500) & (testv>0.5)) )[0]
        if len(kk) > 0:
            newv = np.repeat(4,len(kk))
            fp['s_qc'][kk+1] = 4#np.max( np.concatenate( (fp['s_qc'][kk+1], newv) ) )
            fp['testsfailed'][10] = 1
    
    # Test 12: Digit Rollover Test
    if len(fp['t_raw']) > 0:
        jj = np.where( np.diff(fp['t_raw']) > 10. )[0]
        fp['testsperformed'][11] = 1
        if len(jj) > 0:
            newv = np.repeat(4,len(jj))
            fp['t_qc'][jj+1] = 4#np.max( np.concatenate( (fp['t_qc'][jj+1], newv) ) )
            fp['testsfailed'][11] = 1
        
    if len(fp['s_raw']) > 0:
        fp['testsperformed'][11] = 1
        kk = np.where( np.diff(fp['s_raw'])>5. )[0]
        if len(kk):
            newv = np.repeat(4,len(kk))
            fp['s_qc'][kk+1] = 4#np.max( np.concatenate( (fp['s_qc'][kk+1], newv) ) )
            fp['testsfailed'][11] = 1

    # Test 13: Stuck Value Test
    # Test to see all measurements are the same value. If so, the profile is bad
    flds = ['s_raw','t_raw']
    fldsq = ['s_qc','t_qc']
    fp['testsperformed'][12] = 1
    for a in range( len(flds) ):
        if flds[a] in fp:
            if ( len(fp[flds[a]]) > 0 ) & all( fp[flds[a]]==fp[flds[a]][0] ):
                newv = np.repeat( 4,len(fp[flds[a]]) )
                # print('s_qc:',fp['s_qc'], 'newv:',newv)
                # fp[fldsq[a]] = np.max( np.concatenate( (fp['s_qc'], newv) ) )
                fp[fldsq[a]] = 4#np.max( np.append( newv,fp['s_qc']) )
                fp['testsfailed'][12] = 1

    # Test 14: Salinity Gradient Test
    # Check only from 10-m to MXL, see if gradient of salinity > 0.1 psu/m
    if ( len(fp['s_raw']) > 1 ) & ( len(fp['t_raw']) > 1) & ( len(fp['p_raw'])> 1 ):
        fp['testsperformed'][13] = 1
        ds = np.gradient(fp['s_raw'],fp['p_raw'])
        den0 = sw.eos80.pden(fp['s_raw'],fp['t_raw'],fp['p_raw'],0)
        if np.sum( ~np.isnan( den0[ (fp['p_raw']>3) & (fp['p_raw']<10) ] ) ) > 0:
            # print(den0[ (fp['p_raw']>3) & (fp['p_raw']<10) ])
            sfcd = np.nanmean( den0[ (fp['p_raw']>3) & (fp['p_raw']<10) ] )
            diffd = den0-sfcd
            p_diffd = fp['p_raw'][ (diffd<0.1) & (fp['p_raw']>10) ]
            if len(p_diffd) > 1:
                mxld = np.nanmin(p_diffd)
                
                jj = np.where( (np.abs(ds)>0.1) & (fp['p_raw'] > 10) & (fp['p_raw'] < mxld) )[0]
                if len(jj) > 0:
                    newv = np.repeat(4,len(jj))
                    fp['s_qc'][jj] = 4#np.max( np.concatenate( (fp['s_qc'][jj], newv) ) )
                    fp['s_qc'][ np.max([jj-1,0]) ] = 4#np.max( np.concatenate( (fp['s_qc'][np.max([jj-1,0])], newv) ) )
                    fp['s_qc'][ np.min([jj+1,len(fp['p_raw'])-1]) ] = 4#np.max( [fp['s_qc'][ np.min([jj+1,len(fp['p_raw'])-1]) ], newv] )
                    fp['testsfailed'][13] = 1
    
    # Test 15: Density Gradient Test
    # This test is operated to test whether the whole gradient of density in
    # the upper 100-m depth increasing with depth
    # Additional criteria is used if any of the density gradient more than -0.03 kg/m^4
    if ( len(fp['s_raw']) > 1 ) & ( len(fp['t_raw']) > 1 ) & ( len(fp['p_raw']) > 1 ):
        fp['testsperformed'][14] = 1
        den0 = sw.eos80.pden(fp['s_raw'],fp['t_raw'],fp['p_raw'],0);
        # ds = np.gradient(den0,fp['p_raw'])
        ### from top to bottom
        kz = np.where( (fp['p_raw']<100) & (fp['p_raw']>10) & (fp['s_qc']==1) & (fp['t_qc']==1) & (fp['p_qc']==1) )[0]
        if np.sum(kz) > 0:
            ps = np.diff( den0[kz] )
            jj = np.where( ps < -0.03 )[0]
            if len(jj) > 0:
                fp['s_qc'][kz[jj]] = 4
                fp['t_qc'][kz[jj]] = 4
                fp['p_qc'][kz[jj]] = 4
                fp['testsfailed'][14] = 1
        ### from bottom to top
        kz = np.where( (fp['p_raw']<100) & (fp['p_raw']>10) & (fp['s_qc']==1) & (fp['t_qc']==1) & (fp['p_qc']==1) )[0]
        if np.sum(kz) > 0:
            ps = np.diff( den0[kz] )
            jj = np.where( ps < -0.03 )[0]
            if len(jj) > 0:
                fp['s_qc'][kz[jj+1]] = 4
                fp['t_qc'][kz[jj+1]] = 4
                fp['p_qc'][kz[jj+1]] = 4
                fp['testsfailed'][14] = 1
    # Test 16: Air temperature Test
    # Flag the air temperature point in the profile
    if (len(fp['s_raw']) > 0) & (len(fp['t_raw']) > 0) & (len(fp['p_raw']) > 0 ):
        fp['testsperformed'][15] = 1
        jj = np.where( (fp['s_raw']<1) & (fp['p_raw']<1) )[0]
        if len(jj) > 0:
            newv = np.repeat(6,len(jj))
            fp['s_qc'][jj] = newv
            fp['t_qc'][jj] = newv
            fp['p_qc'][jj] = newv
        # else:
            fp['testsfailed'][15] = 1
             
    return fp

def findspike(values, pres, VarName):
    # VarName: 't' for temperature. 's' for salinity
    # values: temperature or salinity profile values
    # pres: pressure profile values
    spikes = np.array([])
    if VarName not in ['t','s']:
        print('invalid input variable name')
        return np.array([])
    else:
        for i in range( 1,len(pres)-2 ):
            test_value = np.abs(values[i]-(values[i+1]+values[i-1])/2) - \
                np.abs(values[i+1]-values[i-1]/2)
            if VarName == 't':
                if pres[i] < 500:
                    if test_value > 6:
                        spikes = np.append(spikes, i)
                else: # pres> 500 dbar
                    if test_value > 2:
                        spikes = np.append(spikes, i)
            if VarName == 's':
                if pres[i] < 500:
                    if test_value > 0.9:
                        spikes = np.append(spikes, i)
                else: # pres> 500 dbar
                    if test_value > 0.3:
                        spikes = np.append(spikes, i)
        spikes = spikes.astype(int)
        return spikes