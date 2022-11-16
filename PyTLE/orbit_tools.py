from sgp4.io import twoline2rv
from sgp4.earth_gravity import wgs72
from sgp4.propagation import sgp4 as sgprop

import numpy as np
import scipy.interpolate

T = twoline2rv( '1 25544U 98067A   22319.90680528  .00010313  00000+0  18817-3 0  9990',
                '2 25544  51.6437 304.9926 0006899  77.7150   6.3579 15.50029035368715',
                wgs72)


offsets = np.arange(0,1440,5)
testoff = np.arange(2.5,1440,5)

eph = np.vstack( [np.hstack( sgprop(T,o) ) for o in offsets] )
inteph = scipy.interpolate.interp1d( offsets, eph.T, kind='cubic')
#inteph = scipy.interpolate.barycentric_interpolate( offsets, eph.T, testoff )


for to in testoff:
    interpolated = np.hstack(inteph(to))
    actual = np.hstack( sgprop(T,to) ) 
    print(interpolated-actual)
    print()
