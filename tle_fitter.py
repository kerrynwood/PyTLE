# ###############################################################################
# MIT License
# 
# Copyright (c) 2023 Kerry Wood
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# ###############################################################################

from datetime import datetime, timedelta
import numpy as np
from PyTLE import TLE
from PyTLE.utils import julian

# common fields
MAP = [
        #( human name, data struct name, range, map-to-range)
        ('mean_motion',     '_mm'       ,[0,20]     ,[0,1] ),
        ('eccentricity',    '_ecc'      ,[1e-15,1]  ,[0,1] ),
        ('inclination',     '_incl'     ,[0,360]    ,[0,1] ),
        ('argp',            '_argp'     ,[0,360]    ,[0,1] ),
        ('raan',            '_raan'     ,[0,360]    ,[0,1] ),
        ('mean_anomaly',    '_ma'       ,[0,360]    ,[0,1] ),
        ('n_dot',           '_ndot'     ,[-1,1]     ,[0,1] ),
        ('n_dot_dot',       '_ndotdot'  ,[-1,1]     ,[0,1] )
        ]

# for type 0/2
MAP_T0 = MAP + [ ('Bstar','_bstar',[-1,1],[0,1] ) ]

# for type 4
MAP_T4 = MAP + [ 
        ('Bterm',           '_B',       [-1,1]      ,[0,1] ),
        ('AGOM',            '_agom',    [1e-15,100] ,[0,1] ) 
        ]


class tle_fitter( TLE ):
    '''
    convenience routines to map internal TLE fields to a range an optimizer can use (generally 0--1)
    '''
    def __init__(self, TLE : TLE ):
        self._tle = TLE

    def get_map( self ):
        if self._tle._type == 0 or self._tle._type == 2: return MAP_T0
        if self._tle._type == 4 : return MAP_T4

    def _val_to_mapval( self, M ):
        human, field, orig_range, new_range = M
        return np.interp( getattr(self._tle,field), orig_range, new_range )

    def to_array( self ):
        return [ self._val_to_mapval(X) for X in self.get_map() ]

    def from_array( self, array ):
        # assume that order is preserved
        MAP = self.get_map()
        for i,M in enumerate(MAP): 
            human, field, orig_range, new_range = M
            # map it back to the range
            setattr(self._tle,field,np.interp( array[i], new_range, orig_range ) )
        return self

    def testme( self, **kwargs ):
        print(kwargs)

    @staticmethod
    def parseLines( L1, L2 ):
        return tle_fitter( TLE.parseLines(L1,L2) )

    def __str__( self ):
        return "\n".join( self._tle.generateLines() )
    
    def generateLine1( self ):
        return self._tle.generateLine1() 

    def generateLine2( self ):
        return self._tle.generateLine2() 

def test() :
    from sgp4.earth_gravity import wgs72
    from sgp4.io import twoline2rv
    from sgp4.propagation import sgp4 as sgprop

    # Aerocube 12A
    L1 = '1 43556U 18046C   22321.55519027  .00025005  00000+0  49749-3 0  9993'
    L2 = '2 43556  51.6329 154.1269 0008144 222.8163 137.2191 15.46745497242947'

    # parse the TLE and wrap a fitter
    TEST = tle_fitter( TLE.parseLines( L1, L2 ) )
    print('Str of fitter')
    print(str(TEST))

    # get array and then parse back into structure
    arr = TEST.to_array()
    print('Array: {}'.format(arr))
    print('Re-injecting...')
    TEST.from_array( arr )
    print(str(TEST))



    # test ephemeris
    tle = twoline2rv( L1, L2, wgs72 )
    print()
    print('tle epoch is ', julian.from_jd( tle.jdsatepoch ) )
    tledate = julian.from_jd( tle.jdsatepoch )
    
    mins   = np.arange(0,1440*5,10)
    jdates = tle.jdsatepoch + mins/1440
    eph = np.vstack( [np.hstack( sgprop(tle,D)) for D in mins ] )

        # init as if we only had a state-vector
    FIT = tle_fitter(TLE.fromPV( epoch=tledate, P=eph[0,0:3], V=eph[0,3:] ) )
    print()
    print( str(FIT) )


# =====================================================================================================
if __name__ == '__main__':
    test()
