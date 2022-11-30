from datetime import datetime, timedelta
import numpy as np

from sgp4.earth_gravity import wgs72
from sgp4.ext import rv2coe
WGS72_MU = wgs72.mu

# alpha numeric lookups
from_alpha = {'A': 10, 'B': 11, 'C': 12, 'D': 13, 'E': 14, 'F': 15, 'G': 16, 'H': 17, 'J': 18, 'K': 19, 'L': 20, 'M': 21, 'N': 22, 'P': 23, 'Q': 24, 'R': 25, 'S': 26, 'T': 27, 'U': 28, 'V': 29, 'W': 30, 'X': 31, 'Y': 32, 'Z': 33}
to_alpha = {10: 'A', 11: 'B', 12: 'C', 13: 'D', 14: 'E', 15: 'F', 16: 'G', 17: 'H', 18: 'J', 19: 'K', 20: 'L', 21: 'M', 22: 'N', 23: 'P', 24: 'Q', 25: 'R', 26: 'S', 27: 'T', 28: 'U', 29: 'V', 30: 'W', 31: 'X', 32: 'Y', 33: 'Z'}

# -----------------------------------------------------------------------------------------------------
def alpha5(s):
    # from Brandon Rhodes' SGP4 library
    ''' compute an INTEGER from a TLE number string'''
    if not s[0].isalpha():
        return int(s)
    c = s[0]
    return (from_alpha[c] * 10000 ) + int(s[1:])

# -----------------------------------------------------------------------------------------------------
def integer5(I):
    ''' compute a string from an integer '''
    I = int(I)
    if I > 339999 : 
        raise Exception('cannot convert integers > 339999 to TLE strings')
    if I < 100000 : return str(I).zfill(5)
    intstr = str(I)
    lkup = intstr[0:2]
    return to_alpha[ int(lkup) ] + intstr[2:][0:5]

# -----------------------------------------------------------------------------------------------------
def generate_expo_format(flt):
    [mant, _, exp] = '{:+4.4e}'.format(flt).partition('e')
    mant = mant.replace('.', '')
    rV = '{:s}{:+1d}'.format(mant, int(exp) + 1)
    if flt == 0: return "+00000-0"
    return rV

# -----------------------------------------------------------------------------------------------------
# this takes the "00000-0" format as specified in TLE's and outputs a float
def process_expo_format(string):
    if string[0] == '-': neg = -1
    else: neg = 1
    mant = string[1:-2]
    exp = string[-2:]
    return neg * float("0.{}".format(mant)) * (10 ** int(exp))

# -----------------------------------------------------------------------------------------------------
def generate_checksum(line):
    digits = sum( [ int(X) for X in filter( lambda c :str.isdigit(c), line) ] )
    minus = line.count('-')
    rV = str(digits + minus)
    return rV[-1]

# -----------------------------------------------------------------------------------------------------
def revs_per_day_from_semimajor( a ):
    return 86400 * np.sqrt( WGS72_MU / (4 * np.pi ** 2 * a ** 3) )

# -----------------------------------------------------------------------------------------------------
class TLE:
    def __init__(self, l1=None, l2=None, vec=None, **kwargs):
        self._userfields = ['satnum','classification','elset_no','revnum']

        if l1 and l2: 
            self.parse_line1( l1 )
            self.parse_line2( l2 )
        else:
            self.parse_line1( '1 99999U 00000A   00001.00000000 +.00000000  00000-0 -00000-0 0  0001')
            self.parse_line2( '2 99999  10.0000 010.0000 0000000 010.0000 010.0000 01.00000000 00017')

    def _clean_line(self, l):
        while len(l) < 70: l+= '0'

    def parse_line1(self,l1):
        assert l1.startswith('1')
        self._clean_line(l1)
        self.satnum     = int(l1[2:7])
        self.classification    = l1[7]
        self.note       = l1[9:17]
        self.mm_dot     = float(l1[33:43])
        self.mm_dot_dot = process_expo_format(l1[44:52])
        self.bstar      = process_expo_format(l1[53:61])
        self.eph_type   = int( l1[62] )
        self.elset_no   = int( l1[64:68] )
        # self.checksum1  = int(l1[69])
        
        # generate epoch
        epoch_year = int(l1[18:20])
        if epoch_year <= 57: epoch_full_year = epoch_year + 2000
        else: epoch_full_year = 1900 + epoch_year
        epoch_day  = float(l1[20:32])
        self.epoch = datetime(epoch_full_year, 1, 1 ) + timedelta( days=epoch_day )

    def parse_line2(self, l2):
        assert l2.startswith('2')
        self._clean_line(l2)
        self.satnum = int(l2[2:7])
        self.inclination = float(l2[8:16])
        self.raan = float(l2[17:25])
        self.eccentricity = float('0.{}'.format(l2[26:33]))
        self.argp = float(l2[34:42])
        self.mean_anomaly = float(l2[43:51])
        self.mean_motion = float(l2[52:63])
        self.revnum = int(l2[63:68])
        # self.checksum2    = int( l2[69] )

    @property
    def epoch_day(self):
        return 1 + (self.epoch - datetime(self.epoch.year,1,1)).total_seconds() / 86400
        #return int( self.epoch.strftime('%j') )

    @property
    def epoch_year(self):
        return self.epoch.year

    @property
    def epoch_year_str(self):
        return self.epoch.strftime('%y')

    @property
    def seconds_per_orbit(self): return 86400. / self.mean_motion

    @property
    def orbit_mean_motion(self):
        return (np.pi * 2) / self.seconds_per_orbit

    @property
    def semimajor_axis(self):
        return ( WGS72_MU / self.orbit_mean_motion ** 2 ) ** (1/3.)

    @property
    def apogee(self):
        return self.semimajor_axis * ( 1 + self.eccentricity )

    @property
    def line1( self ): 
        return self._generate_new_line1()

    @property
    def line2( self ): 
        return self._generate_new_line2()


    def _generate_new_line1( self ):
        rV = '1 {}{:1} {:8} {:2}{:012.8f} '.format(
            integer5(self.satnum),
            self.classification[0],
            self.note[:8],
            self.epoch_year_str,
            self.epoch_day)
        rV += "{:+14.13f}".format( self.mm_dot ).replace('0.','.')[:10]
        rV += ' {}'.format( generate_expo_format( self.mm_dot_dot ) )
        rV += ' {}'.format(generate_expo_format( self.bstar ) )
        rV += ' {:1}'.format( self.eph_type )
        rV += ' {:4}'.format( self.elset_no )
        rV += generate_checksum(rV)
        return rV

    def _generate_new_line2( self ):
        rV = '2 {} {:08.4f} {:08.4f} '.format( integer5(self.satnum), self.inclination, self.raan)
        rV += '{:9.7f}'.format( self.eccentricity ).partition('.')[2]
        rV += ' {:08.4f} {:08.4f} {:011.8f}{:05d}'.format( self.argp, self.mean_anomaly, self.mean_motion, self.revnum )
        rV += generate_checksum(rV)
        return rV

    def getlines(self):
        return self._generate_new_line1(), self._generate_new_line2()

    @classmethod
    def fromPoliastro( cls, pol ):
        import astropy.units as u
        def setDegrees(V): return (V.to_value(u.deg) + 360) % 360.
        a,ecc,inc,raan,argp,nu = pol.classical()
        newcls = cls()
        # NOTE that these values have units in the AstroPy style
        newcls.set_date( pol.epoch.datetime ) # set the epoch
        newcls.inclination = setDegrees( inc )
        newcls.eccentricity = ecc.to_value(u.one)
        newcls.argument = setDegrees(argp)
        newcls.mean_motion = pol.n.value * 86400. / (2*np.pi)
        newcls.ra = setDegrees(raan)
        newcls.mean_anomaly = setDegrees(nu)
        newcls.bstar = 0.
        return newcls

    @classmethod
    def fromPV(cls, P, V, epoch : datetime, **kwargs):
        '''
        fromPV : given state position and velocity (in TEME), build an initial TLE
        note   : this is *not* going to build mean elements
        '''
        # return p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper
        newcls = cls()
        if V[2] == 0:
            raise Exception('cannot init an orbit with perfectly zero inclination (velocity[Z] ~ 1e-5km/s minimum)')
            return newcls
        try: p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper = rv2coe(P, V, wgs72.mu)
        except Exception as e:
            print('could not init TLE from P: {} V: {}'.format( P,V ) )
            return newcls
        # rv2coe in sgp4 returns > 999999 to indicate undefined or infinite... (see code)
        if any( (X > 999999. for X in p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ) ):
            print('rv2coe returned undefined (> 999999) value')
            return newcls
        newcls.inclination = np.degrees(incl)
        newcls.eccentricity = ecc
        newcls.argp = np.degrees(argp)
        newcls.raan = np.degrees(omega)
        newcls.mean_anomaly = np.degrees(m)
        newcls.bstar = 0
        # calculate the mean motion (these are km values, so we get rads/s, convert to TLE units)
        mm = np.sqrt(wgs72.mu / a ** 3)
        newcls.mean_motion = (mm * 86400) / (2 * np.pi)
        newcls.bstar = 0.
        newcls.epoch = epoch
        # override any parameters that are in our user fields list
        for k in kwargs:
            if k in self._userfields: setattr(self,k,kwargs[k] )
        return newcls

    def __str__(self): return "\n".join(self.getlines())
    def __repr__(self): return str(self)


# ========================================================================================================
if __name__=="__main__":
    L1 = '1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927'
    L2 = '2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537'

    Q = TLE.fromPV( [7000,0,0], [0,8,0.1], datetime.utcnow() )

    print(Q)



