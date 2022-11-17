import numpy as np
import PyTLE
import scipy.optimize
from sgp4.io import twoline2rv
from sgp4.earth_gravity import wgs72
from sgp4.propagation import sgp4 as sgp4prop
from sgp4.ext import rv2coe


###########################################################################################
def setDegrees( V ): return ( V.to_value(u.deg) + 360 ) % 360.

###########################################################################################
###########################################################################################
class tle_fitter( PyTLE.TLE ):
    '''
    -  this is a helper class to allow transforming orbits to and from Poliastro and TLE
    -  it also has routines to turn TLE into arrays for optimization
    -  set `override_bstar` to 0 to keep bstar from fitting
    '''
    def __init__( self,  L1=None, L2=None, cat='NA', override_bstar = None ): 
        super( tle_fitter, self).__init__( L1, L2, cat )
        self.fit_range    = [0,1]
        self.circle_range = [0,360]
        self.ecc_range    = [1e-15,1]
        self.mm_range     = [0, 19]
        self.override_bstar = override_bstar

    def fromPoliastro( self, pol ):
        a,ecc,inc,raan,argp,nu = pol.classical()
        # NOTE that these values have units in the AstroPy style
        self.set_date( pol.epoch.datetime ) # set the epoch
        self.inclination = setDegrees( inc )
        self.eccentricity = ecc.to_value(u.one)
        self.argument = setDegrees(argp)
        self.mean_motion = pol.n.value * 86400. / (2*np.pi) 
        self.ra = setDegrees(raan)
        self.mean_anomaly = setDegrees(nu)
        self.bstar = 0.

    def toArray( self ):
        '''
        return an array that represents the TLE orbital parameters
        mean_motion
        ecc
        inclination
        argument
        raan
        mean_anomaly
        bstar  <--- return in raw form
        '''
        # these values should never be negative or wrapped
        if self.override_bstar == None: bstar = self.bstar
        else: bstar = self.override_bstar
        return np.array( [np.interp( self.mean_motion, self.mm_range, self.fit_range ),
                np.interp( self.eccentricity, self.ecc_range, self.fit_range ),
                np.interp( self.inclination, self.circle_range, self.fit_range),
                np.interp( self.argument, self.circle_range, self.fit_range),
                np.interp( self.ra, self.circle_range, self.fit_range),
                np.interp( self.mean_anomaly, self.circle_range, self.fit_range),
                bstar ] )

    def circ( self, val ):
        return (val + self.fit_range[-1]) % self.fit_range[-1]

    def fromArray( self, X ):
        '''
        take the array mapped from toArray and turn it back into this data structure
        mean_motion
        ecc
        inclination
        argument
        raan
        mean_anomaly
        bstar  <--- return in raw form
        '''
        # all these should be circular
        X3 = self.circ(X[3]) # argument
        X4 = self.circ(X[4]) # ra
        X5 = self.circ(X[5]) # mean_anomaly
        self.mean_motion   = np.interp( X[0], self.fit_range, self.mm_range )
        self.eccentricity  = np.interp( X[1], self.fit_range, self.ecc_range )
        self.inclination   = np.interp( X[2], self.fit_range, self.circle_range )
        self.argument      = np.interp( X[3], self.fit_range, self.circle_range )
        self.ra            = np.interp( X[4], self.fit_range, self.circle_range )
        self.mean_anomaly  = np.interp( X[5], self.fit_range, self.circle_range )
        if self.override_bstar != None:  self.bstar = self.override_bstar
        else:                            self.bstar         = X[6]
        return self

    def fromPV( self, P, V, epoch=None) :
        '''
        fromPV : given state position and velocity (in TEME), build an initial TLE
        note   : this is *not* going to build mean elements
        '''
        # return p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper
        p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper = rv2coe( P, V, wgs72.mu )
        self.inclination = np.degrees( incl  )
        self.eccentricity = ecc
        self.argument = np.degrees( argp )
        self.ra = np.degrees( omega )
        self.mean_anomaly = np.degrees( m )
        self.bstar = 0
        # calculate the mean motion (these are km values, so we get rads/s, convert to TLE units)
        mm = np.sqrt( wgs72.mu / a**3 )
        self.mean_motion = ( mm * 86400 ) / (2*np.pi)
        self.bstar = 0.
        if epoch != None: self.set_date( epoch )
        return self
 

###########################################################################################
def rmsval( V ):
    try: return np.sqrt( np.mean( V ** 2 ) ) * 1e3
    except: return np.inf

###########################################################################################
def compare_TLE_ephem( dates, ephem, tle_fitter, rms = False ):
    '''
    dates : a list of astropy.time objects
    ephem : a list of states (TEME, km)
        dates and ephem must be the same length
    TLE   : (L1, L2) text lines

    returns a list of residuals or rms value (depending on flag)
    '''
    N          = len(dates)
    L1,L2      = tle_fitter.getLines()
    try: prop = twoline2rv( L1, L2, wgs72 )
    except:
        if rms: return np.inf(N)
        else: return np.inf
    tle_epoch  = prop.jdsatepoch
    offset     = np.array( [ D.jd - tle_epoch for D in dates ] ) * 1440.  # SGP4 wants minutes
    testpos, testvel  = zip( *[sgp4prop( prop, O ) for O in offset ] )
    rV = np.array( [ np.linalg.norm(X) for X in testpos - ephem ]  )
    print( "{:010.2f}".format(rmsval(rV)), end='\r')
    if rms: return rmsval( rV )
    return rV

###########################################################################################
def fit_fcn( X, dates, given_ephem, tle_fitter, rms=False ):
    tle_fitter.fromArray(X)
    return compare_TLE_ephem( dates, given_ephem, tle_fitter, rms=rms )

###########################################################################################
###########################################################################################
if __name__ == '__main__':
    import astropy.time
    import astropy.units as u
    from datetime import datetime, timedelta
    import numpy as np

    # test TLE
    L1='1 25544U 98067A   21007.32953112  .00000789  00000-0  22255-4 0  9990'
    L2='2 25544  51.6452  58.3686 0000693 182.6292 276.9175 15.49260299263666'
    # vanguard
    #L1='1     5U 58002B   20320.71165500 -.00000004 +00000-0 -15140-4 0  9993'
    #L2='2     5 034.2525 122.2981 1846950 313.0280 032.9388 10.84868096221685'
    # make some ephemeris
    sgpo    = twoline2rv( L1, L2, wgs72 )
    now     = astropy.time.Time( datetime.utcnow() )
    dates   = [now + (u.min*X) for X in range(0,1440*2,5)]
    offsets = [ (D-now).jd * 1440 for D in dates ]
    eph     = [ sgp4prop( sgpo, O ) for O in offsets ]
    P,V     = zip( *eph )
    P       = np.array(P)

    # return p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper
    p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper = rv2coe( P[0], V[0], wgs72.mu )
    T = tle_fitter()
    T.fromPV( P[0], V[0] )
    T.set_date( dates[0].datetime )
    print('-' * 80)
    print('The initial TLE was')
    print(L1)
    print(L2)
    print('-' * 80)

    print()
    print('-' * 80)
    print('The first guess TLE is')
    print(T)
    print('-' * 80)
    print()
    print('optimizing...')
    

    oval = scipy.optimize.minimize( fit_fcn,              # function
                                 T.toArray(),            # initial search argument
                                 method='Nelder-Mead',
                                 args=( dates, P, T, True),
                                 options = { 'disp':True, 'adaptive' : False, 'fatol':1 } )
    if oval.success:
        N = len(dates)
        Z = T.fromArray( oval.x )
        print()
        print('-' * 80)
        print('New TLE found:')
        print()
        print(Z)
        print()
        print('\t fit between {} and {}'.format( dates[0],dates[-1]) )
        print('\t RMS is     : {:7.3f}'.format(oval.fun))
        print('\t average is : {:7.3f}'.format( oval.fun / N ) )
    else:
        print('-------> Could not optimize')
