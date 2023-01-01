from datetime import datetime, timedelta
import numpy as np
from PyTLE import TLE
from PyTLE import julian

class tle_fitter( TLE ):
    def __init__(self, *args, **kwargs):
        super().__init__( *args, **kwargs )
        self.allfields     = ['mean_motion','eccentricity','inclination','argp','raan','mean_anomaly',
                               'mm_dot','mm_dot_dot','bstar']

        self.default_fields = ['mean_motion','eccentricity','inclination','argp','raan','mean_anomaly',
                               'mm_dot','mm_dot_dot','bstar']

        if 'fields' in kwargs: self._optfields = kwargs['fields']
        else: self._optfields = self.default_fields

    def _setup_mappers(self):
        self._tovec = {
            'mean_motion'   : lambda : np.interp( self.mean_motion, [0,20], [0,1] ),
            'eccentricity'  : lambda : np.interp( self.eccentricity, [1e-15,1], [0,1]),
            'inclination'   : lambda : np.interp( self.inclination, [0,360], [0,1] ),
            'argp'          : lambda : np.interp( self.argp, [0,360], [0,1]),
            'raan'          : lambda : np.interp( self.raan, [0,360], [0,1]),
            'mean_anomaly'  : lambda : np.interp( self.mean_anomaly, [0,360], [0,1]),
            'mm_dot'        : lambda : np.interp( self.mm_dot, [-1,1], [0,1]),
            'mm_dot_dot'    : lambda : np.interp( self.mm_dot_dot, [-1,1], [0,1]),
            'bstar'         : lambda : np.interp(self.bstar, [-1, 1], [0, 1])
        }

        self._fromvec = {
            'mean_motion'   : lambda X: np.interp( X % 1, [0,1], [0,20] ),
            'eccentricity'  : lambda X: np.interp( X % 1, [0,1], [1e-15,1]),
            'inclination'   : lambda X: np.interp( X % 1, [0,1], [0,360]),
            'argp'          : lambda X: np.interp( X % 1, [0,1], [0,360]),
            'raan'          : lambda X: np.interp( X % 1, [0,1], [0,360]),
            'mean_anomaly'  : lambda X: np.interp( X % 1, [0,1], [0,360]),
            'mm_dot'        : lambda X: np.interp( X % 1, [0,1], [-1,1]),
            'mm_dot_dot'    : lambda X: np.interp( X % 1, [0,1], [-1,1]),
            'bstar'         : lambda X: np.interp( X % 1, [0,1], [-1,1])
        }

    def to_array(self):
        out = np.zeros( len(self._optfields) )
        for i,f in enumerate(self._optfields): out[i] = getattr( self, f )
        return out

    def from_array(self, vec):
        assert len(vec) == len(self._optfields)
        for i,f in enumerate(self._optfields): setattr(self,f,vec[i])
        return self


    def ephem_fit( self, jds, ephem_matrix ):
        ''' 
        you must initialize this TLE first or it will seed with default values
        ephem_matrix is : jd, temex, temey, temez, temedx, temedy, temedz...
        '''
        # do the imports here so we can use the class as just a parser without these dependencies (if we want to)
        import scipy.optimize
        from sgp4.earth_gravity import wgs72
        from sgp4.io import twoline2rv
        from sgp4.propagation import sgp4 as sgprop

        offsets     = 1440 * (jds - julian.to_jd(self.epoch) )
        
        # opt closure
        def fit_fcn( X, offset_mins, true_eph ):
            self.from_array(X)
            try: 
                testtle = twoline2rv( self.line1, self.line2, wgs72 )
                testeph = np.vstack( [np.hstack( sgprop(testtle,D)) for D in offset_mins] )
            except Exception as e:
                print(e)
                return np.inf
            
            diffpos = testeph[:,0:3] - true_eph[:,0:3]
            # return the rms value
            rmsval = np.sqrt( np.mean( np.linalg.norm(diffpos,axis=1) ** 2 ) )
            print(rmsval,end='\r')
            return rmsval


        opt_result = scipy.optimize.minimize( 
                        fit_fcn, 
                        self.to_array(), 
                        args=( offsets, ephem_matrix) , 
                        method='Nelder-Mead',
                        options = { 
                        'fatol' : 1. } )
        print()
        if opt_result.success:
            self.from_array( opt_result.x )
            return self, opt_result

        

# =====================================================================================================
if __name__ == '__main__':
    from sgp4.earth_gravity import wgs72
    from sgp4.io import twoline2rv
    from sgp4.propagation import sgp4 as sgprop

    # Aerocube 12A
    L1 = '1 43556U 18046C   22321.55519027  .00025005  00000+0  49749-3 0  9993'
    L2 = '2 43556  51.6329 154.1269 0008144 222.8163 137.2191 15.46745497242947'
    # test ephemeris
    tle = twoline2rv( L1, L2, wgs72 )
    print('tle epoch is ', julian.from_jd( tle.jdsatepoch ) )
    tledate = julian.from_jd( tle.jdsatepoch )
    
    mins   = np.arange(0,1440*5,10)
    jdates = tle.jdsatepoch + mins/1440
    eph = np.vstack( [np.hstack( sgprop(tle,D)) for D in mins ] )

    # init as if we only had a state-vector
    FIT = tle_fitter.fromPV( eph[0,0:3], eph[0,3:], tledate )
    print('Init:', FIT )

    # fit
    print('Fitting a new TLE to the test ephem')
    newtle, res = FIT.ephem_fit( jdates, eph )
    print(newtle.line1)
    print(L1)
    print(newtle.line2)
    print(L2)

    print()
    print(L1)
    print(L2)
    print(newtle)


