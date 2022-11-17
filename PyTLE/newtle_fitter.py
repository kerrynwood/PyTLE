import numpy as np
from RefactorTLE_working import TLE


class tle_fitter( TLE ):
    def __init__(self, **kwargs):
        super().__init__( **kwargs )
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

    def to_vec(self):
        out = np.zeros( len(self._optfields) )
        for i,f in enumerate(self._optfields): out = getattr( self, f )
        return out

    def from_vec(self, vec):
        assert len(vec) == len(self._optfields)
        for i,f in enumerate(self._optfields): setattr(self,f,vec[i])
        return self

if __name__ == '__main__':
    A = tle_fitter(fields=['mean_motion'])
    print(A)
    print(A.to_vec())
    print(A.from_vec([10]))