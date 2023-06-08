# https://github.com/dannyzed/julian/blob/master/julian/julian.py
from datetime import datetime
import math
import unittest

def to_jd(dt: datetime) -> float:
    """
    Converts a given datetime object to Julian date.
    Algorithm is copied from https://en.wikipedia.org/wiki/Julian_day
    All variable names are consistent with the notation on the wiki page.
    Parameters
    ----------
    fmt
    dt: datetime
        Datetime object to convert to MJD
    Returns
    -------
    jd: float
    """
    a = math.floor((14-dt.month)/12)
    y = dt.year + 4800 - a
    m = dt.month + 12*a - 3
    jdn = dt.day + math.floor((153*m + 2)/5) + 365*y + math.floor(y/4) - math.floor(y/100) + math.floor(y/400) - 32045
    jd = jdn + (dt.hour - 12) / 24 + dt.minute / 1440 + dt.second / 86400 + dt.microsecond / 86400000000
    return jd


def from_jd(jd: float) -> datetime:
    """
    Converts a Julian Date to a datetime object.
    Algorithm is from Fliegel and van Flandern (1968)
    Parameters
    ----------
    jd: float
        Julian Date as type specified in the string fmt
    fmt: str
    Returns
    -------
    dt: datetime
    """
    jd,jdf =  math.floor(jd + 0.5), jd + 0.5 - math.floor(jd + 0.5)
    l = jd+68569
    n = 4*l//146097
    l = l-(146097*n+3)//4
    i = 4000*(l+1)//1461001
    l = l-1461*i//4+31
    j = 80*l//2447
    k = l-2447*j//80
    l = j//11
    j = j+2-12*l
    i = 100*(n-49)+i+l

    year = int(i)
    month = int(j)
    day = int(k)
    # in microseconds
    frac_component = int(jdf * (1e6*24*3600))
    hours = int(frac_component // (1e6*3600))
    frac_component -= hours * 1e6*3600
    minutes = int(frac_component // (1e6*60))
    frac_component -= minutes * 1e6*60
    seconds = int(frac_component // 1e6)
    frac_component -= seconds*1e6
    frac_component = int(frac_component)
    dt = datetime(year=year, month=month, day=day,
                  hour=hours, minute=minutes, second=seconds, microsecond=frac_component)
    return  dt

# -----------------------------------------------------------------------------------------------------
class Testing( unittest.TestCase ):
    @classmethod
    def setUpClass(self):
        import astropy.time
        import numpy as np
        self.np = np
        self._ref_jds = self.np.arange( astropy.time.Time('2000-01-01T00:00:00.000Z',format='isot').jd,
                            astropy.time.Time.now().jd, 0.1 )
        self._ref_dts = astropy.time.Time( self._ref_jds, format='jd').datetime

    def test_fromjd(self):
        my_dt = [ from_jd(D) for D in self._ref_jds ]
        diff_seconds = [ (ref - test).total_seconds()
                         for ref,test in zip( self._ref_dts, my_dt ) ]
        print('Max delta fromjd (seconds) {:10.8f}'.format( self.np.max(diff_seconds) ) )
        self.assertTrue( self.np.max(diff_seconds) < 1 )

    def test_tojd(self):
        my_jd = [ to_jd(D) for D in self._ref_dts ]
        diff_days = [ ref - test
                         for ref,test in zip( self._ref_jds, my_jd ) ]
        print( 'Max delta tojd (seconds) {:10.8f}'.format( self.np.max(diff_days) * 86400 ) )
        self.assertTrue( self.np.max( diff_days ) < 1/86400. )

# =====================================================================================================
if __name__ == '__main__' :
    unittest.main()
    # now = datetime.utcnow()
    # print(now)
    # print(to_jd(now))
    # print(from_jd(to_jd(now)))
