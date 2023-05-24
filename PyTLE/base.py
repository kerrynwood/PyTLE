from datetime import datetime, timedelta
import re
from alpha import alpha_to_integer, integer_to_alpha

from formatters import generate_expo_format, process_expo_format
from formatters import epoch_str_todatetime, datetime_to_epochstr

# -----------------------------------------------------------------------------------------------------
class TLE:
    def __init__( self ):
        self.clear()

    def clear( self ):
        self._epoch = None
        self._line1 = None
        self._line2 = None
        self._satno = None
        self._class = 'U'
        self._intld = ''
        self._elset = None
        self._incl    = 0.
        self._raan    = 0.
        self._ecc     = 0.
        self._argp    = 0.
        self._ma      = 0.
        self._mm      = 0. 
        # type 2
        self._ndot    = 0.
        self._ndotdot = 0.
        self._type    = 0.
        self._bstar   = 0.
        # type 4
        self._B       = 0.
        self._agom    = 0.

        # human readable / helper
        self._perigee = None
        self._apogee  = None


    @property
    def perigee( self ):
        if self._perigee is None: self._calculate_apogee_perigee()
        return self._perigee
    
    @property
    def apogee( self ):
        if self._apogee is None: self._calculate_apogee_perigee()
        return self._apogee

    def parseDate( self, S ):
        ''' assume this is just the date string '''
        self._epoch = epoch_str_todatetime( S )

    @staticmethod 
    def parseLines( L1, L2 ):
        if L1[62] == '0' or L1[62] == '2' : return TLE_2( L1, L2 )
        if L1[62] == '4' : return TLE_4( L1, L2 )

    @staticmethod
    def get_type0( self ):
        return TLE_2()

    @staticmethod
    def get_type4( self ):
        return TLE_4()


    def _calculate_apogee_perigee( self, earth_rad = 6378.135):
        ''' Default value for earth_rad is taken from space-track.
        space-track : https://www.space-track.org/documentation#/faq
        Additional references: http://www.satobs.org/seesat/Dec-2002/0197.html
        '''
        semi_major = (8681663.653 / self.mean_motion) ** (2.0/3.0)
        self.perigee = ( semi_major * (1 - self.eccentricity) ) - earth_rad
        self.apogee =  ( semi_major * (1 + self.eccentricity) ) - earth_rad

    def _calculate_apogee_perigee( self, earth_rad = 6378.135):
        ''' Default value for earth_rad is taken from space-track.
        space-track : https://www.space-track.org/documentation#/faq
        Additional references: http://www.satobs.org/seesat/Dec-2002/0197.html
        '''
        semi_major = (8681663.653 / self.mm) ** (2.0/3.0)
        self._perigee = ( semi_major * (1 - self.eccentricity) ) - earth_rad
        self._apogee =  ( semi_major * (1 + self.eccentricity) ) - earth_rad


# -----------------------------------------------------------------------------------------------------
class TLE_2( TLE ):
    """
    1 25544U 98067A   23137.83559306  .00011914  00000-0  21418-3 0  9990
    2 25544  51.6409 118.9691 0006630 359.0829  72.4864 15.50282135397083
    """
    def __init__(self, L1=None, L2=None):
        self.clear()
        if L1 and L2: self.parseLines(L1,L2)

    def parseLine1( self, S ):
        if S[0] != '1' : raise Exception('LINE1 must begin with 1')
        self._satno = alpha_to_integer( S[2:7] )
        self._class = S[7]
        self._intld = S[9:17]
        self.parseDate( S[18:32] ) 
        self._ndot    = float( S[33:43] )
        self._ndotdot = process_expo_format( S[44:52] )
        self._bstar  = process_expo_format( S[53:61] )
        self._type   = int(S[62])
        self._elset = int(S[64:68])
    
    def parseLine2( self, S ):
        if S[0] != '2' : raise Exception('LINE2 must begin with 2')
        if alpha_to_integer( S[2:7] ) != self._satno : raise Exception('satno does not match')
        self._incl = float(S[8:16])
        self._raan = float(S[17:25])
        self._ecc  = float( '0.{}'.format( S[26:33] ) )
        self._argp = float(S[34:42])
        self._ma   = float(S[43:51])
        self._mm   = float(S[52:63])

    def parseLines( self, L1, L2 ):
        self.parseLine1( L1 )
        self.parseLine2( L2 )

    def generateLine1( self ):
        #1 25544U 98067A   23137.83559306  .00011914  00000-0  21418-3 0  9990
        L1 = '1 {:5}{:1} {:8} {:14} {} {} {} {} {}0'.format(
                integer_to_alpha( self._satno ).rjust(5,'0'),
                self._class[0],
                self._intld[:8].rjust(8,' '),
                datetime_to_epochstr( self._epoch ),
                "{:09.8f}".format( self._ndot ).ljust(9,' '),
                generate_expo_format( self._ndotdot),
                generate_expo_format( self._bstar ),
                0,  # <--------- TYPE FLAG
                "{}".format( self._elset ).rjust(4,' ')
                )

                #self._elset )
        return L1

    def generateLine2( self ):
        #L2='2 12345   9.7332 113.4837 7006332 206.5371  38.9576 01.00149480000003'
        L2 = '2 {:5} {:8} {:8} {:7} {:8} {:8} {} {}3 '.format(
                integer_to_alpha( self._satno ).rjust(5,'0'),
                "{:>8.4f}".format( self._incl )[:8], 
                "{:>8.4f}".format( self._raan)[:8], 
                "{}".format( self._ecc )[2:10],
                "{:>8.4f}".format( self._argp)[:8], 
                "{:>8.4f}".format( self._ma)[:8], 
                "{:>011.8f}".format( self._mm)[:12], 
                "{:04d}".format( self._elset)[:4]
                )
        return L2
    
    def generateLines( self ):
        return ( self.generateLine1(), self.generateLine2() )



# -----------------------------------------------------------------------------------------------------
class TLE_4( TLE ):
    def __init__(self, L1=None, L2=None):
        self.clear()
        if L1 and L2: self.parseLines(L1,L2)

    def parseLine1( self, S ):
        if S[0] != '1' : raise Exception('LINE1 must begin with 1')
        self._satno = alpha_to_integer( S[2:7] )
        self._class = S[7]
        self._intld = S[9:17]
        self.parseDate( S[18:32] ) 
        self._agom = process_expo_format( S[44:52] )
        self._B    = process_expo_format( S[53:61] )
        self._type = int(S[62])
        self._elset = int(S[64:68])
    
    def parseLine2( self, S ):
        if S[0] != '2' : raise Exception('LINE2 must begin with 2')
        if alpha_to_integer( S[2:7] ) != self._satno : raise Exception('satno does not match')
        self._incl = float(S[8:16])
        self._raan = float(S[17:25])
        self._ecc  = float( '0.{}'.format( S[26:33] ) )
        self._argp = float(S[34:42])
        self._ma   = float(S[43:51])
        self._mm   = float(S[52:63])

    def parseLines( self, L1, L2 ):
        self.parseLine1( L1 )
        self.parseLine2( L2 )

    def generateLine1( self ):
        #L1='1 12345U xyzzyz   23038.45547454 +.00000000 +46171+0 +33000-1 4 99992'
        L1 = '1 {:5}{:1} {:8} {:14} +.00000000 {} {} 4 {}0'.format(
                integer_to_alpha( self._satno ).rjust(5,'0'),
                self._class[0],
                self._intld[:8].rjust(8,' '),
                datetime_to_epochstr( self._epoch ),
                generate_expo_format( self._agom ),
                generate_expo_format( self._B ),
                self._elset )
        return L1

    def generateLine2( self ):
        #L2='2 12345   9.7332 113.4837 7006332 206.5371  38.9576 01.00149480000003'
        L2 = '2 {:5} {:8} {:8} {:7} {:8} {:8} {} {}3 '.format(
                integer_to_alpha( self._satno ).rjust(5,'0'),
                "{:>8.4f}".format( self._incl )[:8], 
                "{:>8.4f}".format( self._raan)[:8], 
                "{}".format( self._ecc )[2:10],
                "{:>8.4f}".format( self._argp)[:8], 
                "{:>8.4f}".format( self._ma)[:8], 
                "{:>011.8f}".format( self._mm)[:12], 
                "{:04d}".format( self._elset)[:4]
                )
        return L2
    
    def generateLines( self ):
        return ( self.generateLine1(), self.generateLine2() )



# =====================================================================================================
if __name__ == "__main__":
    L1='1 12345U xyzzyz   23038.45547454 +.00000000 +46171+0 +33000-1 4 99992'
    L2='2 12345   9.7332 113.4837 7006332 206.5371  38.9576 01.00149480000003'
    print('Original TLE')
    print(L1)
    print(L2)
    print('....parsing...')
    Q = TLE_4( L1, L2 )
    print('original   ', L1)
    print('generated  ',Q.generateLine1())
    print('original   ',L2)
    print('generated  ',Q.generateLine2())

    ISTR = 'hitherefromkerry'
    print('Setting international designator to {}'.format( ISTR ))
    Q._intld = 'hitherefromkerry'
    RAAN = 15000000
    print('Setting RAAN to {}'.format(RAAN))
    Q.raan = RAAN
    print("\n".join(Q.generateLines()))
    AGOM = 10e6
    print('Setting AGOM to {}'.format(AGOM))
    Q.agom = AGOM
    print("\n".join(Q.generateLines()))


    L1 = '1 25544U 98067A   23137.83559306  .00011914  00000-0  21418-3 0  9990'
    L2 = '2 25544  51.6409 118.9691 0006630 359.0829  72.4864 15.50282135397083'
    Z = TLE_2( L1, L2 )
    print(Z._bstar)
    print(L1)
    print(Z.generateLine1())
    print(L2)
    print(Z.generateLine2())
    Z._raan = 999999999
    Z._incl = 999999999
    print("\n".join(Z.generateLines()))

    Y = TLE.parseLines( L1, L2 )
    print(Y.__dict__)
