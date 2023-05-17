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

    def parseDate( self, S ):
        ''' assume this is just the date string '''
        self._epoch = epoch_str_todatetime( S )

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
        self.ndot    = float( S[33:43] )
        self.ndotdot = process_expo_format( S[44:52] )
        self.bstar  = process_expo_format( S[53:61] )
        self.type   = int(S[62])
        self._elset = int(S[64:68])
    
    def parseLine2( self, S ):
        if S[0] != '2' : raise Exception('LINE2 must begin with 2')
        if alpha_to_integer( S[2:7] ) != self._satno : raise Exception('satno does not match')
        self.incl = float(S[8:16])
        self.raan = float(S[17:25])
        self.ecc  = float( '0.{}'.format( S[26:33] ) )
        self.argp = float(S[34:42])
        self.ma   = float(S[43:51])
        self.mm   = float(S[52:63])

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
                "{:09.8f}".format( self.ndot ).ljust(9,' '),
                generate_expo_format( self.ndotdot),
                generate_expo_format( self.bstar ),
                0,  # <--------- TYPE FLAG
                "{}".format( self._elset ).rjust(4,' ')
                )

                #self._elset )
        return L1

    def generateLine2( self ):
        #L2='2 12345   9.7332 113.4837 7006332 206.5371  38.9576 01.00149480000003'
        L2 = '2 {:5} {:8} {:8} {:7} {:8} {:8} {} {}3 '.format(
                integer_to_alpha( self._satno ).rjust(5,'0'),
                "{:>8.4f}".format( self.incl )[:8], 
                "{:>8.4f}".format( self.raan)[:8], 
                "{}".format( self.ecc )[2:10],
                "{:>8.4f}".format( self.argp)[:8], 
                "{:>8.4f}".format( self.ma)[:8], 
                "{:>011.8f}".format( self.mm)[:12], 
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

        self.agom = process_expo_format( S[44:52] )
        self.B    = process_expo_format( S[53:61] )
        self.type = int(S[62])
        self._elset = int(S[64:68])
    
    def parseLine2( self, S ):
        if S[0] != '2' : raise Exception('LINE2 must begin with 2')
        if alpha_to_integer( S[2:7] ) != self._satno : raise Exception('satno does not match')
        self.incl = float(S[8:16])
        self.raan = float(S[17:25])
        self.ecc  = float( '0.{}'.format( S[26:33] ) )
        self.argp = float(S[34:42])
        self.ma   = float(S[43:51])
        self.mm   = float(S[52:63])

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
                generate_expo_format( self.agom ),
                generate_expo_format( self.B ),
                self._elset )
        return L1

    def generateLine2( self ):
        #L2='2 12345   9.7332 113.4837 7006332 206.5371  38.9576 01.00149480000003'
        L2 = '2 {:5} {:8} {:8} {:7} {:8} {:8} {} {}3 '.format(
                integer_to_alpha( self._satno ).rjust(5,'0'),
                "{:>8.4f}".format( self.incl )[:8], 
                "{:>8.4f}".format( self.raan)[:8], 
                "{}".format( self.ecc )[2:10],
                "{:>8.4f}".format( self.argp)[:8], 
                "{:>8.4f}".format( self.ma)[:8], 
                "{:>011.8f}".format( self.mm)[:12], 
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
    print(Z.bstar)
    print(L1)
    print(Z.generateLine1())
    print(L2)
    print(Z.generateLine2())
