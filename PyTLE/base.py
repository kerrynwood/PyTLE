from datetime import datetime, timedelta
import re
from alpha import alpha_to_integer, integer_to_alpha

# -----------------------------------------------------------------------------------------------------
def generate_expo_format(flt):
    [mant, crap, exp] = '{:+4.4e}'.format(flt).partition('e')
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
    #return neg * float('0.%s' % mant) * (10 ** int(exp))

# -----------------------------------------------------------------------------------------------------
def epoch_str_todatetime( S ):
    year = int(S[0:2])
    if year > 57: tyear = datetime( year=1900+year, month=1, day=1 )
    if year < 57: tyear = datetime( year=2000+year, month=1, day=1 )
    frac = float( S[2:] )
    return tyear + timedelta( days=frac )

# -----------------------------------------------------------------------------------------------------
def datetime_to_epochstr( dt ):
    tyear = datetime(year=dt.year, day=1, month=1 )
    frac = (dt - tyear).total_seconds() / 86400.
    year = dt.strftime('%y')
    ifrac = int(frac)
    lfrac = '{}'.format( frac - ifrac )
    return '{}{:03d}.{}'.format( year, ifrac, lfrac[2:] )[:14]


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
                self._class,
                self._intld,
                datetime_to_epochstr( self._epoch ),
                generate_expo_format( self.agom ),
                generate_expo_format( 0 ), #self.B ),
                self._elset )
        return L1


if __name__ == "__main__":
    L1='1 12345U xyzzyz   23038.45547454 +.00000000 +46171+0 +33000-1 4 99992'
    L2='2 12345   9.7332 113.4837 7006332 206.5371  38.9576 01.00149480000003'
    Q = TLE_4( L1, L2 )
    print(L1)
    print(Q.generateLine1())

