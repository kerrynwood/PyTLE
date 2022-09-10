import re
import datetime
import math
import sys
from . import formatExceptionInfo
from . import tle_epoch
from astropy.time import Time as astro_date

# KNW changes for alphanumeric TLE
# alpha numeric lookups
from_alpha = {'A': 10, 'B': 11, 'C': 12, 'D': 13, 'E': 14, 'F': 15, 'G': 16, 'H': 17, 'J': 18, 'K': 19, 'L': 20, 'M': 21, 'N': 22, 'P': 23, 'Q': 24, 'R': 25, 'S': 26, 'T': 27, 'U': 28, 'V': 29, 'W': 30, 'X': 31, 'Y': 32, 'Z': 33}
to_alpha = {10: 'A', 11: 'B', 12: 'C', 13: 'D', 14: 'E', 15: 'F', 16: 'G', 17: 'H', 18: 'J', 19: 'K', 20: 'L', 21: 'M', 22: 'N', 23: 'P', 24: 'Q', 25: 'R', 26: 'S', 27: 'T', 28: 'U', 29: 'V', 30: 'W', 31: 'X', 32: 'Y', 33: 'Z'}

def alpha5(s):
    # from Brandon Rhodes' SGP4 library
    ''' compute an INTEGER from a TLE number string'''
    if not s[0].isalpha():
        return int(s)
    c = s[0]
    return (from_alpha[c] * 10000 ) + int(s[1:])

def integer5(I):
    ''' compute a string from an integer '''
    I = int(I)
    if I > 339999 : 
        raise Exception('cannot convert integers > 339999 to TLE strings')
    if I < 100000 : return str(I).zfill(5)
    intstr = str(I)
    lkup = intstr[0:2]
    return to_alpha[ int(lkup) ] + intstr[2:][0:5]


# NOTE! As of 2013/08/06, SpaceTrack is outputting records that are somewhat malformed.  Sometimes they have a non-standard NORAD identifier,
# and do not include a checksum.  Since we don't check checksum anyways,
# ignore it.  (Commented out in fields dictionary below).

# for SGP4 satrec
deg2rad  =   math.pi / 180.0         #    0.0174532925199433
xpdotp   =  1440.0 / (2.0 *math.pi)  #  229.1831180523293

# expoRE = re.compile("([\+\-\ ]?[0-9]{4,6})([\+\-]?[0-9]{1,2})")
launch_year_re = re.compile(r'(\d{2})\d{2,3} *')
launch_number_re = re.compile(r'\d{2}(\d{2,3}) *')
launch_piece_re = re.compile(r'\d{2}\d{2,3}.*([A-Z]{1,3})')

##########################################################################

def generate_checksum(line):
    digits = sum(map(lambda x: int(x), filter(
        lambda y: re.match('[0-9]{1}', y), [z for z in line])))
    minus = line.count('-')
    rV = str(digits + minus)
    return rV[-1]
##########################################################################
##########################################################################

def generate_expo_format(flt):
    [mant, crap, exp] = '{:+4.4e}'.format(flt).partition('e')
    mant = mant.replace('.', '')
    rV = '{:s}{:+1d}'.format(mant, int(exp) + 1)
    if flt == 0: return "+00000-0"
    return rV
##########################################################################
##########################################################################
# this takes the "00000-0" format as specified in TLE's and outputs a float

def process_expo_format(string):
    if string[0] == '-': neg = -1
    else: neg = 1
    mant = string[1:-2]
    exp = string[-2:]
    return neg * float("0.{}".format(mant)) * (10 ** int(exp))
    #return neg * float('0.%s' % mant) * (10 ** int(exp))


##########################################################################
catDate = re.compile("([0-9]{4})_([0-9]{3})")

def get_catalog_date(STR):
    R = catDate.search(STR)
    if not R or not len(R.groups()) == 2:
        print("get_catalog_date: could not pull info from", STR)
        return int(0)
    try: retVal = int(R.groups()[0] + R.groups()[1])
    except:
        print("get_catalog_date: could not reformat info:", Y, D)
        return int(0)
    return retVal


##########################################################################
def get_tle_datatype(VAL, TYPE):
    if TYPE == 'float':
        try: return float(VAL)
        except: print("Could not convert", VAL, "to type:", TYPE)

    if TYPE == 'int':
        try: return int(VAL)
        except: print("Could not convert", VAL, "to type:", TYPE)

    if TYPE == 'expo':
        return process_expo_format(VAL)
    # except: print "Could not convert", VAL, "to type:", TYPE

    if TYPE == 'string':
        try: return str(VAL)
        except: print("Could not convert", VAL, "to type:", TYPE)

    if TYPE == 'eccentricity':
        # put in the implicit floating point notation
            try: return float(str('0.%s' % VAL))
            except: print("Could not convert", VAL, "to type:", TYPE)

    if TYPE == "launch_year":
        try: return int(launch_year_re.match(VAL).groups()[0])
        except:
            print("Could not convert", VAL, "to launch_year")

    if TYPE == "launch_number":
        try: return int(launch_number_re.match(VAL).groups()[0])
        except: print("Could not convert", VAL, "to launch_number")

    if TYPE == "launch_piece":
        return VAL
    # try: return launch_piece_re.match( VAL ).groups()[0]
            # except: print "Could not convert", VAL, "to launch_piece"

    return None

##########################################################################
# list the fields that will populate the struct here.  This also lists the
# beginning spot (B), ending spot (E), and data type to feed into above.
tle_fields = {'line1': {
                'sat_no': {'B': 2, 'E': 7, 'T': 'alpha5', 'default' : 99999},
                'classification': {'B': 7, 'E': 8, 'T': 'string','default' : 'U'},
                'launch_year': {'B': 9, 'E': 17, 'T': 'launch_year', 'default' : 0 },
                'launch_number': {'B': 9, 'E': 17, 'T': 'launch_number', 'default' : 1},
                'launch_piece': {'B': 14, 'E': 17, 'T': 'launch_piece','default' : 'A'},
                'epoch_year': {'B': 18, 'E': 20, 'T': 'int','default' : 0},
                'epoch_day': {'B': 20, 'E': 32, 'T': 'float','default' : 1.},
                'mean_motion_1': {'B': 33, 'E': 43, 'T': 'float','default' : 0.},
                'mean_motion_2': {'B': 44, 'E': 52, 'T': 'expo','default' : 0.},
                'bstar': {'B': 53, 'E': 61, 'T': 'expo','default' : 0.},
                'ephem_type': {'B': 62, 'E': 63, 'T': 'int','default' : 0},
                'element_number': {'B': 64, 'E': 68, 'T': 'int','default' : 1},
                # 'line1_checksum' : {'B':68,'E':69,'T':'int'}
    },
            'line2': {
                'sat_no2': {'B': 2, 'E': 7, 'T': 'alpha5','default' : 99999},
                'inclination': {'B': 8, 'E': 16, 'T': 'float','default' : 25.},
                'ra': {'B': 17, 'E': 25, 'T': 'float','default' : 60.},
                'eccentricity': {'B': 26, 'E': 33, 'T': 'eccentricity','default' : 0.1},
                'argument': {'B': 34, 'E': 42, 'T': 'float','default' : 0.1},
                'mean_anomaly': {'B': 43, 'E': 51, 'T': 'float','default' : 0.1},
                'mean_motion': {'B': 52, 'E': 63, 'T': 'float','default' : 0.1},
                'rev_number': {'B': 63, 'E': 68, 'T': 'int','default' : 1},
                # 'line2_checksum' : {'B':68,'E':69,'T':'int'}

                }
    }
##########################################################################
# NOTE: pre-pend all non-data fields with "_".  Any class field that does NOT have this
# as a starting char will be automatically exported as data to the DB
class tle_class:
    '''
    A class to parse and store two-line element set data from the NORAD files
    This currently does NO error checking (even checksum).  Next version

    All field units are native to TLE's (degrees, minutes)
    '''

    def _set_tle_fields(self):
        for key in tle_fields['line1'] : setattr( self, key, None )
        for key in tle_fields['line2'] : setattr( self, key, None )

    def _set_defaults( self ):
        for key in tle_fields['line1'] : setattr( self, key, tle_fields['line1'][key]['default'])
        for key in tle_fields['line2'] : setattr( self, key, tle_fields['line2'][key]['default'])

    def _clear(self):
        self.error = 0
        self._set_tle_fields()
        self.catalog_file = None
        self._astro_date = None

    def _parse_line(self, LINE, DATA):
        tLine = DATA.strip()
        setattr(self, LINE, DATA.strip())
        for Z in tle_fields[LINE].keys():
            field_spec = tle_fields[LINE][Z]
            begin,end,dtype = field_spec['B'], field_spec['E'], field_spec['T']
            try:
                newVal = get_tle_datatype(DATA[begin:end],dtype)
            except:
                self.error = 1
                self._err_msg("exception on line: {} field: {}".format(LINE,Z) )
                return
            if newVal == None:
                 self._err_msg('get_tle_datatype failed on field {}, data: {}, type: {}'.format( Z, begin, end, dtype ) )
                 self.error = 1
                 return
            setattr(self, Z, newVal)

    def _derive_dates(self):
        try:
            self.epoch_int_day = int(self.epoch_day)
            self.epoch_frac_day = self.epoch_day - self.epoch_int_day
        except:
            self._err_msg( 'unexpected error in _derive_dates')
            formatExceptionInfo()
            self.error = 1
        if self.epoch_year <= 57: self.epoch_full_year = self.epoch_year + 2000
        else:   self.epoch_full_year = self.epoch_year + 1900
        self._astro_date = astro_date(tle_epoch.year_day_to_datetime( self.epoch_year, self.epoch_day) )

    def set_date( self , DT ):
        try:
            self.epoch_year = int( DT.strftime('%y'))
            self.epoch_day  = int( DT.strftime('%j'))
            TT = DT.utctimetuple()  
            FRAC = ((1.0/24.0) * TT.tm_hour) + ((1.0/1440.0) * TT.tm_min) + ((1.0/86400.0) * TT.tm_sec) + ((1.0/86400.0e6) * DT.microsecond)
            self.epoch_day += FRAC
        except: 
            self._err_msg('set_date: could not handle {}'.format(DT))
            return

    def printFields( self ):
        A = str()
        if self.catalog_file != None: A = "tle_class (%s)" % self.catalog_file
        else: A = ''
        for I in self.__dict__.keys():
            A = A + str(I) + " ---> " + str(self.__dict__[I]) + "\n"
        return A

    def __str__( self ):
        return "{}\n{}".format( self._generate_new_line1(), self._generate_new_line2() )


    def _struct_print( self ):
        print (self.sat_no, self.classification, self.launch_year, self.launch_number, self.launch_place, self.epoch_year, self.epoch_day, self.mean_motion_1, self.mean_motion_2, self.bstar, self.ephem_type, self.element_number, self.line1_checksum )

    def __lt__( self, other ):
        if not isinstance( other, tle_class ): return NotImplemented
        return self._astro_date.to_jday() < other._astro_date.to_jday() 

    def __gt__( self, other ):
        if not isinstance( other, tle_class ): return NotImplemented
        return self._astro_date.to_jday() > other._astro_date.to_jday() 

    def __ge__( self, other ):
        if not isinstance( other, tle_class ): return NotImplemented
        return self._astro_date.to_jday() >= other._astro_date.to_jday() 

    def __le__( self, other ):
        if not isinstance( other, tle_class ): return NotImplemented
        return self._astro_date.to_jday() <= other._astro_date.to_jday() 

    def __eq__( self, other ):
        if not isinstance( other, tle_class ): return NotImplemented
        return self._astro_date.to_jday() == other._astro_date.to_jday() 

    def _mongoDB( self ):
        retVal = dict()
        for F in dir(self): 
            if not F[0]=="_" and not callable( F ): retVal[F] = getattr(self,F)
        return retVal

    def _check_lines( self ):
        while len(self.line1) < 69 : self.line1+='0'
        while len(self.line2) < 69 : self.line2+='0'

    def __init__(self, line1=None, line2=None, catalogName='N/A'):
        self._clear()
        if line1 == None or line2 == None: 
            self._set_defaults()
            self.catalog_file = catalogName
            return
        self.catalog_file = catalogName
        self.from_lines( line1, line2 )

    def _generate_new_tle( self ):
        ''' generate_new_tle(): we can modify the TLE values stored in this class, and generate a representative TLE.
                    Note, that the original TLE lines are still stored in line1, line2
                    --> NOTE: we are currently kluging checksums by just appending a zero'''
        try: return self._generate_new_line1() + '\n' + self._generate_new_line2()
        except: return None

    def __str__( self ):
        return self._generate_new_tle( )

    def getLines( self ):
        try: return self._generate_new_line1(), self._generate_new_line2()
        except: return None

    def _generate_new_line1( self ):
        rV = '1 {}U {:02d}{:03d}{:3s} {:02d}{:012.8f} '.format(integer5(self.sat_no),self.launch_year,self.launch_number,self.launch_piece,self.epoch_year,self.epoch_day)

        mm1 = [x for x in "{:+14.13f}".format( self.mean_motion_1 )]
        rV += str(mm1[0]) + ''.join(mm1[ mm1.index('.') : mm1.index('.')+9 ] )

        rV += ' %s' % generate_expo_format( self.mean_motion_2 )
        rV += ' %s' % generate_expo_format( self.bstar )

        rV += ' 0 %04d' % self.element_number
        rV += '0' # fake checksum
        return rV

    def _generate_new_line2( self ):
        rV = '2 {} {:08.4f} {:08.4f} '.format( integer5(self.sat_no), self.inclination, self.ra) 
        rV += '{:9.7f}'.format( self.eccentricity ).partition('.')[2]
        rV += ' {:08.4f} {:08.4f} {:011.8f}{:05d}'.format( self.argument, self.mean_anomaly, self.mean_motion, self.rev_number )
        rV += '0' # fake checksum
        return rV

    def _err_msg( self, S ):
        sys.stderr.write('%s %s\n' % ('tle_class', S) )

    def _calculate_apogee_perigee( self, earth_rad = 6378.135):
        ''' Default value for earth_rad is taken from space-track.
        space-track : https://www.space-track.org/documentation#/faq
        Additional references: http://www.satobs.org/seesat/Dec-2002/0197.html
        '''
        semi_major = (8681663.653 / self.mean_motion) ** (2.0/3.0)
        self.perigee = ( semi_major * (1 - self.eccentricity) ) - earth_rad
        self.apogee =  ( semi_major * (1 + self.eccentricity) ) - earth_rad

    def from_lines( self, line1, line2 ):
        self._parse_line( 'line1', line1 )
        self._parse_line( 'line2', line2 )
        self._derive_dates()
        self._check_lines()
        self._calculate_apogee_perigee()

    def from_sgp4_satrec( self, satrec ):
        ''' 
        Brandon Rhode's SGP4 Python code is a direct port from Vallado's code, and outputs a 'satrec' structure after
        parsing a TLE.  Let's see if we can use that to serialize as well

        Note that some fields are converted immediately by SGP4 on ingest, so we'll convert them back (this serializer
        does *not* keep units in the code, just raw fields

        #  ---- find no, ndot, nddot ----
        satrec.no_kozai = satrec.no_kozai / xpdotp; #   rad/min
        satrec.nddot= satrec.nddot * pow(10.0, nexp);
        satrec.bstar= satrec.bstar * pow(10.0, ibexp);

        #  ---- convert to sgp4 units ----
        satrec.ndot = satrec.ndot  / (xpdotp*1440.0);  #   ? * minperday
        satrec.nddot= satrec.nddot / (xpdotp*1440.0*1440);

        #  ---- find standard orbital elements ----
        satrec.inclo = satrec.inclo  * deg2rad;
        satrec.nodeo = satrec.nodeo  * deg2rad;
        satrec.argpo = satrec.argpo  * deg2rad;
        satrec.mo    = satrec.mo     * deg2rad;
        '''
        self.sat_no = satrec.satnum
        try: self.classification = satrec.classification
        except: self.classification = 'U'
        try: self.launch_year = launch_year( satrec.intldesg) 
        except: self.launch_year = 99
        try: self.launch_number = launch_number( satrec.intldesg )
        except: self.launch_number = 1
        try: self.launch_piece = launch_piece( satrec.intldesg )
        except: self.launch_piece = 'A'
        self.epoch_full_year =  satrec.epochyr
        if self.epoch_full_year >= 2000: self.epoch_year = self.epoch_full_year - 2000
        if self.epoch_full_year < 2000: self.epoch_year = self.epoch_full_year - 1900
        self.epoch_day = satrec.epochdays
        self.mean_motion_1 = satrec.ndot * (xpdotp*1440.0)          # <--- converted
        self.mean_motion_2 = satrec.nddot *(xpdotp*1440.*1440.)     # <--- converted
        self.bstar = satrec.bstar
        try: self.ephem_type = int( satrec.ephtype )
        except: self.ephem_type = 0
        try: self.element_number = int( satrec.elnum )
        except: self.element_number = 1
        self.inclination = math.degrees( satrec.inclo )
        self.ra = math.degrees( satrec.nodeo )
        self.eccentricity = satrec.ecco
        self.argument = math.degrees( satrec.argpo )
        self.mean_anomaly = math.degrees( satrec.mo )
        self.mean_motion = satrec.mdot * xpdotp                 # <--- converted
        try: self.rev_number = int(satrec.revnum)
        except: self.rev_number = 0
        return self

def testTLE():
    from sgp4.io import twoline2rv 
    import sgp4.earth_gravity as EG

    L1 = '1 25544U 98067A   20125.70112759  .00000857  00000-0  23422-4 0  9999'
    L2 = '2 25544  51.6442 203.6971 0001138 234.0065 260.8716 15.49346924225270'

    # have SGP4 parse that TLE
    satrec = twoline2rv( L1, L2, EG.wgs72 )
    # build our local class from that line
    A = tle_class( L1, L2, catalogName='blah')
    # build a local class from the satrec
    B = tle_class()
    B.from_sgp4_satrec( satrec )

    # print the different TLEs
    print(' -- Original -- ')
    print(L1)
    print(L2)
    print()
    print(' -- tle_class -- ')
    print( str(A) )
    print()
    print(' -- from satrec -- ')
    print( str(B) )
    print()

    # make some modifications
    B.inclination = 90
    B.eccentricity = 0.75
    B.argument = 0.1
    print(' -- modified TLE (inc: 90, eccen: 0.75, arg: 0.1) --')
    print(' NOTE: this is not a valid TLE, just field mashing ' )
    print( str(B) )


# =====================================================================================================
# =====================================================================================================
if __name__=="__main__":
    testTLE()
