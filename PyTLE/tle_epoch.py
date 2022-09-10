from datetime import datetime,timedelta
from math import floor

def tle_epoch_to_datetime( EP ):
	'''tle_epoch_to_datetime( TLE epoch field ):
		-> this prefers a string, but will try to convert.  
		-> note that if you pass in a float, you'll get some floating point round-off (tle epoch is a big number)'''
	if not isinstance(EP,str): 
		try: EP = str(EP)
		except: return None
	# fix any JSpOC blanks
	EP = EP.replace(' ','0')
	try: 
		(YD,crap,frac) = EP.partition('.')
		DAY = datetime.strptime( YD, "%y%j" )
		FRAC = float( '0.%s' % frac )
		return DAY + timedelta( seconds = ((86400.0) *  FRAC ) )
	except: return None

def year_day_to_datetime( YEAR, DAY ):
	'''year_day_to_datetime( year, day ):
		take TLE formatted date YYDDD.DDDDDD and return a datetime object
		year : int
		day  : float

            # ------------- from stackoverflow ---- https://stackoverflow.com/questions/34849083/convert-tle-times-decimal-days-to-seconds-after-epoch
            # get year 2 digit and floating seconds days 
            y_d, nbs = "16012.375".split('.') 

            # parse to datetime (since midnight and add the seconds) %j Day of the year as a zero-padded decimal number.
            d = datetime.datetime.strptime(y_d, "%y%j") + datetime.timedelta(seconds=float("." + nbs) * 24 * 60 * 60)
            # 1.0 => 1 day
            # from time tuple get epoch time. 
            time.mktime(d.timetuple())r
                '''
	if 57 <= YEAR <= 99: YEAR += 1900
	if 0 <= YEAR <= 56: YEAR += 2000

	S = '{:02d}{:03d}'.format( YEAR, int(DAY) )	
	FRAC = DAY - floor(DAY)
	return datetime.strptime( S, '%Y%j' ) +   timedelta( seconds=((86400.0) *  FRAC ) )

def datetime_to_tle_epoch( DT ):
	TT = DT.utctimetuple()
	FRAC = ((1.0/24.0) * TT.tm_hour) + ((1.0/1440.0) * TT.tm_min) + ((1.0/86400.0) * TT.tm_sec) + ((1.0/86400.0e6)*DT.microsecond)
	RV = int( DT.strftime('%y%j') )
	return RV+FRAC
	#except: return None

if __name__ == "__main__":
	IN = '70090.03500497'
	print (IN)
	print (tle_epoch_to_datetime( IN ))
	print (year_day_to_datetime( 70, 90.03500497 ))
	print (year_day_to_datetime( 1970, 90.03500497 ))
	A = tle_epoch_to_datetime( float(IN) )
	print (A)
	print (datetime_to_tle_epoch( A ) )
