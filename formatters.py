# ###############################################################################
# MIT License
#
# Copyright (c) 2023 Kerry Wood
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# ###############################################################################

from datetime import datetime, timedelta
import numpy as np

# -----------------------------------------------------------------------------------------------------
def generate_checksum(line):
    digits = sum(map(lambda x: int(x), filter(
        lambda y: re.match('[0-9]{1}', y), [z for z in line])))
    minus = line.count('-')
    rV = str(digits + minus)
    return rV[-1]


# -----------------------------------------------------------------------------------------------------
def generate_expo_format(flt):
    [mant, crap, exp] = '{:+4.4e}'.format(flt).partition('e')
    mant = mant.replace('.', '')
    if abs(int(exp)) > 9 : raise Exception('exponent is too large to express')
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
    days = int(frac)
    decimals = frac - days 
    decimals = '{}'.format(np.round( decimals, 8 )).split('.')[1]
    year = dt.strftime('%y')[:2]
    days = '{:03d}'.format(days)[:3]
    return '{}{}.{}'.format( year, days, decimals )
