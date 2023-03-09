# KNW changes for alphanumeric TLE
# alpha numeric lookups
from_alpha = {'A': 10, 'B': 11, 'C': 12, 'D': 13, 'E': 14, 'F': 15, 'G': 16, 'H': 17, 'J': 18, 'K': 19, 'L': 20, 'M': 21, 'N': 22, 'P': 23, 'Q': 24, 'R': 25, 'S': 26, 'T': 27, 'U': 28, 'V': 29, 'W': 30, 'X': 31, 'Y': 32, 'Z': 33}

to_alpha = {10: 'A', 11: 'B', 12: 'C', 13: 'D', 14: 'E', 15: 'F', 16: 'G', 17: 'H', 18: 'J', 19: 'K', 20: 'L', 21: 'M', 22: 'N', 23: 'P', 24: 'Q', 25: 'R', 26: 'S', 27: 'T', 28: 'U', 29: 'V', 30: 'W', 31: 'X', 32: 'Y', 33: 'Z'}


def alpha_to_integer(s):
    # from Brandon Rhodes' SGP4 library
    ''' compute an INTEGER from a TLE number string'''
    if isinstance( s, int) : return int(s)
    if not s[0].isalpha(): return int(s)
    c = s[0]
    return (from_alpha[c] * 10000 ) + int(s[1:])

def integer_to_alpha(I):
    ''' compute a string from an integer '''
    I = int(I)
    if I > 339999 : 
        raise Exception('cannot convert integers > 339999 to TLE strings')
    if I < 100000 : return str(I).zfill(5)
    intstr = str(I)
    lkup = intstr[0:2]
    return to_alpha[ int(lkup) ] + intstr[2:][0:5]


