import math


#Function to calculate the TM based on the nearest neighbors method according sigmaaldrich.com
def tm_calc(seq):

    deltaH = 0
    deltaS = 0

    for i in range(len(seq) - 1):
        pair = seq[i] + seq[i + 1]
        if 'AA' in pair or 'TT' in pair:
            deltaH += -9.1
            deltaS += -0.0240
        elif 'AT' in pair:
            deltaH += -8.6
            deltaS += -0.0239
        elif 'TA' in pair:
            deltaH += -6.0
            deltaS += -0.0169
        elif 'CA' in pair or 'TG' in pair:
            deltaH += -5.8
            deltaS += -0.0129
        elif 'GT' in pair or 'AC' in pair:
            deltaH += -6.5
            deltaS += -0.0173
        elif 'CT' in pair or 'AG' in pair:
            deltaH += -7.8
            deltaS += -0.0208
        elif 'GA' in pair or 'TC' in pair:
            deltaH += -5.6
            deltaS += -0.0135
        elif 'CG' in pair:
            deltaH += -11.9
            deltaS += -0.0278
        elif 'GC' in pair:
            deltaH += -11.1
            deltaS += -0.0267
        elif 'GG' in pair or 'CC' in pair:
            deltaH += -11.0
            deltaS += -0.0266
    A = float(-0.0108)
    R = float(0.00199)
    C = float(0.0000005)
    T = float(-273.15)
    CNa = float(0.05)
    term1 = deltaH / (A + deltaS + R * math.log(C / 4))
    Tm_i = term1 +T
    Tm= Tm_i +16.6 * math.log10(CNa)
    return Tm
seq=str(input("Your sequence: ))
result=tm_calc(seq)
print(result)
