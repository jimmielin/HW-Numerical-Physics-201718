####################################################################
# Computational Physics, 2017-18 Sem1
# HW-1 Ex-1
#
# (c) 2017 Haipeng Lin <linhaipeng@pku.edu.cn>
# All Rights Reserved.
#
# This program is written as Homework for the Computational Physics
# Course in Peking University, School of Physics.
# NO WARRANTIES, EXPRESS OR IMPLIED, ARE OFFERED FOR THIS PRODUCT.
# Copying is strictly prohibited and its usage is solely restricted
# to uses permitted by the author and homework grading.
#
# This program is Python 3 (3.6.2 on MSYS_NT)
# compatible, and does not support Python 2.
#
# Now Playing: 山丘 - 李宗盛
####################################################################

####################################################################
# IEEE Any-Precision Floating Point Format to Decimal
####################################################################
# string ieee2dec(string strI)
# Accepts i in "string" format, returning the decimal result
# in a "string" format. (By default as an example we give double)
# The function is freely configurable.
def ieee2dec(strN):
    # ~ * Configuration Parameters Here * ~
    # ieee2dec supports any (theoretically) IEEE-specified precision spec,
    # simply configure the parameters below to your liking.
    #
    # Beta = 2 (fixed for binary, this you can't change, sorry)
    # Precision: Number of precision bits (can be ARBITRARY), d0 ~ d(p-1)
    P = 53
    # L (lower exponential bound), E >= -L
    expL = 1022
    # U (upper exponential bound), E <= +U
    expU = 1023
    # Exponential Bits
    E = 11
    # The one sign bit is contained in P.
    # All numbers generated are SIGNED (for UNSIGNED it is trivial to modify the below code.)

    # ~ * NO USER CONFIGURABLE VARIABLES BELOW * ~
    result = ""

    # Split the strings according to configuration.
    signBit = strN[0]
    if(signBit == "1"):
        # is negative
        result = result + "-"

    expBits = strN[1:(E+1)]
    valBits = strN[(E+1):]

    # Calculate the exponent by deducting the exponent bias
    # (IEEE defined as 2**(E-1) - 1)
    expVal = 1 - 2**(E-1)
    i = E - 1
    while(i != -1):
        expVal += 2**(E-1-i) * (1 if expBits[i] == "1" else 0)
        i -= 1

    # Isolate ValBits through exponent.
    # If exponent is non-negative, then simply
    # 1.xxxxxxxx -> move . to the right expVal bits
    wholeBits = ""
    fractionBits = ""

    if(expVal >= 0):
        wholeBits = "1" + valBits[:(expVal)]
        if(len(wholeBits) < (1 + expVal)):
            wholeBits = wholeBits.ljust(1 + expVal, "0")
        fractionBits = valBits[(expVal):]
        if(len(fractionBits) == 0):
            fractionBits = "0"
    else:
        wholeBits = "0"
        # If its, e.g. raw result
        # 0.0001100101... -> 1.100101... x 2^(-4) (base 2)
        # then the fractional result is 0001100101..
        # Using 4-1=3 padded left-zeroes, then a 1, then the given component
        # (100101...) encoded in the IEEE standard
        fractionBits = "0" * ((-1) * expVal - 1) + "1" + valBits

        # print("* ieee2dec diag: expVal=", expVal, " is neg, wholeBits=0, fractionBits=", fractionBits)
    
    # print("* ieee2dec diag: expVal conversion", expBits, "->", expVal)

    # For the whole number component it is easy to recreate
    # by simply multiplying by 2 where it counts
    # Up from expU down to 0
    # We refrain from left-padding zeroes; instead we simply loop on the string
    # Reverse the bits
    wholeResult = 0
    # print("* ieee2dec diag wholeBits:", wholeBits)
    wholeBitsR = wholeBits[::-1]
    for i in range(len(wholeBitsR)):
        wholeResult += 2**i * (1 if wholeBitsR[i] == "1" else 0)

    fracResult = 0
    # For higher-precision, non-Python double restrained results
    # we need to do as dec2ieee and pre-build a decimal table,
    # then later adjust precision as needed.
    # ## PRECOMPILATION OF INTEGER-CODED NEGATIVE TWO POWERS FOR IEEE ##
    # For double, configure as follows: decPowerOfTenSeed = 1074 (1022+52),
    # minNegTwoPtr = 1075. minNegTwoPtr can be arbitrarily large, but don't
    # be wasteful.
    decPowerOfTenSeed = expL + P - 1
    minNegTwoPtr = expL + P + 1

    # Build this table (does not take long)
    negTwoIntResultsTable = []
    negTwoPtr  = 2
    negTwoCurr = int("5" + "0" * decPowerOfTenSeed)
    negTwoIntResultsTable.append(negTwoCurr)
    while(negTwoPtr != minNegTwoPtr):
        negTwoCurr = negTwoCurr // 2
        negTwoIntResultsTable.append(negTwoCurr)
        negTwoPtr += 1

    # Now simply generate the number. No padding is required as the
    # fracResult is already post-decimal divider component.
    for i in range(len(fractionBits)):
        # note this begins with 0, 1, 2, ...
        # but the actual index meaning should be -1, -2, -3, ...
        # luckily, negTwoIntResultsTable is also 0, 1, 2, ... indexed
        # print("* ieee2dec diag: adding", negTwoIntResultsTable[i])
        fracResult += negTwoIntResultsTable[i] * (1 if fractionBits[i] == "1" else 0)

    # Padd fracResult to desired precision
    fracResult = fracResult.__str__()
    if(len(fracResult) < (1 + decPowerOfTenSeed)):
        fracResult = fracResult.rjust(1 + decPowerOfTenSeed, "0")

    # Chop chop
    fracResult = fracResult.rstrip("0")

    result = result + wholeResult.__str__() + "." + fracResult.__str__()

    return result

print(ieee2dec(input("> ieee: ")))