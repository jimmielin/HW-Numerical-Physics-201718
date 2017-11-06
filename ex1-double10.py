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
# Now Playing: Into Your Arms - The Maine
####################################################################

####################################################################
# Decimal to IEEE Any-Precision Floating Point Format
####################################################################
# string dec2ieee(string strN)
# Accepts n in "string" format, returning the IEEE any-precision result
# in a "string" format. (By default as an example we give double)
# The function is freely configurable.
def dec2ieee(strN):
    # ~ * Configuration Parameters Here * ~
    # dec2ieee supports any (theoretically) IEEE-specified precision spec,
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
    # TODO: Sanity check if strN is actually a number
    if(strN[0] == "-"):
        sign = 1
        strN = strN[1:]
    elif(strN[0] == "+"):
        sign = 0
        strN = strN[1:]
    else:
        sign = 0

    components = strN.split('.')
    if(len(components) > 2):
        raise ValueError("The decimal inputted has more than a whole and fractional component.")

    # "Whole" component
    # accept e.g. ".572" input.
    whole = (int(components[0]) if len(components[0]) > 0 else 0) 

    # "Fractional" component
    # accept e.g "17." "17" input.
    frac  = (components[1] if len(components) == 2 and len(components[1]) > 0 else "0") 

    # For the "Whole" component, decompose it as a set of exponents
    # ranging from 2^1024 all the way down to 2^0 = 1
    # For 0, it is a special case and we'll handle it differently.
    # Store the "whole part" (pre-dot) binary result in the wholeResult string.
    wholeResult = ""
    if(whole == 0):
        wholeResult = "0"
    else:
        # Gradually write the result to wholeResult, from 2^1023 to 0
        i = expU
        while(i != -1):
            # print("* dec2double diag: [whole] whole = ", whole, " checking i =", i, flush=True)
            if(whole >= 2**i):
                whole -= 2**i
                wholeResult = wholeResult + "1"
            else:
                wholeResult = wholeResult + "0"
            i -= 1

    # Strip the trailing left-zeroes
    wholeResult = wholeResult.lstrip("0")
    if(wholeResult == ""):
        wholeResult = "0"

    # For the "Fractional" component, it is trickier as you have to
    # check for the *negative* exponents of 2, which may of course
    # underflow at any time, also at the cost of precision loss.
    #
    # In order not to be constrained by the Python's default "double"
    # precision (which results in increasingly lower precision as one
    # goes below and below), we had to devise a runtime compiled set of
    # "Integer" coded negative powers of 2 with a bunch of zeroes,
    # and compare the result with that instead,
    # after padding both the "frac" part and all these results.
    #
    # A reference, Python-double dependent implementation is commented out below
    # for your viewing pleasure.

    # fracResult = ""
    # doubleFrac = float("." + frac.__str__()) # I feel so ashamed for using python double
    # i = -1
    # while(i != -1075):
    #     if(doubleFrac >= 2**i):
    #         doubleFrac -= 2**i
    #         fracResult = fracResult + "1"
    #     else:
    #         fracResult = fracResult + "0"
    #     i -= 1
    # fracResult = fracResult.rstrip("0")

    # 2**(-1074) is approximately 4.6e-625, 2**(-1073) is 9.22e-625
    # We need to ensure precision of 2**(-1074). A naive try would be 1e626,
    # but we need much more than that to prevent precision loss beyond the point
    # of 2**-623. A experiment reveals we need to use 1e1074 as the starting seed,
    # so we can get all the 751 significant digits of 2**(-1074)

    # Starting from 2^(-1)...
    # ## PRECOMPILATION OF INTEGER-CODED NEGATIVE TWO POWERS FOR IEEE ##
    # For double, configure as follows: decPowerOfTenSeed = 1074 (1022+52),
    # minNegTwoPtr = 1075. minNegTwoPtr can be arbitrarily large, but don't
    # be wasteful.
    decPowerOfTenSeed = expL + P - 1
    minNegTwoPtr = expL + P + 1

    # Build this table (does not take long)
    negTwoIntResultsTable = [None]
    negTwoPtr  = 2
    negTwoCurr = int("5" + "0" * decPowerOfTenSeed)
    negTwoIntResultsTable.append(negTwoCurr)
    while(negTwoPtr != minNegTwoPtr):
        negTwoCurr = negTwoCurr // 2
        negTwoIntResultsTable.append(negTwoCurr)
        # print(negTwoPtr, ":", negTwoCurr) # diag only
        negTwoPtr += 1

    # Now we can do the calculations necessary
    fracResult = ""

    # First, padd the number to our desired Integer-coded precision
    # Given 0.5 = 5 + 0 * decPowerOfTenSeed (strlen = decPowerOfTenSeed + 1)
    frac = frac.ljust(1 + decPowerOfTenSeed, '0')
    intFrac = int(frac) # mutable
    i = -1
    while(i != (-1) * (minNegTwoPtr)):
        # print("* dec2double diag: i =", i, " intFracRemain =", intFrac, flush=True)
        if(intFrac >= negTwoIntResultsTable[(-1) * i]):
            fracResult = fracResult + "1"
            intFrac -= negTwoIntResultsTable[(-1) * i]
        else:
            fracResult = fracResult + "0"
        i -= 1

    # Strip the trailing right-zeroes
    fracResult = fracResult.rstrip("0")
    # print("* dec2double diag: out of fracResult=", fracResult, flush=True)

    # Align the dots (not stars) properly
    if(wholeResult == "0"):
        offsetTwoPower = (-1) * (len(fracResult) - len(fracResult.lstrip("0")) + 1)
        irBinaryRepTail = fracResult.lstrip("0")[1:]
    else:
        offsetTwoPower = len(wholeResult) - 1
        # irBinaryRep = wholeResult[:1] + "." + wholeResult[1:] + fracResult
        irBinaryRepTail = wholeResult[1:] + fracResult
    # print("* dec2double diag: intermediate result " + irBinaryRepTail + " x 2^", offsetTwoPower)

    # Calculate the exponential offset (+1023 for IEEE Double Float)
    # 2^(e-1) - 1, e is the length of the exponential bit storage size
    # e = 11 for double-precision
    decExpOffset = 2**(E - 1) - 1 + offsetTwoPower
    # print("* dec2double diag: decExpOffsetted=", decExpOffset)

    # Convert to binary (copied from above)
    binExpOffset = ""
    i = E - 1
    while(i != -1):
        # print("* dec2double diag: [whole] whole = ", whole, " checking i =", i, flush=True)
        if(decExpOffset >= 2**i):
            decExpOffset -= 2**i
            binExpOffset = binExpOffset + "1"
        else:
            binExpOffset = binExpOffset + "0"
        i -= 1

    # print("* dec2double diag: expOffset (bin) =", binExpOffset)

    # Generate the final IEEE double-precision representation,
    # chopping off any non-significant bits (sorry precision!)
    result = sign.__str__()
    result = result + binExpOffset
    # print("* dec2ieee diag: ieee tail raw result =", irBinaryRepTail)
    if(len(irBinaryRepTail) < P - 1):
        result = result + irBinaryRepTail.ljust(P - 1, '0')
    else:
        # Handle rounding errors... round up if necessary
        if(len(irBinaryRepTail) >= P and irBinaryRepTail[P-1] == "1"):
            lookAt = P-2
            while(irBinaryRepTail[lookAt] == "1"):
                # Also note that we need to decrement the 1s on this way...
                irBinaryRepTail = irBinaryRepTail[0:lookAt] + "0" + irBinaryRepTail[(lookAt+1):]
                lookAt -= 1
            irBinaryRepTail = irBinaryRepTail[0:lookAt] + "1" + irBinaryRepTail[(lookAt+1):]

        result = result + irBinaryRepTail[0:(P-1)]

    return result
    
print(dec2ieee(input("> dec: ")))
