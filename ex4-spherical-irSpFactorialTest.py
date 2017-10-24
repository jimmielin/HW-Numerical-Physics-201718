def _irSpFactorial(l):
    # (2l-3)!!(2l-3)!!/(2l-1)! = ((2l-3)(2l-5)...)^2/(2l-1)(2l-2)(2l-3)...
    # u*u/(u+2)
    result = 1
    ptr2 = 2*l - 3
    ptr1 = 2*l - 1
    while(ptr2 != -1):
        result = result * ptr2 / ptr1 
        if(ptr1 % 2 == 0):
            ptr2 = ptr2 - 2
        ptr1 = ptr1 - 1
    return result

# 5!!5!!/7!
print(_irSpFactorial(100000))