import math

class BoundsError(Exception):
    pass

def brent(a, b, tol, fun, *args):
    '''Uses Brent's Method to find the root; an improved secant search that also falls
back to bisection if the secant search fails.
a, b are bounds (which is low and which is high doesn't matter)
tol is the tolerance
fun is the function to find the root of
*args are other arguments required by that function.'''

    fa = fun(a,*args)
    fb = fun(b,*args)
    if fa*fb > 0:
        raise BoundsError('a and b do not bracket the root of fun(x)')

    if abs(fa) < abs(fb): #swap a and b
        a, b = b, a
        fa, fb = fb, fa

    c, fc = a, fa
    d = c
    modify = True
    err = abs(b - a)
    while err > tol:
        

        if fc not in [fa, fb]: #parabolic interpolation
            s = (a*fb*fc)/((fa-fb)*(fa-fc)) + (b*fa*fc)/((fb - fa)*(fb-fc)) + (c*fa*fb)/((fc-fa)*(fc-fb))

        else: #linear interpolation, aka secant method
            s = b - fb*(b-a)/(fb-fa)

        #we define some booleans to make the next if statement not so monstrous
        c1 = s < (3*a + b)/4 or s > b
        c2 = modify and abs(s-b) >= abs(b-c)/2
        c3 = not modify and abs(s-b) >= abs(c-d)/2
        c4 = modify and abs(b-c) < 1e-8
        c5 = not modify and abs(c-d) < 1e-8

        if c1 or c2 or c3 or c4 or c5:
            s = (a + b)/2 #bisection
            modify = True
        else:
            modify = False

        fs = fun(s, *args)
        d = c
        c = b

        if fa*fs < 0:
            b, fb = s, fs
        else:
            a, fa = s, fs

        if abs(fa) < abs(fb):
            a, b = b, a
            fa, fb = fb, fa

        err = abs(b-a)

    return s

def bisect(high, low, tol, fun, *args):
    '''Given a function handle fun that takes one variable as an argument,
a value high such that fun(high) > 0, a value low such that fun(low) < 0,
and a tolerance tol, find the root of the function within tol.
    *args is a tuple containing the other, nonvariable necessary arguments to the function.'''

    x1 = high
    x2 = low
    err = 1
    while err > tol:
        x = 0.5*(x1 + x2)
        f = fun(x, *args)

        #now actually re-assign the bisection
        if f > 0:
            #err = abs(x1 - x)
            x1 = x
        elif f < 0:
            #err = abs(x2 - x)
            x2 = x
        else: #f1 = 0; we got very lucky
            return x

        err = abs(f)
    return x


