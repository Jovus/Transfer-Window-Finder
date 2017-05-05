import math, vector, roots
from scipy.constants import G, pi
from scipy.optimize import fsolve #just in case

M0 = 1.989e30 #mass of the Sun
MU = G*M0 #units: m**3/s**2

#find machine epsilon for 1.0

MEPSILON = 1
E0 = MEPSILON/10
while 1.0-E0 !=1.0:
    MEPSILON = E0
    E0 = MEPSILON/10

#print(MEPSILON)

days_convert = 86400 # num of seconds in a (Julian) day
AU_convert = 149.597870e9 # num of meters in AU

#testing data
R1 = (0.473265, -0.899215, 0) #in AUs...don't forget to convert
R2 = (0.066842, 1.561256, 0.030948)
DT = 207 #days
 
def _phi(x, fun='acot'):
    '''The x or y part of F(x) = phi(x) - phi(y) - tau; or the root equation for Sun's method.
fun can be either 'acot' or 'atan' and determines which one we use (acot is for x, atan is for y)'''

    #these two should actually give exactly the same answer, but we include the decision tree to be safe
    g = math.sqrt(1 - x*x)
    if fun == 'acot': #phix
        out = acot(x/g) - (2 + x*x) * g / (3*x)
    elif fun == 'atan': #phiy
        out = math.atan(g/x) - (2 + x*x) * g / (3*x)
    else:
        raise InvalidArgumentError('Fun must be either acot or atan')

    return out
    
def acot(x):
    '''Given some x, return the angle corresponding to the inverse cotangent of x'''

    return (pi/2) - math.atan(x)

def acoth(x):
    '''Given some x, return the angle corresponding to the inverse hyperbolic cotangent of x'''
        
    return 0.5*math.log((x+1)/(x-1)) #clearly, undefined at x=1

        
def _getY(x, sigma):
    '''Given x and the angle parameter, return y for that x'''

    y = math.sqrt(1 - sigma**2 * (1 - x*x))
    if sigma < 0:
        y = -y

    return y

def _Ftau(x, angleParameter, parabolicNormTime, normalTime, N):
    '''Equation we want to find the root of in order to tell us what x is for this problem.'''
    y = _getY(x, angleParameter)

    if abs(1-x)<1e-8: #close enough to parabolic
        tofdiff = parabolicNormTime - normalTime #this had better be zero

    elif x > 1: #hyperbolic trajectory
       # print('x = {}'.format(x))
        g = math.sqrt(x*x - 1)
        h = math.sqrt(y*y - 1)
        tofdiff = (-acoth(x/g) + acoth(y/h) + x*g - y*h) / (g**3) - normalTime

    else: #elliptical case; -1 < x < 1
        g = math.sqrt(1 - x*x)
        h = math.sqrt(1 - y*y)
        tofdiff = (acot(x/g) - math.atan(h/y) - x*g + y*h + N*pi) / (g**3) - normalTime

    return tofdiff

def _solutions(x, y, n, R1, R2, angle):
    '''Gives solutions to the Lambert problem, in the form of (R1, R2, V1, V2, transferAngle, #revs)'''
    #because Python suffers strict scope, we get to recalculate known quantities
    #but this saves clarity in reading over just passing all these as args
    r1 = vector.mag(R1)
    r2 = vector.mag(R2)
    deltaR = vector.sub(R2, R1)
    c = vector.mag(deltaR)
    m = r1 + r2 + c
    n = r1 + r2 - c
    sqrtMu = math.sqrt(MU)
    invSqrtM = 1/math.sqrt(m)
    invSqrtN = 1/math.sqrt(n)

    vc = sqrtMu * (y*invSqrtN + x*invSqrtM) #scalar
    vr = sqrtMu * (y*invSqrtN - x*invSqrtM) #scalar
    ec = [elem*(vc/c) for elem in deltaR]
    V1 = tuple(vector.add(ec, [elem*(vr/r1) for elem in R1]))
    V2 = tuple(vector.sub(ec, [elem*(vr/r2) for elem in R2]))

    return [tuple(R1), tuple(R2), V1, V2, n*2*pi + angle]
    

def sunMethod(mu, R1, R2, dt, maxRevs=0): 
    '''Given the standard gravitational parameter mu, the starting position R1,
ending position R2, and time of flight dt, return parameter and semi-major axis of a keplerian orbit.
    R1, R2 are 3-element vectors for position, stored as tuples in the form (x,y,z)
    units are AU
    dt is in days
    maxRevs tells us how many total revolutions about the Sun are allowed

    We always want a prograde orbit

Implementation of Sun, F.T. "On the Minium Time Trajectory and Multiple Solutions of Lambert's Problem"
AAS/AIAA Astrodynamics Conference, Provincetown, Massachusetts, AAS 79-164, June 25-27, 1979'''

    #first, let's change to meters and seconds so our units are consistent

    R1 = [i*AU_convert for i in R1]
    R2 = [i*AU_convert for i in R2]

    dt = days_convert*dt

    #vector magnitudes
    r1 = vector.mag(R1)
    r2 = vector.mag(R2)

    deltaR = vector.sub(R2, R1) #change in position

    c = vector.mag(deltaR) #chord length
    m = r1 + r2 + c #augmented chord length
    n = r1 + r2 - c #diminished chord length

    solutions = []
    
    transferAngle = math.acos(vector.dot(R1, R2)/(r1*r2))
    if vector.cross(R1,R2)[2] < 0:
        transferAngle = 2*pi - transferAngle #we always want to go counterclockwise around +z

    #now we get to some parameters defined by Sun
    angleParameter = math.sqrt(n/m) #what the paper calls 'sigma'. Bounds -1 <= sigma <= 1
    if transferAngle > pi:
        angleParameter = -angleParameter

    normTime = 4*dt*math.sqrt(mu/(m**3)) #normalized dt; what the paper calls 'tau'
    parabolicNormTime = (2/3)*(1-angleParameter**3) #the normalized time taken for a parabolic trajectory from R1 to R2

    unwrapparabolicNormTime = (parabolicNormTime/4)/(math.sqrt(mu/(m**3))) #just a check variable

    parabolicDays = unwrapparabolicNormTime/days_convert # days taken by a parabolic trajectory

   #print('Parabolic time: {} days'.format(parabolicDays))

    #now we make an intelligent division of the solution space

    if abs(normTime - parabolicNormTime) < 1e-8: #parabolic solution
        x = 1
        if angleParameter < 0:
            y = -1
        else:
            y = 1

        solutions.append(_solutions(x, y, 0, R1, R2, transferAngle))

    elif normTime < parabolicNormTime: #hyperbolic solution
        x1 = 1
        x2 = 2

        #now we do a quick explosion to find the bounds of the hyperbolic solution
        ans = _Ftau(x2, angleParameter, parabolicNormTime, normTime, maxRevs)

        while ans > 0:
            print('ftau(x) = {}'.format(ans))
            x1 = x2
            x2 = x2*2
            ans = _Ftau(x2, angleParameter, parabolicNormTime, normTime, maxRevs)
        
        #and now that we've found bounds that straddle the zero, do a bisection search
        x = roots.brent(x1, x2, 1e-8, _Ftau, angleParameter, parabolicNormTime, normTime, maxRevs) #slow, but it works, and hyperbolic solutions should be rare
        y = _getY(x, angleParameter)
        solutions.append(_solutions(x, y, 0, R1, R2, transferAngle))
        
    else: #one of multiple possible elliptic solutions. Most of the rest of the algorithm rests here, in this else statement
        maxRevs = min(maxRevs, math.floor(normTime/pi)) #figure out how many revolutions we're allowed to go through; the larger ellipses will intersect on the return leg
        minENormTime = math.acos(angleParameter) + angleParameter*math.sqrt(1 - angleParameter**2) #the time of the minimum energy ellipse
        
        for n in range(maxRevs+1):
            if n > 0 and n == maxRevs:

                if angleParameter == 1:
                    x = 0
                    minNormTime = minENormTime
                elif angleParameter == 0:
                    augphi = lambda x: _phi(x) + n*pi
                    try:
                        x = roots.brent(0, 1-MEPSILON, 1e-6, augphi)
                    except:
                        return [0, 0, 0, 0, 0]
                    minNormTime = 2/(3*x)
                else:
                    diffphi = lambda x: _phi(x) - _phi(_getY(x, angleParameter), fun='atan') + n*pi

                    try:
                        x = roots.brent(0, 1-MEPSILON, 1e-6, diffphi)
                    except:
                        return [0, 0, 0, 0, 0]
                    minNormTime = 2/3 * (1/x - angleParameter**3 / abs(_getY(x, angleParameter)))

                #now we know what minNormTime is, we can check against normTime. Should be equal, but they might not be!

                if abs(normTime - minNormTime) < 1e-6:
                    #we found the unique solution! Yay
                    solutions.append(_solutions(x, _getY(x, angleParameter), (n+1)*2*pi - transferAngle, R1, R2, transferAngle))
                    break
                
                elif normTime < minNormTime:
                    #No solutions for n revolutions!
                    break #feel free to change this

                elif normTime < minENormTime:
                    #Two 'low path' solutions for n revolutions. Do a simple test, then just pick one.
                    try:
                        x_one = roots.brent(0, x, 1e-4, _Ftau, angleParameter, parabolicNormTime, normTime, n)
                        solutions.append(_solutions(x_one, _getY(x_one, angleParameter), n, R1, R2, transferAngle))
                    except:
                        pass
                    try:
                        x_two = roots.brent(x, 1-1e-50, 1e-4, _Ftau, angleParameter, parabolicNormTime, normTime, n)
                        solutions.append(_solutions(x_two, _getY(x_one, angleParameter), n, R1, R2, transferAngle))
                    except:
                        pass #one of these two will work; probably both

                    break

            if abs(normTime - minENormTime) < 1e-6: #the minimum energy ellipse
                    solutions.append(_solutions(0, _getY(0, angleParameter), n, R1, R2, transferAngle))

            else:
                if n > 0 or normTime > minENormTime: #high path solution - i.e., the transfer goes to apoapsis and intersects coming back down
                    try:
                        x = roots.brent(-1+MEPSILON, 0, 1e-4, _Ftau, angleParameter, parabolicNormTime, normTime, n)
                        solutions.append(_solutions(x, _getY(x, angleParameter), n, R1, R2, transferAngle))
                    except:
                        pass

                if n > 0 or normTime < minENormTime: #low path solution
                    try:
                        x = roots.brent(0, 1-MEPSILON, 1e-4, _Ftau, angleParameter, parabolicNormTime, normTime, n)
                        solutions.append(_solutions(x, _getY(x, angleParameter), n, R1, R2, transferAngle))
                    except:
                        pass
            
            minENormTime = pi*minENormTime
        
    return solutions
