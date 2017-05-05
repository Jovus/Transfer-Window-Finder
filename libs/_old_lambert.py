#THIS DOES NOT WORK and is here for posterity/ideas/editing/whatever


import numpy as np
import math
from scipy.constants import G, pi
from scipy.optimize import fsolve
import vector
M = 1.989e30 #mass of the sun
MU = G*M
days_conversion_factor = 86400
AU_conversion_factor = 149.597870e9

print(MU)

R1 = (0.473265, -0.899215, 0) #in AUs...should be fine?
R2 = (0.066842, 1.561256, 0.030948)

def relativeError(a, b):
    return abs(1 - a/b)

def sec2Days(seconds):
    '''Given a time in seconds, return a time in days. Assumes exact 24-hr days'''
    days = seconds/86400
    return days

def acot(x):
    '''Return the inverse cotangent of x'''
    return(0.5*pi - math.atan(x))

def acoth(x):
    '''Return the inverse hyperbolic cotangent of x'''
    return(0.5*math.log((x+1)/(x-1)))

def getY(x, sigma):
    '''Given the Lambert parameters x and sigma, return parameter y. From Sun's notation.'''

    y = math.sqrt(1 - sigma*sigma*(1 - x*x))
    if sigma < 0:
        y = -y

    return y

def LagrangeToF(x, y, tau, sigma):
    '''Expression of the Lagrange Time of Flight equation for the scaled space of the Sun solution method.
x is the independent variable, tau is the desired normalized time of flight, sigma is the angle parameter'''

    if abs(x - 1) < 1e-8: #use the parabolic expresison
        print('x = {}, using parabolic F(x), which corresponds to x = 1'.format(x))
        tof = (2/3)*(1-sigma**3) - tau

    elif x > 1: #hyperbolic solution
        print('x = {}, using hyperbolic F(x), which corresponds to x > 1'.format(x))
        g = math.sqrt(x**2 - 1)
        h = math.sqrt(y**2 - 1)
        tof = (-acoth(x/g) + acoth(y/h) + x*g - y*h)/(g**3) - tau

    else: #elliptic solution
        print('x = {}, using elliptic F(x), which corresponds to |x|<1'.format(x))
        g = math.sqrt(1-x**2)
        h = math.sqrt(1-y**2)
        tof = (acot(x/g) - acot(y/h) - x*g + y*h)/(g**3) - tau

    return tof


def LPrime(x, y, tau, sigma):
    '''First derivative of Time of Flight equation. For use in iteration.'''

    F = LagrangeToF(x, y, tau, sigma)
    top = 3*x*(F + tau) - 2*(1 - (sigma**3 * x/abs(y)))
    bottom = 1 - x**2

    return top/bottom
    
def LPrime2(x, y, tau, sigma):
    '''Second derivative of time of flight equation.'''

    Fp = LPrime(x, y, tau, sigma)

    top = (1 + 4*x**2)*Fp + 2*(1 - (sigma**5 * (x**3/abs(y)**3)))
    bottom = x*(1 - x**2)

    return (top/bottom)

def LagrangeToFwithYCalc(x, tau, sigma):
    '''Erase thise after testing'''
    y = getY(x, sigma)

    return LagrangeToF(x, y, tau, sigma)

##def halley(x0, f, fp, f2p, *args, tol=1e-4, maxiters=1000):
##    '''Given an initial guess x0, a function f -> f(x),
##its first and second derivative functions fp and f2p, and a tolerance tol,
##perform Halley iteration until convergence, defined by tol.
##maxiters: max number of iterations. Defaults to 1000, which should be way more than enough'''
##    #cubic convergence
##    err = 1
##    iterations = 0
##    while err > tol:
##        F = f(x0, *args)
##        FP = fp(x0, *args)
##        F2P = f2p(x0, *args)
##        x = x0 - (2*F*FP)/(2*(FP*FP) - (FP*F2P))
##        err = relativeError(x0, x)
##        
##        iterations +=1
##        x0 = x
##        if iterations >= maxiters:
##            print("Reached maximum number of iterations in Halley's method, which probably means something's wrong.")
##            break
##
##    return x0


##def ftau(x, tau, sigma, taup):
### returns the diff of desired normalTime and the calculated normalTime for our given x (given sigma)
###should be defined over (-1, inf)
##
##    if abs(x-1) <= 1e-8: #parabolic, feel free to tighten this to ==
##        tof = taup - tau
##    else:
##        y = getY(x, sigma)
##
##        if x > 1: #hyperbolic
##            g = math.sqrt(x*x - 1)
##            h = math.sqrt(y*y - 1)
##            tof = (-acoth(x/g) + acoth(y/h) + x*g - y*h)/(g**3) - tau


##    '''Given x, normalized time of flight tau, and angle parameter sigma, return F(x) (the Lagrange time of flight equation)
##This effectively returns the difference between our desired normalized time of flight and our normalized time of flight on the current iteration of x'''

##    y = getY(x, sigma)
##    if abs(x-1) <=1e-8: #parabolic case; included for completeness
##        tof =  ((2/3)*math.sqrt(1-sigma**3)-tau)
##
##    elif x > 1: #hyperbolic case
##        g = math.sqrt(x*x - 1)
##        h = math.sqrt(y*y - 1)
##        tof = (-acoth(x/g) + acoth(y/h) + x*g - y*h)/g**3 - tau #possibly wrong, should be 1-x*2
##
##    else: #elliptic case, -1 > x > 1
##        g = math.sqrt(1 - x*x)
##        h = math.sqrt(1 - y*y)
##        tof = (acot(x/g) - acot(y/h) - x*g + y*h)/g**3 - tau
##
##    return tof

##def ftaup(x, tau, sigma):
##    '''First derivative of function ftau'''
##
##    y = getY(x, sigma)
##    f = ftau(x, tau, sigma)
##    denom = 1 - x*x
##    num1 = 3*x*(f + tau)
##    num2 = 2 - 2*sigma**3 * x/abs(y)
##    num = num1 - num2
##    fp = num/denom
##    #fp = (1/(1-x*x)) * (3*x*(f+tau) - 2*(1-(sigma**3 * (x/abs(y)))))
##
##    return fp
##
##def ftau2p(x, tau, sigma):
##    '''Second derivative of function ftau'''
##
##    y = getY(x, sigma)
##    fp = ftaup(x, tau, sigma)
##
##    denom = x - x**3
##    num1 = 1 + 4*x*x*fp
##    num2 = 2 - 2*sigma**5 * x**3/abs(y)**3
##    num = num1 + num2
##    fpp = num/denom
##    #fpp = (1/(x*(1 - x*x))) * ((1 + 4*x*x)*fp + 2*(1 - (sigma**5 * (x**3/abs(y)**3))))
##    return fpp
    
def sunSolve(mu, R1, R2, dt, start, target, maxRevs=0): 
    '''Given the standard gravitational parameter mu, the starting position r1,
ending position r2, and time of flight dt, return parameter and semi-major axis of a keplerian orbit.
    R1, R2 are 3-element vectors for position, stored as tuples in the form (x,y,z) or (r,theta,z)
    units are AU
    dt is in days
    start and target are both body objects that contain orbital elements and various other information for their respective bodies
    maxRevs tells us how many total revolutions about the Sun are allowed
    prograde is a boolean that tells us whether to go prograde (the short way)
        or retrograde (the long way)

Based on Sun, F.T. "On the Minium Time Trajectory and Multiple Solutions of Lambert's Problem"
AAS/AIAA Astrodynamics Conference, Provincetown, Massachusetts, AAS 79-164, June 25-27, 1979

Very helpfully laid out in "Automated trajectory design for impulsive and low thrust
interplanetary mission analysis" by Samuel Arthur Wagner, 2014 Iowa College, included in docs'''
    #first, convert from AUs and days into meters and seconds
    R1 = [i*AU_conversion_factor for i in R1]
    R2 = [i*AU_conversion_factor for i in R2]
    dt = days_conversion_factor*dt

    r1 = vector.mag(R1)
    r2 = vector.mag(R2)
    deltaR = vector.sub(R2, R1) #change in position
    #radMu = math.sqrt(mu)

    c = vector.mag(deltaR) #chord of the section
    s = (r1 + r2 + c)/2 #semiperimeter
    m = r1 + r2 + c #aka 2*s
    n = r1 + r2 -c
    #iradM = 1/math.sqrt(m)
    #iradN = 1/math.sqrt(n)

    #we only want prograde orbits, because retrograde orbits are way expensive
    transferAngle = math.acos(vector.dot(R1, R2)/(r1*r2))
    if vector.cross(R1, R2)[2] < 0:
        print('Transforming transfer angle')
        transferAngle = 2*pi - transferAngle #keep the prograde despite the fact that it's longer

    normalTime = 4*dt*math.sqrt(mu/(m*m*m)) #Sun's method time normalization tau
    #for clarity: we want to converge on an orbit that takes this unpacked normalTime to get from R1 to R2

    
#    sigma = math.sqrt(((4*r1*r2)/(m*m))*math.cos(transferAngle/2))
    sigma = math.sqrt(n/m)
    if transferAngle > pi:
        sigma = -sigma


    angleparameter = math.sqrt(n/m)

    print('Sigma = {}, aparam = {}'.format(sigma, angleparameter))
    #independent iteration variable x: x**2 = 1 - (m/4a)

    parabolicNormTime = (2/3)*(1-(sigma**3)) #normalized time for x = 1
    minENormTime = math.acos(sigma) + sigma*(math.sqrt(1 - sigma**2)) #normalized time for x = 0

    #determine starting x based on normalTime comparisons
    print('taup = {}'.format(parabolicNormTime))
    print('tau = {}.'.format(normalTime))
    print('tauminE = {}'.format(minENormTime))
    
    if abs(normalTime - parabolicNormTime) < 1e-8: #parabolic
        x0 = 1

    elif normalTime < parabolicNormTime:
        x0 = 3

    elif abs(normalTime - minENormTime) < 1e-8: #elliptical case close to minimum energy
        x0 = 0

    elif normalTime < minENormTime:
        x0 = 0.5

    else:
        x0 = -0.5

    print('x0 is {}'.format(x0))


    #fun = lambda x: LagrangeToFwithYCalc(x0, normalTime, sigma)
    #x = fsolve(fun, x0)
    
    tol = 1e-10 #relative error in x
    res = 1
    iterations = 1
    maxiters = 500
    while res > tol:

        y = getY(x0, sigma)
        F = LagrangeToF(x0, y, normalTime, sigma)
        Fp = LPrime(x0, y, normalTime, sigma)
        Fpp = LPrime2(x0, y, normalTime, sigma)

        #x = x0 - (2*F*Fp)/((2*Fp**2) - F*Fpp) #Halley's method
        #let's use newton instead and see what happens
        x = x0 - F/Fp
        print("Just got through {} runs of Halley's method, and x is now {}".format(iterations, x))

        res = abs(x - x0)/abs(x0) #relative change in x due to iteration step
        x0 = x

        iterations += 1
        if iterations > maxiters:
            print('Dumping out of convergence loop because maxiters hit.')
            break


    
    
    y = getY(x0, sigma)
    tof = LagrangeToF(x0, y, normalTime, sigma)

    #ok, now find the velocity vector solutions
    sqrtMu = math.sqrt(mu)
    invSqrtN = 1/math.sqrt(n)
    invSqrtM = 1/math.sqrt(m)

    vc = sqrtMu + (y*invSqrtN + x*invSqrtM)
    vr = sqrtMu + (y*invSqrtN - x*invSqrtM)
    ec = [elem*(vc/c) for elem in deltaR]
    v1 = vector.add(ec, [elem*(vr/r1) for elem in R1])
    v2 = vector.sub(ec, [elem*(vr/r2) for elem in R2])

    invA = (2/r1) - ((vector.mag(v1)**2)/mu)

    a = 1/invA
    
    print('ToF = {}'.format(tof))
    print('NormTime = {}'.format(normalTime))
    print(v1, v2)
    return x, a, m, sigma

        

    































        

##    parabolicNormalTime = (2/3)*(1-sigma**3) #should be sqrt?
##
##    minENormalTime = math.acos(sigma) + sigma*math.sqrt(1-sigma*sigma)
##
##    #now we're ready to start actually solving for x, the independent variable
##    
##    #first, choose a sensible x0
##    
##    if abs(normalTime - parabolicNormalTime) <= 1e-8: #parabolic
##        x0 = 1
##    elif normalTime > parabolicNormalTime: #hyperbolic
##        x0 = 3
##    elif abs(normalTime - minENormalTime) <= 1e-8: #still elliptical
##        x0 = 0
##    elif normalTime < minENormalTime: #elliptical
##        x0 = 0.5
##    else: #yet another elliptical solution case
##        x0 = -0.5
##    func = lambda x: ftau(x, normalTime, sigma)
##    #x = fsolve(func, x0)
##    x = halley(x0, ftau, ftaup, ftau2p, normalTime, sigma, tol=1e-2)
##    a = (1-s)/2*x*x
##    return x, a
    
    

































##
##    #this algorithm works for at least some cases, but it's fragile. Try to improve
##    
##    #it's a bit unfortunate that R1, R2 are capitalized, but I need some way to distinguish them from their magnitudes
##    r1 = vector.mag(R1) #magnitude of R1
##    r2 = vector.mag(R2) #magnitude of R2
##
##    deltaR = vector.sub(R2, R1) #change in position
##    #c = vector.mag(deltaR) #magnitude of the change in position vector
##
##
###totally start over, because I can't seem to properly debug that code
##
##    #we want a transfer orbit going counter-clockwise around the Z axis, viz. prograde
##    transferAngle = math.acos(vector.dot(R1, R2)/(r1 * r2))
##    if (R1[0]*R2[1]) - (R1[1]*R2[0]) < 0: 
##        transferAngle = 2*pi - transferAngle
##
##
##    cosNu = math.cos(transferAngle)
##    sinNu = math.sin(transferAngle)
##    
##    k = r1*r2*(1-cosNu)
##    l = r1 + r2
##    m = r1*r2*(1+cosNu)
##
##    c = math.sqrt(2*m)
##    #these parameters correspond to parabolic orbits; one in each direction
##    p_min = k/(l+c)
##    p_max = k/(l-c)
##
##    #now, the rule is:
##    #if transferAngle > PI, 0 > p > p_max
##    #if transferAngle < PI, p_min > p > inf
##
##    #this rule is always satisfied if we pick p such that p_min > p > p_max
##    if transferAngle >= pi:
##        print('p cannot go above {}'.format(p_max))
##    if transferAngle < pi:
##        print('p cannot go below {}'.format(p_min))
##    p_span = p_max - p_min
##    p0 = p_min + p_span/3
##    p1 = p_min + p_span/2
##    print('Starting guesses: p0 = {}, p1 = {}'.format(p0, p1))
##
##    #given these initial guesses, calculate their corresponding times of flight
##    #p0 = 1.2
##    #p1 = 1.3
##    p = [p0, p1, 0]
##    t = [0, 0, 0]
##
##    for i in range(2):
##
##        a = m*k*p[i]/((2*m-l*l)*p[i]*p[i] + 2*k*l*p[i] - k*k)
##        f = 1 - r2/p[i]*(1-cosNu)
##        g = r1*r2*sinNu/math.sqrt(mu*p[i])
##        fdot = math.sqrt(mu/p[i])*math.tan(transferAngle/2)*((1-cosNu)/p[i] - 1/r1 - 1/r2)
##
##        if a >=0: #elliptic
##            pushterm = (r1/a)*(1-f)
##            comp = 1 - pushterm
##            dE = math.acos(comp)
##            seconds = g + math.sqrt(a**3/MU)*(dE-math.sin(dE))
##        else: #hyperbolic
##            print("Orbit has gone hyperbolic.")
##            dE = math.acosh(1 - (r1/a)*(1-f))
##            seconds = g + math.sqrt((-a)**3/MU)*(math.sinh(dE) - dE)
##
##        t[i] = sec2Days(seconds)
##
##    tol = 1e-2
##    res = 1
##    n = 1
##    print('Seed times: t0 = {}, t1 = {}'.format(t[0], t[1]))
##    while tol < res and n < 5000:
##
##        p[2] = p[1] + (dt - t[1])*(p[1] - p[0])/(t[1] - t[0]) #fragile
##
##        a = m*k*p[2]/((2*m-l*l)*p[2]*p[2] + 2*k*l*p[2] - k*k)
##        f = 1 - r2/p[2]*(1-cosNu)
##        g = r1*r2*sinNu/math.sqrt(mu*p[2])
##        fdot = math.sqrt(mu/p[2])*math.tan(transferAngle/2)*((1-cosNu)/p[2] - 1/r1 - 1/r2)
##
##        if a >=0: #elliptic
##            dE = math.acos(1 - (r1/a)*(1-f))
##            seconds = g + math.sqrt(a**3/MU)*(dE-math.sin(dE))
##        else: #hyperbolic
##            print("Orbit has gone hyperbolic")
##            dE = math.acosh(1 - (r1/a)*(1-f))
##            seconds = g + math.sqrt((-a)**3/MU)*(math.sinh(dE) - dE)
##
##        t[2] = sec2Days(seconds)
##
##        #clean up for next loop
##        res = abs(t[2] - t[1])
##       # p[0] = p[1]
##        p[1] = p[2]
##       # t[0] = t[1]
##        t[1] = t[2]
##        print('n = {0:<4} p = {1:.6f} AU a = {2:0.6f} AU t = {3:0.4f} days'.format(n, p[2], a, t[2]))
##        n +=1
##        if n == 5000:
##            print('Failed to converge; dropping out because of max iterations.')
##
##    return p, a, t
##        
##    
######    more = r1 + r2 + c
######    less = r1 + r2 - c
####    transferAngle = math.acos(vector.dot(R1, R2)/(r1 * r2)) #delta nu, or change in true anomaly (difference between starting angle and ending angle of the spacecraft). In a sense, this is the most difficult bit to understand.
####
####    if (R1[0]*R2[1]) - (R1[1]*R2[0]) < 0: #retrograde is shorter, but we want prograde anyway because retrograde transfers are beaucoup expensive
####        transferAngle = 2*pi - transferAngle
####
####    
####    #using the cosine and sine regularly, so store them
####    nuCos = math.cos(transferAngle)
####    nuSin = math.sin(transferAngle)
####    
####    #intermediate values to make reading easier. These don't have any specific meaning; they're just values we're going to use a lot
####    k = r1*r2*(1 - nuCos)
####    l = r1 + r2 #major axis of Hohmann transfer
####    m = r1*r2*(1 + nuCos)
####    
####    #now find the limits of our possible guess for the parameter of the orbit
####    p_low = k / (l + math.sqrt(2*m))
####    p_hi = k / (l - math.sqrt(2*m))
####    
####    #those are the parameters of the parabolic intersection orbits
####    #for nu < pi, p_low < p < inf
####    #for nu > pi, 0 > p > p_hi
####
####    #for the sake of simplicity, our initial p guess will be either half or double the corresponding parabolic p
####    if transferAngle <= pi:
####        print('Transfer angle less than pi, using p_low')
####        print('Transfer angle is {}'.format(math.degrees(transferAngle)))
####        p1 = p_low + (p_low+p_hi)/3
####        p0 = p_low + (p_low+p_hi)/2
####        pass
####    else:
####        print('using p_hi')
####        p1 = p_hi/2
####        p0 = p_hi - (p_hi - p1)/2 #average of guess and max
####
####    #p1 = 1.3
####    #p0 = 1.2    
####
####    #if a is positive, trial orbit is elliptic. If a is negative, trial orbit is hyperbolic (which will never be lowest delta-v transfer unless we aerobrake, but let's keep it anyway)
#### 
####    t = [0,0,0]
####    p = [p0,p1,0]
####    #for the linear interpolation method, we need two starting values. So we need to initialize some t0 as well
####    for j in range(2): #two iterations
####       
####        a = m*k*p[j]/((2*m - l*l)*p[j]*p[j] + (2*k*l*p[j] - k*k)) #semi-major axis for guess orbit
####        
####        #some intermediate values that make reading easier
####        f = 1 - (r2/p[j])*(1 - nuCos)
####        fdot = math.sqrt(mu/p[j])*math.tan(transferAngle/2)*(((1-nuCos)/p[j]) - (1/r1) - (1/r2))
####        g = (r1*r2*nuSin)/math.sqrt(mu*p[j])
####        gdot = 1 - (r1/p[j])*(1 - nuCos)
####
####        #Calculating eccentric anomaly
####
####
####        if a >= 0: #elliptic
####            costerm = 1 - (r1/a)*(1-f)
####            sinterm = (-r1)*r2*fdot/math.sqrt(mu*a)
####            try:
####                dE = math.acos(costerm)
####            except ValueError:
####                dE = math.asin(sinterm)
####            #print('Delta E = {}'.format(dE))
####
####            seconds = g + math.sqrt(math.pow(a,3)/mu)*(dE - math.sin(dE))
####
####        if a < 0: #hyperbolic
####            dE = math.acosh(1 - (r1/a)*(1 - f))
####            seconds = g + math.sqrt((-a)**3/mu)*(math.sinh(dE) - dE)
####
####        t[j] = sec2Days(seconds)
####        
####
####    #now that we've found two previous values, we can select a new p with linear interpolation and calculate its corresponding t until we're within our tolerance
####    trial = 1
####    residue = 1
####    tol = 1e-2 #time is in days; this gives a time-of-flight tolerance of ~15 minutes, which is plenty close enough
####    while residue > tol:
####
####        p[2] = p[1] + ((dt - t[1])*(p[1] - p[0]))/(t[1] - t[0])
####
####        #now calculate our latest time of flight, given our latest parameter
####
####        a = m*k*p[2]/((2*m - l*l)*p[2]*p[2] + 2*k*l*p[2] - k*k)
####        
####        f = 1 - (r2/p[2])*(1 - nuCos)
####        #fdot = math.sqrt(mu/p[2])*math.tan(transferAngle/2)*(((1-nuCos)/p[2]) - (1/r1) - (1/r2))
####        g = (r1*r2*nuSin)/math.sqrt(mu*p[2])
####        gdot = 1 - (r1/p[2])*(1 - nuCos)
####        
####
####        #this looks familiar, doesn't it? Calculating eccentric anomaly, but with p instead of p0
####        if a >= 0: #elliptic
####            costerm = 1 - (r1/a)*(1-f)
####            sinterm = (-r1)*r2*fdot/math.sqrt(mu*a)
####            try:
####                dE = math.acos(costerm)
####            except ValueError:
####                dE = math.asin(sinterm)
####                #print('Delta E = {}'.format(dE))
####        
####            seconds = g + math.sqrt(a**3/mu)*(dE - math.sin(dE))
####
####        else: #hyperbolic
####            dE = math.acosh(1 - (r1/a)*(1 - f))
####            seconds = g + math.sqrt((-a)**3/mu)*(math.sinh(dE) - dE)
####
####        t[2] = sec2Days(seconds)
####        
####        #clean up for the next iteration and check exit condition
####        
####        residue = abs(dt - t[2])
####        print('Iter = {0:<4} p = {1:.6f} AU a = {2:0.6f} AU t = {3:0.4f} days'.format(trial, p[2], a, t[2]))
####
####        t[0] = t[1]
####        t[1] = t[2]
####        p[0] = p[1]
####        p[1] = p[2]
####        trial +=1
####        if trial > 1000:
####            print("Brent's method has failed to converge.")
####            break
####    
####    return p[2], t[2], a
##
###def brentsMethod(p, t, dt):

def transferOrbitElements(p, t, a):
    '''Given the parameter, time of flight, and semi-major axis of a transfer orbit,
find its other orbital parameters.'''
    pass


        
        
