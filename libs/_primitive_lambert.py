#NEVER USE. Technically the algorithm 'lambertSun' does work for some cases, but it has problems and singularities. Kept for posterity



import numpy as np
import math
from scipy.constants import G, pi
import vector
sun = 1.989e30
MU = G*sun
AU_conversion_factor = 149.597870e9
MU = MU/(AU_conversion_factor**3)
print(MU)

R1 = (0.473265, -0.899215, 0)
R2 = (0.066842, 1.561256, 0.030948)

def sec2Days(seconds):
    '''Given a time in seconds, return a time in days. Assumes exact 24-hr days'''
    days = seconds/86400
    return days

def lambertSun(mu, R1, R2, dt, start, target, maxRevs=0): #remove a to shove in a class later, this is just for testing
    '''Given the standard gravitational parameter mu, the starting position r1,
ending position r2, and time of flight dt, return a keplerian orbit.
    R1, R2 are 3-element vectors for position, stored as tuples in the form (x,y,z) or (r,theta,z)
    units are AU
    dt is in days
    start and target are both body objects that contain orbital elements and various other information for their respective bodies
    maxRevs tells us how many total revolutions about the Sun are allowed
    prograde is a boolean that tells us whether to go prograde (the short way)
        or retrograde (the long way)

Based on Sun, F.T. "On the Minium Time Trajectory and Multiple Solutions of Lambert's Problem"
AAS/AIAA Astrodynamics Conference, Provincetown, Massachusetts, AAS 79-164, June 25-27, 1979'''

    #this algorithm works for at least some cases, but it's fragile. Try to improve
    
    #it's a bit unfortunate that R1, R2 are capitalized, but I need some way to distinguish them from their magnitudes
    r1 = vector.mag(R1) #magnitude of R1
    r2 = vector.mag(R2) #magnitude of R2

    deltaR = vector.sub(R2, R1) #change in position
    #c = vector.mag(deltaR) #magnitude of the change in position vector


#totally start over, because I can't seem to properly debug that code

    transferAngle = math.acos(vector.dot(R1, R2)/(r1 * r2))
    if (R1[0]*R2[1]) - (R1[1]*R2[0]) < 0: #retrograde is shorter, but we want prograde anyway because retrograde transfers are beaucoup expensive
        transferAngle = 2*pi - transferAngle


    cosNu = math.cos(transferAngle)
    sinNu = math.sin(transferAngle)
    
    k = r1*r2*(1-cosNu)
    l = r1 + r2
    m = r1*r2*(1+cosNu)

    c = math.sqrt(2*m)
    #these parameters correspond to parabolic orbits; one in each direction
    p_min = k/(l+c)
    p_max = k/(l-c)

    #now, the rule is:
    #if transferAngle > PI, 0 > p > p_max
    #if transferAngle < PI, p_min > p > inf

    #this rule is always satisfied if we pick p such that p_min > p > p_max
    p_span = p_max - p_min
    p0 = p_min + p_span/3
    p1 = p_min + p_span/2

    #given these initial guesses, calculate their corresponding times of flight
    #p0 = 1.2
    #p1 = 1.3
    p = [p0, p1, 0]
    t = [0, 0, 0]

    for i in range(2):

        a = m*k*p[i]/((2*m-l*l)*p[i]*p[i] + 2*k*l*p[i] - k*k)
        f = 1 - r2/p[i]*(1-cosNu)
        g = r1*r2*sinNu/math.sqrt(mu*p[i])
        fdot = math.sqrt(mu/p[i])*math.tan(transferAngle/2)*((1-cosNu)/p[i] - 1/r1 - 1/r2)

        if a >=0: #elliptic
            pushterm = (r1/a)*(1-f)
            comp = 1 - pushterm
            dE = math.acos(comp)
            seconds = g + math.sqrt(a**3/MU)*(dE-math.sin(dE))
        else: #hyperbolic
            dE = math.acosh(1 - (r1/a)*(1-f))
            seconds = g + math.sqrt((-a)**3/MU)*(math.sinh(dE) - dE)

        t[i] = sec2Days(seconds)

    tol = 1e-2
    res = 1
    n = 1
    while tol < res and n < 5000:

        p[2] = p[1] + (dt - t[1])*(p[1] - p[0])/(t[1] - t[0]) #fragile

        a = m*k*p[2]/((2*m-l*l)*p[2]*p[2] + 2*k*l*p[2] - k*k)
        f = 1 - r2/p[2]*(1-cosNu)
        g = r1*r2*sinNu/math.sqrt(mu*p[2])
        fdot = math.sqrt(mu/p[2])*math.tan(transferAngle/2)*((1-cosNu)/p[2] - 1/r1 - 1/r2)

        if a >=0: #elliptic
            dE = math.acos(1 - (r1/a)*(1-f))
            seconds = g + math.sqrt(a**3/MU)*(dE-math.sin(dE))
        else: #hyperbolic
            dE = math.acosh(1 - (r1/a)*(1-f))
            seconds = g + math.sqrt((-a)**3/MU)*(math.sinh(dE) - dE)

        t[2] = sec2Days(seconds)

        #clean up for next loop
        res = abs(t[2] - t[1])
        p[0] = p[1]
        p[1] = p[2]
        t[0] = t[1]
        t[1] = t[2]
        print('n = {0:<4} p = {1:.6f} AU a = {2:0.6f} AU t = {3:0.4f} days'.format(n, p[2], a, t[2]))
        n +=1
        if n == 5000:
            print('Failed to converge; dropping out because of max iterations.')

    return p, a, t
        
    
####    more = r1 + r2 + c
####    less = r1 + r2 - c
##    transferAngle = math.acos(vector.dot(R1, R2)/(r1 * r2)) #delta nu, or change in true anomaly (difference between starting angle and ending angle of the spacecraft). In a sense, this is the most difficult bit to understand.
##
##    if (R1[0]*R2[1]) - (R1[1]*R2[0]) < 0: #retrograde is shorter, but we want prograde anyway because retrograde transfers are beaucoup expensive
##        transferAngle = 2*pi - transferAngle
##
##    
##    #using the cosine and sine regularly, so store them
##    nuCos = math.cos(transferAngle)
##    nuSin = math.sin(transferAngle)
##    
##    #intermediate values to make reading easier. These don't have any specific meaning; they're just values we're going to use a lot
##    k = r1*r2*(1 - nuCos)
##    l = r1 + r2 #major axis of Hohmann transfer
##    m = r1*r2*(1 + nuCos)
##    
##    #now find the limits of our possible guess for the parameter of the orbit
##    p_low = k / (l + math.sqrt(2*m))
##    p_hi = k / (l - math.sqrt(2*m))
##    
##    #those are the parameters of the parabolic intersection orbits
##    #for nu < pi, p_low < p < inf
##    #for nu > pi, 0 > p > p_hi
##
##    #for the sake of simplicity, our initial p guess will be either half or double the corresponding parabolic p
##    if transferAngle <= pi:
##        print('Transfer angle less than pi, using p_low')
##        print('Transfer angle is {}'.format(math.degrees(transferAngle)))
##        p1 = p_low + (p_low+p_hi)/3
##        p0 = p_low + (p_low+p_hi)/2
##        pass
##    else:
##        print('using p_hi')
##        p1 = p_hi/2
##        p0 = p_hi - (p_hi - p1)/2 #average of guess and max
##
##    #p1 = 1.3
##    #p0 = 1.2    
##
##    #if a is positive, trial orbit is elliptic. If a is negative, trial orbit is hyperbolic (which will never be lowest delta-v transfer unless we aerobrake, but let's keep it anyway)
## 
##    t = [0,0,0]
##    p = [p0,p1,0]
##    #for the linear interpolation method, we need two starting values. So we need to initialize some t0 as well
##    for j in range(2): #two iterations
##       
##        a = m*k*p[j]/((2*m - l*l)*p[j]*p[j] + (2*k*l*p[j] - k*k)) #semi-major axis for guess orbit
##        
##        #some intermediate values that make reading easier
##        f = 1 - (r2/p[j])*(1 - nuCos)
##        fdot = math.sqrt(mu/p[j])*math.tan(transferAngle/2)*(((1-nuCos)/p[j]) - (1/r1) - (1/r2))
##        g = (r1*r2*nuSin)/math.sqrt(mu*p[j])
##        gdot = 1 - (r1/p[j])*(1 - nuCos)
##
##        #Calculating eccentric anomaly
##
##
##        if a >= 0: #elliptic
##            costerm = 1 - (r1/a)*(1-f)
##            sinterm = (-r1)*r2*fdot/math.sqrt(mu*a)
##            try:
##                dE = math.acos(costerm)
##            except ValueError:
##                dE = math.asin(sinterm)
##            #print('Delta E = {}'.format(dE))
##
##            seconds = g + math.sqrt(math.pow(a,3)/mu)*(dE - math.sin(dE))
##
##        if a < 0: #hyperbolic
##            dE = math.acosh(1 - (r1/a)*(1 - f))
##            seconds = g + math.sqrt((-a)**3/mu)*(math.sinh(dE) - dE)
##
##        t[j] = sec2Days(seconds)
##        
##
##    #now that we've found two previous values, we can select a new p with linear interpolation and calculate its corresponding t until we're within our tolerance
##    trial = 1
##    residue = 1
##    tol = 1e-2 #time is in days; this gives a time-of-flight tolerance of ~15 minutes, which is plenty close enough
##    while residue > tol:
##
##        p[2] = p[1] + ((dt - t[1])*(p[1] - p[0]))/(t[1] - t[0])
##
##        #now calculate our latest time of flight, given our latest parameter
##
##        a = m*k*p[2]/((2*m - l*l)*p[2]*p[2] + 2*k*l*p[2] - k*k)
##        
##        f = 1 - (r2/p[2])*(1 - nuCos)
##        #fdot = math.sqrt(mu/p[2])*math.tan(transferAngle/2)*(((1-nuCos)/p[2]) - (1/r1) - (1/r2))
##        g = (r1*r2*nuSin)/math.sqrt(mu*p[2])
##        gdot = 1 - (r1/p[2])*(1 - nuCos)
##        
##
##        #this looks familiar, doesn't it? Calculating eccentric anomaly, but with p instead of p0
##        if a >= 0: #elliptic
##            costerm = 1 - (r1/a)*(1-f)
##            sinterm = (-r1)*r2*fdot/math.sqrt(mu*a)
##            try:
##                dE = math.acos(costerm)
##            except ValueError:
##                dE = math.asin(sinterm)
##                #print('Delta E = {}'.format(dE))
##        
##            seconds = g + math.sqrt(a**3/mu)*(dE - math.sin(dE))
##
##        else: #hyperbolic
##            dE = math.acosh(1 - (r1/a)*(1 - f))
##            seconds = g + math.sqrt((-a)**3/mu)*(math.sinh(dE) - dE)
##
##        t[2] = sec2Days(seconds)
##        
##        #clean up for the next iteration and check exit condition
##        
##        residue = abs(dt - t[2])
##        print('Iter = {0:<4} p = {1:.6f} AU a = {2:0.6f} AU t = {3:0.4f} days'.format(trial, p[2], a, t[2]))
##
##        t[0] = t[1]
##        t[1] = t[2]
##        p[0] = p[1]
##        p[1] = p[2]
##        trial +=1
##        if trial > 1000:
##            print("Brent's method has failed to converge.")
##            break
##    
##    return p[2], t[2], a

#def brentsMethod(p, t, dt):

def transferOrbitElements(p, t, a):
    '''Given the parameter, time of flight, and semi-major axis of a transfer orbit,
find its other orbital parameters.'''
    pass


        
        
