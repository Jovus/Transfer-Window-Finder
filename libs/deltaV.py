import math, vector, body, lambert
from scipy.constants import pi, G

##M0 = 1.989e30 #mass of the Sun
##earth = body.Body('Earth',0)
##MU = G*earth.mass #units: m**3/s**2
#testing data

##R1 = (0.47326, -0.899215, 0)
##R1 = [elem*149.597870e9 for elem in R1]
##VP = (25876.6, 13759.5, 0)
##VS = (28996.2, 15232.7, 1289.2)
AU_convert = 149.597870e9 # num of meters in AU

def injection(Vs, planet, parkRadius, ignorePlaneChange=False):
    '''Given the injection velocity V1 for the transfer, the body around which we start,
and the radius, in km, of the parking orbit around which we start, return the injection
delta-V for the desired orbit (i.e. the injection delta-V to reach v1) and the injection angle gamma.

Injection delta-V is to leave the starting body. Insertion delta-V is to arrive at the ending body.
This function gives injection delta-V (and injection angle).'''

    if ignorePlaneChange: #probably because we're going to do a mid-course correction burn
        Vs = (Vs[0], Vs[1], 0)

    mu = G*planet.mass
        
    #first, make sure our units are clear
    r0 = planet.radius + (parkRadius*1000) #yes, we assume a circular parking orbit. It's close enough to true, most likely
    Rp = planet.cartesianPos(units='m')

    Vp = planet.cartesianV(units='m') #remember to add this function to Class Body
    Vdiff = vector.sub(Vs, Vp)

    v_inf = vector.mag(Vdiff) #the desired velocity at infinity - that is, the amount of hyperbolic excess velocity
    v0 = math.sqrt(v_inf**2 + (2*mu/r0))

    deltaV = v0 - math.sqrt(mu/r0)

    rp = vector.mag(Rp)

    gamma = math.acos(vector.dot(Rp,Vdiff)/(rp*v_inf))

    return(deltaV, gamma)

def insertion(Vs, planet, parkPeri, parkApo, ignorePlaneChange=False):
    '''Given the insertion velocity Vs and insertion position R2 for the transfer, the destination body,
and the periapsis and apoapsis of the ending orbit in km, calculate the insertion delta-V'''
    
    #we don't bother to change the Lambert problem to account for a miss distance because it's approximate anyway
    #instead, we'll just calculate the delta-V it takes to slow from the hyperbolic trajectory
    #to the desired orbital trajectory, assuming we can just thrust retrograde from our
    #desired periapsis and not hit the planet (which is probably true)

    #This does a little violence to the actual parameters of the ending orbit (because the periapsis we desire
    #may not actually end up as the real periapsis of the orbit) but we don't care - it'll be close enough
    #given this is all approximate anyway


    #basically, we run injection in reverse - instead of starting in an orbit and getting hyperbolic, we start hyperbolic and find an orbit
    mu = G*planet.mass

    parkApo = parkApo*1000 + planet.radius
    parkPeri = parkPeri*1000 + planet.radius #okay, true, we could do away with planet.radius since it's added to both sides, but let's just be careful

    diameter = parkPeri + parkApo #periapsis and apoapsis are around the focus, not the centre
    a = diameter/2 #semi-major axis of our desired ellipse; our elliptical orbit has the same energy as a circular orbit of radius a

    
    Vp = planet.cartesianV('m')
    if ignorePlaneChange:
        Vs = list(Vs)
        Vp = list(Vp)
        Vs[2] = 0
        Vp[2] = 0
        
    Vdiff = vector.sub(Vs, Vp)
    v_inf = vector.mag(Vdiff) #magnitude of velocity at infinity, including extra transfer velocity

    
    #now we need the necessary velocity of an orbit of radius a - which is the same as the average velocity of an orbit with semi-major axis a
    v0 = math.sqrt(v_inf**2 + (2*mu/a)) #note we assume our hyperbolic periapsis coincides with our desired periapsis; it's probably pretty close to true as long as we don't want a periapsis WAAAY far away from the planet

    deltaV = v0 - math.sqrt(mu/a)

    if deltaV <= 0:
        print('Whoops! Insertion delta-V is wrong!')

    return deltaV

def HohmannTime(mu, peri, apo, wantDays=True):
    '''Given a gravitational parameter, periapsis, and apoapsis, calculate the time required for a Hohmann transfer between the two.'''
    a = (peri+apo)/2 #a is in AU, so let's quickly convert
    a = a*AU_convert

    time = pi*math.sqrt(a**3/mu)
    days = time/86400

    if wantDays:
        return days
    else:
        return time
