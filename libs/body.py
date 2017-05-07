from scipy.constants import pi
import configparser as cp
import vector, math



def utcJulian(date, to='Julian'):
    '''Provided a date in UTC in form dd/mm/yyyy, convert to days since J2000. Alternatively, given a Julian Day, convert to a date in UTC.
    Please note this is horribly, horribly inaccurate and isn't intended to have any real semblance to reality.'''
    if to == 'Julian':
        day, month, year = date.split('/')
        day = int(day)
        month = int(month)
        year = int(year) - 2000
        jules = day-1 + 30*(month-1) + year*365 #every month has exactly 30 days, right?
        return jules
    #if to == 'UTC':

class Body():
    '''Name is a heading in a cfg/ini file given in path'''

    def __init__(self, name, jules, path='/configs/bodies.cfg'):
        self.define(name, jules, path)
        self.E = self.eccentricAnomaly(self.e, self.mean_anomaly)
        self.jules = jules #the Julian date for this body definition
        pass #shove in the other Body functions once I've defined them

    def define(self, name, jules, path): #jules is julian days from J2000
        cfg = cp.ConfigParser()
        cfg.read(path)
        self.mass = float(cfg[name]['mass'])
        self.radius = float(cfg[name]['radius'])
        self.period = float(cfg[name]['period'])
        self.a = float(cfg[name]['semimajor axis'])
        self.e = float(cfg[name]['eccentricity'])
        self.inc = math.radians(float(cfg[name]['inclination']))
        self.L2000 = math.radians(float(cfg[name]['mean longitude J2000']))
        self.L = self.mean_long(self.L2000, jules, self.period)
        self.wbar = math.radians(float(cfg[name]['lperi']))
        self.lan = math.radians(float(cfg[name]['lan']))
        self.w = self.wbar - self.lan #argument of periapsis
        self.mean_anomaly = self.M(self.L, self.wbar)
        
    def eccentricAnomaly(self, e, M):
        '''Find the eccentric anomaly from the eccentricity and the mean anomaly. Iterative.'''
        tol = 10e-6
        #solving Kepler's equation, M = E - esin(E) by Newton's Method

        E = M #starting guess
        dE = 1
        while dE > tol:
            dE = (E - e*math.sin(E) - M)/(1 - e*math.cos(E))
            E -= dE

        return E
        
    def mean_long(self, L2000, jules, period):
        '''Given a mean longitude at J2000, julian days since that date, and period, return
current mean longitude'''

        if jules > 0:
            while jules >=period:
                jules -=period
        elif jules < 0:
            while jules <= period:
                jules += period

        L = L2000 + 2*pi*(jules/period)
        return L

    def M(self, L, wbar):
        mean_anomaly = L - wbar
        while mean_anomaly < -pi:
            mean_anomaly += 2*pi
        while mean_anomaly > pi:
            mean_anomaly -= 2*pi

        return mean_anomaly

    def update(self, jules):
        '''Given Julian days, update the various changing parameters. Useful so you don't have to instanstiate a whole other body'''

        self.L = self.mean_long(self.L2000, jules, self.period)
        self.mean_anomaly = self.M(self.L, self.wbar)
        self.E = self.eccentricAnomaly(self.e, self.mean_anomaly)

    def cartesianPos(self,units='AU'):
        '''Return the position of the body in X, Y, Z coordinates'''

        #let's make it easier to read than having 'self' everywhere
        a = self.a
        E = self.E
        e = self.e
        w = self.w
        inc = self.inc
        lan = self.lan
        
        #we're going to be tricksy and use geometric rotations

        #I and J form a 2d coordinate system in the plane of the orbit
        I = a * (math.cos(E) - e)
        J = a * math.sin(E) * math.sqrt(1 - e**2)

        #now we rotate these into the plane of the ecliptic and record z

        x = math.cos(w)*I - math.sin(w)*J
        y = math.sin(w)*I - math.cos(w)*J #converting from polar to cartesian using the argument of periapsis

        #now by inclination

        z = math.sin(inc)*x
        x = math.cos(inc)*x

        #now by the longitude of the ascending node
        x0 = x
        x = math.cos(lan)*x0 - math.sin(lan)*y
        y = math.sin(lan)*x0 + math.cos(lan)*y

        if units=='m':
            x,y,z = 149.597870e9*x, 149.597870e9*y, 149.597870e9*z
        return x, y, z

    def cartesianV(self, units='m'):
        '''Return the velocity of the body in X, Y, Z coordinates'''

        #we could do things the hard and correct way to calculate velocity, but instead we're going to cheat
        #central difference approximation
        jules = self.jules
        self.update(jules-1)
        pos1 = self.cartesianPos(units)

        self.update(jules+1)
        pos2 = self.cartesianPos(units)

        deltaPos = vector.sub(pos2, pos1)
        velocity = [elem/(2*86400) for elem in deltaPos]

        #leave things the way we found them
        self.update(jules)

        return velocity
