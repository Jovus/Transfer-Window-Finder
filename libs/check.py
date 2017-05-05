import math
from scipy.constants import G, pi
import vector
sun = 1.989e30
MU = G*sun
AU_conversion_factor = 149.597870e9
MU = MU/(AU_conversion_factor**3)

def lambert (mu, pos1, pos2, dt, maxRevs = 0):
  # Based on Sun, F.T. "On the Minium Time Trajectory and Multiple Solutions of Lambert's Problem"
  # AAS/AIAA Astrodynamics Conference, Provincetown, Massachusetts, AAS 79-164, June 25-27, 1979
    r1 = vector.mag(pos1)
    r2 = vector.mag(pos2)

    deltaPos = vector.sub(pos2, pos1)
    c = vector.mag(deltaPos)
    m = r1 + r2 + c
    n = r1 + r2 - c

    #we want a prograde orbit
    transferAngle = math.acos(vector.dot(pos1, pos2))/ (r1*r2)

    if pos1[0]*pos2[1] - pos1[1]*pos2[0] < 0:
        transferAngle = 2*pi - transferAngle
        print('Forcing transfer angle to prograde')

    angleParameter = math.sqrt(n/m)
    if transferAngle > pi:
        angleParameter = -angleParameter
        print('Angle parameter goes negative')

    return(transferAngle, angleParameter)
