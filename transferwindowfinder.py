import sys, os, matplotlib

try:
    sys.path.append('libs/') #add the libs to the system path so we can import them
    #import matrix_solver as ms #add custom modules to this line
    import vector, math, deltaV, lambert
    from body import Body, utcJulian
    from lambert import MU as MU_SUN
except ModuleNotFoundError:
    print('Please only run this from the project_transfer_windows/ directory.')
    raise ModuleNotFoundError
    sys.exit()


import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import pi, G

PATH = './configs/bodies.cfg'
AU_convert = 149.597870e9 # num of meters in AU

def calcAgain():
    calc = input('Would you like to calculate another transfer?').lower().startswith('y')
    return calc

def validateInput(prompt, category, responses=(), default='', fail="I didn't understand your input."):
    '''Validate user input. Prompt is what is shown to the user. Category is
type of response, and determines how the function validates
Valid types include:
    bool, which looks for a y/n answer and defaults to n
    float, which returns a float
    int, which returns an int
    list, which checks against possible responses
    date, which accepts a date in mm/dd/yyyy format (though it doesn't check except to see that it has a str of length 10 that consists of 3 integers separated by /)

responses is a tuple of allowed responses, and is only used by type list.
fail is a failure string to show the user. '''
    #TODO: add default behaviour handling so the user doesn't have to enter the value
    loop = True
    while loop:
        check = input(prompt)

        if category == 'list':
            if check in responses:
                loop = False
            else:
                print(fail)

        elif category == 'bool':
            check = check.lower().startswith('y')
            loop = False

        elif category == 'float':
            try:
                check = float(check)
                loop = False
            except ValueError:
                print(fail)

        elif category == 'int':
            try:
                check = int(check)
                loop = False
            except ValueError:
                print(fail)

        elif category == 'date':
            date = check.split('/')
            if len(date) == 3 and len(date[0]) == 2 and len(date[1]) == 2 and len(date[2]) == 4:
                try:
                    [int(elem) for elem in date] #we don't care about keeping it, we're just checking
                    loop = False
                except ValueError:
                    print(fail)
            else:
                print(fail)

        else: #something went wrong, probably that 'category' was not valid
            raise ValueError('Cannot validate with that category type.')

    return check
            
            
            

calc = True

while calc:

    date = validateInput('Please enter an earliest departure date in the form of dd/mm/yyyy: ', category='date')
    departureWindow = validateInput('Please enter the latest departure date as a number of days from the earliest departure date. (For example, for six months enter "180" without quotes): ', category='int')
    jules = utcJulian(date)

    validPlanets = ['Earth', 'Mars']

    origin = validateInput('Please enter a starting body. Valid entries include Earth and Mars: ', category='list', responses=validPlanets)
    startRad = validateInput('Please enter the altitude of the starting circular parking orbit in (km): ', category='float')
    #startRad = startRad
    ignorePlaneChange = validateInput('Should plane change delta-V be ignored? For example, ignore plane change delta-V if you launched into the plane of the target (y/n): ', category='bool')
    
    destination = validateInput('Please enter a destination body. Valid entries include Earth and Mars: ', category='list', responses=validPlanets)
    flyby = validateInput('Are you calculating a flyby or have some other reason to ignore insertion delta-V? (For example, you plan to aerocapture) (y/n): ', category='bool')
    if not flyby:
        destApo = validateInput('Please enter the apoapsis of the desired capture orbit, in km: ', category='float')
        #destApo = destApo
        destPeri = validateInput('Please enter the periapsis of the desired capture orbit, in km: ', category='float')
        #destPeri = destPeri
        
    
    start = Body(origin, jules, PATH)
    dest = Body(destination, jules, PATH)

    #now figure out the Hohmann time for a transfer from start to dest and use that as the centre of the transfer window
    hohmann = math.floor(deltaV.HohmannTime(MU_SUN, start.a, dest.a))
    

    #that's the centre, what's the window?

    #our minimum time of flight is the larger of either the hohmann time minus the destination's period, or half the hohmann time
    begin = max(hohmann - math.floor(dest.period), math.floor(hohmann/2))
    end = begin + min(2*math.floor(dest.period), hohmann) 
    #our maximum time of flight is the minimum time of flight plus the smaller of twice the destination period or the hohmann time
    #this accounts nicely for a destination that has both a longer period than a starting body, and a shorter period

    flightRange = end-begin + 1 #number of days between end and begin

    timesOfFlight = [begin+i for i in range(flightRange)] #the range of days we'll pass to our lambert solver for time of flight

    departureDates = [i for i in range(departureWindow + 1)] #list of departure dates, as julian days since J2000
    print('Departure dates run from {} to {}'.format(min(departureDates), max(departureDates)))
    print('Times of flight run from {} to {}'.format(min(timesOfFlight), max(timesOfFlight)))

    #now for each departure date, lambert solve for every time of flight and plug the output of the lambert solver
    #through our injection and insertion dV finders. Skip insertion dV if we're doing a flyby.
    #store each total dV in a nested list, where inner list index corresponds to time of flight, and outer list index corresponds to departure date\

    transferDVs = []
    
    for day in departureDates:
        deltaVRangeforDeparture = []
        start.update(jules+day)
        for time in timesOfFlight:
            dest.update(jules+day+time)
            solutions = lambert.sunMethod(MU_SUN, start.cartesianPos(), dest.cartesianPos(), time)
            solutionDeltaVs = []
            for array in solutions: #each element in the solutions list is a list of R1, R2, V1, V2, angle. We care about R1, V1. There may be multiple solutions for a given time-of-flight. Probably needs some work, because there are situations where it will return [0,0,0,0,0]
                if array[0] == 0: #degenerate solution
                    totalDeltaV = float('nan') #might get spot discontinuities in delta V plot...
                else:
                    injectDeltaV = deltaV.injection(array[2], start, startRad, ignorePlaneChange)
                    if not flyby:
                        insertDeltaV = deltaV.insertion(array[3], dest, destPeri, destApo, ignorePlaneChange)
                    else:
                        insertDeltaV = 0
                    
                    totalDeltaV = (injectDeltaV[0] + insertDeltaV)/1000 #we want km/s

                solutionDeltaVs.append(totalDeltaV)

            deltaVRangeforDeparture.append(min(solutionDeltaVs))

        transferDVs.append(deltaVRangeforDeparture) #now transferDVs should be a nested list, with each 'row' corresponding to a delta V for the time of flight range starting at a new departure date

    #plotting section
    
    minDV = min([item for sublist in transferDVs for item in sublist])
    maxDV = max([item for sublist in transferDVs for item in sublist])
    print('Minimum transfer delta V (km/s): {}'.format(minDV))
    print('Maximum transfer delta V (km/s): {}'.format(maxDV))

    #xx,yy = np.meshgrid(timesOfFlight, departureDates)
    dateWidth = len(departureDates)
    tofWidth = len(timesOfFlight)
    aspect = tofWidth/dateWidth

    
    
    #plt.pcolormesh(xx, yy, transferDVs, vmin=minDV, vmax=3*minDV)
    plt.imshow(transferDVs, norm=matplotlib.colors.LogNorm(), cmap=plt.cm.coolwarm, vmin=minDV, vmax=min(maxDV, minDV*10), origin='lower', extent=[timesOfFlight[0], timesOfFlight[-1], departureDates[0], departureDates[-1]])#, vmin=minDV, vmax=min(maxDV, minDV**2)
    plt.axes().set_aspect(aspect)
    plt.xlabel('Time of flight in days')
    plt.ylabel('Departure date as days from {}'.format(date))
    #formatter = matplotlib.ticker.LogFormatter(base=10, labelOnlyBase=False)
    if flyby:
        titlestring = '{}/{} transfer from {:0g}km parking orbit, flyby'.format(origin, destination, startRad)
    else:
        titlestring = '{}/{} transfer from {:0g}km parking orbit, to \n{:0g}x{:0g}km ending orbit'.format(origin, destination, startRad, destApo, destPeri)
    plt.suptitle(titlestring)
    plt.title('Minimum delta V: {:3.2}km/s'.format(minDV))

    
    cbar = plt.colorbar()
    cbar.set_label('Delta-V, km/s')
    plt.show()

    #plt.sur

    calc = calcAgain()
                
                        

                


    #TODO: figure the Hohmann time for a transfer from one circular orbit to the other
    #THEN: make a list with a range of days around that transfer time, such that:
        #[MAX, MIN] inclusive between
        #MIN = hohmann - dest.period or hohmann/2, whichever is less
        #MAX = MIN + ( 2*dest.period or hohmann, whichever is less)

    #THEN: make an array to hold delta V
    #THEN: solve the Lambert problem with each time-of-flight, leaving on the departure date
        #WHILE: reduce each solution to its delta V to save space
        #WHILE: shove those delta V into a nested list whose index represents days from the earliest departure date
    #THEN: solve and reduce for dates up to ceil(longer body's orbital period), shoving them into said nested list with index for day from earliest departure
        
    #THEN: you have a huuuuge nested list - if we're lucky
    #THEN: using matplotlib/numpy/whatever, make coordinate arrays for a colormap
    #THEN: Graph the color map. Axes are y(Time of Flight in days), x(departure date in days from earliest departure)
    #THEN: you're done
    
