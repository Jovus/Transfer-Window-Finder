from math import sqrt

def checkLength(a, b):
    '''Given two vectors a and b, check to see if they're the same length. If not, raise an error.
If they are, return the length.'''
    #technically you can do all the operations I want to with vectors of different lengths,
    #by subbing in 0 for the elements of the shorter vector, but it's always better
    #to explicitly fail and make the user be clear instead of assuming what
    #he might want to do.

    la = len(a)
    lb = len(b)
    if lb is not la:
        raise IndexError('Vectors are not the same length!')

    return len(a)

def add(a, b, op='add'):
    '''Given two vectors a and b, return their sum (or difference)'''
    l = checkLength(a,b)

    if op=='add':
        vect = [a[i]+b[i] for i in range(l)]
    elif op=='sub':
        vect = [a[i]-b[i] for i in range(l)]
    else:
        raise ValueError('Inappropriate value to "op" argument. Valid ops are add or subtract.')

    return tuple(vect)

def sub(a, b):
    '''Given two vectors a and b, return their difference.'''
    #this is just a wrapper to make code easier to read. I suppose I could use a decorator, but I'm not yet comfortable with those
    return add(a, b, op='sub')

def dot(a, b):
    '''Given two vectors a and b, return their dot product.'''
    l = checkLength(a, b)
    return sum(a[i]*b[i] for i in range(l))

def cross(a, b):
    '''a and b are three-dimensional vectors. Returns the cross product of AxB'''
    l = checkLength(a, b)
    if l is not 3:
        raise ValueError('I can only work with 3-dimensional vectors.')


    #might do this more easily by just calling a determinant function, but this works for now
    
    i = a[1]*b[2] - a[2]*b[1] #first term of the cross product
    j = a[2]*b[0] - a[0]*b[2] #second
    k = a[0]*b[1] - a[1]*b[0] #third

    return (i, j, k)

def mag(a):
    '''Given a vector a, return its magnitude.'''
    return sqrt(sum([elem**2 for elem in a]))
