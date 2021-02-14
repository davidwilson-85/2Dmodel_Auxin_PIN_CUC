#!/usr/bin/env python

import math

def vector_to_degrees(vector_x, vector_y):

    '''
    Convert vector information (x and y pair) to degrees (0-360).
    
    1. Calculate the quadrant in which the vector lies and then
    convert to radians.
    2. Calculate the sine or cosine of the vector to know the
    distance from beginning of quadrant, and then convert to
    radians.
    3. Convert from radians to degrees.
    
    There are Python functions to covert from different units but I
    still need to perform the stept above to get correct degree values.
    
    Degrees are later used to rotate the arrow images that show the
    direction of auxin flux in images/videos.
    '''

    x, y = vector_x, vector_y

    # If there is no flux (no vector) return float '366.0'
    if x == 0 and y == 0:
        deg = '361.0'
        return deg

    # Hypotenuse
    h = math.sqrt(x**2 + y**2)

    if x >= 0 and y >= 0:
        radians_quadrant = 0
        sin = x / h

    elif x >= 0 and y < 0:
        cos = y / h
        radians_quadrant = math.pi * (1/2)

    elif x < 0 and y <= 0:
        sin = x / h
        radians_quadrant = math.pi

    elif x < 0 and y > 0:
        cos = y / h
        radians_quadrant = math.pi * (3/2)

    try:
        radians_vector = abs(math.asin(sin))
    except:
        radians_vector = abs(math.asin(cos))
    
    degrees = math.degrees(radians_vector + radians_quadrant)

    #print("sin: " + str(sin))
    print("rdn_vec: " + str(radians_vector))
    print("rdn_quad: " + str(radians_quadrant))
    print("degrees: " + str(degrees))

    return degrees


vector_to_degrees(0,-1)
