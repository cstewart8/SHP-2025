"""
Cameron Stewart - SHP
End of week 2 progress, 26/1/25
"""

import numpy as np
import math

from astroquery.jplhorizons import Horizons

def get_coordinates(file, file_size):
    """
    Function to obtain the coordinates of the centre of plate, converted to degrees.
    file = Input file, the plate catalogue
    file_size = Number of lines in the file
    coords = Array of plate coordinates in RA, Dec.
    """
    # Empty array to store coordinates
    coords = np.zeros((file_size, 2))
    # Counter, used for indexing
    x = 0

    # Resets file pointer to start, needed since file_size already used
    file.seek(0)
    # Loop over every plate
    for line in file:
        x+= 1
        
        # Converts RA of plate to degrees and stores in array
        coords[x-1,0] = (int(line[20:22]))*15 + (int(line[22:24])*15/60) + (int(line[24])*15/600)

        # Used to account for the range -90 < Dec < 90, controls addition of minutes
        if line[25] == "-":
            #Converts Dec of plate to degrees and stores in array
            coords[x-1,1] = (int(line[25:28])) - (int(line[28:30])*1/60)
        else:
            coords[x-1,1] = (int(line[25:28])) + (int(line[28:30])*1/60)

    return coords

def julian_date(date):
    """
    Function used to convert the date from yyyy/mm/dd format to Julian date.
    Follows procedure detailed in Duffet-Smith (pg. 9)
    date = Input array storing regular date in form [yyyy,mm,dd]
    JD = Output Julian date, floating point number
    """
    # Extract information from input
    y, m, d = date[0], date[1], date[2]

    # Conditions used to alter parameters, if month is early in year
    if m == (1 or 2):
        Y = y - 1
        M = m + 12
    else:
        Y = y
        M = m

    # Corrects for Gregorian calendar differences
    if y >= 1582:
        A = round(Y/100)
        B = 2 - A + round(A/4)
    else:
        B = 0

    # Keeps year as a positive number
    if Y < 0:
        C = round((365.25 * Y) - 0.75)
    else:
        C = round(365.25 * Y)

    D = round(30.6001 * (M + 1))

    # Calculates final Julian date
    JD = B + C + D + d + 1720994.5

    return JD

def convert_LST(LST, date):
    """
    Function to convert the Local Sidereal Time to Greenwich Sidereal Time
    then to Universal Time.
    LST = Input array of the local sidereal time
    date = Input array storing regular date in form [yyyy,mm,dd]
    UT = Universal Time
    """
    # LST to GST
    # Define Longitude of observation, for UKST: 149° 03' 24''
    longitude = 149 + 3/60 + 24/3600 # degrees
    # Convert longitude to hours
    lon = longitude/15

    #Convert LST to decimal hours 
    hours = LST[0] + LST[1]/60 + LST[2]/3600
    # Calculate GST by subtracting longitudinal correction from LST
    GST = hours - lon

    # Ensures that variable stays within range 0-24
    while GST < 0:
        GST += 24

    while GST > 24:
        GST -= 24
    
    # GST to UT
    # Convert current date to JD
    JD = julian_date(date)
    
    S = JD - 2451545
    T = S/36525
    To = 6.697374558 + (2400.051336 * T) + (0.000025862 * T**2)

    # Ensures that variable stays within range 0-24
    while To < 0:
        To += 24
    while To > 24:
        To -= 24

    A = GST - To

    # Ensures that variable stays within range 0-24
    while A < 0:
            A += 24
    while A > 24:
        A -= 24

    # Calculate UT as decimal hours
    UT = A * 0.9972695663

    return UT

def get_date(file, file_size):
    """
    Function used to extract the date, LST and exposure time of each plate.
    file = Input file, the plate catalogue
    file_size = Number of lines in the file
    dates = Array used to store desired information, with first element an array containing
    exposure date as [yyyy,mm,dd], second element is LST (hours) and third element the exposure
    time (hours)
    """
    # Empty array to store date of observation, LST and exposure time
    dates = np.zeros((file_size, 3), dtype=object)
    # Counter, used for indexing
    x = 0

    # Resets file pointer to start, needed since file_size already used
    file.seek(0)
    # Loop over every plate
    for line in file:
        x+= 1
        
        # Stores observation date in array as [yyyy,mm,dd]
        # Correct for fact catalog only gives year as XX, eg. not 19XX or 20XX
        if 10 <= (int(line[30:32])) <= 99:
            dates[x-1,0] = [1900 + int(line[30:32]), int(line[32:34]), int(line[34:36])]
        else:
            dates[x-1,0] = [2000 + int(line[30:32]), int(line[32:34]), int(line[34:36])]

        # Stores LST in hours, code used to deal with blank spaces in file, for hours of LST
        try:
            dates[x-1,1] = int(line[36:38]) + int(line[38:40])/60
        except ValueError:
            dates[x-1,1] = int(line[38:40])/60

        # Stores exposure time in hours, code used to deal with blank spaces in file
        try:
            dates[x-1,2] = (int(line[52:55]) + int(line[55:56])/10)/60
        except ValueError:
            dates[x-1,2] =  (int(line[55:56])/10)/60

    return dates

def convert_date():
    return test


def main():
    # Opens desired plate catalogue
    with open("catlog_ukstu.lis", "r") as file:
        # Returns number of lines in the file, used for intitialising arrays
        num_lines = sum(1 for x in file)

        date = [2009,6,19.75] # test line
        LST = [0,24, 5.23] # test line
        # Get coordinates of centre of each plate in degrees
        coords = get_coordinates(file, num_lines)
        JD = julian_date(date) # test line
        UT = convert_LST(LST, date) # test line
        dates = get_date(file, num_lines)

    #test coords are ok
    print(coords[0:3])
    print(JD)
    print(UT)
    print(dates[0])









if __name__ == "__main__":
    main()