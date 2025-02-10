"""
Cameron Stewart - SHP
"""

import numpy as np
import math
import time

from astroquery.jplhorizons import Horizons
from astropy.coordinates import SkyCoord, FK5, FK4
import astropy.units as u

# Functions to calculate sin/cos with degree input, not currently used
def cosd(degrees):
    cosx = math.cos(math.radians(degrees))
    return cosx
def sind(degrees):
    sinx = math.sin(math.radians(degrees))
    return sinx

def get_coordinates(file, file_size):
    """
    Function to obtain the number and coordinates of the centre of the plate, converted to radians.
    file = Input file, the plate catalogue
    file_size = Number of lines in the file
    coords = Array of plate coordinates in RA, Dec.
    """
    # Empty array to store coordinates
    coords = np.zeros((file_size, 3))
    # Counter, used for indexing
    x = 0

    # Resets file pointer to start, needed since file_size already used
    file.seek(0)
    # Loop over every plate
    for line in file:
        x+= 1
        
        # Store plate number in first element of array
        coords[x-1,0] = (int(line[2:7]))
        
        # Converts RA of plate to radians and stores in array
        coords[x-1,1] = ((int(line[20:22]))*15 + (int(line[22:24])*15/60) + (int(line[24])*15/600)) * np.pi/180

        # Used to account for the range -90 < Dec < 90, controls addition of minutes
        if line[25] == "-":
            #Converts Dec of plate to radians and stores in array
            coords[x-1,2] = ((int(line[25:28])) - (int(line[28:30])*1/60)) * np.pi/180
        else:
            coords[x-1,2] = ((int(line[25:28])) + (int(line[28:30])*1/60)) * np.pi/180

    return coords

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
    if 1<= m <= 2:
        Y = y - 1
        M = m + 12
    else:
        Y = y
        M = m
    
    # Corrects for Gregorian calendar differences
    if y >= 1582: # (dont need this extra because all dates after 1582) and m >= 10 and d > 18:
        A = int(Y/100)
        B = 2 - A + int(A/4)
    else:
        B = 0
    
    # Keeps year as a positive number
    if Y < 0:
        C = int((365.25 * Y) - 0.75)
    else:
        C = int(365.25 * Y)
    
    D = int(30.6001 * (M + 1))
    
    # Calculates final Julian date
    JD = B + C + D + d + 1720994.5
    
    return JD

def convert_LST(date, LST):
    """
    Function to convert the Local Sidereal Time to Greenwich Sidereal Time
    then to Universal Time.
    date = Input array storing regular date in form [yyyy,mm,dd]
    LST = Input of the local sidereal time
    UT = Universal Time
    """
    # LST to GST
    # Define Longitude of observation, for UKST: 149° 03' 24''
    longitude = 149 + 3/60 + 24/3600 # degrees
    # Convert longitude to hours
    lon = longitude/15
    
    # Calculate GST by subtracting longitudinal correction from LST
    GST = LST - lon
    
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

def convert_date(dates, file_size):
    """
    Function to convert the date of mid-observation to the corresponding Julian date
    dates = Array containing observation date, LST and exposure time parsed from plates.
    file_size = Number of lines in the file.
    julian_dates = List of Modified Julian dates of observation times of each plate
    """
    # Array to store the combined values of Date + LST + 1/2 exp_time
    combined_dates = np.zeros((file_size, 3))
    # Output array to store calculated Julian dates
    julian_dates = np.zeros(file_size)

    # Loop over each plate
    for i in range(len(dates)):
        # Convert each LST in dates array to the UT
        UT = convert_LST(dates[i,0], dates[i,1])
        dates[i,1] = UT
        
        # Assign year to first element
        combined_dates[i][0] = dates[i,0][0]
        # Assign month to second element
        combined_dates[i][1] = dates[i,0][1]
        # Assign decimal days to third element by calculating day
        # plus decimal hours of LST + half the exposure time
        combined_dates[i][2] = dates[i,0][2] + (dates[i,1] + 0.5*dates[i,2]) / 24

        # Convert to Modified Julian date and store in array
        julian_dates[i] = julian_date(combined_dates[i]) - 2400000.5

    return julian_dates

def call_horizons(object, julian_dates, file_size):
    """
    Function used to call the JPL Horizons interface, which can return the astrometric
    coordinates of Solar System bodies at a given time.
    object = Desired body that the user wants coordinates of
    julian_dates = Array of the modified Julian dates of each plate
    file_size = Number of lines in the file
    horizon_coords = Output array of returned coordinates (RA, Dec), in radians, of the object
    """
    # Initialise array to store output coordinates
    horizon_coords = np.zeros((file_size, 2))
    """
    # Loop over every Julian date
    for i in range(len(julian_dates)):
        # Call Horizons to create an object of the target body at a given date
        obj = Horizons(id=f'{object}', location='260', epochs=julian_dates[i])
        # Return the ephemerides of this body at the desired time
        eph = obj.ephemerides()
        # Append the RA and Dec of the object to the coordinate array
        horizon_coords[i,0] = eph['RA']
        horizon_coords[i,1] = eph['DEC']
    """
    chunk_size = 50

    for i in range(0, len(julian_dates), chunk_size):
        chunk = julian_dates[i:i+chunk_size]
        obj = Horizons(id=f'{object}', location='260', epochs=chunk)

        # Return ephemeris of object
        eph = obj.ephemerides()
        #print(eph['datetime_jd', 'RA', 'DEC'])
        horizon_coords[i:i+chunk_size,0] = (eph["RA"]) * np.pi/180
        horizon_coords[i:i+chunk_size,1] = (eph["DEC"]) * np.pi/180

    return horizon_coords

def convert_coords(coords):
    """
    Function to convert the returned coordinates from Horizons from the J2000 reference to B1950,
    the same reference as the UKST plate collection.
    coords = input 2D array containing RA, Dec values (radians)
    transformed_coords = Output 2D array of coordinates converted to B1950 reference
    """
    # Separate the input array into lists of RA, Dec
    ra_list = coords[:, 0]
    dec_list = coords[:, 1]

    # Initialise a SkyCoord object, defining the coordinates in the J2000 reference frame
    J2000_coords = SkyCoord(ra=ra_list*u.rad, dec=dec_list*u.rad, frame=FK5(equinox="J2000"))#, equinox="J2000")
    # Convert these coords to the B19500 reference
    B1950_coords = J2000_coords.transform_to(FK4(equinox="B1950"))
    # Returns the transformed coordinates to the input size and shape
    transformed_coords = np.column_stack([B1950_coords.ra.rad, B1950_coords.dec.rad])

    return transformed_coords

def projection(plate_coords, horizon_coords):
    """
    Function to perform a tangent plane projection of the spherical coordinates returned by Horizons, onto
    the plate. Follows the process defined in Dick 1991 paper.
    plate_coords = Array of RA, Dec of centre of plate
    horizon_coords = Returned object coordinates from Horizons, converted to B1950 reference
    projected_coords = Angular coordinate projections, angular distances from centre of plate in radians
    """
    # Initialise array to store projected coordinates
    projected_coords = np.zeros((len(plate_coords), 2))

    # Loop over each set of plate coords and horizon coords
    for i in range(len(plate_coords)):
        # Define the angles used in calculations
        A = plate_coords[i,1]
        D = plate_coords[i,2]
        alpha = horizon_coords[i,0]
        delta = horizon_coords[i,1]

        # Calculate cos(theta), useful variable
        theta = np.sin(delta)*np.sin(D) + np.cos(delta)*np.cos(D)*np.cos(alpha-A)
        
        # Calculate the (x,y) projections on the plane
        xi = (np.cos(delta)*np.sin(alpha-A))/theta
        eta = (np.sin(delta)*np.cos(D) - np.cos(delta)*np.sin(D)*np.cos(alpha-A))/theta
        
        # Add to array
        projected_coords[i,0] = xi
        projected_coords[i,1] = eta

    return projected_coords

def coordinate_comparison(plate_coords, projected_coords):
    """
    Function to test if projection coords are in plate range and hence if object expected to appear on plate.
    plate_coords = Array of RA, Dec of centre of plate
    projected_coords = Angular coordinate projections, angular distances from centre of plate in radians
    """
    high = 3.3 * np.pi/180
    low = -3.3 * np.pi/180

    count = 0
    # Loop over every plate and corresponding object coordinates
    for i in range(len(plate_coords)):
        """
        # Define upper and lower bounds on plate, defining a 6.6 x 6.6 degree area.
        x_high = plate_coords[i,1] + (3.3 * np.pi/180)
        x_low = plate_coords[i,1] - (3.3 * np.pi/180)
        
        y_high = plate_coords[i,2] + (3.3 * np.pi/180)
        y_low = plate_coords[i,2] - (3.3 * np.pi/180)

        # Check if both coordinates are in the range
        if x_low <= projected_coords[i,0] <= x_high and y_low <= projected_coords[i,1] <= y_high:
            print(f"Plate number: {int(plate_coords[i,0])}")
        """

        # Define upper and lower bounds on plate, defining a 6.6 x 6.6 degree area.

        # Check if both coordinates are in the range
        if low <= projected_coords[i,0] <= high and low <= projected_coords[i,1] <= high:
            print(f"Plate number: {int(plate_coords[i,0])}")
        else:
            count += 1
        
    if count == len(plate_coords):
        print(f"No plate matches")


def main():
    # Time how long the code takes to run for curiosity
    start_time = time.time()
    # Opens desired plate catalogue
    with open("catlog_ukstu.lis", "r") as file:
        # Returns number of lines in the file, used for intitialising arrays
        num_lines = sum(1 for x in file)
    
        # Get coordinates of centre of each plate in radians
        coords = get_coordinates(file, num_lines)
    
        dates = get_date(file, num_lines)

        JDs = convert_date(dates, num_lines)
        #JDs = [2441875.27757215, 2441878.27580884, 2441884.92134529]
        horizon_coords = call_horizons('499', JDs, num_lines)
        horizon_1950 = convert_coords(horizon_coords)

        proj = projection(coords, horizon_1950)

        coordinate_comparison(coords, proj)

    #test coords are ok
    print(coords[0:3])
    #print(JD)
    #print(UT)
    print(dates[0:3])
    print(JDs[0:3])
    print(horizon_coords[0:3])
    print(horizon_1950[0:3])
    print(proj[0:3])


    end_time = time.time()
    print(f"Execution time is {end_time - start_time} seconds")







if __name__ == "__main__":
    main()