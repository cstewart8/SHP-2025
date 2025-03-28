"""
Cameron Stewart - SHP
Project to analyse UKSTU plate collection and test for the presence
of potentially hazardous asteroids
"""

import numpy as np
import pandas as pd
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

def RA_rad_to_hms(radians):
    """
    Function to convert RA from radians to hh:mm:ss
    radians = Float, angle in radians
    RA = Output of RA in array of shape [hh,mm,ss]
    """
    # Convert to degrees
    degrees = radians * 180/np.pi
    
    # Calculate hours component
    h = degrees/15
    hours = int(h)

    # Calculate minutes component
    m = (h-hours)*60
    mins = int(m)

    # Calculate seconds component
    secs = (m - mins)*60
    # Output RA as array in desired format
    RA = [hours, mins, secs]

    return RA

def DEC_rad_to_dms(radians):
    """
    Function to convert Declination from radians to dd:mm:ss
    radians = Float, angle in radians
    DEC = Output of DEC in array of shape [dd,mm,ss]
    """
    # Calculate degrees component
    degrees = radians * 180/np.pi
    deg = int(degrees)

    # Calculate minutes component
    m = abs(degrees-deg) * 60
    mins = int(m)

    # Calculate seconds component
    secs = (m-mins)*60

    # Output RA as array in desired format
    DEC = [deg, mins, secs]

    return DEC

def get_coordinates(file, file_size):
    """
    Function to obtain the number and coordinates of the centre of the plate, converted to radians.
    file = Input file, the plate catalogue
    file_size = Number of lines in the file
    coords = Output, Array of plate coordinates in RA, Dec, plate number and minimum magnitude
    """
    # Empty array to store coordinates
    coords = np.zeros((file_size, 4))
    # Counter, used for indexing
    x = 0

    # Resets file pointer to start, needed since file_size already used
    file.seek(0)
    # Loop over every plate
    for line in file:
        # Excludes prism/grating plates
        if line[7] == "P" or line[7] == "S" or line[7] == "G" or line[7] == "T":
            continue
        # Excludes non-sidereal tracked plates
        elif "T" in line[57:60]:
            continue
        else:
            # Store plate number in first element of array
            coords[x,0] = (int(line[2:7]))
        
            # Converts RA of plate to radians and stores in array
            coords[x,1] = ((int(line[20:22]))*15 + (int(line[22:24])*15/60) + (int(line[24])*15/600)) * np.pi/180

            # Used to account for the range -90 < Dec < 90, controls addition of minutes
            if line[25] == "-":
                #Converts Dec of plate to radians and stores in array
                coords[x,2] = ((int(line[25:28])) - (int(line[28:30])*1/60)) * np.pi/180
            else:
                coords[x,2] = ((int(line[25:28])) + (int(line[28:30])*1/60)) * np.pi/180
            
            # Blank list to store plate prefix and exposure time
            plate_data = []
            plate_data.append(line[0:2])
            try:
                plate_data.append((int(line[52:55]) + int(line[55:56])/10))
            except ValueError:
                plate_data.append((int(line[55:56])/10))
            
            # Calculate minimium magnitude and store in array
            coords[x,3] = min_magnitude(plate_data)
            
            # Add to count
            x+= 1
    
    coords = coords[:x]

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
        
        # Excludes prism/grating plates
        if line[7] == "P" or line[7] == "S" or line[7] == "G" or line[7] == "T":
            continue
        # Excludes non-sidereal tracked plates
        elif "T" in line[57:60]:
            continue
        else:
            # Stores observation date in array as [yyyy,mm,dd]
            # Correct for fact catalog only gives year as XX, eg. not 19XX or 20XX
            if 10 <= (int(line[30:32])) <= 99:
                dates[x,0] = [1900 + int(line[30:32]), int(line[32:34]), int(line[34:36])]
            else:
                dates[x,0] = [2000 + int(line[30:32]), int(line[32:34]), int(line[34:36])]

            # Stores LST in hours, code used to deal with blank spaces in file, for hours of LST
            try:
                dates[x,1] = int(line[36:38]) + int(line[38:40])/60
            except ValueError:
                dates[x,1] = int(line[38:40])/60

            # Stores exposure time in hours, code used to deal with blank spaces in file
            try:
                dates[x,2] = (int(line[52:55]) + int(line[55:56])/10)/60
            except ValueError:
                dates[x,2] =  (int(line[55:56])/10)/60

            # Add to count
            x+= 1

    dates = dates[:x]

    return dates

def julian_date(date):
    """
    Function used to convert the date from yyyy/mm/dd format to Julian date.
    Follows procedure detailed in Duffet-Smith et al (pg. 9)
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

def convert_date(dates):
    """
    Function to convert the date of mid-observation to the corresponding Julian date
    dates = Array containing observation date, LST and exposure time parsed from plates.
    julian_dates = List of Modified Julian dates of observation times of each plate
    """
    # Array to store the combined values of Date + LST + 1/2 exp_time
    combined_dates = np.zeros((len(dates), 3))
    # Output array to store calculated Julian dates
    julian_dates = np.zeros(len(dates))

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

def call_horizons(object, julian_dates):
    """
    Function used to call the JPL Horizons interface, which can return the astrometric
    coordinates of Solar System bodies at a given time.
    object = Desired body that the user wants coordinates of
    julian_dates = Array of the modified Julian dates of each plate
    horizon_coords = Output array of returned coordinates (RA, Dec), in radians, of the object
    horizon_stats = Ouptut array storing the associated object properties from Horizons, errors(rads) and visual magnitude
    """
    if object == "499":
        print(f"Mars")
    else:
        print(f"{object}")
    # Initialise array to store output coordinates
    horizon_coords = np.zeros((len(julian_dates), 2))
    horizon_stats = np.zeros((len(julian_dates), 3))
    # Define number of dates to be sent to horizon at once.
    chunk_size = 50
    
    # Loop over every JD to return a Horizons object
    for i in range(0, len(julian_dates), chunk_size):
        chunk = julian_dates[i:i+chunk_size]
        obj = Horizons(id=f'{object}', location='260', epochs=chunk)

        # Return ephemeris of object
        eph = obj.ephemerides(quantities='1,9,36')
        # Convert coords from degrees to radians and add to array
        horizon_coords[i:i+chunk_size,0] = (eph["RA"]) * np.pi/180
        horizon_coords[i:i+chunk_size,1] = (eph["DEC"]) * np.pi/180

        # Convert error in horzion coordinates from arcseconds to radians and store in array
        horizon_stats[i:i+chunk_size,0] = (eph["RA_3sigma"]) * (np.pi/(180*3600))
        horizon_stats[i:i+chunk_size,1] = (eph["DEC_3sigma"]) * (np.pi/(180*3600))
        horizon_stats[i:i+chunk_size,2] = (eph["V"])
    
    return horizon_coords, horizon_stats

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
    J2000_coords = SkyCoord(ra=ra_list*u.rad, dec=dec_list*u.rad, frame=FK5, equinox="J2000")#, equinox="J2000")
    # Convert these coords to the B19500 reference
    B1950_coords = J2000_coords.transform_to(FK4(equinox="B1950"))
    # Returns the transformed coordinates to the input size and shape
    transformed_coords = np.column_stack([B1950_coords.ra.rad, B1950_coords.dec.rad])

    return transformed_coords

def gcd(coord1, coord2):
    """
    Function used to calculate the Great Circle Distance between two points on the celestial sphere.
    coord1 = Point 1 in form [RA1, Dec1]
    coord2 = Point 2 in form [RA2, Dec2]
    gcd = Great circle distance
    """
    # Calculate gcd between two points
    gcd = np.arccos(np.sin(coord1[2])*np.sin(coord2[1]) + np.cos(coord1[2])*np.cos(coord2[1])*np.cos(coord1[1] - coord2[0]))

    return gcd

def projection(plate_coords, horizon_coords, horizon_stats):
    """
    Function to perform a tangent plane projection of the spherical coordinates returned by Horizons, onto
    the plate. Follows the process defined in Dick 1991 paper.
    plate_coords = Array of RA, Dec of centre of plate
    horizon_coords = Returned object coordinates from Horizons, converted to B1950 reference
    horizon_stats = Array storing the associated object properties from Horizons, errors(rads) and visual magnitude
    projected_coords = Angular coordinate projections, angular distances from centre of plate in radians. Also stores
                       plate data such as number and limiting magnitude, with horizon data
    horizon_coords = Masked array of object coordinates from Horizons, converted to B1950 reference
    """
    # Array to store great circle separations
    circle_distances = np.zeros(len(plate_coords))

    # Calculate GCD of plate coord and horizon coord
    for i in range(len(plate_coords)):
        circle_distances[i] = gcd(plate_coords[i], horizon_coords[i])
        
    # Define a mask to filter out all instances where coords separated by angular distance greater than 10 degrees
    mask = (circle_distances) <= 10*np.pi/180
    
    # Filter coord arrays to only include closer points
    plate_coords_mask = plate_coords[mask]
    horizon_coords_mask = horizon_coords[mask]
    horizon_stats_mask = horizon_stats[mask]

    # Initialise array to store projected coordinates
    projected_coords = np.zeros((len(plate_coords_mask), 7))

    # Loop over each set of plate coords and horizon coords
    for i in range(len(plate_coords_mask)):
        # Define the angles used in calculations
        A = plate_coords_mask[i,1]
        D = plate_coords_mask[i,2]
        alpha = horizon_coords_mask[i,0]
        delta = horizon_coords_mask[i,1]
        
        # Calculate cos(theta), useful variable
        theta = np.sin(delta)*np.sin(D) + np.cos(delta)*np.cos(D)*np.cos(alpha-A)

        # Calculate the (x,y) projections on the plane
        xi = (np.cos(delta)*np.sin(alpha-A))/theta
        eta = (np.sin(delta)*np.cos(D) - np.cos(delta)*np.sin(D)*np.cos(alpha-A))/theta

        # Add plate number and projected coords, with errors, to array
        projected_coords[i,0] = plate_coords_mask[i,0]
        projected_coords[i,1] = xi
        projected_coords[i,2] = eta
        projected_coords[i,3] = horizon_stats_mask[i,0]
        projected_coords[i,4] = horizon_stats_mask[i,1]
        # Add visual magnitude and limiting magnitude to same array
        projected_coords[i,5] = horizon_stats_mask[i,2]
        projected_coords[i,6] = plate_coords_mask[i,3]
    
    return projected_coords, horizon_coords_mask

def coordinate_comparison(projected_coords, object_coords):
    """
    Function to test if projection coords are in plate range and hence if object expected to appear on plate.
    plate_coords = Array of RA, Dec of centre of plate
    projected_coords = Angular coordinate projections, angular distances from centre of plate in radians
    """
    # Catch if there are no projected coordinates available, prevents error of trying to loop over no array
    if projected_coords.shape[0] == 0:
        print("No plate matches")
        return
    
    # Define upper and lower bounds on plate, using a 6.6 x 6.6 degree area.
    lim = 3.2 * np.pi/180
    
    # Counter used to check if there are no plate hits, increases by 1 with each 'miss'
    count = 0
    # Counter used to display total number of hits
    hits = 0
    # Loop over every plate and corresponding object coordinates
    for i in range(len(projected_coords)):
        # Check if both coordinates are in the desired range
        if -lim <= projected_coords[i,1] <= lim and -lim <= projected_coords[i,2] <= lim:
            hits += 1
            
            # Tests what catalogue is being used
            if int(projected_coords[i,0]) > 20000:
                # Tests first number of identifier, used to decode the correction made to distinguish between surveys
                survey, correction, shift, plate_scale = plate_type(str(projected_coords[i,0]))
            else:
                # Standard UKST catalogue
                survey, correction, shift, plate_scale = "UKSTU", 0, 177.5, 14.9e-3

            # Calculate physical plate postions of object, using plate scale converted to mm per radian
            x_mm = (projected_coords[i,1])*plate_scale*((3600*180)/np.pi) + shift
            y_mm = (projected_coords[i,2])*plate_scale*((3600*180)/np.pi) + shift
            # Calculate error on physical plate postions
            x_err = (projected_coords[i,3])*plate_scale*((3600*180)/np.pi)
            y_err = (projected_coords[i,4])*plate_scale*((3600*180)/np.pi)
            
            # Convert the object's coordinates to desired output
            object_RA = RA_rad_to_hms(object_coords[i,0])
            object_DEC = DEC_rad_to_dms(object_coords[i,1])
            # Print details of object hit
            print(f"\n{survey} plate: {int(projected_coords[i,0])-correction} at x = {x_mm:.3f} ± {x_err:.3f} mm, y = {y_mm:.3f} ± {y_err:.3f}mm")
            print(f"Object at RA = {object_RA[0]}:{object_RA[1]}:{object_RA[2]:.3f}, DEC = {object_DEC[0]}:{object_DEC[1]}:{object_DEC[2]:.3f} in B1950 reference")
            print(f"Expected Visual Magnitude: {projected_coords[i,5]:.3f} on plate of Minimum Magnitude: {projected_coords[i,6]:.3f}")
            
        else:
            count += 1
    
    # Determine if there are any matches
    if count == len(projected_coords):
        print(f"\nNo plate matches")
    else:
        print(f"\nObject found on {hits} plates")

def min_magnitude(plate_data):
    """
    Function used to calculate the minimum magnitude that can be detected by a plate.
    Decided by it's filter, emulsion and exposure time.
    plate_data = Input list in format: [Plate Prefix, exposure time]
    magnitude = Output, minimum magnitude
    """
    # Dictionary used to store limiting magnitude and reference times for each filter prefix
    # Referenced in s5.6.3 of UKSTU Handbook
    prefix = {
        'U': (21, 180),
        'B': (21, 60),
        'J': (22.5, 60),
        'V': (21, 60),
        'O': (22, 60),
        'R': (21.5, 75),
        'I': (19.5, 90),
        "W": (19.5, 90),
        "Z": (19.5, 90)
    }

    # Assign limiting mag and reference time, based on first letter of prefix
    # (21, plate_data[1]) is the 'anything else' catch
    # .strip used to remove blank spaces and only use first letter of prefix.
    mag_lim, t_ref = prefix.get(((plate_data[0]).strip()[0]), (21, plate_data[1]))
    
    # Calculate minimum magnitude of plate
    if plate_data[1] != 0:
        # Add magnitude shift to reference limiting magnitude
        magnitude = mag_lim + 1.25*np.log10(plate_data[1]/t_ref)
    # For cases where exposure time = zero
    else:
        magnitude = 0

    return magnitude

def plate_type(identifier):
    """
    Function used to decode the correction made to the plate number
    Used to determine which survey the plate is from when using full ROE holdings
    identifier = Plate identifier, the uncorrected plate number
    """
    # Dictionary used to match identifier to survey type, correction,
    # shift to move coordinates to corner and plate scale.
    platetype = {
        '1': ('UKSTU', 1000000, 177.5, 14.9e-3),
        '2': ('POSS-I (E)', 2000000, 177.5, 14.9e-3),
        '3': ('POSS-I (O)', 3000000, 177.5, 14.9e-3),
        '4': ('POSS-II (B)', 4000000, 177.5, 14.9e-3),
        '5': ('POSS-II (R)', 5000000, 177.5, 14.9e-3),
        '6': ('POSS-II (I)', 6000000, 177.5, 14.9e-3),
        '7': ('ESO-B', 7000000, 150, 14.8e-3),
        '8': ('ESO-R', 8000000, 150, 14.8e-3)
    }
    # Matches first digit of identifier to survey type
    survey, correction, shift, plate_scale = platetype.get(identifier[0], ('None', 0, 0, 0))

    return survey, correction, shift, plate_scale

def get_csvdata(dataframe):
    """
    Function to parse in the relevant information from a plate catalogue in .csv format
    Obtains the plate identifier, exposure time and plate grade to calculate minimum magnitude,
    plate suffix to skip plates, and the coordinates of the centre of the plate.
    Dataframe = Input data, in pandas dataframe format
    coordinates = Output, array of plate identifier, coordinates and minimum magnitude
    """
    # Create array from dataframe of suffix and plate grade, fill any empty entries with zeros
    suffix = dataframe["suffix"].fillna(0).to_numpy()
    plate_grade = dataframe["grade"].fillna(0).to_numpy()
    # Convert array to str type, so common mask can be used
    plate_grade = plate_grade.astype(str)

    # Mask to exclude prism/grating plates
    mask1 = ~np.isin(suffix, ["P", "S", "G", "T"])
    # Mask to exclude non-sidereal tracked plates
    mask2 = (np.char.find(plate_grade, "T") == -1)

    # Create total mask
    mask = mask1 & mask2

    # Read in coordinates of plate centre in degrees, J2000 reference
    RA_arr = dataframe["raj2000.0deg"].to_numpy()
    Dec_arr = dataframe["dcj2000.0deg"].to_numpy()
    # Read in plate identifier
    identifier = dataframe["identifier"].to_numpy()
    # Read in prefix and exposure time
    prefix = dataframe["prefix"].fillna(0).to_numpy()
    exp_time = dataframe["expmins"].fillna(0).to_numpy()

    # Mask each array
    RA_arr = RA_arr[mask]
    Dec_arr = Dec_arr[mask]
    identifier = identifier[mask]
    prefix = prefix[mask]
    exp_time = exp_time[mask]

    # Convert plate coordinates to B1950 reference
    combined_coords = np.column_stack([RA_arr, Dec_arr]) * np.pi/180
    plate_coords = convert_coords(combined_coords)
    
    # Define array to store coordinates and plate data
    coords = np.zeros((len(RA_arr), 4))

    # Loop to add values to coordinate array
    for i in range(len(coords)):
        # Add plate identifier
        coords[i,0] = identifier[i]
        # Add coordinates to array
        coords[i,1] = plate_coords[i,0]
        coords[i,2] = plate_coords[i,1]
        
        plate_data = []
        # Append prefix and exposure time to empty list
        plate_data.append(prefix[i])
        plate_data.append(exp_time[i])
        # Calculate limiting magnitude of plate
        coords[i,3] = min_magnitude(plate_data)

    coordinates = coords[::-1]

    return coordinates

def get_MJdates(dataframe):
    """
    Function used to obtain the Modified Julian Dates of mid exposure, from a .csv file
    Dataframe = Input data, in pandas dataframe format
    julian_dates = Output, array of modified julian dates
    """
    # Create array from dataframe of suffix and plate grade
    suffix = dataframe["suffix"].fillna(0).to_numpy()
    plate_grade = dataframe["grade"].fillna(0).to_numpy()
    # Convert array to str type, so common mask can be used
    plate_grade = plate_grade.astype(str)

    # Mask to exclude prism/grating plates
    mask1 = ~np.isin(suffix, ["P", "S", "G", "T"])
    # Excludes non-sidereal tracked plates
    mask2 = (np.char.find(plate_grade, "T") == -1)
    
    # Create large mask
    mask = mask1 & mask2

    # Read in mjds from dataframe, then mask this array
    mjds = dataframe["mjdmidexp"].to_numpy()
    mjds = mjds[mask]

    julian_dates = mjds[::-1]

    return julian_dates

def main():
    # Time how long the code takes to run
    start_time = time.time()
    """
    # Code for running UKSTU catalogue only
    # Opens desired plate catalogue
    with open("catlog_ukstu.lis", "r") as file:
        # Returns number of lines in the file, used for intitialising arrays
        num_lines = sum(1 for x in file)
    
        # Get coordinates of centre of each plate in radians
        coords = get_coordinates(file, num_lines)
        # Get dates and times of plate exposure
        dates = get_date(file, num_lines)
        # Convert dates to Julian dates
        JDs = convert_date(dates)
        # Call Horizons and return object coordinates + errors
        horizon_coords, horizon_errors = call_horizons('499', JDs)
        # Convert Horzion coordinates to B1950 reference system
        horizon_1950 = convert_coords(horizon_coords)
        # Perform tangent plane projection
        proj, object_coords = projection(coords, horizon_1950, horizon_errors)
        # Test if objects are expected to appear on plate
        coordinate_comparison(proj, object_coords)
    
    """
    # Code for running full ROE holdings
    # Create dataframe from csv input file
    df = pd.read_csv("ROE-plate-holdings-srt.csv")
    # Get coordinates of centre of each plate in radians
    coords = get_csvdata(df)
    # Get array of modified julian dates of each plate
    mjds = get_MJdates(df)
    # Call Horizons and return object coordinates + errors
    horizon_coords, horizon_errors = call_horizons('499', mjds)
    # Convert Horzion coordinates to B1950 reference system
    horizon_1950 = convert_coords(horizon_coords)
    # Perform tangent plane projection
    proj, object_coords = projection(coords, horizon_1950, horizon_errors)
    # Test if objects are expected to appear on plate
    coordinate_comparison(proj, object_coords)


    end_time = time.time()
    print(f"Execution time is {end_time - start_time} seconds")







if __name__ == "__main__":
    main()