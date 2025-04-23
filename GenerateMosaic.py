#
# © Copyright Scott Sumner-Moore 2025 All Rights Reserved
#
# To run this file, you need to have the following packages installed:
# pip install geopy
# pip install geocoder
# pip install astropy
# pip install requests
# pip install urllib3
# pip install matplotlib


import math
from math import ceil
from datetime import datetime, timedelta
import requests
import ssl
from requests.adapters import HTTPAdapter
from urllib3.poolmanager import PoolManager
import json
import xml.etree.ElementTree as ET
from astropy import units as u
from astropy.coordinates import SkyCoord, Distance, EarthLocation, AltAz
from astropy.time import Time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import geocoder
from geopy.geocoders import Nominatim

# Convert RA (hms) to decimal degrees
def dms_to_deg(d, m, s):
    return round((abs(d) + m/60 + s/3600), 4) * (-1 if d < 0 else 1)

# Convert DMS to decimal degrees
def hms_to_dec(h, m, s):
    dechour = round((h + m/60 + s/3600), 4)
    if (dechour < 0):
        dechour += 24
    return dechour

# Convert RA (hms) and Dec (dms) to Unistellar link
def ra_dec_to_unistellar_link(ra_h, ra_m, ra_s, dec_d, dec_m, dec_s, scitag):
    # unistellar://science/defense?ra=12&dec=-30&fov=15&fovunit=arcmin&tags=NGC+6946
    ra = hms_to_dec(ra_h, ra_m, ra_s)
    dec = dms_to_deg(dec_d, dec_m, dec_s)
    return f"unistellar://science/defense?ra={ra}&dec={dec}&sci={scitag}"

# Convert decimal RA/Dec to Unistellar link
def generate_unistellar_link(ra, dec, tag):
    # unistellar://science/defense?ra=12&dec=-30&fov=15&fovunit=arcmin&tags=NGC+6946
    return f"unistellar://science/defense?ra={ra}&dec={dec}&sci={tag}"

# Convert decimal degrees to HMS. Return hours, minutes, seconds
def degrees_to_hms(degrees):
    degSign = 1
    if degrees < 0:
        degSign = -1
        degrees = abs(degrees)
    hours = int(degrees / 15)
    minutes = int((degrees - (hours * 15)) * 4)
    seconds = (degrees - (hours * 15) - (minutes / 4)) * 240
    return hours*degSign, minutes, seconds

# Convert decimal degrees to DMS. Return degrees, minutes, seconds
def degrees_to_dms(degrees):
    degSign = 1
    if degrees < 0:
        degSign = -1
        degrees = abs(degrees)
    deg = int(degrees)
    minutes = int((degrees - deg) * 60)
    seconds = (degrees - deg - minutes / 60) * 3600
    return deg*degSign, minutes, seconds

# Get the user's location, returning latitude and longitude
def get_location():
    # First, ask the user if they want to provide an address, a specific latitude and longitude, or use their current location
    while True:
        print("Where will you be observing from? Please choose one of the following:")
        print("  1) Provide an address")
        print("  2) Enter a latitude and longitude")
        print("  3) Use my current location (default)")
        choice = input()
        if choice == '1':
            break
        elif choice == '2':
            break
        elif choice == '3':
            break
        else:
            print("Assuming current location")
            break
    if choice == '1':
        geolocator = Nominatim(user_agent="GenerateMosaic")
        # Loop until a valid address is provided
        while True:
            try:
                # Ask the user for an address
                address = input("Enter the address: ")
                g = geolocator.geocode(address)
                latitude = g.latitude
                longitude = g.longitude
                break
            except Exception as e:
                print("Error:", e)
        
    elif choice == '2':
        while True:
            try:
                latitude = float(input("Enter your latitude: "))
                longitude = float(input("Enter your longitude: "))
                break
            except ValueError:
                print("Invalid latitude or longitude. Please enter numeric values.")
    else: # choice == '3':
        g = geocoder.ip('me')
        if g.ok:
            latitude = g.latlng[0]
            longitude = g.latlng[1]
            #address = g.address
            #city = g.city
            #state = g.state
            #country = g.country
        else:
            print("Unable to retrieve location.")

    return latitude, longitude

# Convert decimal RA/Dec to decimal Alt/Az
def convert_radec_to_altAz(ra_obj, dec_obj, altaz):
    object_location = SkyCoord(ra_obj, dec_obj, frame='icrs', unit='deg')
    object_location_alt_az = object_location.transform_to(altaz)
    altitude = object_location_alt_az.alt.degree
    azimuth = object_location_alt_az.az.degree
    return (altitude, azimuth)

# Convert decimal Alt/Az to decimal RA/Dec
def convert_altAz_to_radec(alt_obj, az_obj, altaz):
    object_location = SkyCoord(alt=alt_obj, az=az_obj, frame='altaz', unit='deg', obstime=altaz.obstime, location=altaz.location)
    object_radec = object_location.transform_to('icrs')
    ra_degrees = object_radec.ra.degree
    dec_degrees = object_radec.dec.degree
    return (ra_degrees, dec_degrees)

# Get the center of the mosaic in Alt/Az (return altitude and azimuth)
def get_center_of_mosiac(alt_az_list):
    # alt_az_list is a list of (alt, az) pairs
    # Get center of altitude and azimuth
    altitudes = [alt for alt, az in alt_az_list]
    azimuths = [az for alt, az in alt_az_list]
    altitude_center = (max(altitudes) + min(altitudes)) / 2
    azimuth_center = (max(azimuths) + min(azimuths)) / 2
    return (altitude_center, azimuth_center)

# Get the width of the mosaic in degrees azimuth
def get_width_of_mosaic(alt_az_list):
    # alt_az_list is a list of (alt, az) pairs
    # Get width of azimuth
    azimuths = [az for alt, az in alt_az_list]
    azimuth_width = max(azimuths) - min(azimuths)
    return azimuth_width

# Get the height of the mosaic in degrees altitude
def get_height_of_mosaic(alt_az_list):
    # alt_az_list is a list of (alt, az) pairs
    # Get width of altitude
    altitudes = [alt for alt, az in alt_az_list]
    altitude_width = max(altitudes) - min(altitudes)
    return altitude_width

# Generage the mosaic links given the following parameters:
# object_list: list of object names
# ra_list: list of RA values in degrees
# dec_list: list of Dec values in degrees
# size_list: list of object sizes in degrees
# fov_width_asec: field of view width in arcseconds
# fov_height_asec: field of view height in arcseconds
# overlap: overlap percentage between tiles (0-1)
# border: border percentage around objects (0-1)
# tag: tag to identify the mosaic
# latitude: latitude of the observing location in degrees
# longitude: longitude of the observing location in degrees
# observingDate: date and time of the observation
# ax: matplotlib axis for RA/Dec plot
# ax_altaz: matplotlib axis for Alt/Az plot
# output_file: name of the file to save the links to (optional)
# output_type: type of output file (text, csv or html) (optional)
def generate_mosaic(object_list, ra_list, dec_list, size_list, fov_width_asec, fov_height_asec, overlap, border, tag, latitude, longitude, observingDate, ax, ax_altaz, output_file=None, output_type='text'):

    # Convert each coordinate in the ra_list / dec_list to alt / az
    observerLocation = EarthLocation.from_geodetic(lat=latitude*u.deg, lon=-longitude*u.deg)
    altaz = AltAz(obstime=observingDate, location=observerLocation)

    alt_az_list = []
    for i in range(len(object_list)):
        alt_az_list.append(convert_radec_to_altAz(ra_list[i], dec_list[i], altaz))
        #print(f"Object {i+1}: Alt: {alt_az_list[i][0]}, Az: {alt_az_list[i][1]}, RA: {ra_list[i]}, Dec: {dec_list[i]}")

    # Determine the center of the objects in Alt/Az
    centerOfMosaicAltAz = get_center_of_mosiac(alt_az_list)
    center_alt = centerOfMosaicAltAz[0]
    center_az  = centerOfMosaicAltAz[1]

    # Get the center of the mosaic in RA/Dec
    centerOfMosaicRaDec = convert_altAz_to_radec(center_alt, center_az, altaz)
    ra_deg = centerOfMosaicRaDec[0]
    dec_deg = centerOfMosaicRaDec[1]

    # Calculate the width of the mosaic based on the object's Alt Az coordinates
    obj_width_degrees = get_width_of_mosaic(alt_az_list)

    # Calculate the height of the mosaic based on the object's Alt Az coordinates
    obj_height_degrees = get_height_of_mosaic(alt_az_list)

    fov_height_degrees = fov_height_asec / 3600.0
    fov_width_degrees = fov_width_asec / 3600.0

    # If only one object in the list, use the size of the object
    if len(alt_az_list) == 1:
        try:
            obj_width_degrees = size_list[0][0]
            obj_height_degrees = size_list[0][1]
        except TypeError:
            obj_width_degrees = 0
            obj_height_degrees = 0

        # If the object size is not provided, use the default size
        if obj_width_degrees == 0:
            obj_width_degrees = fov_width_degrees / 4
        if obj_height_degrees == 0:
            obj_height_degrees = fov_height_degrees / 4

    ## Calculate the effective size of each tile accounting for the overlap
    tile_width_degrees = fov_width_degrees * (1 - overlap)
    tile_height_degrees = fov_height_degrees * (1 - overlap)

    # Calculate the size of the object in arcseconds with the border
    obj_width_degrees += obj_width_degrees * border * 2
    obj_height_degrees += obj_height_degrees * border * 2

    # Calculate the number of tiles needed to cover the adjusted object size
    az_tiles = ceil(obj_width_degrees / tile_width_degrees)
    alt_tiles = ceil(obj_height_degrees / tile_height_degrees)

    # Get the center of the mosaic alt / az
    centerOfMosaicAlt = center_alt #centerOfMosaicAltAz.alt.degree
    centerOfMosaicAz = center_az #centerOfMosaicAltAz.az.degree

    # Determine the center of the upper left tile in degrees
    ul_center_alt = centerOfMosaicAlt - ((alt_tiles-1) * tile_height_degrees / 2)
    ul_center_az = centerOfMosaicAz - ((az_tiles-1) * tile_width_degrees / 2)
    
    ra_graph = []
    dec_graph = []
    alt_graph = []
    az_graph = []

    if (output_file):
        f = open(output_file, 'a')
        if output_type == 'html':
            # Write the header for HTML output
            f.write("<html><head></head><body>\n")
            f.write(f"<p>Mosaic for {tag} (centered on {ra_deg}, {dec_deg}) ({az_tiles}x{alt_tiles}) tiles of size {tile_width_degrees * 3600}' x {tile_height_degrees * 3600} arcsec</p>")
            f.write("<table border='1'>\n")
            f.write("<tr><th></th><th>Link</th><th>RA (hms)</th><th>Dec (dms)</th></tr>\n")
        elif output_type == 'csv':
            # Write the header for CSV output
            f.write(f"Mosaic for {tag} (centered on {ra_deg}, {dec_deg}) ({az_tiles}x{alt_tiles}) tiles of size {tile_width_degrees * 3600}' x {tile_height_degrees * 3600} arcsec\n")
            f.write("Row,Col,Link,RA (hms),Dec (dms)\n")
        else:
            # Write the header for text output
            f.write(f"Mosaic for {tag} (centered on {ra_deg}, {dec_deg}) ({az_tiles}x{alt_tiles}) tiles of size {tile_width_degrees * 3600}' x {tile_height_degrees * 3600} arcsec\n")

    else:
        f = None


    # First, we determine the centers of each tile in Alt/Az, then convert to RA/Dec
    # The list of RA/Dec centers will be used to build the links for the observations, but
    # as the observation goes on, the RA/Dec positions will shift due to the rotation of the Earth.
    # How do we correct for that?

    # Determine the center of each tile in degrees
    # Save centers in two lists for RA/Dec and Alt/Az
    old_rotation_angle = 0
    xtile_centers_altaz = []
    xtile_centers_radec = []
    for j in range(alt_tiles - 1, -1, -1):
        for i in range(az_tiles):
            xaz_deg = ul_center_az + (i * tile_width_degrees)
            xalt_deg = ul_center_alt + (j * tile_height_degrees)
            xtile_centers_altaz.append((xalt_deg, xaz_deg))
            xtile_ra, xtile_dec = convert_altAz_to_radec(xalt_deg, xaz_deg, altaz)
            xtile_centers_radec.append((xtile_ra, xtile_dec))
            xlink = generate_unistellar_link(xtile_ra, xtile_dec, tag)

            # Output the link for this tile based on the output type and file
            xtile_rah, xtile_ram, xtile_ras = degrees_to_hms(xtile_ra)
            xtile_rah_string = f"{xtile_rah}h{xtile_ram}m{xtile_ras}s"
            xtile_decd, xtile_decm, xtile_decs = degrees_to_dms(xtile_dec)
            xtile_decd_string = f"{xtile_decd}d{xtile_decm}m{xtile_decs}s"
            if f:
                if output_type == 'html':
                    f.write(f"<tr><td><input type='checkbox'></input></td><td><a href='{xlink}'>Row {alt_tiles - j} Col {i+1}</a></td><td>{xtile_rah_string}</td><td>{xtile_decd_string}</td></tr>\n")
                elif output_type == 'csv':
                    f.write(f"{alt_tiles - j},{i+1},{xlink},{xtile_rah_string},{xtile_decd_string}\n")
                else:
                    f.write(f"{xlink} {xtile_rah_string} {xtile_decd_string}\n")
            else:
                print(f"Row {alt_tiles - j} Col {i+1}: {xlink} ra: {xtile_rah_string} dec: {xtile_decd_string}")

            # Now add the tile center to the graph arrays
            ra_graph.append(xtile_ra)
            dec_graph.append(xtile_dec)
            alt_graph.append(xalt_deg)
            az_graph.append(xaz_deg)

            # Calculate the lower left corner of the tile in alt-az
            xlower_left_alt = xalt_deg - (fov_height_degrees/2)
            xlower_left_az = xaz_deg - (fov_width_degrees/2)

            # Calculate the upper left corner of the tile in Alt/Az
            xupper_left_alt = xlower_left_alt + fov_height_degrees
            xupper_left_az = xlower_left_az

            # Calculate the upper right corner of the tile in Alt/Az
            xupper_right_alt = xlower_left_alt + fov_height_degrees
            xupper_right_az = xlower_left_az + fov_width_degrees

            # Calculate the lower right corner of the tile in Alt/Az
            xlower_right_alt = xlower_left_alt
            xlower_right_az = xlower_left_az + fov_width_degrees

            # Calculate the lower left corner of the tile in RA/Dec
            xlower_left_ra, xlower_left_dec = convert_altAz_to_radec(xlower_left_alt, xlower_left_az, altaz)

            # Calculate the upper left corner of the tile in RA/Dec
            xupper_left_ra, xupper_left_dec = convert_altAz_to_radec(xupper_left_alt, xupper_left_az, altaz)

            # Calculate the upper right corner of the tile in RA/Dec
            xupper_right_ra, xupper_right_dec = convert_altAz_to_radec(xupper_right_alt, xupper_right_az, altaz)

            # Calculate the lower right corner of the tile in RA/Dec
            xlower_right_ra, xlower_right_dec = convert_altAz_to_radec(xlower_right_alt, xlower_right_az, altaz)

            # Now add the tile to the graph by drawing a closed polygon
            xpolygon = patches.Polygon([[xlower_left_ra, xlower_left_dec], [xupper_left_ra, xupper_left_dec], [xupper_right_ra, xupper_right_dec], [xlower_right_ra, xlower_right_dec], [xlower_left_ra, xlower_left_dec]], closed=True, edgecolor='blue', facecolor='none', linewidth=1)
            ax.add_patch(xpolygon)
            ax.annotate(f"({alt_tiles - j},{i+1})", (xlower_left_ra, xlower_left_dec), textcoords="offset points", xytext=(0, 10), ha='center')

            xpolygon_altaz = patches.Polygon([[xlower_left_az, xlower_left_alt], [xupper_left_az, xupper_left_alt], [xupper_right_az, xupper_right_alt], [xlower_right_az, xlower_right_alt], [xlower_left_az, xlower_left_alt]], closed=True, edgecolor='blue', facecolor='none', linewidth=1)
            ax_altaz.add_patch(xpolygon_altaz)
            ax_altaz.annotate(f"({alt_tiles - j},{i+1})", (xlower_left_az, xlower_left_alt), textcoords="offset points", xytext=(0, 10), ha='center')

     # Close the file
    if f:
        if output_type == 'html':
            f.write(f"</body></html>\n")
        f.close()
    

    xcoordinates = SkyCoord(ra=ra_graph*u.degree, dec=dec_graph*u.degree, frame='icrs')
    xcoordinates_altaz = SkyCoord(alt=alt_graph*u.degree, az=az_graph*u.degree, frame='altaz', obstime=altaz.obstime, location=altaz.location)
    
    ## Plotting
    ax.scatter(xcoordinates.ra, xcoordinates.dec, marker='o', color='blue', label='Centers of Tiles')
    ax.set_xlabel('Right Ascension (degrees)')
    ax.set_ylabel('Declination (degrees)')
    ax.set_title('RA and Dec Plot')
    #ax.legend()
    ax.invert_xaxis()
    ax.grid(True)

    ax_altaz.scatter(xcoordinates_altaz.az, xcoordinates_altaz.alt, marker='o', color='blue', label='Centers of Tiles')
    ax_altaz.set_ylabel('Altitude (degrees)')
    ax_altaz.set_xlabel('Azimuth (degrees)')
    ax_altaz.set_title('Alt Az Plot')
    #ax_altaz.legend(bbox_to_anchor=(-0.05, 1), loc='upper left', borderaxespad=0.5)
    ax_altaz.grid(True)

    xra_center_list = [ra_deg]
    xdec_center_list = [dec_deg]
    ax.scatter(xra_center_list, xdec_center_list, marker='x', color='green', label='Center of Mosaic')

    xalt_center_list = [center_alt]
    xaz_center_list = [center_az]
    ax_altaz.scatter(xaz_center_list, xalt_center_list, marker='x', color='green', label='Center of Mosaic')

    return

# Invoke the NED database to get the coordinates and size of the object
class LibreSSLAdapter(HTTPAdapter):
    """
    A custom HTTPAdapter to use a specific SSL context compatible with LibreSSL.
    """
    def __init__(self, ssl_context=None, **kwargs):
        self.ssl_context = ssl_context or ssl.create_default_context()
        super().__init__(**kwargs)

    def init_poolmanager(self, *args, **kwargs):
        kwargs['ssl_context'] = self.ssl_context
        return super().init_poolmanager(*args, **kwargs)

def query_ned_object(object_name):
    """
    Query the NED Object Lookup web service for an astronomical object.

    Args:
        object_name (str): The name of the astronomical object to look up.

    Returns:
        dict: A dictionary containing the RA, Dec, and size of the object if found.
    """
    base_url = "https://ned.ipac.caltech.edu/srs/ObjectLookup"
    #params = {
    #    "of": "xml_main",  # Request XML response
    #    "objname": object_name
    #}
    #params = 'json={"name":{"v":"M31"}}'
    #params = 'json={"name":{"v":"Stephan\'s Quintet"}}'
    params = 'json={"name":{"v":"' + object_name + '"}}'

    # Create a custom SSL context compatible with LibreSSL
    ssl_context = ssl.create_default_context()

    # Create a session and mount the custom adapter
    session = requests.Session()
    session.mount("https://", LibreSSLAdapter(ssl_context=ssl_context))

    try:
        response = session.get(base_url, params=params)
        response.raise_for_status()  # Raise an exception for HTTP errors

        # Parse the JSON string
        json_data = response.text
        data = json.loads(json_data)

        # Extract RA, Dec, and size from the XML response
        query_time = data.get("QueryTime")
        copyright_info = data.get("Copyright")
        version = data.get("Version")
        supplied_name = data.get("Supplied")
        interpreted_name = data.get("Interpreted", {}).get("Name")
        preferred_name = data.get("Preferred", {}).get("Name")
        position = data.get("Preferred", {}).get("Position", {})
        ra = position.get("RA")
        dec = position.get("Dec")
        unc_semi_major = position.get("UncSemiMajor")
        unc_semi_minor = position.get("UncSemiMinor")
        pos_angle = position.get("PosAngle")
        ref_code_position = position.get("RefCode")
        obj_type = data.get("Preferred", {}).get("ObjType", {}).get("Value")
        redshift = data.get("Preferred", {}).get("Redshift", {})
        redshift_value = redshift.get("Value")
        redshift_uncertainty = redshift.get("Uncertainty")
        redshift_ref_code = redshift.get("RefCode")
        redshift_quality_flag = redshift.get("QualityFlag")
        status_code = data.get("StatusCode")
        result_code = data.get("ResultCode")

        # Print the extracted information
        #print(f"Query Time: {query_time}")
        #print(f"Copyright: {copyright_info}")
        #print(f"Version: {version}")
        #print(f"Supplied Name: {supplied_name}")
        #print(f"Interpreted Name: {interpreted_name}")
        #print(f"Preferred Name: {preferred_name}")
        #print(f"Position:")
        #print(f"  RA: {ra}")
        #print(f"  Dec: {dec}")
        #print(f"  Uncertainty (Semi-Major): {unc_semi_major}")
        #print(f"  Uncertainty (Semi-Minor): {unc_semi_minor}")
        #print(f"  Position Angle: {pos_angle}")
        #print(f"  Reference Code (Position): {ref_code_position}")
        #print(f"Object Type: {obj_type}")
        #print(f"Redshift:")
        #print(f"  Value: {redshift_value}")
        #print(f"  Uncertainty: {redshift_uncertainty}")
        #print(f"  Reference Code: {redshift_ref_code}")
        #print(f"  Quality Flag: {redshift_quality_flag}")
        #print(f"Status Code: {status_code}")
        #print(f"Result Code: {result_code}")
        return {
            "RA": ra,
            "Dec": dec,
            "Size": unc_semi_major
        }

    except requests.exceptions.RequestException as e:
        print(f"Error querying NED: {e}")
        return None
    except ET.ParseError as e:
        print(f"Error parsing NED response: {e}")
        return None




# This script generates a mosaic of astronomical objects using the Unistellar telescope.
# It will prompt the user for the following information:
# 1. The names of the objects to include in the mosaic
# 2. The type of telescope being used
# 3. The amount of overlap between tiles
# 4. The amount of border around the objects
# 5. The desired aspect ratio of the mosaic
# 6. The name of the file (if any) to save the links to, along with the file format (html, csv, or text)
# 7. The user's location (address, specific latitude and longitude, or current location)
# 8. The date and time of the observation
#
# The script will then use the NED API to get information about each object and generate a mosaic based on that information.
# The script will then generate a list of links that will allow the user to create a mosaic of the objects using the Unistellar telescope.
# This script is designed to be run from the command line


# Need to get a list of objects and their sizes and generate the mosaic to cover them all
object_list = []  # List of object names
size_list = []  # List of object sizes
ra_list = []  # List of RA values
dec_list = []  # List of Dec values
alt_list = []  # List of Alt values
az_list = []  # List of Az values

# Ask the user for the object names, one at a time, until they are done
while True:
    object_name = input("Enter the name of an object (or 'done' or just hit enter to finish): ")
    if object_name.lower() == "done":
        break
    # if object name is empty, done
    if not object_name:
        break
    object_list.append(object_name)
    # Query NED for each object and get its RA, Dec, and size
    result = query_ned_object(object_name)
    if result:
        if (result["RA"] is None) or (result["Dec"] is None) or (result["Size"] is None):
            print(f"Failed to retrieve data for object: {object_name}")
            # remove object name from list
            object_list.pop()
        else:
            ra_list.append(result["RA"])
            dec_list.append(result["Dec"])
            size_list.append(result["Size"])
            print(f"Found object {object_name}: RA: {result['RA']}, Dec: {result['Dec']}, Size: {result['Size']}")
    else:
        print(f"Failed to retrieve data for object: {object_name}")
        # remove object name from list
        object_list.pop()
        continue

# Build a string with the object names
tag = "+".join(object_list)

# Get the type of telescope
while True:
    telescope_type = input("Enter the type of telescope you are using ('EV1', 'EV2', 'Equinox1', 'Equinox2', 'Odyssey', 'Other'): ")
    if telescope_type.lower() in ["ev1", "ev2", "equinox1", "equinox2", "odyssey", "other"]:
        if telescope_type.lower() == "ev1":
            fov_width = 0.61*60*60 # based on https://astronomynow.com/2022/02/22/unistellar-evscope-equinox/
            fov_height = 0.46*60*60
        elif telescope_type.lower() == "ev2":
            fov_width = 45.6*60 # based on https://www.unistellar.com/expert/
            fov_height = 34.2*60
        elif telescope_type.lower() == "equinox1":
            fov_width = 0.61*60*60 # based on https://astronomynow.com/2022/02/22/unistellar-evscope-equinox/
            fov_height = 0.46*60*60
        elif telescope_type.lower() == "equinox2":
            fov_width = 45.6*60 # based on https://www.unistellar.com/expert/
            fov_height = 34.2*60
        elif telescope_type.lower() == "odyssey":
            fov_width = 45*60 # based on https://www.unistellar.com/discovery/
            fov_height = 33.6*60
        else:
            fov_width = float(input("Enter the field of view width of the telescope (in arcminutes): "))
            fov_height = float(input("Enter the field of view height of the telescope (in arcminutes): "))
            # convert width and height to arcseconds
            fov_width = fov_width * 60
            fov_height = fov_height * 60
        break
    else:
        print("Invalid telescope type. Assuming 'EV2' telescope.")
        telescope_type = "EV2"
        fov_width = 45.6*60 # based on https://www.unistellar.com/expert/
        fov_height = 34.2*60
        break

# Get the overlap percentage
overlap = float(input("Enter the overlap percentage between tiles (0-100) - default is 40%: ") or 40) / 100
# Get the size of the border as a percentage of the object's size
border_percentage = float(input("Enter the size of the border as a percentage of the mosaic's size (0-100) - default is 20%: ") or 20) / 100
desired_aspect_ratio = input("Enter the decimal value of the aspect ratio of the desired mosaic (width/height) - default is 16/9, enter min for smallest mosaic size: ")
if desired_aspect_ratio.lower() == "min":
    aspect_ratio = 0
else:
    try: 
        aspect_ratio = float(desired_aspect_ratio)
    except ValueError:
        print("Invalid aspect ratio. Using default value of 16/9.")
        aspect_ratio = 16/9

# Get the name of the file to save the links to
output_file = input("Enter the name of the file to save the links to (or just hit enter to skip): ")
if output_file:
    output_file = output_file.strip()
    # Specify the type of output to generate (default is html)
    output_type = input("Enter the type of output to generate (html, csv, text) - default is html: ")
    if output_type.lower() == "csv":
        output_type = "csv"
    elif output_type.lower() == "text":
        output_type = "text"
    else:
        # Warn the user if they enter something other than html, csv, or text
        print("Invalid output type. Using default value of html.")
        output_type = "html"
else:
    output_file = None
    output_type = "text"

# Get the user's location
latitude, longitude = get_location()
# Get the observing location from the latitude and longitude
observerLocation = EarthLocation.from_geodetic(lat=latitude*u.deg, lon=longitude*u.deg, height=0*u.m)

# Get the observing date and time from the user
while True:
    try:
        timestring = input("Enter the observing date and time (YYYY-MM-DD HH:MM) -- default is tonight at 9pm: ")
        if not timestring:
            timestring = datetime.now().replace(hour=21, minute=0, second=0, microsecond=0).strftime("%Y-%m-%d %H:%M")
        # Parse the date and time string
        observingDate = Time(timestring)
        break
    except ValueError:
        #print("Invalid date and time format. Please enter the observing date and time in the format YYYY-MM-DD HH:MM.")
        print("Invalid date and time format. Assuming tonight at 9pm.")
        observingDate = Time(datetime.now().replace(hour=21, minute=0, second=0, microsecond=0))
        break

## Get the expected minutes between observations
#while True:
#    try:
#        minutes_string = input("Enter the expected number of minutes between observations (default 10): ")
#        if not minutes_string:
#            minutes_string = "10"
#        minutesBetweenObservations = int(minutes_string)
#        break
#    except ValueError:
#        print("Invalid input. Please enter a positive integer.")

# Print the Alt/Az of the objects input by the user
# Calculate the Alt/Az of the objects input by the user
xobserverLocation = EarthLocation.from_geodetic(lat=latitude*u.deg, lon=-longitude*u.deg)
xaltaz = AltAz(obstime=observingDate, location=observerLocation)
# Loop through the list of objects and calculate their Alt/Az
for xobject_index, xobject_ra in enumerate(ra_list):
    # Get the Alt/Az of the object at the specified observing date and time
    xobject_dec = dec_list[xobject_index]
    # Calculate the Alt/Az of the object
    xobject_alt, xobject_az = convert_radec_to_altAz(xobject_ra, xobject_dec, xaltaz)
    alt_list.append(xobject_alt)
    az_list.append(xobject_az)
    print(f'Object {xobject_index+1}: RA: {xobject_ra}, Dec: {xobject_dec}, Alt: {xobject_alt}, Az: {xobject_az}')

# Get the common subplot for the tile rectangles
fig, ax = plt.subplots(1)
fig_altaz, ax_altaz = plt.subplots(1)

# Convert each coordinate in the ra_list / dec_list to alt / az
main_observerLocation = EarthLocation.from_geodetic(lat=latitude*u.deg, lon=-longitude*u.deg)
main_altaz = AltAz(obstime=observingDate, location=main_observerLocation)

main_alt_az_list = []
for i in range(len(object_list)):
    main_alt_az_list.append(convert_radec_to_altAz(ra_list[i], dec_list[i], main_altaz))

# Generate the mosaic using the above information
generate_mosaic(object_list, ra_list, dec_list, size_list, fov_width, fov_height, overlap, border_percentage, tag, latitude, longitude, observingDate, ax, ax_altaz, output_file, output_type)

# Using object_list, ra_list and dec_list, add the objects to the plot, with different shapes and colors for each object
# Create a scatter plot of the RA and Dec values
ax.scatter(ra_list, dec_list, marker='o', color='red', label='Objects')

# Get list of alt and az values from main_alt_az_list
main_alt_list = []
main_az_list = []
for i in range(len(object_list)):
    main_alt_list.append(main_alt_az_list[i][0])
    main_az_list.append(main_alt_az_list[i][1])
ax_altaz.scatter(main_az_list, main_alt_list, marker='o', color='red', label='Objects')
# Add labels for each object
for i, obj_name in enumerate(object_list):
    ax.annotate(obj_name, (ra_list[i], dec_list[i]), textcoords="offset points", xytext=(0, 10), ha='center')
    ax_altaz.annotate(obj_name, (main_az_list[i], main_alt_list[i]), textcoords="offset points", xytext=(0, 10), ha='center')
# Add a legend
ax.legend()
#ax_altaz.legend()
# Alt-Az legend needs to go outside the graph because the graph generally fills the plot area
# To do that, we need to set the position of the legend and reduce the size of the graph so the legend fits without being cut off
pos = ax_altaz.get_position()
#ax_altaz.set_position([pos.x0, pos.y0, pos.width * 0.8, pos.height * 0.8])
ax_altaz.set_position([pos.x0, pos.y0 + (pos.height * 0.2), pos.width, pos.height * 0.8])
lgd = ax_altaz.legend(bbox_to_anchor=(0., -0.13), loc='upper left', borderaxespad=0.5)


# Show the plot
plt.show()

exit(0)

