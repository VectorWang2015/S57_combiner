import math
from typing import Tuple
from geographiclib.geodesic import Geodesic


# in km
EARTH_RADIUS = 6371
BASE_POINT_NAME = "Haitang"
BASE_POINT_LON = 119.37116703323744
BASE_POINT_LAT = 34.75485890660782

def haversin(phi):
    """phi in radius"""
    return math.pow(math.sin(phi/2), 2)


def reverse_haversin(hs):
    """hs for haversin value"""
    sin_phi_divide_2 = math.sqrt(hs)
    phi_divide_2 = math.asin(sin_phi_divide_2)
    return phi_divide_2 * 2


def haversin_distance(
        coord_1: Tuple[float, float],
        coord_2: Tuple[float, float],
        ) -> float:
    """coords in angles"""
    lon_1, lat_1 = coord_1
    lon_2, lat_2 = coord_2
    # coords to radius
    lon_1 = lon_1 / 180 * math.pi
    lon_2 = lon_2 / 180 * math.pi
    lat_1 = lat_1 / 180 * math.pi
    lat_2 = lat_2 / 180 * math.pi
    # apply haversin equation
    arc_radius_haversin = haversin(lat_2 - lat_1) + math.cos(lat_1) * math.cos(lat_2) * haversin(lon_1 - lon_2)
    arc_radius = reverse_haversin(arc_radius_haversin)
    arc_distance = arc_radius * EARTH_RADIUS
    return arc_distance


def latlon_to_enu(ref_lon, ref_lat, lon, lat):
    """
    enu: east-north-up
    """
    geod = Geodesic.WGS84
    result = geod.Inverse(ref_lat, ref_lon, lat, lon)
    distance = result['s12']    # in meters
    azimuth = math.radians(result['azi1'])

    x = distance * math.sin(azimuth)
    y = distance * math.cos(azimuth)
    return x, y


def enu_to_latlon(ref_lon, ref_lat, x, y):
    """
    reversed latlon_to_enu
    """
    azimuth = math.degrees(math.atan2(x, y)) 
    distance = math.sqrt(x**2 + y**2)  # in meters

    geod = Geodesic.WGS84
    result = geod.Direct(ref_lat, ref_lon, azimuth, distance)

    return result['lon2'], result['lat2']