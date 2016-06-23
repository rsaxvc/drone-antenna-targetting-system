#Ported from http://cosinekitty.com/compass.html by Don Cross
import math

def ParseAngle(angle, limit):
	if (math.isnan (angle) or (angle < -limit) or (angle > limit)):
		return None
	else:
		return angle

def ParseElevation (angle):
	if (math.isnan (angle)):
		return None
	else:
		return angle

def parseLocation (lat,lon,elv):
	lat = ParseAngle (lat, 90.0)
	if (lat != None):
		lon = ParseAngle (lon, 180.0)
		if (lon != None):
			elv = ParseElevation (elv)
			if (elv != None):
				return { 'lat':lat, 'lon':lon, 'elv':elv}
	return None

def EarthRadiusInMeters(latitudeRadians): # latitude is geodetic, i.e. that reported by GPS
	# http:#en.wikipedia.org/wiki/Earth_radius
	a = 6378137.0  # equatorial radius in meters
	b = 6356752.3  # polar radius in meters
	cos = math.cos (latitudeRadians)
	sin = math.sin (latitudeRadians)
	t1 = a * a * cos
	t2 = b * b * sin
	t3 = a * cos
	t4 = b * sin
	return math.sqrt ((t1*t1 + t2*t2) / (t3*t3 + t4*t4))

def GeocentricLatitude(lat):
	# Convert geodetic latitude 'lat' to a geocentric latitude 'clat'.
	# Geodetic latitude is the latitude as given by GPS.
	# Geocentric latitude is the angle measured from center of Earth between a point and the equator.
	# https:#en.wikipedia.org/wiki/Latitude#Geocentric_latitude
	e2 = 0.00669437999014
	return math.atan((1.0 - e2) * math.tan(lat))

def LocationToPoint (latDeg, lonDeg, elvM, oblate):
	# Convert (lat, lon, elv) to (x, y, z).
	lat = latDeg * math.pi / 180.0
	lon = lonDeg * math.pi / 180.0
	clat = lat
	radius = 6371009.0
	if( oblate ):
		radius = EarthRadiusInMeters(lat)
		clat = GeocentricLatitude(lat)

	cosLon = math.cos(lon)
	sinLon = math.sin(lon)
	cosLat = math.cos(clat)
	sinLat = math.sin(clat)
	x = radius * cosLon * cosLat
	y = radius * sinLon * cosLat
	z = radius * sinLat

	# We used geocentric latitude to calculate (x,y,z) on the Earth's ellipsoid.
	# Now we use geodetic latitude to calculate normal vector from the surface, to correct for elevation.
	cosGlat = math.cos(lat)
	sinGlat = math.sin(lat)

	nx = cosGlat * cosLon
	ny = cosGlat * sinLon
	nz = sinGlat

	x += elvM * nx
	y += elvM * ny
	z += elvM * nz

	return { 'x':x, 'y':y, 'z':z, 'radius':radius, 'nx':nx, 'ny':ny, 'nz':nz}

def Distance (ap, bp):
	dx = ap['x'] - bp['x']
	dy = ap['y'] - bp['y']
	dz = ap['z'] - bp['z']
	return math.sqrt (dx*dx + dy*dy + dz*dz)

def RotateGlobe (latB, lonB, elvB, latA, lonA, elvA, bradius, aradius, oblate):
	# Get modified coordinates of 'b' by rotating the globe so that 'a' is at lat=0, lon=0.
	brp = LocationToPoint(latA, lonB - lonA, elvA, oblate)

	# Rotate brp cartesian coordinates around the z-axis by lonA degrees,
	# then around the y-axis by latA degrees.
	# Though we are decreasing by latA degrees, as seen above the y-axis,
	# this is a positive (counterclockwise) rotation (if B's longitude is east of A's).
	# However, from this point of view the x-axis is pointing left.
	# So we will look the other way making the x-axis pointing right, the z-axis
	# pointing up, and the rotation treated as negative.

	alat = -latA * math.pi / 180.0
	if (oblate):
		alat = GeocentricLatitude(alat)
	acos = math.cos(alat)
	asin = math.sin(alat)

	bx = (brp['x'] * acos) - (brp['z'] * asin)
	by = brp['y']
	bz = (brp['x'] * asin) + (brp['z'] * acos)

	return { 'x':bx, 'y':by, 'z':bz, 'radius':bradius}

def NormalizeVectorDiff(b,a):
	# Calculate norm(b-a), where norm divides a vector by its length to produce a unit vector.
	dx = b['x'] - a['x']
	dy = b['y'] - a['y']
	dz = b['z'] - a['z']
	dist2 = dx*dx + dy*dy + dz*dz
	if (dist2 == 0):
		return None
	dist = math.sqrt(dist2)
	return { 'x':(dx/dist), 'y':(dy/dist), 'z':(dz/dist) }
	#return { 'x':(dx/dist), 'y':(dy/dist), 'z':(dz/dist), 'radius':1.0 }

def Calculate(latA, lonA, elvA, latB, lonB, elvB, oblate):
	ap = LocationToPoint(latA, lonA, elvA, oblate)
	bp = LocationToPoint(latB, lonB, elvB, oblate)

	azimuth = 0.0
	altitude = 0.0
	distance = Distance(ap,bp)

	# Let's use a trick to calculate azimuth:
	# Rotate the globe so that point A looks like latitude 0, longitude 0.
	# We keep the actual radii calculated based on the oblate geoid,
	# but use angles based on subtraction.
	# Point A will be at x=radius, y=0, z=0.
	# Vector difference B-A will have dz = N/S component, dy = E/W component.
	br = RotateGlobe (latB, lonB, elvB, latA, lonA, elvA, bp['radius'], ap['radius'], oblate)
	if (br['z']*br['z'] + br['y']*br['y'] > 1.0e-6):
		theta = math.atan2(br['z'], br['y']) * 180.0 / math.pi
		azimuth = 90.0 - theta
		if (azimuth < 0.0):
			azimuth += 360.0
		if (azimuth > 360.0):
			azimuth -= 360.0

	bma = NormalizeVectorDiff(bp, ap)
	if (bma != None) :
		# Calculate altitude, which is the angle above the horizon of B as seen from A.
		# Almost always, B will actually be below the horizon, so the altitude will be negative.
		# The dot product of bma and norm = cos(zenith_angle), and zenith_angle = (90 deg) - altitude.
		# So altitude = 90 - acos(dotprod).
		altitude = 90.0 - (180.0 / math.pi)*math.acos(bma['x']*ap['nx'] + bma['y']*ap['ny'] + bma['z']*ap['nz'])
	return {"distance":distance, "azimuth":azimuth, "altitude":altitude}

print Calculate(40,38,0,40,38.001,100,True)
