"""
pyart.map.

"""

import numpy as np

from scipy import spatial
from mpl_toolkits.basemap import pyproj

from ..config import get_fillvalue, get_field_name
from ..core.grid import Grid

def _radar_coords_to_cartesian(radar, debug=False):
	"""
	"""

	# Effective radius of Earth
	Re = 6371.0 * 4.0 / 3.0 * 1000.0

	# Get radar scan geometry and convert scan angles to radians
	r = radar.range['data']
	theta_e = np.radians(radar.elevation['data'])
	theta_a = np.radians(radar.azimuth['data'])

	# Create radar gate mesh
	theta_a = np.meshgrid(r, theta_a)[1]
	r, theta_e = np.meshgrid(r, theta_e)
	r = r.flatten()
	theta_a = theta_a.flatten()
	theta_e = theta_e.flatten()

	# Compute vertical height (z), arc length (s), eastward distance (x),
	# and northward distance (y)
	z = np.sqrt(r**2 + 2.0 * r * Re * np.sin(theta_e) + Re**2) - Re
	s = Re * np.arcsin(r * np.cos(theta_e) / (z + Re))
	x = s * np.sin(theta_a)
	y = s * np.cos(theta_a)

	return z, y, x

def map_to_grid(radar, grid_coords, grid_origin=None, fields=None,
				weighting_function='Cressman', toa=17000.0, leafsize=10.0, 
				k=1, eps=0.0, roi_func='constant', constant_roi=1000.0, 
				cutoff_radius=5000.0, min_radius=250.0, x_factor=0.01,
				y_factor=0.01, z_factor=0.01, map_roi=False, map_dist=True, 
				proj='lcc', datum='NAD83', ellps='GRS80', 
				fill_value=None, debug=False):
	"""
	Map one or more radars to a Cartesian analysis grid.

	Parameters
	----------
	radar : Radar
		A radar object to be mapped to the Cartesian analysis grid.
	grid_coords : tuple
		The (z, y, x) coordinates of the grid in meters. These can describe
		either a uniform or non-uniform grid.
	k : int
		The number of nearest neighbors to return. When 'weighting_function'
		is 'Barnes' or 'Cressman', k should be larger than 1.
	eps : float
		Return approximate nearest neighbors. The k-th returned value is
		guaranteed to be no further than (1 + eps) times the distance to
		the real k-th nearest neighbor.
	cutoff_radius : float
		The largest radius in meters to search for points. This is only valid 
		when 'weighting_function' is 'Barnes', and should be large enough to
		capture all points where the Barnes weight is nonzero.
	"""

	# Parse missing value
	if fill_value is None:
		fill_value = get_fillvalue()

	# Check weight function parameters
	if weighting_function not in ['Cressman', 'Barnes', 'nearest']:
		raise ValueError('Unsupported weighting_function')
	if weighting_function in ['Cressman', 'Barnes'] and k == 1:
		raise ValueError('For Cressman or Barnes weights, k > 1')

	# Get grid origin if not given
	if grid_origin is None:
		lat0 = radar.latitude['data'][0]
		lon0 = radar.longitude['data'][0]
		grid_origin = (lat0, lon0)

	# Get fields which should be mapped
	# If no fields are given, then all fields will be mapped
	if fields is None:
		fields = radar.fields.keys()
	else:
		if not (set(fields) & set(radar.fields.keys())):
			raise ValueError('One or more specified fields do not exist')


	# Calculate radar offset from the origin
	radar_alt = radar.altitude['data'][0]
	radar_lat = radar.latitude['data'][0]
	radar_lon = radar.longitude['data'][0]
	pj = pyproj.Proj(proj=proj, lat_0=grid_origin[0], lon_0=grid_origin[1],
					 x_0=0.0, y_0=0.0, datum=datum, ellps=ellps)
	radar_x, radar_y = pj(radar_lon, radar_lat)
	offset = (0.0, radar_y, radar_x)

	# Calculate Cartesian locations of radar gates relative to grid origin
	z_g, y_g, x_g = _radar_coords_to_cartesian(radar, debug=debug)
	z_g = z_g + offset[0]
	y_g = y_g + offset[1]
	x_g = x_g + offset[2]

	# Parse Cartesian locations of analysis grid
	z_a, y_a, x_a = grid_coords
	nz = len(z_a)
	ny = len(y_a)
	nx = len(x_a)
	if debug:
		print 'Grid shape is nz = %i, ny = %i, nx = %i' % (nz, ny, nx)

	# Create analysis grid mesh
	z_a, y_a, x_a = np.meshgrid(z_a, y_a, x_a, indexing='ij')
	z_a = z_a.flatten()
	y_a = y_a.flatten()
	x_a = x_a.flatten()

	if debug:
		print 'Grid array has shape %s' % (z_a.shape,)

	# Compute the radius of influence for each analysis point, if necessary
	if not hasattr(roi_func, '__call__'):
		if roi_func == 'constant':
			roi = constant_roi
			max_roi = constant_roi

		elif roi_func == 'dist':
			roi = default_roi_func_dist(
				    x_a, y_a, z_a, offset, x_factor=x_factor,
				    y_factor=y_factor, z_factor=z_factor,
				    min_radius=min_radius)
			max_roi = roi.max()

		else:
			raise ValueError('Unsupported roi_func')

		if debug:
			print 'Minimum ROI is %.2f m' % np.min(roi)
			print 'Maximum ROI is %.2f m' % np.min(roi)
			print 'ROI array has shape %s' % (np.shape(roi),)

	# Create k-d tree object for radar gate locations
	# Depending on the number of radar gates this can be resource intensive
	# but should nonetheless take on the order of seconds to create
	tree_g = spatial.cKDTree(zip(z_g, y_g, x_g), leafsize=leafsize)

	if debug:
		print 'tree.m = %i, tree.n = %i' % (tree_g.m, tree_g.n)

	# Query k-d tree
	if weighting_function == 'nearest':
		dist, ind = tree_g.query(zip(z_a, y_a, x_a), k=1, p=2.0, eps=eps,
			                     distance_upper_bound=np.inf)
	elif weighting_function == 'Barnes':
		dist, ind = tree_g.query(zip(z_a, y_a, x_a), k=k, p=2.0, eps=eps,
			                     distance_upper_bound=cutoff_radius)
	else:
		dist, ind = tree_g.query(zip(z_a, y_a, x_a), k=k, p=2.0, eps=eps,
			                     distance_upper_bound=max_roi)

	if debug:
		print 'Distance array has shape %s' % (dist.shape,)
		print 'Minimum distance is %.2f m' % dist.min()
		print 'Maximum distance is %.2f m' % dist.max()
		print 'Index array has shape %s' % (ind.shape,)
		print 'Minimum index is %i' % ind.min()
		print 'Maximum index is %i' % ind.max()

	# Missing neighbors are indicated with an index set to tree.n
	# This condition will not be met for the nearest neighbor scheme, but
	# it can be met for the Cressman and Barnes schemes
	# We can safely set the index of missing neighbors to 0 since
	# its weight can also be set to 0 later and thus not affect the
	# interpolation (averaging)
	bad_index = ind == tree_g.n
	ind[bad_index] = 0

	# Interpolate the radar data onto the analysis grid and populate the
	# mapped fields dictionary
	map_fields = {}
	for field in fields:
		radar_data = radar.fields[field]['data'].flatten()

		if weighting_function == 'nearest':
			fq = radar_data[ind]
			
		elif weighting_function == 'Barnes':
			wq = np.exp(-dist**2 / kappa)
			fq = np.ma.average(radar_data[ind], weights=wq, axis=1)

		else:
			roi_stack = np.repeat(roi, k).reshape(roi.size, k)
			wq = (roi_stack**2 - dist**2) / (roi_stack**2 + dist**2)
			not_close = dist > roi_stack
			wq[not_close] = 0.0
			fq = np.ma.average(radar_data[ind], weights=wq, axis=1)

		map_fields[field] = {'data': fq.reshape(nz, ny, nx)}

		# Populate mapped field metadata
		[map_fields[field].update({meta: value}) for meta, value in 
		 radar.fields[field].iteritems() if meta != 'data']


	# Map the nearest neighbor distances, if necessary
	if map_dist:
		field = 'nearest_neighbor_distance'
		map_fields[field] = {
			'units': 'meters',
			'long_name': 'Distance to closest radar gate'}
		if weighting_function == 'nearest':
			map_fields[field]['data'] = dist.reshape(nz, ny, nx)
		else:
			map_fields[field]['data'] = dist.min(axis=1).reshape(nz, ny, nx)

	# Map the radius of influence, if necessary
	if map_roi and 'roi' in locals():
		field = 'radius_of_influence'
		map_fields[field] = {
			'data': roi.reshape(nz, ny, nx),
			'units': 'meters',
			'long_name': 'Radius of influence used for mapping'}

	# Populate Cartesian axes dictionaries
	x_disp = {
		'data': grid_coords[2],
		'units': 'meters',
		'axis': 'X',
		'long_name': 'x-coordinate in Cartesian system'}

	y_disp = {
		'data': grid_coords[1],
		'units': 'meters',
		'axis': 'Y',
		'long_name': 'y-coordinate in Cartesian system'}

	z_disp = {
		'data': grid_coords[0],
		'units': 'meters',
		'axis': 'Z',
		'positive': 'up',
		'long_name': 'z-coordinate in Cartesian system'}

	# Populate grid origin dictionaries
	altitude = {
		'data': 0.0,
		'units': 'meters',
		'standard_name': 'altitude',
		'long_name': 'Altitude at grid origin above mean sea level'}

	latitude = {
		'data': grid_origin[0],
		'units': 'degrees_N',
		'standard_name': 'latitude',
		'valid_min': -90.0,
		'valid_max': 90.0,
		'long_name': 'Latitude at grid origin'}

	longitude = {
		'data': grid_origin[1],
		'units': 'degrees_E',
		'standard_name': 'longitude',
		'valid_min': -180.0,
		'valid_max': 180.0,
		'long_name': 'Longitude at grid origin'}


	# Populate time dictionaries
	time = {
		'data': radar.time['data'].min(),
		'units': radar.time['units'],
		'calendar': radar.time['calendar'],
		'standard_name': radar.time['standard_name'],
		'long_name': 'Time in seconds since volume start'}

	time_start = {
		'data': radar.time['data'].min(),
		'units': radar.time['units'],
		'calendar': radar.time['calendar'],
		'standard_name': radar.time['standard_name'],
		'long_name': 'Time in seconds since volume start'}

	time_end = {
		'data': radar.time['data'].max(),
		'units': radar.time['units'],
		'calendar': radar.time['calendar'],
		'standard_name': radar.time['standard_name'],
		'long_name': 'Time in seconds since volume end'}

	# Create axes
	axes = {
		'time': time,
		'time_start': time_start,
		'time_end': time_end,
		'x_disp': x_disp,
		'y_disp': y_disp,
		'z_disp': z_disp,
		'altitude': altitude,
		'latitude': latitude,
		'longitude': longitude}

	# Create metadata
	metadata = {
		'reference': '',
		'Conventions': '',
		'site': '',
		'facility_id': '',
		'project': '',
		'state': ''}

	return Grid(map_fields, axes, metadata)


def default_roi_func_dist(x_a, y_a, z_a, offset, x_factor=1.0, y_factor=1.0, 
	                      z_factor=1.0, min_radius=250.0):
	"""
	"""

	# Apply the offset to the analysis grid such that it is now
	# radar-centric 
	z_a = z_a - offset[0]
	y_a = y_a - offset[1]
	x_a = x_a - offset[2]

	roi = np.sqrt(x_factor * x_a**2 + y_factor * y_a**2 + 
		          z_factor * z_a**2) + min_radius

	return roi