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
	z = np.sqrt(r**2 + 2.0 * r * Re * np.sin(theta_e) + Re**2)
	s = Re * np.arcsin(r * np.cos(theta_e) / (z + Re))
	x = s * np.sin(theta_a)
	y = s * np.cos(theta_a)

	return z, y, x

def map_to_grid(radar, grid_dimensions, grid_origin=None, fields=None,
				map_roi=True, weighting_function='Cressman',
				toa=17000.0, leafsize=10.0, roi_func='dist_beam',
				proj='lcc', datum='NAD83', ellps='GRS80', 
				fill_value=None):
	"""
	"""

	# Parse missing value
	if fill_value is None:
		fill_value = get_fillvalue()

	# Check weight function parameters
	if weighting_function.upper() not in ['CRESSMAN', 'BARNES']:
		raise ValueError('Unsupported weighting_function')

	# Get grid origin if not given
	if grid_origin is None:
		lat0 = radar.latitude['data']
		lon0 = radar.longitude['data']
		grid_origin = (lat0, lon0)

	# Get fields which should be mapped
	if fields is None:
		fields = radar.fields.keys()
	nfields = len(fields)

	# Calculate radar offset from the origin
	radar_alt = radar.altitude['data']
	radar_lat = radar.latitude['data']
	radar_lon = radar.longitude['data']
	pj = pyproj.Proj(proj=proj, lat_0=grid_origin[0], lon_0=grid_origin[1],
					 x_0=0.0, y_0=0.0, datum=datum, ellps=ellps)
	radar_x, radar_y = pj(radar_lon, radar_lat)
	offset = (radar_alt, radar_x, radar_y)

	# Calculate Cartesian locations of radar gates
	z_g, y_g, x_g = _radar_coords_to_cartesian(radar, debug=debug)

	# Parse Cartesian locations of analysis grid
	z_a, y_a, x_a = grid_dimensions

	# Create k-d tree object from analysis grid
	tree_a = spatial.cKDTree(zip(z_a, y_a, x_a), leafsize=leafsize)

	return grid