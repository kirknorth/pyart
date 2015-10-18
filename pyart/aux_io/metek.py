"""
pyart.aux_io.metek
==================

Routines for reading METEK GmbH vertically pointing and scanning cloud radars,
e.g., the MIRA-35 cloud radar.

"""

import netCDF4
import numpy as np

from scipy import constants

from ..io.cfradial import _ncvar_to_dict
from ..config import FileMetadata
from ..core.radar import Radar


def read_mira35(filename, field_names=None, additional_metadata=None,
                file_field_names=False, exclude_fields=None, debug=False,
                verbose=False):
    """
    Read METEK MIRA-35 vertically pointing or scanning cloud radar data. Radar
    volumes with variable sweep modes are supported.

    Parameters
    ----------
    filename : str
        Name of NetCDF file to read data from.
    field_names : dict, optional
        Dictionary mapping field names in the file names to radar field names.
        Unlike other read functions, fields not in this dictionary or having a
        value of None are still included in the radar.fields dictionary, to
        exclude them use the `exclude_fields` parameter. Fields which are
        mapped by this dictionary will be renamed from key to value.
    additional_metadata : dict of dicts, optional
        This parameter is not used, it is included for uniformity.
    file_field_names : bool, optional
        True to force the use of the field names from the file in which
        case the `field_names` parameter is ignored. False will use to
        `field_names` parameter to rename fields.
    exclude_fields : list or None, optional
        List of fields to exclude from the radar object. This is applied
        after the `file_field_names` and `field_names` parameters.
    debug : bool, optional
        True to print debugging information, False to suppress.
    verbose : bool, optional
        True to print progress and identification information, False to
        suppress.

    Returns
    -------
    radar : Radar
        Radar containing available axes information, data fields, and metadata.

    """


    # create metadata retrieval object
    filemetadata = FileMetadata(
        'cfradial', field_names, additional_metadata, file_field_names,
        exclude_fields)

    # read the data
    ncobj = netCDF4.Dataset(filename)

    # parse variables and attributes
    ncatts = ncobj.ncattrs
    ncvars = ncobj.variables

    # 4.1 Global attributes -> move to metadata dictionary
    # hard code mobile platform and varying gates attributes
    metadata = dict([(k, getattr(ncobj, k)) for k in ncatts])
    metadata['platform_is_mobile'] = 'false'
    metadata['n_gates_vary'] = 'false'
    metadata['scan_id'] = 0

    # 4.2 Dimensions (do nothing)
    # TODO: address frequency dimension

    # 4.3 Global variables -> move to metadata dictionary
    # ignore time_* global variables, these are calculated from the time
    # variable when the file is written
    # hard code most of these variables
    metadata['volume_number'] = 0
    metadata['platform_type'] = 'fixed'
    metadata['instrument_type'] = 'radar'
    metadata['primary_axis'] = 'axis_z'

    # 4.4 coordinate variables -> create attribute dictionaries
    # hard code units of time variable
    time = _ncvar_to_dict(ncvars['time'])
    time['units'] = 'seconds since 1970-01-01 00:00:00 UTC'
    _range = _ncvar_to_dict(ncvars['range'])

    # 4.5 Ray dimension variables
    # check for constant gate spacing
    delta_range = np.diff(_range['data'])
    if np.allclose(delta_range.std(), 0.0, atol=0.1):
        _range['spacing_is_constant'] = 'true'
        _range['meters_between_gates'] = delta_range.mean()

    # 4.6 Location variables -> create attribute dictionaries
    # from the data files I have looked at so far, I suspect that the latitude,
    # longitude, and altitude data is recorded incorrectly
    latitude = filemetadata('latitude')
    latitude['data'] = float(ncobj.Latitude.split()[0])

    longitude = filemetadata('longitude')
    longitude['data'] = float(ncobj.Longitude.split()[0])

    altitude = filemetadata('altitude')
    altitude['data'] = float(ncobj.Altitude.split()[0])

    # 4.7 Sweep variables -> create attribute dictionaries
    # note that no sweep information is available in the original files
    sweep_number = _sweep_number(
        ncobj, filemetadata, debug=debug, verbose=verbose)

    sweep_mode = _sweep_mode(ncobj, filemetadata, debug=debug, verbose=verbose)

    sweep_start_idx, sweep_end_idx = _sweep_idx(
        ncobj, filemetadata, debug=debug, verbose=verbose)

    fixed_angle = _fixed_angle(
        ncobj, sweep_mode, sweep_start_idx, sweep_end_idx, filemetadata,
        debug=debug, verbose=verbose)

    # first sweep mode determines scan type
    scan_type = sweep_mode['data'][0]


    # 4.8 Sensor pointing variables -> create attribute dictionaries
    azimuth = _ncvar_to_dict(ncvars['azi'])
    elevation = _ncvar_to_dict(ncvars['elv'])

    # 4.9 Moving platform geo-reference variables -> not applicable

    # 4.10 Moments field data variables -> field attribute dictionary
    # all variables with dimensions of 'time', 'range' are fields
    fields = {}
    for var, data in ncvars.iteritems():
        if exclude_fields is not None and var in exclude_fields:
            continue
        if data.dimensions == ('time', 'range'):
            field = filemetadata.get_field_name(var)
            if field is None:
                field = var
            fields[field] = _ncvar_to_dict(data)

    # 4.5 instrument_parameters sub-convention -> instrument_parameters dict
    wavelength = _ncvar_to_dict(ncvars['lambda'])
    prf = _ncvar_to_dict(ncvars['prf'])
    nyquist = _ncvar_to_dict(ncvars['NyquistVelocity'])
    samples = _ncvar_to_dict(ncvars['nave'])


    # 4.6 radar_parameters sub-convention -> instrument_parameters dict
    instrument_parameters = {}

    # 4.7 lidar_parameters sub-convention -> not applicable

    # 4.8 radar_calibration sub-convention -> skip

    # close NetCDF object
    ncobj.close()

    return Radar(
        time, _range, fields, metadata, scan_type, latitude, longitude,
        altitude, sweep_number, sweep_mode, fixed_angle, sweep_start_idx,
        sweep_end_idx, azimuth, elevation, antenna_transition,
        instrument_parameters=instrument_parameters)


def _sweep_number(ncobj, filemetadata, debug=False, verbose=False):
    """
    Determine the number of sweeps in the radar volume.

    Parameters
    ----------
    ncobj : netCDF4.Dataset

    filemetadata : FileMetadata

    debug : bool, optional
        True to print debugging information, False to suppress.
    verbose : bool, optional
        True to print progress and identification information, False to
        suppress.

    Returns
    -------
    sweep_number : dict
        Sweep number dictionary compatible with Radar.

    """

    if verbose:
        print 'Processing sweep number'

    # parse sweep number metadata
    sweep_number = filemetadata('sweep_number')

    # parse azimuth and elevation angle data
    azi = ncobj.variables['azi'][:]
    elv = ncobj.variables['elv'][:]

    # compute fixed angles, indices, and counts of azimuth and elevation data
    azi_fixed, azi_idx, azi_counts = np.unique(
        azi, return_index=True, return_counts=True)
    elv_fixed, elv_idx, elv_counts = np.unique(
        elv, return_index=True, return_counts=True)

    if debug:
        nsweeps = sweep_number['data'].size
        print 'Number of sweeps in volume: {}'.format(nsweeps)

    return sweep_number


def _sweep_mode(ncobj, filemetadata, debug=False, verbose=False):
    """
    Detect the scan modes of the radar by investigating the relative changes
    between azimuth and elevation angle.

    The CF/Radial convention defines the possible sweep modes: 'sector',
    'coplane', 'rhi', 'vertical_pointing', 'idle', 'azimuth_surveillance',
    'elevation_surveillance', 'sunscan', 'pointing', 'manual_ppi',
    'manual_rhi'.

    Parameters
    ----------
    ncobj : netCDF4.Dataset
        Radar Dataset containing valid azimuth and elevation data.
    debug : bool, optional
        True to print debugging information, False to suppress.
    verbose : bool, optional
        True to print progress and identification information, False to
        suppress.


    Returns
    -------
    sweep_mode : dict
        Sweep mode dictionary compatible with Radar.

    """

    # parse sweep mode metadata
    sweep_mode = filemetadata('sweep_mode')

    # parse azimuth and elevation angle data
    azi = ncobj.variables['azi'][:]
    elv = ncobj.variables['elv'][:]

    # pompute changes in azimuth and elevation angle
    delta_azi = np.diff(azi, n=1)
    delta_elv = np.diff(elv, n=1)



    return sweep_mode
