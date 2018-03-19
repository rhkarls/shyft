from shyft import api
from netCDF4 import Dataset
from .time_conversion import convert_netcdf_time
from shyft.repository.interfaces import GeoTsRepository, ForecastSelectionCriteria
from shyft.repository.netcdf.concat_data_repository import ConcatDataRepository
from shyft.repository.netcdf.met_netcdf_data_repository import MetNetcdfDataRepository


class WXRepositoryError(Exception):
    pass

class WXRepository(GeoTsRepository):

    def __init__(self, epsg, file_name, flattened=False):
        """
        Construct the netCDF4 dataset reader for concatenated gridded forecasts and initialize data retrieval.

        Parameters
        ----------
        epsg: string
            Unique coordinate system id for result coordinates. Currently "32632" and "32633" are supported.
        filename: string
            Path to netcdf file containing concatenated forecasts
        flattened: bool
            Flags whether grid_points are flattened
        """
        if flattened:
            self.wx_repo = ConcatDataRepository(epsg, file_name)
        elif not flattened:
            self.wx_repo = MetNetcdfDataRepository(epsg, None, file_name)
            with Dataset(file_name) as dataset:
                time = dataset.variables.get("time", None)
                time = convert_netcdf_time(time.units, time)
                self.wx_repo.time = time

    def get_timeseries_ensemble(self, input_source_types, utc_period, geo_location_criteria=None):
        """
        Get ensemble of shyft source vectors of time series covering utc_period
        for input_source_types.

        Time series are constructed by concatenating values from forecasts
        according to fc_periodicity whose lead period
        (nb_lead_intervals_to_drop, nb_lead_intervals_to_drop + fc_len_to_concat)
        intersect the utc_period. See _get_time_structure_from_dataset for details
        on fc_len_to_concat.


        Parameters
        ----------
        see interfaces.GeoTsRepository

        Returns
        -------
        see interfaces.GeoTsRepository
        """
        wx_repo = self.wx_repo
        diff_years = int((utc_period.start - wx_repo.time[0])//(365 * 24 * 3600))
        utc_start_shifted = utc_period.start - (diff_years) * 365 * 24 * 3600
        utc_end_shifted = utc_period.end - (diff_years) * 365 * 24 * 3600
        utc_period_shifted = api.UtcPeriod(utc_start_shifted, utc_end_shifted)
        return wx_repo.get_timeseries_ensemble(input_source_types, utc_period_shifted, geo_location_criteria)

