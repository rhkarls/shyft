"""
Tests GeoTsRepositoryCollection
"""
from __future__ import print_function
from __future__ import absolute_import
from os import path
# import random
import unittest
# import numpy as np

from shyft import api
from shyft import shyftdata_dir
from shyft.repository.geo_ts_repository_collection import GeoTsRepositoryCollection
from shyft.repository.geo_ts_repository_collection import GeoTsRepositoryCollectionError
from shyft.repository.netcdf.met_netcdf_data_repository import MetNetcdfDataRepository
from shyft.repository.netcdf.met_netcdf_data_repository import MetNetcdfDataRepositoryError
from shapely.geometry import box


class GeoTsRepositoryCollectionTestCase(unittest.TestCase):
    @property
    def arome_epsg_bbox(self):
        """A slice of test-data located in shyft-data repository/arome."""
        EPSG = 32632
        x0 = 436100.0  # lower left
        y0 = 6823000.0  # lower right
        nx = 74
        ny = 24
        dx = 1000.0
        dy = 1000.0
        return EPSG, ([x0, x0 + nx * dx, x0 + nx * dx, x0], [y0, y0, y0 + ny * dy, y0 + ny * dy]), box(x0, y0, x0 + dx * nx, y0 + dy * ny)

    def test_get_timeseries_collection(self):
        tc= api.YMDhms(2015, 8, 24, 6)
        n_hours = 30
        dt = api.deltahours(1)
        utc = api.Calendar()  # No offset gives Utc
        t0 = utc.time(tc)
        period = api.UtcPeriod(t0, t0 + api.deltahours(n_hours))
        date_str = "{}{:02}{:02}_{:02}".format(tc.year, tc.month, tc.day, tc.hour)

        epsg, bbox, bpoly = self.arome_epsg_bbox

        base_dir = path.join(shyftdata_dir, "repository", "arome_data_repository")
        f1 = "arome_metcoop_red_default2_5km_{}.nc".format(date_str)
        f2 = "arome_metcoop_red_test2_5km_{}.nc".format(date_str)

        ar1 = MetNetcdfDataRepository(epsg, base_dir, filename=f1, allow_subset=True)
        ar2 = MetNetcdfDataRepository(epsg, base_dir, filename=f2, elevation_file=f1, allow_subset=True)

        geo_ts_repository = GeoTsRepositoryCollection([ar1, ar2])
        sources = geo_ts_repository.get_timeseries(("temperature", "radiation"),
                                                   period, geo_location_criteria=bpoly)

        with self.assertRaises(GeoTsRepositoryCollectionError) as context:
            GeoTsRepositoryCollection([ar1, ar2], reduce_type="foo")

        geo_ts_repository = GeoTsRepositoryCollection([ar1, ar2], reduce_type="add")
        with self.assertRaises(GeoTsRepositoryCollectionError) as context:
            sources = geo_ts_repository.get_timeseries(("temperature", "radiation"),
                                                       period, geo_location_criteria=bpoly)

    def test_get_forecast_collection(self):
        n_hours = 30
        dt = api.deltahours(1)
        utc = api.Calendar()  # No offset gives Utc
        tc = api.YMDhms(2015, 8, 24, 6)
        t0 = utc.time(tc)
        period = api.UtcPeriod(t0, t0 + api.deltahours(n_hours))
        date_str = "{}{:02}{:02}_{:02}".format(tc.year, tc.month, tc.day, tc.hour)

        epsg, bbox, bpoly = self.arome_epsg_bbox

        base_dir = path.join(shyftdata_dir, "repository", "arome_data_repository")
        f1_elev = "arome_metcoop_red_default2_5km_{}.nc".format(date_str)
        f1 = "arome_metcoop_red_default2_5km_(\d{4})(\d{2})(\d{2})[T_](\d{2})Z?.nc$"
        #f2 = "arome_metcoop_red_test2_5km_{}.nc".format(date_str)
        f2 = "arome_metcoop_red_test2_5km_(\d{4})(\d{2})(\d{2})[T_](\d{2})Z?.nc$"

        ar1 = MetNetcdfDataRepository(epsg, base_dir, filename=f1, allow_subset=True)
        ar2 = MetNetcdfDataRepository(epsg, base_dir, filename=f2, elevation_file=f1_elev, allow_subset=True)

        geo_ts_repository = GeoTsRepositoryCollection([ar1, ar2])
        source_names = ("temperature", "radiation")
        sources = geo_ts_repository.get_forecast(source_names, period, t0,
                                                 geo_location_criteria=bpoly)
        self.assertTrue(all([x in source_names for x in sources]))

        geo_ts_repository = GeoTsRepositoryCollection([ar1, ar2], reduce_type="add")
        with self.assertRaises(GeoTsRepositoryCollectionError) as context:
            sources = geo_ts_repository.get_forecast(("temperature", "radiation"),
                                                     period, t0, geo_location_criteria=bpoly)

    def test_get_ensemble_forecast_collection(self):
        EPSG = 32633
        upper_left_x = 436100.0
        upper_left_y = 7417800.0
        nx = 74
        ny = 94
        dx = 1000.0
        dy = 1000.0
        t0 = api.YMDhms(2015, 7, 26, 0)
        n_hours = 30
        utc = api.Calendar()  # No offset gives Utc
        period = api.UtcPeriod(utc.time(t0), utc.time(t0) + api.deltahours(n_hours))
        t_c = utc.time(t0) + api.deltahours(1)

        base_dir = path.join(shyftdata_dir, "netcdf", "arome")
        # pattern = "fc*.nc"
        pattern = "fc_(\d{4})(\d{2})(\d{2})[T_](\d{2})Z?.nc$"
        bpoly = box(upper_left_x, upper_left_y - ny * dy, upper_left_x + nx * dx, upper_left_y)
        try:
            ar1 = MetNetcdfDataRepository(EPSG, base_dir, filename=pattern)
            ar2 = MetNetcdfDataRepository(EPSG, base_dir, filename=pattern)
            repos = GeoTsRepositoryCollection([ar1, ar2])
            data_names = ("temperature", "wind_speed", "relative_humidity")
            ensemble = repos.get_forecast_ensemble(data_names, period, t_c, None)
            self.assertTrue(isinstance(ensemble, list))
            self.assertEqual(len(ensemble), 10)
            with self.assertRaises(GeoTsRepositoryCollectionError) as context:
                repos = GeoTsRepositoryCollection([ar1, ar2], reduce_type="add")
                repos.get_forecast_ensemble(data_names, period, t_c, geo_location_criteria=bpoly)
            self.assertEqual("Only replace is supported yet", context.exception.args[0])
        except MetNetcdfDataRepositoryError as adre:
            self.skipTest("(test inconclusive- missing arome-data {0})".format(adre))


if __name__ == '__main__':
    unittest.main()
