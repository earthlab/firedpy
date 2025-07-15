import os
import re
import sys
import requests
import urllib.parse
from tqdm import tqdm
from typing import List, Dict
from http.cookiejar import CookieJar
import rasterio
from rasterio.merge import merge

# Assume LPDAAC and LandCoverType are defined elsewhere in your codebase
class LandCover(LPDAAC):
    def __init__(self, out_dir: str, n_cores: int = None, username: str = None, password: str = None):
        super().__init__(out_dir)
        self._lp_daac_url = 'https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/MCD12Q1.061/'
        self._parallel_cores = n_cores if n_cores is not None else os.cpu_count() - 1
        self._username = username
        self._password = password
        self._file_regex = r'MCD12Q1\.A\d{7}\.h\d{2}v\d{2}\.061\.\d{13}\.hdf'

    def _generate_local_hdf_path(self, year: str, remote_name: str) -> str:
        return os.path.join(self._land_cover_dir, year, remote_name)

    def _generate_local_hdf_dir(self, year: str) -> str:
        return os.path.join(self._land_cover_dir, year)

    def _query_cmr_granules(self, tiles: List[str]) -> Dict[str, List[str]]:
        cmr_url = "https://cmr.earthdata.nasa.gov/search/granules.json"
        short_name = "MCD12Q1"
        version = "061"

        tile_patterns = [f"h{tile[1:3]}v{tile[4:]}" for tile in tiles]  # e.g., 'h10v09'

        params = {
            "short_name": short_name,
            "version": version,
            "provider": "LPCLOUD",
            "page_size": 2000,
            "page_num": 1,
        }

        granules_by_year = {}
        while True:
            response = requests.get(cmr_url, params=params)
            response.raise_for_status()
            items = response.json()["feed"]["entry"]
            if not items:
                break
            for item in items:
                granule_id = item["title"]
                if any(tile in granule_id for tile in tile_patterns):
                    year = granule_id.split(".")[1][:4]
                    granules_by_year.setdefault(year, []).append(granule_id)
            params["page_num"] += 1

        return granules_by_year

    def _create_requests(self, granules: List[str]):
        download_requests = []
        for granule in granules:
            year = granule.split('.')[1][:4]
            remote_file = f"{granule}.hdf"
            local_file_path = self._generate_local_hdf_path(year, remote_file)
            os.makedirs(os.path.dirname(local_file_path), exist_ok=True)

            if os.path.exists(local_file_path):
                if not self._verify_hdf_file(local_file_path):
                    print('Removing corrupt file:', local_file_path)
                    os.remove(local_file_path)
                else:
                    continue

            url = f"{self._lp_daac_url}{granule}/{remote_file}"
            download_requests.append((url, local_file_path))

        return download_requests

    def _create_annual_mosaic(self, year: str, land_cover_type: LandCoverType = LandCoverType.IGBP):
        output_file = f"lc_mosaic_{land_cover_type.value}_{year}.tif"
        lc_files = [self._generate_local_hdf_path(year, f) for f in os.listdir(self._generate_local_hdf_dir(year))
                    if re.match(self._file_regex, os.path.basename(f)) is not None]

        datasets = []
        for lc_file_path in lc_files:
            with rasterio.open(lc_file_path) as lc_file:
                datasets.append([sd for sd in lc_file.subdatasets if str(land_cover_type.value) in sd.lower()][0])

        tiles = [rasterio.open(ds) for ds in datasets]
        mosaic, transform = merge(tiles)

        crs = tiles[0].meta.copy()
        crs.update({"driver": "GTiff", "height": mosaic.shape[1], "width": mosaic.shape[2], "transform": transform})

        with rasterio.open(os.path.join(self._mosaics_dir, output_file), "w+", **crs) as dst:
            dst.write(mosaic)

    def get_land_cover(self, tiles: List[str] = None, land_cover_type: LandCoverType = LandCoverType.IGBP):
        if tiles is None:
            tiles = self.get_all_available_tiles()

        print("Querying granules across all years...")
        granules_by_year = self._query_cmr_granules(tiles)

        for year, granules in granules_by_year.items():
            output_file = f"lc_mosaic_{land_cover_type.value}_{year}.tif"
            if os.path.exists(os.path.join(self._mosaics_dir, output_file)):
                continue

            print(f"Downloading data for {year}...")
            download_requests = self._create_requests(granules)
            self._download_files(download_requests)

            print(f"Mosaicking land cover tiles for {year}...")
            self._create_annual_mosaic(str(year), land_cover_type)

        print(f"Land cover data saved to {self._mosaics_dir}")
