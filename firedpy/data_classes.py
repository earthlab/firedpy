import os
from typing import List, Union
from datetime import datetime
import geopandas as gpd
import shutil


PROJECT_DIR = os.path.dirname(os.path.dirname(__file__))


class Base:
    """
    """
    MODIS_CRS = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
    MODIS_SINUSOIDAL_PATH = os.path.join(PROJECT_DIR, "ref", "modis_grid.gpkg")
    CONUS_SHAPEFILE_PATH = os.path.join(PROJECT_DIR, 'ref', 'boundaries', 'conus.gpkg')
    DEFAULT_TILES = ["h08v04", "h09v04", "h10v04", "h11v04", "h12v04",
                     "h13v04", "h08v05", "h09v05", "h10v05", "h11v05",
                     "h12v05", "h08v06", "h09v06", "h10v06", "h11v06"]

    def __init__(self, out_dir: str, start_year: Union[None, int], end_year: Union[None, int], username: str,
                 password: str, tiles: List[str] = None):
        self._out_dir = out_dir
        self._start_year = start_year
        self._end_year = end_year
        self._username = username
        self._password = password
        self._date = datetime.today().strftime('%m-%d-%Y')
        self._cpus = os.cpu_count()
        self._modis_template_path = os.path.join(out_dir, 'rasters')
        self._modis_template_file_root = 'mosaic_template.tif'
        self._land_cover_path = os.path.join(out_dir, 'rasters', 'landcover')
        self._land_cover_file_root = 'lc_mosaic_'
        self._nc_path = os.path.join(out_dir, 'rasters', 'burn_area', 'netcdfs')
        self._hdf_path = os.path.join(out_dir, 'rasters', 'burn_area', 'hdfs')
        self._tiles = self.DEFAULT_TILES if tiles is None else tiles

        # Initialize output directory folders and files
        self._create_paths()
        self._get_shape_files()
        print("Project Folder: " + out_dir)

    def _create_paths(self):
        sub_folders = [
            os.path.join('rasters', 'burn_area'),
            os.path.join('rasters', 'burn_area', 'hdfs'),
            os.path.join('rasters', 'ecoregion'),
            os.path.join('rasters', 'landcover'),
            os.path.join('rasters', 'landcover', 'mosaics'),
            os.path.join('shapefiles', 'ecoregion'),
            os.path.join('tables')
        ]

        for f in [os.path.join(self._out_dir, sf) for sf in sub_folders]:
            os.makedirs(f, exist_ok=True)

    def _get_shape_files(self):
        """
        Just to grab some basic shapefiles needed for calculating statistics.
        """
        shape_dir = os.path.join(self._out_dir, "shapefiles")
        os.makedirs(shape_dir, exist_ok=True)

        files_to_copy = {
            'modis_sinusoidal_grid_world.shp': self.MODIS_SINUSOIDAL_PATH,
            'conus.shp': self.CONUS_SHAPEFILE_PATH
        }

        for dest_name, source_path in files_to_copy.items():
            dest_path = os.path.join(shape_dir, dest_name)
            if not os.path.exists(dest_path):
                shutil.copy(source_path, dest_path)

        run_conus_modis_path = os.path.join(shape_dir, 'conus_modis.shp')
        if not os.path.exists(run_conus_modis_path):
            print("Reprojecting state shapefile to MODIS Sinusoidal...")
            conus_path = os.path.join(shape_dir, 'conus.shp')
            conus = gpd.read_file(conus_path)
            modis_conus = conus.to_crs(self.MODIS_CRS)
            modis_conus.to_file(run_conus_modis_path)


class BurnData(Base):
    def __init__(self, out_dir: str, start_year: Union[None, int], end_year: Union[None, int], username: str,
                 password: str, tiles: List[str] = None):
        super().__init__(out_dir, start_year, end_year, username, password, tiles)

    def get_burns(self):
        """
        This will download the MODIS burn event data set tiles and create a
        singular mosaic to use as a template file for coordinate reference
        information and geometries.

        User manual:
            http://modis-fire.umd.edu/files/MODIS_C6_BA_User_Guide_1.2.pdf

        Update 02/2021 -> fuoco server transitioned to SFTP Dec 2020
            Update firedpy to use Paramiko SSHClient / SFTPClient
            Server-side changes are described in the user manual linked above

        SFTP:
            sftp://fire:burnt@fuoco.geog.umd.edu/gfed4/MCD64A1/C6/
            username: fire
            password: burnt

        """
        # Check into the UMD SFTP fuoco server using Paramiko
        ssh_client = paramiko.SSHClient()
        ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh_client.connect(hostname="fuoco.geog.umd.edu", username="fire", password="burnt")
        print("Connected to 'fuoco.geog.umd.edu' ...")

        # Open the connection to the SFTP
        sftp_client = ssh_client.open_sftp()

        # Download the available files and catch failed downloads
        base_sftp_folder = os.path.join('data', 'MODIS', 'C6', 'MCD64A1', 'HDF')

        def get_tile_range(start_year, end_year):
            if start_year is None or end_year is None:
                return []
            return ["MCD64A1.A" + str(yr) for yr in range(start_year, end_year + 1)]

        def download_files(hdfs, sftp_client, hdf_tile_folder, tile_range):
            for h in tqdm(hdfs):
                if not tile_range or any(name in h for name in tile_range):
                    try:
                        remote = os.path.join(sftp_folder, h)
                        os.chdir(hdf_tile_folder)
                        sftp_client.get(remote, h)
                    except Exception as e:
                        print(e)

        def check_and_redownload_missing_files(hdfs, sftp_client, hdf_tile_folder, tile_range, sftp_folder):
            missings = []

            for hdf in hdfs:
                if not tile_range or any(name in hdf for name in tile_range):
                    trgt = os.path.join(hdf_tile_folder, hdf)
                    if not os.path.exists(trgt) or not is_file_valid(trgt):
                        missings.append(os.path.join(sftp_folder, hdf))

            if not missings:
                return

            print(f"Missed Files: {missings}")
            print("trying again...")
            download_files(missings, sftp_client, hdf_tile_folder, tile_range)

        def is_file_valid(filepath):
            try:
                gdal.Open(filepath).GetSubDatasets()[0][0]
                return True
            except Exception:
                print("Bad file detected, removing to try again...")
                os.remove(filepath)
                return False

        tile_range = get_tile_range(self._start_year, self._end_year)

        for tile in self._tiles:
            nc_file = os.path.join(self._out_dir, 'rasters', 'burn_area', 'netcdfs', f"{tile}.nc")

            if os.path.exists(nc_file):
                continue

            sftp_folder = os.path.join(base_sftp_folder, tile)
            hdf_tile_folder = os.path.join(self._hdf_path, tile)
            os.makedirs(hdf_tile_folder, exist_ok=True)

            try:
                print(f"Downloading/Checking HDF files for: {tile}")
                sftp_client.chdir(sftp_folder)
                hdfs = [h for h in sftp_client.listdir() if ".hdf" in h]

                download_files(hdfs, sftp_client, hdf_tile_folder, tile_range)
                check_and_redownload_missing_files(hdfs, sftp_client, hdf_tile_folder, tile_range, sftp_folder)
            except Exception as e:
                print(f"No MCD64A1 Product for tile: {tile}, skipping...")

        # Close new SFTP connection
        ssh_client.close()
        sftp_client.close()
        print("Disconnected from 'fuoco.geog.umd.edu' ...")

        # Build the netcdfs here
        tile_files = {}
        for tid in tiles:
            files = glob(os.path.join(self._hdf_path, tid, "*hdf"))
            tile_files[tid] = files

        # Merge one year into a reference mosaic
        if not os.path.exists(self._modis_template_path):
            print("Creating reference mosaic ...")
            folders = glob(os.path.join(self._hdf_path, "*"))
            file_groups = [glob(os.path.join(f, "*hdf")) for f in folders]
            for f in file_groups:
                f.sort()
            files = [f[0] for f in file_groups]
            dss = [rasterio.open(f).subdatasets[0] for f in files]
            tiles = [rasterio.open(d) for d in dss]
            mosaic, transform = merge(tiles)
            crs = tiles[0].meta.copy()
            template_path = os.path.join(self._modis_template_path,
                                         self._modis_template_file_root)
            crs.update({"driver": "GTIFF",
                        "height": mosaic.shape[1],
                        "width": mosaic.shape[2],
                        "transform": transform})

            with rasterio.open(template_path, "w+", **crs) as dst:
                dst.write(mosaic)

        # Build one netcdf per tile
        for tid in tiles:
            files = tile_files[tid]
            if len(files) > 0:
                try:
                    self.buildNCs(files)
                except Exception as e:
                    file_name = os.path.join(self._nc_path, tid + ".nc")
                    print("Error on tile " + tid + ": " + str(e))
                    print("Removing " + file_name + " and moving on.")
                    os.remove(file_name)