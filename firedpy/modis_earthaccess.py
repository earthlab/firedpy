"""Simplified EarthAccess MODIS Data Wrapper for Firedpy."""
import earthaccess
import os
import requests

import logging

logger = logging.getLogger(__name__)


class MODISEarthAccess:
    """Simplified EarthAccess-based MODIS data access for Firedpy."""

    def __init__(self):
        """Initialize MODISEarthAccess object."""
        self._earthaccess = self._authenticate()
        self.username = self._earthaccess.username
        self.password = self._earthaccess.password

    def _authenticate(self):
        # This will use a config file or a prompt if one isn't available
        auth = earthaccess.login(strategy="all")
        if not auth:
            raise RuntimeError("EarthAccess authentication failed")
        return auth

    def download_file(self, url, dest_path):
        """Download file with earthaccess authentication.

        Parameters
        ----------
        url : str
            Target Earthaccess file URL.
        dest_path : str
            Target local file path.

        Returns
        -------
        bool : True if download successful, False if request not authenticated
               or if down load fails.
        """
        try:
            os.makedirs(os.path.dirname(dest_path), exist_ok=True)

            # Use requests with earthaccess authentication
            response = requests.get(
                url,
                auth=(self.username, self.password),
                stream=True
            )
            response.raise_for_status()

            with open(dest_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)

            logger.info(f"✅ Downloaded: {os.path.basename(dest_path)}")
            return True

        except Exception as e:
            logger.error(f"Download failed for {url}: {e}")
            return False
