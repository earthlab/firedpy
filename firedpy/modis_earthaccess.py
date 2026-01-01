"""Simplified EarthAccess MODIS Data Wrapper for Firedpy."""
import earthaccess
import os
import requests

import logging

logger = logging.getLogger(__name__)


# Global instance for backward compatibility
MODIS_ACCESS = None


class MODISEarthAccess:
    """Simplified EarthAccess-based MODIS data access for Firedpy."""

    def __init__(self, username=None, password=None):
        """Initialize MODISEarthAccess object.

        Paramters
        ---------
        username : str
            Earth Access Username. Defaults to user prompt, Optional.
        """
        self.username = username
        self.password = password
        self._authenticated = False
        self._setup_authentication()

    def _setup_authentication(self):
        """Setup earthaccess authentication."""
        try:
            if self.username and self.password:
                self.auth = earthaccess.login(self.username, self.password)
            else:
                self.auth = earthaccess.login()

            if self.auth:
                self._authenticated = True
                logger.info("✅ EarthAccess authentication successful")
            else:
                raise RuntimeError("Authentication failed")

        except Exception as e:
            logger.error(f"EarthAccess authentication failed: {e}")
            raise

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
        if not self._authenticated:
            return False

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


def setup_modis_earthaccess(username, password):
    """Setup global MODIS earthaccess instance.

    Parameters
    ----------
    username : str
        Earthaccess account username.
    password : str
        Earthaccess account password.

    Returns
    -------
    firedpy.modis_earthaccess.MODISEarthAccess : A Firedpy MODIS Earth Access
        object.
    """

    global MODIS_ACCESS

    try:
        MODIS_ACCESS = MODISEarthAccess(username, password)
        return MODIS_ACCESS
    except Exception as e:
        print(f"EarthAccess setup failed: {e}")
        return None


def get_modis_earthaccess():
    """Get the global MODIS earthaccess instance.

    Returns
    -------
    firedpy.modis_earthaccess.MODISEarthAccess : A Firedpy MODIS Earth Access
        object.
    """
    if MODIS_ACCESS is None:
        raise RuntimeError("Call setup_modis_earthaccess() first")
    return MODIS_ACCESS
