"""
Simplified EarthAccess MODIS Data Wrapper for FireDPy
====================================================
"""

import earthaccess
import os
import requests
from typing import List, Optional, Tuple, Union
import logging

logger = logging.getLogger(__name__)

class MODISEarthAccess:
    """Simplified EarthAccess-based MODIS data access for FireDPy."""

    def __init__(self, username: Optional[str] = None, password: Optional[str] = None):
        self.username = username
        self.password = password
        self._authenticated = False
        self._setup_authentication()

    def _setup_authentication(self):
        """Setup earthaccess authentication."""
        try:
            if self.username and self.password:
                self.auth = earthaccess.login(username=self.username, password=self.password)
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

    def download_file(self, url: str, dest_path: str) -> bool:
        """Download a file using standard requests with earthaccess authentication."""

        if not self._authenticated:
            return False

        try:
            os.makedirs(os.path.dirname(dest_path), exist_ok=True)

            # Use requests with earthaccess authentication
            response = requests.get(url, auth=(self.username, self.password), stream=True)
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

# Global instance for backward compatibility
_modis_access = None

def setup_modis_earthaccess(username: str, password: str) -> MODISEarthAccess:
    """Setup global MODIS earthaccess instance."""
    global _modis_access
    try:
        _modis_access = MODISEarthAccess(username, password)
        return _modis_access
    except Exception as e:
        print(f"EarthAccess setup failed: {e}")
        return None

def get_modis_earthaccess() -> MODISEarthAccess:
    """Get the global MODIS earthaccess instance."""
    if _modis_access is None:
        raise RuntimeError("Call setup_modis_earthaccess() first")
    return _modis_access
