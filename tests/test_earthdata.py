"""
Testing access to NASA's earth data
"""

from http.cookiejar import CookieJar
import urllib.request

def test_earthdata_credentials(username: str, password: str) -> None:
    # Earthdata Login
    # test url for correct user/password
    root_url = 'https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-public/MCD12Q1.061/'
    url = root_url + 'MCD12Q1.A2019001.h10v09.061.2022169160720/BROWSE.MCD12Q1.A2019001.h10v09.061.2022169160720.1.jpg'
    # the below url does not work as of July 2025
    # url = "https://e4ftl01.cr.usgs.gov/MOTA/MCD12Q1.061/2019.01.01/BROWSE.MCD12Q1.A2019001.h10v09.061.2022169160720.1.jpg"

    password_manager = urllib.request.HTTPPasswordMgrWithDefaultRealm()
    password_manager.add_password(None, "https://urs.earthdata.nasa.gov", username, password)
    # Create a cookie jar for storing cookies. This is used to store and return
    # the session cookie given to use by the data server (otherwise it will just
    # keep sending us back to Earthdata Login to authenticate).  Ideally, we
    # should use a file based cookie jar to preserve cookies between runs. This
    # will make it much more efficient.
    cookie_jar = CookieJar()
    # Install all the handlers.
    opener = urllib.request.build_opener(
        urllib.request.HTTPBasicAuthHandler(password_manager),
        # urllib.request.HTTPHandler(debuglevel=1),    # Uncomment these two lines to see
        # urllib.request.HTTPSHandler(debuglevel=1),   # details of the requests/responses
        urllib.request.HTTPCookieProcessor(cookie_jar))
    urllib.request.install_opener(opener)

    request = urllib.request.Request(url)
    urllib.request.urlopen(request)

# test the function
un = 'maco4303'
pw = '5Stringbanjo1!'
test_earthdata_credentials(un, pw)