import enum


class ShapeType(enum.Enum):
    SHP = "shp"
    GPKG = "gpkg"
    BOTH = "both"
    NONE = "none"


class EcoRegionType(enum.Enum):
    NA = "na"
    WORLD = "world"


class TileChoice(enum.Enum):
    A = "a"
    B = "b"
    C = "c"
    D = "d"
    E = "e"


class LandCoverType(enum.IntEnum):
    NONE = 0
    IGBP = 1
    UMD = 2
    MODIS_LAI = 3
    MODIS_BGC = 4
    PFT = 5
