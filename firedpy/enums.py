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
    """MODIS Land Cover Types.

    LC_Type1: Annual International Geosphere-Biosphere Programme (IGBP)
    LC_Type2: Annual University of Maryland (UMD)
    LC_Type3: Annual Leaf Area Index (LAI)
    LC_Type4: Annual BIOME-Biogeochemical Cycles (BGC)
    LC_Type5: Annual Plant Functional Types (PFT)
    """
    NONE = 0
    IGBP = 1
    UMD = 2
    MODIS_LAI = 3
    MODIS_BGC = 4
    PFT = 5
