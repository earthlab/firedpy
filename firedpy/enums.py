import enum


class ShapeType(enum.Enum):
    SHP = "shp"
    GPKG = "gpkg"
    BOTH = "both"
    NONE = "none"


class EcoRegionType(enum.Enum):
    """Ecoregion Data Sources/types

    world : World Terrestrial Ecoregions (World Wildlife Fund)
    na : North American ecoregions (Omernick, 1987)
    """
    NA = "na"
    WORLD = "world"
    NONE = None


class TileChoice(enum.Enum):
    A = "a"
    B = "b"
    C = "c"
    D = "d"
    E = "e"


class LandCoverType(enum.IntEnum):
    """MODIS Land Cover Types.

    LC_Type1 : International Geosphere-Biosphere Programme (IGBP) scheme
    LC_Type2 : University of Maryland (UMD) scheme
    LC_Type3 : MODIS-derived Leaf Area Index (LAI/fPAR) scheme
    LC_Type4 : MODIS-derived Net Primary Production (NPP) scheme
    LC_Type5 : Plant Functional Type (PFT) scheme
    """
    NONE = 0
    IGBP = 1
    UMD = 2
    MODIS_LAI = 3
    MODIS_BGC = 4
    PFT = 5
