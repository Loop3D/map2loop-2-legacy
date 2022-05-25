from enum import IntEnum, Enum

class Datatype(IntEnum):
    GEOLOGY = 0
    STRUCTURE = 1
    FAULT = 2
    FOLD = 3
    MINERAL_DEPOSIT = 4
    DTM = 5
    SECTION = 6
    DRILLHOLE = 7
    DTB_GRID = 8
    COVER_MAP = 9
    METADATA = 10

class Datastate(IntEnum):
    UNNAMED = 0
    UNLOADED = 1
    LOADED = 2
    REPROJECTED = 3
    CLIPPED = 4
    CONVERTED = 5
    COMPLETE = 6
    MERGED = 7
    ERRORED = 9

class ErrorState(IntEnum):
    NONE = 0
    URLERROR = 1
    CONFIGERROR = 2

class VerboseLevel(IntEnum):
    NONE = 0
    TEXTONLY = 1
    ALL = 2


