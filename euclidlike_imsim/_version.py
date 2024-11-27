# Read version number from root directory
with open("../VERSION") as version_file:
    __version__ = version_file.read().strip()
    
__version_info__ = tuple(map(lambda x:int(x.split('-')[0]), __version__.split('.')))[:3]