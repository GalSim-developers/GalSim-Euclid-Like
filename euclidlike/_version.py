# Read version number from root directory
current_dir = os.path.dirname(os.path.abspath(__file__))
version_path = os.path.join(current_dir, "..", "VERSION")

with open(version_path, "r") as version_file:
    __version__ = version_file.read().strip()
    
__version_info__ = tuple(map(lambda x:int(x.split('-')[0]), __version__.split('.')))[:3]