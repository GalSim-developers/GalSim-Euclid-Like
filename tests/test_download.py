## Testing code was adapted from: https://github.com/GalSim-developers/GalSim/blob/f6499baeadf433a34f98f5c918b2580e6d28bf7b/tests/test_download.py

import sys
import logging
import builtins
import shutil
from unittest import mock
import galsim
import euclidlike
import scripts.download_psf as download_psf  # Not imported automatically
from galsim.utilities import  timer
from numpy.testing import assert_raises

# This file tests the euclidlike_download_psf executable.

@timer
def test_args():
    """Test the argument parsing
    """
    # -h prints description and exits
    print('This should print the description...')  # Can only check by inspection, not assert
    with assert_raises(SystemExit):
        args = download_psf.parse_args(['-h'])

    # Check defaults:
    args = download_psf.parse_args([])
    assert args.verbosity == 2
    assert args.force is False
    assert args.quiet is False
    assert args.unpack is False
    assert args.save is False
    assert args.dir is None
    assert args.nolink is False

    # Check setting each of those to be something else
    args = download_psf.parse_args(['-v', '3'])
    assert args.verbosity == 3
    args = download_psf.parse_args(['-v', '0'])
    assert args.verbosity == 0

    args = download_psf.parse_args(['-f'])
    assert args.force is True

    args = download_psf.parse_args(['-q'])
    assert args.quiet is True

    args = download_psf.parse_args(['-u'])
    assert args.unpack is True

    args = download_psf.parse_args(['--save'])
    assert args.save is True


    args = download_psf.parse_args(['--nolink'])
    assert args.nolink is True

    # Some invalid parameters
    # To avoid ugly text output during nosetests runs, redirect stderr to stdout for a moment.
    sys_stderr = sys.stderr
    sys.stderr = sys.stdout
    with assert_raises(SystemExit):
        download_psf.parse_args(['-v', '-1'])
    with assert_raises(SystemExit):
        download_psf.parse_args(['-v', '4'])
    sys.stderr = sys_stderr

# global for the bleh and delay functions
count = 0

@timer
def test_query():
    """Test the query_yes_no function

    Need to mock the input function for this
    """
    from scripts.download_psf import query_yes_no

    def bleh():
        global count
        count += 1
        return 'y' if count % 5 == 0 else 'bleh'

    def delay_y():
        global count
        count += 1
        return 'y' if count % 5 == 0 else ''

    def delay_n():
        global count
        count += 1
        return 'n' if count % 5 == 0 else ''

    with mock.patch('scripts.download_psf.get_input', return_value='y'):
        assert query_yes_no('', 'yes') == 'yes'
    with mock.patch('scripts.download_psf.get_input', return_value='yes'):
        assert query_yes_no('', 'yes') == 'yes'
    with mock.patch('scripts.download_psf.get_input', return_value='n'):
        assert query_yes_no('', 'yes') == 'no'
    with mock.patch('scripts.download_psf.get_input', return_value='no'):
        assert query_yes_no('', 'yes') == 'no'
    with mock.patch('scripts.download_psf.get_input', return_value=''):
        assert query_yes_no('', 'yes') == 'yes'
    with mock.patch('scripts.download_psf.get_input', bleh):
        assert query_yes_no('', 'yes') == 'yes'

    with mock.patch('scripts.download_psf.get_input', return_value='y'):
        assert query_yes_no('', 'no') == 'yes'
    with mock.patch('scripts.download_psf.get_input', return_value='yes'):
        assert query_yes_no('', 'no') == 'yes'
    with mock.patch('scripts.download_psf.get_input', return_value='n'):
        assert query_yes_no('', 'no') == 'no'
    with mock.patch('scripts.download_psf.get_input', return_value='no'):
        assert query_yes_no('', 'no') == 'no'
    with mock.patch('scripts.download_psf.get_input', return_value=''):
        assert query_yes_no('', 'no') == 'no'
    with mock.patch('scripts.download_psf.get_input', bleh):
        assert query_yes_no('', 'yes') == 'yes'

    with mock.patch('scripts.download_psf.get_input', return_value='y'):
        assert query_yes_no('', None) == 'yes'
    with mock.patch('scripts.download_psf.get_input', return_value='yes'):
        assert query_yes_no('', None) == 'yes'
    with mock.patch('scripts.download_psf.get_input', return_value='n'):
        assert query_yes_no('', None) == 'no'
    with mock.patch('scripts.download_psf.get_input', return_value='no'):
        assert query_yes_no('', None) == 'no'
    with mock.patch('scripts.download_psf.get_input', delay_n):
        assert query_yes_no('', None) == 'no'
    with mock.patch('scripts.download_psf.get_input', delay_y):
        assert query_yes_no('', None) == 'yes'
    with mock.patch('scripts.download_psf.get_input', bleh):
        assert query_yes_no('', None) == 'yes'

    with assert_raises(ValueError):
        query_yes_no('', 'invalid')

# Need to call these before each time make_logger is repeated.  Else duplicate handles.
def remove_handler():
    logger = logging.getLogger('euclidlike')
    for handler in logger.handlers:
        logger.removeHandler(handler)

@timer
def test_names():
    """Test the get_names function
    """
    from scripts.download_psf import get_names

    args =download_psf.parse_args([])
    remove_handler()
    logger = download_psf.make_logger(args)
    url, target, target_dir, link_dir, unpack_dir, do_link = get_names(args, logger)
    assert url == 'https://cosmo.vera.psc.edu/Euclid-like-PSFs/monopsfs_euclidlike.zip'


class fake_urlopen:
    err = None

    # We don't want to actually check things on the internet.  So this class fakes
    # up the real urlopen function.
    def __init__(self, url):
        self.n = 10
        pass
    def info(self):
        return {
            "Server": "nginx/1.16.1",
            "Content-Type": "application/octet-stream",
            "Content-Length": "728",
            "Connection": "close",
            "Content-MD5": "e05cfe60c037c645d61ac70545cc2a99",
            "Content-Security-Policy": "default-src 'none';",
            "X-Content-Type-Options": "nosniff",
            "X-Download-Options": "noopen",
            "X-Permitted-Cross-Domain-Policies": "none",
            "X-Frame-Options": "sameorigin",
            "X-XSS-Protection": "1; mode=block",
            "Content-Disposition": "attachment; filename= monopsfs_5_4.fits.gz",
            "ETag": "\"md5:e05cfe60c037c645d61ac70545cc2a99\"",
            "Last-Modified": "Sun, 31 May 2024 02:19:18 GMT",
            "Date": "Thu, 11 Jun 2024 16:06:07 GMT",
            "Accept-Ranges": "none",
            "X-RateLimit-Limit": "60",
            "X-RateLimit-Remaining": "59",
            "X-RateLimit-Reset": "1591891628",
            "Retry-After": "60",
            "Strict-Transport-Security": "max-age=0",
            "Referrer-Policy": "strict-origin-when-cross-origin",
            "Set-Cookie": "session=3765c5a1d5211e53_5ee2566f",
            "X-Session-ID": "3765c5a1d5211e53_5ee2566f",
            "X-Request-ID": "2b720f14bdd71a29031a5cb415b391f8"
        }
    def read(self, block_sz):
        if self.err:
            raise OSError(self.err)
        if self.n:
            self.n -= 1
            return b'x' * 80  # ignore block_sz
        else:
            return b''


if __name__ == "__main__":
    testfns = [v for k, v in vars().items() if k[:5] == 'test_' and callable(v)]
    for testfn in testfns:
        testfn()