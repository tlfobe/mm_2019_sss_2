"""
mm_2019_sss_2
A short description of the project.
"""

# Add imports here
from .system import *
import mm_2019_sss_2.energy
import mm_2019_sss_2.monte_carlo

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
