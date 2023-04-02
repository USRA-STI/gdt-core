__version__ = '2.0.0'

import os
from pathlib import Path

suite_path = Path(__file__).parent.parent

if 'GDT_BASE' in os.environ:
    base_path = Path(os.environ['GDT_BASE'])
else:
    base_path = Path.home().joinpath('.gammaray_data_tools', __version__)

cache_path = base_path.joinpath('cache')
# test_data = {instrument: test_data_path.joinpath(instrument) \
#              for instrument in os.listdir(test_data_path)}

if 'GDT_DATA' in os.environ:
    data_path = Path(os.environ['GDT_DATA'])
else:
    data_path = base_path.joinpath('test_data')
