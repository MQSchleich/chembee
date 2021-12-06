from chembee.graphics import polar_plot
from chembee.processing import load_data
import os
import sys
import pytest

from pathlib import Path


@pytest.fixture(scope="module")
def script_loc(request):
    '''Return the directory of the currently running test script'''

    return Path(request.fspath).parent 
def test_polar_plot(script_loc):
    
    frame =  load_data(script_loc.joinpath("data/Biodeg.sdf"))
    polar_plot(frame, save_path=script_loc.joinpath("plots/"))