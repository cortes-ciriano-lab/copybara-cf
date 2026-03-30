"""
Testing module for COPYBARA
Created: 30/03/2026
Python 3.9.6
Carolin Sauer
"""
#!/usr/bin/env python3

import subprocess

from pathlib import Path

ROOTDIR = Path(__file__).parent.parent

def test_install():
    """ test COPYBARA is installed successfully """
    cmd = f"python3 {ROOTDIR}/copybara/copybara.py --help"
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    _, err = p.communicate()
    if p.returncode != 0:
        raise RuntimeError(f"FAILED: {cmd}\n{err}")
