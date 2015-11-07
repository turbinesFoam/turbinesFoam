#!/usr/bin/env python

from subprocess import check_output

def ale():
    """Test actuatorLineElement."""
    output = check_output("./testActuatorLineElement")
