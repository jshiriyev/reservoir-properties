import importlib

import os

import sys

module = importlib.import_module('phaseo._standings_correlation')

print(getattr(module,'StandingsCorrelation'))