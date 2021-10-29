import os, sys, importlib

home_dir = os.path.expanduser("~")
proj_dir = os.path.join(home_dir, '1000Genomes')
sys.path.insert(0, proj_dir)

# __all__ = ["alleles",
#            "pool",
#            "pybcf"]
