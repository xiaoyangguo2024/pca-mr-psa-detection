"""Shared paths and helpers. Scripts are run as `python3 -I code/<script>.py` from the project root.

In isolated mode the user site-packages directory is not on sys.path, so it is added here explicitly;
nothing is ever imported from the data directories.
"""
import sys, os, site
_us = site.getusersitepackages()
if _us not in sys.path:
    sys.path.append(_us)
import json, math
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw")
REF = os.path.join(ROOT, "data", "ref")
PROC = os.path.join(ROOT, "data", "proc")
WORK = os.path.join(ROOT, "work")
RES = os.path.join(ROOT, "results")
TAB = os.path.join(RES, "tables")
FIG = os.path.join(RES, "figures")
ENV_BIN = os.path.join(ROOT, "env", "bin")
for _d in (PROC, WORK, TAB, FIG):
    os.makedirs(_d, exist_ok=True)


def save_json(obj, path):
    with open(path, "w") as f:
        json.dump(obj, f, indent=1, sort_keys=True, default=float)


def load_json(path):
    with open(path) as f:
        return json.load(f)


COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def is_palindromic(a1, a2):
    return COMP.get(a1) == a2


# display names of exposures (figures and tables)
NICE = {"age of sexual initiation": "Age at first sexual intercourse", "education attainment": "Educational attainment",
        "Class. Alphaproteobacteria": "Class Alphaproteobacteria", "Genus. Adlercreutzia": "Genus Adlercreutzia",
        "LDL": "LDL cholesterol", "LTL": "Leukocyte telomere length", "SLE": "Systemic lupus erythematosus"}


def nice(x):
    return NICE.get(x, x if x[:1].isupper() else x[:1].upper() + x[1:])
