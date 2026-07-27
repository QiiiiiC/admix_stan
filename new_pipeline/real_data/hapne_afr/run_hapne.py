from configparser import ConfigParser
import argparse, logging
from hapne.ibd import build_hist
from hapne import hapne_ibd

logging.basicConfig(level=logging.INFO)
ap = argparse.ArgumentParser()
ap.add_argument("--config_file", required=True)
a = ap.parse_args()
cfg = ConfigParser(); cfg.read(a.config_file)
build_hist(cfg)
hapne_ibd(cfg)
