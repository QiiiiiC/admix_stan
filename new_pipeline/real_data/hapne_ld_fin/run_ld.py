from configparser import ConfigParser
import argparse, logging
from hapne.ld import compute_ld, compute_ccld
from hapne import hapne_ld
logging.basicConfig(level=logging.INFO)
ap=argparse.ArgumentParser(); ap.add_argument("--config_file", required=True)
a=ap.parse_args(); cfg=ConfigParser(); cfg.read(a.config_file)
print("stage 2: compute_ld", flush=True);  compute_ld(cfg)
print("stage 2b: compute_ccld", flush=True); compute_ccld(cfg)
print("stage 3: hapne_ld", flush=True);   hapne_ld(cfg)
