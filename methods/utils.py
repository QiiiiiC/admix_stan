"""Small shared utilities: silencing, atomic JSON, hashing, the run-manifest guard."""
from __future__ import annotations

import contextlib
import hashlib
import io
import json
import os
from pathlib import Path


def quiet(fn, *args, **kwargs):
    """Call fn with stdout discarded (msprime/cmdstanpy chatter)."""
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*args, **kwargs)


def save_json(path, value):
    """Atomic write.  The temp name carries the PID: parallel workers write
    identical provenance/manifest files, and a shared temp path let one worker's
    rename delete another's file mid-flight (FileNotFoundError on replace)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(f"{path.suffix}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    tmp.replace(path)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def file_hashes(paths, root):
    """sha256 of each file, keyed by its path relative to root."""
    return {str(Path(p).resolve().relative_to(Path(root).resolve())):
            hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in paths}


def open_run(folder, manifest):
    """Create or reopen a run folder.  The first call saves the manifest (config,
    code hashes, candidate set, block selections); later calls must match it
    exactly, so cached results from different code or config never mix."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    saved = folder / "manifest.json"
    if saved.exists():
        if digest(json.loads(saved.read_text())) != digest(manifest):
            raise ValueError("Configuration or source code changed. Use a new --out directory; "
                             "cached results must not be mixed.")
    else:
        save_json(saved, manifest)
    return manifest
