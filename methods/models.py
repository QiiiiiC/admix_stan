"""Stan model registry, and normalised copies whose lp__ is a true log density.

Graphs with different numbers of events or admixtures are different models, and
comparing their evidence (ELBO or importance-sampled logZ) needs every constant.
Whether `~` statements keep theirs depends on how the engine evaluates the
density (Pathfinder was measured to keep them: raw vs normalised ELBO differ by
< 0.5 nats on snp_Nfixed, 2026-09-28), but no engine supplies the truncation
constants of `times >= 1` or of half-Normal scales.  `normalized_model`
rewrites every `~` as an explicit `target += *_lpdf(...)` and adds the
truncation constants, which it derives from the model source, so lp__ is the
full log density under every engine.  Posteriors are unchanged; only additive
constants move.

Stan files live in methods/stan/:
  *_Nfixed         one shared effective_N, gamma(4, 4/15000) prior; normal IBD
                   likelihood on ibd_hat/ibd_se with the Poisson zero term for
                   empty bins (needs an IBD SE; unreliable when that SE comes
                   from few blocks)
  *_Nfixed_poisson the same with the Poisson count likelihood below
  *_Nvarying       independent Ne per branch: log-normal around a shared level
                   (mu_log, sigma_log), no smoothing between neighbouring branches
  *_poisson        (IBD-only / mixed) derived from the normal-likelihood file by
                   stan/make_poisson.py
  *_Nsmooth        per-branch Ne on a tree-structured log-normal random walk
  *_poisson        IBD likelihood Poisson on raw segment counts, every bin
  *_separate_ne    independent Ne trajectories for the IBD and SNP components
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

STAN = Path(__file__).resolve().parent / "stan"


@dataclass(frozen=True)
class ModelSpec:
    name: str
    file: str
    data: str             # "ibd", "snp" or "mixed"
    ne: str               # "fixed" (effective_N), "varying" (independent per branch) or "smooth" (random walk)
    likelihood: str       # IBD likelihood: "normal" or "poisson"; "none" for SNP-only
    trajectories: int = 1  # Nsmooth: 2 = separate IBD and SNP Ne

    @property
    def path(self):
        return STAN / self.file

    @property
    def data_names(self):
        """Variables the model's `data {}` block declares."""
        block = re.sub(r"//[^\n]*", "", _block(self.path.read_text(), "data"))
        return set(re.findall(r"\b(\w+)\s*;", block))


MODELS = {m.name: m for m in [
    ModelSpec("ibd_Nfixed", "ibd_model_Nfixed.stan", "ibd", "fixed", "normal"),
    ModelSpec("snp_Nfixed", "snp_model_Nfixed.stan", "snp", "fixed", "none"),
    ModelSpec("mixed_Nfixed", "mixed_model_Nfixed.stan", "mixed", "fixed", "normal"),
    ModelSpec("ibd_Nfixed_poisson", "ibd_model_Nfixed_poisson.stan", "ibd", "fixed", "poisson"),
    ModelSpec("mixed_Nfixed_poisson", "mixed_model_Nfixed_poisson.stan", "mixed", "fixed", "poisson"),
    ModelSpec("ibd_Nvarying_poisson", "ibd_model_Nvarying_poisson.stan", "ibd", "varying", "poisson"),
    ModelSpec("snp_Nvarying", "snp_model_Nvarying.stan", "snp", "varying", "none"),
    ModelSpec("mixed_Nvarying_poisson", "mixed_model_Nvarying_poisson.stan", "mixed", "varying", "poisson"),
    ModelSpec("ibd_Nsmooth", "ibd_model_Nsmooth.stan", "ibd", "smooth", "normal"),
    ModelSpec("ibd_Nsmooth_poisson", "ibd_model_Nsmooth_poisson.stan", "ibd", "smooth", "poisson"),
    ModelSpec("snp_Nsmooth", "snp_model_Nsmooth.stan", "snp", "smooth", "none"),
    ModelSpec("mixed_Nsmooth_poisson", "mixed_model_Nsmooth_poisson.stan", "mixed", "smooth", "poisson"),
    ModelSpec("mixed_Nsmooth_poisson_separate", "mixed_model_Nsmooth_poisson_separate_ne.stan",
              "mixed", "smooth", "poisson", trajectories=2),
]}

_SAMPLING = re.compile(
    r"(?m)^(\s*)([\w]+(?:\[[^\]\n]+\])?)\s*~\s*(normal|std_normal|exponential|beta|gamma)\(([^;\n]*)\);")


def _block(source, name):
    """Text of a top-level block (`parameters {...}`), brace-matched."""
    m = re.search(rf"(?m)^{name}\s*\{{", source)
    if not m:
        return ""
    depth, i = 1, m.end()
    while depth:
        depth += {"{": 1, "}": -1}.get(source[i], 0)
        i += 1
    return source[m.end():i - 1]


def truncation_constants(source):
    """Stan expression for the normalisers `~` cannot supply for bounded parameters.

    times: vector<lower=1> with exponential(r) prior -> P(T > 1) = exp(-r), so
    +r per event.  Scalars declared <lower=0> with a normal(0, s) prior are
    half-Normals -> +log 2 each.  Other bounded priors (beta on [0,1], gamma on
    [0, inf)) are already normalised on their support.
    """
    params = _block(source, "parameters")
    terms = []
    m = re.search(r"(?m)^\s*times\s*~\s*exponential\(([^;\n]*)\);", source)
    if m and re.search(r"vector<\s*lower\s*=\s*1\s*>\s*\[n_events\]\s*times\s*;", params):
        rate = float(eval(m.group(1), {"__builtins__": {}}))   # a numeric literal like 0.01 or 1.0/100
        terms.append(f"n_events * {rate!r}")
    half = [name for name in re.findall(r"(?m)^\s*real<\s*lower\s*=\s*0\s*>\s*(\w+)\s*;", params)
            if re.search(rf"(?m)^\s*{name}\s*~\s*normal\(\s*0(\.0*)?\s*,", source)]
    if half:
        terms.append(f"{len(half)} * log(2.0)")
    return " + ".join(terms), half


def normalized_source(spec):
    source = spec.path.read_text()
    extra, half = truncation_constants(source)

    def replace(m):
        indent, lhs, distribution, args = m.groups()
        return f"{indent}target += {distribution}_lpdf({lhs}" + (f" | {args}" if args.strip() else "") + ");"

    source, n = _SAMPLING.subn(replace, source)
    if re.search(r"(?m)^\s*[\w\[\], ]+~", _block(source, "model")):
        raise ValueError(f"{spec.file}: a sampling statement the normaliser does not know survived")
    if extra:
        source = source.replace(
            "\nmodel {", "\nmodel {\n    // Normalisers for truncated priors: times >= 1"
            + (f" and half-Normal {', '.join(half)}" if half else "") + ".\n"
            f"    target += {extra};", 1)
    return source, n


def normalized_model(spec, folder):
    """Compile the normalised copy of `spec` into folder/ (reused if unchanged)."""
    from cmdstanpy import CmdStanModel
    if isinstance(spec, str):
        spec = MODELS[spec]
    source, _ = normalized_source(spec)
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / spec.file
    if not path.exists() or path.read_text() != source:
        path.write_text(source)
    return CmdStanModel(stan_file=str(path))


def model_files():
    """Every Stan source, for manifest hashing."""
    return sorted(STAN.glob("*.stan"))


