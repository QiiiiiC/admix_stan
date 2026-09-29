"""Derive the *_poisson IBD-only / mixed models from their normal-likelihood originals.

    python make_poisson.py            # from methods/stan/

The normal IBD likelihood needs an IBD standard error, and any estimate of it
from few blocks or few segments can collapse far below the counting error.
The derived model is the original with three changes and nothing else:

  * data: `ibd_hat` / `ibd_se` replaced by the raw counts `ibd_count`
  * model: the IBD likelihood replaced by the Poisson count likelihood of the
    *_Nsmooth_poisson models (lambda = segments per pair per cM * cM * pairs)
  * final IBD epoch anchored to the root (`t_end = t_root + T_max`), as in the
    Nsmooth models, instead of the absolute `t_end = T_max`

plus component `lp_ibd`/`chi2_ibd` (and `lp_snp`/`chi2_snp`) generated quantities.
An original that already anchors the final epoch keeps its own; an original whose
generated quantities only report the normal likelihood (they read `ibd_hat`) has
them replaced, since those quantities no longer exist in the derived model.
"""
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCES = [("ibd_model_Nfixed", False), ("mixed_model_Nfixed", True),
           ("ibd_model_Nvarying", False), ("mixed_model_Nvarying", True),
           ("ibd_model_Nsmooth", False)]

POISSON = '''    // ---- IBD likelihood: Poisson on the raw segment counts, every bin ----
    // Same likelihood as the *_Nsmooth_poisson models: lambda = (expected
    // segments per pair per cM) * (genome cM) * (number of haplotype pairs), so
    // the variance is the expected count and no IBD standard error is needed.
    // (The normal likelihood of *_Nfixed needs an SE, and any estimate of it
    // from few blocks or few segments can collapse far below the counting error.)
    // Full lpmf via target+=, so the -log(k!) terms are kept.
    for (i in 1:n_leaves) {
        for (j in i:n_leaves) {
            real n_pairs = (i == j)
                ? n_samples[i] * (n_samples[i] - 1) / 2.0
                : n_samples[i] * 1.0 * n_samples[j];
            for (b in 1:n_bins) {
                target += poisson_lpmf(ibd_count[b, i, j] |
                                       fmax(ibd_number[b][i, j] * cm * n_pairs, 1e-12));
            }
        }
    }'''
GQ_IBD = '''
generated quantities {
    // Component log-likelihood and Pearson chi2 of the Poisson counts, over every
    // bin (empty bins are observations under a Poisson).
    real lp_ibd = 0;
    real chi2_ibd = 0;
    for (i in 1:n_leaves) {
        for (j in i:n_leaves) {
            real n_pairs = (i == j)
                ? n_samples[i] * (n_samples[i] - 1) / 2.0
                : n_samples[i] * 1.0 * n_samples[j];
            for (b in 1:n_bins) {
                real lam = fmax(ibd_number[b][i, j] * cm * n_pairs, 1e-12);
                lp_ibd += poisson_lpmf(ibd_count[b, i, j] | lam);
                chi2_ibd += square(ibd_count[b, i, j] - lam) / lam;
            }
        }
    }
'''
GQ_SNP = '''    real lp_snp = 0;
    real chi2_snp = 0;
    for (i in 1:n_leaves) {
        for (j in i:n_leaves) {
            lp_snp += normal_lpdf(w_hat[i, j] | W_centered[i, j], w_se[i, j]);
            chi2_snp += square((w_hat[i, j] - W_centered[i, j]) / w_se[i, j]);
        }
    }
'''
DECL = ("    array[n_bins] matrix<lower=0>[n_leaves, n_leaves] ibd_hat;\n"
        "    array[n_bins] matrix<lower=0>[n_leaves, n_leaves] ibd_se;\n")
COUNTS = ("    // Raw segment counts per bin per leaf pair: what the Poisson likelihood conditions on.\n"
          "    array[n_bins, n_leaves, n_leaves] int<lower=0> ibd_count;\n")
TAIL_OLD = "                t_end = T_max;\n"
TAIL_NEW = ("                // T_max is the tail length PAST the root (as in the Nsmooth models), not an\n"
            "                // absolute cap: `t_end = T_max` inverts the final epoch once the root is\n"
            "                // older than T_max.\n"
            "                t_end = cumulative_times[n_events] + T_max;\n")
SNP_COMMENTS = ["    // ---- SNP likelihood (composite, with learned overdispersion) ----\n"]
SNP_COMMENT_NEW = "    // ---- SNP likelihood (composite: added to the IBD term as independent) ----\n"


def loop_end(src, start):
    """End index of the brace-delimited statement that begins at `start`."""
    i = src.index("{", start) + 1
    depth = 1
    while depth:
        depth += {"{": 1, "}": -1}.get(src[i], 0)
        i += 1
    return i


ANCHORED = "t_end = cumulative_times[n_events] + T_max;"


def convert(src, mixed):
    if src.count(DECL) != 1:
        raise ValueError(f"expected exactly one {DECL.strip()!r}")
    src = src.replace(DECL, COUNTS)
    if ANCHORED not in src:
        if src.count(TAIL_OLD) != 1:
            raise ValueError(f"expected exactly one {TAIL_OLD.strip()!r}")
        src = src.replace(TAIL_OLD, TAIL_NEW)
    gq = src.find("\ngenerated quantities {")
    if gq >= 0:
        if "ibd_hat" not in src[gq:]:
            raise ValueError("original has generated quantities not tied to the normal likelihood; merge by hand")
        src = src[:gq + 1]
    model = src.index("\nmodel {")
    start = src.index("    for (i in 1:n_leaves) {", model)
    end = loop_end(src, start)
    if "ibd_hat" not in src[start:end] or "ibd_number" not in src[start:end]:
        raise ValueError("first leaf-pair loop of the model block is not the IBD likelihood")
    line = src.rfind("\n", model, start - 1)
    if src[line + 1:start].strip().startswith("// ---- IBD likelihood"):
        start = line + 1
    src = src[:start] + POISSON + src[end:]
    for old in SNP_COMMENTS:
        src = src.replace(old, SNP_COMMENT_NEW)
    if "ibd_hat" in src or "ibd_se" in src:
        raise ValueError("normal-likelihood IBD data left in the derived model")
    return src.rstrip() + "\n" + GQ_IBD + (GQ_SNP if mixed else "") + "}\n"


if __name__ == "__main__":
    for name, mixed in SOURCES:
        out = HERE / f"{name}_poisson.stan"
        out.write_text(convert((HERE / f"{name}.stan").read_text(), mixed))
        print(f"{name}.stan -> {out.name}")
