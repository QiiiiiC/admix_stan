"""Compile the Stan models of this folder (they #include stan/graph.stanfunctions)."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
STAN = HERE / "stan"


def compile_model(name):
    from cmdstanpy import CmdStanModel
    return CmdStanModel(stan_file=str(STAN / f"{name}.stan"), stanc_options={"include-paths": [str(STAN)]})


def at(model, data, params, folder=None):
    """Transformed parameters and generated quantities at fixed parameter values."""
    fit = model.sample(data=data, inits=params, fixed_param=True, chains=1, iter_sampling=1, seed=1,
                       show_progress=False, show_console=False, output_dir=folder, sig_figs=17)
    return {k: v[0] for k, v in fit.stan_variables().items()}
