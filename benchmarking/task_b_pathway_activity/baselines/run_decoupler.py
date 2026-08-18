"""decoupleR baseline(s) for Task B: standard non-spatial pathway-activity scoring.

Serves as the "does spatial information actually help" control -- these methods score
each spot's pathway activity purely from its own expression vector, no spatial graph, no
generative/latent model at all. Uses the same pathways.gmt (and hence the same gene-pathway
mask) as GRAPHIST/VEGA for a fair, apples-to-apples comparison.

Runs 2 representative decoupleR methods: ULM (fast, linear-model-based, close to the
field's current default) and GSVA (the classic single-sample gene-set enrichment method).
AUCell is intentionally excluded from the default set: it hits a numba JIT compilation
error ("multiple values for argument 'stop'") in this environment's numba/numpy version
pairing, unrelated to our code or to AUCell's actual algorithm -- not worth debugging an
unrelated environment issue when ULM+GSVA already cover both major method families
(linear-model-based and rank/enrichment-based) decoupleR represents.
"""
import argparse
import os

import decoupler as dc
import pandas as pd


METHODS = {
    "ulm": dc.mt.ulm,
    "gsva": dc.mt.gsva,
    "aucell": dc.mt.aucell,  # excluded by default -- see module docstring
}
DEFAULT_METHODS = ["ulm", "gsva"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--methods", nargs="+", default=DEFAULT_METHODS, choices=list(METHODS.keys()))
    args = parser.parse_args()

    expr = pd.read_csv(os.path.join(args.data_dir, "st_expression.csv"), index_col=0)
    net = dc.pp.read_gmt(os.path.join(args.data_dir, "pathways.gmt"))

    os.makedirs(args.out_dir, exist_ok=True)
    for name in args.methods:
        fn = METHODS[name]
        print(f"Running decoupleR method: {name}...")
        scores, _ = fn(data=expr, net=net, tmin=3, verbose=False)
        out_path = os.path.join(args.out_dir, f"decoupler_{name}_predictions.csv")
        scores.to_csv(out_path)
        print(f"  wrote {scores.shape} to {out_path}")


if __name__ == "__main__":
    main()
