"""Score Task B pathway-activity predictions against known ground truth.

Metrics:
- per-pathway correlation (Pearson/Spearman) between true and inferred activity
- recall@k: does each spot's top-k inferred pathways overlap its top-k truly-active ones
- DE-pathway F1: run a group A/B test on the inferred activity per pathway (Welch t-test),
  compare the predicted-DE pathway set against the ground-truth injected DE set -- the
  metric that maps most directly to GRAPHIST's actual use case (finding pathway biomarkers)
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr, ttest_ind


def per_pathway_correlation(true_activity: pd.DataFrame, pred: pd.DataFrame) -> dict:
    common = [p for p in true_activity.columns if p in pred.columns]
    pearsons, spearmans = [], []
    for p in common:
        t = true_activity[p].values
        pr = pred.loc[true_activity.index, p].values
        if np.std(pr) < 1e-10:
            continue
        pearsons.append(pearsonr(t, pr)[0])
        spearmans.append(spearmanr(t, pr)[0])
    return {
        "n_pathways_matched": len(common),
        "mean_pearson": float(np.nanmean(pearsons)) if pearsons else np.nan,
        "median_pearson": float(np.nanmedian(pearsons)) if pearsons else np.nan,
        "mean_spearman": float(np.nanmean(spearmans)) if spearmans else np.nan,
    }


def recall_at_k(true_activity: pd.DataFrame, pred: pd.DataFrame, k: int = 5) -> float:
    common = [p for p in true_activity.columns if p in pred.columns]
    true_sub = true_activity[common]
    pred_sub = pred.loc[true_activity.index, common]
    recalls = []
    for spot in true_activity.index:
        true_top = set(true_sub.loc[spot].nlargest(k).index)
        pred_top = set(pred_sub.loc[spot].nlargest(k).index)
        recalls.append(len(true_top & pred_top) / k)
    return float(np.mean(recalls))


def benjamini_hochberg(pvals: dict, alpha: float) -> set:
    """Standard BH step-up procedure. With n=hundreds of spots per group, raw p<alpha
    alone flags nearly everything as significant (huge power to detect trivial spatial
    autocorrelation differences between contiguous halves, not just the injected DE
    effect) -- multiple-testing correction is not optional here.
    """
    names = list(pvals.keys())
    p = np.array([pvals[n] for n in names])
    m = len(p)
    order = np.argsort(p)
    ranked = p[order]
    thresh = (np.arange(1, m + 1) / m) * alpha
    below = ranked <= thresh
    if not below.any():
        return set()
    max_rank = np.max(np.where(below)[0])
    return {names[i] for i in order[: max_rank + 1]}


def de_pathway_f1(
    pred: pd.DataFrame, groups: pd.Series, true_de: list, alpha: float = 0.05, min_cohens_d: float = 0.3
) -> dict:
    """A DE call requires BOTH BH-significance AND a minimum effect size (Cohen's d).
    With hundreds of spots per group, p-value alone flags near-everything as significant
    once model estimation noise leaks a whisper of group signal into a null pathway --
    exactly the reasoning behind GRAPHIST's own differential-activity test combining a
    probability threshold with an effect-size delta, not p-value alone.
    """
    a_idx = groups[groups == "A"].index
    b_idx = groups[groups == "B"].index
    pvals, cohens_d = {}, {}
    for p in pred.columns:
        a = pred.loc[a_idx, p].values
        b = pred.loc[b_idx, p].values
        pvals[p] = ttest_ind(a, b, equal_var=False).pvalue
        pooled_sd = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
        cohens_d[p] = (a.mean() - b.mean()) / (pooled_sd + 1e-12)
    sig = benjamini_hochberg(pvals, alpha)
    predicted_de = {p for p in sig if abs(cohens_d[p]) > min_cohens_d}
    true_de_set = set(true_de)
    all_pathways = set(pred.columns)
    tp = len(predicted_de & true_de_set)
    fp = len(predicted_de - true_de_set)
    fn = len(true_de_set - predicted_de)
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
    return {"de_precision": precision, "de_recall": recall, "de_f1": f1,
            "de_n_predicted": len(predicted_de), "de_n_true": len(true_de_set)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--sim-name", required=True)
    parser.add_argument("--out-csv", default=None)
    parser.add_argument("--k", type=int, default=5)
    args = parser.parse_args()

    true_activity = pd.read_csv(os.path.join(args.data_dir, "true_activity.csv"), index_col=0)
    groups = pd.read_csv(os.path.join(args.data_dir, "st_groups.csv"), index_col=0)["group"]
    true_de = pd.read_csv(os.path.join(args.data_dir, "true_de_pathways.csv"))["pathway"].tolist()

    rows = []
    for pred_path in sorted(glob.glob(os.path.join(args.data_dir, "*_predictions.csv"))):
        method = os.path.basename(pred_path).replace("_predictions.csv", "")
        pred = pd.read_csv(pred_path, index_col=0)
        metrics = {"sim": args.sim_name, "method": method}
        metrics.update(per_pathway_correlation(true_activity, pred))
        metrics[f"recall_at_{args.k}"] = recall_at_k(true_activity, pred, k=args.k)
        metrics.update(de_pathway_f1(pred, groups, true_de))
        rows.append(metrics)

    results = pd.DataFrame(rows).set_index(["sim", "method"])
    pd.set_option("display.width", 200)
    print(results.round(3))

    if args.out_csv:
        os.makedirs(os.path.dirname(args.out_csv), exist_ok=True)
        header = not os.path.exists(args.out_csv)
        results.to_csv(args.out_csv, mode="a", header=header)


if __name__ == "__main__":
    main()
