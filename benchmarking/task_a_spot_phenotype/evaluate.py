"""Score Task A predictions (Graphist(+/-)/Background spot calls) against ground truth.

For each method's predictions CSV (cell_id, predicted in {-1,0,1}, score = raw regression
coefficient), computes precision/recall/F1 for the positive-group and negative-group calls
separately (score used for PR-AUC, since it's a continuous ranking signal), plus a combined
"any phenotype association" metric (predicted != 0 vs. ground_truth != 0) -- these are the
metrics SpaPheno and SpaLinker themselves report, chosen so results are directly comparable.
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, f1_score, precision_score, recall_score


def score_predictions(ground_truth: pd.Series, pred_df: pd.DataFrame) -> dict:
    used_cutoff = pred_df["used_cutoff"].iloc[0] if "used_cutoff" in pred_df.columns else np.nan
    pred_df = pred_df.set_index("cell_id").reindex(ground_truth.index).fillna(0)
    predicted = pred_df["predicted"].astype(int)
    score = pred_df["score"].astype(float) if "score" in pred_df.columns else predicted.astype(float)

    metrics = {"used_cutoff": used_cutoff}
    for group, label in [("positive", 1), ("negative", -1)]:
        y_true = (ground_truth == label).astype(int)
        y_pred = (predicted == label).astype(int)
        rank_score = score if label == 1 else -score
        metrics[f"{group}_precision"] = precision_score(y_true, y_pred, zero_division=0)
        metrics[f"{group}_recall"] = recall_score(y_true, y_pred, zero_division=0)
        metrics[f"{group}_f1"] = f1_score(y_true, y_pred, zero_division=0)
        metrics[f"{group}_pr_auc"] = (
            average_precision_score(y_true, rank_score) if y_true.sum() > 0 else np.nan
        )

    y_true_any = (ground_truth != 0).astype(int)
    y_pred_any = (predicted != 0).astype(int)
    metrics["any_precision"] = precision_score(y_true_any, y_pred_any, zero_division=0)
    metrics["any_recall"] = recall_score(y_true_any, y_pred_any, zero_division=0)
    metrics["any_f1"] = f1_score(y_true_any, y_pred_any, zero_division=0)
    return metrics


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", required=True, help="Directory with st_ground_truth.csv + *_predictions.csv")
    parser.add_argument("--scenario", required=True)
    parser.add_argument("--out-csv", default=None, help="Append results here instead of just printing")
    args = parser.parse_args()

    ground_truth = pd.read_csv(os.path.join(args.data_dir, "st_ground_truth.csv"), index_col=0)["ground_truth"]

    rows = []
    for pred_path in sorted(glob.glob(os.path.join(args.data_dir, "*_predictions.csv"))):
        method = os.path.basename(pred_path).replace("_predictions.csv", "")
        pred_df = pd.read_csv(pred_path)
        metrics = score_predictions(ground_truth, pred_df)
        metrics["scenario"] = args.scenario
        metrics["method"] = method
        rows.append(metrics)

    results = pd.DataFrame(rows).set_index(["scenario", "method"])
    pd.set_option("display.width", 160)
    print(results.round(3))

    if args.out_csv:
        os.makedirs(os.path.dirname(args.out_csv), exist_ok=True)
        header = not os.path.exists(args.out_csv)
        results.to_csv(args.out_csv, mode="a", header=header)


if __name__ == "__main__":
    main()
