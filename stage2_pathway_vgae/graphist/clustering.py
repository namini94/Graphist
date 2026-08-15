"""Post-hoc clustering comparisons (Mclust / KMeans vs. ground-truth annotation).

Mclust itself runs in R (outside this package, see ``stage1_bulk_regression_R``
and the original scripts' comments); this module just loads its label output
and scores it, plus provides an optional KMeans baseline (used only by the
Maynard dataset in the original scripts).
"""
import pandas as pd
from sklearn import metrics
from sklearn.cluster import KMeans


def load_mclust_labels(path: str) -> pd.Series:
    """Load a single-column Mclust label CSV (as exported by the R Stage-1/Stage-2 clustering)."""
    labels = pd.read_csv(path, header=0, index_col=0)
    return pd.Categorical(labels.iloc[:, 0])


def adjusted_rand_index(ground_truth, predicted_labels) -> float:
    return metrics.adjusted_rand_score(ground_truth, predicted_labels)


def kmeans_labels(embedding, n_clusters: int, random_state: int = 0):
    return KMeans(n_clusters=n_clusters, random_state=random_state, n_init="auto").fit_predict(embedding)
