import numpy as np
import umap
import hdbscan

def hdbscan_cluster(X):
  print("Embedding in 10 UMAP dimensions")
  clusterable_embedder = umap.UMAP(
    n_neighbors=30,
    min_dist=0.0,
    n_components=10,
    random_state=42
  )

  print("Finding the clusters using HDBSCAN")
  emb = clusterable_embedder.fit_transform(X)

  labels = hdbscan.HDBSCAN(
    min_samples=15,
    min_cluster_size=20).fit_predict(emb)

  return(labels)
