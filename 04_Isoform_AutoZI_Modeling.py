# =============================================
# AUTOZI SCRIPT FOR ISOFORM-LEVEL ANALYSIS
# =============================================

import os
import numpy as np
import pandas as pd
import scanpy as sc
import muon
from anndata import read_h5ad
from scipy.stats import beta
from scvi.model import AUTOZI
import scipy.sparse

# === Load and Prepare Gene-Level Data ===
Int_folder = "Intermediate_Files/QC"

# Load concatenated, QC-filtered isoform-level data (PBMC1 + PBMC2)
adata_g_filter = read_h5ad(os.path.join(Int_folder, "Concatenated_Iso_Data.h5ad"))

# Quick sanity check: batches should be present and balanced
print(adata_g_filter.obs["batch"].value_counts())

# Ensure batch column is categorical
adata_g_filter.obs["batch"] = adata_g_filter.obs["batch"].astype("category")

# Store raw counts in dedicated layer ('counts') and convert to dense
adata_g_filter.layers["counts"] = (
    adata_g_filter.X.toarray() if scipy.sparse.issparse(adata_g_filter.X) else adata_g_filter.X.copy()
)

# Wrap into MuData for consistency
mdata_g_std = muon.MuData({
    "rna": adata_g_filter.copy(),
    "log_norm_rna": adata_g_filter.copy()
})


# Make sure the active matrix for 'rna' is the counts layer
mdata_g_std.mod["rna"].X = mdata_g_std.mod["rna"].layers["counts"]

# Drop .raw if present to avoid ambiguity
if mdata_g_std.mod["rna"].raw is not None:
    del mdata_g_std.mod["rna"].raw

rna_adata_g_std = mdata_g_std.mod['rna']
rna_adata_g_std.layers["counts"] = rna_adata_g_std.X.copy()

# === Setup and Train AUTOZI ===
AUTOZI.setup_anndata(rna_adata_g_std, batch_key="batch")

# Initialize the AutoZI model with tuned parameters
autozi_g_std = AUTOZI(
    adata=rna_adata_g_std,
    n_latent=30,
    dropout_rate=0.4
)

# Train with early stopping and tuned hyperparameters
autozi_g_std.train(
    max_epochs=300,
    train_size=0.9,
    early_stopping=True,
    early_stopping_patience=30,
    batch_size=64,
    plan_kwargs={'lr': 5e-3, 'weight_decay': 1e-2}
)

# === Save Latents and ZI Output ===
latent = autozi_g_std.get_latent_representation()
posterior = autozi_g_std.get_alphas_betas()
alpha = posterior['alpha_posterior']
beta_ = posterior['beta_posterior']

rna_adata_g_std.obsm['X_AutoZI'] = latent

# Zero-inflation probabilities
zi_probs = beta.cdf(0.5, alpha, beta_)
is_zi_pred = zi_probs > 0.5

print("Fraction of predicted ZI genes:", is_zi_pred.mean())

# Expression filter (for interpretability stats)
mask_expr = (rna_adata_g_std.X.mean(axis=0) > 1).reshape(-1)
print("Fraction expressed >1:", mask_expr.mean())
print("Fraction ZI & expressed >1:", is_zi_pred[mask_expr].mean())

# === Save Model and Annotated Data ===
model_dir = "Modeling/PBMCs_AutoZI/Isoform"
os.makedirs(model_dir, exist_ok=True)

# Save trained model weights/config (reloadable with AUTOZI.load)
autozi_g_std.save(model_dir, overwrite=True)

# Save MuData with latent embedding for easy reuse downstream
mdata_g_std.mod['rna'] = rna_adata_g_std
mdata_g_std.write("Modeling/PBMCs_AutoZI/Isoform/mdata_iso_with_latent_PBMCs_AutoZI.h5mu")

print("✅ AutoZI model and isoform-level MuData saved successfully.")
