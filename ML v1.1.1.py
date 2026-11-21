import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import sys
from collections import Counter
from scipy.stats import rankdata

# ---------------------------------------------------------
# IMPORT SCIKIT-LEARN & SHAP TOOLS
# ---------------------------------------------------------
print("\n🔍 Running GBM ML Pipeline v1.1.1 from:", sys.executable, "\n")
try:
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.svm import SVC
    from sklearn.feature_selection import RFE
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import StratifiedKFold
    from sklearn.metrics import roc_curve, auc, confusion_matrix, accuracy_score
    import shap
except ImportError:
    print("\n❌ CRITICAL ERROR: Library not found. Please run: pip install scikit-learn shap pandas numpy matplotlib scipy")
    exit()

# ==============================================================================
# CONFIGURATION & TOPOLOGICAL WEIGHTS
# ==============================================================================
script_dir = os.path.dirname(os.path.abspath(__file__))
MASTER_FILE = os.path.join(script_dir, "MASTER_GBM_DATASET.csv")
OUTPUT_DIR = os.path.join(script_dir, "ML_Outputs")
os.makedirs(OUTPUT_DIR, exist_ok=True)
NUM_TOP_FEATURES = 4
ITERATIONS = 10  # Number of balanced random undersampling iterations
FINAL_PANEL = ['hsa-miR-769-3p', 'hsa-miR-485-3p', 'hsa-miR-137', 'hsa-miR-765']

# Network degree centrality (C_degree) weights for target genes
TOPOLOGY_DICT = {
    'hsa-miR-137': 6.543,
    'hsa-miR-485-3p': 3.727,
    'hsa-miR-765': 3.771,
    'hsa-miR-769-3p': 2.273
}

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================
def bootstrap_auc_ci(y_true, y_score, n_bootstraps=1000, seed=42):
    """Compute 95% bootstrap confidence interval for ROC-AUC."""
    rng = np.random.RandomState(seed)
    y_true = np.asarray(y_true)
    y_score = np.asarray(y_score)
    bootstrapped_aucs = []

    for _ in range(n_bootstraps):
        indices = rng.randint(0, len(y_true), len(y_true))
        if len(np.unique(y_true[indices])) < 2:
            continue
        fpr_b, tpr_b, _ = roc_curve(y_true[indices], y_score[indices])
        bootstrapped_aucs.append(auc(fpr_b, tpr_b))

    if not bootstrapped_aucs:
        return np.nan, np.nan
    return np.percentile(bootstrapped_aucs, 2.5), np.percentile(bootstrapped_aucs, 97.5)


def run_feature_selection(X_train, y_train):
    """Nested consensus + topologically-weighted feature selection."""
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)

    # 1. Random Forest Rank
    rf_fs = RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=-1)
    rf_fs.fit(X_train, y_train)
    rf_ranks = rankdata(-rf_fs.feature_importances_)

    # 2. LASSO Rank
    lasso = LogisticRegression(penalty='l1', solver='liblinear', random_state=42, max_iter=1000)
    lasso.fit(X_train_scaled, y_train)
    lasso_ranks = rankdata(-np.abs(lasso.coef_[0]))

    # 3. SVM-RFE Rank
    svm = SVC(kernel='linear', random_state=42)
    rfe = RFE(estimator=svm, n_features_to_select=1, step=0.1)
    rfe.fit(X_train_scaled, y_train)
    svm_ranks = rfe.ranking_

    # Consensus & Topologically Weighted Ranks
    base_consensus_rank = (rf_ranks + lasso_ranks + svm_ranks) / 3.0
    features = X_train.columns
    topo_weights = np.array([TOPOLOGY_DICT.get(f, 1.0) for f in features])
    weighted_scores = base_consensus_rank / topo_weights

    return {
        'RF_Alone': features[np.argsort(rf_ranks)[:NUM_TOP_FEATURES]],
        'LASSO_Alone': features[np.argsort(lasso_ranks)[:NUM_TOP_FEATURES]],
        'SVM_Alone': features[np.argsort(svm_ranks)[:NUM_TOP_FEATURES]],
        'Consensus_Unweighted': features[np.argsort(base_consensus_rank)[:NUM_TOP_FEATURES]],
        'Consensus_Weighted': features[np.argsort(weighted_scores)[:NUM_TOP_FEATURES]]
    }


# ==============================================================================
# PART 1: RECURSIVE MERGE OF ALL *_ML_Ready.csv DATASETS
# ==============================================================================
print("==========================================================================")
print("📂 PART 1: MERGING INDIVIDUAL DATASETS")
print("==========================================================================")

search_pattern = os.path.join(script_dir, "**", "*_ML_Ready.csv")
all_files = sorted(glob.glob(search_pattern, recursive=True))

if not all_files:
    print("\n❌ NO *_ML_Ready.csv FILES FOUND!")
    print("Ensure you have run the R scripts for all 7 datasets first.")
    exit()

print(f"✅ Found {len(all_files)} dataset(s):")
dataframes = []

for filename in all_files:
    rel_path = os.path.relpath(filename, script_dir)
    dataset_name = os.path.basename(os.path.dirname(filename))
    if not dataset_name or dataset_name == os.path.basename(script_dir):
        dataset_name = os.path.basename(filename).replace("_ML_Ready.csv", "")

    print(f"   - {rel_path}  [Dataset: {dataset_name}]")
    df_temp = pd.read_csv(filename)
    df_temp['Sample_ID'] = df_temp['Sample_ID'].astype(str)
    df_temp['Dataset'] = dataset_name
    dataframes.append(df_temp)

print("\n🔄 Merging datasets (Inner Join for universally detected miRNAs)...")
df = pd.concat(dataframes, axis=0, join='inner', ignore_index=True)

cols = [c for c in df.columns if c not in ['Sample_ID', 'Dataset', 'Target']]
df = df[['Sample_ID', 'Dataset'] + cols + ['Target']]

df = df.sample(frac=1, random_state=42).reset_index(drop=True)
df.to_csv(MASTER_FILE, index=False)

print(f"🎉 Master dataset saved to: {MASTER_FILE}")
print(f"📊 Total Samples: {df.shape[0]} | 🧬 Common miRNAs: {len(cols)} | 📁 Datasets: {df['Dataset'].nunique()}\n")

if len(cols) < 5:
    print("⚠️  WARNING: Very few common miRNAs found. Check miRNA naming conventions.")

# ==============================================================================
# PART 2: INTERNAL 5-FOLD CV & ABLATION STUDY
# ==============================================================================
print("==========================================================================")
print("⚙️  PART 2: INTERNAL 5-FOLD CV & ABLATION STUDY")
print("==========================================================================")

df_normal = df[df['Target'] == 0]
df_tumor = df[df['Target'] == 1]
num_normal = len(df_normal)
print(f"📊 Pooled Data: {num_normal} Normal samples, {len(df_tumor)} Tumor samples.")

ablation_metrics = {
    name: {'acc': [], 'sens': [], 'spec': [], 'auc': [], 'selected': []}
    for name in ['RF_Alone', 'LASSO_Alone', 'SVM_Alone', 'Consensus_Unweighted', 'Consensus_Weighted']
}

for i in range(ITERATIONS):
    print(f"   Running Internal CV Iteration {i+1}/{ITERATIONS}...")
    df_tumor_sampled = df_tumor.sample(n=num_normal, random_state=42 + i)
    df_balanced = pd.concat([df_normal, df_tumor_sampled]).sample(frac=1, random_state=42 + i).reset_index(drop=True)

    X_bal = df_balanced.drop(columns=['Sample_ID', 'Dataset', 'Target'])
    y_bal = df_balanced['Target']

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42 + i)

    for train_idx, test_idx in cv.split(X_bal, y_bal):
        X_train, X_test = X_bal.iloc[train_idx], X_bal.iloc[test_idx]
        y_train, y_test = y_bal.iloc[train_idx], y_bal.iloc[test_idx]

        panels = run_feature_selection(X_train, y_train)

        for model_name, selected_features in panels.items():
            ablation_metrics[model_name]['selected'].extend(selected_features)

            eval_model = RandomForestClassifier(n_estimators=100, random_state=42)
            eval_model.fit(X_train[selected_features], y_train)
            y_proba = eval_model.predict_proba(X_test[selected_features])[:, 1]

            if len(np.unique(y_test)) > 1:
                fpr, tpr, thresholds = roc_curve(y_test, y_proba)
                optimal_threshold = thresholds[np.argmax(tpr - fpr)]
                roc_auc = auc(fpr, tpr)
            else:
                optimal_threshold = 0.5
                roc_auc = 0.5

            y_pred = (y_proba >= optimal_threshold).astype(int)
            tn, fp, fn, tp = confusion_matrix(y_test, y_pred, labels=[0, 1]).ravel()

            ablation_metrics[model_name]['acc'].append(accuracy_score(y_test, y_pred))
            ablation_metrics[model_name]['sens'].append(tp / (tp + fn) if (tp + fn) > 0 else 0)
            ablation_metrics[model_name]['spec'].append(tn / (tn + fp) if (tn + fp) > 0 else 0)
            ablation_metrics[model_name]['auc'].append(roc_auc)

print("\n📊 ABLATION STUDY RESULTS (Internal 5-Fold CV):")
print(f"{'Model Configuration':<25} | {'AUC':<7} | {'Accuracy':<9} | {'Sensitivity':<11} | {'Specificity':<11}")
print("-" * 74)
for model_name in ablation_metrics:
    print(f"{model_name:<25} | {np.mean(ablation_metrics[model_name]['auc']):.4f}  | "
          f"{np.mean(ablation_metrics[model_name]['acc'])*100:.2f}%    | "
          f"{np.mean(ablation_metrics[model_name]['sens'])*100:.2f}%      | "
          f"{np.mean(ablation_metrics[model_name]['spec'])*100:.2f}%")

top_overall = Counter(ablation_metrics['Consensus_Weighted']['selected']).most_common(10)
print("\n🏆 TOP STABLE miRNA BIOMARKERS (Topologically Weighted Consensus):")
for miRNA, count in top_overall:
    print(f"   - {miRNA}: Selected {count} times")

# ==============================================================================
# PART 3: LEAVE-ONE-DATASET-OUT (LODO) & POOLED LODO VALIDATION
# ==============================================================================
print("\n==========================================================================")
print("🌍 PART 3: LEAVE-ONE-DATASET-OUT (LODO) & POOLED LODO VALIDATION")
print("==========================================================================")

datasets = sorted(df['Dataset'].unique())
total_samples_all = len(df)
lodo_summary = []
pooled_y_true = []
pooled_y_proba = []

plt.close('all')
fig_lodo, ax_lodo = plt.subplots(figsize=(9, 7))

for test_ds in datasets:
    train_df = df[df['Dataset'] != test_ds]
    test_df = df[df['Dataset'] == test_ds]

    X_test_lodo = test_df.drop(columns=['Sample_ID', 'Dataset', 'Target'])
    y_test_lodo = test_df['Target'].values

    train_norm = train_df[train_df['Target'] == 0]
    train_tum = train_df[train_df['Target'] == 1]
    n_bal = min(len(train_norm), len(train_tum))

    proba_iterations = []
    iter_accs, iter_sens, iter_specs = [], [], []

    for it in range(ITERATIONS):
        tr_n_samp = train_norm.sample(n=n_bal, random_state=42 + it)
        tr_t_samp = train_tum.sample(n=n_bal, random_state=42 + it)
        tr_bal = pd.concat([tr_n_samp, tr_t_samp]).sample(frac=1, random_state=42 + it)

        X_tr = tr_bal.drop(columns=['Sample_ID', 'Dataset', 'Target'])
        y_tr = tr_bal['Target']

        # Nested feature selection strictly on the 6 training datasets
        panels = run_feature_selection(X_tr, y_tr)
        selected_feats = panels['Consensus_Weighted']

        clf = RandomForestClassifier(n_estimators=100, random_state=42 + it)
        clf.fit(X_tr[selected_feats], y_tr)

        y_proba = clf.predict_proba(X_test_lodo[selected_feats])[:, 1]
        proba_iterations.append(y_proba)

        # Strict 0.5 threshold (no test-set threshold tuning to prevent leakage)
        y_pred = (y_proba >= 0.5).astype(int)
        tn, fp, fn, tp = confusion_matrix(y_test_lodo, y_pred, labels=[0, 1]).ravel()

        iter_accs.append(accuracy_score(y_test_lodo, y_pred))
        iter_sens.append(tp / (tp + fn) if (tp + fn) > 0 else np.nan)
        iter_specs.append(tn / (tn + fp) if (tn + fp) > 0 else np.nan)

    # Ensemble out-of-fold probability for this held-out dataset
    mean_y_proba = np.mean(proba_iterations, axis=0)
    pooled_y_true.extend(y_test_lodo)
    pooled_y_proba.extend(mean_y_proba)

    if len(np.unique(y_test_lodo)) > 1:
        fpr_ds, tpr_ds, _ = roc_curve(y_test_lodo, mean_y_proba)
        ds_auc = auc(fpr_ds, tpr_ds)
        ci_low, ci_high = bootstrap_auc_ci(y_test_lodo, mean_y_proba)
        ax_lodo.plot(fpr_ds, tpr_ds, lw=1.6, alpha=0.75,
                     label=f"{test_ds} (n={len(y_test_lodo)}, AUC={ds_auc:.3f})")
    else:
        ds_auc, ci_low, ci_high = np.nan, np.nan, np.nan

    ds_weight = (len(y_test_lodo) / total_samples_all) * 100

    lodo_summary.append({
        'Held_Out_Dataset': test_ds,
        'N_Total': len(y_test_lodo),
        'Weight_%': round(ds_weight, 2),
        'N_Normal': int(np.sum(y_test_lodo == 0)),
        'N_GBM': int(np.sum(y_test_lodo == 1)),
        'AUC': round(ds_auc, 4) if not np.isnan(ds_auc) else np.nan,
        'AUC_95_CI': f"[{ci_low:.3f}, {ci_high:.3f}]" if not np.isnan(ci_low) else "N/A",
        'Accuracy_%': round(np.nanmean(iter_accs) * 100, 2),
        'Sensitivity_%': round(np.nanmean(iter_sens) * 100, 2),
        'Specificity_%': round(np.nanmean(iter_specs) * 100, 2)
    })

# Pooled LODO Metrics (all out-of-fold predictions combined)
pooled_y_true = np.array(pooled_y_true)
pooled_y_proba = np.array(pooled_y_proba)
fpr_pooled, tpr_pooled, _ = roc_curve(pooled_y_true, pooled_y_proba)
pooled_auc = auc(fpr_pooled, tpr_pooled)
p_ci_low, p_ci_high = bootstrap_auc_ci(pooled_y_true, pooled_y_proba)

pooled_y_pred = (pooled_y_proba >= 0.5).astype(int)
tn_p, fp_p, fn_p, tp_p = confusion_matrix(pooled_y_true, pooled_y_pred, labels=[0, 1]).ravel()
pooled_acc = accuracy_score(pooled_y_true, pooled_y_pred) * 100
pooled_sens = (tp_p / (tp_p + fn_p)) * 100
pooled_spec = (tn_p / (tn_p + fp_p)) * 100

lodo_df = pd.DataFrame(lodo_summary)

# Sample-size-weighted mean AUC
valid_lodo = lodo_df.dropna(subset=['AUC'])
weighted_mean_auc = np.average(valid_lodo['AUC'], weights=valid_lodo['N_Total'])

# Append Pooled LODO summary row to CSV
pooled_row = pd.DataFrame([{
    'Held_Out_Dataset': 'POOLED_LODO_ALL',
    'N_Total': len(pooled_y_true),
    'Weight_%': 100.00,
    'N_Normal': int(np.sum(pooled_y_true == 0)),
    'N_GBM': int(np.sum(pooled_y_true == 1)),
    'AUC': round(pooled_auc, 4),
    'AUC_95_CI': f"[{p_ci_low:.3f}, {p_ci_high:.3f}]",
    'Accuracy_%': round(pooled_acc, 2),
    'Sensitivity_%': round(pooled_sens, 2),
    'Specificity_%': round(pooled_spec, 2)
}])
lodo_df_export = pd.concat([lodo_df, pooled_row], ignore_index=True)

lodo_csv_path = os.path.join(script_dir, "LODO_Validation_Results.csv")
lodo_df_export.to_csv(lodo_csv_path, index=False)

print(lodo_df_export.to_string(index=False))
print("-" * 98)
print(f"📌 UNWEIGHTED MEAN LODO AUC : {lodo_df['AUC'].mean():.4f} ± {lodo_df['AUC'].std():.4f}")
print(f"📌 WEIGHTED MEAN LODO AUC   : {weighted_mean_auc:.4f}")
print(f"📌 POOLED LODO AUC (95% CI) : {pooled_auc:.4f} [{p_ci_low:.3f}, {p_ci_high:.3f}]")
print(f"✅ Saved LODO results table to: {lodo_csv_path}")

# Finalize LODO ROC Plot (with bold Pooled LODO ROC)
ax_lodo.plot([0, 1], [0, 1], linestyle="--", lw=1.5, color="black", alpha=0.6, label="Chance (AUC = 0.500)")
ax_lodo.plot(fpr_pooled, tpr_pooled, color="navy", lw=3.2,
             label=f"Pooled LODO ROC (AUC = {pooled_auc:.3f}, 95% CI: {p_ci_low:.3f}-{p_ci_high:.3f})")

ax_lodo.set_xlim([-0.05, 1.05])
ax_lodo.set_ylim([-0.05, 1.05])
ax_lodo.set_xlabel("False Positive Rate (1 - Specificity)", fontsize=12, fontweight='bold')
ax_lodo.set_ylabel("True Positive Rate (Sensitivity)", fontsize=12, fontweight='bold')
ax_lodo.set_title("Leave-One-Dataset-Out (LODO) & Pooled Validation ROC Curves", fontsize=13, fontweight='bold', pad=15)
ax_lodo.legend(loc="lower right", fontsize=8.5)
plt.tight_layout()

plt.savefig(os.path.join(OUTPUT_DIR, "LODO_ROC_Curves.png"), dpi=600, bbox_inches='tight')
plt.savefig(os.path.join(OUTPUT_DIR, "LODO_ROC_Curves.svg"), bbox_inches='tight')
plt.close()
print("✅ Saved LODO_ROC_Curves.png & .svg")

# ==============================================================================
# PART 4: EXPLAINABLE AI (SHAP) & HIGH-DEFINITION INTERNAL ROC CURVE
# ==============================================================================
print("\n==========================================================================")
print("🧠 PART 4: SHAP EXPLAINABILITY & HIGH-DEFINITION ROC")
print("==========================================================================")

df_tumor_shap = df_tumor.sample(n=num_normal, random_state=42)
df_balanced_shap = pd.concat([df_normal, df_tumor_shap]).sample(frac=1, random_state=42).reset_index(drop=True)

X_final = df_balanced_shap[FINAL_PANEL]
y_final = df_balanced_shap['Target']

final_model = RandomForestClassifier(n_estimators=100, random_state=42)
final_model.fit(X_final, y_final)

explainer = shap.TreeExplainer(final_model)
shap_values = explainer.shap_values(X_final)

if isinstance(shap_values, list):
    shap_values_gbm = shap_values[1]
elif isinstance(shap_values, np.ndarray) and shap_values.ndim == 3:
    shap_values_gbm = shap_values[:, :, 1]
else:
    shap_values_gbm = shap_values

# SHAP Summary Plot
plt.close('all')
plt.figure(figsize=(9, 4.5))
shap.summary_plot(shap_values_gbm, X_final, show=False)
plt.title("SHAP Summary Plot: Directional Impact on GBM Prediction", fontsize=12, fontweight='bold', pad=15)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "SHAP_Summary_Plot.png"), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(OUTPUT_DIR, "SHAP_Summary_Plot.svg"), bbox_inches='tight')
plt.close()
print("✅ Saved SHAP_Summary_Plot.png & .svg")

# SHAP Dependence Plots
for feature_name in FINAL_PANEL:
    plt.close('all')
    fig, ax = plt.subplots(figsize=(7, 5))
    shap.dependence_plot(feature_name, shap_values_gbm, X_final, interaction_index=None, ax=ax, show=False)
    plt.title(f"SHAP Feature Dependence: {feature_name}", fontsize=12, fontweight='bold', pad=15)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f"SHAP_Dependence_{feature_name}.png"), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(OUTPUT_DIR, f"SHAP_Dependence_{feature_name}.svg"), bbox_inches='tight')
    plt.close()
    print(f"✅ Saved SHAP_Dependence_{feature_name}.png & .svg")

# High-Definition Internal 5-Fold ROC Curve
cv_roc = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
classifier_roc = RandomForestClassifier(n_estimators=100, random_state=42)

tprs, aucs = [], []
mean_fpr = np.linspace(0, 1, 100)

plt.close('all')
fig, ax = plt.subplots(figsize=(8, 6))

for i, (train_idx, test_idx) in enumerate(cv_roc.split(X_final, y_final)):
    classifier_roc.fit(X_final.iloc[train_idx], y_final.iloc[train_idx])
    y_proba_roc = classifier_roc.predict_proba(X_final.iloc[test_idx])[:, 1]

    fpr_fold, tpr_fold, _ = roc_curve(y_final.iloc[test_idx], y_proba_roc)
    roc_auc_fold = auc(fpr_fold, tpr_fold)
    aucs.append(roc_auc_fold)
    ax.plot(fpr_fold, tpr_fold, lw=1, alpha=0.3, label=f"ROC fold {i+1} (AUC = {roc_auc_fold:.3f})")

    interp_tpr = np.interp(mean_fpr, fpr_fold, tpr_fold)
    interp_tpr[0] = 0.0
    tprs.append(interp_tpr)

ax.plot([0, 1], [0, 1], linestyle="--", lw=2, color="red", label="Chance (AUC = 0.500)", alpha=0.8)
mean_tpr = np.mean(tprs, axis=0)
mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)
std_auc = np.std(aucs)

ax.plot(mean_fpr, mean_tpr, color="blue", label=f"Mean ROC (AUC = {mean_auc:.4f} ± {std_auc:.3f})", lw=2.5, alpha=0.8)
std_tpr = np.std(tprs, axis=0)
ax.fill_between(mean_fpr, np.maximum(mean_tpr - std_tpr, 0), np.minimum(mean_tpr + std_tpr, 1),
                color="grey", alpha=0.2, label=r"$\pm$ 1 std. dev.")

ax.set_xlim([-0.05, 1.05])
ax.set_ylim([-0.05, 1.05])
ax.set_xlabel("False Positive Rate (1 - Specificity)", fontsize=12, fontweight='bold')
ax.set_ylabel("True Positive Rate (Sensitivity)", fontsize=12, fontweight='bold')
ax.set_title("High-Definition ROC Curve: Topologically-Weighted 4-miRNA Panel", fontsize=14, fontweight='bold', pad=15)
ax.legend(loc="lower right", fontsize=10)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "High_Def_ROC_Curve.png"), dpi=600, bbox_inches='tight')
plt.savefig(os.path.join(OUTPUT_DIR, "High_Def_ROC_Curve.svg"), bbox_inches='tight')
plt.close()
print("✅ Saved High_Def_ROC_Curve.png & .svg")

print("\n🎉 Complete Pipeline v1.1.1 Execution Finished!")
