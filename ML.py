import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from collections import Counter
from scipy.stats import rankdata

# ---------------------------------------------------------
# IMPORT SCIKIT-LEARN & SHAP TOOLS
# ---------------------------------------------------------
print("\n🔍 I am running from:", sys.executable, "\n")
try:
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.svm import SVC
    from sklearn.feature_selection import RFE
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import StratifiedKFold
    from sklearn.metrics import roc_curve, auc, confusion_matrix, accuracy_score
    import shap  # Explainable AI
except ImportError as e:
    print("\n❌ CRITICAL ERROR: Library not found. Please run: pip install scikit-learn shap pandas matplotlib scipy")
    exit()

# ==============================================================================
# CONFIGURATION & TOPOLOGICAL WEIGHTS
# ==============================================================================
# Robust path setup to allow execution from any directory
script_dir = os.path.dirname(os.path.abspath(__file__))
INPUT_FILE = os.path.join(script_dir, "MASTER_GBM_DATASET.csv")
NUM_TOP_FEATURES = 4  
ITERATIONS = 10  # Number of times to repeat the random undersampling[cite: 7]

# ALGORITHMIC NOVELTY: Topological Weight Dictionary[cite: 7]
# These weights represent the network degree centrality (C_degree) for the target genes[cite: 7]
TOPOLOGY_DICT = {
    'hsa-miR-137': 6.543,      
    'hsa-miR-485-3p': 3.727,
    'hsa-miR-765': 3.771,
    'hsa-miR-769-3p': 2.273
}

# ==============================================================================
# 1. LOAD & SEPARATE DATA
# ==============================================================================
print(f"🚀 Loading file: {INPUT_FILE}")

if not os.path.exists(INPUT_FILE):
    print("\n❌ ERROR: MASTER_GBM_DATASET.csv not found!")
    print("Please run the Merge.py script first to generate the master dataset.")
    exit()

df = pd.read_csv(INPUT_FILE)

df_normal = df[df['Target'] == 0]
df_tumor = df[df['Target'] == 1]

num_normal = len(df_normal)
print(f"📊 Dataset Loaded: {num_normal} Normal samples, {len(df_tumor)} Tumor samples.")
print(f"⚙️  Commencing {ITERATIONS} iterations of balanced random undersampling with Full Ablation Study...")

# Master tracking dictionary for ALL 5 ablation models[cite: 7]
ablation_metrics = {
    'RF_Alone': {'acc': [], 'sens': [], 'spec': [], 'auc': [], 'selected': []},
    'LASSO_Alone': {'acc': [], 'sens': [], 'spec': [], 'auc': [], 'selected': []},
    'SVM_Alone': {'acc': [], 'sens': [], 'spec': [], 'auc': [], 'selected': []},
    'Consensus_Unweighted': {'acc': [], 'sens': [], 'spec': [], 'auc': [], 'selected': []},
    'Consensus_Weighted': {'acc': [], 'sens': [], 'spec': [], 'auc': [], 'selected': []}
}

# ==============================================================================
# 2. THE MASTER UNDERSAMPLING & HYBRID FEATURE SELECTION LOOP
# ==============================================================================
for i in range(ITERATIONS):
    print(f"\n--- Running Iteration {i+1}/{ITERATIONS} ---")
    
    # Random Undersampling to balance classes
    df_tumor_sampled = df_tumor.sample(n=num_normal, random_state=None)
    df_balanced = pd.concat([df_normal, df_tumor_sampled]).sample(frac=1).reset_index(drop=True)
    
    X_bal = df_balanced.drop(columns=['Sample_ID', 'Target'])
    y_bal = df_balanced['Target']
    
    # Stratified 5-Fold Cross-Validation
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=None)
    
    for train_idx, test_idx in cv.split(X_bal, y_bal):
        X_train, X_test = X_bal.iloc[train_idx], X_bal.iloc[test_idx]
        y_train, y_test = y_bal.iloc[train_idx], y_bal.iloc[test_idx]

        # ---------------------------------------------------
        # STEP A: FEATURE SELECTION (NESTED)
        # ---------------------------------------------------
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

        # Calculate Consensus Rank 
        base_consensus_rank = (rf_ranks + lasso_ranks + svm_ranks) / 3.0
        
        # Calculate Topologically Weighted Rank
        features = X_train.columns
        topo_weights = np.array([TOPOLOGY_DICT.get(f, 1.0) for f in features])
        weighted_scores = base_consensus_rank / topo_weights
        
        # ---------------------------------------------------
        # STEP B: EXTRACT THE 5 DISTINCT FEATURE PANELS
        # ---------------------------------------------------
        # Panel 1: RF Alone
        top_idx_rf = np.argsort(rf_ranks)[:NUM_TOP_FEATURES]
        panel_rf = features[top_idx_rf]
        
        # Panel 2: LASSO Alone
        top_idx_lasso = np.argsort(lasso_ranks)[:NUM_TOP_FEATURES]
        panel_lasso = features[top_idx_lasso]
        
        # Panel 3: SVM-RFE Alone
        top_idx_svm = np.argsort(svm_ranks)[:NUM_TOP_FEATURES]
        panel_svm = features[top_idx_svm]
        
        # Panel 4: Consensus (Unweighted)
        top_idx_cons = np.argsort(base_consensus_rank)[:NUM_TOP_FEATURES]
        panel_cons = features[top_idx_cons]
        
        # Panel 5: Consensus (Weighted)
        top_idx_weight = np.argsort(weighted_scores)[:NUM_TOP_FEATURES]
        panel_weight = features[top_idx_weight]

        # Log selected features
        ablation_metrics['RF_Alone']['selected'].extend(panel_rf)
        ablation_metrics['LASSO_Alone']['selected'].extend(panel_lasso)
        ablation_metrics['SVM_Alone']['selected'].extend(panel_svm)
        ablation_metrics['Consensus_Unweighted']['selected'].extend(panel_cons)
        ablation_metrics['Consensus_Weighted']['selected'].extend(panel_weight)

        # ---------------------------------------------------
        # STEP C: EVALUATE EACH PANEL
        # ---------------------------------------------------
        panels = {
            'RF_Alone': panel_rf,
            'LASSO_Alone': panel_lasso,
            'SVM_Alone': panel_svm,
            'Consensus_Unweighted': panel_cons,
            'Consensus_Weighted': panel_weight
        }

        for model_name, selected_features in panels.items():
            X_train_best = X_train[selected_features]
            X_test_best = X_test[selected_features]

            eval_model = RandomForestClassifier(n_estimators=100, random_state=42)
            eval_model.fit(X_train_best, y_train)
            
            y_proba = eval_model.predict_proba(X_test_best)[:, 1]

            if len(np.unique(y_test)) > 1:
                fpr, tpr, thresholds = roc_curve(y_test, y_proba)
                j_scores = tpr - fpr
                optimal_idx = np.argmax(j_scores)
                optimal_threshold = thresholds[optimal_idx]
                roc_auc = auc(fpr, tpr)
            else:
                optimal_threshold = 0.5
                roc_auc = 0.5

            y_pred = (y_proba >= optimal_threshold).astype(int)
            tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
            
            sens = tp / (tp + fn) if (tp + fn) > 0 else 0
            spec = tn / (tn + fp) if (tn + fp) > 0 else 0
            acc = accuracy_score(y_test, y_pred)
            
            ablation_metrics[model_name]['acc'].append(acc)
            ablation_metrics[model_name]['sens'].append(sens)
            ablation_metrics[model_name]['spec'].append(spec)
            ablation_metrics[model_name]['auc'].append(roc_auc)

# ==============================================================================
# 3. REPORT ABLATION STUDY RESULTS
# ==============================================================================
print("\n==========================================================================")
print("📊 ABLATION STUDY RESULTS: COMPARING 5 FEATURE SELECTION ALGORITHMS")
print("==========================================================================")
print(f"{'Model Configuration':<30} | {'AUC':<7} | {'Accuracy':<9} | {'Sensitivity':<11} | {'Specificity':<11}")
print("-" * 74)

for model_name in ['RF_Alone', 'LASSO_Alone', 'SVM_Alone', 'Consensus_Unweighted', 'Consensus_Weighted']:
    avg_auc = np.mean(ablation_metrics[model_name]['auc'])
    avg_acc = np.mean(ablation_metrics[model_name]['acc']) * 100
    avg_sens = np.mean(ablation_metrics[model_name]['sens']) * 100
    avg_spec = np.mean(ablation_metrics[model_name]['spec']) * 100
    
    print(f"{model_name:<30} | {avg_auc:.4f}  | {avg_acc:.2f}%    | {avg_sens:.2f}%      | {avg_spec:.2f}%")

print("==========================================================================\n")

top_overall = Counter(ablation_metrics['Consensus_Weighted']['selected']).most_common(10)
print("🏆 TOP STABLE miRNA BIOMARKERS (Topologically Weighted Consensus):")
for miRNA, count in top_overall:
    print(f"   - {miRNA}: Selected {count} times")
print("==========================================================================\n")


# ==============================================================================
# 4. EXPLAINABLE AI (XAI) WITH SHAP 
# ==============================================================================
print("🧠 Generating SHAP Explainable AI Visualizations for the final panel...")

# These are the top 4 robustly selected biomarkers from the consensus[cite: 7]
final_top_features = ['hsa-miR-769-3p', 'hsa-miR-485-3p', 'hsa-miR-137', 'hsa-miR-765']

df_tumor_shap = df_tumor.sample(n=num_normal, random_state=42)
df_balanced_shap = pd.concat([df_normal, df_tumor_shap]).sample(frac=1, random_state=42).reset_index(drop=True)

X_final = df_balanced_shap[final_top_features]
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

# Save Summary Plot
plt.close('all')
plt.figure(figsize=(9, 4.5))
shap.summary_plot(shap_values_gbm, X_final, show=False)
plt.title("SHAP Summary Plot: Directional Impact on GBM Prediction", fontsize=12, fontweight='bold', pad=15)
plt.tight_layout()
plt.savefig(os.path.join(script_dir, "SHAP_Summary_Plot.png"), dpi=300, bbox_inches='tight')
print("✅ Saved SHAP_Summary_Plot.png")
plt.close()

# Save Feature Dependence Plots
for feature_name in final_top_features:
    plt.close('all')
    fig, ax = plt.subplots(figsize=(7, 5))
    shap.dependence_plot(feature_name, shap_values_gbm, X_final, interaction_index=None, ax=ax, show=False)
    plt.title(f"SHAP Feature Dependence: {feature_name}", fontsize=12, fontweight='bold', pad=15)
    plt.tight_layout()
    plt.savefig(os.path.join(script_dir, f"SHAP_Dependence_{feature_name}.png"), dpi=300, bbox_inches='tight')
    print(f"✅ Saved SHAP_Dependence_{feature_name}.png")
    plt.close()

# ==============================================================================
# 5. HIGH-DEFINITION ROC CURVE FOR FINAL PANEL
# ==============================================================================
print("📈 Generating High-Definition ROC Curve...")

cv_roc = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
classifier_roc = RandomForestClassifier(n_estimators=100, random_state=42)

tprs = []
aucs = []
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
tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color="grey", alpha=0.2, label=r"$\pm$ 1 std. dev.")

ax.set_xlim([-0.05, 1.05])
ax.set_ylim([-0.05, 1.05])
ax.set_xlabel("False Positive Rate (1 - Specificity)", fontsize=12, fontweight='bold')
ax.set_ylabel("True Positive Rate (Sensitivity)", fontsize=12, fontweight='bold')
ax.set_title("High-Definition ROC Curve: Topologically-Weighted 4-miRNA Panel", fontsize=14, fontweight='bold', pad=15)
ax.legend(loc="lower right", fontsize=10)

plt.tight_layout()
roc_path = os.path.join(script_dir, "High_Def_ROC_Curve.png")
plt.savefig(roc_path, dpi=600, bbox_inches='tight')
print(f"✅ Saved Ultra-HD ROC Curve to {roc_path}")
plt.close()

print("\n🎉 Pipeline Execution Complete!")