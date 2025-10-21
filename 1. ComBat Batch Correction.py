import pandas as pd
import numpy as np
import statsmodels.api as sm
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import confusion_matrix, classification_report, accuracy_score


# ==========================================================
# COMBAT BATCH CORRECTION
# Target: Group
# Predictors: MIR132, MIR137, MIR9, MIR941, MIR34A
# Batch: BATCH
# ==========================================================

import pandas as pd
import numpy as np
import statsmodels.api as sm
from neuroHarmonize import harmonizationLearn, harmonizationApply
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import confusion_matrix, classification_report, accuracy_score, roc_curve, auc
import matplotlib.pyplot as plt

# === 1️⃣ LOAD DATA ===
file_path = "/Users/mac/Desktop/Eugei vf 0810 final BATCH 091025.xlsx"
df = pd.read_excel(file_path, engine="openpyxl")

# Clean column names
df.columns = df.columns.str.strip()

print("✅ Data loaded successfully!")
print("Columns:", df.columns.tolist())

# === 2️⃣ DEFINE VARIABLES ===
target = "Group"
batch_col = "BATCH"
predictors = ["MIR132", "MIR137", "MIR9", "MIR941", "MIR34A"]

# Drop rows with missing data
df = df.dropna(subset=[target, batch_col] + predictors).reset_index(drop=True)

# === 3️⃣ SETUP DATA FOR COMBAT ===
# === 🧬 COMBAT BATCH CORRECTION (neuroCombat, compatible version) ===
# === 🧬 COMBAT BATCH CORRECTION (final, neuroCombat style) ===
from neuroCombat import neuroCombat

# Numeric features (samples × features)
X = df[predictors].apply(pd.to_numeric, errors="coerce")
X = X.to_numpy().T  # neuroCombat expects features × samples

# Build a covars DataFrame that INCLUDES the batch column
covars = pd.DataFrame({
    "batch": df[batch_col].astype(str).values,
    "Group": df[target].astype("category").cat.codes
})

print("\n🧬 Applying ComBat batch correction (neuroCombat)...")

# ✅ Specify which column in covars is the batch variable
combat_data = neuroCombat(
    dat=X,                 # features × samples
    covars=covars,         # must include batch column
    batch_col="batch",     # name of the batch column in covars
    parametric=True
)

# Extract corrected data
X_corrected = pd.DataFrame(combat_data["data"].T, columns=predictors)

print("✅ Batch correction complete!")

