import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset, random_split
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import os

# ------------------------------
# Configuration
# ------------------------------
RAW_DIR = "data/raw/"
MODEL_DIR = "models"
BATCH_SIZE = 32
EPOCHS = 2000
LEARNING_RATE = 1e-3
SEED = 42
VALID_SPLIT = 0.2

os.makedirs(MODEL_DIR, exist_ok=True)
torch.manual_seed(SEED)
np.random.seed(SEED)

# ------------------------------
# Load and merge all CSVs
# ------------------------------
all_files = glob.glob(RAW_DIR + "config_kappa_*.csv")
if len(all_files) == 0:
    raise RuntimeError(f"No files found in {RAW_DIR}")

dfs = [pd.read_csv(f, delim_whitespace=True, header=None, on_bad_lines='skip') for f in all_files]
df = pd.concat(dfs, ignore_index=True)

dataset = df.values.astype(np.float32)
X = dataset[:, :-1]  # angles
Y = dataset[:, -1:]  # kappa

# ------------------------------
# PyTorch dataset and split
# ------------------------------
X_tensor = torch.tensor(X)
Y_tensor = torch.tensor(Y)

full_dataset = TensorDataset(X_tensor, Y_tensor)
val_size = int(len(full_dataset) * VALID_SPLIT)
train_size = len(full_dataset) - val_size
train_ds, val_ds = random_split(full_dataset, [train_size, val_size])

train_loader = DataLoader(train_ds, batch_size=BATCH_SIZE, shuffle=True)
val_loader = DataLoader(val_ds, batch_size=BATCH_SIZE)

# ------------------------------
# Define simple feedforward model
# ------------------------------
class SimpleNN(nn.Module):
    def __init__(self, input_dim):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, 20),
            nn.ReLU(),
            nn.Linear(20, 10),
            nn.ReLU(),
            nn.Linear(10, 1)
        )

    def forward(self, x):
        return self.net(x)

model = SimpleNN(X.shape[1])

# ------------------------------
# Loss & optimizer
# ------------------------------
criterion = nn.MSELoss()
optimizer = optim.Adam(model.parameters(), lr=LEARNING_RATE)

# ------------------------------
# Training loop
# ------------------------------
train_losses, val_losses = [], []

for epoch in range(1, EPOCHS+1):
    model.train()
    epoch_loss = 0.0
    for xb, yb in train_loader:
        optimizer.zero_grad()
        out = model(xb)
        loss = criterion(out, yb)
        loss.backward()
        optimizer.step()
        epoch_loss += loss.item() * xb.size(0)
    train_losses.append(epoch_loss / train_size)

    model.eval()
    val_loss = 0.0
    with torch.no_grad():
        for xb, yb in val_loader:
            out = model(xb)
            val_loss += criterion(out, yb).item() * xb.size(0)
    val_losses.append(val_loss / val_size)

    if epoch % 100 == 0 or epoch == 1:
        print(f"Epoch {epoch}/{EPOCHS} | Train MSE: {train_losses[-1]:.5f} | Val MSE: {val_losses[-1]:.5f}")

# ------------------------------
# Save model & training history
# ------------------------------
torch.save(model.state_dict(), os.path.join(MODEL_DIR, "simple_nn.pth"))
np.savetxt(os.path.join(MODEL_DIR, "train_loss.csv"), np.array(train_losses), delimiter=",")
np.savetxt(os.path.join(MODEL_DIR, "val_loss.csv"), np.array(val_losses), delimiter=",")

# ------------------------------
# Plot training history
# ------------------------------
plt.figure(figsize=(6,3))
plt.plot(train_losses, label="Train MSE")
plt.plot(val_losses, label="Validation MSE")
plt.xlabel("Epoch")
plt.ylabel("MSE")
plt.yscale("log")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(MODEL_DIR, "training_curve.png"))
plt.show()
