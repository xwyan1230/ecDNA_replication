import muon as mu
import numpy as np
from muon import MuData
from anndata import AnnData
import tempfile

np.random.seed(1)

n, d, k = 1000, 100, 10

z = np.random.normal(loc=np.arange(k), scale=np.arange(k)*2, size=(n, k))
w = np.random.normal(size=(d, k))
y = np.dot(z, w.T)

adata = AnnData(y)
adata.obs_names = [f'obs_{i+1}' for i in range(n)]
adata.var_names = [f'var_{j+1}' for j in range(d)]

d2 = 50
w2 = np.random.normal(size=(d2, k))
y2 = np.dot(z, w2.T)

adata2 = AnnData(y2)
adata2.obs_names = [f"obs_{i+1}" for i in range(n)]
adata2.var_names = [f"var2_{j+1}" for j in range(d2)]
adata2

mdata = MuData({"A": adata, "B": adata2})
print(np.sum(mdata.obsm['A']))
print(mdata.varm['A'])

# Only keep variables with value > 1 in obs_1
# with in-place filtering for the variables
mu.pp.filter_var(adata, adata["obs_1",:].X.flatten() > 1)
print(adata)
# Modalities can be accessed within the .mod attributes
print(mdata.mod["A"])

print(f"Outdated size:", mdata.varm["A"].sum())
mdata.update()
print(f"Updated size:", mdata.varm["A"].sum())

# Throw away the last sample in the modality 'B'
# with in-place filtering for the observations
mu.pp.filter_obs(mdata.mod["B"], [True for _ in range(n-1)] + [False])

# adata2 object has also changed
assert mdata.mod["B"].shape == adata2.shape
mdata.update()
print(mdata)

# drop the observations that are not present in all the modalities
mu.pp.intersect_obs(mdata)
print(mdata)

# Create a temporary file
temp_file = tempfile.NamedTemporaryFile(mode="w", suffix=".h5mu", prefix="muon_getting_started_")

mdata.write(temp_file.name)
mdata_r = mu.read(temp_file.name, backed=True)
mdata_r