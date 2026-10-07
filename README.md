# OT12 All-Cell OT1/OT2 Workflow Summary

## Folder Layout

- `code_snapshot/`: scripts used by the full-cell OT1/OT2 workflow, copied with the original directory structure.
- `result_summaries/`: key json/csv result summaries. Large h5ad/parquet outputs are not copied here.
- `README_en.md`: English workflow summary.
- `MANIFEST.md`: script-by-script description.

## Final Workflow

The final workflow used for reporting is:

1. Full-cell OT1: PCA-expression FGW, low-rank rank=20, alpha=0.9, epsilon=0.
2. Full-cell OT2: molecule-neighborhood hard-count OT2.
3. Full-data second iteration: use the first OT2 log-normalized hard-count output as the next OT1 spatial input, then run OT2 again.

## OT1 Logic

Scripts:

- `code_snapshot/ot1/ot1_fgw_pcaexpr_allcells_lowrank20_eps0_a0p9.py`
- `code_snapshot/ot1/run_ot1_allcells_lowrank20_eps0_a0p9.sbatch`

Main settings:

- spatial cells: 142,272
- scRNA cells: 100,064
- shared genes: 275
- alpha: 0.9
- epsilon: 0
- low-rank rank: 20
- tau_a: 0.8
- tau_b: 0.9
- expression representation for FGW cost: log-normalized expression -> scale(max_value=10) -> PCA(50)
- structure representation: unweighted one-hot label structure
- marginal: type-balanced marginal
- post-FGW expression output: unbalanced blend

Important details:

- PCA coordinates are used only for the FGW expression cost.
- The final barycentric output is written in the original log-normalized gene space, not in PCA space.
- The low-rank coupling is not materialized as a dense `142272 x 100064` transport matrix. The code uses low-rank factors `q/r/g` to compute row mass and barycentric expression, avoiding a dense matrix of about 57GB.

### OT1 Method Details

The purpose of OT1 is to align spatial-cell expression to the scRNA reference while avoiding a complete overwrite of the original spatial signal. The workflow has three layers: solve FGW, compute scRNA barycentric expression, then use unbalanced row mass to blend the scRNA barycentric expression with the original spatial expression.

#### 1. Expression cost: lognorm -> scale -> PCA(50)

Code:

```text
code_snapshot/ot1/ot1_fgw_pcaexpr_allcells_lowrank20_eps0_a0p9.py
build_pca_expression_embedding()
```

Procedure:

1. Keep shared genes and concatenate spatial and scRNA cells.
2. The input expression is already log-normalized gene expression.
3. Apply `sc.pp.scale(max_value=10)` on the joint matrix.
4. Run PCA with 50 PCs.
5. Use the PCA coordinates in `fgw.prepare(joint_attr="X_expr")` for the FGW expression cost.

Rationale:

- Raw gene-space distance is high-dimensional and sensitive to extreme gene values.
- `scale(max_value=10)` limits extreme values.
- PCA creates a more stable low-dimensional continuous expression space.
- This only changes the coordinate system used by the FGW solver. The final output remains in the original gene-expression space.

In short, PCA is the coordinate system for solving OT1, not the final expression output.

#### 2. Structure cost: unweighted one-hot label structure

Code:

```text
make_binary_structural_embedding()
adata.obsm["X_struct"]
fgw.prepare(x_attr="X_struct", y_attr="X_struct")
```

Procedure:

- Take the union of spatial labels and scRNA labels.
- Represent each cell by a one-hot vector for its label.
- The active one-hot entry is `1 / sqrt(2)`.
- No extra class weighting is applied, hence "unweighted" structure.

Why `1 / sqrt(2)`:

- If two cells have the same label, the one-hot squared Euclidean distance is 0.
- If two cells have different labels, the squared Euclidean distance is 1.
- The structure cost is therefore easy to interpret: same label is 0, different label is about 1.

Why keep this term:

- The expression cost captures expression similarity.
- The structure cost adds a label-level constraint so expression alone does not over-mix clearly different cell types.
- alpha=0.9 gives strong weight to this structure term, which was the stable setting in our experiments.

#### 3. Type-balanced marginal

Code:

```text
make_type_balanced_marginal()
fgw.prepare(a="fgw_typebalanced_marginal", b="fgw_typebalanced_marginal")
```

A standard uniform marginal gives each cell the same mass:

```text
cell mass = 1 / N
```

This gives larger cell types more total mass and can make small cell types under-represented.

The type-balanced marginal instead gives each observed label the same total mass:

```text
total mass per observed label = 1 / number of labels
within each label, split that mass uniformly across cells
```

If a domain has K labels and label k has n_k cells, then each cell in that label gets:

```text
a_i = (1 / K) / n_k
```

This is computed independently for spatial and scRNA:

- spatial uses spatial cell-type labels;
- scRNA uses scRNA cell-type labels;
- each side is normalized to total mass 1.

Why this matters:

- Each cell type has the same total opportunity to transport mass.
- Large cell types do not dominate simply because they contain more cells.
- Small cell types are not ignored because of low abundance.

#### 4. Low-rank coupling and row mass

The low-rank FGW solver returns three factors instead of a dense transport matrix:

```text
q: spatial cells x rank
r: scRNA cells x rank
g: rank
```

Mathematically these factors still represent the coupling matrix `P`, but the code does not materialize `P`. It directly computes:

```text
row_mass_i = sum_j P_ij
```

This avoids constructing a dense `142272 x 100064` matrix. A dense float32 version alone would be about 57GB, and temporary arrays can make memory usage even higher. Using `q/r/g` is an equivalent computation, not a change to the OT objective.

#### 5. Barycentric expression

After FGW, each spatial cell gets a scRNA barycentric expression:

```text
X_scrna_fgw_i = sum_j normalized(P_ij) * X_scrna_j
```

Here `normalized(P_ij)` is the row-normalized transport for one spatial cell. It is a weighted average over scRNA reference cells.

#### 6. Unbalanced blend

Code:

```text
unbalanced_scrna_blend_fraction()
lambda_sc = scrna_blend_fraction[:, None]
X_soft_fgw = lambda_sc * X_scrna_fgw + (1 - lambda_sc) * X_sp
```

Because FGW is unbalanced, the actual transported row mass of a spatial cell may be smaller than its source marginal. The code converts this row mass into a blending fraction:

```text
lambda_i = row_mass_i / source_marginal_i
lambda_i = clip(lambda_i, 0, 1)
```

Then:

```text
X_final_i = lambda_i * X_scrna_barycentric_i + (1 - lambda_i) * X_spatial_original_i
```

Interpretation:

- If `row_mass_i` is close to `source_marginal_i`, the cell is well explained by the scRNA reference. `lambda_i` is close to 1, so the output mostly uses scRNA barycentric expression.
- If `row_mass_i` is much smaller than `source_marginal_i`, the unbalanced solver did not fully match this cell. `lambda_i` is smaller, so more original spatial expression is kept.
- This compensates for unbalanced mass loss and avoids forcibly replacing weakly matched spatial cells with scRNA expression.

In the first-round OT1 result, the median blend fraction is about 0.984, so most cells mainly use scRNA barycentric expression while still retaining the unbalanced residual component from the original spatial data.

First-round OT1 output:

```text
OT12_all_exper/results/ot1/allcells_pcaexpr_lowrank20_a0p9_eps0_pca50/
```

Key h5ad:

```text
spatial_patch_soft_barycentric_fgw_labelstruct_unweighted_typebalancedmarginal_unbalancedblend.h5ad
```

## OT2 Logic

The final OT2 implementation is molecule-neighborhood hard-count OT2.

Core scripts:

- `code_snapshot/ot2/run_submit_ot2_hardcount_molneigh_from_existing.sbatch`
- `code_snapshot/ot2/run_prepare_ot2_hardcount_molneigh_from_existing.sbatch`
- `code_snapshot/ot2/submit_ot2_hardcount_molneigh_windows.sh`
- `code_snapshot/ot2/run_ot2_hardcount_molneigh_windows.sbatch`
- `code_snapshot/ot2/ot2_window_hardcount_molneigh_worker.py`
- `code_snapshot/ot2/run_merge_ot2_hardcount_molneigh_allcells.sbatch`
- `code_snapshot/ot2/merge_window_hardcount_assignments.py`
- `code_snapshot/ot2/make_lognorm_from_counts.py`
- `code_snapshot/ot2/summarize_ot2_convergence.py`

OT2 steps:

1. Use the OT1 soft-expression h5ad as the cell-by-gene reference.
2. Use full molecule/transcript tables and cell coordinates.
3. Run OT2 independently across spatial windows.
4. For each molecule/gene, restrict candidate cells to those within `hard_cut_um=12`.
5. Run Sinkhorn OT per window/gene with epsilon=0.05, max_iters=2000, thresh=1e-3.
6. Assign each molecule to the cell with the highest confidence.
7. Merge window outputs. For overlapping molecules, keep the assignment with the highest confidence.
8. Build an integer hard-count cell-by-gene h5ad.
9. Build a log-normalized h5ad from hard counts using `normalize_total(target_sum=1e4) + log1p`.

### OT2 Method Details

The purpose of OT2 is to reassign molecules/transcripts to cells and produce a real count matrix. Unlike OT1, OT2 ends with hard assignments, not a soft expression matrix.

#### 1. Inputs

OT2 uses three inputs:

- molecule/transcript table: `mol_id`, `gene`, `x`, `y`;
- cell coordinate table: cell center coordinates;
- OT1 soft h5ad: cell-by-gene soft expression used as the cell-side target distribution.

For each gene, the code extracts the OT1 soft expression across all cells:

```text
b_full = X_soft[:, gene]
b_full = max(b_full, 0) + EPS
b_full = b_full / sum(b_full)
```

This distribution tells OT2 which cells should be more likely to receive molecules of that gene.

#### 2. Window-level parallelism

The full dataset is too large for one global OT problem, so OT2 is run across spatial windows:

- each Slurm array task processes one window;
- molecules are selected from the window bbox;
- windows may overlap, so the same molecule can be processed more than once;
- the merge step keeps only the highest-confidence assignment for duplicated molecules.

#### 3. Molecule-neighborhood candidate cells

Code:

```text
cell_tree.query_ball_point(mol_xy, r=hard_cut_um)
candidate_rule = "per_gene_union_cells_within_hard_cut_um_of_molecules"
```

For each window and gene:

1. Collect molecules of that gene inside the window.
2. Use a KDTree to find cells within `hard_cut_um=12` um of each molecule.
3. Take the union of those nearby cells as candidate target cells.

OT2 therefore does not allow a molecule to choose from every cell in the tissue. It only solves a local assignment problem. This is both biologically more reasonable and computationally smaller.

If a gene-window has too few molecules or too few candidate cells, the worker skips that subproblem.

#### 4. Distance cost: hinge exponential + hard-cut penalty

Code:

```text
HardCutHingeExpDistanceCost
```

For distance `d`, the normal cost is:

```text
base = (d / r0) ^ p
extra = alpha * (exp(max(d - r0, 0) / sigma) - 1)
normal_cost = base + extra
```

Parameters:

```text
r0 = 15
p = 2
sigma = 5
alpha = 10
hard_cut_um = 12
big_penalty = 1e6
```

Final cost:

```text
if d <= hard_cut_um:
    cost = normal_cost
else:
    cost = big_penalty + normal_cost
```

Interpretation:

- Close molecules/cells receive a smooth distance penalty.
- Beyond `r0`, an exponential penalty grows quickly.
- Beyond `hard_cut_um=12`, a large penalty of `1e6` is added, effectively forbidding long-distance assignments.
- Candidate construction already restricts cells to 12 um; the hard-cut cost is an additional safeguard.

#### 5. Per-gene/per-window OT problem

For one gene in one window:

- source points are molecules, count M;
- target points are candidate cells, count N;
- source marginal:

```text
a_m = 1 / M
```

- target marginal from OT1 soft expression:

```text
b_c = X_soft[c, gene] / sum_candidate_cells X_soft[c, gene]
```

Then the code runs unbalanced Sinkhorn:

```text
tau_a = 0.8
tau_b = 0.9
epsilon = 0.05
max_iters = 2000
thresh = 1e-3
```

Intuition:

- molecule-side mass encourages every molecule to be assigned;
- cell-side mass follows the OT1 soft expression profile;
- unbalanced OT allows local molecule counts and predicted soft expression mass to disagree.

#### 6. From soft OT matrix to hard assignment

Sinkhorn returns a transport matrix `P`:

```text
M molecules x N candidate cells
```

The code normalizes each molecule row:

```text
P_row = P / row_sum(P)
```

Then each molecule is assigned to the highest-probability cell:

```text
assigned_cell = argmax_j P_row[molecule, j]
assign_conf = max_j P_row[molecule, j]
```

Thus OT2 produces one final cell assignment per molecule. `assign_conf` is the confidence of that hard assignment.

In other words, each RNA/molecule contributes to exactly one cell:

```text
count[assigned_cell, gene] += 1
```

The molecule is not split across cells according to the soft OT probabilities. The row-normalized OT matrix is used only to choose the highest-probability cell, and then that molecule is added as one hard count to that cell-gene entry.

#### 7. Merge overlapping windows

Code:

```text
merge_window_hardcount_assignments.py
```

After all windows are combined, the same molecule may appear in multiple windows. The merge rule is:

```text
group by mol_id and keep the row with the highest assign_conf
```

This ensures that each molecule is counted exactly once.

#### 8. How the hard-count cell-by-gene h5ad is built

After merge, each molecule has:

```text
mol_id, gene, assigned_cell
```

The code counts molecules per cell and gene:

```text
count(cell, gene) = number of molecules assigned to that cell with that gene
```

Then it builds a sparse matrix:

```text
X[cell_index, gene_index] = count
```

The AnnData object is:

```text
adata.X = integer count matrix
adata.obs = obs from OT1 soft h5ad
adata.var = var from OT1 soft h5ad
adata.obsm["spatial"] = spatial coordinates from OT1 soft h5ad
```

This is the `*_cell_by_gene.h5ad` file. Its `X` is molecule count, not log-normalized expression.

Therefore the order is:

```text
each molecule -> argmax cell -> add 1 count to count(cell, gene)
```

The full hard-count matrix is created first. Normalization happens only afterward.

#### 9. Purpose of the lognorm h5ad

From the hard-count h5ad, the pipeline also creates:

```text
*_lognorm_cell_by_gene.h5ad
```

Procedure:

```text
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
```

Thus normalization is not applied before molecule assignment. It is applied after all molecules have been hard assigned and accumulated into a cell-by-gene count matrix.

Uses:

- UMAP and silhouette metrics are computed on log-normalized expression.
- Round-2 OT1 reads this lognorm h5ad using `existing_lognorm` mode.
- Raw molecule counts remain available in `.layers["counts"]`.

First-round OT2 output:

```text
OT12_all_exper/results/ot2/hardcount_molneigh_allcells_from_allcells_pcaexpr_lowrank20_a0p9_eps0_pca50/
```

Key outputs:

```text
ot2_allcells_hardcount_molneigh_cell_by_gene.h5ad
ot2_allcells_hardcount_molneigh_lognorm_cell_by_gene.h5ad
ot2_allcells_hardcount_molneigh_molecule_assignments.parquet
```

## Second Iteration

Scripts:

- `code_snapshot/iter_round2/run_submit_iter_round2_pipeline.sbatch`
- `code_snapshot/iter_round2/ot1/ot1_fgw_pcaexpr_iter2_from_ot2.py`
- `code_snapshot/iter_round2/ot1/run_ot1_iter2_from_ot2_molneigh.sbatch`
- `code_snapshot/iter_round2/ot2/run_prepare_iter2_ot2_inputs.sbatch`
- `code_snapshot/iter_round2/ot2/run_ot2_iter2_windows.sbatch`
- `code_snapshot/iter_round2/ot2/run_merge_iter2_ot2.sbatch`
- `code_snapshot/iter_round2/eval/run_umap_ot1_iter2.sbatch`
- `code_snapshot/iter_round2/eval/run_umap_ot2_iter2.sbatch`
- `code_snapshot/iter_round2/eval/collect_iter_round2_summary.py`

Second-round logic:

1. Read the first OT2 log-normalized hard-count h5ad:
   `ot2_allcells_hardcount_molneigh_lognorm_cell_by_gene.h5ad`
2. Use it as the spatial expression input for the second OT1.
3. Use `existing_lognorm` mode. No additional gene-sum mass matching is applied.
4. Keep OT1 parameters unchanged: alpha=0.9, epsilon=0, rank=20, PCA(50) expression cost.
5. Run the same molecule-neighborhood OT2 from the second OT1 soft h5ad.
6. Save UMAP plots, metrics, molecule assignment summary, and Sinkhorn convergence summary.

Second-round outputs:

```text
OT12_all_exper/iter_round2/results/
```

## Run Commands

First-round OT1:

```bash
cd /scratch/beul/jzhang79/projects/Moscot/code
sbatch OT12_all_exper/ot1/run_ot1_allcells_lowrank20_eps0_a0p9.sbatch
```

First-round OT1 UMAP:

```bash
sbatch OT12_all_exper/eval/run_umap_ot1_allcells.sbatch
```

First-round molecule-neighborhood OT2:

```bash
sbatch OT12_all_exper/ot2/run_submit_ot2_hardcount_molneigh_from_existing.sbatch 80
```

Second full-data iteration:

```bash
sbatch OT12_all_exper/iter_round2/run_submit_iter_round2_pipeline.sbatch 80
```

`80` is the maximum concurrent Slurm array throttle. If the Slurm controller is unstable, use a smaller value such as `40` or `20`.

## Key Results

### OT1

| stage | coupling sum | gene var ratio | cell total std ratio | median blend |
|---|---:|---:|---:|---:|
| iter1 OT1 | 0.980987 | 0.277457 | 0.607772 | 0.983954 |
| iter2 OT1 | 0.982815 | 0.261178 | 0.576796 | 0.985811 |

### OT2 Assignment

| stage | molecule loss | final molecules | cells used | genes |
|---|---:|---:|---:|---:|
| iter1 OT2 | 1.7997% | 15,361,526 | 142,053 | 275 |
| iter2 OT2 | 1.7997% | 15,361,526 | 141,995 | 275 |

### UMAP / Expression-Space Metrics

| stage | separate UMAP silhouette | separate PCA silhouette | ASIS lognorm | scaled lognorm |
|---|---:|---:|---:|---:|
| iter1 OT1 | -0.067 | 0.395 | 0.379 | 0.375 |
| iter2 OT1 | -0.123 | 0.501 | 0.475 | 0.469 |
| iter1 OT2 | 0.178 | 0.122 | 0.035 | 0.026 |
| iter2 OT2 | 0.166 | 0.106 | 0.029 | 0.021 |

Interpretation:

- The second OT1 strengthens separation in PCA/expression space.
- After the second OT2 hard-count step, the final OT2 scores do not improve over the first OT2 and are slightly lower.
- The more stable final result is therefore the first OT1 -> first molecule-neighborhood OT2 result.

### OT2 Sinkhorn Convergence

Second OT2:

- gene-window problems: 526,623
- converged rate: 100%
- reached max_iters: 0
- median iterations: 20
- mean iterations: 22.29
- p99 iterations: 30
- max iterations: 30

Although `max_iters=2000`, all OT2 subproblems converged within 30 Sinkhorn iterations.

## Slurm Notes

OT2 launches many window-array jobs. We observed Slurm-side instability:

```text
Slurm temporarily unable to accept job
Resource temporarily unavailable
DependencyNeverSatisfied
```

This is usually not a code-logic failure. If a few array windows are missing, the merge job may remain blocked because the `afterok` dependency is not satisfied. The repair strategy is:

1. Detect missing windows.
2. Re-run the same `run_ot2*_windows.sbatch` only for missing windows.
3. Re-submit merge, UMAP, and summary jobs.

`code_snapshot/ot2/run_repair_ot2_hardcount_molneigh_missing_then_merge.sbatch` is kept for this repair case.
