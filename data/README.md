# Data

This folder is not tracked by git. Put the Cell Ranger `per_sample_outs`
directory here, or point `paths.cellranger_root` in `config/config.yml`
somewhere else (an absolute path is fine).

The pipeline expects this layout:

```
data/per_sample_outs/
├── MCC8209A3_Male_BM_3/
│   └── sample_filtered_feature_bc_matrix/
└── MCC8209A4_Female_BM_4/
    └── sample_filtered_feature_bc_matrix/
```

Only the filtered (cell-containing) matrix is read. The raw matrix is not
used: it was only needed for ambient-RNA estimation, which the pipeline no
longer does.

A symlink works just as well as a copy:

```bash
ln -s "/path/to/Aging Mouse Data/CITE-seq/aging_BM/per_sample_outs" data/per_sample_outs
```
