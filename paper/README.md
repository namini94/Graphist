# Paper materials

- `GRAPHIST_RECOMB2025.pdf` — the paper draft.
- `working_deck/Total.pptx` — the working slide deck; each finalized figure below is a curated subset of
  its slides.
- `Fig2/` — PDAC main results (region annotation, DDX60L/MYEOV knockdown phenotype-selected spots, region
  composition bar charts, marker genes, pathway activity spatial maps + boxplot). Source: `Total.pptx`
  slides 3–4.
- `Fig4/` — BRCA-PACSI main results (region annotation, survival-phenotype selected spots, histology
  morphology clustering, immune-subset pathway activity, ligand-receptor co-expression). Source:
  `Total.pptx` slides 6–7.
- `SFig1/` — PDAC MYEOV-knockdown supplementary results (volcano plot + marker genes). Source: `Total.pptx`
  slide 5.
- `misc/PhenoGraphAgent_Nami.pptx` — related slide deck (multi-agent phenotype-discovery framework, see
  `[20]` in the paper's references).

## Figure → code mapping

Most panels trace back to scripts in `downstream_analysis/` and `stage2_pathway_vgae/`:

| Panel type | Script |
|---|---|
| Region annotation / Graphist(+/-) spatial plots | `downstream_analysis/preprocessing/{PDAC,BRCA-PACSI}-preprocess.py` |
| Region-composition bar/stacked-bar charts | `downstream_analysis/figures/{Barplot,Barplot_PACSI,Stacked_Barplot,Stacked_Barplot_Compare}.py` |
| Volcano plot (gene-level DE) | `downstream_analysis/figures/Volcano.py` |
| Pathway lollipop plots | `downstream_analysis/figures/Lollipop_PDAC.py` |
| Pathway spatial maps + "Pathway Activity Across Annotations" boxplot | `stage2_pathway_vgae/legacy/VGAE-PA-ST-{PDAC,PACSI}.py` (downstream section) |
| Image-feature (ResNet50) clustering | `downstream_analysis/morphology/PACSI-MORPH-V2.py` |
| Ligand-receptor dotplot / co-expression | `downstream_analysis/preprocessing/BRCA-PACSI-preprocess.py` |

The bulk-sample clustering heatmaps and survival scatter plots (Fig2 panel B, Fig4 panel B) come from the
Stage 1 R code in `stage1_bulk_regression_R/`, not from any Python script here.
