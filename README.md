# Code to: A transcriptional map <br> of human tonsil architecture

## setup

- R version and library have are specified in the `config.yaml` file  
  (e.g., `R: "R_LIBS_USER=/path/to/library /path/to/R/executable"`)
- `.Rprofile` is used for handling command line arguments
- `logs/` capture `.Rout` files from `R CMD BATCH` executions
- intermediate results are written to `outs/` (as *.rds*)
- visualizations are written to `plts/` (as *.pdf*)

## steps

- `01-raw.R`
  - read flat files as `SingleCellExperiment`
  - stash non-RNA targets as `altExps`
  - write out as *.h5*-backed object
  
- `02-fil.R`
  - exclude cells too close to any FOV border
  - exclude cells with low counts, low counts per area,  
  and high negative probe or false code count
  
- `03-ist.R`: semi-supervised `InSituType` clustering  
using `data/raf/mtx.rds` as reference profiles
  
- `04-sub.R`:  
    - subset cells into epithelia (epi), stromal (str),  
    myeloid (mye), B cells (bcs), T cells (tcs)  
    - subpopulations to pool for each subset  
    are specified in `meta/lab/sub.json`
  
- `05-jst.R`:  
`InSituType` subclustering (analogous to above),  
using distinct references profiles for each subset
  
- `00-lab.R`: used to relabel  
`03-ist.R` outputs (automated; using `meta/lab/lv1.json`), and  
`05-jst.R` outputs (manual; using `meta/lab/lv2,<sub>.json`)
  
- downstream analyses
    - `03-pro.R` and `05-rep.R`: PCA using feature subset  
    according to `03-ist.R` and `05-jst.R`, respectively
    - `04-ccc.R`: cell cell-communication
    - `06-ctx.R`: spatial contexts/niches

- `10-plt__<by1>__<by2>-<out1>,<out2>,<plt>.R`
  - pool outputs `<out1/2>` according to `<by1/2>`  
  - generates `plts/<out1>,<out2>,<plt>.pdf`  
  (depending on `<by1/2>`, name may also include  
  subset (`sub`) or section (`sid`) identifier)
    - `sid` = one section
    - `all_sid` = all sections
    - `all_sid_all_sub` = all sections and subsets
    