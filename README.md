# scRNASeqAnalysisAndModelling
This is a repository containing all my efforts to learn tools that are related to the analysis and use of single-cell RNA-seq (scRNA-seq) data. In particular, we are focussed on tools that can be used to infer cell-cell communication pathways Packages used and their purposes include:
	- [Seurat](https://satijalab.org/seurat/), for the analysis and clustering of scRNA-seq data.
	- [Scanpy](https://scanpy.readthedocs.io/en/stable/), for the analysis and clustering of scRNA-seq data.
    - [scVelo](https://scvelo.readthedocs.io/scVelo), for inferring cell state trajectories via dynamical models of RNA velocity.
	- [py-pde](https://py-pde.readthedocs.io/en/latest/), for the numerical solution of reaction-diffusion PDE models that
    we use to model ligand-receptor binding across spatial domains.
	- [CellChat](https://github.com/sqjin/CellChat), a package by Suoqin Jin (fellow postdoc) that infers cell-cell communication by calculating ligand-receptor probabilities and communication patterns from scRNA-seq data.
	- [NicheNet](https://github.com/saeyslab/nichenetr), a package that infers interactions between ligands and their target genes, which occur as a downstream response to the ligand-receptor bindings.
