# branchiostoma-orthofinder-pipeline
OrthoFinder pipeline for ortholog inference in Branchiostoma lanceolatum

Pipeline for assigning human gene names to predicted proteins of *Branchiostoma lanceolatum* (assembly klBraLanc5.hap2) via OrthoFinder ortholog inference across five deuterostome proteomes.

## Species used
- *Branchiostoma lanceolatum* GCF_035083965.1 (query)
- *Homo sapiens* GCF_000001405.40
- *Mus musculus* GCF_000001635.27
- *Danio rerio* GCF_000002035.6
- *Strongylocentrotus purpuratus* GCF_000002235.5

## Results
- 24,423 protein-coding genes analyzed
- 11,961 genes (49.0%) assigned a human ortholog gene name

## Dependencies
```bash
conda create -n orthofinder_env -c bioconda -c conda-forge \
    orthofinder python=3.10 ncbi-datasets-cli pandas biopython -y
conda activate orthofinder_env
```

## Usage
```bash
python ortholog_pipeline_zenodo.py --step all
```

## How to cite this pipeline
If you use this pipeline in your research, please cite:
- This repository: https://github.com/bozzom9/branchiostoma-orthofinder-pipeline

## Citation
If you use this pipeline, please cite:
- Emms D.M. et al. (2025). OrthoFinder. bioRxiv https://doi.org/10.1101/2025.07.15.664860
- Emms D.M. & Kelly S. (2019). Genome Biology 20:238.
- Buchfink B. et al. (2015). Nature Methods 12:59-60.

## Extended annotation pipeline (v1.1)

In addition to OrthoFinder ortholog inference, the pipeline now includes
a Swiss-Prot homology search for genes without a human ortholog.

### Additional scripts

- `build_NP_mapping.py` — builds NP_ accession → gene symbol mapping from NCBI files
- `annotate_unannotated.py` — runs DIAMOND blastp against Swiss-Prot for unannotated genes
- `merge_annotations_v2.py` — merges OrthoFinder and Swiss-Prot annotations into a single file
- `add_gene_names_v2.py` — adds gene_name and annotation_source columns to count matrices

### Annotation results (B. lanceolatum klBraLanc5.hap2)

| Source | Genes | % |
|--------|-------|---|
| OrthoFinder (human ortholog) | 11,961 | 49.0% |
| Swiss-Prot homolog (metazoa) | 5,717 | 23.4% |
| Unknown | 6,745 | 27.6% |
| **Total** | **24,423** | **100%** |

### annotation_source column
Output files include an `annotation_source` column indicating the origin:
- `OrthoFinder` — formal human ortholog inferred by OrthoFinder
- `SwissProt_homolog` — best homolog from metazoa via DIAMOND vs Swiss-Prot
- `unknown` — no annotation found
