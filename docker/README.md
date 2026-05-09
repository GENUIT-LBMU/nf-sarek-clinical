# Docker images for nf-sarek-clinical

Three custom images supplement the public biocontainers used by the pipeline:

| Image | Purpose | Source |
|---|---|---|
| `civicpy:5.2.0` | CIViC variant annotation + Python data parsing (pandas/numpy) | `Dockerfile.civicpy` |
| `r-annotation:1.0.0` | R + GenomicRanges + tidyverse for gene-coverage / MANE annotation | `Dockerfile.r_annotation` |
| `sigprofiler:0.1.7` | SigProfilerAssignment for mutational signature decomposition (COSMIC v3.4) | `Dockerfile.sigprofiler` |

The other public images (snpsift, bcftools, GATK, HRProfiler) are pulled directly from quay.io / Docker Hub by Nextflow at run time.

## Build locally (default)

```bash
bash docker/build_local.sh
```

This tags the images as `nf-sarek-clinical/<name>:<version>`, which matches the pipeline default (`params.docker_namespace = 'nf-sarek-clinical'`). Just run the pipeline afterwards and Nextflow will pick them up.

## Build under a custom namespace

```bash
NAMESPACE=ghcr.io/myorg bash docker/build_local.sh
docker push ghcr.io/myorg/civicpy:5.2.0
docker push ghcr.io/myorg/r-annotation:1.0.0
docker push ghcr.io/myorg/sigprofiler:0.1.7
```

Then at run time:

```bash
nextflow run ... --docker_namespace ghcr.io/myorg
```

## Use a pre-published registry

If the GENUIT-LBMU lab (or another fork) has published the images:

```bash
nextflow run ... --docker_namespace ghcr.io/GENUIT-LBMU/nf-sarek-clinical
```

No local build needed — Nextflow will pull the images on first use.

## Image notes

- **`Dockerfile.civicpy`** is also used by the `withLabel:civicpy` (CIViC annotation) **and** by `SOMATIC_CNV_INTERPRET` (pandas-based CNV parsing). The image carries pandas as a transitive dependency of civicpy, so a single image serves both labels.
- **`Dockerfile.r_annotation`** uses a biocontainer mulled-v2 base. For a slimmer image, consider migrating to `rocker/r-base` + manual `install.packages` (future work).
- **`Dockerfile.sigprofiler`** currently pins Python 3.10. SigProfilerExtractor 1.1.21+ supports Python 3.11 — unifying to 3.11 is a TODO (see Phase C in the audit).
