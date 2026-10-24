# Docker images for nf-sarek-clinical

Three custom images supplement the public biocontainers used by the pipeline:

| Image | Purpose | Source |
|---|---|---|
| `civicpy:5.4.0` | CIViC variant annotation + Python data parsing (pandas/numpy) | `Dockerfile.civicpy` |
| `r-annotation:1.0.0` | R + GenomicRanges + tidyverse for gene-coverage / MANE annotation | `Dockerfile.r_annotation` |
| `sigprofiler:0.1.8` | SigProfilerAssignment for mutational signature decomposition (COSMIC v3.4) | `Dockerfile.sigprofiler` |

The other public images (snpsift, bcftools, GATK, HRProfiler) are pulled directly from quay.io / Docker Hub by Nextflow at run time.

## Build locally (default)

```bash
bash docker/build_local.sh
```

This tags the images as `nf-sarek-clinical/<name>:<version>`, which matches the pipeline default (`params.docker_namespace = 'nf-sarek-clinical'`). Just run the pipeline afterwards and Nextflow will pick them up.

## Build under a custom namespace

```bash
NAMESPACE=ghcr.io/myorg bash docker/build_local.sh
docker push ghcr.io/myorg/civicpy:5.4.0
docker push ghcr.io/myorg/r-annotation:1.0.0
docker push ghcr.io/myorg/sigprofiler:0.1.8
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

## Camino B — mirror all 7 images, digest-pinned (clinical reproducibility)

For air-gapped or accredited-diagnostic environments, mirror **every** image
(the 3 custom + the 4 public) into the org registry and pin each by its
immutable `sha256` digest. This removes any dependency on Docker Hub / quay.io
at run time and guarantees byte-for-byte identical images forever.

```bash
# 1. Authenticate to the destination registry (write access required)
echo "$GHCR_TOKEN" | docker login ghcr.io -u <user> --password-stdin

# 2. Build+push the 3 custom images, copy the 4 public ones, and pin all 7.
#    Regenerates conf/containers.mirror.config with the captured digests.
MIRROR_REGISTRY=ghcr.io/genuit-lbmu bash docker/mirror_images.sh

# 3. Review the regenerated conf/containers.mirror.config, then commit it.

# 4. Run against the pinned mirror (works with docker OR singularity):
nextflow run ... -profile docker,mirror
nextflow run ... -profile singularity,mirror
```

The `mirror` profile (in `nextflow.config`) overrides the container assignments
with `conf/containers.mirror.config`. Without `,mirror`, the pipeline keeps
using the upstream public tags (`conf/containers.config`).

`docker/mirror_images.sh` uses `docker buildx imagetools` (no `skopeo` needed):
public images are copied registry-to-registry without a local pull, preserving
their manifests. To update a tool, bump its tag in the script, re-run it, and
commit the new digest.

> The committed `conf/containers.mirror.config` ships with **tag** references as
> a starting point so `-profile mirror` works right after a push — but tags are
> mutable. Run the script to replace them with `@sha256:…` digests, which is the
> actual Camino B guarantee.

The custom images can also be built+pushed automatically by the
`.github/workflows/build-images.yml` Action (on a `v*` tag or manual dispatch);
you then only need the script to mirror the 4 public images and pin digests.

## Reproducible builds (dependency pinning)

All three custom images build deterministically:

- **Base images are pinned by digest** (`FROM python@sha256:…`, `FROM …mulled…@sha256:…`), not by mutable tags. Refresh a digest with
  `docker buildx imagetools inspect <ref> --format '{{.Manifest.Digest}}'`.
- **Python deps are fully locked** in `requirements.civicpy.txt` and
  `requirements.sigprofiler.txt` (complete transitive closure with `==` pins).
  Each file's header has the one-liner to regenerate it.
- **R/CRAN deps** resolve against a **dated Posit Package Manager snapshot**
  (`packagemanager.posit.co/cran/<date>` in `Dockerfile.r_annotation`), so
  `install.packages()` is frozen instead of tracking CRAN HEAD. Bump the date to update.

This is what makes the Camino B `@sha256` digests meaningful: the frozen digest
corresponds to a build you can actually recreate bit-for-bit.

## Image notes

- **`Dockerfile.civicpy`** is also used by the `withLabel:civicpy` (CIViC annotation) **and** by `SOMATIC_CNV_INTERPRET` (pandas-based CNV parsing). The image carries pandas as a transitive dependency of civicpy, so a single image serves both labels.
- **`Dockerfile.r_annotation`** uses a biocontainer mulled-v2 base. For a slimmer image, consider migrating to `rocker/r-base` + manual `install.packages` (future work).
- **`Dockerfile.sigprofiler`** currently pins Python 3.10. SigProfilerExtractor 1.1.21+ supports Python 3.11 — unifying to 3.11 is a TODO (see Phase C in the audit).
- **apt packages** (bcftools, tabix, compilers) are not version-pinned; the base-image digest fixes the Debian release, but `apt-get` still tracks the live mirror. Full apt determinism would need `snapshot.debian.org` — residual gap, low impact for these tools.
