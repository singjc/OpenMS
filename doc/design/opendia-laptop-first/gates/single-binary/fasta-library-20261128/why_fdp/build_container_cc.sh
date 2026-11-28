#!/usr/bin/env bash
# Container for PLAN_FASTA arms. Same recipe as runtime-round-20261127/hela-evidence: the evshare tools image plus the
# OpenDIA binary, lib/ and the PeptDeep models (build-tree share/OpenMS/models, where the default model lookup looks).
set -euo pipefail
here=/home/singjc/.jcode/scratch/fasta-library-20261128/container_cc
name=opendia-fasta-cc-20261128
repo=/home/singjc/Documents/github/OpenMS-r1-sparse-experiment
build=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build
mkdir -p "$here"
cd "$here"
test ! -e "$name.sif"
test -z "$(git -C "$repo" status --porcelain -- src/topp/OpenDIA)"
git -C "$repo" rev-parse HEAD > source_commit.txt
sha256sum "$build/bin/OpenDIA" > binary.sha256
rm -rf lib models && mkdir -p models
cp "$build/bin/OpenDIA" OpenDIA
cp -a "$build/lib" lib
cp "$build"/share/OpenMS/models/peptdeep_*_dynamic.onnx models/
sha256sum models/* > models.sha256
cat > Dockerfile <<EOF
FROM opendia-evshare-20261126:tools
COPY OpenDIA $build/bin/OpenDIA
COPY lib/ $build/lib/
COPY models/ $build/share/OpenMS/models/
EOF
docker build -t "$name:tools" .
docker run --rm "$name:tools" bash -c 'set -euo pipefail; sha256sum "$(which OpenDIA)"; OpenDIA --help 2>&1 | grep -c -- "-fasta "' > image_check.txt
grep -q "^$(cut -d' ' -f1 binary.sha256) " image_check.txt
singularity build "$name.sif" "docker-daemon://$name:tools" > sif_build.log 2>&1
sha256sum "$name.sif" > container.sha256
singularity exec --cleanenv "$name.sif" bash -c 'sha256sum "$(which OpenDIA)"' > sif_binary.sha256
grep -q "^$(cut -d' ' -f1 binary.sha256) " sif_binary.sha256
rm -rf lib OpenDIA
cat source_commit.txt binary.sha256 container.sha256
