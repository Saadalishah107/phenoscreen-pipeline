#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
docker build -t phenoscreen/genomics:2.0 -f docker/Dockerfile.genomics docker
docker build -t phenoscreen/structure:2.0 -f docker/Dockerfile.structure docker
docker build -t phenoscreen/cheminformatics-test:2.0 -f docker/Dockerfile.cheminformatics-test docker
if [[ "${1:-}" == "--full" ]]; then
  docker build -t phenoscreen/cheminformatics:2.0 -f docker/Dockerfile.cheminformatics docker
fi
echo "Images built."
