#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=160G
#SBATCH -t 2:00:00
cd /nfs/research/juan/singj/opendia/fasta-20261128
ls -la cache_v4/Pep/ | head; ls survey_v3/V0-D/
uv run -q --with pyarrow --with numpy python diag_pool.py diag_pool.json
