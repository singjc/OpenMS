# OpenDIA module layout, 27 November 2026

Implementation commit: `1aeda81e11`, branch `feature/opendia-clean-rewrite` on `singjc/OpenMS`.

## Scope

The private TOPP implementation keeps separate `include/` and `source/` trees. Both now group files under `workflow/`, `library/`, `io/`, `calibration/`, `extraction/`, `features/` and `scoring/`. Redundant `OpenDIA` filename prefixes are removed. Qualified includes, for example `<extraction/Traces.h>`, resolve from the private include root.

The executable entry point remains `OpenDIA.cpp`. The namespace stays `OpenMS::OpenDIA`. The implementation stays one static target, without per-module CMake libraries. Tests stay flat in `test/`, with shorter filenames, while executable target names and CTest names retain their OpenDIA prefixes. Source order, compiler flags and link dependencies are unchanged.

There are no algorithm, numeric type, option or output-schema changes. The source README describes the new layout. Historical experiment documents retain their original source-path references as provenance.

## Validation

- Reconfigured and rebuilt OpenDIA and all 14 core test executables successfully in the existing full OpenMS Ninja build.
- All 20 OpenDIA CTests passed, including TOPP INI/CTD and CLI rejection checks.
- Moved C++ headers, sources and test bodies were checked against the parent revision: only include directives differ.
- The no-OpenSWATH architectural guard now scans recursively through module folders and fails on an empty scan. Isolated negative fixtures confirmed that a forbidden nested include and an absent source tree are both rejected.
- Representative Astral raw01 seed0 run passed: all three output files are SHA256-identical to the read-only pre-refactor `out_base`, the file set is identical, and the spectrum cache directory was removed. Exit code 0. The loaded-laptop run took 10:03.48 with peak RSS 8,105,328 KiB. This is an identity check, not a timing comparison. Logs, binary hash, output hashes and GNU time evidence are in [evidence/](evidence/).

Exception diagnostics naturally report the new source filenames through `__FILE__`. Successful analysis output must remain byte-identical. This is a source-layout change, not a new performance experiment.
