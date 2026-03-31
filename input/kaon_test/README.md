Use this directory for small test manifests that should not be auto-discovered
from the production `input/kaon` tree.

Typical workflow:

1. Copy one or two production manifests into `input/kaon_test/`.
2. Point their `runs_file` entries at short test run lists.
3. Run farm helpers with `-m input/kaon_test`.
4. Run Jasmine with `--manifest-dir input/kaon_test`.

Examples:

```bash
./run_farm.sh -m input/kaon_test 3p0 3p14
./run_farm.sh -c -m input/kaon_test 3p0 3p14
python farm_env/jasmine_put_from_manifest.py center high 3p0 3p14 lh2 --manifest-dir input/kaon_test --product-kind replay
```
