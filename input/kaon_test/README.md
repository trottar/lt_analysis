Use this directory for small test manifests that should not be auto-discovered
from the production `input/kaon` tree.

Typical workflow:

1. Copy one or two production manifests into `input/kaon_test/`.
2. Point their `runs_file` entries at short test run lists.
3. Run farm helpers with `-m input/kaon_test`.
4. After the SWIF workflow succeeds, run Jasmine from an interactive ifarm session with `-j -m input/kaon_test`.

Examples:

```bash
./run_farm.sh -m input/kaon_test 3p0 3p14
./run_farm.sh -c -m input/kaon_test 3p0 3p14
./run_farm.sh -j -m input/kaon_test -s 3p0 3p14
./run_farm.sh -j -c -m input/kaon_test -s 3p0 3p14
```
