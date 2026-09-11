# Runtime evidence records

Store one concise record per farm/runtime validation gate. Each record must
identify the evaluated commit, setting/gate, fresh artifact paths, relevant
artifact hashes when available, inspection performed, and conclusion.

Do not store source-only test output here as runtime evidence. The bootstrap
scan identified validation-bundle tooling in `testing/`, but no completed
repository-owned farm result record; tooling is not evidence that a run
occurred or passed. Add historical records only from an authoritative handoff
or supplied artifacts.
