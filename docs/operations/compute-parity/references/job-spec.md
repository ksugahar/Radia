# Job specification and recovery commands

The local `pack --spec` input uses this structure. Paths in `sources.files`
are Git-tracked paths at the declared commit; keys are destination bundle paths.
Private CAD/mesh/material/reference files belong in `inputs`, not the shared
repository. All source/input paths are local controller paths.

```json
{
  "sources": [{
    "repository": "C:/path/to/reviewed/repository",
    "commit": "full-commit-sha",
    "files": {
      "tools/compute_smoke.py": "tools/compute_smoke.py",
      "validation_test/feec/test_hdiv_vim_demag_solve.py": "validation_test/feec/test_hdiv_vim_demag_solve.py",
      "validation_test/feec/test_hdiv_vim_energy_newton.py": "validation_test/feec/test_hdiv_vim_energy_newton.py"
    }
  }],
  "inputs": {},
  "entry": "tools/compute_smoke.py",
  "args": [],
  "workload": "smoke",
  "minimum_free_memory_bytes": 4294967296,
  "physics_contract": {
    "purpose": "Installed-wheel acceptance, not customer-model validation",
    "units": "SI",
    "route": "Existing public HDiv linear and energy-Newton tests",
    "mesh_material_and_acceptance": "Exactly as committed in the two test modules"
  }
}
```

For a customer solve use `workload: analysis`; bundle a committed public-API
case entrypoint, label-checker contract and all inputs. The physical contract
must name actual settings, not merely copy the smoke example. Set the memory
budget from a small measured case. All arguments are an array, never shell text.

`run` preserves `job.json`, `runtime.json`, `run.log` and `execution.json`; the
case entrypoint must write its own `result.json`. `execution.status: finished`
means the process returned zero and does not itself certify numerical accuracy.

On the compute host after the process exits:

```powershell
<runtime-python> -I <runtime-helper> recover --job-root C:/temp/<owned-job> --output C:/temp/<owned-job>-recovery.zip
```

After transferring that exact archive to durable private controller storage:

```powershell
python tools/compute_runtime.py verify-recovery --archive <durable-archive> --sha256 <remote-reported-hash> --destination <durable-recovered-directory> --report <private-recovery-json>
```

Inspect/commit the intended evidence, then execute the cleanup helper on the
compute host with absolute literal paths and that verified commit:

```powershell
pwsh -NoProfile -File C:/temp/<owned-cleanup-helper>.ps1 -JobRoot C:/temp/<owned-job> -RecoveryArchive C:/temp/<owned-job>-recovery.zip -RecoveredArchiveSHA256 <verified-hash> -EvidenceCommit <full-private-commit> -StagingFiles C:/temp/<owned-job>.zip
```

The helper is staged outside the directory it removes. After recovering its
JSON certificate, remove that helper and any deployment staging using checked
literal-path operations. Do not recursively delete a shared C:\temp directory.
