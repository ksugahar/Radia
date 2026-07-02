# loop_learning validation

This directory contains public-safe validation utilities for CAE-AI Lab loop
learning.  The utilities take queue/result paths as command-line arguments; the
repository does not hard-code private source-tool paths.

Run an autonomous basic-learning pass:

```powershell
python packages/radia-mcp/validation/loop_learning/autonomous_basic_learning.py --queue-json <queue.json> --out-dir <out-dir>
```
