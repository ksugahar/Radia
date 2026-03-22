Synchronize knowledge across all documentation layers after code changes or new discoveries.

## Steps

1. **Detect changes**: Check recent git commits (`git log --oneline -10`) and unstaged changes for new features, bug fixes, API changes, or research findings.

2. **CLAUDE.md**: Update policies, architecture, API references. Add new sections for new features. Remove references to deleted code.

3. **MCP knowledge** (`src/radia/mcp/`):
   - `radia/radia_knowledge.py`: Radia API, conventions, best practices
   - `ngsolve/ngsolve_knowledge.py`: NGSolve formulations, solver tips
   - `ngsolve/sparsesolv_knowledge.py`: Compact AMS, COCR, IC usage
   - `cubit/export_knowledge.py`: Cubit mesh export workflows
   - `cubit/rules.py`: Lint rules for deleted/deprecated APIs

4. **docs/research/**: Add mathematical derivations, formulation documents for new theoretical results.

5. **memory/**: Update project memories for decisions, discoveries, or status changes that affect future conversations.

6. **Verify**: Run `python -m pytest tests/ -q --ignore=tests/cubit --ignore=tests/mcp` to ensure nothing broke.

7. **Commit**: Stage all documentation changes and commit with descriptive message.

## Guidelines

- Only update files where changes are actually needed (don't touch unchanged areas)
- MCP knowledge updates should include code examples that users can copy
- CLAUDE.md updates should be policies or architecture, not tutorials
- docs/research/ is for mathematical derivations and theoretical background
- memory/ is for decisions and context that help in future conversations
- Always verify code examples in MCP knowledge are syntactically correct
