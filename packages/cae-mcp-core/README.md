# CAE MCP Core

Shared Python support for independently distributed CAE MCP servers. It owns
tool registration/contracts, status and source provenance, editable reload,
failure logs, resource/search helpers, and portable mesh artifact inspection.
It does not depend on Radia, radia-mcp, cubit-mesh-export, native solver binaries
or a Coreform Cubit installation. Optional libraries are imported only by the
operations that require them.

`cubit-mesh-export` owns its own MCP server and depends on this foundation.
`radia-mcp` also depends on this foundation; it does not own the Cubit server.
Radia topology optimization may still use Cubit as a mesh-generation backend.
Independent distribution does not eliminate that explicitly selected workflow
dependency or its compatibility validation.

This package is extracted from the existing BSD-3-Clause implementation, not a
second framework. Old `radia_mcp.common` and reload forwarding modules are not
shipped. Consumers import `cae_mcp_core` directly. Release this foundation
before publishing packages that declare its version as a dependency.
