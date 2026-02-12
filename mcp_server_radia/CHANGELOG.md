# Changelog

All notable changes to the Radia MCP Server will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0] - 2026-02-12

### Added
- **Diagnostic and Debugging Tools** (4 new tools)
  - `radia_server_info`: Get server version and status information
  - `radia_list_objects`: List all Radia objects in server state
  - `radia_get_object_info`: Get detailed information about specific objects
  - `radia_clear_state`: Reset server state (clear all objects)

- **NGSolve Integration Tools** (5 new tools)
  - `radia_ngsolve_create_field`: Create field data for NGSolve CoefficientFunction
  - `radia_ngsolve_batch_evaluate`: Batch field evaluation with H-matrix acceleration
  - `radia_ngsolve_enable_hmatrix`: Enable H-matrix acceleration for field evaluation
  - `radia_ngsolve_export_field_grid`: Export structured field grid to NPZ format
  - `radia_ngsolve_mesh_field_sample`: Sample field at NGSolve mesh vertices

- **Documentation Enhancements**
  - Development policy and branch management guidelines
  - Detailed usage examples for NGSolve integration
  - H-matrix acceleration performance characteristics (10-100x speedup)
  - Complete workflow examples

### Improved
- Server architecture with modular tool organization
- Error handling with detailed traceback information
- Logging infrastructure for debugging and monitoring

### Technical Details
- Total tools: 35 (up from 26)
- H-matrix acceleration support for large-scale field evaluation
- NPZ field grid export for scipy.interpolate integration
- State management and inspection capabilities

## [1.0.0] - 2026-02-11

### Added
- Initial release of Radia MCP Server
- **Geometry Tools** (6 tools)
  - Rectangular magnets, hexahedrons, tetrahedrons
  - Race-track coils, containers
  - Unit system management

- **Material Tools** (6 tools)
  - Linear and anisotropic materials
  - Nonlinear B-H curve materials
  - Permanent magnets (fixed and linear)

- **Solver Tools** (5 tools)
  - LU, BiCGSTAB, and HACApK solvers
  - H-matrix parameter configuration
  - Background field support

- **Field Calculation Tools** (5 tools)
  - Point, line, and grid field evaluation
  - Force and torque calculations

- **Export and Workspace Tools** (4 tools)
  - VTK export for visualization
  - Workspace session management
  - Object export for NGSolve coupling

- **Infrastructure**
  - MCP protocol implementation with stdio transport
  - Asynchronous tool execution
  - State management for Radia objects
  - Comprehensive error handling and logging

### Documentation
- Installation instructions
- Tool reference documentation
- Usage examples for Claude Desktop integration
- Standalone and coupled workflow examples

## Release Notes

### Version 1.1.0 Highlights

This release significantly enhances the Radia MCP server with:

1. **NGSolve Integration**: Direct support for Radia-NGSolve coupling with H-matrix acceleration, enabling efficient field evaluation for large-scale FEM problems.

2. **Diagnostic Tools**: New debugging and monitoring capabilities for production deployments, including server health checks and state inspection.

3. **Performance**: H-matrix batch evaluation provides 10-100x speedup for large point clouds (>1000 points), critical for NGSolve mesh integration.

4. **Documentation**: Comprehensive usage examples and best practices for electromagnetic simulation workflows.

### Migration Notes

No breaking changes. All existing tools maintain backward compatibility.

### Known Issues

None reported.

### Contributors

- Research Lab Team
- Co-Authored-By: Claude Sonnet 4.5 <noreply@anthropic.com>
