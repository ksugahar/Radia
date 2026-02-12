"""
Shared Workspace Manager for Radia-NGSolve MCP Servers

Manages shared workspace for inter-server communication:
- Session creation and management
- Object serialization and export
- File-based data exchange
"""

import json
import os
import shutil
from pathlib import Path
from typing import Any, Dict, List, Optional
import uuid
from datetime import datetime


class SharedWorkspace:
    """Manage shared workspace for Radia-NGSolve communication."""

    def __init__(self, base_path: Optional[str] = None):
        """
        Initialize shared workspace.

        Args:
            base_path: Base directory for workspace (default: ~/.radia-ngsolve-workspace)
        """
        if base_path is None:
            base_path = os.path.expanduser("~/.radia-ngsolve-workspace")
        self.base_path = Path(base_path)
        self.base_path.mkdir(parents=True, exist_ok=True)

        # Create config directory
        config_dir = self.base_path / "config"
        config_dir.mkdir(exist_ok=True)

    def create_session(self, session_name: Optional[str] = None) -> str:
        """
        Create a new session and return session ID.

        Args:
            session_name: Optional human-readable session name

        Returns:
            Session ID (UUID)
        """
        session_id = str(uuid.uuid4())
        session_dir = self.base_path / "sessions" / session_id

        # Create directory structure
        (session_dir / "radia" / "objects").mkdir(parents=True, exist_ok=True)
        (session_dir / "radia" / "geometry").mkdir(parents=True, exist_ok=True)
        (session_dir / "radia" / "fields").mkdir(parents=True, exist_ok=True)
        (session_dir / "ngsolve" / "meshes").mkdir(parents=True, exist_ok=True)
        (session_dir / "ngsolve" / "solutions").mkdir(parents=True, exist_ok=True)

        # Create session metadata
        metadata = {
            "session_id": session_id,
            "session_name": session_name or f"Session {session_id[:8]}",
            "created_at": datetime.now().isoformat(),
            "radia_objects": [],
            "ngsolve_objects": []
        }
        self._write_session_metadata(session_id, metadata)

        return session_id

    def export_radia_object(
        self,
        session_id: str,
        obj_name: str,
        radia_obj_id: int,
        geometry_vtk_path: Optional[str] = None,
        additional_metadata: Optional[Dict[str, Any]] = None
    ) -> str:
        """
        Export Radia object to shared workspace.

        Args:
            session_id: Session ID
            obj_name: Object name
            radia_obj_id: Radia object ID
            geometry_vtk_path: Path to VTK geometry file (optional)
            additional_metadata: Additional metadata to store

        Returns:
            Path to exported object metadata file
        """
        session_dir = self.base_path / "sessions" / session_id / "radia"

        # Export object metadata
        obj_metadata = {
            "name": obj_name,
            "object_id": radia_obj_id,
            "type": "radia_object",
            "exported_at": datetime.now().isoformat()
        }

        if additional_metadata:
            obj_metadata.update(additional_metadata)

        obj_file = session_dir / "objects" / f"{obj_name}.json"
        with open(obj_file, 'w') as f:
            json.dump(obj_metadata, f, indent=2)

        # Export geometry if provided
        vtk_file = None
        if geometry_vtk_path and os.path.exists(geometry_vtk_path):
            vtk_file = session_dir / "geometry" / f"{obj_name}.vtk"
            shutil.copy(geometry_vtk_path, vtk_file)

        # Update session metadata
        metadata = self._read_session_metadata(session_id)
        metadata["radia_objects"].append({
            "name": obj_name,
            "object_file": str(obj_file),
            "geometry_file": str(vtk_file) if vtk_file else None,
            "exported_at": datetime.now().isoformat()
        })
        self._write_session_metadata(session_id, metadata)

        return str(obj_file)

    def import_radia_object(self, session_id: str, obj_name: str) -> Dict[str, Any]:
        """
        Import Radia object from shared workspace.

        Args:
            session_id: Session ID
            obj_name: Object name

        Returns:
            Dictionary with metadata and file paths
        """
        session_dir = self.base_path / "sessions" / session_id / "radia"
        obj_file = session_dir / "objects" / f"{obj_name}.json"

        if not obj_file.exists():
            raise FileNotFoundError(
                f"Radia object '{obj_name}' not found in session {session_id}"
            )

        with open(obj_file, 'r') as f:
            metadata = json.load(f)

        # Check for associated files
        vtk_file = session_dir / "geometry" / f"{obj_name}.vtk"
        field_file = session_dir / "fields" / f"{obj_name}_field.npz"

        return {
            "metadata": metadata,
            "geometry_file": str(vtk_file) if vtk_file.exists() else None,
            "field_file": str(field_file) if field_file.exists() else None
        }

    def export_field_data(
        self,
        session_id: str,
        obj_name: str,
        points: List[List[float]],
        field_values: List[List[float]],
        field_type: str = "b"
    ) -> str:
        """
        Export pre-computed field data.

        Args:
            session_id: Session ID
            obj_name: Object name
            points: List of evaluation points [[x,y,z], ...]
            field_values: List of field values corresponding to points
            field_type: Type of field ('b', 'h', 'a', 'phi')

        Returns:
            Path to exported field data file
        """
        import numpy as np
        session_dir = self.base_path / "sessions" / session_id / "radia" / "fields"

        field_file = session_dir / f"{obj_name}_field.npz"
        np.savez_compressed(
            field_file,
            points=np.array(points),
            field_values=np.array(field_values),
            field_type=field_type,
            exported_at=datetime.now().isoformat()
        )

        # Update session metadata
        metadata = self._read_session_metadata(session_id)
        for obj in metadata["radia_objects"]:
            if obj["name"] == obj_name:
                obj["field_file"] = str(field_file)
                obj["field_type"] = field_type
                break
        self._write_session_metadata(session_id, metadata)

        return str(field_file)

    def export_ngsolve_mesh(
        self,
        session_id: str,
        mesh_name: str,
        mesh_file_path: str
    ) -> str:
        """
        Export NGSolve mesh to shared workspace.

        Args:
            session_id: Session ID
            mesh_name: Mesh name
            mesh_file_path: Path to mesh file

        Returns:
            Path to exported mesh file
        """
        session_dir = self.base_path / "sessions" / session_id / "ngsolve" / "meshes"

        # Copy mesh file
        mesh_file = session_dir / f"{mesh_name}.vol.gz"
        shutil.copy(mesh_file_path, mesh_file)

        # Update session metadata
        metadata = self._read_session_metadata(session_id)
        metadata["ngsolve_objects"].append({
            "name": mesh_name,
            "type": "mesh",
            "mesh_file": str(mesh_file),
            "exported_at": datetime.now().isoformat()
        })
        self._write_session_metadata(session_id, metadata)

        return str(mesh_file)

    def list_sessions(self) -> List[Dict[str, Any]]:
        """
        List all available sessions.

        Returns:
            List of session metadata dictionaries
        """
        sessions_dir = self.base_path / "sessions"
        if not sessions_dir.exists():
            return []

        sessions = []
        for session_dir in sessions_dir.iterdir():
            if session_dir.is_dir():
                try:
                    metadata = self._read_session_metadata(session_dir.name)
                    sessions.append(metadata)
                except Exception:
                    # Skip corrupted sessions
                    continue

        return sessions

    def get_session_info(self, session_id: str) -> Dict[str, Any]:
        """
        Get detailed session information.

        Args:
            session_id: Session ID

        Returns:
            Session metadata dictionary
        """
        metadata = self._read_session_metadata(session_id)

        # Add statistics
        metadata["statistics"] = {
            "num_radia_objects": len(metadata.get("radia_objects", [])),
            "num_ngsolve_objects": len(metadata.get("ngsolve_objects", []))
        }

        return metadata

    def delete_session(self, session_id: str):
        """
        Delete a session and all its data.

        Args:
            session_id: Session ID
        """
        session_dir = self.base_path / "sessions" / session_id
        if session_dir.exists():
            shutil.rmtree(session_dir)

    def cleanup_old_sessions(self, days: int = 30):
        """
        Delete sessions older than specified days.

        Args:
            days: Number of days to keep
        """
        from datetime import timedelta

        cutoff_date = datetime.now() - timedelta(days=days)
        sessions = self.list_sessions()

        for session in sessions:
            created_at = datetime.fromisoformat(session["created_at"])
            if created_at < cutoff_date:
                self.delete_session(session["session_id"])

    def _read_session_metadata(self, session_id: str) -> Dict[str, Any]:
        """Read session metadata."""
        metadata_file = self.base_path / "sessions" / session_id / "metadata.json"
        if metadata_file.exists():
            with open(metadata_file, 'r') as f:
                return json.load(f)
        return {}

    def _write_session_metadata(self, session_id: str, metadata: Dict[str, Any]):
        """Write session metadata."""
        metadata_file = self.base_path / "sessions" / session_id / "metadata.json"
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)


# Convenience function for testing
def create_test_session() -> str:
    """Create a test session for development."""
    workspace = SharedWorkspace()
    session_id = workspace.create_session("Test Session")
    print(f"Created test session: {session_id}")
    return session_id


if __name__ == "__main__":
    # Test workspace functionality
    workspace = SharedWorkspace()
    print(f"Workspace location: {workspace.base_path}")

    # Create test session
    session_id = workspace.create_session("Test Session")
    print(f"Created session: {session_id}")

    # List sessions
    sessions = workspace.list_sessions()
    print(f"Total sessions: {len(sessions)}")

    # Get session info
    info = workspace.get_session_info(session_id)
    print(f"Session info: {json.dumps(info, indent=2)}")
