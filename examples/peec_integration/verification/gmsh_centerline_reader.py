"""
GMSH v2.2 Centerline Reader for PEEC

Reads 1D edge elements (wire centerlines) from GMSH v2.2 ASCII files.

Author: Radia Development Team
Date: 2026-02-13
"""

import numpy as np


class GMSHCenterlineReader:
    """Read 1D edge elements from GMSH v2.2 ASCII files"""

    def __init__(self, filename):
        """
        Args:
            filename: Path to GMSH v2.2 ASCII file
        """
        self.filename = filename
        self.nodes = {}          # node_id -> [x, y, z]
        self.edge_elements = []  # List of [node1_id, node2_id]
        self.physical_groups = {}  # group_id -> name

    def read(self):
        """Read GMSH file and extract nodes and 1D edge elements"""
        with open(self.filename, 'r') as f:
            lines = f.readlines()

        i = 0
        while i < len(lines):
            line = lines[i].strip()

            # Read physical names
            if line == '$PhysicalNames':
                i = self._read_physical_names(lines, i + 1)

            # Read nodes
            elif line == '$Nodes':
                i = self._read_nodes(lines, i + 1)

            # Read elements
            elif line == '$Elements':
                i = self._read_elements(lines, i + 1)

            i += 1

        return self.nodes, self.edge_elements

    def _read_physical_names(self, lines, start_idx):
        """Read $PhysicalNames section"""
        num_names = int(lines[start_idx].strip())
        idx = start_idx + 1

        for _ in range(num_names):
            parts = lines[idx].strip().split()
            dimension = int(parts[0])
            group_id = int(parts[1])
            name = parts[2].strip('"')
            self.physical_groups[group_id] = name
            idx += 1

        # Skip to $EndPhysicalNames
        while lines[idx].strip() != '$EndPhysicalNames':
            idx += 1

        return idx

    def _read_nodes(self, lines, start_idx):
        """Read $Nodes section"""
        num_nodes = int(lines[start_idx].strip())
        idx = start_idx + 1

        for _ in range(num_nodes):
            parts = lines[idx].strip().split()
            node_id = int(parts[0])
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
            self.nodes[node_id] = [x, y, z]
            idx += 1

        # Skip to $EndNodes
        while lines[idx].strip() != '$EndNodes':
            idx += 1

        return idx

    def _read_elements(self, lines, start_idx):
        """Read $Elements section (only 1D edge elements)"""
        num_elements = int(lines[start_idx].strip())
        idx = start_idx + 1

        for _ in range(num_elements):
            parts = lines[idx].strip().split()
            element_id = int(parts[0])
            element_type = int(parts[1])

            # Element type 1 = 2-node line (1D edge)
            if element_type == 1:
                num_tags = int(parts[2])
                # Tags: physical_group, elementary_entity, ...
                # Node IDs start after tags
                node_start_idx = 3 + num_tags
                node1 = int(parts[node_start_idx])
                node2 = int(parts[node_start_idx + 1])
                self.edge_elements.append([node1, node2])

            idx += 1

        # Skip to $EndElements
        while lines[idx].strip() != '$EndElements':
            idx += 1

        return idx

    def get_segments(self):
        """
        Get wire segments as [start_point, end_point] pairs.

        Returns:
            List of segments: [[start, end], ...]
            where start/end are [x, y, z] coordinates
        """
        segments = []
        for node1_id, node2_id in self.edge_elements:
            start = self.nodes[node1_id]
            end = self.nodes[node2_id]
            segments.append([start, end])
        return segments

    def get_total_length(self):
        """Compute total wire length"""
        total_length = 0.0
        for segment in self.get_segments():
            start, end = segment
            length = np.linalg.norm(np.array(end) - np.array(start))
            total_length += length
        return total_length

    def print_info(self):
        """Print mesh information"""
        print(f"GMSH Centerline Mesh Info:")
        print(f"  Nodes: {len(self.nodes)}")
        print(f"  Edge elements: {len(self.edge_elements)}")
        print(f"  Total length: {self.get_total_length():.6f} m")
        if self.physical_groups:
            print(f"  Physical groups: {self.physical_groups}")


if __name__ == "__main__":
    # Test with example file
    reader = GMSHCenterlineReader("circular_wire_centerline.msh")
    nodes, edges = reader.read()
    reader.print_info()

    # Print first few segments
    segments = reader.get_segments()
    print(f"\nFirst 3 segments:")
    for i, seg in enumerate(segments[:3]):
        start, end = seg
        length = np.linalg.norm(np.array(end) - np.array(start))
        print(f"  Segment {i}: {start} -> {end} (length: {length:.6f} m)")
