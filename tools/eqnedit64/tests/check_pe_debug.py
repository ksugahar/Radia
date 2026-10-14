"""Verify public PE CodeView records contain only the expected PDB basename."""
import struct
import sys
from pathlib import Path


def check(path):
    data = Path(path).read_bytes()
    pe = struct.unpack_from("<I", data, 0x3C)[0]
    assert data[pe:pe + 4] == b"PE\0\0", "Not PE"
    sections, optional_size = struct.unpack_from("<H12xH", data, pe + 6)
    optional = pe + 24
    magic = struct.unpack_from("<H", data, optional)[0]
    assert magic == 0x20B, "Expected PE32+"
    rva, size = struct.unpack_from("<II", data, optional + 112 + 6 * 8)
    table = optional + optional_size
    offset = None
    for index in range(sections):
        virtual_size, virtual_address, raw_size, raw_offset = struct.unpack_from("<IIII", data, table + index * 40 + 8)
        if virtual_address <= rva < virtual_address + max(virtual_size, raw_size):
            offset = raw_offset + rva - virtual_address
    assert offset is not None and size % 28 == 0, "Missing debug directory"
    names = []
    for entry in range(offset, offset + size, 28):
        kind, length, _, pointer = struct.unpack_from("<IIII", data, entry + 12)
        if kind != 2:
            continue
        record = data[pointer:pointer + length]
        assert record[:4] == b"RSDS", "Unsupported CodeView record"
        name = record[24:].split(b"\0", 1)[0].decode("utf-8")
        assert name == "Eqnedit64.pdb", "Public PDB path must be the product basename"
        names.append(name)
    assert names, "No CodeView record"
    print("PASS: PE debug basename")


if __name__ == "__main__":
    check(sys.argv[1])
