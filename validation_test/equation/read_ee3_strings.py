"""Equation Editor's whole string table.

Its status bar describes whatever the pointer is over, so if the palettes are
described cell by cell, the strings enumerate their contents -- which is the
cheap way to compare them against ours without driving the window.
"""
import struct
import pefile

EXE = (r"C:\Program Files\Microsoft Office\root\vfs"
       r"\ProgramFilesCommonX64\Microsoft Shared\EQUATION\EQNEDT32.EXE")

def leaves(e, path=()):
    if hasattr(e, "directory"):
        for c in e.directory.entries:
            i = c.name.string.decode("utf-8", "replace") if c.name else c.id
            yield from leaves(c, path + (i,))
    else:
        yield path, e.data.struct.OffsetToData, e.data.struct.Size

pe = pefile.PE(EXE, fast_load=False)
out = {}
for top in pe.DIRECTORY_ENTRY_RESOURCE.entries:
    if top.id != 6:
        continue
    for path, rva, size in leaves(top):
        raw = pe.get_data(rva, size)
        off = 0
        for i in range(16):
            if off + 2 > len(raw):
                break
            n = struct.unpack_from("<H", raw, off)[0]
            off += 2
            if n:
                out[(path[0] - 1) * 16 + i] = raw[off:off + n * 2].decode(
                    "utf-16-le", "replace")
                off += n * 2
for k in sorted(out):
    t = out[k].replace("\n", " | ")
    print("%6d  %s" % (k, t))
