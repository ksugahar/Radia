"""Equation Editor's own menus, with the accelerator text it shows the user.

Written to settle whether PageUp, PageDown and Ctrl+Tab are commands a user
could find and name -- they are not: no menu mentions them (see the handover
section on the keys deliberately left unbound).  Kept because "is this in a
menu, and what does the program call it?" is the first question to ask of any
chord found in the key table.

A menu item names what a command DOES, in the words the program itself uses --
which is what decides whether a key may be bound here at all.  Reading a
resource to learn behaviour is fine; it is the artwork that must never be
copied.
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

def wstr(buf, off):
    out = []
    while off + 1 < len(buf):
        (ch,) = struct.unpack_from("<H", buf, off)
        off += 2
        if ch == 0:
            break
        out.append(chr(ch))
    return "".join(out), off

def walk(buf, off, depth, seen):
    """MENU template: each item is flags, [id], text."""
    while off + 2 <= len(buf):
        (flags,) = struct.unpack_from("<H", buf, off)
        off += 2
        ident = None
        if not (flags & 0x10):                      # not a popup
            if off + 2 > len(buf):
                return off
            (ident,) = struct.unpack_from("<H", buf, off)
            off += 2
        text, off = wstr(buf, off)
        label = "  " * depth + (text if text else "-" * 8)
        print(f"{label:<52} {'' if ident is None else ident}")
        if flags & 0x10:
            off = walk(buf, off, depth + 1, seen)
        if flags & 0x80:                            # last item of this popup
            return off
    return off

pe = pefile.PE(EXE, fast_load=False)
for top in pe.DIRECTORY_ENTRY_RESOURCE.entries:
    if top.id != 4:                                 # RT_MENU
        continue
    for path, rva, size in leaves(top):
        buf = pe.get_data(rva, size)
        print(f"==== MENU {path} ({size} bytes) ====")
        ver, hdr = struct.unpack_from("<HH", buf, 0)
        start = 4 if ver == 0 else 4
        try:
            walk(buf, start, 0, set())
        except Exception as ex:
            print("  (parse stopped:", ex, ")")
