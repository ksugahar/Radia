"""Every record of Equation Editor's key table, one line each.

The companion dumper groups by the command pair so that aliases show up
together -- but records with the same (g, c) and a DIFFERENT kind are not
aliases at all, and merging them put Ctrl+Z beside Down and made the table
unreadable.  This prints what is there.
"""
import struct
import pefile

EXE = (r"C:\Program Files\Microsoft Office\root\vfs"
       r"\ProgramFilesCommonX64\Microsoft Shared\EQUATION\EQNEDT32.EXE")

VK = {0x08:"Backspace",0x09:"Tab",0x0D:"Enter",0x1B:"Escape",0x20:"Space",
      0x21:"PageUp",0x22:"PageDown",0x23:"End",0x24:"Home",0x25:"Left",
      0x26:"Up",0x27:"Right",0x28:"Down",0x2D:"Insert",0x2E:"Delete",
      0xBA:";",0xBB:"=",0xBC:",",0xBD:"-",0xBE:".",0xBF:"/",0xC0:"`",
      0xDB:"[",0xDC:"\\",0xDD:"]",0xDE:"'"}
VK.update({0x70+i:"F%d"%(i+1) for i in range(24)})
VK.update({0x30+i:chr(0x30+i) for i in range(10)})
VK.update({0x41+i:chr(0x41+i) for i in range(26)})

def key(v): return VK.get(v, "VK 0x%02X" % v)

def mods(m):
    return "+".join(p for b, p in ((1,"Ctrl"),(2,"Shift"),(4,"mod4"),(8,"Alt")) if m & b)

def chord(vk, m):
    p = [x for x in (mods(m),) if x]; p.append(key(vk)); return "+".join(p)

def leaves(e, path=()):
    if hasattr(e, "directory"):
        for c in e.directory.entries:
            i = c.name.string.decode("utf-8","replace") if c.name else c.id
            yield from leaves(c, path + (i,))
    else:
        yield path, e.data.struct.OffsetToData, e.data.struct.Size

pe = pefile.PE(EXE, fast_load=False)
for top in pe.DIRECTORY_ENTRY_RESOURCE.entries:
    if top.id != 10:
        continue
    for path, rva, size in leaves(top):
        raw = pe.get_data(rva, size)
        print("%-22s %-5s %-6s %s" % ("chord", "kind", "g,c", "note"))
        for i in range(size // 12):
            live, vk, m, kind, g, c = struct.unpack_from("<6h", raw, i * 12)
            if live == 0:
                break
            note = ""
            if kind == 4:
                note = "insert U+%04X" % (c & 0xFFFF)
            print("%-22s %-5d %-6s %s" % (chord(vk, m), kind, "%d,%d" % (g, c), note))
