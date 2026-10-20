"""Standard #line mappings keep native diagnostic filenames independent of scratch."""
from pathlib import Path
import hashlib

C_FAMILY={'.h','.hpp','.hxx','.cuh','.inc','.c','.cc','.cpp','.cxx'}


def stage_logical_sources(paths, root, destination, label):
    root=Path(root).resolve();destination=Path(destination);manifest={}
    for source in dict.fromkeys(Path(p) for p in paths):
        source=source.resolve();relative=source.relative_to(root)
        logical=(Path(label)/relative).as_posix()
        if any(c in logical for c in ('"','\n','\r')):
            raise ValueError('Unsupported logical source filename')
        original=source.read_bytes();compiled=original
        if source.suffix.lower() in C_FAMILY:
            bom=b'\xef\xbb\xbf' if original.startswith(b'\xef\xbb\xbf') else b''
            body=original[len(bom):]
            compiled=bom+('#line 1 "'+logical+'"\n').encode('utf-8')+body
        target=destination/relative;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(compiled)
        manifest[relative.as_posix()]=dict(original_sha256=hashlib.sha256(original).hexdigest(),
            compiled_sha256=hashlib.sha256(compiled).hexdigest(),logical_filename=logical)
    return manifest
