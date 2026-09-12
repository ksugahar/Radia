"""Offline outline-conversion checks: never load a native module or register fonts."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from fontTools.pens.areaPen import AreaPen
from fontTools.pens.boundsPen import BoundsPen
from fontTools.ttLib import TTFont

ROOT = Path(__file__).resolve().parents[1]


def bounds(glyph_set, name):
    pen = BoundsPen(glyph_set)
    glyph_set[name].draw(pen)
    return pen.bounds


class ConversionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.scratch = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.scratch.cleanup)
        cls.output = Path(cls.scratch.name) / "diagnostic.ttf"
        subprocess.run([sys.executable, str(ROOT / "build/make_diagnostic_truetype.py"),
                        str(ROOT / "assets/latinmodern-math.otf"), str(cls.output)],
                       check=True, capture_output=True, text=True)
        cls.source = TTFont(ROOT / "assets/latinmodern-math.otf")
        cls.converted = TTFont(cls.output)
        cls.addClassCleanup(cls.source.close)
        cls.addClassCleanup(cls.converted.close)

    def test_complete_geometry_and_winding(self):
        a, b = self.source.getGlyphSet(), self.converted.getGlyphSet()
        self.assertEqual(self.source.getGlyphOrder(), self.converted.getGlyphOrder())
        for name in self.source.getGlyphOrder():
            with self.subTest(glyph=name):
                x, y = bounds(a, name), bounds(b, name)
                self.assertEqual(x is None, y is None)
                if x is None:
                    continue
                self.assertLessEqual(max(abs(u-v) for u, v in zip(x, y)), 1.5)
                pa, pb = AreaPen(a), AreaPen(b)
                a[name].draw(pa)
                b[name].draw(pb)
                # cu2qu reverses CFF winding for the TrueType convention.
                self.assertLessEqual(pa.value * pb.value, 0)
                self.assertEqual(self.source['hmtx'].metrics[name][0],
                                 self.converted['hmtx'].metrics[name][0])

    def test_original_bearing_mutation_is_detected(self):
        name = 'u1D664'
        original = self.converted['hmtx'].metrics[name]
        try:
            self.converted['hmtx'].metrics[name] = self.source['hmtx'].metrics[name]
            x = bounds(self.source.getGlyphSet(), name)
            y = bounds(self.converted.getGlyphSet(), name)
            self.assertGreater(max(abs(u-v) for u, v in zip(x, y)), 16)
        finally:
            self.converted['hmtx'].metrics[name] = original

    def test_format_tables_and_diagnostic_contract(self):
        self.assertIn('glyf', self.converted)
        self.assertNotIn('CFF ', self.converted)
        for tag in ('cmap', 'MATH', 'GSUB', 'GPOS'):
            self.assertEqual(self.source.getTableData(tag), self.converted.getTableData(tag))
        report = json.loads(self.output.with_suffix('.json').read_text(encoding='utf-8'))
        self.assertTrue(report['diagnostic_only'])
        self.assertTrue(report['advance_widths_preserved'])


if __name__ == '__main__':
    unittest.main()
