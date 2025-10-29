"""
Simple Radia Test
Quick test to verify basic functionality
"""

import sys
import os

# Set UTF-8 encoding for output
if sys.platform == 'win32':
	import codecs
	sys.stdout = codecs.getwriter('utf-8')(sys.stdout.buffer, 'strict')
	sys.stderr = codecs.getwriter('utf-8')(sys.stderr.buffer, 'strict')

# Add build output directory to path
build_dir = os.path.join(os.path.dirname(__file__), 'build', 'lib', 'Release')
sys.path.insert(0, build_dir)

print("=" * 60)
print("Simple Radia Module Test")
print("=" * 60)

try:
	print("\n1. Importing radia module...")
	import radia as rad
	print(f"   [OK] SUCCESS")
	print(f"   Module: {rad.__file__}")

	print("\n2. Getting version...")
	version = rad.UtiVer()
	print(f"   [OK] Version: {version}")

	print("\n3. Creating a rectangular magnet (10x10x10 mm)...")
	magnet = rad.ObjRecMag([0, 0, 0], [10, 10, 10])
	print(f"   [OK] Magnet ID: {magnet}")

	print("\n4. Setting magnetization (1000 A/m in Z direction)...")
	rad.ObjSetM(magnet, [0, 0, 1000])
	print(f"   [OK] Magnetization set")

	print("\n5. Calculating field at point (0, 0, 20) mm...")
	field = rad.Fld(magnet, 'b', [0, 0, 20])
	print(f"   [OK] Field (Tesla):")
	print(f"      Bx = {field[0]:.6f}")
	print(f"      By = {field[1]:.6f}")
	print(f"      Bz = {field[2]:.6f}")

	print("\n6. Cleaning up...")
	rad.UtiDelAll()
	print(f"   [OK] All objects deleted")

	print("\n" + "=" * 60)
	print("[OK] ALL TESTS PASSED!")
	print("=" * 60)

except Exception as e:
	print(f"\n[FAIL] ERROR: {e}")
	import traceback
	traceback.print_exc()
	sys.exit(1)
