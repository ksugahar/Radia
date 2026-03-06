#!/usr/bin/env python3
"""
Apply security fixes to Radia source code
"""

import re

def fix_radpy_cpp():
	"""Fix buffer overflow in radpy.cpp CombErStr function"""
	file_path = r'S:\Visual_Studio\02_Visual_Studio_2022_コマンドライン_コンパイル\04_Radia\src\python\radpy.cpp'

	with open(file_path, 'r', encoding='utf-8') as f:
		content = f.read()

	# Fix CombErStr function
	old_code = r'static char\* CombErStr\(const char\* s1, const char\* s2\)\s*\{\s*return strcat\(strcpy\(g_strErTot, s1\), s2\);\s*\}'

	new_code = '''static char* CombErStr(const char* s1, const char* s2)
{
	// Use safe string operations to prevent buffer overflow
	size_t len1 = strlen(s1);
	size_t len2 = strlen(s2);
	size_t total_len = len1 + len2;

	if(total_len >= sizeof(g_strErTot)) {
		// Truncate if necessary to fit in buffer
		strncpy(g_strErTot, s1, sizeof(g_strErTot) - 1);
		size_t remaining = sizeof(g_strErTot) - strlen(g_strErTot) - 1;
		if(remaining > 0) {
			strncat(g_strErTot, s2, remaining);
		}
		g_strErTot[sizeof(g_strErTot) - 1] = '\\0';
	} else {
		strcpy(g_strErTot, s1);
		strcat(g_strErTot, s2);
	}
	return g_strErTot;
}'''

	content = re.sub(old_code, new_code, content, flags=re.MULTILINE)

	with open(file_path, 'w', encoding='utf-8') as f:
		f.write(content)

	print(f"[OK] Fixed CombErStr in {file_path}")

def fix_pyparse_h():
	"""Fix buffer overflow in pyparse.h CopyPyStringToC function"""
	file_path = r'S:\Visual_Studio\02_Visual_Studio_2022_コマンドライン_コンパイル\04_Radia\src\python\pyparse.h'

	with open(file_path, 'r', encoding='utf-8') as f:
		content = f.read()

	# Fix array bounds checking
	old_code = 'if(len > maxLenStr) len = maxLenStr;'
	new_code = 'if(len >= maxLenStr) len = maxLenStr - 1;'

	content = content.replace(old_code, new_code)

	with open(file_path, 'w', encoding='utf-8') as f:
		f.write(content)

	print(f"[OK] Fixed CopyPyStringToC in {file_path}")

def remove_unnecessary_incref():
	"""Remove unnecessary Py_XINCREF calls after Py_BuildValue"""
	file_path = r'S:\Visual_Studio\02_Visual_Studio_2022_コマンドライン_コンパイル\04_Radia\src\python\radpy.cpp'

	with open(file_path, 'r', encoding='utf-8') as f:
		lines = f.readlines()

	new_lines = []
	count = 0
	i = 0
	while i < len(lines):
		line = lines[i]
		# Check if this line is Py_BuildValue and next line is Py_XINCREF
		if 'Py_BuildValue' in line and i + 1 < len(lines):
			next_line = lines[i + 1]
			if 'Py_XINCREF' in next_line and 'oResInd' in next_line:
				# Keep current line, skip next line
				new_lines.append(line)
				new_lines.append('\t\t// Py_XINCREF removed - Py_BuildValue already returns new reference\n')
				i += 2
				count += 1
				continue
		new_lines.append(line)
		i += 1

	with open(file_path, 'w', encoding='utf-8') as f:
		f.writelines(new_lines)

	print(f"[OK] Removed {count} unnecessary Py_XINCREF calls")

if __name__ == '__main__':
	print("Applying security fixes...")
	print()

	fix_radpy_cpp()
	fix_pyparse_h()
	remove_unnecessary_incref()

	print()
	print("[OK] All security fixes applied successfully!")
