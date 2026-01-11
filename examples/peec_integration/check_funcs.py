import sys
sys.path.insert(0, '../../src/radia')
import radia as rad

# List all functions containing "Div"
funcs = [f for f in dir(rad) if "Div" in f or "div" in f]
print("Division functions:", funcs)
