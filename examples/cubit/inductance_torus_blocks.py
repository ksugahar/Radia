import cubit
surfaces = cubit.get_relatives("volume", 1, "surface")
for sid in surfaces:
	s = cubit.surface(sid)
	if s.is_planar():
		cy = s.center_point()[1]
		if cy > 0:
			cubit.cmd(f'block 1 add tri in surface {sid}')
			cubit.cmd('block 1 name "source"')
			print(f"Source: surface {sid} (y={cy:.6f})")
		else:
			cubit.cmd(f'block 2 add tri in surface {sid}')
			cubit.cmd('block 2 name "sink"')
			print(f"Sink:   surface {sid} (y={cy:.6f})")
