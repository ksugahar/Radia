* PRIMA-reduced SPICE netlist
* Wire: L=100mm, w=1.0mm, h=1.0mm
* Original segments: 10, PRIMA stages: 5

.SUBCKT WIRE_PRIMA port_in port_out

* === PRIMA Ladder (Tridiagonal) ===
R1 port_in port_ina 1.724138e-04
L1 port_ina n1 5.791936e-09
R2 n1 n1a 1.724138e-04
L2 n1a n2 7.922445e-09
R3 n2 n2a 1.724138e-04
L3 n2a n3 6.711010e-09
R4 n3 n3a 1.724138e-04
L4 n3a n4 6.027782e-09
R5 n4 n4a 1.724138e-04
L5 n4a port_out 5.400067e-09

* === Tridiagonal Coupling ===
K1_2 L1 L2 2.637116e-01
K2_3 L2 L3 1.974179e-01
K3_4 L3 L4 1.145372e-01
K4_5 L4 L5 1.828444e-01

.ENDS WIRE_PRIMA

* === Test Circuit ===
Xwire port_in port_out WIRE_PRIMA
Vin port_in 0 AC 1

* AC analysis: 1kHz to 100MHz
.AC DEC 10 1k 100MEG
.PRINT AC V(port_out) I(Vin)
.END