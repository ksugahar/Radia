* PRIMA + Dowell Skin Effect SPICE netlist
* Wire: L=100mm, w=1.0mm, h=1.0mm
* PRIMA stages: 5, Dowell stages: 5
* Skin effect: Dowell continued fraction ladder
* Skin effect significant above: ~8.7 kHz

.SUBCKT WIRE_PRIMA_SKIN port_in port_out

* === Dowell Skin Effect Ladder ===
Rskin0 port_in skin_n0 1.724138e-06
Rskin1 skin_n0 skin_n1 5.747126e-07
Lskin1 skin_n1 0 2.094395e-11
Rskin2 skin_n1 skin_n2 3.448276e-07
Lskin2 skin_n2 0 1.256637e-11
Rskin3 skin_n2 skin_n3 2.463054e-07
Lskin3 skin_n3 0 8.975979e-12
Rskin4 skin_n3 skin_n4 1.915709e-07
Lskin4 skin_n4 0 6.981317e-12
Rskin5 skin_n4 skin_out 1.567398e-07

* === PRIMA External Inductance Ladder ===
Lext1 skin_out ext_n1 5.562684e-09
Lext2 ext_n1 ext_n2 8.019373e-09
Lext3 ext_n2 ext_n3 6.607106e-09
Lext4 ext_n3 ext_n4 5.885469e-09
Lext5 ext_n4 port_out 5.747754e-09

* === PRIMA Tridiagonal Coupling ===
Kext1_2 Lext1 Lext2 2.725819e-01
Kext2_3 Lext2 Lext3 1.872218e-01
Kext3_4 Lext3 Lext4 1.566764e-01
Kext4_5 Lext4 Lext5 1.311712e-01

.ENDS WIRE_PRIMA_SKIN

* === Test Circuit ===
Xwire port_in port_out WIRE_PRIMA_SKIN
Vin port_in 0 AC 1

* AC analysis: 1kHz to 100MHz
.AC DEC 20 1k 100MEG
.PRINT AC V(port_out) I(Vin)
.END