.include "stage2.inc"
.subckt rc_stage in out R=1k C=1u
R1 in out {R}
C1 out 0 {C}
.ends rc_stage
