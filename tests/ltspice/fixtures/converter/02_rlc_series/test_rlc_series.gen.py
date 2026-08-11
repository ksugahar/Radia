#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Auto-generated schemdraw script from test_rlc_series.asc"""

import schemdraw
import schemdraw.elements as elm

with schemdraw.Drawing(show=False) as d:
    d.config(unit=3, font='Times New Roman')

    I1 = d.add(elm.SourceI().up().label('I1\nAC 1', loc='left'))
    d.add(elm.Dot())
    d.push()
    R1 = d.add(elm.Resistor().at(I1.end).right().label('R1\n100'))
    L1 = d.add(elm.Inductor2().at(R1.end).right().label('L1\n10m'))
    d.add(elm.Dot())
    d.push()
    C1 = d.add(elm.Capacitor().at(L1.end).down().label('C1\n1u', loc='left'))
    d.pop()
    d.pop()
    d.push()
    L1 = d.add(elm.Inductor2().at(I1.end).right().label('L1\n10m'))
    d.pop()
    R1 = d.add(elm.Resistor().at(I1.end).right().label('R1\n100'))
    d.add(elm.Line().left().to(I1.start))
    d.add(elm.Ground())
    d.add(elm.Dot())
    d.add(elm.Label().label('out', loc='right'))
    # SPICE: .ac dec 100 10 100k

    d.save('test_rlc_series.pdf')
    print('Saved: test_rlc_series.pdf')