#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Auto-generated schemdraw script from test05_voltage_rc.asc"""

import schemdraw
import schemdraw.elements as elm

with schemdraw.Drawing(show=False) as d:
    d.config(unit=3, font='Times New Roman')

    V1 = d.add(elm.SourceV().up().label('V1\nAC 1', loc='left'))
    R1 = d.add(elm.Resistor().at(V1.end).right().label('R1\n10k'))
    d.add(elm.Dot())
    d.push()
    C1 = d.add(elm.Capacitor().at(R1.end).down().label('C1\n10n', loc='left'))
    d.pop()
    d.add(elm.Line().left().to(V1.start))
    d.add(elm.Ground())
    d.add(elm.Dot())
    d.add(elm.Label().label('in', loc='right'))
    d.add(elm.Dot())
    d.add(elm.Label().label('out', loc='right'))
    # SPICE: .ac dec 100 100 1Meg

    d.save('test05_voltage_rc.pdf')
    print('Saved: test05_voltage_rc.pdf')