#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Auto-generated schemdraw script from test_pi_filter.asc"""

import schemdraw
import schemdraw.elements as elm

with schemdraw.Drawing(show=False) as d:
    d.config(unit=3, font='Times New Roman')

    I1 = d.add(elm.SourceI().up().label('I1\nAC 1', loc='left'))
    d.add(elm.Dot())
    d.push()
    L1 = d.add(elm.Inductor2().at(I1.end).right().label('L1\n1m'))
    C1 = d.add(elm.Capacitor().at(L1.end).down().label('C1\n100n', loc='left'))
    d.add(elm.Dot())
    d.push()
    C2 = d.add(elm.Capacitor().at(C1.end).down().label('C2\n100n', loc='left'))
    d.add(elm.Dot())
    d.push()
    RL = d.add(elm.Resistor().at(C2.end).down().label('RL\n50', loc='left'))
    d.pop()
    RL = d.add(elm.Resistor().at(C2.end).down().label('RL\n50', loc='left'))
    d.pop()
    d.push()
    RL = d.add(elm.Resistor().at(C1.end).down().label('RL\n50', loc='left'))
    d.pop()
    d.push()
    RL = d.add(elm.Resistor().at(C1.end).down().label('RL\n50', loc='left'))
    d.pop()
    C2 = d.add(elm.Capacitor().at(C1.end).down().label('C2\n100n', loc='left'))
    d.pop()
    d.push()
    C2 = d.add(elm.Capacitor().at(I1.end).down().label('C2\n100n', loc='left'))
    d.pop()
    d.push()
    C2 = d.add(elm.Capacitor().at(I1.end).down().label('C2\n100n', loc='left'))
    d.pop()
    d.push()
    RL = d.add(elm.Resistor().at(I1.end).down().label('RL\n50', loc='left'))
    d.pop()
    d.push()
    RL = d.add(elm.Resistor().at(I1.end).down().label('RL\n50', loc='left'))
    d.pop()
    C1 = d.add(elm.Capacitor().at(I1.end).down().label('C1\n100n', loc='left'))
    d.add(elm.Line().left().to(I1.start))
    d.add(elm.Ground())
    d.add(elm.Dot())
    d.add(elm.Label().label('out', loc='right'))
    # SPICE: .ac dec 100 1k 1Meg

    d.save('test_pi_filter.pdf')
    print('Saved: test_pi_filter.pdf')