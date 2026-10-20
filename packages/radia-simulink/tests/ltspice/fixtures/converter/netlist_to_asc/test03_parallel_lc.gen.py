#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Auto-generated schemdraw script from test03_parallel_lc.asc"""

import schemdraw
import schemdraw.elements as elm

with schemdraw.Drawing(show=False) as d:
    d.config(unit=3, font='Times New Roman')

    I1 = d.add(elm.SourceI().up().label('I1\nAC 1', loc='left'))
    d.add(elm.Dot())
    d.push()
    L1 = d.add(elm.Inductor2().at(I1.end).down().label('L1\n1m', loc='left'))
    d.add(elm.Dot())
    d.push()
    C1 = d.add(elm.Capacitor().at(L1.end).down().label('C1\n100n', loc='left'))
    d.pop()
    C1 = d.add(elm.Capacitor().at(L1.end).down().label('C1\n100n', loc='left'))
    d.pop()
    d.push()
    C1 = d.add(elm.Capacitor().at(I1.end).down().label('C1\n100n', loc='left'))
    d.pop()
    d.push()
    C1 = d.add(elm.Capacitor().at(I1.end).down().label('C1\n100n', loc='left'))
    d.pop()
    L1 = d.add(elm.Inductor2().at(I1.end).down().label('L1\n1m', loc='left'))
    d.add(elm.Line().left().to(I1.start))
    d.add(elm.Ground())
    d.add(elm.Dot())
    d.add(elm.Label().label('top', loc='right'))
    # SPICE: .ac dec 100 1k 100k

    d.save('test03_parallel_lc.pdf')
    print('Saved: test03_parallel_lc.pdf')