# PyLTSpice / LTspice / SPICE API and Format Reference

_This document is a structured reformat of the public PyLTSpice / spicelib API and the LTspice / SPICE file formats. Provided for orientation; consult upstream projects ([PyLTSpice](https://github.com/nunobrum/PyLTSpice), [spicelib](https://github.com/nunobrum/spicelib), [LTspice documentation](https://www.analog.com/en/resources/technical-articles/ltspice-getting-started-guide.html)) for the canonical reference._

## Contents
- [PyLTSpice / spicelib Python API](#circuit-pyltspice-api)
- [LTspice symbol anchor points](#circuit-ltspice-anchor-points)
- [LTspice `.asc` file format](#circuit-ltspice-asc-format)
- [LTspice advanced component reference](#circuit-ltspice-advanced-components)
- [SPICE netlist (`.cir`) format](#circuit-spice-netlist-format)

## PyLTSpice / spicelib Python API

<a id="circuit-pyltspice-api"></a>

```
PyLTSpice / spicelib API Reference
====================================

VERSION: PyLTSpice v5.4+ wraps spicelib v1.4+.
The actual classes live in spicelib; PyLTSpice re-exports them.

IMPORT
------
from PyLTSpice import AscEditor, SpiceEditor, SimRunner, RawRead
# or equivalently:
from spicelib import AscEditor, SpiceEditor, SimRunner, RawRead

KEY DATA CLASSES (from spicelib.editor.base_schematic)
-----------------------------------------------------
Point(X, Y)              -- coordinate on schematic canvas
Line(V1: Point, V2: Point)  -- wire or graphical line
Text(coord, text, size, type)  -- label, directive, or comment
TextTypeEnum             -- .DIRECTIVE, .COMMENT, .LABEL, .ATTRIBUTE
ERotation                -- R0, R90, R180, R270, M0, M90, M180, M270
SchematicComponent       -- .reference, .symbol, .position, .rotation, .attributes

1. EDITING .ASC FILES WITH AscEditor
=====================================

Opening an existing file:
    editor = AscEditor("./circuit.asc")

IMPORTANT: AscEditor requires an existing file. To create from scratch,
write the minimal ASC text first:
    with open('new.asc', 'w') as f:
        f.write('Version 4\\nSHEET 1 880 680\\n')
    editor = AscEditor('new.asc')

Modifying component values:
    editor.set_component_value('R1', '4.7k')
    editor.set_component_value('C1', 100e-9)     # auto-converts to eng notation
    editor.set_component_value('R2', 2000)

Getting component info:
    val = editor.get_component_value('R1')         # string '4.7k'
    fval = editor.get_component_floatvalue('R1')   # float
    info = editor.get_component_info('R1')         # dict
    params = editor.get_component_parameters('R1')

Listing components:
    all_refs = editor.get_components()             # ['R1', 'C1', 'V1', ...]
    resistors = editor.get_components('R')         # only R-prefix
    caps_res = editor.get_components('RC')         # R and C prefix

Changing component models:
    editor.set_element_model('D1', '1N4148')
    editor.set_element_model('V3', "SINE(0 1 3k 0 0 0)")

Component position/rotation:
    from spicelib.editor.base_schematic import Point, ERotation
    pos, rot = editor.get_component_position('R1')
    editor.set_component_position('R1', Point(200, 300), ERotation.R90)

Adding a component programmatically:
    from spicelib.editor.base_schematic import SchematicComponent, Point, ERotation
    comp = SchematicComponent(editor, "SYMBOL res 200 300 R0")
    comp.symbol = 'res'
    comp.position = Point(200, 300)
    comp.rotation = ERotation.R0
    comp.reference = 'R2'
    comp.attributes['Value'] = '10k'
    editor.add_component(comp)

Removing a component:
    editor.remove_component('R2')

Adding wires (no add_wire() method -- use wires list directly):
    from spicelib.editor.base_schematic import Line, Point
    wire = Line(Point(0, 192), Point(128, 192))
    editor.wires.append(wire)
    editor.updated = True

Adding net labels / flags:
    from spicelib.editor.base_schematic import Text, TextTypeEnum
    flag = Text(coord=Point(192, 192), text='out', type=TextTypeEnum.LABEL)
    editor.labels.append(flag)
    editor.updated = True

Simulation directives:
    editor.add_instruction('.tran 10m')
    editor.add_instruction('.ac dec 100 1 100k')
    editor.add_instructions('.meas TRAN vout_max MAX V(out)',
                            '.step param R1 1k 10k 1k')

Parameters:
    editor.set_parameter('freq', '1k')
    editor.set_parameters(R_val=1000, C_val=1e-6)
    # NOTE: use set_parameter(), NOT add_instruction('.param ...')

Removing instructions:
    editor.remove_instruction('.tran 10m')
    editor.remove_Xinstruction(r'\\.step.*')   # regex-based

Saving:
    editor.save_netlist('./output.asc')    # must be .asc extension

Scaling/transforming:
    editor.scale(offset_x=100, offset_y=50, scale_x=2.0, scale_y=2.0)

Subcircuit access:
    sub = editor.get_subcircuit('XU1')
    editor.set_component_value('XU1:C2', 20e-12)  # ':' separator

2. RUNNING SIMULATIONS WITH SimRunner
=======================================

from PyLTSpice import SimRunner
from spicelib.simulators.ltspice_simulator import LTspice

runner = SimRunner(output_folder='./temp')

Run from .asc file:
    runner.run('./circuit.asc')

Create netlist, modify, run:
    LTspice.create_netlist('./circuit.asc')    # produces .net
    netlist = SpiceEditor('./circuit.net')
    netlist.set_component_value('R1', '10k')
    runner.run(netlist)

Run with AscEditor:
    editor = AscEditor('./circuit.asc')
    editor.set_component_value('R1', '10k')
    runner.run(editor)

Wait and iterate:
    runner.wait_completion(timeout=120)
    for raw_file, log_file in runner:
        print(f"Raw: {raw_file}, Log: {log_file}")

Synchronous run:
    raw_path, log_path = runner.run_now('./circuit.asc', timeout=60)

Batch sweep:
    for r_val in ['1k', '4.7k', '10k']:
        netlist.set_component_value('R1', r_val)
        runner.run(netlist)
    runner.wait_completion()
    runner.file_cleanup()

LTspice executable: auto-detected at C:/Program Files/ADI/LTspice/LTspice.exe.
Custom: LTspice.create_from(path)

3. READING RESULTS WITH RawRead
=================================

from PyLTSpice import RawRead

raw = RawRead('./simulation.raw')

List traces:
    raw.get_trace_names()   # ['time', 'V(out)', 'I(R1)', ...]

Get waveform data:
    wave = raw.get_wave('V(out)', step=0)   # numpy array
    time = raw.get_axis(step=0)
    time = raw.get_time_axis(step=0)

Stepped simulations:
    steps = raw.get_steps()
    for step_idx in range(len(steps)):
        t = raw.get_axis(step=step_idx)
        v = raw.get_wave('V(out)', step=step_idx)

Export:
    raw.to_csv('./results.csv')
    raw.to_excel('./results.xlsx')
    df = raw.to_dataframe()

KEY NOTES
---------
- AscEditor requires an existing file -- cannot create blank .asc from scratch
- No add_wire() method -- append Line objects to editor.wires directly
- Labels append to editor.labels directly
- set_parameter() for .param, NOT add_instruction('.param ...')
- Unique sim instructions (.tran, .ac, .dc) auto-replaced on add_instruction()
- LTspice 17+ uses UTF-16 LE encoding; older uses ASCII/UTF-8
- File encoding is auto-detected
```

## LTspice symbol anchor points

<a id="circuit-ltspice-anchor-points"></a>

```
LTSpice Symbol Anchor Points (Measured)
=========================================

All coordinates in LTSpice pixels. Grid size = 16px.
sym(x, y) = SYMBOL placement coordinate in ASC file.

RESISTOR (res)
--------------
R0  (vertical):   sym(x,y) -> upper terminal (x+16, y),    lower terminal (x+16, y+64)
R90 (horizontal): sym(x,y) -> left terminal  (x-64, y+16), right terminal (x, y+16)
                  NOTE: R90 anchor is at the RIGHT terminal side
R180:             sym(x,y) -> lower terminal (x+16, y+64), upper terminal (x+16, y)
R270:             sym(x,y) -> right terminal (x+64, y+16), left terminal  (x, y+16)

INDUCTOR (ind)
--------------
R0  (vertical):   sym(x,y) -> upper terminal (x+16, y-16), lower terminal (x+16, y+96)
R90 (horizontal): sym(x,y) -> left terminal  (x-96, y+16), right terminal (x, y+16)
                  NOTE: R90 anchor is at the RIGHT terminal side
R180:             sym(x,y) -> upper terminal (x+16, y+112), lower terminal (x+16, y)
R270:             sym(x,y) -> right terminal (x+96, y+16), left terminal  (x, y+16)

CAPACITOR (cap)
---------------
R0  (vertical):   sym(x,y) -> upper terminal (x+16, y),    lower terminal (x+16, y+64)
R90 (horizontal): sym(x,y) -> left terminal  (x-64, y+16), right terminal (x, y+16)
                  NOTE: R90 anchor is at the RIGHT terminal side
R180:             sym(x,y) -> lower terminal (x+16, y+64), upper terminal (x+16, y)
R270:             sym(x,y) -> right terminal (x+64, y+16), left terminal  (x, y+16)

CURRENT SOURCE (current)
-------------------------
R0:   sym(x,y) -> upper terminal (x+16, y),     lower terminal (x+16, y+96)
R180: sym(x,y) -> upper terminal (x, y-64),     lower terminal (x, y+64)

VOLTAGE SOURCE (voltage)
-------------------------
R0:   sym(x,y) -> upper terminal(+) (x, y),     lower terminal(-) (x, y+96)
R180: sym(x,y) -> upper terminal(+) (x, y+96),  lower terminal(-) (x, y)

SYMBOL SIZES (width x height in LTSpice px)
--------------------------------------------
res:     R0 (32x64),  R90 (64x32)
ind:     R0 (32x112), R90 (112x32)
cap:     R0 (32x64),  R90 (64x32)
current: R0 (32x96),  R180 (32x96)
voltage: R0 (32x96),  R180 (32x96)

COORDINATE SYSTEM
-----------------
- LTSpice Y-axis points DOWNWARD (opposite of schemdraw)
- Grid snapping: round to nearest multiple of 16
- Scale factor: schemdraw 1 unit = 64 LTSpice pixels
```

## LTspice `.asc` file format

<a id="circuit-ltspice-asc-format"></a>

```
LTSpice .asc File Format Specification
========================================

An ASC file is plain-text (UTF-8 or UTF-16 LE for LTspice 17+).
Version can be "Version 4" or "Version 4.1" (newer LTspice).

HEADER
------
Version 4
SHEET 1 <width> <height>

Typical sheet sizes: 880x680 (small), 1600x680 (medium), 2844x1336 (large)

ELEMENTS
--------

WIRE x1 y1 x2 y2
  Connects two points. Coordinates in LTSpice pixels.
  Grid: 16px. All coordinates should be multiples of 16.
  Wires are always axis-aligned (horizontal or vertical).
  Wires connecting at the same coordinate are electrically joined.
  Negative coordinates are valid and commonly used.

FLAG x y <netname>
  Net label at position (x,y). netname "0" = ground.
  Multiple FLAG 0 entries allowed for multiple ground connections.
  Any other name creates a named net (used for power: +V, -V, Vdd;
  signals: OUT, IN, etc.)
  Example: FLAG 192 192 0       (ground)
           FLAG 192 64 out      (net label "out")
           FLAG 144 -1168 +V    (positive supply rail)

IOPIN x y <direction>
  I/O port marker. Must follow a FLAG line.
  direction: In, Out, BiDir
  Used mainly in hierarchical designs.
  Example: IOPIN 192 64 Out

SYMBOL <symbolname> x y <rotation>
  Component placement. symbolname is relative to sym/ directory.
  Built-in symbols (no path prefix): res, cap, ind, current, voltage
  Subdirectory symbols use backslash:
    opamps\\AD711, References\\AD590, SpecialFunctions\\LTC6905
  Top-level custom symbols: LT1021-7, AD820

  Rotation: R0, R90, R180, R270 (normal)
            M0, M90, M180, M270 (mirrored)
  R0 = vertical (default), R90 = horizontal (rotated 90 deg clockwise)
  M = mirrored version of corresponding R rotation

WINDOW <num> x y <alignment> <size>
  Display positioning for component attributes.
  num: 0=InstName, 3=Value, 123=Value2/SpiceLine, 39=SpiceLine2
  alignment: Left, Right, Top, VBottom, VTop
  size: font size. 0 = hidden, 2 = normal
  Setting position to (0,0) with size 0 effectively hides the window.
  Negative offsets are valid for repositioning.
  Each WINDOW belongs to the preceding SYMBOL.

  Common patterns from real LTSpice files:
    WINDOW 0 0 56 VBottom 2     (for R90 resistors - InstName)
    WINDOW 3 32 56 VTop 2       (for R90 resistors - Value)
    WINDOW 123 0 0 Left 0       (hide SpiceLine)
    WINDOW 39 0 0 Left 0        (hide SpiceLine2)
    WINDOW 0 24 80 Left 2       (for R180 current source - InstName)
    WINDOW 3 24 0 Left 2        (for R180 current source - Value)

SYMATTR <attribute> <value>
  Component attribute. Must follow SYMBOL line (after WINDOW lines).
  Attributes: InstName, Value, Value2, SpiceModel, SpiceLine,
              SpiceLine2, Prefix, Def_Sub
  Value2 used for continuation of long parameters:
    SYMATTR Value PULSE(-1 1 5u
    SYMATTR Value2 1n 1n 5u 10u)
  SpiceLine for subcircuit parameters:
    SYMATTR SpiceLine Gain=1
  Examples:
    SYMATTR InstName R1
    SYMATTR Value 10k
    SYMATTR SpiceModel 1N4148

TEXT x y <alignment> <size> !<directive>
TEXT x y <alignment> <size> ;<comment>
  Simulation directive (! prefix) or comment (; prefix).
  alignment: Left, Right, Center, Top, Bottom
  size: font size (typically 2)
  \\n encodes literal newlines for multi-line text.
  Examples:
    TEXT 0 512 Left 2 !.ac dec 100 1 1Meg
    TEXT 0 544 Left 2 !.tran 10m
    TEXT 584 80 Left 2 !.dc V1 0 10 1m\\n.temp -55 25 150
    TEXT 0 576 Left 2 ;This is a comment

  Common SPICE directives:
    .tran <tstop> [startup]
    .ac dec <npoints> <fstart> <fstop>
    .dc <source> <start> <stop> <step>
    .noise V(out) V1 dec 100 1 100k
    .op
    .meas TRAN <name> <function> <expression>
    .step param <name> <start> <stop> <step>
    .param <name>=<value>
    .temp <temp1> [temp2] ...

  Common source value formats:
    AC 1                              (AC analysis source)
    SINE(offset amplitude frequency)
    PULSE(V1 V2 Tdelay Trise Tfall Ton Tperiod)
    PWL(t1 v1 t2 v2 ...)

GRAPHICAL ELEMENTS (optional)
-----------------------------
LINE Normal x1 y1 x2 y2 [style]
RECTANGLE Normal x1 y1 x2 y2 [style]
CIRCLE Normal x1 y1 x2 y2 [style]
ARC Normal x1 y1 x2 y2 x3 y3 x4 y4 [style]

COMPLETE MINIMAL EXAMPLE (RC lowpass from converter)
-----------------------------------------------------
Version 4
SHEET 1 1600 680
WIRE 0 192 0 384
WIRE 0 192 128 192
WIRE 128 192 192 192
WIRE 192 192 192 256
WIRE 192 256 192 384
WIRE 192 384 0 384
FLAG 0 384 0
FLAG 192 192 out
IOPIN 192 192 Out
SYMBOL current 0 320 R180
WINDOW 0 24 80 Left 2
WINDOW 3 24 0 Left 2
SYMATTR InstName I1
SYMATTR Value AC 1
SYMBOL res 192 176 R90
WINDOW 0 0 56 VBottom 2
WINDOW 3 32 56 VTop 2
SYMATTR InstName R1
SYMATTR Value 1k
SYMBOL cap 176 192 R0
SYMATTR InstName C1
SYMATTR Value 1u
TEXT 0 512 Left 2 !.ac dec 100 1 100k

COMPLETE EXAMPLE (inverting op-amp from AD820.asc)
---------------------------------------------------
Version 4
SHEET 1 1240 700
WIRE 432 -1216 416 -1216
WIRE 528 -1216 512 -1216
WIRE 544 -1216 528 -1216
WIRE 640 -1216 624 -1216
WIRE 528 -1120 528 -1216
WIRE 544 -1120 528 -1120
WIRE 640 -1104 640 -1216
WIRE 640 -1104 608 -1104
WIRE 704 -1104 640 -1104
WIRE 544 -1088 432 -1088
WIRE 432 -1072 432 -1088
WIRE 432 -976 432 -992
SYMBOL voltage 144 -1168 R0
WINDOW 123 0 0 Left 2
WINDOW 39 0 0 Left 2
SYMATTR InstName V1
SYMATTR Value 15
SYMBOL voltage 432 -1088 R0
WINDOW 123 24 146 Left 2
WINDOW 39 24 125 Left 2
SYMATTR InstName Vin
SYMATTR Value SINE(0 1 10K)
SYMBOL res 528 -1232 R90
WINDOW 0 0 56 VBottom 2
WINDOW 3 32 56 VTop 2
SYMATTR InstName R1
SYMATTR Value 10K
SYMBOL res 640 -1232 R90
WINDOW 0 0 56 VBottom 2
WINDOW 3 32 56 VTop 2
SYMATTR InstName R2
SYMATTR Value 10K
SYMBOL AD820 576 -1168 R0
SYMATTR InstName U1
FLAG 144 -1168 +V
FLAG 432 -976 0
FLAG 704 -1104 OUT
FLAG 576 -1136 +V
FLAG 576 -1072 -V
TEXT 688 -1000 Left 2 !.tran 1m

NOTE: Op-amp power pins connected via FLAG net labels (+V, -V),
      not by explicit wires -- this is the standard LTSpice pattern.

DIRECTION MAPPING (schemdraw -> LTSpice)
-----------------------------------------
schemdraw right (dx>0) -> LTSpice R90
schemdraw left  (dx<0) -> LTSpice R270
schemdraw up    (dy>0) -> LTSpice R0   (Y-axis inverted!)
schemdraw down  (dy<0) -> LTSpice R180

VALUE FORMATTING (SI prefixes)
-------------------------------
1e12 -> 1T       1e-3  -> 1m
1e9  -> 1G       1e-6  -> 1u
1e6  -> 1Meg     1e-9  -> 1n
1e3  -> 1k       1e-12 -> 1p
1    -> 1        1e-15 -> 1f

ASC FILE ORDERING CONVENTION
------------------------------
Standard element ordering in real LTSpice files:
  1. Version header
  2. SHEET declaration
  3. All WIRE lines
  4. All FLAG lines (ground and net labels)
  5. IOPIN lines (if any, after corresponding FLAG)
  6. SYMBOL + WINDOW + SYMATTR blocks (each component together)
  7. TEXT lines (directives and comments)
```

## LTspice advanced component reference

<a id="circuit-ltspice-advanced-components"></a>

```
LTSpice Advanced Components & Patterns (from Educational Examples)
===================================================================

TRANSISTORS
-----------
BJT NPN:
  SYMBOL npn x y R0           (or NPN -- case-insensitive)
  SYMATTR InstName Q1
  SYMATTR Value 2N3904

BJT PNP:
  SYMBOL pnp x y M180         (common: M180 for standard orientation)
  WINDOW 0 60 68 Left 2
  WINDOW 3 64 28 Left 2
  SYMATTR InstName Q2
  SYMATTR Value 2N3906

JFET N-channel:
  SYMBOL njf x y R0           (also NJF)
  SYMATTR InstName J1
  SYMATTR Value 2N5484

MOSFET:
  SYMBOL nmos x y R0          (N-channel enhancement)
  SYMBOL pmos x y M180        (P-channel)
  SYMATTR InstName M1
  SYMATTR Value IRFP240

IGBT:
  SYMBOL misc\\nigbt x y R0
  SYMATTR InstName Z1
  SYMATTR Prefix Z

DIODES
------
Standard diode:
  SYMBOL diode x y R0
  SYMATTR InstName D1
  SYMATTR Value 1N4148

  R180 (flipped): WINDOW 0 24 72 Left 2 / WINDOW 3 24 0 Left 2

Zener diode:
  SYMBOL zener x y M180
  SYMATTR InstName D1
  SYMATTR Value 6.3V

Schottky diode:
  SYMBOL schottky x y R0
  SYMATTR InstName D1
  SYMATTR Value 1N5818

OP-AMPS
-------
Ideal (no supply pins):
  SYMBOL OPAMPS\\OPAMP x y R0
  SYMATTR InstName U1
  (Requires .include opamp.sub)

UniversalOpamp2 (built-in, no supply needed):
  SYMBOL opamps\\UniversalOpamp2 x y R0
  SYMATTR InstName U1

Real models (need +V/-V supply via FLAG nets):
  SYMBOL opamps\\LT1001 x y R0        (or Opamps\\LT1001)
  SYMBOL AD820 x y R0
  Power: FLAG sx sy+32 +V / FLAG sx sy-32 -V  (offsets from symbol pos)

BEHAVIORAL / CONTROLLED SOURCES
---------------------------------
VCVS (Voltage-Controlled Voltage Source):
  SYMBOL e x y R0
  SYMATTR InstName E1
  SYMATTR Value 1                    (gain)
  SYMATTR Value Laplace=1./(1+.0005*s)**3   (Laplace transfer function)

VCCS (Voltage-Controlled Current Source):
  SYMBOL g x y R0
  SYMATTR InstName G1
  SYMATTR Value {2/R1}               (parameterized)

Behavioral Voltage Source:
  SYMBOL bv x y R0
  SYMATTR InstName B1
  SYMATTR Value V=exp(time-7)        (arbitrary expression)

Behavioral Current Source:
  SYMBOL bi2 x y R0
  SYMATTR InstName B1
  SYMATTR Value I={Cjo}/(1+max(V(bias),-.5*{Vj})/{Vj})**{m}

COUPLED INDUCTORS / TRANSFORMERS
---------------------------------
Use ind2 elements with K coupling statement:

  SYMBOL ind2 x1 y1 R0              (winding 1)
  SYMATTR InstName L1
  SYMATTR Value 100u
  SYMATTR Type ind

  SYMBOL ind2 x2 y2 M0              (winding 2)
  SYMATTR InstName L2
  SYMATTR Value 900u
  SYMATTR Type ind

  TEXT x y alignment 2 !K1 L1 L2 1   (coupling coefficient = 1)

For 3+ windings: !K1 L1 L2 L3 1

CRYSTAL:
  SYMBOL MISC\\XTAL x y R90
  SYMATTR InstName Y1
  SYMATTR Value 0.25p
  SYMATTR SpiceLine Rser=0.1 Lser=0.001 Cpar=5e-011

SWITCHES
--------
Voltage-controlled switch:
  SYMBOL sw x y M180
  SYMATTR InstName S1
  SYMATTR Value MYSW
  TEXT x y Left 2 !.model MYSW SW(Ron=1 Roff=1Meg Vt=.5 Vh=-.4)

DIGITAL ELEMENTS
-----------------
  SYMBOL Digital\\XOR x y R0
  SYMBOL DIGITAL\\SCHMTBUF x y R0
  SYMBOL Digital\\dflop x y M0

SPECIAL / MISC
--------------
  SYMBOL Misc\\jumper x y R0          (test point / jumper)
  SYMBOL SpecialFunctions\\sample x y R0  (sample & hold)
  SYMBOL SpecialFunctions\\MODULATE x y R0 (modulator)
  SYMBOL POWERPRODUCTS\\LT1184F x y R0    (power IC)

GRAPHICAL ANNOTATIONS
----------------------
  RECTANGLE Normal x1 y1 x2 y2 2     (box around circuit section)
  LINE Normal x1 y1 x2 y2            (annotation line)
  DATAFLAG x y ""                     (data readout point)

ADVANCED ANALYSIS DIRECTIVES
------------------------------
DC sweep:
  .dc V1 0 10 1m
  .dc V1 0 15 10m I1 20u 100u 20u     (nested sweep)

Noise analysis:
  .noise V(out) V1 oct 10 1K 100K

S-parameter / Network analysis:
  .net V(out) V1 Rout=50 Rin=50
  .net I(Rout) V4                      (auto-detect impedances)

Fourier analysis:
  .four 1K V(out)

Monte Carlo:
  .step param X 0 20 1                 (cycle MC runs)
  Component: {mc(1n, tol)}             (random tolerance)
  .param tol=.05

Parameter sweep:
  .step param X list 1 10 100 1K       (list values)
  .step oct param V 1m 1.44 2          (octave steps)

Simulation options:
  .options method=trap
  .options maxstep=.0125u
  .options plotwinsize=0 numdgt=15

Include / Model:
  .include opamp.sub
  .model NP NPN(BF=125 Cje=.5p Cjc=.5p Rb=500)
  .model PN LPNP(BF=25 Cje=.3p Cjc=1.5p Rb=250)

Subcircuit (.subckt) in TEXT directive:
  TEXT x y Left 2 !.subckt MYCOMP T1 T2\\n...\\n.ends MYCOMP

PARAMETERIZATION
-----------------
Use {} braces for expressions:
  SYMATTR Value {6*R}
  SYMATTR Value {mc(1n, tol)}
  .param f0=1k Q=0.5
  .param L1=R1*Q/(2*pi*f0)
  .param C1=1/(L1*(2*pi*f0)**2)

Functions: mc(val,tol), flat(x), gauss(x)

WINDOW PATTERNS BY ROTATION (comprehensive)
---------------------------------------------
R0  (vertical, default):   usually no WINDOW override needed
R90 (horizontal):
  res:  WINDOW 0 0 56 VBottom 2 / WINDOW 3 32 56 VTop 2
  cap:  WINDOW 0 0 32 VBottom 2 / WINDOW 3 32 32 VTop 2
  ind:  (uses R270 convention below)
R270 (horizontal inductor):
        WINDOW 0 32 56 VTop 2 / WINDOW 3 5 56 VBottom 2
  or:   WINDOW 0 32 56 VTop 2 / WINDOW 3 4 56 VBottom 2
M0  (mirrored vertical):
  g:    WINDOW 0 -10 9 Right 2 / WINDOW 3 -15 96 Right 2
  ind:  WINDOW 0 -2 30 Right 2 / WINDOW 3 -2 59 Right 2
M90 (mirrored horizontal):
  current: WINDOW 0 -32 40 VBottom 2 / WINDOW 3 32 40 VTop 2
M180 (mirrored, flipped):
  res:  WINDOW 0 36 76 Left 2 / WINDOW 3 36 40 Left 2
  pnp:  WINDOW 0 60 68 Left 2 / WINDOW 3 64 28 Left 2
  cap:  WINDOW 0 24 56 Left 2 / WINDOW 3 24 8 Left 2
R180 (flipped):
  current: WINDOW 0 24 80 Left 2 / WINDOW 3 24 0 Left 2
  diode:   WINDOW 0 24 72 Left 2 / WINDOW 3 24 0 Left 2

Hiding windows: WINDOW 123 0 0 Left 0 / WINDOW 39 0 0 Left 0
```

## SPICE netlist (`.cir`) format

<a id="circuit-spice-netlist-format"></a>

```
SPICE Netlist (.cir/.net/.sp) Format Reference
================================================

A SPICE netlist is a plain-text file describing a circuit.

BASIC STRUCTURE
---------------
* Title line (first line is always the title)
R1 node1 node2 value           ; Resistor
C1 node1 node2 value           ; Capacitor
L1 node1 node2 value           ; Inductor
V1 node+ node- value           ; Voltage source
I1 node+ node- value           ; Current source
D1 anode cathode model         ; Diode
Q1 C B E model                 ; BJT (Collector Base Emitter)
M1 D G S B model               ; MOSFET (Drain Gate Source Bulk)
J1 D G S model                 ; JFET
X1 node1 node2 ... subckt_name ; Subcircuit instance
.model name type(params)        ; Model definition
.subckt name node1 node2 ...    ; Subcircuit definition
.ends                           ; End subcircuit
.end                            ; End of netlist

COMPONENT LINE FORMAT
---------------------
First character determines type:
  R = Resistor       R<name> <n+> <n-> <value>
  C = Capacitor      C<name> <n+> <n-> <value> [IC=<v>]
  L = Inductor       L<name> <n+> <n-> <value> [IC=<i>]
  V = Voltage src    V<name> <n+> <n-> <value_or_spec>
  I = Current src    I<name> <n+> <n-> <value_or_spec>
  D = Diode          D<name> <anode> <cathode> <model>
  Q = BJT            Q<name> <C> <B> <E> [<S>] <model>
  M = MOSFET         M<name> <D> <G> <S> <B> <model>
  J = JFET           J<name> <D> <G> <S> <model>
  X = Subcircuit     X<name> <nodes...> <subckt_name>
  E = VCVS           E<name> <n+> <n-> <nc+> <nc-> <gain>
  F = CCCS           F<name> <n+> <n-> <vname> <gain>
  G = VCCS           G<name> <n+> <n-> <nc+> <nc-> <gain>
  H = CCVS           H<name> <n+> <n-> <vname> <gain>
  K = Coupling       K<name> L1 L2 <coefficient>
  * = Comment line

Node "0" or "GND" = ground reference.

SOURCE SPECIFICATIONS
---------------------
DC:     V1 n+ n- 5
AC:     V1 n+ n- AC 1
SINE:   V1 n+ n- SINE(offset amplitude frequency)
PULSE:  V1 n+ n- PULSE(V1 V2 Tdelay Trise Tfall Ton Tperiod)
PWL:    V1 n+ n- PWL(t1 v1 t2 v2 ...)

ANALYSIS COMMANDS
-----------------
.tran <tstop>
.ac dec <npts> <fstart> <fstop>
.dc <src> <start> <stop> <step>
.op
.noise V(out) Vin dec <npts> <fstart> <fstop>
.param <name>=<value>
.include <filename>

MAPPING: NETLIST -> ASC
========================

Component type -> ASC SYMBOL:
  R -> res
  C -> cap
  L -> ind (or ind2 for coupled)
  V -> voltage
  I -> current
  D -> diode
  Q (NPN) -> npn
  Q (PNP) -> pnp
  M (NMOS) -> nmos
  M (PMOS) -> pmos
  J (NJF) -> njf
  J (PJF) -> pjf
  X -> subcircuit symbol name (from .subckt)
  E -> e (VCVS)
  G -> g (VCCS)
  K -> TEXT directive (!K L1 L2 coefficient)

Node -> coordinate mapping strategy:
  1. Assign each unique node a (x, y) coordinate
  2. Place components between their nodes
  3. Generate WIRE elements to connect
  4. Add FLAG 0 for ground nodes
  5. Add FLAG name for labeled nodes

AUTOMATIC PLACEMENT ALGORITHM
-------------------------------
For .cir -> .asc conversion, a placement strategy is needed:

1. Parse netlist into components and nodes
2. Build adjacency graph (which nodes connect to which)
3. Assign grid positions to nodes:
   - Ground node (0) at bottom center
   - Source nodes at left
   - Output nodes at right
   - Internal nodes placed to minimize wire crossings
4. Place components between their nodes:
   - Determine orientation (R0/R90) from node positions
   - Calculate symbol position from terminal positions
5. Generate wires connecting components to nodes
6. Add FLAGS for ground and named nets
7. Add SPICE directives as TEXT elements

EXAMPLE: Simple RC netlist -> ASC
----------------------------------
Input (.cir):
  * RC Lowpass Filter
  V1 in 0 AC 1
  R1 in out 1k
  C1 out 0 1u
  .ac dec 100 1 100k
  .end

Output (.asc):
  Version 4
  SHEET 1 880 680
  WIRE 0 192 0 384
  WIRE 0 192 192 192
  WIRE 192 192 192 256
  WIRE 192 384 0 384
  FLAG 0 384 0
  FLAG 0 192 in
  FLAG 192 192 out
  IOPIN 192 192 Out
  SYMBOL voltage 0 192 R0
  SYMATTR InstName V1
  SYMATTR Value AC 1
  SYMBOL res 192 176 R90
  WINDOW 0 0 56 VBottom 2
  WINDOW 3 32 56 VTop 2
  SYMATTR InstName R1
  SYMATTR Value 1k
  SYMBOL cap 176 192 R0
  SYMATTR InstName C1
  SYMATTR Value 1u
  TEXT 0 512 Left 2 !.ac dec 100 1 100k
```
