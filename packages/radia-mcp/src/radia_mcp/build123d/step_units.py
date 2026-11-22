"""Declared STEP lengths and import into a millimetre coordinate system."""


_CONVERSION_LENGTH_UNITS = {
    "metre", "meter", "centimetre", "centimeter", "millimetre", "millimeter",
    "micrometre", "micrometer", "inch", "inches", "foot", "feet", "yard", "mil",
    "thou", "microinch", "mile", "mm", "cm", "m", "in", "ft",
}


def _declared_units(model):
    from OCP.StepBasic import (
        StepBasic_ConversionBasedUnitAndLengthUnit,
        StepBasic_SiUnitAndLengthUnit,
        StepBasic_SiUnitName,
    )
    from OCP.StepShape import StepShape_ShapeRepresentation

    names = {}
    for index in range(1, model.NbEntities() + 1):
        entity = model.Value(index)
        if not isinstance(entity, StepShape_ShapeRepresentation):
            continue
        context = entity.ContextOfItems()
        if not hasattr(context, "NbUnits"):
            raise ValueError(f"STEP shape representation #{index} has no declared length unit")
        declared = []
        for position in range(1, context.NbUnits() + 1):
            unit = context.UnitsValue(position)
            if isinstance(unit, StepBasic_SiUnitAndLengthUnit):
                if unit.Name() != StepBasic_SiUnitName.StepBasic_sunMetre:
                    raise ValueError(f"Unknown STEP declared length unit in representation #{index}")
                prefix = unit.Prefix().name.removeprefix("StepBasic_sp").lower() if unit.HasPrefix() else ""
                declared.append(prefix + "metre")
            elif isinstance(unit, StepBasic_ConversionBasedUnitAndLengthUnit):
                name = unit.Name().ToCString().strip() if unit.Name() is not None else ""
                if name.casefold() not in _CONVERSION_LENGTH_UNITS:
                    raise ValueError(f"Unknown STEP declared length unit: {name!r}")
                declared.append(name)
        if len(declared) != 1:
            raise ValueError(f"STEP shape representation #{index} requires one declared length unit")
        names.setdefault(declared[0].casefold(), declared[0])
    if not names:
        raise ValueError("STEP has no declared length unit in its shape representations")
    return list(names.values())


def read_step_length_units(path):
    from OCP.IFSelect import IFSelect_ReturnStatus
    from OCP.STEPControl import STEPControl_Reader

    if path.suffix.lower() not in {".step", ".stp"}:
        raise ValueError("STEP readback requires a .step or .stp file")
    reader = STEPControl_Reader()
    if reader.ReadFile(str(path)) != IFSelect_ReturnStatus.IFSelect_RetDone:
        raise ValueError("Could not read STEP length units")
    return _declared_units(reader.StepModel())


def load_step_mm(path):
    from build123d import import_step
    from OCP.Interface import Interface_Static

    names = read_step_length_units(path)

    # build123d's XCAF import reads the OCCT destination unit. Set and restore
    # it for the import. MCP readback, render and diff run in isolated workers.
    setting = "xstep.cascade.unit"
    previous = Interface_Static.CVal_s(setting)
    if not Interface_Static.SetCVal_s(setting, "MM"):
        raise RuntimeError("Could not set the OCCT STEP measurement unit to mm")
    try:
        shape = import_step(path)
    finally:
        Interface_Static.SetCVal_s(setting, previous)
    return shape, {"declared_length_units": names, "measurement_length_unit": "mm"}
