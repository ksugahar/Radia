"""Solver-neutral air-gap corrections for electrical machines."""

from __future__ import annotations


def carter_coefficient(slot_pitch, gap, slot_opening):
    """Return Carter's coefficient for a slotted air gap.

    The coefficient is the factor by which slot openings make the gap behave
    magnetically larger::

        gamma = (b_o/g)^2 / (5 + b_o/g)
        k_C = tau_s / (tau_s - gamma * g)

    Here ``tau_s`` is the slot pitch, ``g`` the physical gap, and ``b_o`` the
    slot opening.  The result is at least one for a valid machine geometry.
    """
    gamma = (slot_opening / gap) ** 2 / (5.0 + slot_opening / gap)
    return slot_pitch / (slot_pitch - gamma * gap)


def effective_air_gap(slot_pitch, gap, slot_opening):
    """Return the smooth-gap equivalent ``g_eff = k_C * g``.

    Use this in air-gap permeance, magnetizing-inductance, and no-load-flux
    calculations when the physical machine has slot openings.
    """
    return carter_coefficient(slot_pitch, gap, slot_opening) * gap


def slotted_air_gap_permeance_factor(slot_pitch, gap, slot_opening):
    """Return the slotted-to-smooth mean permeance ratio ``1 / k_C``.

    This is the companion ratio ``P_slot / P_smooth`` and is no greater than
    one for a valid machine geometry.
    """
    return 1.0 / carter_coefficient(slot_pitch, gap, slot_opening)
