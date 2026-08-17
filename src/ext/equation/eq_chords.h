/*
 * eq_chords.h -- turning the published shortcut table into key presses
 *
 * Equation::shortcuts() states the chords as text ("Ctrl+F", "Ctrl+T, S") so
 * that one table serves the window, the documentation and the user.  This is
 * the only place that text becomes virtual-key codes and modifier flags.
 *
 * It lives apart from the window so it can be tested without one.  The window
 * is a WIN32 executable that no test can drive, and a silent mistake here
 * costs the user every keyboard shortcut in the editor -- which is exactly
 * what happened once.
 */
#ifndef MTEF_EQ_CHORDS_H
#define MTEF_EQ_CHORDS_H

#include <string>
#include <vector>

namespace mtef {

/* One key press: a virtual-key code and the modifiers held down with it. */
struct Step {
    unsigned vk = 0;
    bool ctrl = false, shift = false, alt = false;

    bool operator==(const Step& o) const {
        return vk == o.vk && ctrl == o.ctrl && shift == o.shift && alt == o.alt;
    }
    bool operator!=(const Step& o) const { return !(*this == o); }
};

/* A chord is one or more presses in order: "Ctrl+T, S" is two. */
struct Chord {
    std::vector<Step> steps;
    std::string command;
};

/* Parse one chord's text.  False if it names a key this does not know. */
bool parse_chord(const std::string& text, Chord& out);

/* The whole published table, parsed once. */
const std::vector<Chord>& chords();

}  // namespace mtef

#endif  /* MTEF_EQ_CHORDS_H */
