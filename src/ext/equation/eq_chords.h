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

/* What a key press means, decided without a window.
 *
 * The window used to read the modifiers from GetKeyState -- the calling
 * thread's input state -- so the only way to exercise this was to take over
 * somebody's keyboard and actually press the keys.  That made the one layer
 * with a history of silent breakage the one layer with no test, and made
 * checking it an imposition on whoever was at the machine.
 *
 * Passing the modifiers in instead costs the window nothing and makes the
 * whole dispatch, two-step chords included, an ordinary function call. */
struct KeyState {
    std::vector<Step> pending;   /* a chord waiting for its second key */
    bool swallow_char = false;   /* this press was a chord, not typing */
};

enum class KeyResult {
    Ignored,     /* not a chord; the caller should treat it as typing */
    Consumed,    /* a command ran */
    Pending,     /* the first key of a two-step chord */
};

/* Look one key press up in the table and, if it names a command, run it. */
KeyResult press_key(class Equation& eq, KeyState& st, unsigned vk,
                    bool ctrl, bool shift, bool alt);

}  // namespace mtef

#endif  /* MTEF_EQ_CHORDS_H */
