#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>

#include "eq_chords.h"

#include "eq_edit.h"

namespace mtef {
namespace {

unsigned named_key(const std::string& name) {
    if (name == "Tab")       return VK_TAB;
    if (name == "Left")      return VK_LEFT;
    if (name == "Right")     return VK_RIGHT;
    if (name == "Up")        return VK_UP;
    if (name == "Down")      return VK_DOWN;
    if (name == "Home")      return VK_HOME;
    if (name == "End")       return VK_END;
    if (name == "Backspace") return VK_BACK;
    if (name == "Delete")    return VK_DELETE;
    if (name == "Insert")    return VK_INSERT;
    if (name == "PageUp")    return VK_PRIOR;
    if (name == "PageDown")  return VK_NEXT;
    if (name == "Escape")    return VK_ESCAPE;
    if (name == "Enter")     return VK_RETURN;
    return 0;
}

}  // namespace

bool parse_chord(const std::string& text, Chord& out) {
    size_t pos = 0;
    while (pos <= text.size()) {
        size_t comma = text.find(", ", pos);
        std::string part = text.substr(
            pos, comma == std::string::npos ? std::string::npos : comma - pos);

        Step st;
        size_t i = 0;
        for (;;) {
            size_t plus = part.find('+', i);
            /* A lone '+' at the end is the key itself, not a separator. */
            if (plus == std::string::npos || plus + 1 >= part.size()) break;
            std::string mod = part.substr(i, plus - i);
            if (mod == "Ctrl")       st.ctrl = true;
            else if (mod == "Shift") st.shift = true;
            else if (mod == "Alt")   st.alt = true;
            else return false;
            i = plus + 1;
        }
        std::string key = part.substr(i);
        if (key.empty()) return false;

        if (unsigned vk = named_key(key)) {
            st.vk = vk;
        } else if (key.size() == 1 && key[0] >= 'A' && key[0] <= 'Z') {
            /* A letter is the NAME of a key, not a character to type: "Ctrl+F"
             * is the F key, and only an explicit "Shift+" in the table asks for
             * shift.  Reading it as a character instead would make every
             * Ctrl+<letter> chord unreachable, since typing a capital F does
             * take shift. */
            st.vk = unsigned((unsigned char)key[0]);
        } else if (key.size() == 1) {
            /* Punctuation IS written as the character, so let the layout
             * decide: "[" and "{" are one physical key, told apart by shift. */
            SHORT r = VkKeyScanW(wchar_t((unsigned char)key[0]));
            if (r == -1) return false;
            st.vk = LOBYTE(r);
            if (HIBYTE(r) & 1) st.shift = true;
        } else {
            return false;
        }
        out.steps.push_back(st);

        if (comma == std::string::npos) break;
        pos = comma + 2;
    }
    return !out.steps.empty();
}

const std::vector<Chord>& chords() {
    static const std::vector<Chord> table = [] {
        std::vector<Chord> v;
        for (const Equation::Binding& b : Equation::shortcuts()) {
            Chord c;
            c.command = b.command;
            if (parse_chord(b.chord, c)) v.push_back(c);
        }
        return v;
    }();
    return table;
}

}  // namespace mtef
