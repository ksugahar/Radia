/* Shared semantic codes used by the TeX syntax tree and renderer. */
#ifndef EQUATION_TYPES_H
#define EQUATION_TYPES_H

#include <string>

namespace eqnedit {

enum Typeface {
    TF_TEXT = 1,
    TF_FUNCTION = 2,
    TF_VARIABLE = 3,
    TF_LCGREEK = 4,
    TF_UCGREEK = 5,
    TF_SYMBOL = 6,
    TF_VECTOR = 7,
    TF_NUMBER = 8,
    TF_USER1 = 9,
    TF_USER2 = 10,
    TF_MATH_EXTRA = 11,
    /* Explicit spacing: a character that is an advance and nothing else.
     * Modelled as a character rather than its own node so the caret, undo,
     * selection, and Backspace treat it as one thing without any new cases;
     * the command it came from is kept in CharNode::latex. */
    TF_SPACE = 12,
    /* Explicit math alphabets used by the always-visible style palette.
     * Keep these distinct from TF_TEXT (\text) and TF_VARIABLE (automatic
     * math variables) so saving teaches the command the user actually chose. */
    TF_ROMAN = 13,       /* \mathrm{...} */
    TF_MATH_ITALIC = 14, /* \mathit{...} */
    TF_MATH_SANS = 15,   /* \mathsf{...} */
    TF_MATH_MONO = 16,   /* \mathtt{...} */
    TF_MATH_SCRIPT = 17, /* \mathcal{...} */
    TF_MATH_DOUBLE = 18, /* \mathbb{...} */
    TF_MATH_FRAKTUR = 19,/* \mathfrak{...} */
    TF_BOLD_SYMBOL = 20, /* \bm{...}; accepts \boldsymbol input */
    TF_DISPLAY = 0x96,
};

/* Width of an explicit space, in ems.  TeX's own values: \thinmuskip 3mu,
 * \medmuskip 4mu, \thickmuskip 5mu (18mu = 1em), and \quad = 1em. */
inline double space_width_em(const char* command) {
    if (!command) return 0.0;
    const std::string c(command);
    if (c == "\\!") return -3.0 / 18.0;
    if (c == "\\,") return 3.0 / 18.0;
    if (c == "\\:") return 4.0 / 18.0;
    if (c == "\\;") return 5.0 / 18.0;
    if (c == "\\ ") return 1.0 / 3.0;
    /* TeX's active `~` is a non-breaking interword space. */
    if (c == "~") return 1.0 / 3.0;
    if (c == "\\quad") return 1.0;
    if (c == "\\qquad") return 2.0;
    return 0.0;
}

/* Template selectors.  The values are internal identifiers; files store TeX,
 * never these numbers. */
enum TemplateSelector {
    tmANGLE = 0, tmPAREN = 1, tmBRACE = 2, tmBRACK = 3,
    tmBAR = 4, tmDBAR = 5, tmFLOOR = 6, tmCEIL = 7,
    tmLBLB = 8, tmRBRB = 9, tmRBLB = 10, tmLBRP = 11, tmLPRB = 12,
    tmROOT = 13, tmFRACT = 14, tmSCRIPT = 15,
    tmUBAR = 16, tmOBAR = 17, tmLARROW = 18, tmRARROW = 19,
    tmBARROW = 20,
    tmSINT = 21, tmDINT = 22, tmTINT = 23,
    tmSSINT = 24, tmDSINT = 25, tmTSINT = 26,
    tmUHBRACE = 27, tmLHBRACE = 28,
    tmSUM = 29, tmISUM = 30, tmPROD = 31, tmIPROD = 32,
    tmCOPROD = 33, tmICOPROD = 34, tmUNION = 35, tmIUNION = 36,
    tmINTER = 37, tmIINTER = 38,
    tmLIM = 39, tmLDIV = 40, tmSLFRACT = 41, tmINTOP = 42,
    tmSUMOP = 43, tmLSCRIPT = 44, tmDIRAC = 45,
    tmUARROW = 46, tmOARROW = 47, tmOARC = 48,
};

enum Embellishment {
    EM_DOT = 2, EM_DDOT = 3, EM_TDOT = 4,
    EM_PRIME = 5, EM_DPRIME = 6, EM_BPRIME = 7,
    EM_TILDE = 8, EM_HAT = 9, EM_NOT = 10,
    EM_RARROW = 11, EM_LARROW = 12, EM_BARROW = 13,
    EM_R1ARROW = 14, EM_L1ARROW = 15, EM_MBAR = 16,
    EM_OBAR = 17, EM_TPRIME = 18, EM_FROWN = 19, EM_SMILE = 20,
};

enum EquationSize {
    SIZETYPE_FULL = 0,
    SIZETYPE_SUB = 1,
    SIZETYPE_SUB2 = 2,
    SIZETYPE_SYM = 3,
    SIZETYPE_SUBSYM = 4,
};

}  // namespace eqnedit

#endif
