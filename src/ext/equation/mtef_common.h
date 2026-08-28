/*
 * mtef_common.h -- shared constants for the historical equation tree
 *
 * The namespace and filename are retained to avoid a mechanical rewrite of
 * the mature TeX layout engine. Binary MTEF record/header constants are gone:
 * no MTEF or .eqn reader/writer remains in the supported source.
 */
#ifndef MTEF_COMMON_H
#define MTEF_COMMON_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================
 * Typeface codes (stored in CHAR record, with 0x80 flag for 16-bit code)
 * ============================================================ */
/* Typeface categories used by the TeX parser and renderers. */
enum {
    TF_TEXT     = 1,
    TF_FUNCTION = 2,
    TF_VARIABLE = 3,
    TF_LCGREEK  = 4,
    TF_UCGREEK  = 5,
    TF_SYMBOL   = 6,
    TF_VECTOR   = 7,
    TF_NUMBER   = 8,
    TF_USER1    = 9,
    TF_USER2    = 10,
    TF_MTEXTRA  = 11,
    /* Upright bold, \mathbf -- a matrix name, NOT a vector.  The lab sets a
     * vector as \vec\bm, which is BOLD ITALIC, and that is TF_VECTOR above.
     * The two were one code until the vector rule was written down, so
     * \mathbf{A} and \bm{A} drew the same and one of them was wrong.
     * MTEF has no code of its own for it and folds it onto TF_VECTOR. */
    TF_USER3    = 12,
    TF_DISPLAY  = 22,
};

/* ============================================================
 * Template selectors (0..48)
 *
 * Values are retained as an internal structural-tree ABI.
 * ============================================================ */
enum {
    /* Fences (0-12) */
    tmANGLE   = 0,
    tmPAREN   = 1,
    tmBRACE   = 2,
    tmBRACK   = 3,
    tmBAR     = 4,
    tmDBAR    = 5,
    tmFLOOR   = 6,
    tmCEIL    = 7,
    tmLBLB    = 8,
    tmRBRB    = 9,
    tmRBLB    = 10,
    tmLBRP    = 11,
    tmLPRB    = 12,

    /* Structural (13-15) */
    tmROOT    = 13,
    tmFRACT   = 14,
    tmSCRIPT  = 15,

    /* Decorations (16-20) */
    tmUBAR    = 16,
    tmOBAR    = 17,
    tmLARROW  = 18,
    tmRARROW  = 19,
    tmBARROW  = 20,

    /* Integrals (21-26) */
    tmSINT    = 21,
    tmDINT    = 22,
    tmTINT    = 23,
    tmSSINT   = 24,
    tmDSINT   = 25,
    tmTSINT   = 26,

    /* Horizontal braces (27-28) */
    tmUHBRACE = 27,
    tmLHBRACE = 28,

    /* BigOps (29-38) */
    tmSUM     = 29,
    tmISUM    = 30,
    tmPROD    = 31,
    tmIPROD   = 32,
    tmCOPROD  = 33,
    tmICOPROD = 34,
    tmUNION   = 35,
    tmIUNION  = 36,
    tmINTER   = 37,
    tmIINTER  = 38,

    /* Special (39-48) */
    tmLIM     = 39,
    tmLDIV    = 40,
    tmSLFRACT = 41,
    tmINTOP   = 42,
    tmSUMOP   = 43,
    tmLSCRIPT = 44,
    tmDIRAC   = 45,
    tmUARROW  = 46,
    tmOARROW  = 47,
    tmOARC    = 48,
};

/* ============================================================
 * Embellishment types (used in EMBELL records)
 * ============================================================ */
enum {
    EM_DOT     = 2,
    EM_DDOT    = 3,
    EM_TDOT    = 4,
    EM_PRIME   = 5,
    EM_DPRIME  = 6,
    EM_BPRIME  = 7,    /* backward prime / triple prime */
    EM_TILDE   = 8,
    EM_HAT     = 9,
    EM_NOT     = 10,
    EM_RARROW  = 11,   /* right arrow over */
    EM_LARROW  = 12,   /* left arrow over */
    EM_BARROW  = 13,   /* bidirectional arrow over */
    EM_R1ARROW = 14,   /* right harpoon */
    EM_L1ARROW = 15,   /* left harpoon */
    EM_MBAR    = 16,   /* middle bar (strikethrough) */
    EM_OBAR    = 17,   /* overbar */
    EM_TPRIME  = 18,   /* triple prime */
    EM_FROWN   = 19,
    EM_SMILE   = 20,
};

/* ============================================================
 * Size types (used in SIZE records, stored in high nibble of tag)
 * ============================================================ */
enum {
    SIZETYPE_FULL   = 0,
    SIZETYPE_SUB    = 1,
    SIZETYPE_SUB2   = 2,
    SIZETYPE_SYM    = 3,
    SIZETYPE_SUBSYM = 4,
};

/* ============================================================
 * Selector classification helpers
 * ============================================================ */
static inline int mtef_is_fence(int sel) {
    return sel >= tmANGLE && sel <= tmLPRB;
}
static inline int mtef_is_integral(int sel) {
    return sel >= tmSINT && sel <= tmTSINT;
}
static inline int mtef_is_bigop(int sel) {
    return (sel >= tmSINT && sel <= tmTSINT) ||
           (sel >= tmSUM  && sel <= tmIINTER);
}
static inline int mtef_is_display_tmpl(int sel) {
    return mtef_is_bigop(sel) || sel == tmUHBRACE || sel == tmLHBRACE;
}
static inline int mtef_is_decoration(int sel) {
    return sel == tmUBAR || sel == tmOBAR ||
           sel == tmLARROW || sel == tmRARROW || sel == tmBARROW;
}

#ifdef __cplusplus
} /* extern "C" */

#include <string>

/* One code point as UTF-8.
 *
 * There were five copies of this, in the SVG, MathML, RTF, OMML and editing
 * layers, and every one of them stopped at three bytes.  That was invisible
 * for as long as nothing past U+FFFF was set -- and then the letters started
 * coming from the maths font, whose italic alphabet lives at U+1D400, and
 * every variable in every output format became three bytes of nonsense at
 * once.  Five copies of a function is five places for the same bug to hide,
 * so there is now one. */
/* Function names that take LIMITS: the subscript goes UNDER the name rather
 * than beside it, in display style.  TeX's list, from plain.tex's \limits
 * group -- \sin and \log are the other kind and keep their scripts beside.
 *
 * One predicate, consulted by the layout AND by every writer, because the
 * picture and what gets pasted must not disagree about it. */
/* Explicit space, as a fraction of an em -- 0 for anything that is not one.
 *
 * TeX's four, carried as the Unicode spaces that mean the same thing, so they
 * are ordinary characters everywhere: the tree holds one, the writers emit
 * one, Word receives one.  A dedicated node would have needed a case in every
 * writer and a visit() in the interface, to say something a character already
 * says.
 *
 * They were being DROPPED -- a\,b came out the same width as ab -- which for
 * a lab that writes "5\,mm" all day is not a small thing. */
/* The double-struck and script alphabets.
 *
 * Unicode put most of each in the Mathematical Alphanumeric block and left
 * HOLES where a letter already had a code point of its own in Letterlike
 * Symbols -- so double-struck R is U+211D and not U+1D549, and script L is
 * U+2112 and not U+1D4C1.  A font has the glyph at the assigned place and
 * nowhere else, so arithmetic on the block start alone draws nothing for
 * exactly the letters people use most: R for the reals, L for a Lagrangian.
 *
 * The model holds the real character, which is what Word wants, and what
 * reads back as \mathbb{R} on the way out. */
inline unsigned int mtef_double_struck_of(unsigned int cp) {
    switch (cp) {
        case 'C': return 0x2102; case 'H': return 0x210D;
        case 'N': return 0x2115; case 'P': return 0x2119;
        case 'Q': return 0x211A; case 'R': return 0x211D;
        case 'Z': return 0x2124;
        default: break;
    }
    if (cp >= 'A' && cp <= 'Z') return 0x1D538 + (cp - 'A');
    if (cp >= 'a' && cp <= 'z') return 0x1D552 + (cp - 'a');
    if (cp >= '0' && cp <= '9') return 0x1D7D8 + (cp - '0');
    return 0;
}

inline unsigned int mtef_script_of(unsigned int cp) {
    switch (cp) {
        case 'B': return 0x212C; case 'E': return 0x2130;
        case 'F': return 0x2131; case 'H': return 0x210B;
        case 'I': return 0x2110; case 'L': return 0x2112;
        case 'M': return 0x2133; case 'R': return 0x211B;
        case 'e': return 0x212F; case 'g': return 0x210A;
        case 'o': return 0x2134;
        default: break;
    }
    if (cp >= 'A' && cp <= 'Z') return 0x1D49C + (cp - 'A');
    if (cp >= 'a' && cp <= 'z') return 0x1D4B6 + (cp - 'a');
    return 0;
}

/* And back, for writing the command out again.  Returns 0 when the character
 * belongs to neither alphabet. */
inline unsigned int mtef_plain_of_alphabet(unsigned int cp, bool* doubleStruck) {
    for (unsigned int c = 'A'; c <= 'Z'; ++c) {
        if (mtef_double_struck_of(c) == cp) { *doubleStruck = true;  return c; }
        if (mtef_script_of(c) == cp)        { *doubleStruck = false; return c; }
    }
    for (unsigned int c = 'a'; c <= 'z'; ++c) {
        if (mtef_double_struck_of(c) == cp) { *doubleStruck = true;  return c; }
        if (mtef_script_of(c) == cp)        { *doubleStruck = false; return c; }
    }
    for (unsigned int c = '0'; c <= '9'; ++c)
        if (mtef_double_struck_of(c) == cp) { *doubleStruck = true;  return c; }
    return 0;
}

/* Bold italic -- the lab's VECTOR face.
 *
 * \vec\bm{B} is the base spelling of a vector here, so the letter under the
 * arrow is bold italic.  Unicode has no bold-italic digits (a bold digit is
 * the same glyph in both), so those come from the bold set. */
inline unsigned int mtef_bold_italic_of(unsigned int cp) {
    if (cp >= 'A' && cp <= 'Z')     return 0x1D468 + (cp - 'A');
    if (cp >= 'a' && cp <= 'z')     return 0x1D482 + (cp - 'a');
    if (cp >= 0x391 && cp <= 0x3A9) return 0x1D71C + (cp - 0x391);
    if (cp >= 0x3B1 && cp <= 0x3C9) return 0x1D736 + (cp - 0x3B1);
    if (cp == 0x2207)               return 0x1D735;   /* nabla   */
    if (cp == 0x2202)               return 0x1D74F;   /* partial */
    if (cp >= '0' && cp <= '9')     return 0x1D7CE + (cp - '0');
    return 0;
}

/* Upright bold -- \mathbf. */
inline unsigned int mtef_bold_upright_of(unsigned int cp) {
    if (cp >= 'A' && cp <= 'Z')     return 0x1D400 + (cp - 'A');
    if (cp >= 'a' && cp <= 'z')     return 0x1D41A + (cp - 'a');
    if (cp >= '0' && cp <= '9')     return 0x1D7CE + (cp - '0');
    if (cp >= 0x391 && cp <= 0x3A9) return 0x1D6A8 + (cp - 0x391);
    if (cp >= 0x3B1 && cp <= 0x3C9) return 0x1D6C2 + (cp - 0x3B1);
    return 0;
}

inline double mtef_space_em(unsigned int cp) {
    switch (cp) {
        case 0x2006: return 3.0 / 18.0;   /* \, thin        (six-per-em)  */
        case 0x205F: return 4.0 / 18.0;   /* \: medium      (medium math) */
        case 0x2005: return 5.0 / 18.0;   /* \; thick       (four-per-em) */
        case 0x2003: return 1.0;          /* \quad          (em)          */
        default:     return 0.0;
    }
}

inline bool mtef_name_takes_limits(const std::string& name) {
    static const char* kNames[] = {
        "det", "gcd", "inf", "lim", "liminf", "limsup", "max", "min",
        "Pr", "sup", "injlim", "projlim", nullptr
    };
    for (int i = 0; kNames[i]; ++i)
        if (name == kNames[i]) return true;
    return false;
}

/* Which cut of Latin Modern Roman a size wants.
 *
 * Latin Modern is drawn at eight optical sizes, and they are not the same
 * letter scaled: the small ones are wider and lighter-cut so they hold up.
 * TeX picks between them by the size being set -- the ranges are lmodern.sty's
 * -- and the difference is not decorative.  Setting everything from the 10 pt
 * cut left \sin 0.31 pt wide at 12 point, which looked like kerning inside
 * the name until the two cuts were measured side by side: lmroman12 gives
 * 14.4240, which is TeX's number exactly, and lmroman10 gives 14.7360. */
inline int mtef_optical_size(double pt) {
    if (pt < 5.5)  return 5;
    if (pt < 6.5)  return 6;
    if (pt < 7.5)  return 7;
    if (pt < 8.5)  return 8;
    if (pt < 9.5)  return 9;
    if (pt < 11.0) return 10;
    if (pt < 15.0) return 12;
    return 17;
}

/* The GDI family name for that cut: "LM Roman 12", and so on. */
inline std::string mtef_roman_face(double pt) {
    return "LM Roman " + std::to_string(mtef_optical_size(pt));
}

inline std::string mtef_utf8_of(unsigned int cp) {
    std::string s;
    if (cp < 0x80) {
        s += char(cp);
    } else if (cp < 0x800) {
        s += char(0xC0 | (cp >> 6));
        s += char(0x80 | (cp & 0x3F));
    } else if (cp < 0x10000) {
        s += char(0xE0 | (cp >> 12));
        s += char(0x80 | ((cp >> 6) & 0x3F));
        s += char(0x80 | (cp & 0x3F));
    } else {
        s += char(0xF0 | (cp >> 18));
        s += char(0x80 | ((cp >> 12) & 0x3F));
        s += char(0x80 | ((cp >> 6) & 0x3F));
        s += char(0x80 | (cp & 0x3F));
    }
    return s;
}
#endif

#endif /* MTEF_COMMON_H */
