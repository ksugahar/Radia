/*
 * mtef_common.h -- Shared constants for MTEF↔LaTeX bidirectional converter
 *
 * Single source of truth for all enum values used by both mtef2tex and tex2mtef.
 * Replaces duplicated #define constants across mtef2tex.c and tex2mtef.c.
 *
 * Compatible with both C and C++ compilation.
 */
#ifndef MTEF_COMMON_H
#define MTEF_COMMON_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================
 * MTEF Record Types (low 4 bits of tag byte)
 * ============================================================ */
enum {
    REC_END     = 0,
    REC_LINE    = 1,
    REC_CHAR    = 2,
    REC_TMPL    = 3,
    REC_PILE    = 4,
    REC_MATRIX  = 5,
    REC_EMBELL  = 6,
    REC_RULER   = 7,
    REC_FONT    = 8,
    REC_SIZE    = 9,
    REC_FULL    = 10,
    REC_SUB     = 11,
    REC_SUB2    = 12,
    REC_SYM     = 13,
    REC_SUBSYM  = 14,
};

/* ============================================================
 * Tag byte option bits (high 4 bits of tag byte)
 * ============================================================ */
enum {
    OPT_LINE_NULL   = 0x01,
    OPT_CHAR_EMBELL = 0x02,
    OPT_LINE_LSPACE = 0x04,
    OPT_NUDGE       = 0x08,
};

/* ============================================================
 * Typeface codes (stored in CHAR record, with 0x80 flag for 16-bit code)
 * ============================================================ */
/* Typeface: logical index (used after parser subtracts 0x80 from wire byte).
 * mtef2tex.c compares against these values. */
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
};

/* Typeface: wire format (logical | 0x80, indicates 16-bit char code).
 * tex2mtef.c emits these values in CHAR records. */
enum {
    TFW_TEXT     = 0x81,
    TFW_FUNCTION = 0x82,
    TFW_VAR      = 0x83,
    TFW_LCGREEK  = 0x84,
    TFW_UCGREEK  = 0x85,
    TFW_SYMBOL   = 0x86,
    TFW_VECTOR   = 0x87,
    TFW_NUMBER   = 0x88,
    TFW_EXTRA    = 0x8b,
    TFW_DISPLAY  = 0x96,

    /* Backward-compatible aliases (tex2mtef.c uses these directly) */
    TF_VAR       = 0x83,
    TF_FUNC      = 0x82,
    TF_EXTRA     = 0x8b,
    TF_DISPLAY   = 0x96,
};

/* ============================================================
 * Template selectors (0..48)
 *
 * Single unified enum replacing tmXXX (mtef2tex.c) and TM_XXX (tex2mtef.c).
 * Both prefixes are provided as aliases for backward compatibility.
 * ============================================================ */
enum {
    /* Fences (0-12) */
    tmANGLE   = 0,   TM_ANGLE   = 0,
    tmPAREN   = 1,   TM_PAREN   = 1,
    tmBRACE   = 2,   TM_BRACE   = 2,
    tmBRACK   = 3,   TM_BRACK   = 3,
    tmBAR     = 4,   TM_BAR     = 4,
    tmDBAR    = 5,   TM_DBAR    = 5,
    tmFLOOR   = 6,   TM_FLOOR   = 6,
    tmCEIL    = 7,   TM_CEIL    = 7,
    tmLBLB    = 8,   TM_LBLB    = 8,
    tmRBRB    = 9,   TM_RBRB    = 9,
    tmRBLB    = 10,  TM_RBLB    = 10,
    tmLBRP    = 11,  TM_LBRP    = 11,
    tmLPRB    = 12,  TM_LPRB    = 12,

    /* Structural (13-15) */
    tmROOT    = 13,  TM_ROOT    = 13,
    tmFRACT   = 14,  TM_FRACT   = 14,
    tmSCRIPT  = 15,  TM_SCRIPT  = 15,

    /* Decorations (16-20) */
    tmUBAR    = 16,  TM_UBAR    = 16,
    tmOBAR    = 17,  TM_OBAR    = 17,
    tmLARROW  = 18,  TM_LARROW  = 18,
    tmRARROW  = 19,  TM_RARROW  = 19,
    tmBARROW  = 20,  TM_BARROW  = 20,

    /* Integrals (21-26) */
    tmSINT    = 21,  TM_SINT    = 21,
    tmDINT    = 22,  TM_DINT    = 22,
    tmTINT    = 23,  TM_TINT    = 23,
    tmSSINT   = 24,  TM_SSINT   = 24,
    tmDSINT   = 25,  TM_DSINT   = 25,
    tmTSINT   = 26,  TM_TSINT   = 26,

    /* Horizontal braces (27-28) */
    tmUHBRACE = 27,  TM_UHBRACE = 27,
    tmLHBRACE = 28,  TM_LHBRACE = 28,

    /* BigOps (29-38) */
    tmSUM     = 29,  TM_SUM     = 29,
    tmISUM    = 30,  TM_ISUM    = 30,
    tmPROD    = 31,  TM_PRODUCT = 31,  TM_IPROD = 32,
    tmIPROD   = 32,  TM_IPRODUCT = 32,
    tmCOPROD  = 33,  TM_COPRODUCT = 33,
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
 * MTEF binary wire format byte constants
 * ============================================================ */
enum {
    MTEF_HDR_VER     = 0x03,  /* MTEF version 3 */
    MTEF_HDR_PLAT    = 0x01,  /* platform: Windows */
    MTEF_HDR_PROD    = 0x01,  /* product: Equation Editor */
    MTEF_HDR_PRODMAJ = 0x03,  /* product major version */
    MTEF_HDR_PRODMIN = 0x0a,  /* product minor version */

    REC_END_B    = 0x00,
    REC_LINE_B   = 0x01,
    SIZE_FULL_B  = 0x0a,
    SIZE_SUB_B   = 0x0b,
    SIZE_SYM_B   = 0x0d,
    NULL_LINE_B  = 0x11,
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
#endif

#endif /* MTEF_COMMON_H */
