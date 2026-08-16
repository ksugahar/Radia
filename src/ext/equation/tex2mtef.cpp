/*
 * tex2mtef.cpp -- LaTeX -> MTEF v3 binary
 *
 * Parses the LaTeX subset mtef2tex.cpp emits and writes MTEF v3 that
 * Equation Editor accepts.  Only the .eqn round trip needs this direction;
 * OMML and SVG are produced from LaTeX without going near MTEF.
 */

#include "tex2mtef.h"
#include "mtef_common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* ============================================================
 * Constants
 * ============================================================ */

/* MTEF record types */


/* Template selectors */

/* Embellishment types (MTEF spec — must match mtef2tex.cpp EMBELL_MAP indices) */
#define EM_TDOT     4
#define EM_NOT      10
#define EM_VEC      11
#define EM_BAR      16
#define EM_FROWN    19
#define EM_SMILE    20

/* Typeface bytes (tag byte includes options) */
#define TF_LCGREEK  0x84
#define TF_UCGREEK  0x85
#define TF_EXTRA    0x8b

/* ============================================================
 * ByteBuffer — growing byte array
 * ============================================================ */

typedef struct {
    uint8_t *data;
    int      len;
    int      cap;
} ByteBuf;

static void bb_init(ByteBuf *bb) {
    bb->data = NULL; bb->len = 0; bb->cap = 0;
}
static void bb_ensure(ByteBuf *bb, int need) {
    if (bb->len + need > bb->cap) {
        int newcap = bb->cap ? bb->cap * 2 : 256;
        while (newcap < bb->len + need) newcap *= 2;
        bb->data = (uint8_t *)realloc(bb->data, newcap);
        bb->cap = newcap;
    }
}
static void bb_byte(ByteBuf *bb, uint8_t b) {
    bb_ensure(bb, 1); bb->data[bb->len++] = b;
}
static void bb_bytes(ByteBuf *bb, const uint8_t *p, int n) {
    bb_ensure(bb, n); memcpy(bb->data + bb->len, p, n); bb->len += n;
}
static void bb_append(ByteBuf *bb, const ByteBuf *src) {
    if (src->len > 0) bb_bytes(bb, src->data, src->len);
}
static void bb_free(ByteBuf *bb) {
    free(bb->data); bb->data = NULL; bb->len = bb->cap = 0;
}
/* Detach ownership */
static uint8_t *bb_detach(ByteBuf *bb, int *outLen) {
    uint8_t *p = bb->data;
    if (outLen) *outLen = bb->len;
    bb->data = NULL; bb->len = bb->cap = 0;
    return p;
}

/* ============================================================
 * MTEF builder helpers (ported from generate_db.py)
 * ============================================================ */

static void emit_char_rec(ByteBuf *bb, uint8_t tag, uint8_t tf, uint16_t code) {
    bb_byte(bb, tag);
    bb_byte(bb, tf);
    bb_byte(bb, (uint8_t)(code & 0xFF));
    bb_byte(bb, (uint8_t)((code >> 8) & 0xFF));
}

/* char_italic: 12 83 XX 00 (native EQNEDT32 uses tag=0x12 for italic vars) */
static void emit_italic(ByteBuf *bb, char c) {
    emit_char_rec(bb, 0x12, TF_VAR, (uint8_t)c);
}
/* char_number: 02 88 XX 00 */
static void emit_number(ByteBuf *bb, char c) {
    emit_char_rec(bb, 0x02, TFW_NUMBER, (uint8_t)c);
}
/* char_symbol: 02 86 lo hi */
static void emit_symbol(ByteBuf *bb, uint16_t code) {
    emit_char_rec(bb, 0x02, TFW_SYMBOL, code);
}
/* char_text: 02 81 XX 00 */
static void emit_text_char(ByteBuf *bb, char c) {
    emit_char_rec(bb, 0x02, TFW_TEXT, (uint8_t)c);
}
/* char_function: 02 82 XX 00 */
static void emit_func_char(ByteBuf *bb, char c) {
    emit_char_rec(bb, 0x02, TF_FUNC, (uint8_t)c);
}
/* char_greek_lc: 02 84 lo hi (Greek typeface, NOT italic)
 * Native EQNEDT32 uses tf=0x84 (TF_GREEK_LC) for lowercase Greek. */
static void emit_greek_lc(ByteBuf *bb, uint16_t code) {
    emit_char_rec(bb, 0x02, 0x84, code);
}
/* char_greek_uc: 02 85 lo hi (Greek UC typeface)
 * Native EQNEDT32 uses tf=0x85 for uppercase Greek. */
static void emit_greek_uc(ByteBuf *bb, uint16_t code) {
    emit_char_rec(bb, 0x02, 0x85, code);
}
/* char_display: 02 96 lo hi */
static void emit_display(ByteBuf *bb, uint16_t code) {
    emit_char_rec(bb, 0x02, TF_DISPLAY, code);
}
/* char_extra_math: 02 8b XX 00 */
static void emit_extra_math(ByteBuf *bb, uint8_t code) {
    emit_char_rec(bb, 0x02, TF_EXTRA, code);
}
/* char_vector: 02 87 XX 00 (tf=0x87 = TF_VECTOR|0x80 for 16-bit code) */
static void emit_vector_char(ByteBuf *bb, char c) {
    emit_char_rec(bb, 0x02, TFW_VECTOR, (uint8_t)c);
}

#define EMIT_END(bb)      bb_byte(bb, REC_END_B)
#define EMIT_LINE(bb)     bb_byte(bb, REC_LINE_B)
#define EMIT_SIZE_FULL(bb) bb_byte(bb, SIZE_FULL_B)
#define EMIT_SIZE_SUB(bb)  bb_byte(bb, SIZE_SUB_B)
#define EMIT_SIZE_SYM(bb)  bb_byte(bb, SIZE_SYM_B)
#define EMIT_NULL_LINE(bb) bb_byte(bb, NULL_LINE_B)

/* Display chars for BigOp — use TF_SYMBOL (0x86), NOT TF_DISPLAY (0x96).
 * TF_DISPLAY is for fence bracket display only. */
static void emit_display_int(ByteBuf *bb) { emit_char_rec(bb, 0x02, TFW_SYMBOL, 0x222B); }
static void emit_display_sum(ByteBuf *bb) { emit_char_rec(bb, 0x02, TFW_SYMBOL, 0x2211); }
static void emit_display_prod(ByteBuf *bb) { emit_char_rec(bb, 0x02, TFW_SYMBOL, 0x220F); }

/* Fence display character pairs */
static void fence_display_chars(int selector, uint16_t *left, uint16_t *right) {
    switch (selector) {
    case TM_ANGLE: *left = 0x2329; *right = 0x232A; break;
    case TM_PAREN: *left = 0x0028; *right = 0x0029; break;
    case TM_BRACE: *left = 0x007B; *right = 0x007D; break;
    case TM_BRACK: *left = 0x005B; *right = 0x005D; break;
    case TM_BAR:   *left = 0x007C; *right = 0x007C; break;
    case TM_DBAR:  *left = 0xEC09; *right = 0xEC0A; break;
    case TM_FLOOR: *left = 0xF8F0; *right = 0xF8FB; break;
    case TM_CEIL:  *left = 0xF8EE; *right = 0xF8F9; break;
    default:       *left = 0x0028; *right = 0x0029; break;
    }
}

/* ============================================================
 * LaTeX → Unicode reverse lookup table
 * Sorted by command string for binary search.
 * ============================================================ */

typedef struct { const char *cmd; uint16_t code; } CmdEntry;

static const CmdEntry LATEX_TO_UNICODE[] = {
    {"\\Delta", 0x0394}, {"\\Gamma", 0x0393}, {"\\Lambda", 0x039B},
    {"\\Leftarrow", 0x21D0}, {"\\Leftrightarrow", 0x21D4},
    {"\\Longleftarrow", 0x27F8}, {"\\Longleftrightarrow", 0x27FA},
    {"\\Longrightarrow", 0x27F9},
    {"\\Omega", 0x03A9}, {"\\Phi", 0x03A6}, {"\\Pi", 0x03A0},
    {"\\Psi", 0x03A8}, {"\\Rightarrow", 0x21D2},
    {"\\Sigma", 0x03A3}, {"\\Theta", 0x0398},
    {"\\Uparrow", 0x21D1}, {"\\Upsilon", 0x03A5}, {"\\Xi", 0x039E},
    {"\\aleph", 0x2135}, {"\\alpha", 0x03B1}, {"\\approx", 0x2248},
    {"\\ast", 0x2217},
    {"\\because", 0x2235}, {"\\beta", 0x03B2}, {"\\bigcap", 0x22C2},
    {"\\bigcup", 0x22C3}, {"\\bigtriangledown", 0x25BD},
    {"\\bullet", 0x2022},
    {"\\cap", 0x2229}, {"\\cdot", 0x22C5}, {"\\cdots", 0x2026},  /* EQNEDT32 uses U+2026 (HORIZONTAL ELLIPSIS), not U+22EF */
    {"\\chi", 0x03C7}, {"\\clubsuit", 0x2663}, {"\\cong", 0x2245},
    {"\\coprod", 0x2210}, {"\\cup", 0x222A},
    {"\\dag", 0x2020}, {"\\dagger", 0x2020}, {"\\ddagger", 0x2021},
    {"\\ddots", 0x22F1},
    {"\\delta", 0x03B4}, {"\\diamond", 0x22C4},
    {"\\diamondsuit", 0x2666}, {"\\div", 0x00F7},
    {"\\downarrow", 0x2193},
    {"\\ell", 0x2113}, {"\\emptyset", 0x2205},
    {"\\epsilon", 0x03B5}, {"\\equiv", 0x2261}, {"\\eta", 0x03B7},
    {"\\exists", 0x2203},
    {"\\forall", 0x2200}, {"\\frown", 0x2322},
    {"\\gamma", 0x03B3}, {"\\geq", 0x2265}, {"\\gg", 0x226B},
    {"\\heartsuit", 0x2665},
    {"\\hookrightarrow", 0x21AA}, {"\\hookleftarrow", 0x21A9},
    {"\\imath", 0x0131}, {"\\in", 0x2208}, {"\\infty", 0x221E},
    {"\\iota", 0x03B9},
    {"\\kappa", 0x03BA},
    {"\\lambda", 0x03BB}, {"\\langle", 0x2329},
    {"\\lceil", 0x2308}, {"\\ldots", 0x2026},
    {"\\leftarrow", 0x2190}, {"\\leftharpoondown", 0x21BD},
    {"\\leftharpoonup", 0x21BC}, {"\\leftrightarrow", 0x2194},
    {"\\leq", 0x2264}, {"\\lfloor", 0x230A}, {"\\ll", 0x226A},
    {"\\mapsto", 0x21A6}, {"\\mid", 0x2223}, {"\\models", 0x22A8},
    {"\\mp", 0x2213}, {"\\mu", 0x03BC},
    {"\\nabla", 0x2207}, {"\\nearrow", 0x2197}, {"\\neg", 0x00AC},
    {"\\neq", 0x2260}, {"\\nexists", 0x2204}, {"\\ni", 0x220B},
    {"\\not", 0x0338}, {"\\notin", 0x2209}, {"\\nu", 0x03BD},
    {"\\nwarrow", 0x2196},
    {"\\odot", 0x2299}, {"\\omega", 0x03C9}, {"\\ominus", 0x2296},
    {"\\oplus", 0x2295}, {"\\otimes", 0x2297},
    {"\\parallel", 0x2225}, {"\\partial", 0x2202},
    {"\\perp", 0x22A5}, {"\\phi", 0x03C6}, {"\\pi", 0x03C0},
    {"\\pm", 0x00B1}, {"\\prec", 0x227A}, {"\\prime", 0x2032},
    {"\\propto", 0x221D}, {"\\psi", 0x03C8},
    {"\\rangle", 0x232A}, {"\\rceil", 0x2309},
    {"\\rfloor", 0x230B}, {"\\rho", 0x03C1},
    {"\\rightarrow", 0x2192}, {"\\rightharpoondown", 0x21C1},
    {"\\rightharpoonup", 0x21C0},
    {"\\searrow", 0x2198}, {"\\setminus", 0x2216},
    {"\\sigma", 0x03C3}, {"\\sim", 0x223C}, {"\\simeq", 0x2243},
    {"\\smile", 0x2323}, {"\\spadesuit", 0x2660},
    {"\\sqcap", 0x2293}, {"\\sqcup", 0x2294},
    {"\\sqsubseteq", 0x2291}, {"\\sqsupseteq", 0x2292},
    {"\\star", 0x22C6}, {"\\subset", 0x2282},
    {"\\subseteq", 0x2286}, {"\\succ", 0x227B},
    {"\\supset", 0x2283}, {"\\supseteq", 0x2287},
    {"\\surd", 0x221A}, {"\\swarrow", 0x2199},
    {"\\tau", 0x03C4}, {"\\therefore", 0x2234},
    {"\\theta", 0x03B8}, {"\\times", 0x00D7}, {"\\to", 0x2192},
    {"\\top", 0x22A4}, {"\\triangleleft", 0x25C1},
    {"\\triangleright", 0x25B7},
    {"\\uparrow", 0x2191}, {"\\updownarrow", 0x2195},
    {"\\upsilon", 0x03C5},
    {"\\varepsilon", 0x03B5},
    /* \var* Greek: U+03D1/03D5/03D6/03F0/03F1 are "SYMBOL" codepoints absent in EQNEDT32 fonts.
     * Map each to its base-letter equivalent so the glyph appears. */
    {"\\varkappa", 0x03BA},  /* ϰ U+03F0 → κ U+03BA */
    {"\\varphi",   0x03D5},  /* ϕ U+03D5 (EQNEDT32 renders correctly) */
    {"\\varpi",    0x03C0},  /* ϖ U+03D6 → π U+03C0 */
    {"\\varrho",   0x03C1},  /* ϱ U+03F1 → ρ U+03C1 */
    {"\\varsigma", 0x03C2},  /* ς U+03C2 (standard Greek block, fine as-is) */
    {"\\vartheta", 0x03B8},  /* ϑ U+03D1 → θ U+03B8 */
    {"\\vdash", 0x22A2}, {"\\vdots", 0x22EE}, {"\\vee", 0x2228},
    {"\\wedge", 0x2227}, {"\\wp", 0x2118},
    {"\\xi", 0x03BE},
    {"\\zeta", 0x03B6},
};
#define LATEX_TO_UNICODE_N (sizeof(LATEX_TO_UNICODE)/sizeof(LATEX_TO_UNICODE[0]))

static int cmp_cmd(const void *a, const void *b) {
    return strcmp(((const CmdEntry *)a)->cmd, ((const CmdEntry *)b)->cmd);
}

static int lookup_latex_unicode(const char *cmd) {
    CmdEntry key = { cmd, 0 };
    CmdEntry *found = (CmdEntry *)bsearch(&key, LATEX_TO_UNICODE,
        LATEX_TO_UNICODE_N, sizeof(CmdEntry), cmp_cmd);
    return found ? (int)found->code : -1;
}

int tex_command_to_unicode(const char *cmd) {
    return lookup_latex_unicode(cmd);
}

/* Greek LC code set */
static int is_greek_lc(uint16_t code) {
    /* All \var* letters now map to standard Greek block (03B1-03C9), so no extra cases needed */
    return (code >= 0x03B1 && code <= 0x03C9);
}
static int is_greek_uc(uint16_t code) {
    return code >= 0x0391 && code <= 0x03A9;
}

/* Commands that are binary/relational operators (for body termination) */
static int is_symbol_operator(const char *cmd) {
    static const char *ops[] = {
        "\\approx", "\\cap", "\\cdot", "\\cong", "\\cup",
        "\\equiv", "\\geq", "\\gg", "\\in", "\\leq", "\\ll",
        "\\neq", "\\oplus", "\\ominus", "\\otimes", "\\odot",
        "\\perp", "\\prec", "\\sim", "\\subset", "\\subseteq",
        "\\succ", "\\supset", "\\supseteq", "\\times",
        "\\wedge", "\\vee",
        NULL
    };
    for (int i = 0; ops[i]; i++)
        if (strcmp(cmd, ops[i]) == 0) return 1;
    return 0;
}

/* Named functions: sin, cos, etc. */
static int is_named_function(const char *name) {
    static const char *fns[] = {
        "arccos","arcsin","arctan","arg","cos","cosh","cot","coth","csc",
        "deg","det","dim","exp","gcd","hom","inf","ker","lg","lim",
        "liminf","limsup","ln","log","max","min","sec","sin","sinh",
        "sup","tan","tanh", NULL
    };
    for (int i = 0; fns[i]; i++)
        if (strcmp(name, fns[i]) == 0) return 1;
    return 0;
}

/* Named operators (display with \operatorname) */
static int is_named_operator(const char *name) {
    static const char *ops[] = {
        "curl","div","grad","rot","tr", NULL
    };
    for (int i = 0; ops[i]; i++)
        if (strcmp(name, ops[i]) == 0) return 1;
    return 0;
}

/* ============================================================
 * Tokenizer
 * ============================================================ */

enum TokenType {
    TOK_COMMAND, TOK_LBRACE, TOK_RBRACE, TOK_SUBSCRIPT,
    TOK_SUPERSCRIPT, TOK_AMPERSAND, TOK_NEWLINE, TOK_CHAR,
    TOK_SPACE, TOK_EOF
};

typedef struct {
    int  type;
    char value[128]; /* command name or char */
    int  pos;        /* position in source */
} Token;

typedef struct {
    Token *tokens;
    int    count;
    int    cap;
} TokenList;

static void tl_init(TokenList *tl) {
    tl->tokens = NULL; tl->count = 0; tl->cap = 0;
}
static void tl_push(TokenList *tl, int type, const char *val, int pos) {
    if (tl->count >= tl->cap) {
        tl->cap = tl->cap ? tl->cap * 2 : 64;
        tl->tokens = (Token *)realloc(tl->tokens, tl->cap * sizeof(Token));
    }
    Token *t = &tl->tokens[tl->count++];
    t->type = type;
    if (val) {
        int len = (int)strlen(val);
        if (len >= (int)sizeof(t->value)) len = sizeof(t->value) - 1;
        memcpy(t->value, val, len);
        t->value[len] = '\0';
    } else {
        t->value[0] = '\0';
    }
    t->pos = pos;
}
static void tl_free(TokenList *tl) {
    free(tl->tokens); tl->tokens = NULL; tl->count = tl->cap = 0;
}

/* Known command check — for longest-prefix matching in tokenizer */
static int is_known_command(const char *cmd) {
    /* Check lookup table first */
    if (lookup_latex_unicode(cmd) >= 0) return 1;
    /* Structural commands */
    static const char *known[] = {
        "\\dfrac", "\\tfrac", "\\frac", "\\sqrt", "\\left", "\\right", "\\middle",
        "\\begin", "\\end", "\\operatorname", "\\text", "\\mathbf", "\\mathop",
        "\\overset", "\\limits", "\\nolimits",
        "\\overbrace", "\\underbrace", "\\overline", "\\underline",
        "\\overleftarrow", "\\overrightarrow", "\\overleftrightarrow",
        "\\oiint", "\\oiiint",
        "\\dot", "\\ddot", "\\dddot", "\\tilde", "\\hat", "\\vec", "\\bar",
        "\\not", "\\frown", "\\smile", "\\rm", "\\symbol",
        "\\int", "\\iint", "\\iiint", "\\oint",
        "\\sum", "\\prod", "\\coprod", "\\bigcup", "\\bigcap",
        NULL
    };
    for (int i = 0; known[i]; i++)
        if (strcmp(cmd, known[i]) == 0) return 1;
    /* Named functions */
    if (cmd[0] == '\\' && is_named_function(cmd + 1)) return 1;
    if (cmd[0] == '\\' && is_named_operator(cmd + 1)) return 1;
    return 0;
}

/* Read a backslash command from position i */
static int read_command(const char *s, int i, int len, char *out, int outmax) {
    int start = i;
    i++; /* skip backslash */
    if (i >= len) { out[0] = '\\'; out[1] = '\0'; return 1; }

    /* \\  (double backslash = newline) */
    if (s[i] == '\\') { strcpy(out, "\\\\"); return 2; }

    /* Escaped special chars: \{ \} \# \$ \% \& \_ \  — emit as CHAR */
    if (s[i] == '{' || s[i] == '}' || s[i] == '#' || s[i] == '$' ||
        s[i] == '%' || s[i] == '&' || s[i] == '_' || s[i] == ' ') {
        out[0] = s[i]; out[1] = '\0'; /* just the char, no backslash */
        return 2; /* consumed 2 source chars */
    }

    /* single special char commands: \, \; \! \| etc. */
    if (!isalpha((unsigned char)s[i])) {
        out[0] = '\\'; out[1] = s[i]; out[2] = '\0';
        return 2;
    }

    /* alpha command: \cmd */
    int j = 0;
    out[j++] = '\\';
    while (i < len && isalpha((unsigned char)s[i]) && j < outmax - 1)
        out[j++] = s[i++];
    out[j] = '\0';
    return i - start;
}

static void tokenize(const char *latex, TokenList *tl) {
    int len = (int)strlen(latex);
    int i = 0;
    char cmd[128];

    tl_init(tl);

    while (i < len) {
        char c = latex[i];

        if (c == '{')  { tl_push(tl, TOK_LBRACE, "{", i); i++; }
        else if (c == '}')  { tl_push(tl, TOK_RBRACE, "}", i); i++; }
        else if (c == '_')  { tl_push(tl, TOK_SUBSCRIPT, "_", i); i++; }
        else if (c == '^')  { tl_push(tl, TOK_SUPERSCRIPT, "^", i); i++; }
        else if (c == '&')  { tl_push(tl, TOK_AMPERSAND, "&", i); i++; }
        else if (c == '\\') {
            int consumed = read_command(latex, i, len, cmd, sizeof(cmd));
            if (strcmp(cmd, "\\\\") == 0) {
                tl_push(tl, TOK_NEWLINE, "\\\\", i);
            } else if (strcmp(cmd, "\\quad") == 0 || strcmp(cmd, "\\qquad") == 0 ||
                       strcmp(cmd, "\\,") == 0 || strcmp(cmd, "\\;") == 0 ||
                       strcmp(cmd, "\\:") == 0 || strcmp(cmd, "\\!") == 0 ||
                       strcmp(cmd, "\\thinspace") == 0 || strcmp(cmd, "\\medspace") == 0 ||
                       strcmp(cmd, "\\thickspace") == 0 || strcmp(cmd, "\\enspace") == 0) {
                /* Spacing commands — silently drop (MTEF has no explicit spacing) */
            } else if (cmd[0] != '\\') {
                /* Escaped char (e.g., \{ → '{') — emit as CHAR */
                tl_push(tl, TOK_CHAR, cmd, i);
            } else if (!is_known_command(cmd) && strlen(cmd) > 2) {
                /* Unknown command — try longest prefix matching */
                /* e.g., \logp → \log + p, \sincos → \sin + cos */
                char prefix[128];
                int best_len = 0;
                for (int end = (int)strlen(cmd) - 1; end > 1; end--) {
                    strncpy(prefix, cmd, end); prefix[end] = '\0';
                    if (is_known_command(prefix)) { best_len = end; break; }
                }
                if (best_len > 0) {
                    strncpy(prefix, cmd, best_len); prefix[best_len] = '\0';
                    tl_push(tl, TOK_COMMAND, prefix, i);
                    /* Remaining chars as individual CHAR tokens */
                    for (int k = best_len; cmd[k]; k++) {
                        char buf[2] = { cmd[k], '\0' };
                        tl_push(tl, TOK_CHAR, buf, i + k);
                    }
                } else {
                    tl_push(tl, TOK_COMMAND, cmd, i);
                }
            } else {
                tl_push(tl, TOK_COMMAND, cmd, i);
            }
            i += consumed;
        }
        else if (c == ' ' || c == '\t') {
            /* skip spaces (LaTeX math ignores them) */
            i++;
        }
        else if (c == '\n' || c == '\r') {
            /* skip newlines in source */
            i++;
        }
        else {
            char buf[2] = { c, '\0' };
            tl_push(tl, TOK_CHAR, buf, i);
            i++;
        }
    }

    tl_push(tl, TOK_EOF, "", i);
}

/* ============================================================
 * AST Node types
 * ============================================================ */

typedef enum {
    ND_CHAR, ND_SYMBOL, ND_RM, ND_FRAC, ND_SQRT,
    ND_SCRIPT, ND_FENCE, ND_DIRAC, ND_INTEGRAL, ND_BIGOP,
    ND_DECORATION, ND_BRACE_DECO, ND_ENVIRONMENT,
    ND_FUNCTION, ND_TEXT, ND_MATHBF, ND_EMBELL, ND_GROUP,
    ND_PRIME, ND_DEGREE, ND_OVERSET
} NodeType;

/* Forward declaration */
typedef struct AstNode AstNode;

/* Node list (growable) */
typedef struct {
    AstNode **items;
    int       count;
    int       cap;
} NodeList;

struct AstNode {
    NodeType type;
    union {
        struct { char ch; int typeface; uint16_t code; } chr;    /* ND_CHAR */
        struct { char latex[64]; uint16_t code; } sym;           /* ND_SYMBOL */
        struct { char ch; } rm;                                  /* ND_RM */
        struct { NodeList numer; NodeList denom; int display; } frac; /* ND_FRAC */
        struct { NodeList content; NodeList index; int has_index; } sq; /* ND_SQRT */
        struct { NodeList base; NodeList sub; NodeList sup;
                 int has_sub; int has_sup; } script;             /* ND_SCRIPT */
        struct { int selector; int variation; NodeList content; } fence; /* ND_FENCE */
        struct { NodeList bra; NodeList ket; int variation; } dirac; /* ND_DIRAC */
        struct { int selector; int variation; NodeList lower; NodeList upper;
                 NodeList body; int has_lower; int has_upper;
                 int has_limits; int is_echo; } integ;           /* ND_INTEGRAL */
        struct { int selector; NodeList lower; NodeList upper; NodeList body;
                 int has_lower; int has_upper; int has_limits; } bigop; /* ND_BIGOP */
        struct { int selector; int variation; NodeList content; } deco; /* ND_DECORATION */
        struct { int selector; NodeList content; NodeList label; } bdeco; /* ND_BRACE_DECO */
        struct { int kind; NodeList *lines; int nlines; int ncols; } env; /* ND_ENVIRONMENT 0=gathered 1=aligned 2=matrix */
        struct { char name[64]; int is_operator; } func;         /* ND_FUNCTION */
        struct { char text[256]; } text;                         /* ND_TEXT */
        struct { NodeList content; } mathbf;                     /* ND_MATHBF */
        struct { int embell_type; NodeList content; } embell;    /* ND_EMBELL */
        struct { NodeList children; } group;                     /* ND_GROUP */
        struct { int count; } prime;                             /* ND_PRIME */
        struct { int dummy; } degree;                            /* ND_DEGREE */
        struct { NodeList over; NodeList base; } overset;        /* ND_OVERSET */
    } u;
};

/* ============================================================
 * Node allocator (simple pool — no individual frees)
 * ============================================================ */

#define NODE_POOL_SIZE 4096

typedef struct NodePool {
    AstNode nodes[NODE_POOL_SIZE];
    int     used;
    struct NodePool *next;
} NodePool;

static NodePool *pool_new(void) {
    NodePool *p = (NodePool *)calloc(1, sizeof(NodePool));
    return p;
}

typedef struct {
    NodePool *first;
    NodePool *current;
} Allocator;

static void alloc_init(Allocator *a) {
    a->first = a->current = pool_new();
}
static void alloc_free(Allocator *a) {
    NodePool *p = a->first;
    while (p) { NodePool *next = p->next; free(p); p = next; }
    a->first = a->current = NULL;
}
static AstNode *alloc_node(Allocator *a, NodeType type) {
    if (a->current->used >= NODE_POOL_SIZE) {
        NodePool *np = pool_new();
        a->current->next = np;
        a->current = np;
    }
    AstNode *n = &a->current->nodes[a->current->used++];
    memset(n, 0, sizeof(*n));
    n->type = type;
    return n;
}

/* NodeList helpers */
static void nl_init(NodeList *nl) { nl->items = NULL; nl->count = nl->cap = 0; }
static void nl_push(NodeList *nl, AstNode *node) {
    if (nl->count >= nl->cap) {
        nl->cap = nl->cap ? nl->cap * 2 : 8;
        nl->items = (AstNode **)realloc(nl->items, nl->cap * sizeof(AstNode *));
    }
    nl->items[nl->count++] = node;
}
/* NOTE: NodeList items arrays are leaked (freed with the whole pool at end) */

/* ============================================================
 * Greek RM uppercase map (for {\rm X} pattern)
 * ============================================================ */

typedef struct { uint16_t code; char letter; } GreekRmEntry;
static const GreekRmEntry GREEK_RM_UC[] = {
    {0x0391, 'A'}, {0x0392, 'B'}, {0x0395, 'E'}, {0x0396, 'Z'},
    {0x0397, 'H'}, {0x0399, 'I'}, {0x039A, 'K'}, {0x039C, 'M'},
    {0x039D, 'N'}, {0x039F, 'O'}, {0x03A1, 'P'}, {0x03A4, 'T'},
    {0x03A7, 'X'},
};
#define GREEK_RM_UC_N (sizeof(GREEK_RM_UC)/sizeof(GREEK_RM_UC[0]))

static int lookup_greek_rm(char ch) {
    for (int i = 0; i < (int)GREEK_RM_UC_N; i++)
        if (GREEK_RM_UC[i].letter == ch) return GREEK_RM_UC[i].code;
    return -1;
}

/* ============================================================
 * Parser state
 * ============================================================ */

typedef struct {
    TokenList *tl;
    int        pos;
    const char *latex; /* raw source for context checks */
    Allocator  alloc;
} Parser;

/* Forward declarations */
static void parse_expr(Parser *p, NodeList *out, int stop_at);
static AstNode *parse_group(Parser *p);
static void parse_script_arg(Parser *p, NodeList *out);
static void try_attach_scripts(Parser *p, NodeList *nodes);
static int is_binary_op_node(AstNode *node);
static int is_operator_node(AstNode *node);
static void parse_standalone_script(Parser *p, NodeList *out);
static void consume_bigop_body(Parser *p, NodeList *body);

static Token *peek(Parser *p) { return &p->tl->tokens[p->pos]; }
static Token *advance(Parser *p) { return &p->tl->tokens[p->pos++]; }

static int match(Parser *p, int type, const char *val) {
    Token *t = peek(p);
    if (t->type != type) return 0;
    if (val && strcmp(t->value, val) != 0) return 0;
    p->pos++;
    return 1;
}

static Token *expect(Parser *p, int type) {
    Token *t = peek(p);
    if (t->type != type) return NULL; /* parse error — best effort */
    p->pos++;
    return t;
}

/* ============================================================
 * BigOp / Integral command sets
 * ============================================================ */

static int is_bigop_cmd(const char *cmd) {
    return strcmp(cmd, "\\int") == 0 || strcmp(cmd, "\\iint") == 0 ||
           strcmp(cmd, "\\iiint") == 0 || strcmp(cmd, "\\oint") == 0 ||
           strcmp(cmd, "\\oiint") == 0 || strcmp(cmd, "\\oiiint") == 0 ||
           strcmp(cmd, "\\sum") == 0 || strcmp(cmd, "\\prod") == 0 ||
           strcmp(cmd, "\\coprod") == 0 || strcmp(cmd, "\\bigcup") == 0 ||
           strcmp(cmd, "\\bigcap") == 0;
}

static int is_integral_cmd(const char *cmd) {
    return strcmp(cmd, "\\int") == 0 || strcmp(cmd, "\\iint") == 0 ||
           strcmp(cmd, "\\iiint") == 0 || strcmp(cmd, "\\oint") == 0 ||
           strcmp(cmd, "\\oiint") == 0 || strcmp(cmd, "\\oiiint") == 0;
}

static int integral_selector(const char *cmd) {
    if (strcmp(cmd, "\\int") == 0)    return TM_SINT;
    if (strcmp(cmd, "\\iint") == 0)   return TM_DINT;
    if (strcmp(cmd, "\\iiint") == 0)  return TM_TINT;
    if (strcmp(cmd, "\\oint") == 0)   return TM_SINT; /* var=3 */
    if (strcmp(cmd, "\\oiint") == 0)  return TM_DINT; /* var=2 */
    if (strcmp(cmd, "\\oiiint") == 0) return TM_TINT; /* var=2 */
    return TM_SINT;
}

static int integral_variation(const char *cmd) {
    if (strcmp(cmd, "\\oint") == 0)    return 3;
    if (strcmp(cmd, "\\oiint") == 0)   return 2;
    if (strcmp(cmd, "\\oiiint") == 0)  return 2;
    return 0;
}

static int sumop_selector(const char *cmd, int has_limits) {
    if (strcmp(cmd, "\\sum") == 0)    return has_limits ? TM_SUM : TM_ISUM;
    if (strcmp(cmd, "\\prod") == 0)   return has_limits ? TM_PRODUCT : TM_IPRODUCT;
    if (strcmp(cmd, "\\coprod") == 0) return has_limits ? TM_COPRODUCT : 34;
    if (strcmp(cmd, "\\bigcup") == 0) return has_limits ? 35 : 36;
    if (strcmp(cmd, "\\bigcap") == 0) return has_limits ? 37 : 38;
    return TM_ISUM;
}

/* ============================================================
 * Decoration command map
 * ============================================================ */

typedef struct { const char *cmd; int selector; int variation; } DecoEntry;
static const DecoEntry DECO_MAP[] = {
    {"\\overline", TM_OBAR, 0}, {"\\underline", TM_UBAR, 0},
    {"\\overrightarrow", TM_RARROW, 0}, {"\\overleftarrow", TM_LARROW, 0},
    {"\\overleftrightarrow", TM_BARROW, 0},
};
#define DECO_MAP_N (sizeof(DECO_MAP)/sizeof(DECO_MAP[0]))

/* Embellishment command map */
typedef struct { const char *cmd; int embell; } EmbellEntry;
static const EmbellEntry EMBELL_MAP[] = {
    {"\\dot", EM_DOT}, {"\\ddot", EM_DDOT}, {"\\dddot", EM_TDOT},
    {"\\hat", EM_HAT}, {"\\tilde", EM_TILDE}, {"\\vec", EM_VEC},
    {"\\bar", EM_BAR}, {"\\frown", EM_FROWN}, {"\\smile", EM_SMILE},
};
#define EMBELL_MAP_N (sizeof(EMBELL_MAP)/sizeof(EMBELL_MAP[0]))

/* Fence delimiter map */
typedef struct { const char *delim; int selector; } FenceDelim;
static int determine_fence(const char *left, const char *right, int *sel, int *var) {
    *var = 0;
    /* full pair matches */
    if (strcmp(left, "(") == 0 && strcmp(right, ")") == 0) { *sel = TM_PAREN; return 1; }
    if (strcmp(left, "[") == 0 && strcmp(right, "]") == 0) { *sel = TM_BRACK; return 1; }
    if ((strcmp(left, "\\{") == 0 || strcmp(left, "{") == 0) &&
        (strcmp(right, "\\}") == 0 || strcmp(right, "}") == 0)) { *sel = TM_BRACE; return 1; }
    if (strcmp(left, "|") == 0 && strcmp(right, "|") == 0) { *sel = TM_BAR; return 1; }
    if (strcmp(left, "\\|") == 0 && strcmp(right, "\\|") == 0) { *sel = TM_DBAR; return 1; }
    if (strcmp(left, "\\langle") == 0 && strcmp(right, "\\rangle") == 0) { *sel = TM_ANGLE; return 1; }
    if (strcmp(left, "\\lfloor") == 0 && strcmp(right, "\\rfloor") == 0) { *sel = TM_FLOOR; return 1; }
    if (strcmp(left, "\\lceil") == 0 && strcmp(right, "\\rceil") == 0) { *sel = TM_CEIL; return 1; }
    /* left-only (right = ".") */
    if (strcmp(right, ".") == 0 || strcmp(right, "") == 0) {
        *var = 1;
        if (strcmp(left, "(") == 0) { *sel = TM_PAREN; return 1; }
        if (strcmp(left, "[") == 0) { *sel = TM_BRACK; return 1; }
        if (strcmp(left, "\\{") == 0 || strcmp(left, "{") == 0) { *sel = TM_BRACE; return 1; }
        if (strcmp(left, "|") == 0) { *sel = TM_BAR; return 1; }
        if (strcmp(left, "\\|") == 0) { *sel = TM_DBAR; return 1; }
        if (strcmp(left, "\\langle") == 0) { *sel = TM_ANGLE; return 1; }
        *sel = TM_PAREN; return 1;
    }
    /* right-only (left = ".") */
    if (strcmp(left, ".") == 0 || strcmp(left, "") == 0) {
        *var = 2;
        if (strcmp(right, ")") == 0) { *sel = TM_PAREN; return 1; }
        if (strcmp(right, "]") == 0) { *sel = TM_BRACK; return 1; }
        if (strcmp(right, "\\}") == 0 || strcmp(right, "}") == 0) { *sel = TM_BRACE; return 1; }
        if (strcmp(right, "|") == 0) { *sel = TM_BAR; return 1; }
        if (strcmp(right, "\\|") == 0) { *sel = TM_DBAR; return 1; }
        *sel = TM_PAREN; return 1;
    }
    /* Mixed pairs */
    if (strcmp(left, "[") == 0 && strcmp(right, ")") == 0) { *sel = TM_RBLB; return 1; }
    if (strcmp(left, "(") == 0 && strcmp(right, "]") == 0) { *sel = TM_LPRB; return 1; }
    *sel = TM_PAREN; return 1;
}

/* ============================================================
 * Parser: recursive descent
 * ============================================================ */

/* Read delimiter after \left or \right */
static void read_delimiter(Parser *p, char *out, int outmax) {
    Token *t = peek(p);
    if (t->type == TOK_CHAR) {
        strncpy(out, t->value, outmax - 1); out[outmax-1] = '\0';
        advance(p);
    } else if (t->type == TOK_COMMAND) {
        strncpy(out, t->value, outmax - 1); out[outmax-1] = '\0';
        advance(p);
    } else {
        out[0] = '.'; out[1] = '\0'; /* invisible */
    }
}

/* Parse {expr} — returns nodes placed in 'out' list */
static void parse_brace_group(Parser *p, NodeList *out) {
    /* expect LBRACE already consumed */
    /* Check for {\rm X} pattern */
    if (peek(p)->type == TOK_COMMAND && strcmp(peek(p)->value, "\\rm") == 0) {
        advance(p); /* skip \rm */
        if (peek(p)->type == TOK_CHAR) {
            AstNode *n = alloc_node(&p->alloc, ND_RM);
            n->u.rm.ch = peek(p)->value[0];
            advance(p);
            expect(p, TOK_RBRACE);
            nl_push(out, n);
            return;
        }
    }
    parse_expr(p, out, TOK_RBRACE);
    match(p, TOK_RBRACE, NULL);
}

/* Parse a single group */
static AstNode *parse_group(Parser *p) {
    Token *t = peek(p);

    if (t->type == TOK_LBRACE) {
        advance(p);
        NodeList children;
        nl_init(&children);
        parse_brace_group(p, &children);
        if (children.count == 1) return children.items[0];
        AstNode *g = alloc_node(&p->alloc, ND_GROUP);
        g->u.group.children = children;
        return g;
    }

    if (t->type == TOK_COMMAND) {
        const char *cmd = t->value;

        /* Fraction */
        if (strcmp(cmd, "\\dfrac") == 0 || strcmp(cmd, "\\tfrac") == 0 ||
            strcmp(cmd, "\\frac") == 0) {
            advance(p);
            AstNode *n = alloc_node(&p->alloc, ND_FRAC);
            nl_init(&n->u.frac.numer); nl_init(&n->u.frac.denom);
            n->u.frac.display = (strcmp(cmd, "\\tfrac") != 0);
            expect(p, TOK_LBRACE);
            parse_brace_group(p, &n->u.frac.numer);
            expect(p, TOK_LBRACE);
            parse_brace_group(p, &n->u.frac.denom);
            return n;
        }

        /* Square root */
        if (strcmp(cmd, "\\sqrt") == 0) {
            advance(p);
            AstNode *n = alloc_node(&p->alloc, ND_SQRT);
            nl_init(&n->u.sq.content); nl_init(&n->u.sq.index);
            n->u.sq.has_index = 0;
            /* Optional [n] for nth root */
            if (peek(p)->type == TOK_CHAR && peek(p)->value[0] == '[') {
                advance(p); /* skip [ */
                n->u.sq.has_index = 1;
                while (peek(p)->type != TOK_EOF) {
                    if (peek(p)->type == TOK_CHAR && peek(p)->value[0] == ']') {
                        advance(p); break;
                    }
                    AstNode *ch = parse_group(p);
                    if (ch) nl_push(&n->u.sq.index, ch);
                }
            }
            expect(p, TOK_LBRACE);
            parse_brace_group(p, &n->u.sq.content);
            return n;
        }

        /* Fence: \left ... \right */
        if (strcmp(cmd, "\\left") == 0) {
            advance(p);
            char leftd[64], rightd[64];
            read_delimiter(p, leftd, sizeof(leftd));
            NodeList content;
            nl_init(&content);
            NodeList bra_content; /* for Dirac */
            nl_init(&bra_content);
            int middle_found = 0;
            /* Parse until \right */
            while (peek(p)->type != TOK_EOF) {
                if (peek(p)->type == TOK_COMMAND && strcmp(peek(p)->value, "\\right") == 0)
                    break;
                /* Handle \middle| */
                if (peek(p)->type == TOK_COMMAND && strcmp(peek(p)->value, "\\middle") == 0) {
                    advance(p);
                    read_delimiter(p, rightd, sizeof(rightd)); /* consume the | */
                    /* Save current content as bra */
                    middle_found = 1;
                    bra_content = content;
                    nl_init(&content);
                    continue;
                }
                AstNode *ch = parse_group(p);
                if (ch) {
                    nl_push(&content, ch);
                    try_attach_scripts(p, &content);
                    /* Consume body for integrals/bigops inside fence */
                    AstNode *last = content.items[content.count - 1];
                    if (last->type == ND_INTEGRAL || last->type == ND_BIGOP) {
                        NodeList *body_list = (last->type == ND_INTEGRAL) ?
                            &last->u.integ.body : &last->u.bigop.body;
                        consume_bigop_body(p, body_list);
                    }
                }
            }
            if (peek(p)->type == TOK_COMMAND && strcmp(peek(p)->value, "\\right") == 0) {
                advance(p);
                read_delimiter(p, rightd, sizeof(rightd));
            } else {
                strcpy(rightd, ".");
            }

            /* Check for Dirac notation */
            /* \left\langle ... \middle| ... \right\rangle → full Dirac */
            if (middle_found && strcmp(leftd, "\\langle") == 0 && strcmp(rightd, "\\rangle") == 0) {
                AstNode *n = alloc_node(&p->alloc, ND_DIRAC);
                n->u.dirac.bra = bra_content;
                n->u.dirac.ket = content;
                n->u.dirac.variation = 0;
                return n;
            }
            int sel, var;
            determine_fence(leftd, rightd, &sel, &var);
            /* Dirac bra: \left\langle ... \right| */
            if (strcmp(leftd, "\\langle") == 0 && strcmp(rightd, "|") == 0) {
                AstNode *n = alloc_node(&p->alloc, ND_DIRAC);
                n->u.dirac.bra = content;
                nl_init(&n->u.dirac.ket);
                n->u.dirac.variation = 1;
                return n;
            }
            /* Dirac ket: \left| ... \right\rangle */
            if (strcmp(leftd, "|") == 0 && strcmp(rightd, "\\rangle") == 0) {
                AstNode *n = alloc_node(&p->alloc, ND_DIRAC);
                nl_init(&n->u.dirac.bra);
                n->u.dirac.ket = content;
                n->u.dirac.variation = 2;
                return n;
            }
            AstNode *n = alloc_node(&p->alloc, ND_FENCE);
            n->u.fence.selector = sel;
            n->u.fence.variation = var;
            n->u.fence.content = content;
            return n;
        }

        /* Integral commands */
        if (is_integral_cmd(cmd)) {
            /* Track end position for adjacency check:
             * \int_{lo}^{hi} (no space) → limits
             * \int ^{body}   (space)    → body superscript, not limit */
            int last_end = t->pos + (int)strlen(t->value);
            advance(p);
            AstNode *n = alloc_node(&p->alloc, ND_INTEGRAL);
            n->u.integ.selector = integral_selector(cmd);
            n->u.integ.variation = integral_variation(cmd);
            nl_init(&n->u.integ.lower); nl_init(&n->u.integ.upper);
            nl_init(&n->u.integ.body);
            n->u.integ.has_lower = n->u.integ.has_upper = 0;
            n->u.integ.has_limits = 0;
            n->u.integ.is_echo = 0;

            /* Check for \limits — changes selector to display variant */
            if (peek(p)->type == TOK_COMMAND && strcmp(peek(p)->value, "\\limits") == 0) {
                advance(p); n->u.integ.has_limits = 1;
                last_end = p->tl->tokens[p->pos - 1].pos + (int)strlen(p->tl->tokens[p->pos - 1].value);
                /* Switch to display-limits selector */
                if (n->u.integ.selector == TM_SINT)  n->u.integ.selector = TM_SSINT;
                else if (n->u.integ.selector == TM_DINT) n->u.integ.selector = TM_DSINT;
                else if (n->u.integ.selector == TM_TINT) n->u.integ.selector = TM_TSINT;
            }
            /* Lower limit — only if adjacent (no space after \int) */
            if (peek(p)->type == TOK_SUBSCRIPT && peek(p)->pos == last_end) {
                advance(p); n->u.integ.has_lower = 1;
                parse_script_arg(p, &n->u.integ.lower);
                last_end = p->tl->tokens[p->pos - 1].pos + (int)strlen(p->tl->tokens[p->pos - 1].value);
            }
            /* Upper limit — only if adjacent */
            if (peek(p)->type == TOK_SUPERSCRIPT && peek(p)->pos == last_end) {
                advance(p); n->u.integ.has_upper = 1;
                parse_script_arg(p, &n->u.integ.upper);
            }
            /* Promote to display-limits selector if limits present
             * (native EQNEDT32 always uses TM_SSINT/DSINT/TSINT with limits) */
            if ((n->u.integ.has_lower || n->u.integ.has_upper) && !n->u.integ.has_limits) {
                if (n->u.integ.selector == TM_SINT)  n->u.integ.selector = TM_SSINT;
                else if (n->u.integ.selector == TM_DINT) n->u.integ.selector = TM_DSINT;
                else if (n->u.integ.selector == TM_TINT) n->u.integ.selector = TM_TSINT;
            }
            /* Body: parse until stop condition */
            /* For now, gather remaining tokens in this expression */
            /* Body consumed by parse_expr via bigop body logic */
            return n;
        }

        /* Sum/prod/bigop commands */
        if (strcmp(cmd, "\\sum") == 0 || strcmp(cmd, "\\prod") == 0 ||
            strcmp(cmd, "\\coprod") == 0 || strcmp(cmd, "\\bigcup") == 0 ||
            strcmp(cmd, "\\bigcap") == 0) {
            advance(p);
            int has_limits = 0;
            if (peek(p)->type == TOK_COMMAND && strcmp(peek(p)->value, "\\limits") == 0) {
                advance(p); has_limits = 1;
            } else if (peek(p)->type == TOK_COMMAND && strcmp(peek(p)->value, "\\nolimits") == 0) {
                advance(p); has_limits = 0;
            }
            AstNode *n = alloc_node(&p->alloc, ND_BIGOP);
            n->u.bigop.selector = sumop_selector(cmd, has_limits);
            nl_init(&n->u.bigop.lower); nl_init(&n->u.bigop.upper);
            nl_init(&n->u.bigop.body);
            n->u.bigop.has_lower = n->u.bigop.has_upper = 0;
            n->u.bigop.has_limits = has_limits;

            if (peek(p)->type == TOK_SUBSCRIPT) {
                advance(p); n->u.bigop.has_lower = 1;
                parse_script_arg(p, &n->u.bigop.lower);
            }
            if (peek(p)->type == TOK_SUPERSCRIPT) {
                advance(p); n->u.bigop.has_upper = 1;
                parse_script_arg(p, &n->u.bigop.upper);
            }
            /* Handle ^{} before _{} order (e.g. \sum\limits^{}_{i}) */
            if (!n->u.bigop.has_lower && peek(p)->type == TOK_SUBSCRIPT) {
                advance(p); n->u.bigop.has_lower = 1;
                parse_script_arg(p, &n->u.bigop.lower);
            }
            /* Drop empty limits: ^{} with no content → treat as no upper */
            if (n->u.bigop.has_upper && n->u.bigop.upper.count == 0)
                n->u.bigop.has_upper = 0;
            return n;
        }

        /* Decorations: \overline, \underline, etc. */
        for (int i = 0; i < (int)DECO_MAP_N; i++) {
            if (strcmp(cmd, DECO_MAP[i].cmd) == 0) {
                advance(p);
                AstNode *n = alloc_node(&p->alloc, ND_DECORATION);
                n->u.deco.selector = DECO_MAP[i].selector;
                n->u.deco.variation = DECO_MAP[i].variation;
                nl_init(&n->u.deco.content);
                expect(p, TOK_LBRACE);
                parse_brace_group(p, &n->u.deco.content);
                /* Double decoration detection: \overline{\overline{X}} → var=1 */
                if ((n->u.deco.selector == TM_OBAR || n->u.deco.selector == TM_UBAR) &&
                    n->u.deco.content.count == 1 &&
                    n->u.deco.content.items[0]->type == ND_DECORATION &&
                    n->u.deco.content.items[0]->u.deco.selector == n->u.deco.selector) {
                    n->u.deco.variation = 1;
                    n->u.deco.content = n->u.deco.content.items[0]->u.deco.content;
                }
                return n;
            }
        }

        /* Embellishments: \dot, \hat, \vec, etc. */
        for (int i = 0; i < (int)EMBELL_MAP_N; i++) {
            if (strcmp(cmd, EMBELL_MAP[i].cmd) == 0) {
                advance(p);
                AstNode *n = alloc_node(&p->alloc, ND_EMBELL);
                n->u.embell.embell_type = EMBELL_MAP[i].embell;
                nl_init(&n->u.embell.content);
                expect(p, TOK_LBRACE);
                parse_brace_group(p, &n->u.embell.content);
                return n;
            }
        }

        /* Overbrace / underbrace */
        if (strcmp(cmd, "\\overbrace") == 0 || strcmp(cmd, "\\underbrace") == 0) {
            int is_over = (strcmp(cmd, "\\overbrace") == 0);
            advance(p);
            AstNode *n = alloc_node(&p->alloc, ND_BRACE_DECO);
            n->u.bdeco.selector = is_over ? TM_UHBRACE : TM_LHBRACE;
            nl_init(&n->u.bdeco.content); nl_init(&n->u.bdeco.label);
            expect(p, TOK_LBRACE);
            parse_brace_group(p, &n->u.bdeco.content);
            /* label: ^{...} for overbrace, _{...} for underbrace */
            if (is_over && peek(p)->type == TOK_SUPERSCRIPT) {
                advance(p);
                parse_script_arg(p, &n->u.bdeco.label);
            } else if (!is_over && peek(p)->type == TOK_SUBSCRIPT) {
                advance(p);
                parse_script_arg(p, &n->u.bdeco.label);
            }
            return n;
        }

        /* Environments: \begin{...}...\end{...} */
        if (strcmp(cmd, "\\begin") == 0) {
            advance(p);
            expect(p, TOK_LBRACE);
            char envname[64] = {0};
            int ei = 0;
            while (peek(p)->type != TOK_RBRACE && peek(p)->type != TOK_EOF && ei < 63) {
                envname[ei++] = peek(p)->value[0];
                advance(p);
            }
            envname[ei] = '\0';
            match(p, TOK_RBRACE, NULL);

            int kind = 0; /* 0=gathered,1=aligned,2=matrix,3=bmatrix,4=pmatrix,5=vmatrix,6=Vmatrix,7=Bmatrix */
            if (strcmp(envname, "aligned") == 0) kind = 1;
            else if (strcmp(envname, "matrix") == 0) kind = 2;
            else if (strcmp(envname, "bmatrix") == 0) kind = 3;
            else if (strcmp(envname, "pmatrix") == 0) kind = 4;
            else if (strcmp(envname, "vmatrix") == 0) kind = 5;
            else if (strcmp(envname, "Vmatrix") == 0) kind = 6;
            else if (strcmp(envname, "Bmatrix") == 0) kind = 7;

            AstNode *n = alloc_node(&p->alloc, ND_ENVIRONMENT);
            n->u.env.kind = kind;
            n->u.env.lines = NULL;
            n->u.env.nlines = 0;
            n->u.env.ncols = 0;
            int lines_cap = 32;
            n->u.env.lines = (NodeList *)calloc(lines_cap, sizeof(NodeList));

            if (kind == 2) {
                /* Matrix: parse cells with & separators, rows with \\ */
                /* First pass: collect all cells flat, track row lengths */
                int row_lens[256]; /* cols per row */
                int nrows = 0, total_cells = 0, maxcols = 0;
                while (peek(p)->type != TOK_EOF && nrows < 256) {
                    if (peek(p)->type == TOK_COMMAND && strcmp(peek(p)->value, "\\end") == 0)
                        break;
                    int row_cols = 0;
                    do {
                        if (total_cells >= lines_cap) {
                            lines_cap *= 2;
                            n->u.env.lines = (NodeList *)realloc(n->u.env.lines,
                                lines_cap * sizeof(NodeList));
                        }
                        nl_init(&n->u.env.lines[total_cells]);
                        parse_expr(p, &n->u.env.lines[total_cells], TOK_AMPERSAND);
                        total_cells++;
                        row_cols++;
                        if (peek(p)->type == TOK_AMPERSAND) advance(p); else break;
                    } while (1);
                    row_lens[nrows] = row_cols;
                    nrows++;
                    if (row_cols > maxcols) maxcols = row_cols;
                    if (peek(p)->type == TOK_NEWLINE) {
                        advance(p);
                        /* Trailing \\ before \end → add empty row */
                        if (peek(p)->type == TOK_COMMAND &&
                            strcmp(peek(p)->value, "\\end") == 0) {
                            if (nrows < 256) {
                                if (total_cells >= lines_cap) {
                                    lines_cap *= 2;
                                    n->u.env.lines = (NodeList *)realloc(n->u.env.lines,
                                        lines_cap * sizeof(NodeList));
                                }
                                nl_init(&n->u.env.lines[total_cells]);
                                total_cells++;
                                row_lens[nrows] = 1;
                                nrows++;
                            }
                        }
                    }
                }
                /* Pad to nrows*maxcols: shift cells from back to front */
                int padded_total = nrows * maxcols;
                if (padded_total > lines_cap) {
                    n->u.env.lines = (NodeList *)realloc(n->u.env.lines,
                        padded_total * sizeof(NodeList));
                }
                /* Initialize empty cells */
                for (int i = total_cells; i < padded_total; i++)
                    nl_init(&n->u.env.lines[i]);
                /* Copy from back to avoid overwrite: shift rows into padded positions */
                int src = total_cells;
                for (int r = nrows - 1; r >= 0; r--) {
                    int rc = row_lens[r];
                    src -= rc;
                    int dst = r * maxcols;
                    /* Move cells from src..src+rc to dst..dst+rc */
                    if (dst != src) {
                        memmove(&n->u.env.lines[dst], &n->u.env.lines[src],
                                rc * sizeof(NodeList));
                    }
                    /* Initialize padding cells */
                    for (int c = rc; c < maxcols; c++)
                        nl_init(&n->u.env.lines[dst + c]);
                }
                n->u.env.nlines = nrows;
                n->u.env.ncols = maxcols;
                /* Consume \end{...} */
                if (peek(p)->type == TOK_COMMAND && strcmp(peek(p)->value, "\\end") == 0) {
                    advance(p);
                    expect(p, TOK_LBRACE);
                    while (peek(p)->type != TOK_RBRACE && peek(p)->type != TOK_EOF) advance(p);
                    match(p, TOK_RBRACE, NULL);
                }
            } else {
                /* gathered/aligned: parse rows separated by \\ */
                while (peek(p)->type != TOK_EOF) {
                    if (peek(p)->type == TOK_COMMAND && strcmp(peek(p)->value, "\\end") == 0) {
                        advance(p);
                        expect(p, TOK_LBRACE);
                        while (peek(p)->type != TOK_RBRACE && peek(p)->type != TOK_EOF) advance(p);
                        match(p, TOK_RBRACE, NULL);
                        break;
                    }
                    if (n->u.env.nlines >= lines_cap) {
                        lines_cap *= 2;
                        n->u.env.lines = (NodeList *)realloc(n->u.env.lines,
                            lines_cap * sizeof(NodeList));
                    }
                    nl_init(&n->u.env.lines[n->u.env.nlines]);
                    if (kind == 1) {
                        /* aligned: parse left, consume &, parse right, merge */
                        parse_expr(p, &n->u.env.lines[n->u.env.nlines], TOK_AMPERSAND);
                        if (peek(p)->type == TOK_AMPERSAND) {
                            advance(p);
                            parse_expr(p, &n->u.env.lines[n->u.env.nlines], TOK_NEWLINE);
                        }
                    } else {
                        parse_expr(p, &n->u.env.lines[n->u.env.nlines], TOK_NEWLINE);
                    }
                    n->u.env.nlines++;
                    match(p, TOK_NEWLINE, NULL);
                }
            }
            return n;
        }

        /* \text{...} */
        if (strcmp(cmd, "\\text") == 0) {
            advance(p);
            AstNode *n = alloc_node(&p->alloc, ND_TEXT);
            n->u.text.text[0] = '\0';
            if (peek(p)->type == TOK_LBRACE) {
                int bracePos = peek(p)->pos;
                advance(p); /* skip { */
                /* Find matching } */
                int depth = 1;
                while (peek(p)->type != TOK_EOF && depth > 0) {
                    if (peek(p)->type == TOK_LBRACE) depth++;
                    else if (peek(p)->type == TOK_RBRACE) {
                        depth--;
                        if (depth == 0) break;
                    }
                    advance(p);
                }
                /* Extract raw substring from source: after { to before } */
                if (p->latex) {
                    int start = bracePos + 1;
                    int end = peek(p)->pos; /* position of closing } */
                    int len = end - start;
                    if (len > 255) len = 255;
                    if (len > 0) memcpy(n->u.text.text, p->latex + start, len);
                    n->u.text.text[len < 0 ? 0 : len] = '\0';
                }
                if (peek(p)->type == TOK_RBRACE) advance(p); /* consume } */
            }
            return n;
        }

        /* \mathbf{...} */
        if (strcmp(cmd, "\\mathbf") == 0) {
            advance(p);
            AstNode *n = alloc_node(&p->alloc, ND_MATHBF);
            nl_init(&n->u.mathbf.content);
            expect(p, TOK_LBRACE);
            parse_brace_group(p, &n->u.mathbf.content);
            return n;
        }

        /* \overset{top}{base} */
        if (strcmp(cmd, "\\overset") == 0) {
            advance(p);
            AstNode *n = alloc_node(&p->alloc, ND_OVERSET);
            nl_init(&n->u.overset.over); nl_init(&n->u.overset.base);
            expect(p, TOK_LBRACE);
            parse_brace_group(p, &n->u.overset.over);
            expect(p, TOK_LBRACE);
            parse_brace_group(p, &n->u.overset.base);
            return n;
        }

        /* \operatorname{...} */
        if (strcmp(cmd, "\\operatorname") == 0) {
            advance(p);
            AstNode *n = alloc_node(&p->alloc, ND_FUNCTION);
            n->u.func.name[0] = '\0';
            n->u.func.is_operator = 1;
            if (peek(p)->type == TOK_LBRACE) {
                advance(p);
                int fi = 0;
                while (peek(p)->type != TOK_RBRACE && peek(p)->type != TOK_EOF && fi < 63) {
                    n->u.func.name[fi++] = peek(p)->value[0];
                    advance(p);
                }
                n->u.func.name[fi] = '\0';
                match(p, TOK_RBRACE, NULL);
            }
            return n;
        }

        /* \symbol{N} */
        if (strcmp(cmd, "\\symbol") == 0) {
            advance(p);
            expect(p, TOK_LBRACE);
            int code = 0;
            while (peek(p)->type != TOK_RBRACE && peek(p)->type != TOK_EOF) {
                if (peek(p)->type == TOK_CHAR && isdigit((unsigned char)peek(p)->value[0]))
                    code = code * 10 + (peek(p)->value[0] - '0');
                advance(p);
            }
            match(p, TOK_RBRACE, NULL);
            AstNode *n = alloc_node(&p->alloc, ND_SYMBOL);
            n->u.sym.code = (uint16_t)code;
            n->u.sym.latex[0] = '\0';
            return n;
        }

        /* \not — negation slash (prefix for next symbol) */
        if (strcmp(cmd, "\\not") == 0) {
            advance(p);
            /* Next symbol should be modified — simplified: just emit ̸ */
            AstNode *n = alloc_node(&p->alloc, ND_SYMBOL);
            n->u.sym.code = 0x0338;
            strcpy(n->u.sym.latex, "\\not");
            return n;
        }

        /* \mathop{sym}\limits_{lo}^{hi} */
        if (strcmp(cmd, "\\mathop") == 0) {
            advance(p);
            /* Parse the {symbol} */
            expect(p, TOK_LBRACE);
            NodeList content;
            nl_init(&content);
            parse_brace_group(p, &content);
            /* Check for \limits or \nolimits */
            int has_limits = 0;
            if (peek(p)->type == TOK_COMMAND && strcmp(peek(p)->value, "\\limits") == 0) {
                advance(p); has_limits = 1;
            }
            /* Just return the content for now (simplified) */
            if (content.count == 1) return content.items[0];
            AstNode *g = alloc_node(&p->alloc, ND_GROUP);
            g->u.group.children = content;
            return g;
        }

        /* Named functions: \sin, \cos, etc. */
        if (cmd[0] == '\\' && is_named_function(cmd + 1)) {
            advance(p);
            AstNode *n = alloc_node(&p->alloc, ND_FUNCTION);
            strncpy(n->u.func.name, cmd + 1, 63);
            n->u.func.is_operator = 0;
            return n;
        }

        /* Named operators: \curl, \div, \grad, etc. */
        if (cmd[0] == '\\' && is_named_operator(cmd + 1)) {
            advance(p);
            AstNode *n = alloc_node(&p->alloc, ND_FUNCTION);
            strncpy(n->u.func.name, cmd + 1, 63);
            n->u.func.is_operator = 1;
            return n;
        }

        /* Known LaTeX symbol commands */
        if (strcmp(cmd, "\\prime") == 0) {
            advance(p);
            AstNode *n = alloc_node(&p->alloc, ND_SYMBOL);
            n->u.sym.code = 0x2032;
            strcpy(n->u.sym.latex, "\\prime");
            return n;
        }
        if (strcmp(cmd, "\\circ") == 0) {
            advance(p);
            AstNode *n = alloc_node(&p->alloc, ND_DEGREE);
            return n;
        }

        /* General symbol lookup */
        {
            int code = lookup_latex_unicode(cmd);
            if (code >= 0) {
                advance(p);
                AstNode *n = alloc_node(&p->alloc, ND_SYMBOL);
                n->u.sym.code = (uint16_t)code;
                strncpy(n->u.sym.latex, cmd, 63);
                return n;
            }
        }

        /* Unknown command — emit as function-font chars */
        advance(p);
        AstNode *n = alloc_node(&p->alloc, ND_CHAR);
        n->u.chr.ch = 0;
        n->u.chr.typeface = 5; /* 'command' typeface */
        strncpy((char*)&n->u.chr.code, cmd, 2); /* store cmd pointer hack — use ND_FUNCTION instead */
        /* Better: create function node */
        n = alloc_node(&p->alloc, ND_FUNCTION);
        strncpy(n->u.func.name, cmd + 1, 63);
        n->u.func.is_operator = 0;
        return n;
    }

    /* Character token */
    if (t->type == TOK_CHAR) {
        advance(p);
        char c = t->value[0];
        /* Asterisk: context-dependent —
         * with adjacent space → symbol U+2217 (ASTERISK OPERATOR)
         * without → function char (Hodge star etc.) */
        if (c == '*') {
            int tpos = t->pos;
            int has_space = (tpos > 0 && p->latex && p->latex[tpos - 1] == ' ') ||
                            (p->latex && p->latex[tpos + 1] == ' ');
            if (has_space) {
                AstNode *n = alloc_node(&p->alloc, ND_SYMBOL);
                n->u.sym.code = 0x2217;
                return n;
            }
            AstNode *n = alloc_node(&p->alloc, ND_CHAR);
            n->u.chr.ch = '*'; n->u.chr.typeface = 4; /* plain → func */
            return n;
        }
        AstNode *n = alloc_node(&p->alloc, ND_CHAR);
        n->u.chr.ch = c;
        if (isdigit((unsigned char)c)) {
            n->u.chr.typeface = 1; /* number */
        } else if (isalpha((unsigned char)c)) {
            n->u.chr.typeface = 0; /* variable (italic) */
        } else if (c == '+' || c == '-' || c == '=' || c == '<' || c == '>' ||
                   c == '/' || c == '!' || c == ',' || c == ';' || c == ':' ||
                   c == '(' || c == ')' || c == '[' || c == ']' || c == '|' ||
                   c == '.' || c == '\'') {
            n->u.chr.typeface = 4; /* plain */
        } else {
            n->u.chr.typeface = 0; /* default italic */
        }
        return n;
    }

    /* Skip unexpected tokens */
    advance(p);
    return NULL;
}

/* Parse script argument: {expr} or single token */
static void parse_script_arg(Parser *p, NodeList *out) {
    if (peek(p)->type == TOK_LBRACE) {
        advance(p);
        parse_brace_group(p, out);
    } else {
        AstNode *n = parse_group(p);
        if (n) nl_push(out, n);
    }
}

/* Try to attach sub/superscripts to the last node in the list */
static void try_attach_scripts(Parser *p, NodeList *nodes) {
    if (peek(p)->type != TOK_SUBSCRIPT && peek(p)->type != TOK_SUPERSCRIPT)
        return;

    AstNode *base_node = (nodes->count > 0) ? nodes->items[nodes->count - 1] : NULL;

    /* Don't attach scripts to integrals/bigops — body consumer handles them */
    if (base_node && (base_node->type == ND_INTEGRAL || base_node->type == ND_BIGOP))
        return;

    AstNode *n = alloc_node(&p->alloc, ND_SCRIPT);
    nl_init(&n->u.script.base);
    nl_init(&n->u.script.sub);
    nl_init(&n->u.script.sup);
    n->u.script.has_sub = 0;
    n->u.script.has_sup = 0;

    /* Don't attach to binary operators — use empty base */
    if (base_node && !is_binary_op_node(base_node)) {
        nl_push(&n->u.script.base, base_node);
        nodes->count--; /* remove base from list */
    }

    /* Parse sub then sup, or sup then sub */
    for (int pass = 0; pass < 2; pass++) {
        if (peek(p)->type == TOK_SUBSCRIPT && !n->u.script.has_sub) {
            advance(p);
            n->u.script.has_sub = 1;
            parse_script_arg(p, &n->u.script.sub);
        } else if (peek(p)->type == TOK_SUPERSCRIPT && !n->u.script.has_sup) {
            advance(p);
            n->u.script.has_sup = 1;
            parse_script_arg(p, &n->u.script.sup);
        }
    }

    /* Check for prime: ^{\prime} or ^{\prime\prime} (with or without subscript) */
    if (n->u.script.has_sup) {
        int is_prime = 1, count = 0;
        for (int i = 0; i < n->u.script.sup.count; i++) {
            AstNode *s = n->u.script.sup.items[i];
            if (s->type == ND_SYMBOL && strcmp(s->u.sym.latex, "\\prime") == 0)
                count++;
            else { is_prime = 0; break; }
        }
        if (is_prime && count > 0) {
            /* Convert to PrimeNode */
            AstNode *pn = alloc_node(&p->alloc, ND_PRIME);
            pn->u.prime.count = count;
            n->u.script.sup.count = 0;
            nl_push(&n->u.script.sup, pn);
        }
    }

    /* ^{a/b} → ^{\frac{a}{b}} conversion for EQNEDT32 rendering.
     * Slash in superscript/subscript renders poorly as literal text;
     * a tmFRACT template renders as a proper stacked fraction.
     * Only convert when ALL nodes are simple alphanumeric chars (no operators
     * like +, -, =, etc.) to avoid misinterpreting e.g. {n+1/2}. */
    for (int which = 0; which < 2; which++) {
        NodeList *sl = (which == 0) ? &n->u.script.sup : &n->u.script.sub;
        int has = (which == 0) ? n->u.script.has_sup : n->u.script.has_sub;
        if (!has || sl->count < 3) continue;
        /* Check all nodes are simple ND_CHAR (alphanumeric or '/') */
        int all_simple = 1, slash_pos = -1, slash_count = 0;
        for (int i = 0; i < sl->count; i++) {
            AstNode *item = sl->items[i];
            if (item->type != ND_CHAR) { all_simple = 0; break; }
            char c = item->u.chr.ch;
            if (c == '/') { slash_pos = i; slash_count++; }
            else if (!((c >= '0' && c <= '9') || (c >= 'a' && c <= 'z') ||
                       (c >= 'A' && c <= 'Z'))) {
                all_simple = 0; break;
            }
        }
        if (!all_simple || slash_count != 1 ||
            slash_pos == 0 || slash_pos == sl->count - 1)
            continue;
        /* Build ND_FRAC from the split */
        AstNode *frac = alloc_node(&p->alloc, ND_FRAC);
        nl_init(&frac->u.frac.numer);
        nl_init(&frac->u.frac.denom);
        frac->u.frac.display = 0;
        for (int i = 0; i < slash_pos; i++)
            nl_push(&frac->u.frac.numer, sl->items[i]);
        for (int i = slash_pos + 1; i < sl->count; i++)
            nl_push(&frac->u.frac.denom, sl->items[i]);
        sl->count = 0;
        nl_push(sl, frac);
    }

    /* ^{\circ} → DegreeNode */
    if (n->u.script.has_sup && n->u.script.sup.count == 1 &&
        n->u.script.sup.items[0]->type == ND_DEGREE) {
        /* keep as script with degree */
    }

    nl_push(nodes, n);
}

/* Check if node is an operator that shouldn't stop +/- body termination */
static int is_operator_node(AstNode *node) {
    if (!node) return 0;
    if (node->type == ND_SYMBOL) {
        const char *l = node->u.sym.latex;
        return strcmp(l, "\\cdot") == 0 || strcmp(l, "\\times") == 0 ||
               strcmp(l, "\\pm") == 0 || strcmp(l, "\\mp") == 0 ||
               strcmp(l, "\\wedge") == 0 || strcmp(l, "\\vee") == 0 ||
               strcmp(l, "\\cap") == 0 || strcmp(l, "\\cup") == 0 ||
               strcmp(l, "\\oplus") == 0 || strcmp(l, "\\otimes") == 0 ||
               strcmp(l, "\\odot") == 0;
    }
    return 0;
}

/* Check if node is a binary/relational op that shouldn't be a script base */
static int is_binary_op_node(AstNode *node) {
    if (!node) return 0;
    if (node->type == ND_CHAR && node->u.chr.typeface == 4) {
        char c = node->u.chr.ch;
        return c == '+' || c == '-' || c == '=' || c == ',' || c == ';' || c == ':';
    }
    if (node->type == ND_SYMBOL) {
        const char *l = node->u.sym.latex;
        return strcmp(l, "+") == 0 || strcmp(l, "-") == 0 || strcmp(l, "=") == 0 ||
               strcmp(l, "<") == 0 || strcmp(l, ">") == 0 ||
               strcmp(l, "\\cdot") == 0 || strcmp(l, "\\times") == 0 ||
               strcmp(l, "\\pm") == 0 || strcmp(l, "\\mp") == 0 ||
               strcmp(l, "\\div") == 0 || strcmp(l, "\\wedge") == 0 ||
               strcmp(l, "\\vee") == 0 || strcmp(l, "\\cap") == 0 ||
               strcmp(l, "\\cup") == 0 || strcmp(l, "\\leq") == 0 ||
               strcmp(l, "\\geq") == 0 || strcmp(l, "\\neq") == 0 ||
               strcmp(l, "\\approx") == 0 || strcmp(l, "\\equiv") == 0 ||
               strcmp(l, "\\sim") == 0 || strcmp(l, "\\in") == 0 ||
               strcmp(l, "\\subset") == 0 || strcmp(l, "\\supset") == 0;
    }
    return 0;
}

/* Parse standalone _/^ as ScriptNode with empty base */
static void parse_standalone_script(Parser *p, NodeList *out) {
    AstNode *n = alloc_node(&p->alloc, ND_SCRIPT);
    nl_init(&n->u.script.base);
    nl_init(&n->u.script.sub);
    nl_init(&n->u.script.sup);
    n->u.script.has_sub = 0;
    n->u.script.has_sup = 0;

    for (int pass = 0; pass < 2; pass++) {
        if (peek(p)->type == TOK_SUBSCRIPT && !n->u.script.has_sub) {
            advance(p); n->u.script.has_sub = 1;
            parse_script_arg(p, &n->u.script.sub);
        } else if (peek(p)->type == TOK_SUPERSCRIPT && !n->u.script.has_sup) {
            advance(p); n->u.script.has_sup = 1;
            parse_script_arg(p, &n->u.script.sup);
        }
    }
    nl_push(out, n);
}

/* Consume body content after an integral/bigop */
static void consume_bigop_body(Parser *p, NodeList *body) {
    int paren_depth = 0;
    int fence_depth = 0;  /* \left..\right nesting */
    while (peek(p)->type != TOK_EOF) {
        Token *t = peek(p);

        /* Stop at expression boundaries */
        if (t->type == TOK_RBRACE || t->type == TOK_NEWLINE || t->type == TOK_AMPERSAND)
            break;
        if (t->type == TOK_COMMAND && strcmp(t->value, "\\end") == 0)
            break;
        /* Stop before \right only at fence_depth 0 (unmatched \right = outer fence) */
        if (fence_depth == 0 && t->type == TOK_COMMAND && strcmp(t->value, "\\right") == 0)
            break;

        /* Track \left/\right depth */
        if (t->type == TOK_COMMAND && strcmp(t->value, "\\left") == 0)
            fence_depth++;
        else if (t->type == TOK_COMMAND && strcmp(t->value, "\\right") == 0)
            fence_depth = fence_depth > 0 ? fence_depth - 1 : 0;

        /* Track paren depth */
        if (t->type == TOK_CHAR && (t->value[0] == '(' || t->value[0] == '['))
            paren_depth++;
        else if (t->type == TOK_CHAR && (t->value[0] == ')' || t->value[0] == ']'))
            paren_depth = paren_depth > 0 ? paren_depth - 1 : 0;

        /* Only stop at depth 0 */
        if (paren_depth == 0) {
            /* Stop before another integral/bigop.
             * Nested BigOps (e.g. \sum_i \sum_j) crash EQNEDT32 via DS Equation
             * paste, so consecutive BigOps must stay as siblings. */
            if (t->type == TOK_COMMAND && is_bigop_cmd(t->value))
                break;
            /* Stop before '=' */
            if (t->type == TOK_CHAR && t->value[0] == '=')
                break;
            /* Stop before '+'/'-' when preceded by non-operator */
            if (t->type == TOK_CHAR && (t->value[0] == '+' || t->value[0] == '-')) {
                if (body->count > 0 && !is_operator_node(body->items[body->count - 1]))
                    break;
            }
        }

        /* Handle standalone _/^ in body */
        if (t->type == TOK_SUBSCRIPT || t->type == TOK_SUPERSCRIPT) {
            parse_standalone_script(p, body);
            continue;
        }

        AstNode *node = parse_group(p);
        if (node) {
            nl_push(body, node);
            try_attach_scripts(p, body);
        }
    }
}

/* Main expression parser */
static void parse_expr(Parser *p, NodeList *out, int stop_at) {
    while (peek(p)->type != TOK_EOF) {
        Token *t = peek(p);

        /* Stop conditions: always stop at }, \\, & (structure tokens) */
        if (t->type == stop_at) break;
        if (t->type == TOK_RBRACE) break;
        if (t->type == TOK_NEWLINE) break;
        if (t->type == TOK_AMPERSAND) break;

        /* Stop at \end{...} for environments */
        if (t->type == TOK_COMMAND && strcmp(t->value, "\\end") == 0) {
            Token look = p->pos + 1 < p->tl->count ? p->tl->tokens[p->pos + 1] : p->tl->tokens[p->pos];
            if (look.type == TOK_LBRACE) break;
        }

        /* Handle standalone _/^ */
        if (t->type == TOK_SUBSCRIPT || t->type == TOK_SUPERSCRIPT) {
            parse_standalone_script(p, out);
            continue;
        }

        AstNode *node = parse_group(p);
        if (node) {
            nl_push(out, node);
            try_attach_scripts(p, out);
            /* If last node is integral/bigop, consume body */
            AstNode *last = out->items[out->count - 1];
            if (last->type == ND_INTEGRAL || last->type == ND_BIGOP) {
                NodeList *body_list = (last->type == ND_INTEGRAL) ?
                    &last->u.integ.body : &last->u.bigop.body;
                consume_bigop_body(p, body_list);
            }
        }
    }

    /* Display echo detection removed:
     * EQNEDT32 uses TM_DINT/TM_TINT for multi-integrals (not echo TM_SINT records),
     * so this logic was never needed for tex->mtef.  It incorrectly suppressed a
     * standalone \int that follows \oint (e.g. "...dlC\int" pattern). */
}

/* ============================================================
 * MTEF Writer — AST → MTEF binary
 * ============================================================ */

/* Forward declaration.
 * deferred: optional ByteBuf for "promoted" display data (root level).
 * trailing: optional ByteBuf for content at the enclosing fence level
 *           (e.g. tmSCRIPT that should be inside integral body but outside fence). */
static void emit_nodes(ByteBuf *bb, NodeList *nodes, int in_frac, ByteBuf *deferred);
static void emit_node(ByteBuf *bb, AstNode *node, int in_frac, ByteBuf *deferred);
/* Extended version with trailing buffer for BigOp script splitting */
static void emit_node_ext(ByteBuf *bb, AstNode *node, int in_frac,
                          ByteBuf *deferred, ByteBuf *trailing);

/* Add embellishment to last CHAR record in bb */
static void add_embell_to_last(ByteBuf *bb, int embell_type) {
    if (bb->len >= 4) {
        uint8_t tag = bb->data[bb->len - 4];
        if ((tag & 0x0F) == 0x02) {
            bb->data[bb->len - 4] = tag | 0x20; /* set OPT_CHAR_EMBELL */
            bb_byte(bb, 0x06);
            bb_byte(bb, (uint8_t)embell_type);
            /* No 0x00 terminator: MTEF reader stops at non-0x06 byte */
            return;
        }
    }
}

/* Chain embellishment onto existing embell chain or add new one */
static void chain_embell(ByteBuf *bb, int embell_type) {
    /* Check if buffer already has embellishments: last 2 bytes = 0x06 type */
    if (bb->len >= 2 && bb->data[bb->len - 2] == 0x06) {
        bb_byte(bb, 0x06);
        bb_byte(bb, (uint8_t)embell_type);
        return;
    }
    /* No existing embellishment — add from scratch */
    add_embell_to_last(bb, embell_type);
}

static void emit_node(ByteBuf *bb, AstNode *node, int in_frac, ByteBuf *deferred) {
    emit_node_ext(bb, node, in_frac, deferred, NULL);
}
static void emit_node_ext(ByteBuf *bb, AstNode *node, int in_frac,
                           ByteBuf *deferred, ByteBuf *trailing) {
    if (!node) return;

    switch (node->type) {
    case ND_CHAR: {
        char c = node->u.chr.ch;
        switch (node->u.chr.typeface) {
        case 0: emit_italic(bb, c); break;     /* variable */
        case 1: emit_number(bb, c); break;     /* number */
        case 2: emit_text_char(bb, c); break;  /* text */
        case 3: emit_func_char(bb, c); break;  /* function */
        case 4: /* plain */
            if (c == '(' || c == ')' || c == '[' || c == ']' || c == '|') {
                emit_func_char(bb, c);
            } else if (c == '+' || c == '-' || c == '=' || c == '<' || c == '>' ||
                       c == '/' || c == '.' || c == ',' || c == ';' || c == ':' ||
                       c == '!') {
                /* EQNEDT32 uses U+2212 (MINUS SIGN) not U+002D (HYPHEN-MINUS) */
                emit_symbol(bb, (c == '-') ? 0x2212 : (uint16_t)c);
            } else {
                /* Fallback (including ') — use italic to avoid SYMBOL font mismap */
                emit_italic(bb, c);
            }
            break;
        default: emit_italic(bb, c); break;
        }
        break;
    }

    case ND_SYMBOL: {
        uint16_t code = node->u.sym.code;
        if (is_greek_lc(code))
            emit_greek_lc(bb, code);
        else if (is_greek_uc(code))
            emit_greek_uc(bb, code);
        else
            emit_symbol(bb, code);
        break;
    }

    case ND_RM: {
        int code = lookup_greek_rm(node->u.rm.ch);
        if (code >= 0)
            emit_greek_uc(bb, (uint16_t)code);
        else
            emit_text_char(bb, node->u.rm.ch);
        break;
    }

    case ND_FRAC: {
        ByteBuf num, den;
        bb_init(&num); bb_init(&den);
        emit_nodes(&num, &node->u.frac.numer, 1, NULL);
        emit_nodes(&den, &node->u.frac.denom, 1, NULL);
        /* TMPL tmFRACT var=0: NO SIZE_FULL between numer/denom (native pattern) */
        bb_byte(bb, 0x03); bb_byte(bb, TM_FRACT); bb_byte(bb, 0x00);
        EMIT_END(bb); /* slot[0]: empty */
        EMIT_LINE(bb); bb_append(bb, &num); EMIT_END(bb); /* numerator */
        EMIT_LINE(bb); bb_append(bb, &den); EMIT_END(bb); /* denominator */
        EMIT_END(bb); /* slot[1] end */
        bb_free(&num); bb_free(&den);
        break;
    }

    case ND_SQRT: {
        ByteBuf content;
        bb_init(&content);
        /* in_frac=2: fence templates inside tmROOT crash EQNEDT32 DS paste.
         * Convert fences to inline bracket chars instead. */
        emit_nodes(&content, &node->u.sq.content, 2, NULL);
        if (node->u.sq.has_index) {
            ByteBuf idx;
            bb_init(&idx);
            emit_nodes(&idx, &node->u.sq.index, 2, NULL);
            /* tmROOT var=1, 2 slots */
            bb_byte(bb, 0x03); bb_byte(bb, TM_ROOT); bb_byte(bb, 0x01);
            EMIT_END(bb); /* slot[0]: empty */
            EMIT_LINE(bb); bb_append(bb, &content); EMIT_END(bb);
            EMIT_SIZE_SUB(bb); EMIT_LINE(bb); bb_append(bb, &idx); EMIT_END(bb);
            EMIT_END(bb); /* slot[1] end */
            bb_free(&idx);
        } else {
            /* tmROOT var=0: native format is slot[0]=empty, then extra
             * content LINE after it. The binary patch at 0x4423C7 prevents
             * the cold-start NULL crash that this format triggers.
             * SIZE_FULL after tmROOT END restores font size (native pattern). */
            /* Strip trailing SIZE_FULL from content (native places it
             * after the tmROOT scope, not inside content LINE).
             * tmROOT trailing: END + NULL_LINE + END + END + [SIZE_FULL]. */
            int had_trailing_sf = 0;
            if (content.len > 0 && content.data[content.len - 1] == SIZE_FULL_B) {
                content.len--;
                had_trailing_sf = 1;
            }
            bb_byte(bb, 0x03); bb_byte(bb, TM_ROOT); bb_byte(bb, 0x00);
            EMIT_END(bb); /* slot[0]: empty */
            EMIT_LINE(bb); bb_append(bb, &content); EMIT_END(bb);
            /* Trailing NULL_LINE only when content ended with SIZE_FULL
             * (i.e. had tmSCRIPT). When content ends with display data
             * (integral etc.), NULL_LINE is already present. */
            if (had_trailing_sf) EMIT_NULL_LINE(bb);
            EMIT_END(bb);
        }
        bb_free(&content);
        break;
    }

    case ND_SCRIPT: {
        ByteBuf base_bb;
        bb_init(&base_bb);
        emit_nodes(&base_bb, &node->u.script.base, in_frac, NULL);

        /* Check for prime embellishment */
        if (node->u.script.has_sup && node->u.script.sup.count == 1 &&
            node->u.script.sup.items[0]->type == ND_PRIME) {
            int pcount = node->u.script.sup.items[0]->u.prime.count;
            int etype = (pcount == 1) ? EM_PRIME : (pcount == 2) ? EM_DPRIME : EM_BPRIME;
            if (!node->u.script.has_sub) {
                /* Prime only — add embellishment to base */
                bb_append(bb, &base_bb);
                add_embell_to_last(bb, etype);
                bb_free(&base_bb);
                break;
            }
            /* Prime + subscript */
            bb_append(bb, &base_bb);
            chain_embell(bb, etype);
            ByteBuf sub_bb;
            bb_init(&sub_bb);
            emit_nodes(&sub_bb, &node->u.script.sub, in_frac, NULL);
            bb_byte(bb, 0x03); bb_byte(bb, TM_SCRIPT); bb_byte(bb, 0x01);
            EMIT_END(bb); EMIT_SIZE_SUB(bb);
            EMIT_LINE(bb); bb_append(bb, &sub_bb); EMIT_END(bb);
            EMIT_NULL_LINE(bb); EMIT_SIZE_FULL(bb); EMIT_END(bb);
            bb_free(&sub_bb);
            bb_free(&base_bb);
            break;
        }

        /* Degree: ^{\circ} */
        if (node->u.script.has_sup && node->u.script.sup.count == 1 &&
            node->u.script.sup.items[0]->type == ND_DEGREE &&
            !node->u.script.has_sub) {
            bb_append(bb, &base_bb);
            /* Superscript with degree symbol */
            bb_byte(bb, 0x03); bb_byte(bb, TM_SCRIPT); bb_byte(bb, 0x00);
            EMIT_END(bb); EMIT_SIZE_SUB(bb);
            EMIT_NULL_LINE(bb); EMIT_LINE(bb);
            emit_symbol(bb, 0x00B0);
            EMIT_END(bb); EMIT_SIZE_FULL(bb); EMIT_END(bb);
            bb_free(&base_bb);
            break;
        }

        /* Normal scripts */
        bb_append(bb, &base_bb);
        bb_free(&base_bb);


        ByteBuf sub_bb, sup_bb;
        bb_init(&sub_bb); bb_init(&sup_bb);
        if (node->u.script.has_sub) emit_nodes(&sub_bb, &node->u.script.sub, in_frac, NULL);
        if (node->u.script.has_sup) emit_nodes(&sup_bb, &node->u.script.sup, in_frac, NULL);

        if (node->u.script.has_sub && node->u.script.has_sup) {
            bb_byte(bb, 0x03); bb_byte(bb, TM_SCRIPT); bb_byte(bb, 0x02);
            EMIT_END(bb); EMIT_SIZE_SUB(bb);
            EMIT_LINE(bb); bb_append(bb, &sub_bb); EMIT_END(bb);
            EMIT_LINE(bb); bb_append(bb, &sup_bb); EMIT_END(bb);
            EMIT_SIZE_FULL(bb); EMIT_END(bb);
        } else if (node->u.script.has_sub) {
            bb_byte(bb, 0x03); bb_byte(bb, TM_SCRIPT); bb_byte(bb, 0x01);
            EMIT_END(bb); EMIT_SIZE_SUB(bb);
            EMIT_LINE(bb); bb_append(bb, &sub_bb); EMIT_END(bb);
            EMIT_NULL_LINE(bb); EMIT_SIZE_FULL(bb); EMIT_END(bb);
        } else if (node->u.script.has_sup) {
            /* Native format: END SIZE_SUB NULL_LINE LINE(sup) END SIZE_FULL END
             * NULL_LINE returns immediately (no children parsed), so SIZE_FULL
             * is at slot[1] level (not nested inside NULL_LINE).
             * SIZE_FULL restores font size for continuation after the template
             * (EQNEDT32 scopes SIZE within template slots). */
            bb_byte(bb, 0x03); bb_byte(bb, TM_SCRIPT); bb_byte(bb, 0x00);
            EMIT_END(bb); EMIT_SIZE_SUB(bb);
            EMIT_NULL_LINE(bb); EMIT_LINE(bb); bb_append(bb, &sup_bb);
            EMIT_END(bb); EMIT_SIZE_FULL(bb); EMIT_END(bb);
        }
        bb_free(&sub_bb); bb_free(&sup_bb);
        break;
    }

    case ND_FENCE: {
        ByteBuf content;
        bb_init(&content);
        emit_nodes(&content, &node->u.fence.content, in_frac, NULL);
        uint16_t dleft, dright;
        fence_display_chars(node->u.fence.selector, &dleft, &dright);

        /* in_frac==2 (inside sqrt): fence templates crash EQNEDT32 on
         * DS Equation paste when inside tmROOT slots.  Emit inline bracket
         * chars instead.  Brackets won't auto-resize but won't crash. */
        if (in_frac == 2) {
            if (node->u.fence.variation == 0 || node->u.fence.variation == 1) {
                if (dleft) emit_symbol(bb, dleft);
            }
            bb_append(bb, &content);
            if (node->u.fence.variation == 0 || node->u.fence.variation == 2) {
                if (dright) emit_symbol(bb, dright);
            }
            bb_free(&content);
            break;
        }

        /* Check if content is a PILE record (environment inside fence) */
        if (content.len >= 2 && content.data[0] == 0x04) {
            /* PILE as direct child of fence template */
            bb_byte(bb, 0x03); bb_byte(bb, (uint8_t)node->u.fence.selector);
            bb_byte(bb, (uint8_t)node->u.fence.variation);
            EMIT_END(bb); /* slot[0]: empty */
            /* Wrap PILE with fcount=1 header */
            bb_byte(bb, 0x04); /* PILE record */
            bb_byte(bb, content.data[1]); /* halign from original */
            bb_byte(bb, 0x01); /* fcount=1 */
            bb_bytes(bb, content.data + 2, content.len - 2); /* rest of PILE */
            EMIT_END(bb); /* PILE scope end */
            /* Display chars */
            if (node->u.fence.variation == 0) {
                emit_display(bb, dleft); emit_display(bb, dright);
            } else if (node->u.fence.variation == 1) {
                emit_display(bb, dleft);
            } else if (node->u.fence.variation == 2) {
                emit_display(bb, dright);
            }
            EMIT_END(bb);
        } else {
            /* Non-PILE content.
             * If deferred is provided, pass it to inner emit_nodes so BigOp
             * display data gets promoted, and also write fence display chars
             * to deferred (with SIZE_FULL prefix, matching native). */
            ByteBuf fence_content;
            ByteBuf fence_trailing;  /* for tmSCRIPT from BigOp (at fence level) */
            bb_init(&fence_content);
            bb_init(&fence_trailing);
            /* Pass trailing so BigOp can place tmSCRIPT outside content LINE.
             * Use NULL for deferred: BigOp display data inside fences must stay
             * at the fence level, not be promoted to a parent BigOp (e.g. when
             * \sum\limits inside \left( \right) inside \int). */
            for (int i = 0; i < node->u.fence.content.count; i++)
                emit_node_ext(&fence_content, node->u.fence.content.items[i],
                              in_frac, NULL, &fence_trailing);

            bb_byte(bb, 0x03); bb_byte(bb, (uint8_t)node->u.fence.selector);
            bb_byte(bb, (uint8_t)node->u.fence.variation);
            EMIT_END(bb); /* slot[0]: empty */
            EMIT_LINE(bb); bb_append(bb, &fence_content); EMIT_END(bb);
            /* Trailing (tmSCRIPT) at fence level, inside integral body */
            bb_append(bb, &fence_trailing);
            bb_free(&fence_content);
            bb_free(&fence_trailing);

            if (node->u.fence.content.count > 0) {
                ByteBuf *dst = deferred ? deferred : bb;
                if (deferred) bb_byte(dst, SIZE_FULL_B);  /* native: SIZE_FULL before fence display */
                if (node->u.fence.variation == 0) {
                    emit_display(dst, dleft); emit_display(dst, dright);
                } else if (node->u.fence.variation == 1) {
                    emit_display(dst, dleft);
                } else if (node->u.fence.variation == 2) {
                    emit_display(dst, dright);
                }
                if (deferred) bb_byte(dst, REC_END_B);
                else EMIT_END(bb);  /* fence END after display chars */
            }
            if (in_frac && node->u.fence.content.count == 0) EMIT_END(bb);
        }
        bb_free(&content);
        break;
    }

    case ND_INTEGRAL: {
        int sel = node->u.integ.selector;
        int int_count = (sel == TM_SINT || sel == TM_SSINT) ? 1 :
                        (sel == TM_DINT || sel == TM_DSINT) ? 2 : 3;

        /* Check if body contains a fence with BigOp \limits inside.
         * Only this specific pattern needs the deferred display data path
         * (native EQNEDT32 promotes display data to root level). */
        int body_has_fence_bigop = 0;
        for (int i = 0; i < node->u.integ.body.count; i++) {
            AstNode *bi = node->u.integ.body.items[i];
            if (bi->type == ND_FENCE) {
                for (int j = 0; j < bi->u.fence.content.count; j++) {
                    AstNode *fc = bi->u.fence.content.items[j];
                    if (fc->type == ND_BIGOP && fc->u.bigop.has_limits &&
                        (fc->u.bigop.selector == TM_SUM || fc->u.bigop.selector == TM_ISUM)) {
                        body_has_fence_bigop = 1; break;
                    }
                }
                if (body_has_fence_bigop) break;
            }
        }

        if (!node->u.integ.has_lower && !node->u.integ.has_upper) {
            /* Check empty body first */
            if (node->u.integ.body.count == 0) {
                if (node->u.integ.is_echo) { break; }
                if (node->u.integ.variation == 3)
                    emit_symbol(bb, 0x222E);
                else if (node->u.integ.variation == 2 && int_count == 2)
                    emit_symbol(bb, 0x222F);
                else if (node->u.integ.variation == 2 && int_count == 3)
                    emit_symbol(bb, 0x2230);
                else
                    emit_symbol(bb, 0x222B);
                break;
            }

            /* No limits with body */
            if (node->u.integ.variation == 0 && in_frac == 1) {
                ByteBuf body_bb; bb_init(&body_bb);
                emit_nodes(&body_bb, &node->u.integ.body, in_frac, NULL);
                bb_byte(bb, 0x03); bb_byte(bb, (uint8_t)sel); bb_byte(bb, 0);
                EMIT_LINE(bb); bb_append(bb, &body_bb); EMIT_END(bb);
                EMIT_END(bb);
                bb_free(&body_bb);
            } else if (body_has_fence_bigop && in_frac == 0) {
                /* EQNEDT32 native pattern: when body has a fence,
                 * emit body with deferred collection, then promoted data
                 * (BigOp display + fence display + trailing content)
                 * follows at the same level as the tmSINT template. */
                ByteBuf body_deferred;
                bb_init(&body_deferred);

                /* Find the fence and trailing content */
                int fence_idx = -1;
                for (int i = 0; i < node->u.integ.body.count; i++)
                    if (node->u.integ.body.items[i]->type == ND_FENCE) { fence_idx = i; break; }

                bb_byte(bb, 0x03); bb_byte(bb, (uint8_t)sel);
                bb_byte(bb, (uint8_t)node->u.integ.variation);
                EMIT_END(bb);  /* slot[0] empty */

                /* LINE (integral body) — only fence structure inside */
                EMIT_LINE(bb);
                emit_node(bb, node->u.integ.body.items[fence_idx], in_frac, &body_deferred);
                EMIT_END(bb);  /* close integral body LINE */

                /* Promoted data: BigOp display + fence display chars */
                bb_append(bb, &body_deferred);
                bb_free(&body_deferred);

                /* Trailing content after fence (e.g. d, Ω, ^3) — at same level */
                for (int i = fence_idx + 1; i < node->u.integ.body.count; i++)
                    emit_node(bb, node->u.integ.body.items[i], in_frac, NULL);

                /* If trailing content ends with SIZE_FULL (from tmSCRIPT emit),
                 * replace it with END (native pattern) */
                if (bb->len > 0 && bb->data[bb->len - 1] == SIZE_FULL_B)
                    bb->data[bb->len - 1] = REC_END_B;

                /* Integral display data */
                EMIT_NULL_LINE(bb); EMIT_NULL_LINE(bb); EMIT_SIZE_SYM(bb);
                for (int i = 0; i < int_count; i++) emit_display_int(bb);
            } else {
                /* Standard: no fence in body */
                ByteBuf body_bb; bb_init(&body_bb);
                emit_nodes(&body_bb, &node->u.integ.body, in_frac, NULL);
                bb_byte(bb, 0x03); bb_byte(bb, (uint8_t)sel);
                bb_byte(bb, (uint8_t)node->u.integ.variation);
                EMIT_END(bb);
                EMIT_LINE(bb); bb_append(bb, &body_bb); EMIT_END(bb);
                /* Native integral display data pattern:
                 * SIZE_SUB NULL_LINE NULL_LINE SIZE_SYM display_char(s)
                 * END END SIZE_SUB NULL_LINE */
                EMIT_SIZE_SUB(bb);
                EMIT_NULL_LINE(bb); EMIT_NULL_LINE(bb); EMIT_SIZE_SYM(bb);
                for (int i = 0; i < int_count; i++) emit_display_int(bb);
                EMIT_SIZE_SUB(bb); EMIT_NULL_LINE(bb);
                /* v2.9 fix (matrix-loss): removed BOTH trailing ENDs.
                 * Trailing ENDs after integral display chars close the
                 * enclosing LINE prematurely when integral is followed by
                 * '=' + fence-with-matrix at root level. The size markers
                 * and NULL_LINEs do not need explicit ENDs; they are
                 * fixed-format records in the display data block. */
                bb_free(&body_bb);
            }
        } else {
            /* With limits — use bigop pattern */
            ByteBuf body_bb, lo_bb, hi_bb;
            bb_init(&body_bb); bb_init(&lo_bb); bb_init(&hi_bb);
            emit_nodes(&body_bb, &node->u.integ.body, in_frac, NULL);
            if (node->u.integ.has_lower) emit_nodes(&lo_bb, &node->u.integ.lower, in_frac, NULL);
            if (node->u.integ.has_upper) emit_nodes(&hi_bb, &node->u.integ.upper, in_frac, NULL);

            bb_byte(bb, 0x03); bb_byte(bb, (uint8_t)node->u.integ.selector);
            bb_byte(bb, (uint8_t)node->u.integ.variation);
            EMIT_END(bb);
            EMIT_LINE(bb); bb_append(bb, &body_bb); EMIT_END(bb);
            EMIT_SIZE_SUB(bb);
            if (node->u.integ.has_lower) {
                EMIT_LINE(bb); bb_append(bb, &lo_bb); EMIT_END(bb);
            } else { EMIT_NULL_LINE(bb); }
            if (node->u.integ.has_upper) {
                EMIT_LINE(bb); bb_append(bb, &hi_bb); EMIT_END(bb);
            } else { EMIT_NULL_LINE(bb); }
            EMIT_SIZE_SYM(bb); emit_display_int(bb);
            EMIT_SIZE_FULL(bb);

            bb_free(&body_bb); bb_free(&lo_bb); bb_free(&hi_bb);
        }
        break;
    }

    case ND_BIGOP: {
        ByteBuf body_bb, lo_bb, hi_bb;
        bb_init(&body_bb); bb_init(&lo_bb); bb_init(&hi_bb);
        emit_nodes(&body_bb, &node->u.bigop.body, in_frac, NULL);
        if (node->u.bigop.has_lower) emit_nodes(&lo_bb, &node->u.bigop.lower, in_frac, NULL);
        if (node->u.bigop.has_upper) emit_nodes(&hi_bb, &node->u.bigop.upper, in_frac, NULL);

        /* BigOp symbol map — Unicode code points, NOT the legacy Symbol-font
         * slots.  2026-08-16: sum/product used 0xE5 / 0xD5 here, and EQNEDT32
         * draws NOTHING for those, so an inline `\sum _{j}x_{j}` rendered as
         * the limits with no operator glyph.  The display variants
         * (emit_display_sum / emit_display_prod) were already U+2211 / U+220F,
         * as is EQNEDT32's own output (tmISUM operator slot = 02 86 11 22).
         * The roundtrip test could not see it: mtef2tex maps both spellings
         * back to \sum. */
        uint16_t sym_code = 0x2211; /* default: summation */
        int sel = node->u.bigop.selector;
        if (sel == TM_SUM || sel == TM_ISUM) sym_code = 0x2211;
        else if (sel == TM_PRODUCT || sel == TM_IPRODUCT) sym_code = 0x220F;
        else if (sel == TM_COPRODUCT || sel == 34) sym_code = 0x2210;
        else if (sel == 35 || sel == 36) sym_code = 0x22C3;
        else if (sel == 37 || sel == 38) sym_code = 0x22C2;

        if (node->u.bigop.has_limits) {
            /* Mode 1: explicit \limits — EQNEDT32 native "display data" pattern.
             * TMPL with var=0, slot[0]=empty. Body content and display data
             * (limit LINEs, NULL_LINE, SIZE_SYM, display char) follow as siblings
             * at the parent level, NOT inside the template. */
            int disp_sel = sel;
            if (sel == TM_ISUM) disp_sel = TM_SUM;
            else if (sel == TM_IPRODUCT) disp_sel = TM_PRODUCT;
            else if (sel == 34) disp_sel = TM_COPRODUCT;
            else if (sel == 36) disp_sel = 35;
            else if (sel == 38) disp_sel = 37;

            /* var=0: slot[0]=empty, body in LINE.
             * If the last body node is ND_SCRIPT, split: base goes in body LINE,
             * tmSCRIPT goes OUTSIDE body LINE (native places it at fence level).
             * Display data is promoted to deferred if available. */
            /* Use disp_sel for BigOps. Note: tmPROD etc crash EQNEDT32
             * via DS Equation paste, but work fine for .eqn file generation. */
            /* var=1 when limits are present (native EQNEDT32 format) */
            bb_byte(bb, 0x03); bb_byte(bb, disp_sel);
            bb_byte(bb, (node->u.bigop.has_lower || node->u.bigop.has_upper) ? 0x01 : 0x00);
            EMIT_END(bb);  /* slot[0] empty — TMPL fully closed */

            /* Split trailing script from body if trailing buffer available.
             * Only split when inside a fence (trailing != NULL) so the
             * tmSCRIPT can be placed at the fence parent level. */
            {
                AstNode *trailing_script_node = NULL;
                if (trailing && node->u.bigop.body.count > 0) {
                    AstNode *last = node->u.bigop.body.items[node->u.bigop.body.count - 1];
                    /* Only split when base has embellishment (hat/bar/etc).
                     * Native EQNEDT32 uses the embellished CHAR (tag 0x22) as
                     * the tmSCRIPT anchor. Without EMBELL, splitting causes
                     * EQNEDT32 to crash because the CHAR lacks the 0x22 flag. */
                    if (last->type == ND_SCRIPT && last->u.script.base.count > 0) {
                        AstNode *base_last = last->u.script.base.items[last->u.script.base.count - 1];
                        if (base_last->type == ND_EMBELL)
                            trailing_script_node = last;
                    }
                }

                if (trailing_script_node) {
                    /* Rebuild body_bb without trailing script, emit base only */
                    bb_free(&body_bb); bb_init(&body_bb);
                    for (int i = 0; i < node->u.bigop.body.count - 1; i++)
                        emit_node(&body_bb, node->u.bigop.body.items[i], in_frac, NULL);
                    emit_nodes(&body_bb, &trailing_script_node->u.script.base, in_frac, NULL);
                }

                /* Body LINE */
                EMIT_LINE(bb); bb_append(bb, &body_bb); EMIT_END(bb);

                if (trailing_script_node && trailing) {
                    /* Write tmSCRIPT to trailing so it ends up OUTSIDE
                     * the fence content LINE but INSIDE integral body LINE */
                    ByteBuf sub_s, sup_s;
                    bb_init(&sub_s); bb_init(&sup_s);
                    if (trailing_script_node->u.script.has_sub)
                        emit_nodes(&sub_s, &trailing_script_node->u.script.sub, in_frac, NULL);
                    if (trailing_script_node->u.script.has_sup)
                        emit_nodes(&sup_s, &trailing_script_node->u.script.sup, in_frac, NULL);
                    if (trailing_script_node->u.script.has_sub && trailing_script_node->u.script.has_sup) {
                        bb_byte(trailing, 0x03); bb_byte(trailing, TM_SCRIPT); bb_byte(trailing, 0x02);
                        bb_byte(trailing, REC_END_B); bb_byte(trailing, SIZE_SUB_B);
                        bb_byte(trailing, REC_LINE_B); bb_append(trailing, &sub_s); bb_byte(trailing, REC_END_B);
                        bb_byte(trailing, REC_LINE_B); bb_append(trailing, &sup_s); bb_byte(trailing, REC_END_B);
                        bb_byte(trailing, REC_END_B);
                    } else if (trailing_script_node->u.script.has_sub) {
                        bb_byte(trailing, 0x03); bb_byte(trailing, TM_SCRIPT); bb_byte(trailing, 0x01);
                        bb_byte(trailing, REC_END_B); bb_byte(trailing, SIZE_SUB_B);
                        bb_byte(trailing, REC_LINE_B); bb_append(trailing, &sub_s); bb_byte(trailing, REC_END_B);
                        bb_byte(trailing, NULL_LINE_B); bb_byte(trailing, REC_END_B);
                    } else if (trailing_script_node->u.script.has_sup) {
                        bb_byte(trailing, 0x03); bb_byte(trailing, TM_SCRIPT); bb_byte(trailing, 0x00);
                        bb_byte(trailing, REC_END_B); bb_byte(trailing, SIZE_SUB_B);
                        bb_byte(trailing, NULL_LINE_B);
                        bb_byte(trailing, REC_LINE_B); bb_append(trailing, &sup_s); bb_byte(trailing, REC_END_B);
                        bb_byte(trailing, REC_END_B);
                    }
                    bb_free(&sub_s); bb_free(&sup_s);
                } else if (trailing_script_node) {
                    /* No deferred: emit tmSCRIPT inline (standalone BigOp) */
                    emit_node(bb, trailing_script_node, in_frac, NULL);
                }
            }

            /* Display data → promoted to deferred if available */
            {
                ByteBuf *dst = deferred ? deferred : bb;
                if (node->u.bigop.has_lower) {
                    bb_byte(dst, REC_LINE_B); bb_append(dst, &lo_bb); bb_byte(dst, REC_END_B);
                } else { bb_byte(dst, NULL_LINE_B); }
                if (node->u.bigop.has_upper) {
                    bb_byte(dst, REC_LINE_B); bb_append(dst, &hi_bb); bb_byte(dst, REC_END_B);
                } else { bb_byte(dst, NULL_LINE_B); }
                bb_byte(dst, SIZE_SYM_B);
                if (sel == TM_SUM || sel == TM_ISUM) emit_display_sum(dst);
                else if (sel == TM_PRODUCT || sel == TM_IPRODUCT) emit_display_prod(dst);
                else emit_display_sum(dst);
                /* Native has END END after display char (closes display data scope) */
                if (deferred) { bb_byte(dst, REC_END_B); bb_byte(dst, REC_END_B); }
            }
        } else {
            /* Mode 2/3: Always use tmSUM/tmPROD template (same as Mode 1).
             * char_symbol + tmSCRIPT causes garbled rendering in EQNEDT32.
             * Native EQNEDT32 always uses the BigOp template with display data. */
            /* EQNEDT32 only supports tmSUM (29) for CopyMtefData+ParseMtef.
             * tmPROD/tmCOPROD/tmBIGCUP/tmBIGCAP crash or produce empty trees.
             * Use tmSUM for ALL BigOps, with the correct display char (Π, ∐, ∪, ∩). */
            /* var=1 for \limits style BigOps (native EQNEDT32 uses var=1
             * when display data limits are present, required for DS Equation paste) */
            int bigop_var = (node->u.bigop.has_limits &&
                             (node->u.bigop.has_lower || node->u.bigop.has_upper)) ? 1 : 0;
            bb_byte(bb, 0x03); bb_byte(bb, TM_SUM); bb_byte(bb, (uint8_t)bigop_var);
            EMIT_END(bb);  /* slot[0] empty */
            EMIT_LINE(bb); bb_append(bb, &body_bb); EMIT_END(bb);
            /* Display data: lower, upper, SIZE_SYM + display char */
            if (node->u.bigop.has_lower) {
                bb_byte(bb, REC_LINE_B); bb_append(bb, &lo_bb); bb_byte(bb, REC_END_B);
            } else { bb_byte(bb, NULL_LINE_B); }
            if (node->u.bigop.has_upper) {
                bb_byte(bb, REC_LINE_B); bb_append(bb, &hi_bb); bb_byte(bb, REC_END_B);
            } else { bb_byte(bb, NULL_LINE_B); }
            bb_byte(bb, SIZE_SYM_B);
            if (sel == TM_SUM || sel == TM_ISUM) emit_display_sum(bb);
            else if (sel == TM_PRODUCT || sel == TM_IPRODUCT) emit_display_prod(bb);
            else if (sel == TM_COPRODUCT || sel == 34) emit_char_rec(bb, 0x02, TFW_SYMBOL, 0x2210);
            else if (sel == 35 || sel == 36) emit_char_rec(bb, 0x02, TFW_SYMBOL, 0x22C3); /* ∪ */
            else if (sel == 37 || sel == 38) emit_char_rec(bb, 0x02, TFW_SYMBOL, 0x22C2); /* ∩ */
            else emit_display_sum(bb);
        }
        /* Dead code removed — old Mode 2 with char_symbol + tmSCRIPT */
        if (0) {
            if (node->u.bigop.has_lower && node->u.bigop.has_upper) {
                bb_byte(bb, 0x03); bb_byte(bb, TM_SCRIPT); bb_byte(bb, 0x02);
                EMIT_END(bb); EMIT_SIZE_SUB(bb);
                EMIT_LINE(bb); bb_append(bb, &lo_bb); EMIT_END(bb);
                EMIT_LINE(bb); bb_append(bb, &hi_bb); EMIT_END(bb);
                EMIT_END(bb);
                EMIT_SIZE_FULL(bb);
            } else if (node->u.bigop.has_lower) {
                bb_byte(bb, 0x03); bb_byte(bb, TM_SCRIPT); bb_byte(bb, 0x01);
                EMIT_END(bb); EMIT_SIZE_SUB(bb);
                EMIT_LINE(bb); bb_append(bb, &lo_bb); EMIT_END(bb);
                EMIT_NULL_LINE(bb); EMIT_END(bb);
                EMIT_SIZE_FULL(bb);
            } else {
                bb_byte(bb, 0x03); bb_byte(bb, TM_SCRIPT); bb_byte(bb, 0x00);
                EMIT_END(bb); EMIT_SIZE_SUB(bb);
                EMIT_NULL_LINE(bb); EMIT_LINE(bb); bb_append(bb, &hi_bb);
                EMIT_END(bb); EMIT_END(bb);
                EMIT_SIZE_FULL(bb);
            }
            bb_append(bb, &body_bb);
        }
        bb_free(&body_bb); bb_free(&lo_bb); bb_free(&hi_bb);
        break;
    }

    case ND_DECORATION: {
        ByteBuf content;
        bb_init(&content);
        emit_nodes(&content, &node->u.deco.content, in_frac, NULL);
        bb_byte(bb, 0x03); bb_byte(bb, (uint8_t)node->u.deco.selector);
        bb_byte(bb, (uint8_t)node->u.deco.variation);
        EMIT_END(bb);
        EMIT_LINE(bb); bb_append(bb, &content); EMIT_END(bb);
        bb_free(&content);
        break;
    }

    case ND_BRACE_DECO: {
        ByteBuf content, label;
        bb_init(&content); bb_init(&label);
        emit_nodes(&content, &node->u.bdeco.content, in_frac, NULL);
        emit_nodes(&label, &node->u.bdeco.label, in_frac, NULL);
        bb_byte(bb, 0x03); bb_byte(bb, (uint8_t)node->u.bdeco.selector); bb_byte(bb, 0x00);
        EMIT_END(bb);
        EMIT_LINE(bb); bb_append(bb, &content); EMIT_END(bb);
        EMIT_SIZE_SUB(bb); EMIT_LINE(bb); bb_append(bb, &label); EMIT_END(bb);
        EMIT_SIZE_FULL(bb);
        if (node->u.bdeco.selector == TM_UHBRACE)
            emit_display(bb, 0xFE37);
        else
            emit_display(bb, 0xFE38);
        EMIT_END(bb);
        bb_free(&content); bb_free(&label);
        break;
    }

    case ND_ENVIRONMENT: {
        if (node->u.env.kind == 2 && node->u.env.ncols > 0) {
            /* Matrix: MATRIX record (0x05) */
            int rows = node->u.env.nlines;
            int cols = node->u.env.ncols;
            bb_byte(bb, 0x05);     /* record type: MATRIX */
            bb_byte(bb, 0x01);     /* valign: center */
            bb_byte(bb, 0x01);     /* hjust: center */
            bb_byte(bb, 0x01);     /* vjust: center */
            bb_byte(bb, (uint8_t)rows);
            bb_byte(bb, (uint8_t)cols);
            bb_byte(bb, 0x00);     /* reserved */
            bb_byte(bb, 0x00);     /* reserved */
            /* Emit cells row-major */
            for (int i = 0; i < rows * cols; i++) {
                NodeList *cell = &node->u.env.lines[i];
                if (cell->count == 0) {
                    EMIT_END(bb);  /* null element */
                } else {
                    ByteBuf cellbuf;
                    bb_init(&cellbuf);
                    emit_nodes(&cellbuf, cell, in_frac, NULL);
                    EMIT_LINE(bb); bb_append(bb, &cellbuf); EMIT_END(bb);
                    bb_free(&cellbuf);
                }
            }
        } else if (node->u.env.kind == 1) {
            /* Aligned: PILE halign=3 */
            bb_byte(bb, 0x04); bb_byte(bb, 0x03);
            for (int i = 0; i < node->u.env.nlines; i++) {
                ByteBuf line;
                bb_init(&line);
                emit_nodes(&line, &node->u.env.lines[i], in_frac, NULL);
                EMIT_LINE(bb); bb_append(bb, &line); EMIT_END(bb);
                bb_free(&line);
            }
            EMIT_END(bb);
        } else if (node->u.env.kind >= 3) {
            /* bmatrix/pmatrix/vmatrix/Vmatrix/Bmatrix:
             * Emit as fence TMPL (tmBRACK/tmPAREN/etc.) wrapping a PILE.
             * PILE uses a special halign value (20-24) so mtef2tex can
             * reconstruct the correct environment name on decode.
             * halign: 20=bmatrix, 21=pmatrix, 22=vmatrix, 23=Vmatrix, 24=Bmatrix */
            static const struct { uint8_t sel; uint8_t halign; uint16_t dleft, dright; }
                mat_fence[] = {
                    {0,0,0,0},{0,0,0,0},{0,0,0,0},        /* 0,1,2 unused */
                    {TM_BRACK, 20, 0x005B, 0x005D},       /* 3=bmatrix */
                    {TM_PAREN, 21, 0x0028, 0x0029},       /* 4=pmatrix */
                    {TM_BAR,   22, 0x007C, 0x007C},       /* 5=vmatrix */
                    {TM_DBAR,  23, 0xEC09, 0xEC0A},       /* 6=Vmatrix */
                    {TM_BRACE, 24, 0x007B, 0x007D},       /* 7=Bmatrix */
                };
            int k = node->u.env.kind;
            uint8_t sel    = mat_fence[k].sel;
            uint8_t halign = mat_fence[k].halign;
            uint16_t dleft = mat_fence[k].dleft, dright = mat_fence[k].dright;
            /* Fence TMPL header + empty slot[0] */
            bb_byte(bb, 0x03); bb_byte(bb, sel); bb_byte(bb, 0x00);
            EMIT_END(bb);
            /* PILE directly (parse_pile reads rows as LINE records) */
            bb_byte(bb, 0x04); bb_byte(bb, halign);
            for (int i = 0; i < node->u.env.nlines; i++) {
                ByteBuf line; bb_init(&line);
                emit_nodes(&line, &node->u.env.lines[i], in_frac, NULL);
                EMIT_LINE(bb); bb_append(bb, &line); EMIT_END(bb);
                bb_free(&line);
            }
            EMIT_END(bb); /* terminate PILE */
            /* NOTE: display bracket chars intentionally NOT emitted here.
             * parse_pile uses a lookahead heuristic to detect its end; emitting
             * display chars after the PILE causes parse_pile to misidentify them
             * as content of the last row.  mtef2tex identifies the environment
             * name from the TMPL selector + PILE halign (20-24), so display chars
             * are not needed for roundtrip. */
            (void)dleft; (void)dright;
        } else {
            /* Gathered: PILE halign=0 */
            bb_byte(bb, 0x04); bb_byte(bb, 0x00);
            for (int i = 0; i < node->u.env.nlines; i++) {
                ByteBuf line;
                bb_init(&line);
                emit_nodes(&line, &node->u.env.lines[i], in_frac, NULL);
                EMIT_LINE(bb); bb_append(bb, &line); EMIT_END(bb);
                bb_free(&line);
            }
            EMIT_END(bb);
        }
        break;
    }

    case ND_FUNCTION: {
        for (int i = 0; node->u.func.name[i]; i++)
            emit_func_char(bb, node->u.func.name[i]);
        break;
    }

    case ND_TEXT: {
        const unsigned char *s = (const unsigned char *)node->u.text.text;
        int i = 0;
        while (s[i]) {
            uint32_t cp;
            if (s[i] < 0x80) {
                cp = s[i]; i++;
            } else if ((s[i] & 0xE0) == 0xC0 && s[i+1]) {
                cp = ((s[i] & 0x1F) << 6) | (s[i+1] & 0x3F); i += 2;
            } else if ((s[i] & 0xF0) == 0xE0 && s[i+1] && s[i+2]) {
                cp = ((s[i] & 0x0F) << 12) | ((s[i+1] & 0x3F) << 6) | (s[i+2] & 0x3F); i += 3;
            } else if ((s[i] & 0xF8) == 0xF0 && s[i+1] && s[i+2] && s[i+3]) {
                cp = ((s[i] & 0x07) << 18) | ((s[i+1] & 0x3F) << 12) | ((s[i+2] & 0x3F) << 6) | (s[i+3] & 0x3F); i += 4;
            } else {
                cp = s[i]; i++; /* fallback: emit raw byte */
            }
            emit_char_rec(bb, 0x02, TFW_TEXT, (uint16_t)(cp & 0xFFFF));
        }
        break;
    }

    case ND_MATHBF: {
        for (int i = 0; i < node->u.mathbf.content.count; i++) {
            AstNode *child = node->u.mathbf.content.items[i];
            if (child->type == ND_CHAR && child->u.chr.typeface == 0) {
                emit_vector_char(bb, child->u.chr.ch);
            } else {
                emit_node(bb, child, in_frac, NULL);
            }
        }
        break;
    }

    case ND_EMBELL: {
        ByteBuf content;
        bb_init(&content);
        emit_nodes(&content, &node->u.embell.content, in_frac, NULL);
        bb_append(bb, &content);
        /* Add embellishment to last CHAR */
        if (bb->len >= 4) {
            uint8_t tag = bb->data[bb->len - 4];
            if ((tag & 0x0F) == 0x02) {
                bb->data[bb->len - 4] = tag | 0x20;
                bb_byte(bb, 0x06);
                bb_byte(bb, (uint8_t)node->u.embell.embell_type);
            }
        }
        bb_free(&content);
        break;
    }

    case ND_GROUP:
        emit_nodes(bb, &node->u.group.children, in_frac, NULL);
        break;

    case ND_PRIME: {
        int etype = (node->u.prime.count == 1) ? EM_PRIME :
                    (node->u.prime.count == 2) ? EM_DPRIME : EM_BPRIME;
        bb_byte(bb, 0x06); bb_byte(bb, (uint8_t)etype);
        break;
    }

    case ND_DEGREE:
        emit_symbol(bb, 0x00B0);
        break;

    case ND_OVERSET: {
        /* Check if over is \frown → tmOARC */
        if (node->u.overset.over.count == 1 &&
            node->u.overset.over.items[0]->type == ND_SYMBOL &&
            strcmp(node->u.overset.over.items[0]->u.sym.latex, "\\frown") == 0) {
            ByteBuf content;
            bb_init(&content);
            emit_nodes(&content, &node->u.overset.base, in_frac, NULL);
            bb_byte(bb, 0x03); bb_byte(bb, 48); bb_byte(bb, 0x00); /* tmOARC */
            EMIT_END(bb);
            EMIT_LINE(bb); bb_append(bb, &content); EMIT_END(bb);
            bb_free(&content);
        } else {
            /* Fallback: just output base */
            emit_nodes(bb, &node->u.overset.base, in_frac, NULL);
        }
        break;
    }

    case ND_DIRAC: {
        /* tmDIRAC (selector 45) */
        ByteBuf bra_bb, ket_bb;
        bb_init(&bra_bb); bb_init(&ket_bb);
        emit_nodes(&bra_bb, &node->u.dirac.bra, in_frac, NULL);
        emit_nodes(&ket_bb, &node->u.dirac.ket, in_frac, NULL);
        if (node->u.dirac.variation == 0) {
            /* Full: ⟨bra|ket⟩ */
            bb_byte(bb, 0x03); bb_byte(bb, 45); bb_byte(bb, 0x00);
            EMIT_END(bb);
            EMIT_LINE(bb); bb_append(bb, &bra_bb); EMIT_END(bb);
            EMIT_LINE(bb); bb_append(bb, &ket_bb); EMIT_END(bb);
            emit_display(bb, 0x2329); emit_display(bb, 0x007C); emit_display(bb, 0x232A);
            EMIT_END(bb);
        } else if (node->u.dirac.variation == 1) {
            /* Bra: ⟨bra| */
            bb_byte(bb, 0x03); bb_byte(bb, 45); bb_byte(bb, 0x01);
            EMIT_END(bb);
            EMIT_LINE(bb); bb_append(bb, &bra_bb); EMIT_END(bb);
            emit_display(bb, 0x2329); emit_display(bb, 0x007C);
            EMIT_END(bb);
        } else {
            /* Ket: |ket⟩ */
            bb_byte(bb, 0x03); bb_byte(bb, 45); bb_byte(bb, 0x02);
            EMIT_END(bb);
            EMIT_LINE(bb); bb_append(bb, &ket_bb); EMIT_END(bb);
            emit_display(bb, 0x007C); emit_display(bb, 0x232A);
            EMIT_END(bb);
        }
        bb_free(&bra_bb); bb_free(&ket_bb);
        break;
    }

    default:
        break;
    }
}

static void emit_nodes(ByteBuf *bb, NodeList *nodes, int in_frac, ByteBuf *deferred) {
    for (int i = 0; i < nodes->count; i++)
        emit_node(bb, nodes->items[i], in_frac, deferred);
}

/* ============================================================
 * Public API
 * ============================================================ */

uint8_t *tex_to_mtef(const char *latex, int *outLen) {
    if (!latex || !*latex) return NULL;

    /* Tokenize */
    TokenList tl;
    tokenize(latex, &tl);

    /* Parse */
    Parser parser;
    parser.tl = &tl;
    parser.pos = 0;
    parser.latex = latex;
    alloc_init(&parser.alloc);

    NodeList ast;
    nl_init(&ast);
    parse_expr(&parser, &ast, TOK_EOF);

    /* Build MTEF */
    ByteBuf result;
    bb_init(&result);

    /* MTEF header: 03 01 01 03 0a */
    bb_byte(&result, MTEF_HDR_VER);
    bb_byte(&result, MTEF_HDR_PLAT);
    bb_byte(&result, MTEF_HDR_PROD);
    bb_byte(&result, MTEF_HDR_PRODMAJ);
    bb_byte(&result, MTEF_HDR_PRODMIN);
    EMIT_SIZE_FULL(&result);
    EMIT_LINE(&result);

    emit_nodes(&result, &ast, 0, NULL);

    EMIT_END(&result); /* root LINE end */
    EMIT_END(&result); /* top-level end */

    /* Cleanup */
    alloc_free(&parser.alloc);
    tl_free(&tl);

    return bb_detach(&result, outLen);
}

int looks_like_latex(const char *text) {
    if (!text) return 0;
    /* Skip leading whitespace */
    while (*text == ' ' || *text == '\t' || *text == '\n' || *text == '\r')
        text++;
    /* Dollar-sign delimited: $...$ or $$...$$ */
    if (text[0] == '$' && strlen(text) >= 3) return 1;
    /* General heuristic: any \cmd where cmd is 2+ alpha chars is likely LaTeX.
     * This catches \sigma, \omega, \operatorname, etc. without needing
     * an exhaustive command list. */
    for (const char *p = text; *p; p++) {
        if (*p == '\\') {
            int alen = 0;
            for (const char *q = p + 1; isalpha((unsigned char)*q); q++)
                alen++;
            if (alen >= 2) return 1;
        }
    }
    /* Check for ^{ or _{ patterns (subscript/superscript) */
    if (strstr(text, "^{") || strstr(text, "_{"))
        return 1;
    /* Check for ^ or _ (TeX math mode ignores whitespace between tokens) */
    for (const char *p = text; *p; p++) {
        if (*p == '^' || *p == '_') return 1;
    }
    return 0;
}

/* Strip dollar-sign delimiters: $$...$$ → ..., $...$ → ...
 * Whitespace around delimiters is ignored (matching TeX behavior). */
void strip_dollar_delimiters(char *text) {
    if (!text) return;
    /* Trim leading/trailing whitespace first (TeX ignores it) */
    char *s = text;
    while (*s == ' ' || *s == '\t' || *s == '\n' || *s == '\r') s++;
    if (s != text) memmove(text, s, strlen(s) + 1);
    int len = (int)strlen(text);
    while (len > 0 && (text[len-1] == ' ' || text[len-1] == '\t' ||
                        text[len-1] == '\n' || text[len-1] == '\r'))
        text[--len] = '\0';
    /* Strip $$ ... $$ */
    if (len >= 4 && text[0] == '$' && text[1] == '$' &&
        text[len-1] == '$' && text[len-2] == '$') {
        memmove(text, text + 2, len - 4);
        text[len - 4] = '\0';
    /* Strip $ ... $ */
    } else if (len >= 2 && text[0] == '$' && text[len-1] == '$') {
        memmove(text, text + 1, len - 2);
        text[len - 2] = '\0';
    }
    /* Trim whitespace inside delimiters too */
    s = text;
    while (*s == ' ' || *s == '\t') s++;
    if (s != text) memmove(text, s, strlen(s) + 1);
    len = (int)strlen(text);
    while (len > 0 && (text[len-1] == ' ' || text[len-1] == '\t'))
        text[--len] = '\0';
}

