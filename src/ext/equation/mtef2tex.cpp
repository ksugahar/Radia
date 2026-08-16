/*
 * mtef2tex.cpp -- MTEF v3 binary -> LaTeX
 *
 * Reads the MTEF v3 that Equation Editor 3.0 writes.  This is the
 * corpus-validated converter and the one legacy import goes through; the
 * node-tree path next to it is a partial reimplementation.
 */

#include "mtef2tex.h"
#include "mtef_common.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ============================================================
 * Constants (template selectors, record types) live in mtef_common.h
 * ============================================================ */

/* Record types */
/* ============================================================
 * Arena allocator
 * ============================================================ */

#define ARENA_SIZE 131072  /* 128 KB */

typedef struct {
    uint8_t buf[ARENA_SIZE];
    size_t  used;
} Eq2TexArena;

static void arena_init(Eq2TexArena *a)
{
    a->used = 0;
}

static void *arena_alloc(Eq2TexArena *a, size_t n)
{
    /* 4-byte alignment */
    n = (n + 3) & ~(size_t)3;
    if (a->used + n > ARENA_SIZE)
        return NULL;
    void *p = a->buf + a->used;
    a->used += n;
    memset(p, 0, n);
    return p;
}

/* ============================================================
 * StringBuilder
 * ============================================================ */

typedef struct {
    char  *buf;
    size_t len;
    size_t cap;
} StringBuilder;

static void sb_init(StringBuilder *sb)
{
    sb->cap = 512;
    sb->buf = (char *)malloc(sb->cap);
    sb->len = 0;
    if (sb->buf) sb->buf[0] = '\0';
}

static void sb_grow(StringBuilder *sb, size_t need)
{
    if (!sb->buf) return;
    size_t new_cap = sb->cap;
    while (new_cap < sb->len + need + 1)
        new_cap *= 2;
    if (new_cap != sb->cap) {
        char *p = (char *)realloc(sb->buf, new_cap);
        if (!p) return;
        sb->buf = p;
        sb->cap = new_cap;
    }
}

static void sb_append(StringBuilder *sb, const char *s)
{
    if (!s || !sb->buf) return;
    size_t slen = strlen(s);
    if (slen == 0) return;
    sb_grow(sb, slen);
    memcpy(sb->buf + sb->len, s, slen);
    sb->len += slen;
    sb->buf[sb->len] = '\0';
}

static void sb_append_n(StringBuilder *sb, const char *s, int n)
{
    if (!s || !sb->buf || n <= 0) return;
    sb_grow(sb, (size_t)n);
    memcpy(sb->buf + sb->len, s, (size_t)n);
    sb->len += (size_t)n;
    sb->buf[sb->len] = '\0';
}

static void sb_append_char(StringBuilder *sb, char c)
{
    if (!sb->buf) return;
    sb_grow(sb, 1);
    sb->buf[sb->len++] = c;
    sb->buf[sb->len] = '\0';
}

/* Hands over ownership: sb must not be used after this returns. */
static char *sb_detach(StringBuilder *sb)
{
    char *r = sb->buf;
    sb->buf = NULL;
    sb->len = sb->cap = 0;
    return r;
}

static void sb_free(StringBuilder *sb)
{
    free(sb->buf);
    sb->buf = NULL;
    sb->len = sb->cap = 0;
}

/* Is the content nothing but whitespace? */
static int sb_is_blank(const StringBuilder *sb)
{
    if (!sb->buf) return 1;
    for (size_t i = 0; i < sb->len; i++)
        if (sb->buf[i] != ' ' && sb->buf[i] != '\t' && sb->buf[i] != '\n')
            return 0;
    return 1;
}

/* ============================================================
 * Node types
 * ============================================================ */

typedef enum {
    NODE_LINE, NODE_CHAR, NODE_TMPL, NODE_PILE,
    NODE_MATRIX, NODE_EMBELL, NODE_SIZE, NODE_FONT
} NodeType;

typedef struct MtefNode MtefNode;

/* Child node list */
typedef struct {
    MtefNode **items;
    int        count;
    int        cap;
} NodeList;

struct MtefNode {
    NodeType type;
    union {
        struct { NodeList children; int is_null; } line;
        struct { int typeface; uint16_t char_code; NodeList embells; } ch;
        struct {
            int selector; int variation;
            int orig_variation;  /* variation as it was before Pass 2 */
            NodeList slots;  /* each slot is a LINE node */
            MtefNode *display_lower;  /* set by convert_pile */
            MtefNode *display_upper;
        } tmpl;
        struct { int halign; NodeList lines; } pile;
        struct { int rows; int cols; NodeList elements; } matrix;
        struct { int embell_type; } embell;
        struct { int size_type; } size;  /* SIZETYPE_xxx */
        struct { int font_index; int style; char name[64]; } font;
    } u;
};

/* Node construction helpers */
static MtefNode *new_node(Eq2TexArena *a, NodeType type)
{
    MtefNode *n = (MtefNode *)arena_alloc(a, sizeof(MtefNode));
    if (n) n->type = type;
    return n;
}

/* Append to a NodeList */
static void nl_push(Eq2TexArena *a, NodeList *nl, MtefNode *node)
{
    if (!node) return;
    if (nl->count >= nl->cap) {
        int new_cap = nl->cap < 8 ? 8 : nl->cap * 2;
        MtefNode **new_items = (MtefNode **)arena_alloc(a, sizeof(MtefNode *) * new_cap);
        if (!new_items) return;
        if (nl->items && nl->count > 0)
            memcpy(new_items, nl->items, sizeof(MtefNode *) * nl->count);
        nl->items = new_items;
        nl->cap = new_cap;
    }
    nl->items[nl->count++] = node;
}

/* --- Debug tree dump --- */
/* ============================================================
 * Lookup tables
 * ============================================================ */

typedef struct { uint16_t key; const char *val; } MapEntry;

/* Binary search */
static const char *map_lookup(const MapEntry *map, int n, uint16_t key)
{
    int lo = 0, hi = n - 1;
    while (lo <= hi) {
        int mid = (lo + hi) / 2;
        if (map[mid].key == key) return map[mid].val;
        if (map[mid].key < key) lo = mid + 1;
        else hi = mid - 1;
    }
    return NULL;
}

/* Membership test over a uint16 array */
static int set_contains(const uint16_t *set, int n, uint16_t val)
{
    for (int i = 0; i < n; i++)
        if (set[i] == val) return 1;
    return 0;
}

/* --- SYMBOL_MAP (sorted by key) --- */
static const MapEntry SYMBOL_MAP[] = {
    {0x22, "\\forall "}, {0x24, "\\exists "}, {0x25, "\\%"},
    {0x27, "\\ni "}, {0x2A, "\\ast "}, {0x2B, " + "}, {0x2D, " - "},
    {0x2E, "."}, {0x2F, "/"},
    {0x3C, " < "}, {0x3D, " = "}, {0x3E, " > "},
    {0x40, "\\cong "},
    {0x41, "{\\rm A}"}, {0x42, "{\\rm B}"}, {0x43, "\\Chi "},
    {0x44, "\\Delta "}, {0x45, "{\\rm E}"}, {0x46, "\\Phi "},
    {0x47, "\\Gamma "}, {0x48, "{\\rm H}"}, {0x49, "{\\rm I}"},
    {0x4A, "\\vartheta "}, {0x4B, "{\\rm K}"}, {0x4C, "\\Lambda "},
    {0x4D, "{\\rm M}"}, {0x4E, "{\\rm N}"}, {0x4F, "{\\rm O}"},
    {0x50, "\\Pi "}, {0x51, "\\Theta "}, {0x52, "{\\rm P}"},
    {0x53, "\\Sigma "}, {0x54, "{\\rm T}"}, {0x55, "\\Upsilon "},
    {0x56, "\\varsigma "}, {0x57, "\\Omega "}, {0x58, "\\Xi "},
    {0x59, "\\Psi "}, {0x5A, "{\\rm Z}"},
    {0x5B, "["}, {0x5C, "\\therefore "}, {0x5D, "]"},
    {0x5E, "\\perp "}, {0x5F, "\\_"},
    {0x61, "\\alpha "}, {0x62, "\\beta "}, {0x63, "\\chi "},
    {0x64, "\\delta "}, {0x65, "\\varepsilon "}, {0x66, "\\phi "},
    {0x67, "\\gamma "}, {0x68, "\\eta "}, {0x69, "\\iota "},
    {0x6A, "\\varphi "}, {0x6B, "\\kappa "}, {0x6C, "\\lambda "},
    {0x6D, "\\mu "}, {0x6E, "\\nu "}, {0x6F, "o"},
    {0x70, "\\pi "}, {0x71, "\\theta "}, {0x72, "\\rho "},
    {0x73, "\\sigma "}, {0x74, "\\tau "}, {0x75, "\\upsilon "},
    {0x76, "\\varpi "}, {0x77, "\\omega "}, {0x78, "\\xi "},
    {0x79, "\\psi "}, {0x7A, "\\zeta "},
    {0x7B, "\\{"}, {0x7C, "|"}, {0x7D, "\\}"}, {0x7E, "\\sim "},
    {0xA0, " "}, {0xA1, "\\Upsilon "},
    {0xA2, "\\prime "}, {0xA3, "\\leq "}, {0xA4, "/"},
    {0xA5, "\\infty "}, {0xA7, "\\clubsuit "}, {0xA8, "\\diamondsuit "},
    {0xA9, "\\heartsuit "}, {0xAA, "\\spadesuit "},
    {0xAB, "\\leftrightarrow "}, {0xAC, "\\leftarrow "},
    {0xAD, "\\uparrow "}, {0xAE, "\\rightarrow "}, {0xAF, "\\downarrow "},
    {0xB0, "^{\\circ}"}, {0xB1, "\\pm "},
    {0xB2, "\\prime\\prime "}, {0xB3, "\\geq "}, {0xB4, "\\times "},
    {0xB5, "\\propto "}, {0xB6, "\\partial "}, {0xB7, "\\bullet "},
    {0xB8, "\\div "}, {0xB9, "\\neq "}, {0xBA, "\\equiv "},
    {0xBB, "\\approx "}, {0xBC, "\\cdots "},
    {0xC0, "\\aleph "}, {0xC1, "\\Im "}, {0xC2, "\\Re "}, {0xC3, "\\wp "},
    {0xC4, "\\otimes "}, {0xC5, "\\oplus "}, {0xC6, "\\emptyset "},
    {0xC7, "\\cap "}, {0xC8, "\\cup "},
    {0xC9, "\\supset "}, {0xCA, "\\supseteq "},
    {0xCB, "\\not\\subset "}, {0xCC, "\\subset "}, {0xCD, "\\subseteq "},
    {0xCE, "\\in "}, {0xCF, "\\notin "},
    {0xD0, "\\angle "}, {0xD1, "\\nabla "},
    {0xD5, "\\prod "}, {0xD6, "\\surd "}, {0xD7, "\\cdot "},
    {0xD8, "\\neg "}, {0xD9, "\\wedge "}, {0xDA, "\\vee "},
    {0xDB, "\\Leftrightarrow "}, {0xDC, "\\Leftarrow "},
    {0xDD, "\\Uparrow "}, {0xDE, "\\Rightarrow "}, {0xDF, "\\Downarrow "},
    {0xE0, "\\diamond "}, {0xE1, "\\langle "},
    {0xE5, "\\sum "},
    {0xF1, "\\rangle "}, {0xF2, "\\int "},
    {0xF5, "\\lfloor "}, {0xF6, "\\rfloor "},
    {0xF7, "\\lceil "}, {0xF8, "\\rceil "},
};
#define SYMBOL_MAP_N (sizeof(SYMBOL_MAP)/sizeof(SYMBOL_MAP[0]))

/* --- MTEXTRA_MAP (sorted by key) --- */
static const MapEntry MTEXTRA_MAP[] = {
    {0x24, "\\nexists "},
    {0x25, "\\therefore "},
    {0x27, "\\because "},
    {0x2B, "\\oplus "},
    {0x2F, "\\oslash "},
    {0x3B, "\\triangleleft "},
    {0x3C, "\\triangleright "},
    {0x3E, "\\bigtriangledown "},
    {0x43, "\\widehat"},
    {0x44, "\\widetilde"},
    {0x48, "\\overleftarrow"},
    {0x49, "\\overrightarrow"},
    {0x4A, "\\overleftrightarrow"},
    {0x62, "\\Longleftarrow "},
    {0x63, "\\Longrightarrow "},
    {0x64, "\\Longleftrightarrow "},
};
#define MTEXTRA_MAP_N (sizeof(MTEXTRA_MAP)/sizeof(MTEXTRA_MAP[0]))

/* --- UNICODE_MAP (sorted by key) --- */
static const MapEntry UNICODE_MAP[] = {
    {0x00B1, " \\pm "}, {0x00D7, " \\times "}, {0x00F7, " \\div "},
    {0x0391, "{\\rm A}"}, {0x0392, "{\\rm B}"}, {0x0393, "\\Gamma "},
    {0x0394, "\\Delta "}, {0x0395, "{\\rm E}"}, {0x0396, "{\\rm Z}"},
    {0x0397, "{\\rm H}"}, {0x0398, "\\Theta "}, {0x0399, "{\\rm I}"},
    {0x039A, "{\\rm K}"}, {0x039B, "\\Lambda "}, {0x039C, "{\\rm M}"},
    {0x039D, "{\\rm N}"}, {0x039E, "\\Xi "}, {0x039F, "{\\rm O}"},
    {0x03A0, "\\Pi "}, {0x03A1, "{\\rm P}"}, {0x03A3, "\\Sigma "},
    {0x03A4, "{\\rm T}"}, {0x03A5, "\\Upsilon "}, {0x03A6, "\\Phi "},
    {0x03A7, "{\\rm X}"}, {0x03A8, "\\Psi "}, {0x03A9, "\\Omega "},
    {0x03B1, "\\alpha "}, {0x03B2, "\\beta "}, {0x03B3, "\\gamma "},
    {0x03B4, "\\delta "}, {0x03B5, "\\varepsilon "}, {0x03B6, "\\zeta "},
    {0x03B7, "\\eta "}, {0x03B8, "\\theta "}, {0x03B9, "\\iota "},
    {0x03BA, "\\kappa "}, {0x03BB, "\\lambda "}, {0x03BC, "\\mu "},
    {0x03BD, "\\nu "}, {0x03BE, "\\xi "}, {0x03BF, "o"},
    {0x03C0, "\\pi "}, {0x03C1, "\\rho "}, {0x03C2, "\\varsigma "},
    {0x03C3, "\\sigma "}, {0x03C4, "\\tau "}, {0x03C5, "\\upsilon "},
    {0x03C6, "\\phi "}, {0x03C7, "\\chi "}, {0x03C8, "\\psi "},
    {0x03C9, "\\omega "},
    {0x03D1, "\\vartheta "}, {0x03D5, "\\varphi "}, {0x03D6, "\\varpi "},
    {0x03F1, "\\varrho "},
    {0x2013, "--"}, {0x2014, "---"},
    {0x2020, "\\dag "}, {0x2021, "\\ddag "},
    {0x2022, "\\bullet "}, {0x2026, "\\ldots "},
    {0x2032, "\\prime "}, {0x2044, " / "},
    {0x210F, "\\hbar "}, {0x2111, "\\Im "}, {0x2113, "\\ell "},
    {0x2118, "\\wp "}, {0x211C, "\\Re "}, {0x2135, "\\aleph "},
    {0x2190, " \\leftarrow "}, {0x2191, " \\uparrow "},
    {0x2192, " \\to "}, {0x2193, " \\downarrow "},
    {0x2194, " \\leftrightarrow "}, {0x2195, " \\updownarrow "},
    {0x2196, " \\nwarrow "}, {0x2197, " \\nearrow "},
    {0x2198, " \\searrow "}, {0x2199, " \\swarrow "},
    {0x21A6, " \\mapsto "},
    {0x21A9, " \\hookleftarrow "}, {0x21AA, " \\hookrightarrow "},
    {0x21BC, " \\leftharpoonup "}, {0x21BD, " \\leftharpoondown "},
    {0x21C0, " \\rightharpoonup "}, {0x21C1, " \\rightharpoondown "},
    {0x21D0, " \\Leftarrow "}, {0x21D1, " \\Uparrow "},
    {0x21D2, " \\Rightarrow "}, {0x21D3, " \\Downarrow "},
    {0x21D4, " \\Leftrightarrow "},
    {0x2200, "\\forall "}, {0x2202, "\\partial "}, {0x2203, "\\exists "},
    {0x2205, "\\emptyset "}, {0x2206, "\\Delta "}, {0x2207, "\\nabla "},
    {0x2208, " \\in "}, {0x2209, " \\notin "}, {0x220B, " \\ni "},
    {0x220F, "\\prod "}, {0x2210, "\\coprod "}, {0x2211, "\\sum "},
    {0x2212, " - "}, {0x2213, " \\mp "}, {0x2215, " / "},
    {0x2216, " \\setminus "}, {0x2217, " * "}, {0x2218, " \\circ "},
    {0x2219, " \\bullet "}, {0x221A, "\\surd "},
    {0x221D, " \\propto "}, {0x221E, "\\infty "},
    {0x221F, "\\perp "}, {0x2220, "\\angle "},
    {0x2223, "\\mid "}, {0x2225, "\\parallel "},
    {0x2227, " \\wedge "}, {0x2228, " \\vee "},
    {0x2229, " \\cap "}, {0x222A, " \\cup "},
    {0x222B, "\\int "}, {0x222E, "\\oint "},
    {0x2234, "\\therefore "}, {0x2235, "\\because "},
    {0x2243, " \\simeq "}, {0x2245, " \\cong "},
    {0x2248, " \\approx "}, {0x2260, " \\neq "}, {0x2261, " \\equiv "},
    {0x2264, " \\leq "}, {0x2265, " \\geq "},
    {0x226A, " \\ll "}, {0x226B, " \\gg "},
    {0x227A, " \\prec "}, {0x227B, " \\succ "},
    {0x2282, " \\subset "}, {0x2283, " \\supset "},
    {0x2286, " \\subseteq "}, {0x2287, " \\supseteq "},
    {0x2295, " \\oplus "}, {0x2296, " \\ominus "}, {0x2297, " \\otimes "},
    {0x2299, " \\odot "},
    {0x22A5, " \\perp "}, {0x22C5, " \\cdot "}, {0x22EE, "\\vdots "},
    {0x22EF, "\\cdots "}, {0x22F0, "\\iddots "}, {0x22F1, "\\ddots "},
    {0x2329, "\\langle "}, {0x232A, "\\rangle"},
};
#define UNICODE_MAP_N (sizeof(UNICODE_MAP)/sizeof(UNICODE_MAP[0]))

/* --- SYMBOL_HIBYTE_MAP --- */
static const MapEntry SYMBOL_HIBYTE_MAP[] = {
    {0xD1, "\\nabla "}, {0xD6, "\\sqrt"}, {0xE2, "\\infty "},
};
#define SYMBOL_HIBYTE_MAP_N 3

/* --- EMBELL_MAP (indexed by embell_type) --- */
typedef struct { const char *prefix; const char *suffix; } EmbellEntry;
static const EmbellEntry EMBELL_MAP[] = {
    {"", ""},           /* 0 */
    {"", ""},           /* 1 */
    {"\\dot{", "}"},    /* EM_DOT=2 */
    {"\\ddot{", "}"},   /* EM_DDOT=3 */
    {"\\dddot{", "}"},  /* EM_TDOT=4 */
    {"", "^{\\prime}"},          /* EM_PRIME=5 */
    {"", "^{\\prime\\prime}"},   /* EM_DPRIME=6 */
    {"", "^{\\prime\\prime\\prime}"}, /* EM_BPRIME=7 */
    {"\\tilde{", "}"},  /* EM_TILDE=8 */
    {"\\hat{", "}"},    /* EM_HAT=9 */
    {"\\not ", ""},     /* EM_NOT=10 */
    {"\\vec{", "}"},    /* EM_RARROW=11 */
    {"\\overleftarrow{", "}"}, /* EM_LARROW=12 */
    {"\\overleftrightarrow{", "}"}, /* EM_BARROW=13 */
    {"", ""},           /* EM_R1ARROW=14 */
    {"", ""},           /* EM_L1ARROW=15 */
    {"\\bar{", "}"},    /* EM_MBAR=16 */
    {"\\overline{", "}"}, /* EM_OBAR=17 */
    {"", "^{\\prime\\prime\\prime}"}, /* EM_TPRIME=18 */
    {"\\frown{", "}"},  /* EM_FROWN=19 */
    {"\\smile{", "}"},  /* EM_SMILE=20 */
};
#define EMBELL_MAP_N 21

/* --- Greek lowercase (indexed by 'a'-'z') --- */
static const char *GREEK_LOWER[26] = {
    "\\alpha ", "\\beta ", "\\chi ", "\\delta ", "\\varepsilon ",
    "\\phi ", "\\gamma ", "\\eta ", "\\iota ", "\\varphi ",
    "\\kappa ", "\\lambda ", "\\mu ", "\\nu ", "o",
    "\\pi ", "\\theta ", "\\rho ", "\\sigma ", "\\tau ",
    "\\upsilon ", "\\varpi ", "\\omega ", "\\xi ", "\\psi ", "\\zeta ",
};

/* --- Greek uppercase (indexed by 'A'-'Z') --- */
static const char *GREEK_UPPER[26] = {
    "{\\rm A}", "{\\rm B}", "\\Chi ", "\\Delta ", "{\\rm E}",
    "\\Phi ", "\\Gamma ", "{\\rm H}", "{\\rm I}", "\\vartheta ",
    "{\\rm K}", "\\Lambda ", "{\\rm M}", "{\\rm N}", "{\\rm O}",
    "\\Pi ", "\\Theta ", "{\\rm P}", "\\Sigma ", "{\\rm T}",
    "\\Upsilon ", "\\varsigma ", "\\Omega ", "\\Xi ", "\\Psi ", "{\\rm Z}",
};

/* --- BigOp display chars --- */
static const uint16_t BIGOP_DISPLAY_CHARS[] = {
    0xD5, 0xE5, 0xF2,  /* Symbol font: prod, sum, int */
    0x220F, 0x2210, 0x2211,  /* prod, coprod, sum */
    0x222B, 0x222C, 0x222D, 0x222E, 0x222F, 0x2230, /* integrals */
    0x22C2, 0x22C3,  /* intersection, union */
    0xFE37, 0xFE38,  /* overbrace, underbrace decorative */
};
#define BIGOP_DISPLAY_N (sizeof(BIGOP_DISPLAY_CHARS)/sizeof(BIGOP_DISPLAY_CHARS[0]))

/* --- Fence display chars --- */
static const uint16_t FENCE_DISPLAY_CHARS[] = {
    0x28, 0x29, 0x5B, 0x5D, 0x7B, 0x7D, 0x7C,
    0xE9, 0xEB, 0xEC, 0xF9, 0xFB, 0xFC,
    0x2016, 0x2308, 0x2309, 0x230A, 0x230B,
    0x2329, 0x232A, 0x27E8, 0x27E9,
};
#define FENCE_DISPLAY_N (sizeof(FENCE_DISPLAY_CHARS)/sizeof(FENCE_DISPLAY_CHARS[0]))

/* --- Selector set helpers --- */
static int is_fence_selector(int sel) {
    return sel >= tmANGLE && sel <= tmLPRB;
}
static int is_bigop_selector(int sel) {
    return (sel >= tmSINT && sel <= tmTSINT) ||
           (sel >= tmSUM && sel <= tmIINTER);
}
static int is_display_tmpl_selector(int sel) {
    return is_bigop_selector(sel) || sel == tmUHBRACE || sel == tmLHBRACE;
}
static int is_decoration_selector(int sel) {
    return sel == tmUBAR || sel == tmOBAR ||
           sel == tmLARROW || sel == tmRARROW || sel == tmBARROW;
}

/* Global arena pointer for conversion phase (set by mtef_to_latex_c) */
static Eq2TexArena *g_conv_arena = NULL;
/* Depth counter: cursor garbage removal only at depth 0 (top-level LINE) */
static int g_convert_line_depth = 0;
static MtefNode *g_parent_bigop = NULL;  /* Parent BigOp for display data promotion */

/* BigOp display data stack: tmSCRIPT pushes, tmSUM pops (LIFO order) */
#define BIGOP_DISP_STACK_MAX 16
static struct { MtefNode *lo; MtefNode *hi; } g_bigop_disp_stack[BIGOP_DISP_STACK_MAX];
static int g_bigop_disp_stack_n = 0;

/* forward declarations */
typedef struct MtefParser MtefParser;
static MtefNode *parse_record(MtefParser *p, int rec_type, int options);
static void convert_node(MtefNode *node, int prod_ver, StringBuilder *sb);
static void convert_line(MtefNode *node, int prod_ver, Eq2TexArena *arena, StringBuilder *sb);

/* ============================================================
 * MTEF parser
 * ============================================================ */

struct MtefParser {
    const uint8_t *data;
    size_t         len;
    size_t         pos;
    int            prod_ver;
    Eq2TexArena   *arena;
};

static int p_read_byte(MtefParser *p)
{
    if (p->pos >= p->len) return -1;
    return p->data[p->pos++];
}

static int p_read_uint16(MtefParser *p)
{
    if (p->pos + 2 > p->len) return -1;
    uint16_t v = p->data[p->pos] | (p->data[p->pos + 1] << 8);
    p->pos += 2;
    return v;
}

static int p_read_int16(MtefParser *p)
{
    int v = p_read_uint16(p);
    if (v < 0) return 0;
    if (v >= 32768) v -= 65536;
    return v;
}

static void p_read_nudge(MtefParser *p)
{
    int b1 = p_read_byte(p);
    int b2 = p_read_byte(p);
    if (b1 == 128 && b2 == 128) {
        p_read_int16(p);
        p_read_int16(p);
    }
}

static void p_read_null_string(MtefParser *p)
{
    while (p->pos < p->len && p->data[p->pos] != 0)
        p->pos++;
    if (p->pos < p->len) p->pos++; /* skip null terminator */
}

static void p_read_null_string_to(MtefParser *p, char *buf, int buflen)
{
    int i = 0;
    while (p->pos < p->len && p->data[p->pos] != 0) {
        if (i < buflen - 1) buf[i++] = (char)p->data[p->pos];
        p->pos++;
    }
    buf[i] = '\0';
    if (p->pos < p->len) p->pos++;
}

/* --- parse_object_list --- */
static NodeList parse_object_list(MtefParser *p)
{
    NodeList nl = {NULL, 0, 0};
    int safety = 0;
    while (p->pos < p->len && safety < 10000) {
        int tag = p_read_byte(p);
        if (tag < 0) break;
        int rec_type = tag & 0x0F;
        int options = (tag >> 4) & 0x0F;
        if (rec_type == REC_END) break;
        MtefNode *node = parse_record(p, rec_type, options);
        if (node) {
            nl_push(p->arena, &nl, node);
            /* Fence: flatten slot[1] (trailing content absorbed to consume
             * the fence END) back as siblings for convert_line compatibility. */
            if (node->type == NODE_TMPL &&
                is_fence_selector(node->u.tmpl.selector) &&
                node->u.tmpl.slots.count >= 2) {
                MtefNode *s1 = node->u.tmpl.slots.items[1];
                for (int i = 0; i < s1->u.line.children.count; i++)
                    nl_push(p->arena, &nl, s1->u.line.children.items[i]);
                node->u.tmpl.slots.count = 1;
            }
        }
        safety++;
    }
    return nl;
}

/* --- get_template_slot_count --- */
static int get_template_slot_count(int selector, int variation)
{
    /* Fences: 1 slot */
    if (selector >= tmANGLE && selector <= tmLPRB) return 1;

    if (selector == tmROOT) return (variation == 1) ? 2 : 1;
    if (selector == tmFRACT || selector == tmSLFRACT) return 2;
    if (selector == tmSCRIPT || selector == tmLSCRIPT) return 2;
    if (selector == tmUBAR || selector == tmOBAR) return 1;
    if (selector == tmLARROW || selector == tmRARROW || selector == tmBARROW) return 1;

    /* BigOps: two groups with different slot patterns in EQNEDT32.
     *
     * Integrals (tmSINT..tmTSINT): variation bits encode slot structure.
     *   bit 0 = has slot[1] (integrand body), bit 1 = integral type (e.g., ∮).
     *   e.g., tmSINT var=1 → 2 slots: slot[0]=empty, slot[1]=integrand+dx
     *         tmSINT var=3 → 3 slots: ∮ with content in slot[1], slot[2]
     *
     * Sums/Products (tmSUM..tmIINTER): always 1 slot (slot[0]=empty).
     *   Limits are in display data (parsed by convert_line Pass 2), not in slots.
     *   e.g., tmSUM var=1 → has limits in display data, NOT in extra slots. */
    if (is_bigop_selector(selector)) {
        if (selector >= tmSINT && selector <= tmTSINT) {
            /* Integrals: variation-based slot count */
            int has_lower = variation & 0x01;
            int has_upper = variation & 0x02;
            return 1 + (has_lower ? 1 : 0) + (has_upper ? 1 : 0);
        }
        /* Sums/Products: always 1 slot */
        return 1;
    }

    if (selector == tmUHBRACE || selector == tmLHBRACE) return 2;
    if (selector == tmLIM) return (variation == 2) ? 3 : 2;
    if (selector == tmLDIV) return (variation == 0) ? 2 : 1;

    if (selector == tmINTOP || selector == tmSUMOP) {
        int has_lower = variation & 0x01;
        int has_upper = variation & 0x02;
        return (has_lower ? 1 : 0) + (has_upper ? 1 : 0) + 1;
    }

    if (selector == tmDIRAC) return (variation == 0) ? 2 : 1;
    if (selector == tmUARROW || selector == tmOARROW || selector == tmOARC) return 1;

    return 1; /* default */
}

/* --- parse_line --- */
static MtefNode *parse_line(MtefParser *p, int options)
{
    int has_nudge = options & OPT_NUDGE;
    int is_null = options & OPT_LINE_NULL;
    int has_lspace = options & OPT_LINE_LSPACE;

    if (has_nudge) p_read_nudge(p);
    if (has_lspace) p_read_uint16(p);

    MtefNode *n = new_node(p->arena, NODE_LINE);
    if (!n) return NULL;

    if (is_null) {
        n->u.line.is_null = 1;
        return n;
    }

    n->u.line.children = parse_object_list(p);
    return n;
}

/* --- parse_embell --- */
static MtefNode *parse_embell(MtefParser *p, int options)
{
    if (options & OPT_NUDGE) p_read_nudge(p);
    int etype = p_read_byte(p);
    if (etype < 0) return NULL;
    MtefNode *n = new_node(p->arena, NODE_EMBELL);
    if (!n) return NULL;
    n->u.embell.embell_type = etype;
    return n;
}

/* --- parse_embell_list --- */
static NodeList parse_embell_list(MtefParser *p)
{
    NodeList nl = {NULL, 0, 0};
    while (p->pos < p->len) {
        int tag = p_read_byte(p);
        if (tag < 0) break;
        int rec_type = tag & 0x0F;
        int options = (tag >> 4) & 0x0F;
        if (rec_type == REC_END) { p->pos--; break; }
        if (rec_type == REC_EMBELL) {
            MtefNode *node = parse_embell(p, options);
            if (node) nl_push(p->arena, &nl, node);
        } else {
            p->pos--;
            break;
        }
    }
    return nl;
}

/* --- parse_char --- */
static MtefNode *parse_char(MtefParser *p, int options)
{
    if (options & OPT_NUDGE) p_read_nudge(p);
    int tf_byte = p_read_byte(p);
    if (tf_byte < 0) return NULL;
    int typeface = tf_byte - 128;
    int char_code = p_read_uint16(p);
    if (char_code < 0) return NULL;

    MtefNode *n = new_node(p->arena, NODE_CHAR);
    if (!n) return NULL;
    n->u.ch.typeface = typeface;
    n->u.ch.char_code = (uint16_t)char_code;

    if (options & OPT_CHAR_EMBELL)
        n->u.ch.embells = parse_embell_list(p);

    return n;
}

/* --- parse_tmpl --- */
static MtefNode *parse_tmpl(MtefParser *p, int options)
{
    if (options & OPT_NUDGE) p_read_nudge(p);

    int selector = p_read_byte(p);
    if (selector < 0) selector = 0;
    int variation = p_read_byte(p);
    if (variation < 0) variation = 0;

    if (variation >= 128) {
        int v2 = p_read_byte(p);
        if (v2 >= 0) variation = ((variation - 128) << 8) | v2;
    }

    int num_slots = get_template_slot_count(selector, variation);

    MtefNode *n = new_node(p->arena, NODE_TMPL);
    if (!n) return NULL;
    n->u.tmpl.selector = selector;
    n->u.tmpl.variation = variation;
    n->u.tmpl.orig_variation = variation;

    for (int i = 0; i < num_slots; i++) {
        /* EQNEDT32 native tmSCRIPT: if slot[0] has SIZE_SUB, slot[1] contains
         * parent-scope data (fences, BigOp display data) that was consumed by
         * parse_object_list. Save position before reading slot[1] so we can
         * rewind if slot[0] indicates native subscript format. */
        int saved_pos = p->pos;
        int is_script_slot1 = (i == 1 && (selector == tmSCRIPT || selector == tmLSCRIPT) &&
                               n->u.tmpl.slots.count == 1);
        int s0_has_size_sub = 0;
        if (is_script_slot1) {
            MtefNode *s0 = n->u.tmpl.slots.items[0];
            for (int k = 0; k < s0->u.line.children.count; k++) {
                MtefNode *c = s0->u.line.children.items[k];
                if (c->type == NODE_SIZE && c->u.size.size_type == SIZETYPE_SUB) {
                    s0_has_size_sub = 1; break;
                }
            }
        }

        /* Also check: slot[1]'s first record must be SIZE_FULL (= size reset after
         * subscript data, indicating parent-scope overflow, not actual superscript) */
        int next_is_size_full = 0;
        if (is_script_slot1 && s0_has_size_sub && saved_pos < p->len) {
            uint8_t peek = p->data[saved_pos];
            next_is_size_full = ((peek & 0x0F) == REC_FULL ||
                                 (peek == 0x0A));  /* SIZE_FULL shortcut */
        }
        if (is_script_slot1 && s0_has_size_sub && next_is_size_full) {
            /* Don't read slot[1]: restore position so parent gets these records */
            p->pos = saved_pos;
            MtefNode *slot = new_node(p->arena, NODE_LINE);
            if (slot) { slot->u.line.is_null = 1; nl_push(p->arena, &n->u.tmpl.slots, slot); }
            break;
        }

        MtefNode *slot = new_node(p->arena, NODE_LINE);
        if (!slot) break;
        int dbg_pos_before = p->pos;
        slot->u.line.children = parse_object_list(p);
        /* Check if slot is effectively null (empty children) */
        if (slot->u.line.children.count == 0)
            slot->u.line.is_null = 1;
        nl_push(p->arena, &n->u.tmpl.slots, slot);
    }

    /* EQNEDT32 native tmSCRIPT: slot[0] has SIZE_SUB (subscript data),
     * OR slot[0] is empty and slot[1] starts with SIZE_SUB (native superscript).
     * In native format, parse_object_list for slot[1] eats ALL remaining records
     * including content that belongs to the PARENT scope.
     * Fix: detect overflow in slot[1] and move non-display-data back to parent. */
    if ((selector == tmSCRIPT || selector == tmLSCRIPT) &&
        n->u.tmpl.slots.count == 2) {
        MtefNode *s0 = n->u.tmpl.slots.items[0];
        int has_size_sub = 0;
        /* Check slot[0] for SIZE_SUB (native subscript) */
        for (int i = 0; i < s0->u.line.children.count; i++) {
            MtefNode *c = s0->u.line.children.items[i];
            if (c->type == NODE_SIZE && c->u.size.size_type == SIZETYPE_SUB) {
                has_size_sub = 1; break;
            }
        }
        if (has_size_sub) {
            /* slot[1] has parent-scope data: flatten back into parent.
             * We do this by replacing slot[1] with an empty slot and restoring
             * the parse position to before slot[1] was read. */
            MtefNode *s1 = n->u.tmpl.slots.items[1];
            if (s1->u.line.children.count > 0) {
                /* Flatten: move slot[1] children back to parent stream.
                 * Since we can't un-read, we store them in the template and
                 * handle them in convert_tmpl output phase. */
                n->u.tmpl.slots.items[1] = new_node(p->arena, NODE_LINE);
                if (n->u.tmpl.slots.items[1]) {
                    n->u.tmpl.slots.items[1]->u.line.is_null = 1;
                }
                /* Store overflow in a special field: use display_lower as carrier.
                 * convert_tmpl for tmSCRIPT will output these as trailing siblings. */
                n->u.tmpl.display_lower = s1;
            }
        }
    }

    /* Fence templates: after slot[0], read trailing content (content LINE,
     * display chars) until the fence END into slot[1].  The fence END is
     * consumed by parse_object_list.  Slot[1] is then flattened back into
     * the parent scope by parse_object_list (see below).
     *
     * Only keep slot[1] if it actually contains display chars — this
     * confirms the END was a fence END.  Inside fractions (in_frac),
     * fences have no display chars and the END belongs to the fraction. */
    if (is_fence_selector(selector) && n->u.tmpl.slots.count == 1 &&
        p->pos < p->len && (p->data[p->pos] & 0x0F) != REC_END) {
        int saved_pos = p->pos;
        MtefNode *trail = new_node(p->arena, NODE_LINE);
        if (trail) {
            trail->u.line.children = parse_object_list(p);
            /* Verify slot[1] ENDS with display chars (confirms fence END).
             * If the last child is NOT a display char, the END that
             * terminated slot[1] was a scope END, not a fence END. */
            int has_display = 0;
            if (trail->u.line.children.count > 0) {
                MtefNode *last = trail->u.line.children.items[
                    trail->u.line.children.count - 1];
                if (last->type == NODE_CHAR && last->u.ch.typeface == 22)
                    has_display = 1;
            }
            if (has_display) {
                nl_push(p->arena, &n->u.tmpl.slots, trail);
            } else {
                /* No display chars — the END was not a fence END.
                 * Restore position so the parent can re-read these records. */
                p->pos = saved_pos;
            }
        }
    }

    /* EQNEDT32 ROOT pattern: slot[0] empty -> extra slot */
    if (selector == tmROOT && n->u.tmpl.slots.count > 0) {
        MtefNode *first = n->u.tmpl.slots.items[0];
        int is_empty = first->u.line.is_null || first->u.line.children.count == 0;
        if (!is_empty) {
            /* Check all children are null lines */
            is_empty = 1;
            for (int i = 0; i < first->u.line.children.count; i++) {
                MtefNode *c = first->u.line.children.items[i];
                if (c->type != NODE_LINE || !c->u.line.is_null) {
                    is_empty = 0; break;
                }
            }
        }
        if (is_empty && p->pos < p->len) {
            int peek = p->data[p->pos] & 0x0F;
            if (peek != REC_END) {
                MtefNode *extra = new_node(p->arena, NODE_LINE);
                if (extra) {
                    extra->u.line.children = parse_object_list(p);
                    nl_push(p->arena, &n->u.tmpl.slots, extra);
                }
            }
        }
    }

    return n;
}

/* --- parse_pile --- */
static MtefNode *parse_pile(MtefParser *p, int options)
{
    int has_nudge = options & OPT_NUDGE;
    int has_ruler = options & 0x02;

    if (has_nudge) p_read_nudge(p);

    int halign = p_read_byte(p);
    if (halign < 0) halign = 0;

    if (has_ruler) {
        /* skip ruler data */
        int n_stops = p_read_byte(p);
        if (n_stops > 0) {
            for (int i = 0; i < n_stops; i++) {
                p_read_byte(p);   /* tab stop type */
                p_read_uint16(p); /* tab stop offset */
            }
        }
    }

    MtefNode *pile = new_node(p->arena, NODE_PILE);
    if (!pile) return NULL;
    pile->u.pile.halign = halign;

    int saw_end = 0;
    while (p->pos < p->len) {
        int tag = p_read_byte(p);
        if (tag < 0) break;
        int rec_type = tag & 0x0F;
        int opts = (tag >> 4) & 0x0F;

        if (rec_type == REC_END) {
            if (p->pos < p->len) {
                int next_byte = p->data[p->pos];
                int next_rtype = next_byte & 0x0F;
                /* Double-END or outer END: PILE is done */
                if (next_rtype == REC_END) break;
                /* Next is a LINE record: more rows follow */
                if (next_rtype == REC_LINE) { saw_end = 1; continue; }
                /* Content record (CHAR/TMPL/PILE/MATRIX/EMBELL/RULER/FONT)
                 * after an END, with at least one row already parsed:
                 * this END terminates the PILE; outer stream content follows. */
                if (pile->u.pile.lines.count > 0 &&
                    next_rtype >= REC_CHAR && next_rtype <= REC_FONT) {
                    break;
                }
            } else {
                /* At stream end: PILE is done */
                break;
            }
            /* Fall through: SIZE/typesize records or other edge cases —
             * use original 256-byte lookahead heuristic */
            int has_more = 0;
            size_t scan_end = p->pos + 256;
            if (scan_end > p->len) scan_end = p->len;
            for (size_t i = p->pos; i < scan_end; i++) {
                if (p->data[i] != 0x00) { has_more = 1; break; }
            }
            if (!has_more) break;
            saw_end = 1;
            continue;
        }

        if (rec_type == REC_LINE) {
            /* EQNEDT32: template with empty slot -> append next LINE */
            int append_to_current = 0;
            if (pile->u.pile.lines.count > 0) {
                MtefNode *cur = pile->u.pile.lines.items[pile->u.pile.lines.count - 1];
                if (cur->u.line.children.count > 0) {
                    MtefNode *last = cur->u.line.children.items[cur->u.line.children.count - 1];
                    if (last->type == NODE_TMPL &&
                        (is_fence_selector(last->u.tmpl.selector) ||
                         is_bigop_selector(last->u.tmpl.selector)) &&
                        last->u.tmpl.slots.count > 0) {
                        MtefNode *s0 = last->u.tmpl.slots.items[0];
                        if (s0->u.line.is_null || s0->u.line.children.count == 0)
                            append_to_current = 1;
                    }
                }
            }

            MtefNode *line = parse_line(p, opts);
            if (line) {
                if (line->u.line.is_null && !saw_end && pile->u.pile.lines.count > 0) {
                    /* skip null LINE without preceding END */
                } else if (append_to_current) {
                    MtefNode *cur = pile->u.pile.lines.items[pile->u.pile.lines.count - 1];
                    nl_push(p->arena, &cur->u.line.children, line);
                } else {
                    nl_push(p->arena, &pile->u.pile.lines, line);
                }
            }
            saw_end = 0;
        } else {
            MtefNode *node = parse_record(p, rec_type, opts);
            if (node) {
                if (saw_end && pile->u.pile.lines.count > 0) {
                    /* After an END separator, non-LINE records belong to
                     * the NEXT pile line (continuation content like SIZE_FULL, CHAR).
                     * Append to the last line only if it was just created for this segment. */
                    MtefNode *cur = pile->u.pile.lines.items[pile->u.pile.lines.count - 1];
                    nl_push(p->arena, &cur->u.line.children, node);
                } else if (pile->u.pile.lines.count > 0) {
                    MtefNode *cur = pile->u.pile.lines.items[pile->u.pile.lines.count - 1];
                    nl_push(p->arena, &cur->u.line.children, node);
                } else {
                    MtefNode *new_line = new_node(p->arena, NODE_LINE);
                    if (new_line) {
                        nl_push(p->arena, &new_line->u.line.children, node);
                        nl_push(p->arena, &pile->u.pile.lines, new_line);
                    }
                }
                saw_end = 0;
            }
        }
    }
    return pile;
}

/* --- parse_matrix --- */
static MtefNode *parse_matrix(MtefParser *p, int options)
{
    if (options & OPT_NUDGE) p_read_nudge(p);

    p_read_byte(p); /* valign */
    p_read_byte(p); /* hjust */
    p_read_byte(p); /* vjust */
    int rows = p_read_byte(p);
    int cols = p_read_byte(p);
    if (rows < 0) rows = 0;
    if (cols < 0) cols = 0;

    /* EQNEDT32: 2-byte partition */
    p_read_byte(p);
    p_read_byte(p);

    MtefNode *mat = new_node(p->arena, NODE_MATRIX);
    if (!mat) return NULL;
    mat->u.matrix.rows = rows;
    mat->u.matrix.cols = cols;

    int total = rows * cols;
    int elem_idx = 0;
    int size_skips = 0;
    while (elem_idx < total && p->pos < p->len) {
        int tag = p_read_byte(p);
        if (tag < 0) break;
        int rec_type = tag & 0x0F;
        int opts = (tag >> 4) & 0x0F;

        if (rec_type >= REC_FULL && rec_type <= REC_SUBSYM) {
            size_skips++;
            if (size_skips > 100) break;
            continue;
        }
        size_skips = 0;

        MtefNode *elem = NULL;
        if (rec_type == REC_LINE) {
            elem = parse_line(p, opts);
        } else if (rec_type == REC_END) {
            elem = new_node(p->arena, NODE_LINE);
            if (elem) elem->u.line.is_null = 1;
        } else {
            MtefNode *child = parse_record(p, rec_type, opts);
            elem = new_node(p->arena, NODE_LINE);
            if (elem && child) nl_push(p->arena, &elem->u.line.children, child);
        }
        if (!elem) {
            elem = new_node(p->arena, NODE_LINE);
            if (elem) elem->u.line.is_null = 1;
        }
        nl_push(p->arena, &mat->u.matrix.elements, elem);
        elem_idx++;
    }

    /* pad missing elements */
    while (mat->u.matrix.elements.count < total) {
        MtefNode *pad = new_node(p->arena, NODE_LINE);
        if (pad) { pad->u.line.is_null = 1; nl_push(p->arena, &mat->u.matrix.elements, pad); }
        else break;
    }

    return mat;
}

/* --- parse_font --- */
static MtefNode *parse_font(MtefParser *p, int options)
{
    (void)options;
    int font_index = p_read_byte(p);
    int style = p_read_byte(p);

    MtefNode *n = new_node(p->arena, NODE_FONT);
    if (!n) { p_read_null_string(p); return NULL; }
    n->u.font.font_index = font_index;
    n->u.font.style = style;
    p_read_null_string_to(p, n->u.font.name, sizeof(n->u.font.name));
    return n;
}

/* --- parse_size --- */
static MtefNode *parse_size(MtefParser *p, int options)
{
    (void)options;
    int lsize = p_read_byte(p);
    if (lsize < 0) return NULL;
    if (lsize == 101) { p_read_int16(p); }
    else if (lsize == 100) { p_read_byte(p); p_read_int16(p); }
    else { p_read_byte(p); }
    return NULL; /* sizes are ignored in conversion */
}

/* --- parse_typesize --- */
static MtefNode *parse_typesize(MtefParser *p, int rec_type)
{
    (void)p;
    MtefNode *n = new_node(p->arena, NODE_SIZE);
    if (!n) return NULL;
    switch (rec_type) {
        case REC_FULL:   n->u.size.size_type = SIZETYPE_FULL; break;
        case REC_SUB:    n->u.size.size_type = SIZETYPE_SUB; break;
        case REC_SUB2:   n->u.size.size_type = SIZETYPE_SUB2; break;
        case REC_SYM:    n->u.size.size_type = SIZETYPE_SYM; break;
        case REC_SUBSYM: n->u.size.size_type = SIZETYPE_SUBSYM; break;
        default:         n->u.size.size_type = SIZETYPE_FULL; break;
    }
    return n;
}

/* --- parse_ruler --- */
static void parse_ruler_data(MtefParser *p)
{
    int n_stops = p_read_byte(p);
    if (n_stops < 0) return;
    for (int i = 0; i < n_stops; i++) {
        p_read_byte(p);
        p_read_uint16(p);
    }
}

/* --- parse_record dispatch --- */
static MtefNode *parse_record(MtefParser *p, int rec_type, int options)
{
    switch (rec_type) {
    case REC_LINE:   return parse_line(p, options);
    case REC_CHAR:   return parse_char(p, options);
    case REC_TMPL:   return parse_tmpl(p, options);
    case REC_PILE:   return parse_pile(p, options);
    case REC_MATRIX: return parse_matrix(p, options);
    case REC_EMBELL: return parse_embell(p, options);
    case REC_RULER:  parse_ruler_data(p); return NULL;
    case REC_FONT:   return parse_font(p, options);
    case REC_SIZE:   parse_size(p, options); return NULL;
    case REC_FULL: case REC_SUB: case REC_SUB2: case REC_SYM: case REC_SUBSYM:
        return parse_typesize(p, rec_type);
    default: return NULL;
    }
}

/* --- parse (entry point) --- */
static MtefNode *mtef_parse(const uint8_t *data, size_t len, Eq2TexArena *arena, int *out_prod_ver)
{
    if (len < 5) return NULL;

    MtefParser p;
    p.data = data;
    p.len = len;
    p.pos = 0;
    p.arena = arena;

    /* 5-byte MTEF header */
    p_read_byte(&p); /* version */
    p_read_byte(&p); /* platform */
    p_read_byte(&p); /* product */
    int prod_ver = p_read_byte(&p);
    p_read_byte(&p); /* prod_subver */
    p.prod_ver = prod_ver;
    if (out_prod_ver) *out_prod_ver = prod_ver;

    /* Parse all records (skip intermediate REC_END) */
    NodeList children = {NULL, 0, 0};
    while (p.pos < p.len) {
        int tag = p_read_byte(&p);
        if (tag < 0) break;
        int rec_type = tag & 0x0F;
        int options = (tag >> 4) & 0x0F;

        if (rec_type == REC_END) {
            int has_more = 0;
            size_t scan_end = p.pos + 256;
            if (scan_end > p.len) scan_end = p.len;
            for (size_t i = p.pos; i < scan_end; i++) {
                if (p.data[i] != 0x00) { has_more = 1; break; }
            }
            if (!has_more) break;
            continue;
        }

        MtefNode *node = parse_record(&p, rec_type, options);
        if (node) nl_push(arena, &children, node);
    }

    /* EQNEDT32: flatten first LINE into root */
    NodeList flat = {NULL, 0, 0};
    int first_line_done = 0;
    for (int i = 0; i < children.count; i++) {
        MtefNode *c = children.items[i];
        if (!first_line_done && c->type == NODE_LINE && !c->u.line.is_null &&
            c->type != NODE_PILE) {
            for (int j = 0; j < c->u.line.children.count; j++)
                nl_push(arena, &flat, c->u.line.children.items[j]);
            first_line_done = 1;
        } else {
            nl_push(arena, &flat, c);
        }
    }

    MtefNode *root = new_node(arena, NODE_LINE);
    if (!root) return NULL;
    root->u.line.children = flat;
    return root;
}

/* ============================================================
 * LaTeX converter -- helpers
 * ============================================================ */

static int is_slot_empty(MtefNode *slot)
{
    if (!slot) return 1;
    if (slot->type != NODE_LINE) return 0;
    if (slot->u.line.is_null) return 1;
    if (slot->u.line.children.count == 0) return 1;
    /* all children are null lines? */
    for (int i = 0; i < slot->u.line.children.count; i++) {
        MtefNode *c = slot->u.line.children.items[i];
        if (c->type != NODE_LINE || !c->u.line.is_null) return 0;
    }
    return 1;
}

static int is_bigop_display_char(MtefNode *node, int prod_ver)
{
    if (!node || node->type != NODE_CHAR) return 0;
    uint16_t code = node->u.ch.char_code;
    if (set_contains(BIGOP_DISPLAY_CHARS, BIGOP_DISPLAY_N, code)) return 1;
    if (prod_ver >= 10 && code > 0xFF) {
        uint8_t hi = (code >> 8) & 0xFF;
        if (hi >= 0x80 || hi < 0x20) return 1;
    }
    return 0;
}

static int is_fence_display_char(MtefNode *node, int prod_ver)
{
    if (!node || node->type != NODE_CHAR) return 0;
    uint16_t code = node->u.ch.char_code;
    if (set_contains(FENCE_DISPLAY_CHARS, FENCE_DISPLAY_N, code)) return 1;
    if (node->u.ch.typeface == 22 &&
        !set_contains(BIGOP_DISPLAY_CHARS, BIGOP_DISPLAY_N, code)) return 1;
    if (prod_ver >= 10 && code > 0xFF) {
        uint8_t hi = (code >> 8) & 0xFF;
        if (set_contains(FENCE_DISPLAY_CHARS, FENCE_DISPLAY_N, hi)) return 1;
    }
    return 0;
}

static int is_empty_line(MtefNode *node)
{
    return node && node->type == NODE_LINE &&
           (node->u.line.is_null || node->u.line.children.count == 0);
}

/* escape a single ASCII char for math mode */
static void append_escaped_char(StringBuilder *sb, char c)
{
    switch (c) {
    case '#': case '$': case '%': case '&': case '_': case '{': case '}':
        sb_append_char(sb, '\\'); sb_append_char(sb, c); break;
    case '~': sb_append(sb, "\\sim "); break;
    case '^': sb_append(sb, "\\hat{}"); break;
    default:  sb_append_char(sb, c); break;
    }
}

/* ============================================================
 * char_to_latex
 * ============================================================ */

static void char_to_latex(int typeface, uint16_t code, int prod_ver, StringBuilder *sb)
{
    const char *s;

    /* Symbol font */
    if (typeface == TF_SYMBOL) {
        s = map_lookup(UNICODE_MAP, UNICODE_MAP_N, code);
        if (s) { sb_append(sb, s); return; }
        s = map_lookup(SYMBOL_MAP, SYMBOL_MAP_N, code);
        if (s) { sb_append(sb, s); return; }
        if (code >= 0x20 && code <= 0x7E) { sb_append_char(sb, (char)code); return; }
        { char buf[32]; sprintf(buf, "\\symbol{%u}", code); sb_append(sb, buf); return; }
    }

    /* MTEXTRA font */
    if (typeface == TF_MTEXTRA) {
        s = map_lookup(MTEXTRA_MAP, MTEXTRA_MAP_N, code);
        if (s) { sb_append(sb, s); return; }
    }

    /* Unicode lookup */
    s = map_lookup(UNICODE_MAP, UNICODE_MAP_N, code);
    if (s) { sb_append(sb, s); return; }

    /* Greek lowercase */
    if (typeface == TF_LCGREEK) {
        s = map_lookup(SYMBOL_MAP, SYMBOL_MAP_N, code);
        if (s) { sb_append(sb, s); return; }
        s = map_lookup(UNICODE_MAP, UNICODE_MAP_N, code);
        if (s) { sb_append(sb, s); return; }
        if (code >= 'a' && code <= 'z') { sb_append(sb, GREEK_LOWER[code - 'a']); return; }
        if (code >= 0x20 && code < 0x7F) { sb_append_char(sb, (char)code); return; }
        return;
    }

    /* Greek uppercase */
    if (typeface == TF_UCGREEK) {
        s = map_lookup(SYMBOL_MAP, SYMBOL_MAP_N, code);
        if (s) { sb_append(sb, s); return; }
        s = map_lookup(UNICODE_MAP, UNICODE_MAP_N, code);
        if (s) { sb_append(sb, s); return; }
        if (code >= 'A' && code <= 'Z') { sb_append(sb, GREEK_UPPER[code - 'A']); return; }
        if (code >= 0x20 && code < 0x7F) { sb_append_char(sb, (char)code); return; }
        return;
    }

    /* Vector (bold) */
    if (typeface == TF_VECTOR) {
        if ((code >= 'A' && code <= 'Z') || (code >= 'a' && code <= 'z')) {
            sb_append(sb, "\\mathbf{");
            sb_append_char(sb, (char)code);
            sb_append_char(sb, '}');
            return;
        }
        if (code >= 0x20 && code < 0x7F) { sb_append_char(sb, (char)code); return; }
        return;
    }

    /* Function */
    if (typeface == TF_FUNCTION) {
        if (code == 0x5E) { sb_append(sb, " \\wedge "); return; }
        if (code >= 0x20 && code < 0x7F) { sb_append_char(sb, (char)code); return; }
        return;
    }

    /* Text */
    if (typeface == TF_TEXT) {
        if (code >= 0x20 && code < 0x7F) {
            char c = (char)code;
            if (c == '#' || c == '$' || c == '%' || c == '&' ||
                c == '_' || c == '{' || c == '}') {
                sb_append_char(sb, '\\');
                sb_append_char(sb, c);
                return;
            }
            sb_append(sb, "\\text{");
            sb_append_char(sb, c);
            sb_append_char(sb, '}');
            return;
        }
        /* Non-ASCII text: encode Unicode code point as UTF-8 */
        if (code >= 0x80) {
            char utf8[5];
            int len = 0;
            if (code < 0x800) {
                utf8[len++] = (char)(0xC0 | (code >> 6));
                utf8[len++] = (char)(0x80 | (code & 0x3F));
            } else {
                utf8[len++] = (char)(0xE0 | (code >> 12));
                utf8[len++] = (char)(0x80 | ((code >> 6) & 0x3F));
                utf8[len++] = (char)(0x80 | (code & 0x3F));
            }
            utf8[len] = '\0';
            sb_append(sb, "\\text{");
            sb_append(sb, utf8);
            sb_append_char(sb, '}');
            return;
        }
        return;
    }

    /* Variable, Number, etc. - standard ASCII */
    if (code >= 0x20 && code < 0x7F) {
        append_escaped_char(sb, (char)code);
        return;
    }

    /* v10.2 packed char codes */
    if (prod_ver >= 10 && code > 0xFF) {
        uint8_t hi = (code >> 8) & 0xFF;
        uint8_t lo = code & 0xFF;

        if (lo == 0x03) {
            uint16_t candidate = 0x0300 | hi;
            s = map_lookup(UNICODE_MAP, UNICODE_MAP_N, candidate);
            if (s) { sb_append(sb, s); return; }
        }

        if (hi >= 0x20 && hi < 0x7F) {
            append_escaped_char(sb, (char)hi);
            return;
        }

        /* Try Unicode prefix restoration */
        static const uint8_t prefixes[] = {0x22, 0x03, 0x21, 0x00};
        for (int i = 0; i < 4; i++) {
            uint16_t candidate = ((uint16_t)prefixes[i] << 8) | hi;
            s = map_lookup(UNICODE_MAP, UNICODE_MAP_N, candidate);
            if (s) { sb_append(sb, s); return; }
        }

        s = map_lookup(SYMBOL_HIBYTE_MAP, SYMBOL_HIBYTE_MAP_N, hi);
        if (s) { sb_append(sb, s); return; }

        if (hi < 0x20) return; /* control char */
    }
}

/* ============================================================
 * convert_char
 * ============================================================ */

static void convert_char(MtefNode *node, int prod_ver, StringBuilder *sb)
{
    StringBuilder inner;
    sb_init(&inner);
    char_to_latex(node->u.ch.typeface, node->u.ch.char_code, prod_ver, &inner);

    /* Apply embellishments */
    for (int i = 0; i < node->u.ch.embells.count; i++) {
        MtefNode *e = node->u.ch.embells.items[i];
        if (e->type != NODE_EMBELL) continue;
        int et = e->u.embell.embell_type;
        if (et >= 0 && et < EMBELL_MAP_N) {
            const char *pfx = EMBELL_MAP[et].prefix;
            const char *sfx = EMBELL_MAP[et].suffix;
            if (pfx[0] || sfx[0]) {
                StringBuilder wrap;
                sb_init(&wrap);
                sb_append(&wrap, pfx);
                sb_append(&wrap, inner.buf);
                sb_append(&wrap, sfx);
                sb_free(&inner);
                inner = wrap;
            }
        }
    }

    sb_append(sb, inner.buf);
    sb_free(&inner);
}

/* ============================================================
 * convert_tmpl helpers (slot conversion)
 * ============================================================ */

static void convert_slot(MtefNode *tmpl, int index, int prod_ver, StringBuilder *sb)
{
    if (index < tmpl->u.tmpl.slots.count) {
        MtefNode *slot = tmpl->u.tmpl.slots.items[index];
        convert_node(slot, prod_ver, sb);
    }
}

static void convert_slot_to(MtefNode *tmpl, int index, int prod_ver, StringBuilder *out)
{
    StringBuilder tmp;
    sb_init(&tmp);
    convert_slot(tmpl, index, prod_ver, &tmp);
    sb_append(out, tmp.buf);
    sb_free(&tmp);
}

/* Get slot content as a string (caller frees) */
static char *slot_str(MtefNode *tmpl, int index, int prod_ver)
{
    StringBuilder tmp;
    sb_init(&tmp);
    convert_slot(tmpl, index, prod_ver, &tmp);
    return sb_detach(&tmp);
}

/* ============================================================
 * EQNEDT32 slot parsers
 * ============================================================ */

/* _parse_eqnedt32_bigop_slot: Extract integrand, lower, upper from slot[1] */
static void parse_eq_bigop_slot(MtefNode *slot, int prod_ver,
                                char **out_integrand, char **out_lower, char **out_upper)
{
    *out_integrand = *out_lower = *out_upper = NULL;
    if (!slot || slot->type != NODE_LINE || slot->u.line.is_null ||
        slot->u.line.children.count == 0) return;

    MtefNode *integrand_line = NULL;
    MtefNode *limit_lines[8];
    int limit_count = 0;
    int in_limits = 0;

    for (int i = 0; i < slot->u.line.children.count; i++) {
        MtefNode *child = slot->u.line.children.items[i];
        if (child->type == NODE_SIZE) {
            if (child->u.size.size_type == SIZETYPE_SUB) in_limits = 1;
            else if (child->u.size.size_type == SIZETYPE_SYM) break;
            continue;
        }
        if (child->type == NODE_LINE) {
            if (in_limits) {
                if (!child->u.line.is_null && limit_count < 8)
                    limit_lines[limit_count++] = child;
            } else if (!integrand_line)
                integrand_line = child;
            else if (!child->u.line.is_null && limit_count < 8)
                limit_lines[limit_count++] = child;
            continue;
        }
        if (child->type == NODE_CHAR && in_limits) break;
    }

    if (integrand_line) {
        StringBuilder sb; sb_init(&sb);
        convert_node(integrand_line, prod_ver, &sb);
        /* Also convert trailing CHAR/TMPL nodes after integrand (e.g. "dx") */
        int found_integ = 0;
        for (int i = 0; i < slot->u.line.children.count; i++) {
            MtefNode *child = slot->u.line.children.items[i];
            if (child == integrand_line) { found_integ = 1; continue; }
            if (!found_integ) continue;
            if (child->type == NODE_SIZE) break;  /* limits start */
            if (child->type == NODE_CHAR || child->type == NODE_TMPL)
                convert_node(child, prod_ver, &sb);
            else if (child->type == NODE_LINE) break;  /* another LINE = limit */
        }
        *out_integrand = sb_detach(&sb);
    }
    if (limit_count >= 1) {
        StringBuilder sb; sb_init(&sb);
        convert_node(limit_lines[0], prod_ver, &sb);
        *out_lower = sb_detach(&sb);
    }
    if (limit_count >= 2) {
        StringBuilder sb; sb_init(&sb);
        convert_node(limit_lines[1], prod_ver, &sb);
        *out_upper = sb_detach(&sb);
    }
}

/* _parse_eqnedt32_brace_slot: Extract content, label from slot[1] */
static void parse_eq_brace_slot(MtefNode *slot, int prod_ver,
                                char **out_content, char **out_label)
{
    *out_content = *out_label = NULL;
    if (!slot || slot->type != NODE_LINE || slot->u.line.is_null ||
        slot->u.line.children.count == 0) return;

    StringBuilder content_sb, label_sb;
    sb_init(&content_sb);
    sb_init(&label_sb);
    int in_label = 0;

    for (int i = 0; i < slot->u.line.children.count; i++) {
        MtefNode *child = slot->u.line.children.items[i];
        if (child->type == NODE_SIZE) {
            if (child->u.size.size_type == SIZETYPE_SUB) in_label = 1;
            else if (child->u.size.size_type == SIZETYPE_FULL ||
                     child->u.size.size_type == SIZETYPE_SYM) break;
            continue;
        }
        if (in_label) {
            if (child->type == NODE_LINE && child->u.line.is_null) continue;
            convert_node(child, prod_ver, &label_sb);
        } else {
            convert_node(child, prod_ver, &content_sb);
        }
    }
    *out_content = sb_detach(&content_sb);
    *out_label = sb_detach(&label_sb);
}

/* _parse_eqnedt32_lim_slot: Extract op_name, limit from slot[1] */
static void parse_eq_lim_slot(MtefNode *slot, int prod_ver,
                              char **out_opname, char **out_limit)
{
    *out_opname = *out_limit = NULL;
    if (!slot || slot->type != NODE_LINE || slot->u.line.is_null ||
        slot->u.line.children.count == 0) return;

    StringBuilder op_sb, lim_sb;
    sb_init(&op_sb);
    sb_init(&lim_sb);
    int in_limit = 0;

    for (int i = 0; i < slot->u.line.children.count; i++) {
        MtefNode *child = slot->u.line.children.items[i];
        if (child->type == NODE_SIZE) {
            if (child->u.size.size_type == SIZETYPE_SUB) in_limit = 1;
            else if (child->u.size.size_type == SIZETYPE_SYM ||
                     child->u.size.size_type == SIZETYPE_FULL) break;
            continue;
        }
        if (child->type == NODE_LINE && child->u.line.is_null) continue;
        if (in_limit) convert_node(child, prod_ver, &lim_sb);
        else convert_node(child, prod_ver, &op_sb);
    }
    *out_opname = sb_detach(&op_sb);
    *out_limit = sb_detach(&lim_sb);
}

/* _parse_eqnedt32_root_slot: Extract radicand, nth index from slot[1] */
static void parse_eq_root_slot(MtefNode *slot, int prod_ver,
                               char **out_radicand, char **out_nth)
{
    *out_radicand = *out_nth = NULL;
    if (!slot || slot->type != NODE_LINE || slot->u.line.is_null ||
        slot->u.line.children.count == 0) return;

    StringBuilder rad_sb, nth_sb;
    sb_init(&rad_sb);
    sb_init(&nth_sb);
    int in_nth = 0;

    for (int i = 0; i < slot->u.line.children.count; i++) {
        MtefNode *child = slot->u.line.children.items[i];
        if (child->type == NODE_SIZE) {
            if (child->u.size.size_type == SIZETYPE_SUB) in_nth = 1;
            else if (child->u.size.size_type == SIZETYPE_SYM) break;
            continue;
        }
        if (in_nth) {
            if (child->type == NODE_LINE && !child->u.line.is_null) {
                convert_node(child, prod_ver, &nth_sb);
            } else if (child->type == NODE_CHAR) break;
        } else {
            convert_node(child, prod_ver, &rad_sb);
        }
    }
    *out_radicand = sb_detach(&rad_sb);
    *out_nth = sb_detach(&nth_sb);
}

/* _parse_eqnedt32_dirac: Extract bra, ket from slot[1] */
static void parse_eq_dirac(MtefNode *slot, int prod_ver,
                           char **out_bra, char **out_ket)
{
    *out_bra = *out_ket = NULL;
    if (!slot || slot->type != NODE_LINE || slot->u.line.is_null ||
        slot->u.line.children.count == 0) return;

    MtefNode *lines[8];
    int nlines = 0;
    for (int i = 0; i < slot->u.line.children.count; i++) {
        MtefNode *c = slot->u.line.children.items[i];
        if (c->type == NODE_LINE && !c->u.line.is_null && nlines < 8)
            lines[nlines++] = c;
        else if (c->type == NODE_CHAR && c->u.ch.typeface == 22)
            break;
    }
    if (nlines >= 1) { StringBuilder sb; sb_init(&sb); convert_node(lines[0], prod_ver, &sb); *out_bra = sb_detach(&sb); }
    if (nlines >= 2) { StringBuilder sb; sb_init(&sb); convert_node(lines[1], prod_ver, &sb); *out_ket = sb_detach(&sb); }
}

/* is_display_data_slot: slot contains only SizeNode, null LINE, display chars.
 * Requires at least one SIZE or null-LINE marker to distinguish true display
 * data from a real content slot that happens to contain only BIGOP_DISPLAY_CHARs
 * (e.g. a standalone \int that was "stolen" into an extra TMPL slot). */
static int is_display_data_slot(MtefNode *slot)
{
    if (!slot || slot->type != NODE_LINE) return 0;
    if (slot->u.line.is_null || slot->u.line.children.count == 0) return 1;
    int has_size_or_null_line = 0;
    for (int i = 0; i < slot->u.line.children.count; i++) {
        MtefNode *c = slot->u.line.children.items[i];
        if (c->type == NODE_SIZE) { has_size_or_null_line = 1; continue; }
        if (c->type == NODE_LINE && c->u.line.is_null) { has_size_or_null_line = 1; continue; }
        if (c->type == NODE_CHAR && set_contains(BIGOP_DISPLAY_CHARS, BIGOP_DISPLAY_N, c->u.ch.char_code)) continue;
        return 0;
    }
    /* Only treat as display-only when SIZE/null-LINE markers are present */
    return has_size_or_null_line;
}

/* ============================================================
 * convert_tmpl — Template conversion
 * ============================================================ */

static void convert_tmpl(MtefNode *node, int prod_ver, StringBuilder *sb)
{
    int sel = node->u.tmpl.selector;
    int var = node->u.tmpl.variation;

    /* === Fences === */
    /* If slot[0] is a PILE with a matrix-type halign (20-24, set by tex2mtef for
     * bmatrix/pmatrix/vmatrix etc.), output the matrix env directly without
     * \left..\right wrappers — the env name already encodes the bracket style. */
    if (sel >= tmANGLE && sel <= tmLPRB) {
        MtefNode *s0 = (node->u.tmpl.slots.count > 0) ? node->u.tmpl.slots.items[0] : NULL;
        if (s0 && s0->type == NODE_PILE && s0->u.pile.halign >= 20) {
            convert_node(s0, prod_ver, sb);
            return;
        }
    }
    if (sel == tmPAREN) {
        char *s = slot_str(node, 0, prod_ver);
        if (var == 1) { sb_append(sb, "\\left( "); sb_append(sb, s); sb_append(sb, " \\right."); }
        else if (var == 2) { sb_append(sb, "\\left. "); sb_append(sb, s); sb_append(sb, " \\right)"); }
        else { sb_append(sb, "\\left( "); sb_append(sb, s); sb_append(sb, " \\right)"); }
        free(s); return;
    }
    if (sel == tmBRACK) {
        char *s = slot_str(node, 0, prod_ver);
        if (var == 1) { sb_append(sb, "\\left[ "); sb_append(sb, s); sb_append(sb, " \\right."); }
        else if (var == 2) { sb_append(sb, "\\left. "); sb_append(sb, s); sb_append(sb, " \\right]"); }
        else { sb_append(sb, "\\left[ "); sb_append(sb, s); sb_append(sb, " \\right]"); }
        free(s); return;
    }
    if (sel == tmBRACE) {
        char *s = slot_str(node, 0, prod_ver);
        if (var == 1) { sb_append(sb, "\\left\\{ "); sb_append(sb, s); sb_append(sb, " \\right."); }
        else if (var == 2) { sb_append(sb, "\\left. "); sb_append(sb, s); sb_append(sb, " \\right\\}"); }
        else { sb_append(sb, "\\left\\{ "); sb_append(sb, s); sb_append(sb, " \\right\\}"); }
        free(s); return;
    }
    if (sel == tmANGLE) {
        char *s = slot_str(node, 0, prod_ver);
        if (var == 1) { sb_append(sb, "\\left\\langle "); sb_append(sb, s); sb_append(sb, " \\right."); }
        else if (var == 2) { sb_append(sb, "\\left. "); sb_append(sb, s); sb_append(sb, " \\right\\rangle "); }
        else { sb_append(sb, "\\left\\langle "); sb_append(sb, s); sb_append(sb, " \\right\\rangle "); }
        free(s); return;
    }
    if (sel == tmBAR) {
        char *s = slot_str(node, 0, prod_ver);
        if (var == 1) { sb_append(sb, "\\left| "); sb_append(sb, s); sb_append(sb, " \\right."); }
        else if (var == 2) { sb_append(sb, "\\left. "); sb_append(sb, s); sb_append(sb, " \\right|"); }
        else { sb_append(sb, "\\left| "); sb_append(sb, s); sb_append(sb, " \\right|"); }
        free(s); return;
    }
    if (sel == tmDBAR) {
        char *s = slot_str(node, 0, prod_ver);
        if (var == 1) { sb_append(sb, "\\left\\| "); sb_append(sb, s); sb_append(sb, " \\right."); }
        else if (var == 2) { sb_append(sb, "\\left. "); sb_append(sb, s); sb_append(sb, " \\right\\|"); }
        else { sb_append(sb, "\\left\\| "); sb_append(sb, s); sb_append(sb, " \\right\\|"); }
        free(s); return;
    }
    if (sel == tmFLOOR) {
        char *s = slot_str(node, 0, prod_ver);
        sb_append(sb, "\\left\\lfloor "); sb_append(sb, s); sb_append(sb, " \\right\\rfloor ");
        free(s); return;
    }
    if (sel == tmCEIL) {
        char *s = slot_str(node, 0, prod_ver);
        sb_append(sb, "\\left\\lceil "); sb_append(sb, s); sb_append(sb, " \\right\\rceil ");
        free(s); return;
    }
    if (sel == tmLBLB) { char *s = slot_str(node,0,prod_ver); sb_append(sb,"\\left[ "); sb_append(sb,s); sb_append(sb," \\right["); free(s); return; }
    if (sel == tmRBRB) { char *s = slot_str(node,0,prod_ver); sb_append(sb,"\\left] "); sb_append(sb,s); sb_append(sb," \\right]"); free(s); return; }
    if (sel == tmRBLB) { char *s = slot_str(node,0,prod_ver); sb_append(sb,"\\left] "); sb_append(sb,s); sb_append(sb," \\right["); free(s); return; }
    if (sel == tmLBRP) { char *s = slot_str(node,0,prod_ver); sb_append(sb,"\\left[ "); sb_append(sb,s); sb_append(sb," \\right)"); free(s); return; }
    if (sel == tmLPRB) { char *s = slot_str(node,0,prod_ver); sb_append(sb,"\\left( "); sb_append(sb,s); sb_append(sb," \\right]"); free(s); return; }

    /* === Fractions === */
    if (sel == tmFRACT || sel == tmSLFRACT) {
        char *n_str = NULL;
        char *d_str = NULL;

        /* EQNEDT32 packed format: slot[0] empty, slot[1] has num+denom as 2 LINEs.
         * Detect BEFORE calling slot_str because slot_str destructively modifies
         * the tree (Pass 1 fence merging). Re-converting the same LINE after
         * destructive modification produces wrong output (doubled content).
         * Also recover non-LINE children leaked from denom due to fence trailing
         * ENDs — these are siblings in slot[1] after the second LINE. */
        int packed = 0;
        if (node->u.tmpl.slots.count >= 2 && is_slot_empty(node->u.tmpl.slots.items[0])) {
            MtefNode *s1 = node->u.tmpl.slots.items[1];
            MtefNode *lines[8]; int nlines = 0;
            int line_idx[8];
            for (int i = 0; i < s1->u.line.children.count && nlines < 8; i++) {
                MtefNode *c = s1->u.line.children.items[i];
                if (c->type == NODE_LINE && !c->u.line.is_null) {
                    line_idx[nlines] = i;
                    lines[nlines++] = c;
                }
            }
            if (nlines >= 2) {
                packed = 1;
                StringBuilder ns, ds; sb_init(&ns); sb_init(&ds);
                convert_node(lines[0], prod_ver, &ns);
                convert_node(lines[1], prod_ver, &ds);
                /* Append leaked non-LINE children after the denom LINE.
                 * When a fence inside the denom has a trailing END, it truncates
                 * the denom LINE and content after the fence leaks to slot[1]. */
                for (int i = line_idx[1] + 1; i < s1->u.line.children.count; i++) {
                    MtefNode *c = s1->u.line.children.items[i];
                    if (c->type != NODE_LINE || (!c->u.line.is_null && c->u.line.children.count > 0))
                        convert_node(c, prod_ver, &ds);
                }
                n_str = sb_detach(&ns);
                d_str = sb_detach(&ds);
            }
        }
        if (!packed) {
            n_str = slot_str(node, 0, prod_ver);
            d_str = slot_str(node, 1, prod_ver);
        }

        if (sel == tmSLFRACT) {
            sb_append(sb, n_str); sb_append_char(sb, '/'); sb_append(sb, d_str);
        } else if (var == 1) {
            sb_append(sb, "\\tfrac{"); sb_append(sb, n_str); sb_append(sb, "}{"); sb_append(sb, d_str); sb_append_char(sb, '}');
        } else {
            sb_append(sb, "\\dfrac{"); sb_append(sb, n_str); sb_append(sb, "}{"); sb_append(sb, d_str); sb_append_char(sb, '}');
        }
        free(n_str); free(d_str); return;
    }

    /* === Roots === */
    if (sel == tmROOT) {
        char *content = slot_str(node, 0, prod_ver);
        char *nth = NULL;

        if ((!content || !content[0] || content[0] == ' ') && node->u.tmpl.slots.count >= 2) {
            free(content); content = NULL;
            if (var == 1) {
                char *rad = NULL; char *idx = NULL;
                parse_eq_root_slot(node->u.tmpl.slots.items[1], prod_ver, &rad, &idx);
                content = rad; nth = idx;
            } else {
                content = slot_str(node, 1, prod_ver);
            }
            /* suffix from extra slots */
            StringBuilder suffix; sb_init(&suffix);
            for (int i = 2; i < node->u.tmpl.slots.count; i++) {
                char *ex = slot_str(node, i, prod_ver);
                if (ex && ex[0]) sb_append(&suffix, ex);
                free(ex);
            }
            if (nth && nth[0]) {
                sb_append(sb, "\\sqrt["); sb_append(sb, nth); sb_append(sb, "]{");
                sb_append(sb, content ? content : ""); sb_append_char(sb, '}');
            } else {
                sb_append(sb, "\\sqrt{"); sb_append(sb, content ? content : ""); sb_append_char(sb, '}');
            }
            if (suffix.len > 0) sb_append(sb, suffix.buf);
            sb_free(&suffix);
            free(content); free(nth); return;
        }

        if (var == 1 && node->u.tmpl.slots.count >= 2 && !nth) {
            nth = slot_str(node, 1, prod_ver);
        }
        if (nth && nth[0]) {
            sb_append(sb, "\\sqrt["); sb_append(sb, nth); sb_append(sb, "]{");
            sb_append(sb, content ? content : ""); sb_append_char(sb, '}');
        } else {
            sb_append(sb, "\\sqrt{"); sb_append(sb, content ? content : ""); sb_append_char(sb, '}');
        }
        free(content); free(nth); return;
    }

    /* === Scripts === */
    if (sel == tmSCRIPT || sel == tmLSCRIPT) {
        char *sub_s = slot_str(node, 0, prod_ver);
        char *sup_s = NULL;  /* computed lazily to avoid tree mutation before EQNEDT32 check */

        /* EQNEDT32: slot[0] empty — extract sub/sup from slot[1] sub-LINEs.
         * NOTE: sup_s must NOT be computed before this check, because slot_str(node,1,...)
         * calls convert_line which mutates OBAR/fence slot[0] via Pass 1 fence-merging.
         * If that mutation happens first, the subsequent convert_node(lines[0]) call sees
         * the already-filled slot and outputs the content twice. */
        if ((!sub_s || !sub_s[0] || sub_s[0] == ' ') && node->u.tmpl.slots.count >= 2) {
            MtefNode *s1 = node->u.tmpl.slots.items[1];

            MtefNode *lines[8]; int nlines = 0;
            for (int i = 0; i < s1->u.line.children.count && nlines < 8; i++) {
                MtefNode *c = s1->u.line.children.items[i];
                if (c->type == NODE_LINE && !c->u.line.is_null) lines[nlines++] = c;
            }
            free(sub_s); sub_s = NULL;  /* sup_s is already NULL */
            if (var == 1 || var == 2) {
                /* Check if slot[1] has overflow: contains fence/bigop TMPL beyond subscript data */
                int has_overflow = 0;
                if (var == 1 && nlines >= 1) {
                    for (int ii = 0; ii < s1->u.line.children.count; ii++) {
                        MtefNode *c = s1->u.line.children.items[ii];
                        if (c->type == NODE_TMPL &&
                            (is_fence_selector(c->u.tmpl.selector) ||
                             is_display_tmpl_selector(c->u.tmpl.selector))) {
                            has_overflow = 1; break;
                        }
                    }
                }
                if (var == 1 && nlines >= 1 && has_overflow) {
                    /* var=1 (subscript only): first non-null LINE is the subscript.
                     * Remaining content in slot[1] is parent-scope overflow
                     * (fences, BigOp display data, equation continuation).
                     * Process it via slot_str with g_parent_bigop set. */
                    StringBuilder a; sb_init(&a);
                    convert_node(lines[0], prod_ver, &a);
                    sub_s = sb_detach(&a);
                    /* Mark the subscript LINE as null so slot_str won't output it again */
                    lines[0]->u.line.is_null = 1;
                    /* Process rest of slot[1] as overflow: run convert_line
                     * which will execute passes (including g_parent_bigop promotion) */
                    MtefNode *saved_p = g_parent_bigop;
                    sup_s = slot_str(node, 1, prod_ver);
                    g_parent_bigop = saved_p;
                } else if (nlines >= 2) {
                    StringBuilder a, b; sb_init(&a); sb_init(&b);
                    convert_node(lines[0], prod_ver, &a);
                    convert_node(lines[1], prod_ver, &b);
                    sub_s = sb_detach(&a); sup_s = sb_detach(&b);
                } else if (nlines == 1) {
                    StringBuilder a; sb_init(&a);
                    convert_node(lines[0], prod_ver, &a);
                    sub_s = sb_detach(&a); sup_s = _strdup("");
                } else {
                    StringBuilder a; sb_init(&a);
                    convert_node(s1, prod_ver, &a);
                    sub_s = sb_detach(&a); sup_s = _strdup("");
                }
            } else { /* var == 0: superscript only */
                StringBuilder a; sb_init(&a);
                convert_node(s1, prod_ver, &a);
                sup_s = sb_detach(&a); sub_s = _strdup("");
            }
        } else {
            /* Normal path: now safe to compute sup_s (tree not yet mutated) */
            sup_s = slot_str(node, 1, prod_ver);
        }

        const char *prefix = (sel == tmLSCRIPT) ? "{}" : "";
        if (var == 0) {
            sb_append(sb, prefix); sb_append(sb, "^{"); sb_append(sb, sup_s ? sup_s : ""); sb_append_char(sb, '}');
        } else if (var == 1) {
            sb_append(sb, prefix); sb_append(sb, "_{"); sb_append(sb, sub_s ? sub_s : ""); sb_append_char(sb, '}');
            /* sup_s may contain overflow content from EQNEDT32 slot[1] */
            if (sup_s && sup_s[0]) {
                char *t = sup_s; while (*t == ' ') t++;
                if (*t) sb_append(sb, t);
            }
        } else {
            sb_append(sb, prefix);
            sb_append(sb, "_{"); sb_append(sb, sub_s ? sub_s : ""); sb_append(sb, "}^{"); sb_append(sb, sup_s ? sup_s : ""); sb_append_char(sb, '}');
        }
        free(sub_s); free(sup_s);

        /* Output overflow data from EQNEDT32 native tmSCRIPT (display_lower carrier) */
        if (node->u.tmpl.display_lower && node->u.tmpl.display_lower->type == NODE_LINE) {
            convert_line(node->u.tmpl.display_lower, prod_ver, g_conv_arena, sb);
        }
        return;
    }

    /* === Bars === */
    if (sel == tmUBAR) {
        char *s = slot_str(node,0,prod_ver);
        if (var == 1) { sb_append(sb,"\\underline{\\underline{"); sb_append(sb,s); sb_append(sb,"}}"); }
        else { sb_append(sb,"\\underline{"); sb_append(sb,s); sb_append_char(sb,'}'); }
        free(s); return;
    }
    if (sel == tmOBAR) {
        char *s = slot_str(node,0,prod_ver);
        if (var == 1) { sb_append(sb,"\\overline{\\overline{"); sb_append(sb,s); sb_append(sb,"}}"); }
        else { sb_append(sb,"\\overline{"); sb_append(sb,s); sb_append_char(sb,'}'); }
        free(s); return;
    }

    /* === Arrows === */
    if (sel == tmLARROW) { char *s = slot_str(node,0,prod_ver); sb_append(sb,"\\overleftarrow{"); sb_append(sb,s); sb_append_char(sb,'}'); free(s); return; }
    if (sel == tmRARROW) { char *s = slot_str(node,0,prod_ver); sb_append(sb,"\\overrightarrow{"); sb_append(sb,s); sb_append_char(sb,'}'); free(s); return; }
    if (sel == tmBARROW) { char *s = slot_str(node,0,prod_ver); sb_append(sb,"\\overleftrightarrow{"); sb_append(sb,s); sb_append_char(sb,'}'); free(s); return; }

    /* === Horizontal braces === */
    if (sel == tmUHBRACE || sel == tmLHBRACE) {
        char *content = slot_str(node, 0, prod_ver);
        char *label = slot_str(node, 1, prod_ver);
        /* Case 1: slot[0]=empty, slot[1] has content+label (EQNEDT32 packed) */
        if ((!content || !content[0]) && node->u.tmpl.slots.count >= 2) {
            char *c2 = NULL, *l2 = NULL;
            parse_eq_brace_slot(node->u.tmpl.slots.items[1], prod_ver, &c2, &l2);
            if (c2 && c2[0]) { free(content); content = c2; } else free(c2);
            if (l2 && l2[0]) { free(label); label = l2; }
            else {
                free(l2);
                /* parse_eq_brace_slot splits on SIZE_SUB, but EQNEDT32 native
                 * overbrace/underbrace uses LINE boundaries: first non-null LINE
                 * = content, second non-null LINE = label.  Try LINE split. */
                MtefNode *s1 = node->u.tmpl.slots.items[1];
                if (s1 && s1->u.line.children.count >= 2) {
                    int nlines = 0, split_idx = -1;
                    for (int i = 0; i < s1->u.line.children.count; i++) {
                        MtefNode *c = s1->u.line.children.items[i];
                        if (c->type == NODE_LINE && !c->u.line.is_null &&
                            c->u.line.children.count > 0) {
                            nlines++;
                            if (nlines == 2) { split_idx = i; break; }
                        }
                    }
                    if (split_idx >= 0) {
                        /* Rebuild content from first LINE only */
                        StringBuilder csb; sb_init(&csb);
                        for (int i = 0; i < split_idx; i++) {
                            MtefNode *c = s1->u.line.children.items[i];
                            if (c->type == NODE_SIZE) continue;
                            convert_node(c, prod_ver, &csb);
                        }
                        free(content); content = sb_detach(&csb);
                        /* Label from second LINE */
                        StringBuilder lsb; sb_init(&lsb);
                        convert_node(s1->u.line.children.items[split_idx], prod_ver, &lsb);
                        free(label); label = sb_detach(&lsb);
                    }
                }
                if (!label || !label[0]) {
                    if (node->u.tmpl.display_lower) {
                        StringBuilder ds; sb_init(&ds);
                        convert_node(node->u.tmpl.display_lower, prod_ver, &ds);
                        free(label); label = sb_detach(&ds);
                    } else { if (!label) label = _strdup(""); }
                }
            }
        }
        if (sel == tmUHBRACE) {
            sb_append(sb, "\\overbrace{"); sb_append(sb, content); sb_append(sb, "}^{"); sb_append(sb, label); sb_append_char(sb, '}');
        } else {
            sb_append(sb, "\\underbrace{"); sb_append(sb, content); sb_append(sb, "}_{"); sb_append(sb, label); sb_append_char(sb, '}');
        }
        free(content); free(label); return;
    }

    /* === Wide accents === */
    if (sel == tmUARROW || sel == tmOARROW) {
        char *s = slot_str(node,0,prod_ver);
        if (var == 0) { sb_append(sb,"\\overleftarrow{"); sb_append(sb,s); sb_append_char(sb,'}'); }
        else if (var == 1) { sb_append(sb,"\\overrightarrow{"); sb_append(sb,s); sb_append_char(sb,'}'); }
        else { sb_append(sb,"\\overleftrightarrow{"); sb_append(sb,s); sb_append_char(sb,'}'); }
        free(s); return;
    }
    if (sel == tmOARC) {
        char *s = slot_str(node,0,prod_ver);
        sb_append(sb,"\\overset{\\frown}{"); sb_append(sb,s); sb_append_char(sb,'}');
        free(s); return;
    }

    /* === Big Operators — delegate to convert_bigop === */
    /* (inline here to avoid another function) */

    int has_lower = var & 0x01;
    int has_upper = var & 0x02;

    /* tmLIM */
    if (sel == tmLIM) {
        char *base = slot_str(node, 0, prod_ver);
        if ((!base || !base[0] || base[0] == ' ') && node->u.tmpl.slots.count >= 2) {
            char *opname = NULL, *limit = NULL;
            parse_eq_lim_slot(node->u.tmpl.slots.items[1], prod_ver, &opname, &limit);
            free(base);
            if (opname && opname[0]) {
                sb_append(sb, "\\operatorname{"); sb_append(sb, opname); sb_append_char(sb, '}');
            }
            if (limit && limit[0]) {
                sb_append(sb, "\\limits_{"); sb_append(sb, limit); sb_append_char(sb, '}');
            }
            free(opname); free(limit); return;
        }
        if (has_lower && has_upper) {
            sb_append(sb, "\\mathop{"); sb_append(sb, base); sb_append(sb, "}\\limits_{");
            char *lo = slot_str(node,1,prod_ver); char *hi = slot_str(node,2,prod_ver);
            sb_append(sb, lo); sb_append(sb, "}^{"); sb_append(sb, hi); sb_append_char(sb, '}');
            free(lo); free(hi);
        } else if (has_lower) {
            sb_append(sb, "\\mathop{"); sb_append(sb, base); sb_append(sb, "}\\limits_{");
            char *lo = slot_str(node,1,prod_ver); sb_append(sb, lo); sb_append_char(sb, '}'); free(lo);
        } else if (has_upper) {
            sb_append(sb, "\\mathop{"); sb_append(sb, base); sb_append(sb, "}\\limits^{");
            char *hi = slot_str(node,1,prod_ver); sb_append(sb, hi); sb_append_char(sb, '}'); free(hi);
        } else {
            sb_append(sb, "\\mathop{"); sb_append(sb, base); sb_append_char(sb, '}');
        }
        free(base); return;
    }

    /* tmLDIV */
    if (sel == tmLDIV) {
        char *s = slot_str(node,0,prod_ver);
        if (var == 0) { sb_append(sb, "\\overline{"); sb_append(sb, s); sb_append_char(sb, '}'); }
        else sb_append(sb, s);
        free(s); return;
    }

    /* tmDIRAC */
    if (sel == tmDIRAC) {
        int s0_empty = (node->u.tmpl.slots.count > 0) ? is_slot_empty(node->u.tmpl.slots.items[0]) : 1;
        if (s0_empty && var == 0 && node->u.tmpl.slots.count >= 2) {
            char *bra = NULL, *ket = NULL;
            parse_eq_dirac(node->u.tmpl.slots.items[1], prod_ver, &bra, &ket);
            sb_append(sb, "\\left\\langle "); sb_append(sb, bra ? bra : "");
            sb_append(sb, " \\middle| "); sb_append(sb, ket ? ket : "");
            sb_append(sb, " \\right\\rangle ");
            free(bra); free(ket); return;
        }
        if (var == 0) {
            char *a = slot_str(node,0,prod_ver), *b = slot_str(node,1,prod_ver);
            sb_append(sb, "\\left\\langle "); sb_append(sb,a); sb_append(sb, " \\middle| "); sb_append(sb,b); sb_append(sb, " \\right\\rangle ");
            free(a); free(b); return;
        }
        if (var == 1) { char *a = slot_str(node,0,prod_ver); sb_append(sb, "\\left\\langle "); sb_append(sb,a); sb_append(sb, " \\right|"); free(a); return; }
        if (var == 2) { char *a = slot_str(node,0,prod_ver); sb_append(sb, "\\left| "); sb_append(sb,a); sb_append(sb, " \\right\\rangle "); free(a); return; }
        return;
    }

    /* tmINTOP / tmSUMOP */
    if (sel == tmINTOP || sel == tmSUMOP) {
        int idx = 0; char *lo = NULL, *hi = NULL;
        if (has_lower) { lo = slot_str(node, idx, prod_ver); idx++; }
        if (has_upper) { hi = slot_str(node, idx, prod_ver); idx++; }
        char *sym = slot_str(node, idx, prod_ver);
        const char *style = (sel == tmSUMOP) ? "\\limits" : "\\nolimits";
        sb_append(sb, "\\mathop{"); sb_append(sb, sym); sb_append_char(sb, '}');
        if (has_lower && has_upper) {
            sb_append(sb, style); sb_append(sb, "_{"); sb_append(sb, lo); sb_append(sb, "}^{"); sb_append(sb, hi); sb_append_char(sb, '}');
        } else if (has_lower) {
            sb_append(sb, style); sb_append(sb, "_{"); sb_append(sb, lo); sb_append_char(sb, '}');
        } else if (has_upper) {
            sb_append(sb, style); sb_append(sb, "^{"); sb_append(sb, hi); sb_append_char(sb, '}');
        }
        sb_append_char(sb, ' ');
        free(lo); free(hi); free(sym); return;
    }

    /* Standard BigOps */
    {
        const char *op = NULL;
        const char *style = "";

        switch (sel) {
        case tmSINT: op="\\int"; break; case tmDINT: op="\\iint"; break; case tmTINT: op="\\iiint"; break;
        case tmSSINT: op="\\int"; style="\\limits"; break; case tmDSINT: op="\\iint"; style="\\limits"; break; case tmTSINT: op="\\iiint"; style="\\limits"; break;
        case tmSUM: op="\\sum"; style="\\limits"; break; case tmISUM: op="\\sum"; style="\\nolimits"; break;
        case tmPROD: op="\\prod"; style="\\limits"; break; case tmIPROD: op="\\prod"; style="\\nolimits"; break;
        case tmCOPROD: op="\\coprod"; style="\\limits"; break; case tmICOPROD: op="\\coprod"; style="\\nolimits"; break;
        case tmUNION: op="\\bigcup"; style="\\limits"; break; case tmIUNION: op="\\bigcup"; style="\\nolimits"; break;
        case tmINTER: op="\\bigcap"; style="\\limits"; break; case tmIINTER: op="\\bigcap"; style="\\nolimits"; break;
        }
        if (!op) {
            /* fallback: convert all slots */
            for (int i = 0; i < node->u.tmpl.slots.count; i++)
                convert_slot(node, i, prod_ver, sb);
            return;
        }

        /* Special contour integral variants — use orig_variation
         * to avoid collision with Pass 2 display data mutation */
        int ovar = node->u.tmpl.orig_variation;
        if (sel == tmSINT && ovar == 3) op = "\\oint";
        else if (sel == tmDINT && (ovar & 0x02) && !(ovar & 0x01)) op = "\\oiint";
        else if (sel == tmTINT && (ovar & 0x02) && !(ovar & 0x01)) op = "\\oiiint";

        /* Check _display_lower/_display_upper from convert_pile */
        if (node->u.tmpl.display_lower || node->u.tmpl.display_upper) {
            char *lo = NULL, *hi_s = NULL;
            if (node->u.tmpl.display_lower) {
                StringBuilder ts; sb_init(&ts);
                convert_node(node->u.tmpl.display_lower, prod_ver, &ts);
                lo = sb_detach(&ts);
            }
            if (node->u.tmpl.display_upper) {
                StringBuilder ts; sb_init(&ts);
                convert_node(node->u.tmpl.display_upper, prod_ver, &ts);
                hi_s = sb_detach(&ts);
            }
            /* integrand */
            int s0_empty = (node->u.tmpl.slots.count > 0) ? is_slot_empty(node->u.tmpl.slots.items[0]) : 0;
            MtefNode *saved_parent = g_parent_bigop;
            g_parent_bigop = node;
            char *integ;
            if (s0_empty && node->u.tmpl.slots.count >= 2)
                integ = slot_str(node, 1, prod_ver);
            else
                integ = slot_str(node, 0, prod_ver);
            g_parent_bigop = saved_parent;

            sb_append(sb, op);
            if (lo && lo[0] && hi_s && hi_s[0]) {
                if (style[0]) sb_append(sb, style);
                sb_append(sb, "_{"); sb_append(sb, lo); sb_append(sb, "}^{"); sb_append(sb, hi_s); sb_append_char(sb, '}');
            } else if (lo && lo[0]) {
                if (style[0]) sb_append(sb, style);
                sb_append(sb, "_{"); sb_append(sb, lo); sb_append_char(sb, '}');
            } else if (hi_s && hi_s[0]) {
                if (style[0]) sb_append(sb, style);
                sb_append(sb, "^{"); sb_append(sb, hi_s); sb_append_char(sb, '}');
            }
            /* trim integrand */
            if (integ) {
                char *trimmed = integ;
                while (*trimmed == ' ') trimmed++;
                if (*trimmed) { sb_append_char(sb, ' '); sb_append(sb, trimmed); }
                else sb_append_char(sb, ' ');
            }
            free(lo); free(hi_s); free(integ); return;
        }

        /* EQNEDT32: slot[0] empty, slot[1] has integrand+limits */
        int s0_empty = (node->u.tmpl.slots.count > 0) ? is_slot_empty(node->u.tmpl.slots.items[0]) : 0;
        if (s0_empty && node->u.tmpl.slots.count >= 2 && !is_slot_empty(node->u.tmpl.slots.items[1])) {
            char *integ = NULL, *lo = NULL, *hi_s = NULL;
            /* Temporarily clear g_parent_bigop during parse_eq_bigop_slot.
             * parse_eq_bigop_slot handles limit extraction directly from slot[1]
             * children (SIZE_SUB → LINE(lo) → LINE(hi) → SIZE_SYM pattern).
             * If g_parent_bigop is set, Pass 2 inside convert_node(integrand)
             * would incorrectly promote display data from nested BigOps to
             * the WRONG parent (e.g. inner sum limits → outer sum template). */
            MtefNode *saved_parent = g_parent_bigop;
            g_parent_bigop = NULL;
            parse_eq_bigop_slot(node->u.tmpl.slots.items[1], prod_ver, &integ, &lo, &hi_s);
            g_parent_bigop = saved_parent;

            /* Check if nested processing promoted display data to us */
            if (!lo && node->u.tmpl.display_lower) {
                StringBuilder ts; sb_init(&ts);
                convert_node(node->u.tmpl.display_lower, prod_ver, &ts);
                lo = sb_detach(&ts);
            }
            if (!hi_s && node->u.tmpl.display_upper) {
                StringBuilder ts; sb_init(&ts);
                convert_node(node->u.tmpl.display_upper, prod_ver, &ts);
                hi_s = sb_detach(&ts);
            }
            /* Pop from BigOp display data stack (LIFO: innermost tmSUM pops last-pushed = innermost limits) */
            if (!lo && !hi_s && g_bigop_disp_stack_n > 0) {
                g_bigop_disp_stack_n--;
                if (g_bigop_disp_stack[g_bigop_disp_stack_n].lo) {
                    StringBuilder ts; sb_init(&ts);
                    convert_node(g_bigop_disp_stack[g_bigop_disp_stack_n].lo, prod_ver, &ts);
                    lo = sb_detach(&ts);
                }
                if (g_bigop_disp_stack[g_bigop_disp_stack_n].hi) {
                    StringBuilder ts; sb_init(&ts);
                    convert_node(g_bigop_disp_stack[g_bigop_disp_stack_n].hi, prod_ver, &ts);
                    hi_s = sb_detach(&ts);
                }
            }

            sb_append(sb, op);
            if (lo && lo[0] && hi_s && hi_s[0]) {
                if (style[0]) sb_append(sb, style);
                sb_append(sb, "_{"); sb_append(sb, lo); sb_append(sb, "}^{"); sb_append(sb, hi_s); sb_append_char(sb, '}');
            } else if (lo && lo[0]) {
                if (style[0]) sb_append(sb, style);
                sb_append(sb, "_{"); sb_append(sb, lo); sb_append_char(sb, '}');
            } else if (hi_s && hi_s[0]) {
                if (style[0]) sb_append(sb, style);
                sb_append(sb, "^{"); sb_append(sb, hi_s); sb_append_char(sb, '}');
            }
            if (integ) { char *t = integ; while(*t==' ')t++; if(*t){sb_append_char(sb,' ');sb_append(sb,t);}else sb_append_char(sb,' '); }
            else sb_append_char(sb, ' ');

            /* extra slots */
            for (int i = 2; i < node->u.tmpl.slots.count; i++) {
                if (is_display_data_slot(node->u.tmpl.slots.items[i])) continue;
                char *ex = slot_str(node, i, prod_ver);
                if (ex) { char *t = ex; while(*t==' ')t++; if(*t) sb_append(sb,t); free(ex); }
            }
            free(integ); free(lo); free(hi_s); return;
        }

        /* Standard: slots [integrand, lower?, upper?] */
        /* Set parent BigOp so nested convert_line can promote display data */
        MtefNode *saved_parent = g_parent_bigop;
        g_parent_bigop = node;
        char *integ = slot_str(node, 0, prod_ver);
        g_parent_bigop = saved_parent;

        /* Pop from BigOp display data stack if nested processing didn't
         * promote limits.  This handles nested BigOps (∑∑∑) where the
         * outer Pass 2 collected all display data and pushed inner pairs
         * to the stack for inner BigOps to consume here. */
        if (!node->u.tmpl.display_lower && !node->u.tmpl.display_upper &&
            g_bigop_disp_stack_n > 0) {
            g_bigop_disp_stack_n--;
            node->u.tmpl.display_lower = g_bigop_disp_stack[g_bigop_disp_stack_n].lo;
            node->u.tmpl.display_upper = g_bigop_disp_stack[g_bigop_disp_stack_n].hi;
        }

        /* Check if nested convert_line promoted display data to us */
        if (node->u.tmpl.display_lower || node->u.tmpl.display_upper) {
            char *lo2 = NULL, *hi2 = NULL;
            if (node->u.tmpl.display_lower) {
                StringBuilder ts; sb_init(&ts);
                convert_node(node->u.tmpl.display_lower, prod_ver, &ts);
                lo2 = sb_detach(&ts);
            }
            if (node->u.tmpl.display_upper) {
                StringBuilder ts; sb_init(&ts);
                convert_node(node->u.tmpl.display_upper, prod_ver, &ts);
                hi2 = sb_detach(&ts);
            }
            sb_append(sb, op);
            if (lo2 && lo2[0] && hi2 && hi2[0]) {
                if (style[0]) sb_append(sb, style);
                sb_append(sb, "_{"); sb_append(sb, lo2); sb_append(sb, "}^{"); sb_append(sb, hi2); sb_append_char(sb, '}');
            } else if (lo2 && lo2[0]) {
                if (style[0]) sb_append(sb, style);
                sb_append(sb, "_{"); sb_append(sb, lo2); sb_append_char(sb, '}');
            } else if (hi2 && hi2[0]) {
                if (style[0]) sb_append(sb, style);
                sb_append(sb, "^{"); sb_append(sb, hi2); sb_append_char(sb, '}');
            }
            if (integ) { char *t = integ; while(*t==' ')t++; if(*t){sb_append_char(sb,' ');sb_append(sb,t);}else sb_append_char(sb,' '); }
            else sb_append_char(sb, ' ');
            free(lo2); free(hi2); free(integ); return;
        }

        int si = 1;
        char *lo = NULL, *hi_s = NULL;
        if (has_lower) { lo = slot_str(node, si, prod_ver); si++; }
        if (has_upper) { hi_s = slot_str(node, si, prod_ver); si++; }

        sb_append(sb, op);
        if (has_lower && has_upper) {
            if (style[0]) sb_append(sb, style);
            sb_append(sb, "_{"); sb_append(sb, lo); sb_append(sb, "}^{"); sb_append(sb, hi_s); sb_append_char(sb, '}');
        } else if (has_lower) {
            if (style[0]) sb_append(sb, style);
            sb_append(sb, "_{"); sb_append(sb, lo); sb_append_char(sb, '}');
        } else if (has_upper) {
            if (style[0]) sb_append(sb, style);
            sb_append(sb, "^{"); sb_append(sb, hi_s); sb_append_char(sb, '}');
        }
        if (integ) { char *t = integ; while(*t==' ')t++; if(*t){sb_append_char(sb,' ');sb_append(sb,t);}else sb_append_char(sb,' '); }
        else sb_append_char(sb, ' ');
        free(integ); free(lo); free(hi_s);
    }
}

/* Known LaTeX function names for TF_FUNCTION grouping */
static int is_known_latex_func(const char *name) {
    static const char *funcs[] = {
        "arccos", "arcsin", "arctan", "arg",
        "cos", "cosh", "cot", "coth", "csc",
        "deg", "det", "dim", "exp", "gcd",
        "hom", "inf", "ker", "lg", "lim",
        "liminf", "limsup", "ln", "log", "max",
        "min", "sec", "sin", "sinh", "sup",
        "tan", "tanh", NULL
    };
    for (int i = 0; funcs[i]; i++) {
        if (strcmp(name, funcs[i]) == 0) return 1;
    }
    return 0;
}

/* Known operator names (not standard LaTeX commands, use \operatorname{}) */
static int is_known_op_name(const char *name) {
    static const char *ops[] = {
        "rot", "curl", "div", "grad",
        "sgn", "tr", "diag", "rank",
        "supp", "adj", "ord", "card",
        "Re", "Im",
        NULL
    };
    for (int i = 0; ops[i]; i++) {
        if (strcmp(name, ops[i]) == 0) return 1;
    }
    return 0;
}

/* Emit a single operator/function name as LaTeX */
static void emit_func_name(const char *name, StringBuilder *sb) {
    if (is_known_latex_func(name)) {
        sb_append_char(sb, '\\');
        sb_append(sb, name);
    } else {
        sb_append(sb, "\\operatorname{");
        sb_append(sb, name);
        sb_append_char(sb, '}');
    }
}

/* Flush accumulated TF_VARIABLE chars: emit known functions as \det, \sin, etc.
 * ONLY converts when the ENTIRE buffer matches a known function name exactly.
 * This is conservative: "det" → \det, but "sink" → "sink" (not \sin + k).
 * Unknown sequences are emitted as plain italic chars. */
static void flush_var_func_buf(char *buf, int *len, StringBuilder *sb) {
    if (*len == 0) return;
    buf[*len] = '\0';
    /* Exact match only: entire buffer must be a known function/operator */
    if (is_known_latex_func(buf) || is_known_op_name(buf)) {
        emit_func_name(buf, sb);
        /* Add space after function name to prevent \det merging with
         * following letter (e.g., \det Z not \detZ) */
        sb_append_char(sb, ' ');
    } else {
        /* Not a function: emit as regular italic chars */
        for (int i = 0; i < *len; i++)
            sb_append_char(sb, buf[i]);
    }
    *len = 0;
}

/* Flush accumulated TF_FUNCTION chars as \sin, \cos, \operatorname{rot}, etc.
 * Splits concatenated operators: "rotrot" → \operatorname{rot}\operatorname{rot}
 * Uses longest-prefix matching to handle overlaps (e.g. "cosh" vs "cos"+"h"). */
static void flush_func_buf(char *buf, int *len, StringBuilder *sb) {
    if (*len == 0) return;
    buf[*len] = '\0';
    /* Single known name? Emit directly. */
    if (is_known_latex_func(buf) || is_known_op_name(buf)) {
        emit_func_name(buf, sb);
        *len = 0;
        return;
    }
    /* Try splitting via longest-prefix matching */
    int pos = 0, total = *len;
    while (pos < total) {
        int best = 0;
        for (int plen = total - pos; plen >= 2; plen--) {
            char saved = buf[pos + plen];
            buf[pos + plen] = '\0';
            int found = is_known_latex_func(buf + pos) || is_known_op_name(buf + pos);
            buf[pos + plen] = saved;
            if (found) { best = plen; break; }
        }
        if (best > 0) {
            char saved = buf[pos + best];
            buf[pos + best] = '\0';
            emit_func_name(buf + pos, sb);
            buf[pos + best] = saved;
            pos += best;
        } else {
            /* No known prefix — emit entire remainder */
            sb_append(sb, "\\operatorname{");
            sb_append(sb, buf + pos);
            sb_append_char(sb, '}');
            break;
        }
    }
    *len = 0;
}

/* ============================================================
 * convert_line — 3-pass line conversion
 * ============================================================ */

static void convert_line(MtefNode *node, int prod_ver, Eq2TexArena *arena, StringBuilder *sb)
{
    if (!node || node->type != NODE_LINE) return;
    if (node->u.line.is_null && node->u.line.children.count == 0) return;

    int my_depth = g_convert_line_depth++;

    int n = node->u.line.children.count;
    if (n == 0) { g_convert_line_depth--; return; }

    /* skip_set: bit array */
    uint8_t *skip = NULL;
    int skip_bytes = (n + 7) / 8;
    skip = (uint8_t *)calloc(skip_bytes, 1);
    if (!skip) return;

#define SKIP_SET(i) (skip[(i)/8] |= (1 << ((i)%8)))
#define IS_SKIPPED(i) (skip[(i)/8] & (1 << ((i)%8)))

    MtefNode **ch = node->u.line.children.items;

    /* === Pass 0a: Display fraction detection ===
     * Pattern: TMPL(tmFRACT, slot[0]=empty) → numer_content... → SIZE_SUB block → denom_content... → SIZE_SUB block
     * Build synthetic LINE nodes for numerator and denominator, merge into fraction slots. */
    if (arena) {
        for (int idx = 0; idx < n; idx++) {
            if (IS_SKIPPED(idx)) continue;
            if (ch[idx]->type != NODE_TMPL) continue;
            if (ch[idx]->u.tmpl.selector != tmFRACT) continue;
            if (ch[idx]->u.tmpl.slots.count < 2) continue;
            if (!is_slot_empty(ch[idx]->u.tmpl.slots.items[0])) continue;

            /* Standard fraction: slot[1] has 2 sub-LINEs (numerator + denominator).
             * Display fraction: slot[1] has only 1 sub-LINE (partial numerator content). */
            {
                MtefNode *s1 = ch[idx]->u.tmpl.slots.items[1];
                if (s1 && s1->type == NODE_LINE && !s1->u.line.is_null) {
                    int sub_line_count = 0;
                    for (int j = 0; j < s1->u.line.children.count; j++) {
                        MtefNode *c = s1->u.line.children.items[j];
                        if (c->type == NODE_LINE && !c->u.line.is_null) sub_line_count++;
                    }
                    if (sub_line_count >= 2) continue; /* standard fraction, not display */
                }
            }

            /* Find first display data block: SIZE_SUB → ... → SIZE_SYM → CHAR → ... */
            int disp1_start = -1;
            for (int j = idx + 1; j < n; j++) {
                if (IS_SKIPPED(j)) continue;
                /* Non-null LINE before SIZE_SUB → remaining SIZE records belong to
                 * BigOps (integral display data), not this fraction. Stop scanning. */
                if (ch[j]->type == NODE_LINE && !ch[j]->u.line.is_null &&
                    ch[j]->u.line.children.count > 0) break;
                if (ch[j]->type == NODE_SIZE && ch[j]->u.size.size_type == SIZETYPE_SUB) {
                    /* Verify it's followed by SIZE_SYM (actual display data, not just subscript) */
                    int has_sym = 0;
                    for (int k = j + 1; k < n; k++) {
                        if (ch[k]->type == NODE_SIZE && ch[k]->u.size.size_type == SIZETYPE_SYM) { has_sym = 1; break; }
                        if (ch[k]->type != NODE_SIZE && ch[k]->type != NODE_LINE) break;
                    }
                    if (has_sym) { disp1_start = j; break; }
                }
            }
            if (disp1_start < 0) {
                /* Fallback: LINE-based split for display fractions without SIZE markers.
                 * Pattern: FRACT(slot[0]=empty) → CHAR continuation → LINE(denominator)
                 * This occurs when EQNEDT32 stores the fraction across nesting levels. */
                int denom_line_idx = -1;
                for (int j = idx + 1; j < n; j++) {
                    if (IS_SKIPPED(j)) continue;
                    if (ch[j]->type == NODE_LINE && !ch[j]->u.line.is_null &&
                        ch[j]->u.line.children.count > 0) {
                        denom_line_idx = j;
                        break;
                    }
                }
                if (denom_line_idx >= 0) {
                    /* Build numerator: slot[1] children + continuation CHARs/TMPLs before denom LINE */
                    MtefNode *slot1 = ch[idx]->u.tmpl.slots.items[1];
                    int s1c = (slot1 && slot1->type == NODE_LINE && !slot1->u.line.is_null)
                              ? slot1->u.line.children.count : 0;
                    int nc = 0;
                    for (int j = idx + 1; j < denom_line_idx; j++) {
                        if (IS_SKIPPED(j)) continue;
                        nc++;
                    }
                    int total_n = s1c + nc;
                    if (total_n > 0) {
                        MtefNode *numer_line = (MtefNode *)arena_alloc(arena, sizeof(MtefNode));
                        memset(numer_line, 0, sizeof(MtefNode));
                        numer_line->type = NODE_LINE;
                        MtefNode **items = (MtefNode **)arena_alloc(arena, total_n * sizeof(MtefNode*));
                        int ni = 0;
                        for (int j = 0; j < s1c; j++)
                            items[ni++] = slot1->u.line.children.items[j];
                        for (int j = idx + 1; j < denom_line_idx; j++) {
                            if (IS_SKIPPED(j)) continue;
                            items[ni++] = ch[j];
                            SKIP_SET(j);
                        }
                        numer_line->u.line.children.items = items;
                        numer_line->u.line.children.count = ni;
                        ch[idx]->u.tmpl.slots.items[0] = numer_line;
                    }

                    /* Build denominator: denom LINE + continuation (non-LINE, non-SIZE) after it */
                    int denom_end = denom_line_idx;
                    for (int j = denom_line_idx + 1; j < n; j++) {
                        if (IS_SKIPPED(j)) continue;
                        if (ch[j]->type == NODE_LINE) break;
                        if (ch[j]->type == NODE_SIZE) break;
                        denom_end = j;
                    }
                    int dc = 0;
                    for (int j = denom_line_idx; j <= denom_end; j++) {
                        if (IS_SKIPPED(j)) continue;
                        dc++;
                    }
                    if (dc > 0) {
                        MtefNode *denom_line = (MtefNode *)arena_alloc(arena, sizeof(MtefNode));
                        memset(denom_line, 0, sizeof(MtefNode));
                        denom_line->type = NODE_LINE;
                        MtefNode **items = (MtefNode **)arena_alloc(arena, dc * sizeof(MtefNode*));
                        int ni = 0;
                        for (int j = denom_line_idx; j <= denom_end; j++) {
                            if (IS_SKIPPED(j)) continue;
                            items[ni++] = ch[j];
                            SKIP_SET(j);
                        }
                        denom_line->u.line.children.items = items;
                        denom_line->u.line.children.count = ni;
                        ch[idx]->u.tmpl.slots.items[1] = denom_line;
                    }

                    /* Skip remaining BigOp display data (integral limits/symbols) */
                    for (int j = denom_end + 1; j < n; j++) {
                        if (IS_SKIPPED(j)) continue;
                        if (ch[j]->type == NODE_SIZE ||
                            ch[j]->type == NODE_LINE ||
                            (ch[j]->type == NODE_CHAR &&
                             ch[j]->u.ch.char_code == 0x222B)) /* ∫ */
                            SKIP_SET(j);
                    }
                }
                continue;
            }

            /* Find end of display data block 1: scan to SIZE_FULL or last display element */
            int disp1_end = disp1_start;
            {
                int sym_seen = 0;
                for (int j = disp1_start; j < n; j++) {
                    disp1_end = j;
                    if (ch[j]->type == NODE_SIZE) {
                        if (ch[j]->u.size.size_type == SIZETYPE_SYM) sym_seen = 1;
                        else if (ch[j]->u.size.size_type == SIZETYPE_FULL && sym_seen) break;
                    } else if (ch[j]->type == NODE_CHAR && sym_seen) {
                        /* display CHAR after SIZE_SYM — include but check if next is SIZE_FULL */
                    } else if (ch[j]->type == NODE_LINE && ch[j]->u.line.is_null) {
                        /* null LINE — skip */
                    } else if (j > disp1_start) break;
                }
            }

            /* Find second display data block */
            int disp2_start = -1, disp2_end = -1;
            for (int j = disp1_end + 1; j < n; j++) {
                if (IS_SKIPPED(j)) continue;
                if (ch[j]->type == NODE_SIZE && ch[j]->u.size.size_type == SIZETYPE_SUB) {
                    int has_sym = 0;
                    for (int k = j + 1; k < n; k++) {
                        if (ch[k]->type == NODE_SIZE && ch[k]->u.size.size_type == SIZETYPE_SYM) { has_sym = 1; break; }
                        if (ch[k]->type != NODE_SIZE && ch[k]->type != NODE_LINE) break;
                    }
                    if (has_sym) {
                        disp2_start = j;
                        int sym_seen = 0;
                        for (int k = j; k < n; k++) {
                            disp2_end = k;
                            if (ch[k]->type == NODE_SIZE) {
                                if (ch[k]->u.size.size_type == SIZETYPE_SYM) sym_seen = 1;
                                else if (ch[k]->u.size.size_type == SIZETYPE_FULL && sym_seen) break;
                            } else if (ch[k]->type == NODE_CHAR && sym_seen) { /* display char */ }
                            else if (ch[k]->type == NODE_LINE && ch[k]->u.line.is_null) { /* null line */ }
                            else if (k > j) break;
                        }
                    }
                    break;
                }
            }

            /* Build numerator LINE: slot[1] children + continuation from idx+1 to disp1_start-1 */
            MtefNode *slot1 = ch[idx]->u.tmpl.slots.items[1];
            int slot1_count = (slot1 && slot1->type == NODE_LINE && !slot1->u.line.is_null)
                              ? slot1->u.line.children.count : 0;
            int numer_cont = 0;
            for (int j = idx + 1; j < disp1_start; j++) {
                if (IS_SKIPPED(j)) continue;
                numer_cont++;
            }
            int total_numer = slot1_count + numer_cont;
            if (total_numer > 0) {
                MtefNode *numer_line = (MtefNode *)arena_alloc(arena, sizeof(MtefNode));
                memset(numer_line, 0, sizeof(MtefNode));
                numer_line->type = NODE_LINE;
                MtefNode **items = (MtefNode **)arena_alloc(arena, total_numer * sizeof(MtefNode*));
                int ni = 0;
                for (int j = 0; j < slot1_count; j++)
                    items[ni++] = slot1->u.line.children.items[j];
                for (int j = idx + 1; j < disp1_start; j++) {
                    if (IS_SKIPPED(j)) continue;
                    items[ni++] = ch[j];
                    SKIP_SET(j);
                }
                numer_line->u.line.children.items = items;
                numer_line->u.line.children.count = ni;
                ch[idx]->u.tmpl.slots.items[0] = numer_line;
            }

            /* Build denominator LINE: from disp1_end+1 to disp2_start-1 (or n) */
            int denom_limit = (disp2_start >= 0) ? disp2_start : n;
            int denom_count = 0;
            for (int j = disp1_end + 1; j < denom_limit; j++) {
                if (IS_SKIPPED(j)) continue;
                denom_count++;
            }
            if (denom_count > 0) {
                MtefNode *denom_line = (MtefNode *)arena_alloc(arena, sizeof(MtefNode));
                memset(denom_line, 0, sizeof(MtefNode));
                denom_line->type = NODE_LINE;
                MtefNode **items = (MtefNode **)arena_alloc(arena, denom_count * sizeof(MtefNode*));
                int ni = 0;
                for (int j = disp1_end + 1; j < denom_limit; j++) {
                    if (IS_SKIPPED(j)) continue;
                    items[ni++] = ch[j];
                    SKIP_SET(j);
                }
                denom_line->u.line.children.items = items;
                denom_line->u.line.children.count = ni;
                ch[idx]->u.tmpl.slots.items[1] = denom_line;
            }

            /* Skip display data blocks */
            for (int j = disp1_start; j <= disp1_end; j++) SKIP_SET(j);
            if (disp2_start >= 0 && disp2_end >= 0) {
                for (int j = disp2_start; j <= disp2_end; j++) SKIP_SET(j);
            }
        }
    }

    /* === Pass 0: SIZE(sub) -> CHAR(display) -> SIZE(sym) skip === */
    for (int idx = 0; idx < n; idx++) {
        if (ch[idx]->type != NODE_SIZE || ch[idx]->u.size.size_type != SIZETYPE_SUB) continue;
        /* find next non-skipped CHAR (skip NULL_LINEs and SIZE nodes) */
        MtefNode *next_real = NULL;
        int next_j = -1;
        for (int j = idx + 1; j < n; j++) {
            if (IS_SKIPPED(j)) continue;
            MtefNode *cand = ch[j];
            if (cand->type == NODE_LINE && cand->u.line.is_null) continue;
            if (cand->type == NODE_SIZE) continue;
            next_real = cand; next_j = j; break;
        }
        if (!next_real || next_real->type != NODE_CHAR) continue;
        uint16_t code = next_real->u.ch.char_code;
        int is_display = set_contains(BIGOP_DISPLAY_CHARS, BIGOP_DISPLAY_N, code);
        if (!is_display && code > 0xFF) {
            uint8_t hi = (code >> 8) & 0xFF;
            is_display = (hi >= 0x80 || hi < 0x20);
        }
        if (!is_display) continue;
        /* Check for SIZE(sym) after */
        int has_sym = 0;
        for (int j = idx + 1; j < n; j++) {
            MtefNode *c = ch[j];
            if (c->type == NODE_SIZE && c->u.size.size_type == SIZETYPE_SYM) { has_sym = 1; break; }
            if (c->type != NODE_SIZE && c->type != NODE_CHAR) break;
        }
        if (!has_sym) continue;
        /* Skip SIZE sub through SIZE sym + trailing display CHAR (including NULL_LINEs) */
        {
            int past_sym = 0;
            for (int j = idx; j < n; j++) {
                if (IS_SKIPPED(j)) continue;
                MtefNode *c = ch[j];
                if (c->type == NODE_SIZE) {
                    SKIP_SET(j);
                    if (c->u.size.size_type == SIZETYPE_SYM) past_sym = 1;
                } else if (c->type == NODE_CHAR) {
                    SKIP_SET(j);
                    if (past_sym) break; /* display CHAR after SIZE_SYM: skip and stop */
                } else if (c->type == NODE_LINE && c->u.line.is_null) {
                    SKIP_SET(j);
                } else break;
            }
        }
    }

    /* === Pass 1: Fence empty-slot merging === */
#ifdef EQ2TEX_DEBUG
    fprintf(stderr, "[P1] convert_line n=%d\n", n);
#endif
    for (int idx = 0; idx < n; idx++) {
        if (IS_SKIPPED(idx)) continue;

        MtefNode *fence_tmpl = NULL;
        int nested = 0;

#ifdef EQ2TEX_DEBUG
        fprintf(stderr, "[P1] idx=%d type=%d\n", idx, ch[idx]->type);
#endif
        if (ch[idx]->type == NODE_TMPL &&
            (is_fence_selector(ch[idx]->u.tmpl.selector) ||
             is_decoration_selector(ch[idx]->u.tmpl.selector))) {
            fence_tmpl = ch[idx];
#ifdef EQ2TEX_DEBUG
            fprintf(stderr, "[P1] found fence sel=%d\n", ch[idx]->u.tmpl.selector);
#endif
        } else if (ch[idx]->type == NODE_LINE && !ch[idx]->u.line.is_null) {
            /* check sub-LINE for fence */
            for (int si = 0; si < ch[idx]->u.line.children.count; si++) {
                MtefNode *sub = ch[idx]->u.line.children.items[si];
                if (sub->type == NODE_TMPL &&
                    (is_fence_selector(sub->u.tmpl.selector) ||
                     is_decoration_selector(sub->u.tmpl.selector))) {
                    MtefNode *s0 = (sub->u.tmpl.slots.count > 0) ? sub->u.tmpl.slots.items[0] : NULL;
                    if (s0 && (s0->u.line.is_null || s0->u.line.children.count == 0)) {
                        fence_tmpl = sub;
                        nested = 1;
                    }
                    break;
                }
            }
        }

        if (!fence_tmpl) continue;
        int s0_empty = (fence_tmpl->u.tmpl.slots.count > 0) ?
            is_slot_empty(fence_tmpl->u.tmpl.slots.items[0]) : 1;
        if (!s0_empty) continue;

        if (nested) {
            /* If the body LINE already contains fence display chars (from
             * parser fence slot[1] flattening), the fence is self-contained.
             * Don't absorb from the parent level — those siblings belong to
             * a parent BigOp (e.g. \sum\limits_{n} \left( ... \right)). */
            {
                int has_inner_fence_disp = 0;
                for (int ci = 0; ci < ch[idx]->u.line.children.count; ci++) {
                    if (is_fence_display_char(ch[idx]->u.line.children.items[ci], prod_ver)) {
                        has_inner_fence_disp = 1; break;
                    }
                }
                if (has_inner_fence_disp) continue;
            }

            /* Nested fence in sub-LINE: fence content is handled by recursive convert_line.
             * EQNEDT32 native: display data for BigOps inside the fence (limit LINEs,
             * NULL_LINEs, SIZE, bigop display chars, fence display chars, trailing content
             * like "dΩ³") is promoted to the current level. Absorb this data into
             * ch[idx] so recursive passes can process it correctly.
             *
             * Fence display chars are ABSORBED (not just skipped) so the recursive
             * direct-fence scan at depth+1 can find them and set found_fence_display=1,
             * preventing the first-child-fence-extension from over-absorbing.
             *
             * Stop at any LINE that follows a fence display char (= outer BigOp's own
             * display data, e.g. two NULL_LINEs before SIZE+CHAR(∫) for tmSINT). */
            int found_fence_disp2 = 0;
            for (int j = idx + 1; j < n; j++) {
                if (IS_SKIPPED(j)) continue;
                MtefNode *c = ch[j];
                if (is_fence_display_char(c, prod_ver)) {
                    /* Absorb fence display char into sub-LINE so depth+1 can see it */
                    if (arena) { nl_push(arena, &ch[idx]->u.line.children, c); }
                    SKIP_SET(j);
                    found_fence_disp2 = 1;
                } else if (found_fence_disp2 && c->type == NODE_LINE) {
                    /* LINE after fence display chars = outer BigOp's display marker: stop */
                    break;
                } else if (c->type == NODE_LINE && c->u.line.is_null) {
                    /* NULL_LINE before fence display: end-of-limits marker for nested BigOp.
                     * Absorb into sub-LINE so Pass 2b at depth+1 can find it. */
                    if (arena) { nl_push(arena, &ch[idx]->u.line.children, c); SKIP_SET(j); }
                } else if (!found_fence_disp2 && c->type == NODE_LINE && !c->u.line.is_null) {
                    /* Non-null LINE before fence display: limit LINE for nested BigOp.
                     * Absorb only if a NULL_LINE follows (confirms it's display data). */
                    int has_null_ahead = 0;
                    for (int k = j + 1; k < n && k < j + 8; k++) {
                        if (IS_SKIPPED(k)) continue;
                        if (ch[k]->type == NODE_LINE && ch[k]->u.line.is_null) { has_null_ahead = 1; break; }
                        if (ch[k]->type == NODE_LINE || ch[k]->type == NODE_SIZE) continue;
                        break;
                    }
                    if (has_null_ahead && arena) {
                        nl_push(arena, &ch[idx]->u.line.children, c);
                        SKIP_SET(j);
                    } else { break; }
                } else if (c->type == NODE_SIZE) {
                    SKIP_SET(j);
                } else if (is_bigop_display_char(c, prod_ver)) {
                    SKIP_SET(j);
                } else if (found_fence_disp2 && arena && c->type == NODE_CHAR &&
                           !is_fence_display_char(c, prod_ver) &&
                           !is_bigop_display_char(c, prod_ver)) {
                    /* Trailing content CHAR after fence display: absorb into sub-LINE */
                    nl_push(arena, &ch[idx]->u.line.children, c);
                    SKIP_SET(j);
                } else if (found_fence_disp2 && arena && c->type == NODE_TMPL &&
                           !is_display_tmpl_selector(c->u.tmpl.selector)) {
                    /* Trailing non-BigOp TMPL after fence display (e.g. tmSCRIPT^3) */
                    nl_push(arena, &ch[idx]->u.line.children, c);
                    SKIP_SET(j);
                } else {
                    break;
                }
            }
        } else {
            /* Direct fence template: find next content LINE or PILE and merge */
            int content_idx = -1;
            for (int j = idx + 1; j < n; j++) {
                if (IS_SKIPPED(j)) continue;
                MtefNode *c = ch[j];
#ifdef EQ2TEX_DEBUG
                fprintf(stderr, "[P1] fence merge search j=%d type=%d\n", j, c->type);
#endif
                if (c->type == NODE_LINE && !c->u.line.is_null && c->u.line.children.count > 0) {
                    content_idx = j; break;
                }
                if (c->type == NODE_PILE) {
                    content_idx = j; break;
                }
                if (c->type == NODE_SIZE) { SKIP_SET(j); continue; }
                if (c->type == NODE_LINE && c->u.line.is_null) { SKIP_SET(j); continue; }
                break;
            }

#ifdef EQ2TEX_DEBUG
            fprintf(stderr, "[P1] fence merge content_idx=%d\n", content_idx);
#endif
            if (content_idx >= 0) {
                /* Merge content into slot[0] */
                if (fence_tmpl->u.tmpl.slots.count > 0)
                    fence_tmpl->u.tmpl.slots.items[0] = ch[content_idx];
                SKIP_SET(content_idx);

                /* Skip display bracket chars after content.
                 * EQNEDT32 native: BigOp display data (LINEs, NULL_LINEs, non-display TMPLs)
                 * may appear between fence content and fence display chars.
                 * Absorb them into the content LINE so nested passes can process them. */
                int found_fence_display = 0;
                for (int j = content_idx + 1; j < n; j++) {
                    if (IS_SKIPPED(j)) continue;
                    MtefNode *c = ch[j];
                    if (is_fence_display_char(c, prod_ver)) { SKIP_SET(j); found_fence_display = 1; }
                    else if (c->type == NODE_SIZE) {
                        if (found_fence_display) break;  /* SIZE after fence display = BigOp display data, stop */
                        SKIP_SET(j);
                    }
                    else if (!found_fence_display && arena && c->type == NODE_LINE) {
                        /* Absorb both null and non-null LINEs (limits, display data) */
                        nl_push(arena, &ch[content_idx]->u.line.children, c);
                        SKIP_SET(j);
                    }
                    else if (!found_fence_display && arena && c->type == NODE_TMPL &&
                             !is_display_tmpl_selector(c->u.tmpl.selector)) {
                        /* Absorb non-bigop TMPLs (e.g. tmSCRIPT) into content */
                        nl_push(arena, &ch[content_idx]->u.line.children, c);
                        SKIP_SET(j);
                    }
                    else if (!found_fence_display && is_bigop_display_char(c, prod_ver)) {
                        /* Bare bigop display CHAR (e.g. Sigma for tmSUM): skip it.
                         * The empty NULL_LINE absorbed above serves as the display data
                         * marker for Pass 2b; the bare CHAR itself is redundant here. */
                        SKIP_SET(j);
                    }
                    else { break; }
                }

                /* First-child fence extension: if no display chars found at this level
                 * and fence is the first non-SIZE/null element, merge ALL remaining
                 * children into a synthetic LINE as slot[0].
                 * (EQNEDT32 puts display chars at the parent level)
                 * NOTE: Only for fences, NOT decorations (overbar/underbar/arrows).
                 * Decorations merge only ONE sibling LINE, not everything. */
                if (!found_fence_display && arena
                    && is_fence_selector(fence_tmpl->u.tmpl.selector)) {
                    int is_first = 1;
                    for (int j = 0; j < idx; j++) {
                        if (IS_SKIPPED(j)) continue;
                        MtefNode *c = ch[j];
                        if (c->type == NODE_SIZE || (c->type == NODE_LINE && c->u.line.is_null))
                            continue;
                        is_first = 0; break;
                    }
                    if (is_first) {
                        /* Collect remaining: content_line + all non-skipped after it */
                        int rem_count = 0;
                        for (int j = content_idx; j < n; j++) {
                            if (IS_SKIPPED(j) && j != content_idx) continue;
                            rem_count++;
                        }
                        if (rem_count > 1) {
                            MtefNode *synth = new_node(arena, NODE_LINE);
                            if (synth) {
                                synth->u.line.is_null = 0;
                                for (int j = content_idx; j < n; j++) {
                                    if (IS_SKIPPED(j) && j != content_idx) continue;
                                    nl_push(arena, &synth->u.line.children, ch[j]);
                                    SKIP_SET(j);
                                }
                                fence_tmpl->u.tmpl.slots.items[0] = synth;
                            }
                        }
                    }
                }
            }
        }
    }

    /* === Pass 2: BigOp remote display data detection === */
    /* EQNEDT32 pattern: TMPL(BigOp, slot0=empty) → LINE(content) → SIZE(sub) → LINE(null/limit)* → SIZE(sym) → CHAR(display)* → SIZE(full) */
    for (int idx = 0; idx < n; idx++) {
        if (IS_SKIPPED(idx)) continue;
        if (ch[idx]->type != NODE_SIZE || ch[idx]->u.size.size_type != SIZETYPE_SUB) continue;

        /* Collect display data block: SIZE(sub) → LINE* → SIZE(sym) → CHAR(display)* → SIZE(full) */
        int block[128];
        int block_n = 0;
        MtefNode *limit_lines[16];
        int limit_n = 0;
        int sym_found = 0;

        block[block_n++] = idx;

        for (int j = idx + 1; j < n && block_n < 128; j++) {
            if (IS_SKIPPED(j)) continue;
            MtefNode *c = ch[j];
            if (c->type == NODE_LINE) {
                block[block_n++] = j;
                if (!c->u.line.is_null && c->u.line.children.count > 0 && limit_n < 16)
                    limit_lines[limit_n++] = c;
                /* NULL_LINE may contain SIZE_SYM/SUBSYM + bigop display char
                 * (EQNEDT32 native packed display data marker) */
                if (c->u.line.is_null && c->u.line.children.count > 0) {
                    for (int k = 0; k < c->u.line.children.count; k++) {
                        MtefNode *nc = c->u.line.children.items[k];
                        if (nc->type == NODE_SIZE &&
                            (nc->u.size.size_type == SIZETYPE_SYM ||
                             nc->u.size.size_type == SIZETYPE_SUBSYM))
                            sym_found = 1;
                    }
                }
                continue;
            }
            if (c->type == NODE_SIZE) {
                block[block_n++] = j;
                if (c->u.size.size_type == SIZETYPE_SYM ||
                    c->u.size.size_type == SIZETYPE_SUBSYM)
                    sym_found = 1;
                else if ((c->u.size.size_type == SIZETYPE_FULL) && sym_found)
                    break;
                continue;
            }
            if (c->type == NODE_CHAR && sym_found) {
                if (is_bigop_display_char(c, prod_ver)) {
                    block[block_n++] = j;
                    continue;
                }
            }
            break;
        }

        if (!sym_found) {
            continue;
        }

        /* Find nearest preceding DISPLAY_TMPL_SELECTORS template (reverse scan) */
        MtefNode *target_tmpl = NULL;
        int content_idx = -1;
        for (int j = idx - 1; j >= 0; j--) {
            if (IS_SKIPPED(j)) continue;
            MtefNode *c = ch[j];
            if (c->type == NODE_TMPL && is_display_tmpl_selector(c->u.tmpl.selector) &&
                c->u.tmpl.variation == c->u.tmpl.orig_variation) {
                /* Skip BigOps that already have limits assigned (variation changed) */
                target_tmpl = c;
                break;
            }
            if (content_idx < 0 && c->type == NODE_LINE && !c->u.line.is_null)
                content_idx = j;
        }

        /* No BigOp sibling found: promote display data to parent BigOp
         * using display_lower/display_upper (checked post-slot_str). */
        if (!target_tmpl) {
            if (g_parent_bigop && is_display_tmpl_selector(g_parent_bigop->u.tmpl.selector)) {
                /* Set limits directly on parent via display_lower/display_upper */
                if (limit_n >= 1) g_parent_bigop->u.tmpl.display_lower = limit_lines[0];
                if (limit_n >= 2) g_parent_bigop->u.tmpl.display_upper = limit_lines[1];
                /* Skip all display data block indices */
                for (int bi = 0; bi < block_n; bi++)
                    SKIP_SET(block[bi]);
                continue;
            } else {
                continue;
            }
        }

        /* EQNEDT32 packed format: slot[1] already has integrand+limits.
         * Don't merge external content into slot[0] or overwrite slot[1]. */
        int has_packed_slot1 = (target_tmpl->u.tmpl.slots.count >= 2 &&
                                !is_slot_empty(target_tmpl->u.tmpl.slots.items[1]));

        /* Merge content LINE into slot[0] if empty (and not EQNEDT32 packed) */
        if (!has_packed_slot1 &&
            target_tmpl->u.tmpl.slots.count > 0 &&
            is_slot_empty(target_tmpl->u.tmpl.slots.items[0]) &&
            content_idx >= 0 && !IS_SKIPPED(content_idx)) {
            target_tmpl->u.tmpl.slots.items[0] = ch[content_idx];
            SKIP_SET(content_idx);
        }

        /* Set limit lines.  For integral templates (tmSINT–tmTSINT), use
         * display_lower/display_upper instead of modifying variation, because
         * variation has different semantics for integrals (0=standard, 2=oint,
         * 3=contour etc.).  For BigOps (tmSUM etc.), use slot-based approach. */
        if (limit_n > 0 && !has_packed_slot1) {
            int sel = target_tmpl->u.tmpl.selector;
            if (sel >= tmSINT && sel <= tmTSINT) {
                /* Integral: use display_lower/display_upper */
                if (limit_n >= 1) target_tmpl->u.tmpl.display_lower = limit_lines[0];
                if (limit_n >= 2) target_tmpl->u.tmpl.display_upper = limit_lines[1];
            } else if (arena) {
                /* BigOp: use slot-based approach.
                 * For nested BigOps (∑∑∑), display data from multiple levels
                 * appears in REVERSE order (inner first, outer last).
                 * Use the LAST pair of limit LINEs for the current template. */
                if (limit_n >= 2) {
                    target_tmpl->u.tmpl.variation |= 0x03;
                    while (target_tmpl->u.tmpl.slots.count < 3) {
                        MtefNode *nl = (MtefNode *)arena_alloc(arena, sizeof(MtefNode));
                        if (!nl) break;
                        memset(nl, 0, sizeof(MtefNode));
                        nl->type = NODE_LINE; nl->u.line.is_null = 1;
                        nl_push(arena, &target_tmpl->u.tmpl.slots, nl);
                    }
                    if (target_tmpl->u.tmpl.slots.count >= 3) {
                        target_tmpl->u.tmpl.slots.items[1] = limit_lines[limit_n - 2];
                        target_tmpl->u.tmpl.slots.items[2] = limit_lines[limit_n - 1];
                    }
                    /* Nested BigOps (∑∑∑): remaining limit pairs belong to
                     * inner BigOps.  Push them to g_bigop_disp_stack (LIFO)
                     * so that inner BigOps pop their limits during convert_tmpl.
                     * Pairs are in reverse order (inner-first), push in forward
                     * order so that the shallowest inner BigOp pops last-pushed. */
                    for (int pi = 0; pi + 1 < limit_n - 2; pi += 2) {
                        if (g_bigop_disp_stack_n < BIGOP_DISP_STACK_MAX) {
                            g_bigop_disp_stack[g_bigop_disp_stack_n].lo = limit_lines[pi];
                            g_bigop_disp_stack[g_bigop_disp_stack_n].hi = limit_lines[pi + 1];
                            g_bigop_disp_stack_n++;
                        }
                    }
                } else {
                    target_tmpl->u.tmpl.variation |= 0x01;
                    while (target_tmpl->u.tmpl.slots.count < 2) {
                        MtefNode *nl = (MtefNode *)arena_alloc(arena, sizeof(MtefNode));
                        if (!nl) break;
                        memset(nl, 0, sizeof(MtefNode));
                        nl->type = NODE_LINE; nl->u.line.is_null = 1;
                        nl_push(arena, &target_tmpl->u.tmpl.slots, nl);
                    }
                    if (target_tmpl->u.tmpl.slots.count >= 2)
                        target_tmpl->u.tmpl.slots.items[1] = limit_lines[0];
                }
            }
        }

        /* Skip all display data block indices */
        for (int bi = 0; bi < block_n; bi++)
            SKIP_SET(block[bi]);
    }

    /* === Pass 2c: EQNEDT32 native tmSCRIPT sup-LINE absorption ===
     * When TMPL(tmSCRIPT/tmLSCRIPT, var=2) has slot[0]=empty and slot[1] with only one
     * non-null LINE (the sub), the sup LINE is the immediately following non-null sibling.
     * Absorb it into slot[1] so that convert_tmpl sees both sub and sup in one place. */
    if (arena) {
        for (int idx = 0; idx < n; idx++) {
            if (IS_SKIPPED(idx)) continue;
            MtefNode *c = ch[idx];
            if (c->type != NODE_TMPL) continue;
            if (c->u.tmpl.selector != tmSCRIPT && c->u.tmpl.selector != tmLSCRIPT) continue;
            if (c->u.tmpl.orig_variation != 2) continue;  /* only var=2 (both sub+sup) */
            if (c->u.tmpl.slots.count < 2) continue;
            if (!is_slot_empty(c->u.tmpl.slots.items[0])) continue;
            /* Count non-null LINEs in slot[1] */
            MtefNode *s1 = c->u.tmpl.slots.items[1];
            int nlines1 = 0;
            for (int i = 0; i < s1->u.line.children.count; i++) {
                MtefNode *sc = s1->u.line.children.items[i];
                if (sc->type == NODE_LINE && !sc->u.line.is_null) nlines1++;
            }
            if (nlines1 != 1) continue;  /* already has both sub+sup or neither */
            /* Look for the sup LINE as next non-skipped sibling */
            for (int j = idx + 1; j < n; j++) {
                if (IS_SKIPPED(j)) continue;
                MtefNode *nc = ch[j];
                if (nc->type == NODE_SIZE) continue;  /* skip SIZE records */
                if (nc->type == NODE_LINE && !nc->u.line.is_null) {
                    nl_push(arena, &s1->u.line.children, nc);  /* absorb as sup */
                    SKIP_SET(j);
                }
                break;
            }
        }
    }

    /* === Pass 2b: BigOp display data without SIZE_SUB === */
    /* EQNEDT32 native pattern for tmSUM/tmISUM/tmPROD etc. with limits:
     *   TMPL(BigOp, slot0=empty)
     *   → LINE(content)
     *   → [LINE(lo)]?  [LINE(hi)]?   (non-null limit LINEs, before display-data marker)
     *   → NULL_LINE*(with bigop display chars OR SIZE_SYM inside)
     *   → SIZE(SYM)?  SIZE(FULL)?
     *   → CHAR(bigop_display)*
     * "disp_data_found" is set by any of: SIZE_SYM, NULL_LINE containing bigop display char. */
    for (int idx = 0; idx < n; idx++) {
        if (IS_SKIPPED(idx)) continue;
        if (ch[idx]->type != NODE_TMPL) continue;
        if (!is_display_tmpl_selector(ch[idx]->u.tmpl.selector)) continue;
        if (ch[idx]->u.tmpl.slots.count == 0 ||
            !is_slot_empty(ch[idx]->u.tmpl.slots.items[0])) continue;

        /* Find next non-skipped content LINE after this template.
         * Accept empty LINEs (children.count == 0) as body — consecutive
         * BigOps produce empty body when the second BigOp is consumed
         * by consume_bigop_body's stop-at-next-bigop rule. */
        int content_idx = -1;
        for (int j = idx + 1; j < n; j++) {
            if (IS_SKIPPED(j)) continue;
            if (ch[j]->type == NODE_LINE && !ch[j]->u.line.is_null) {
                content_idx = j; break;
            }
            if (ch[j]->type == NODE_SIZE) continue;
            break;
        }
        if (content_idx < 0) continue;

        /* After content, collect limit LINEs and display-data block */
        MtefNode *limit_lines2[4]; int limit2_n = 0;
        int block2[128]; int block2_n = 0;
        int disp_data_found = 0;
        for (int j = content_idx + 1; j < n && block2_n < 128; j++) {
            if (IS_SKIPPED(j)) continue;
            MtefNode *c = ch[j];
            if (c->type == NODE_LINE && c->u.line.is_null) {
                /* NULL_LINE is always the end-of-limits / display-data marker.
                 * In EQNEDT32 native encoding the NULL_LINE may be empty (the
                 * actual SIZE_SYM + bigop CHAR follow as outer siblings and were
                 * already bigop_char_skip'ed by Pass 1 fence scan). */
                disp_data_found = 1;
                block2[block2_n++] = j; continue;
            }
            /* Non-null LINE: before display data = limit; after = stop */
            if (c->type == NODE_LINE && !c->u.line.is_null) {
                if (!disp_data_found && limit2_n < 2) {
                    limit_lines2[limit2_n++] = c;
                    block2[block2_n++] = j; continue;
                }
                break;
            }
            /* Non-display TMPL (e.g. tmSCRIPT) between content and display data:
             * absorb into content LINE so it attaches to the trailing char of content. */
            if (c->type == NODE_TMPL && !is_display_tmpl_selector(c->u.tmpl.selector)) {
                if (arena) nl_push(arena, &ch[content_idx]->u.line.children, c);
                block2[block2_n++] = j; continue;
            }
            if (c->type == NODE_SIZE) {
                block2[block2_n++] = j;
                if (c->u.size.size_type == SIZETYPE_SYM ||
                    c->u.size.size_type == SIZETYPE_SUBSYM) disp_data_found = 1;
                continue;  /* SIZE_FULL no longer terminates — trailing display chars follow */
            }
            if (c->type == NODE_CHAR && disp_data_found &&
                is_bigop_display_char(c, prod_ver)) {
                block2[block2_n++] = j; continue;
            }
            break;
        }
        if (!disp_data_found) continue;

        /* EQNEDT32 packed format: slot[1] already has integrand+limits.
         * Don't merge external content into slot[0] — it's a separate
         * expression AFTER the bigop (e.g. the second "dx" in ∫f(x)dx dx). */
        int has_packed2 = (ch[idx]->u.tmpl.slots.count >= 2 &&
                           !is_slot_empty(ch[idx]->u.tmpl.slots.items[1]));
        if (!has_packed2) {
            ch[idx]->u.tmpl.slots.items[0] = ch[content_idx];
            SKIP_SET(content_idx);
            /* Attach limit LINEs as display_lower/display_upper */
            if (limit2_n >= 1) ch[idx]->u.tmpl.display_lower = limit_lines2[0];
            if (limit2_n >= 2) ch[idx]->u.tmpl.display_upper = limit_lines2[1];
        }

        /* Skip all display data */
        for (int bi = 0; bi < block2_n; bi++)
            SKIP_SET(block2[bi]);
    }

    /* === Pass 2d: tmSCRIPT-adjacent BigOp display data (nested sum pattern) ===
     * EQNEDT32 native nested BigOps (e.g. ∑∑∑) store display data at the
     * continuation LINE level, AFTER tmSCRIPT (not after BigOp templates).
     * Pattern: ..., tmSCRIPT, LINE(lo), LINE(hi), SIZE_SYM, CHAR(∑), ...
     * Push limits to g_bigop_disp_stack for the parent BigOp to pop. */
    if (g_parent_bigop && is_display_tmpl_selector(g_parent_bigop->u.tmpl.selector)) {
        for (int idx = 0; idx < n; idx++) {
            if (IS_SKIPPED(idx)) continue;
            MtefNode *c = ch[idx];
            /* Look for SIZE_SYM as the anchor */
            if (c->type != NODE_SIZE ||
                (c->u.size.size_type != SIZETYPE_SYM &&
                 c->u.size.size_type != SIZETYPE_SUBSYM)) continue;
            /* Verify next record is BigOp display char */
            int disp_idx = -1;
            for (int j = idx + 1; j < n; j++) {
                if (IS_SKIPPED(j)) continue;
                if (ch[j]->type == NODE_CHAR && is_bigop_display_char(ch[j], prod_ver)) {
                    disp_idx = j; break;
                }
                break;
            }
            if (disp_idx < 0) continue;
            /* Check if there's a tmSCRIPT before the display data (with limit LINEs in between) */
            int has_script_before = 0;
            MtefNode *lo_line = NULL, *hi_line = NULL;
            int first_limit_idx = -1;
            for (int j = idx - 1; j >= 0; j--) {
                if (IS_SKIPPED(j)) continue;
                MtefNode *pc = ch[j];
                if (pc->type == NODE_LINE && !pc->u.line.is_null) {
                    if (!hi_line) { hi_line = pc; first_limit_idx = j; }
                    else if (!lo_line) { lo_line = pc; first_limit_idx = j; }
                    continue;
                }
                if (pc->type == NODE_TMPL &&
                    (pc->u.tmpl.selector == tmSCRIPT || pc->u.tmpl.selector == tmLSCRIPT)) {
                    has_script_before = 1;
                    break;
                }
                break;  /* unexpected node type */
            }
            if (!has_script_before || !lo_line) continue;
            /* Swap lo/hi if reverse-scanned (hi was found first) */
            if (lo_line && hi_line) {
                /* lo_line was found second (= earlier in children), hi_line first (= later) */
                /* They're already in correct order: lo before hi */
            }
            /* Push to BigOp display data stack */
            if (g_bigop_disp_stack_n < BIGOP_DISP_STACK_MAX) {
                g_bigop_disp_stack[g_bigop_disp_stack_n].lo = lo_line;
                g_bigop_disp_stack[g_bigop_disp_stack_n].hi = hi_line;
                g_bigop_disp_stack_n++;
                /* Skip limit LINEs + SIZE_SYM + display char */
                if (first_limit_idx >= 0) {
                    for (int j = first_limit_idx; j <= disp_idx; j++)
                        if (!IS_SKIPPED(j)) SKIP_SET(j);
                }
            }
        }
    }

    /* === Pass 3: Cursor garbage removal (trailing ASCII letters) === */
    /* Only at top-level LINE (depth 0).  Nested LINEs (template slots,
     * integrand content, etc.) never have cursor garbage — their trailing
     * chars like "dx" are intentional. */
    if (my_depth == 0) {
        int has_display_tmpl = 0;
        for (int i = 0; i < n && !has_display_tmpl; i++)
            if (ch[i]->type == NODE_TMPL &&
                is_display_tmpl_selector(ch[i]->u.tmpl.selector))
                has_display_tmpl = 1;
        if (has_display_tmpl) {
            int trailing[64]; /* indices of trailing ASCII chars */
            int trail_n = 0;
            for (int i = n - 1; i >= 0 && trail_n < 64; i--) {
                if (IS_SKIPPED(i)) continue;
                if (ch[i]->type == NODE_CHAR) {
                    uint16_t code = ch[i]->u.ch.char_code;
                    if ((code >= 0x41 && code <= 0x5A) || (code >= 0x61 && code <= 0x7A)) {
                        trailing[trail_n++] = i;
                        continue;
                    }
                }
                break;
            }
            if (trail_n >= 2) {
                /* Check if first trailing char follows a non-ASCII operator */
                int first_idx = trailing[trail_n - 1]; /* leftmost trailing ASCII */
                int prev_idx = -1;
                for (int j = first_idx - 1; j >= 0; j--) {
                    if (!IS_SKIPPED(j)) { prev_idx = j; break; }
                }
                /* Don't remove trailing letters after BigOp templates (dx, dt, etc.) */
                int follows_bigop = 0;
                for (int j = first_idx - 1; j >= 0; j--) {
                    if (IS_SKIPPED(j)) continue;
                    if (ch[j]->type == NODE_TMPL &&
                        is_display_tmpl_selector(ch[j]->u.tmpl.selector)) {
                        follows_bigop = 1; break;
                    }
                    if (ch[j]->type == NODE_CHAR || ch[j]->type == NODE_LINE) break;
                }
                if (!follows_bigop) {
                    int keep_first = 0;
                    if (prev_idx >= 0 && ch[prev_idx]->type == NODE_CHAR &&
                        ch[prev_idx]->u.ch.char_code > 0x7F)
                        keep_first = 1;
                    for (int ti = 0; ti < trail_n - keep_first; ti++)
                        SKIP_SET(trailing[ti]);
                }
            }
        }
    }

    /* === Output phase === */
    int in_text = 0; /* grouping TF_TEXT chars into \text{...} */
    char func_buf[32]; int func_len = 0; /* TF_FUNCTION accumulator */
    /* TF_VARIABLE function detection: accumulate lowercase letters,
     * check if they form a known function name (det, sin, cos, ...).
     * Emit as \det etc. if matched, otherwise as regular italic chars. */
    char var_buf[32]; int var_len = 0; int var_start_idx = -1;
    for (int i = 0; i < n; i++) {
        if (IS_SKIPPED(i)) continue;
        MtefNode *c = ch[i];
        if (c->type == NODE_SIZE || c->type == NODE_FONT) continue;

        /* Text char grouping */
        if (c->type == NODE_CHAR && c->u.ch.typeface == TF_TEXT) {
            flush_func_buf(func_buf, &func_len, sb);
            flush_var_func_buf(var_buf, &var_len, sb);
            if (!in_text) { sb_append(sb, "\\text{"); in_text = 1; }
            uint16_t code = c->u.ch.char_code;
            if (code >= 0x20 && code < 0x7F) {
                char cc = (char)code;
                if (cc == '#' || cc == '$' || cc == '%' || cc == '&' ||
                    cc == '_' || cc == '{' || cc == '}') {
                    sb_append_char(sb, '\\');
                }
                sb_append_char(sb, cc);
            } else if (code >= 0x80) {
                /* Non-ASCII: encode as UTF-8 */
                char utf8[5]; int ul = 0;
                if (code < 0x800) {
                    utf8[ul++] = (char)(0xC0 | (code >> 6));
                    utf8[ul++] = (char)(0x80 | (code & 0x3F));
                } else {
                    utf8[ul++] = (char)(0xE0 | (code >> 12));
                    utf8[ul++] = (char)(0x80 | ((code >> 6) & 0x3F));
                    utf8[ul++] = (char)(0x80 | (code & 0x3F));
                }
                utf8[ul] = '\0';
                sb_append(sb, utf8);
            }
            /* Apply embellishments */
            for (int ei = 0; ei < c->u.ch.embells.count; ei++) {
                MtefNode *e = c->u.ch.embells.items[ei];
                if (e->type == NODE_EMBELL && e->u.embell.embell_type >= 0 &&
                    e->u.embell.embell_type < EMBELL_MAP_N) {
                    /* close text group, apply embell, reopen would be complex; skip for text */
                }
            }
            continue;
        }

        /* Function char grouping (TF_FUNCTION → \sin, \cos, etc.) */
        if (c->type == NODE_CHAR && c->u.ch.typeface == TF_FUNCTION) {
            flush_var_func_buf(var_buf, &var_len, sb);
            if (in_text) { sb_append_char(sb, '}'); in_text = 0; }
            uint16_t code = c->u.ch.char_code;
            if (code == 0x5E) {
                flush_func_buf(func_buf, &func_len, sb);
                sb_append(sb, " \\wedge ");
            } else if (code >= 'a' && code <= 'z' && func_len < 30) {
                /* Only lowercase letters accumulate into function names */
                func_buf[func_len++] = (char)code;
            } else if (code >= 0x20 && code < 0x7F) {
                /* Non-letter function chars (parens, etc.): output directly */
                flush_func_buf(func_buf, &func_len, sb);
                sb_append_char(sb, (char)code);
            }
            continue;
        }

        /* TF_VARIABLE function detection (det, sin, cos, Im, Re, etc.)
         * EQNEDT32 sometimes stores known functions as TF_VARIABLE instead
         * of TF_FUNCTION. Accumulate lowercase letters and check on flush. */
        if (c->type == NODE_CHAR && c->u.ch.typeface == TF_VARIABLE &&
            c->u.ch.embells.count == 0) {
            /* Only accumulate chars WITHOUT embellishments.
             * Embellished chars (e.g. \hat{a}) must go through convert_node
             * to preserve their decorations. */
            flush_func_buf(func_buf, &func_len, sb);
            if (in_text) { sb_append_char(sb, '}'); in_text = 0; }
            uint16_t code = c->u.ch.char_code;
            if (code >= 'A' && code <= 'Z') {
                /* Uppercase after lowercase might start a new word (e.g., "detZ").
                 * Flush lowercase first, then start new accumulation for "Re"/"Im". */
                flush_var_func_buf(var_buf, &var_len, sb);
                var_buf[var_len++] = (char)code;
            } else if (code >= 'a' && code <= 'z' && var_len < 30) {
                var_buf[var_len++] = (char)code;
            } else {
                flush_var_func_buf(var_buf, &var_len, sb);
                convert_node(c, prod_ver, sb);
            }
            continue;
        }

        flush_func_buf(func_buf, &func_len, sb);
        flush_var_func_buf(var_buf, &var_len, sb);
        if (in_text) { sb_append_char(sb, '}'); in_text = 0; }
        convert_node(c, prod_ver, sb);
    }
    flush_func_buf(func_buf, &func_len, sb);
    flush_var_func_buf(var_buf, &var_len, sb);
    if (in_text) sb_append_char(sb, '}');

#undef SKIP_SET
#undef IS_SKIPPED
    free(skip);
    g_convert_line_depth--;
}

/* ============================================================
 * convert_pile — Multi-line display handling
 * ============================================================ */

/* Find display_tmpl in line (start or end) */
static MtefNode *find_bigop_in_line(MtefNode *line)
{
    if (!line || line->type != NODE_LINE) return NULL;
    for (int i = 0; i < line->u.line.children.count; i++) {
        MtefNode *c = line->u.line.children.items[i];
        if (c->type == NODE_TMPL && is_display_tmpl_selector(c->u.tmpl.selector))
            return c;
        if (c->type == NODE_SIZE) continue;
        if (c->type == NODE_LINE && c->u.line.is_null) continue;
        if (c->type == NODE_CHAR || c->type == NODE_TMPL) break;
    }
    for (int i = line->u.line.children.count - 1; i >= 0; i--) {
        MtefNode *c = line->u.line.children.items[i];
        if (c->type == NODE_TMPL && is_display_tmpl_selector(c->u.tmpl.selector))
            return c;
        if (c->type == NODE_SIZE) continue;
        break;
    }
    return NULL;
}

static int find_all_bigops_in_line(MtefNode *line, MtefNode **out, int max)
{
    if (!line || line->type != NODE_LINE) return 0;
    int count = 0;
    for (int i = 0; i < line->u.line.children.count && count < max; i++) {
        MtefNode *c = line->u.line.children.items[i];
        if (c->type == NODE_TMPL && is_display_tmpl_selector(c->u.tmpl.selector))
            out[count++] = c;
        else if (c->type == NODE_LINE)
            count += find_all_bigops_in_line(c, out + count, max - count);
    }
    return count;
}

static int is_display_data_line(MtefNode *line, int prod_ver)
{
    if (!line || line->type != NODE_LINE) return 0;
    int has_disp = 0, content_count = 0;
    for (int i = 0; i < line->u.line.children.count; i++) {
        MtefNode *c = line->u.line.children.items[i];
        if (c->type == NODE_TMPL || c->type == NODE_MATRIX) return 0;
        if (is_bigop_display_char(c, prod_ver)) has_disp = 1;
        /* SIZE_SYM inside a NULL_LINE signals display data even if
         * the actual bigop display char appears outside the NULL_LINE. */
        else if (c->type == NODE_SIZE && c->u.size.size_type == SIZETYPE_SYM) has_disp = 1;
        else if (c->type != NODE_SIZE) content_count++;
    }
    if (has_disp) return 1;
    /* Only treat as display data if line has display chars.
     * Don't consume normal content lines (e.g. "epsilon_r = 1") as display data
     * just because they have few items. */
    return 0;
}

/* is cursor garbage line: all ASCII letters */
static int is_cursor_garbage_line(MtefNode *line)
{
    /* EQNEDT32 cursor artifact: a single alphabetic character in TF_VARIABLE
     * typeface appended at the end of a PILE.  Must be exactly 1 char; lines
     * with 2+ chars are genuine content (e.g. a row in a gathered env). */
    if (!line || line->type != NODE_LINE) return 0;
    int nchars = 0;
    for (int i = 0; i < line->u.line.children.count; i++) {
        MtefNode *c = line->u.line.children.items[i];
        if (c->type == NODE_SIZE) continue;   /* skip SIZE records */
        nchars++;
        if (c->type != NODE_CHAR) return 0;
        uint16_t code = c->u.ch.char_code;
        if (!((code >= 'A' && code <= 'Z') || (code >= 'a' && code <= 'z'))) return 0;
    }
    return nchars == 1;   /* single-char only */
}

/* Extract display limits from display lines */
static void extract_display_limits(MtefNode **lines, int nlines, int prod_ver,
                                   Eq2TexArena *arena, MtefNode **out_lower, MtefNode **out_upper)
{
    *out_lower = *out_upper = NULL;
    /* Collect non-display, non-size children from display lines as limits */
    MtefNode *limits[16];
    int nlimits = 0;

    for (int li = 0; li < nlines; li++) {
        MtefNode *dl = lines[li];
        if (!dl || dl->type != NODE_LINE) continue;
        NodeList group = {NULL, 0, 0};
        for (int i = 0; i < dl->u.line.children.count; i++) {
            MtefNode *c = dl->u.line.children.items[i];
            if (is_bigop_display_char(c, prod_ver)) {
                if (group.count > 0 && nlimits < 16) {
                    MtefNode *ln = new_node(arena, NODE_LINE);
                    if (ln) { ln->u.line.children = group; limits[nlimits++] = ln; }
                    group.count = 0; group.items = NULL; group.cap = 0;
                }
            } else if (c->type != NODE_SIZE) {
                nl_push(arena, &group, c);
            }
        }
        if (group.count > 0 && nlimits < 16) {
            MtefNode *ln = new_node(arena, NODE_LINE);
            if (ln) { ln->u.line.children = group; limits[nlimits++] = ln; }
        }
    }
    if (nlimits >= 1) *out_lower = limits[0];
    if (nlimits >= 2) *out_upper = limits[1];
}

static void convert_pile(MtefNode *node, int prod_ver, Eq2TexArena *arena, StringBuilder *sb)
{
    if (!node || node->type != NODE_PILE || node->u.pile.lines.count == 0) return;

    int nlines = node->u.pile.lines.count;
    MtefNode **src = node->u.pile.lines.items;

    /* Build filtered_lines (max 256 lines) */
    MtefNode *filtered[256];
    int fcount = 0;
    int i = 0;

    while (i < nlines && fcount < 256) {
        MtefNode *line = src[i];
        MtefNode *bigops[16];
        int nbigops = find_all_bigops_in_line(line, bigops, 16);
        MtefNode *start_bigop = find_bigop_in_line(line);

        if (nbigops > 0) {
            filtered[fcount++] = line;
            i++;

            /* consume display data lines */
            MtefNode *disp_lines[32];
            int dcount = 0;
            while (i < nlines && dcount < 32) {
                MtefNode *next = src[i];
                if (find_all_bigops_in_line(next, bigops, 16) > 0) break;
                if (is_display_data_line(next, prod_ver)) {
                    disp_lines[dcount++] = next;
                    i++;
                    continue;
                }
                break;
            }

            /* merge limits */
            if (start_bigop && dcount > 0) {
                MtefNode *lo = NULL, *hi = NULL;
                extract_display_limits(disp_lines, dcount, prod_ver, arena, &lo, &hi);
                if (lo) start_bigop->u.tmpl.display_lower = lo;
                if (hi) start_bigop->u.tmpl.display_upper = hi;
            }
        } else {
            filtered[fcount++] = line;
            i++;
        }
    }

    /* Decoration merging */
    i = 0;
    int merge_safe = 0;
    while (i < fcount - 1) {
        MtefNode *line = filtered[i];
        if (line->type == NODE_LINE && line->u.line.children.count > 0) {
            MtefNode *last_tmpl = NULL;
            for (int j = line->u.line.children.count - 1; j >= 0; j--) {
                MtefNode *c = line->u.line.children.items[j];
                if (c->type == NODE_TMPL && is_decoration_selector(c->u.tmpl.selector)) {
                    last_tmpl = c; break;
                }
                if (c->type == NODE_SIZE) continue;
                break;
            }
            if (last_tmpl) {
                int s0e = (last_tmpl->u.tmpl.slots.count > 0) ?
                    is_slot_empty(last_tmpl->u.tmpl.slots.items[0]) : 1;
                if (s0e) {
                    MtefNode *next = filtered[i + 1];
                    if (next->type == NODE_LINE && next->u.line.children.count > 0) {
                        int has_tmpl = 0;
                        for (int j = 0; j < next->u.line.children.count; j++) {
                            if (next->u.line.children.items[j]->type == NODE_TMPL ||
                                next->u.line.children.items[j]->type == NODE_MATRIX)
                                { has_tmpl = 1; break; }
                        }
                        if (!has_tmpl && last_tmpl->u.tmpl.slots.count > 0) {
                            last_tmpl->u.tmpl.slots.items[0] = next;
                            /* remove filtered[i+1] */
                            for (int j = i + 1; j < fcount - 1; j++) filtered[j] = filtered[j+1];
                            fcount--;
                            merge_safe++;
                            if (merge_safe > fcount + 10) break;
                            continue;
                        }
                    }
                }
            }
        }
        merge_safe = 0;
        i++;
    }

    /* Cursor garbage removal.
     * NOTE: Only applies when the last line is a SINGLE char in TF_VARIABLE
     * typeface AND the preceding line contains a bigop (integral etc.) or
     * is non-trivial.  Disabled for now to avoid false positives in PILE
     * environments where the last row genuinely has simple content. */
    /* if (fcount > 1 && is_cursor_garbage_line(filtered[fcount - 1]))
        fcount--; */

    /* Convert lines to LaTeX */
    char *line_strs[256];
    for (int li = 0; li < fcount; li++) {
        StringBuilder ls; sb_init(&ls);
        convert_node(filtered[li], prod_ver, &ls);
        line_strs[li] = sb_detach(&ls);
    }

    /* Single line — but check for multi-LINE children (PILE with 1 outer LINE
     * containing multiple inner LINEs, common in EQNEDT32 PILE encoding). */
    if (fcount == 1) {
        MtefNode *only = filtered[0];
        if (only && only->type == NODE_LINE) {
            int sub_lines = 0;
            for (int si = 0; si < only->u.line.children.count; si++) {
                MtefNode *ch = only->u.line.children.items[si];
                if (ch->type == NODE_LINE && !ch->u.line.is_null) sub_lines++;
            }
            if (sub_lines >= 2) {
                /* Expand sub-LINEs into separate logical lines */
                free(line_strs[0]);
                fcount = 0;
                for (int si = 0; si < only->u.line.children.count && fcount < 256; si++) {
                    MtefNode *ch = only->u.line.children.items[si];
                    if (ch->type == NODE_LINE && !ch->u.line.is_null) {
                        filtered[fcount++] = ch;
                    }
                }
                /* Re-convert expanded lines */
                for (int li = 0; li < fcount; li++) {
                    StringBuilder ls; sb_init(&ls);
                    convert_node(filtered[li], prod_ver, &ls);
                    line_strs[li] = sb_detach(&ls);
                }
                /* fall through to multi-line output below */
                goto multi_line_output;
            }
        }
        sb_append(sb, line_strs[0]);
        free(line_strs[0]);
        return;
    }

multi_line_output:

    /* Alignment */
    static const char *RELATION_OPS[] = {
        "=", "\\leq", "\\geq", "\\neq", "\\approx",
        "\\equiv", "\\sim", "\\propto", "<", ">",
    };

    const char *env;
    if (node->u.pile.halign == 20) env = "bmatrix";
    else if (node->u.pile.halign == 21) env = "pmatrix";
    else if (node->u.pile.halign == 22) env = "vmatrix";
    else if (node->u.pile.halign == 23) env = "Vmatrix";
    else if (node->u.pile.halign == 24) env = "Bmatrix";
    else if (node->u.pile.halign == 3) {
        env = "aligned";
        for (int li = 0; li < fcount; li++) {
            char *tex = line_strs[li];
            int inserted = 0;
            for (int oi = 0; oi < 10 && !inserted; oi++) {
                const char *op = RELATION_OPS[oi];
                char *pos = strstr(tex, op);
                if (pos && pos > tex) {
                    int offset = (int)(pos - tex);
                    int olen = (int)strlen(op);
                    int tlen = (int)strlen(tex);
                    char *new_tex = (char *)malloc(tlen + 2);
                    if (new_tex) {
                        memcpy(new_tex, tex, offset);
                        new_tex[offset] = '&';
                        memcpy(new_tex + offset + 1, tex + offset, tlen - offset + 1);
                        free(tex);
                        line_strs[li] = new_tex;
                        inserted = 1;
                    }
                }
            }
        }
    } else {
        env = "gathered";
    }

    sb_append(sb, "\\begin{"); sb_append(sb, env); sb_append(sb, "}\n");
    for (int li = 0; li < fcount; li++) {
        sb_append(sb, "  "); sb_append(sb, line_strs[li]);
        if (li < fcount - 1) sb_append(sb, " \\\\\n");
        else sb_append_char(sb, '\n');
        free(line_strs[li]);
    }
    sb_append(sb, "\\end{"); sb_append(sb, env); sb_append_char(sb, '}');
}

/* ============================================================
 * convert_matrix
 * ============================================================ */

static void convert_matrix(MtefNode *node, int prod_ver, StringBuilder *sb)
{
    int rows = node->u.matrix.rows;
    int cols = node->u.matrix.cols;
    if (rows == 0 || cols == 0) return;

    sb_append(sb, "\\begin{matrix}\n");
    for (int r = 0; r < rows; r++) {
        sb_append(sb, "  ");
        for (int c = 0; c < cols; c++) {
            int idx = r * cols + c;
            if (idx < node->u.matrix.elements.count) {
                convert_node(node->u.matrix.elements.items[idx], prod_ver, sb);
            }
            if (c < cols - 1) sb_append(sb, " & ");
        }
        if (r < rows - 1) sb_append(sb, " \\\\\n");
        else sb_append_char(sb, '\n');
    }
    sb_append(sb, "\\end{matrix}");
}

/* ============================================================
 * convert_node — main dispatch
 * ============================================================ */

static void convert_node(MtefNode *node, int prod_ver, StringBuilder *sb)
{
    if (!node) return;
    switch (node->type) {
    case NODE_LINE:   convert_line(node, prod_ver, g_conv_arena, sb); break;
    case NODE_CHAR:   convert_char(node, prod_ver, sb); break;
    case NODE_TMPL:   convert_tmpl(node, prod_ver, sb); break;
    case NODE_PILE:   convert_pile(node, prod_ver, g_conv_arena, sb); break;
    case NODE_MATRIX: convert_matrix(node, prod_ver, sb); break;
    case NODE_SIZE:   break; /* ignored */
    case NODE_FONT:   break; /* ignored */
    case NODE_EMBELL: {
        /* Standalone EMBELL in template slots — apply to preceding content */
        int et = node->u.embell.embell_type;
        if (et >= 0 && et < EMBELL_MAP_N) {
            const char *pfx = EMBELL_MAP[et].prefix;
            const char *sfx = EMBELL_MAP[et].suffix;
            if (pfx[0]) sb_append(sb, pfx);
            if (sfx[0]) sb_append(sb, sfx);
        }
        break;
    }
    }
}

/* ============================================================
 * Tree dump, for debugging
 * ============================================================ */

#ifdef EQ2TEX_DEBUG
static void dump_tree(MtefNode *node, int depth)
{
    if (!node) { printf("%*sNULL\n", depth*2, ""); return; }
    switch (node->type) {
    case NODE_LINE:
        printf("%*sLINE%s children=%d\n", depth*2, "",
               node->u.line.is_null ? "(null)" : "",
               node->u.line.children.count);
        for (int i = 0; i < node->u.line.children.count; i++)
            dump_tree(node->u.line.children.items[i], depth+1);
        break;
    case NODE_CHAR:
        printf("%*sCHAR tf=%d code=0x%04x embells=%d\n", depth*2, "",
               node->u.ch.typeface, node->u.ch.char_code, node->u.ch.embells.count);
        break;
    case NODE_TMPL:
        printf("%*sTMPL sel=%d var=%d slots=%d\n", depth*2, "",
               node->u.tmpl.selector, node->u.tmpl.variation, node->u.tmpl.slots.count);
        for (int i = 0; i < node->u.tmpl.slots.count; i++) {
            printf("%*s  slot[%d]:\n", depth*2, "", i);
            dump_tree(node->u.tmpl.slots.items[i], depth+2);
        }
        break;
    case NODE_PILE:
        printf("%*sPILE lines=%d\n", depth*2, "", node->u.pile.lines.count);
        for (int i = 0; i < node->u.pile.lines.count; i++)
            dump_tree(node->u.pile.lines.items[i], depth+1);
        break;
    case NODE_MATRIX:
        printf("%*sMATRIX %dx%d\n", depth*2, "", node->u.matrix.rows, node->u.matrix.cols);
        break;
    case NODE_EMBELL:
        printf("%*sEMBELL type=%d\n", depth*2, "", node->u.embell.embell_type);
        break;
    case NODE_SIZE:
        printf("%*sSIZE type=%d\n", depth*2, "", node->u.size.size_type);
        break;
    case NODE_FONT:
        printf("%*sFONT\n", depth*2, "");
        break;
    }
}
#endif

/* ============================================================
 * Public API
 * ============================================================ */

char *mtef_to_latex_c(const uint8_t *data, size_t len)
{
    if (!data || len < 5) return NULL;

    /* Equation Native header detection (28-byte header) */
    if (len >= 30) {
        uint16_t cb_hdr = data[0] | (data[1] << 8);
        if (cb_hdr == 28 && (size_t)cb_hdr < len) {
            /* Check if data after header looks like MTEF (version 3, platform 0-2) */
            if (data[28] <= 5 && data[29] <= 2) {
                data += cb_hdr;
                len -= cb_hdr;
            }
        }
    }

    /* Stack-allocate arena (large!) — use heap instead for safety */
    Eq2TexArena *arena = (Eq2TexArena *)malloc(sizeof(Eq2TexArena));
    if (!arena) return NULL;
    arena_init(arena);

    int prod_ver = 0;
    MtefNode *root = mtef_parse(data, len, arena, &prod_ver);
    if (!root) { free(arena); return NULL; }

#ifdef EQ2TEX_DEBUG
    printf("=== Parse Tree ===\n");
    dump_tree(root, 0);
    printf("=== End Tree ===\n");
#endif

    StringBuilder sb;
    sb_init(&sb);

    /* Set global arena for nested convert_line calls */
    g_conv_arena = arena;
    g_convert_line_depth = 0;

    /* Check if root LINE contains PILE children */
    if (root->type == NODE_LINE) {
        int has_pile = 0, has_fence_or_tmpl = 0;
        for (int i = 0; i < root->u.line.children.count; i++) {
            MtefNode *c = root->u.line.children.items[i];
            if (c->type == NODE_PILE) has_pile = 1;
            if (c->type == NODE_TMPL &&
                (is_fence_selector(c->u.tmpl.selector) ||
                 is_decoration_selector(c->u.tmpl.selector)))
                has_fence_or_tmpl = 1;
        }

        if (has_pile && has_fence_or_tmpl) {
            /* PILE + fence/decoration: use convert_line for 3-pass processing
             * so fence merge can combine PILE into fence slot[0]. */
            convert_line(root, prod_ver, arena, &sb);
        } else if (has_pile) {
            /* PILE equation without fence: iterate children, special-case PILEs */
            for (int i = 0; i < root->u.line.children.count; i++) {
                MtefNode *c = root->u.line.children.items[i];
                if (c->type == NODE_PILE) {
                    convert_pile(c, prod_ver, arena, &sb);
                } else if (c->type != NODE_SIZE && c->type != NODE_FONT) {
                    convert_node(c, prod_ver, &sb);
                }
            }
        } else {
            /* Single-line equation: use convert_line for 3-pass processing */
            convert_line(root, prod_ver, arena, &sb);
        }
    } else {
        convert_node(root, prod_ver, &sb);
    }

    g_conv_arena = NULL;
    free(arena);

    /* Post-processing: normalize LaTeX output spacing */
    {
        char *result = sb_detach(&sb);

        /* In-place substitution: replace all occurrences of FROM with TO.
         * Requires strlen(TO) <= strlen(FROM). */
#define INPLACE_SUB(str, from, to) do {                              \
            const char *_f = (from), *_t = (to);                     \
            size_t _fl = strlen(_f), _tl = strlen(_t);               \
            char *_p = (str);                                         \
            while ((_p = strstr(_p, _f)) != NULL) {                  \
                memmove(_p + _tl, _p + _fl, strlen(_p + _fl) + 1);  \
                memcpy(_p, _t, _tl);                                  \
                _p += _tl;                                            \
            }                                                         \
        } while (0)

        /* Pass 1: merge triple \cdot → \cdots */
        INPLACE_SUB(result, " \\cdot  \\cdot  \\cdot ", "\\cdots ");

        /* Pass 2: remove spurious space before subscript / superscript brace
         *   e.g.  "\alpha _{n}"  →  "\alpha_{n}"
         *         "\sum _{n=1}^" →  "\sum_{n=1}^"  */
        INPLACE_SUB(result, " _{", "_{");
        INPLACE_SUB(result, " ^{", "^{");

        /* Pass 3: unary minus at the start of a brace group
         *   e.g.  "e^{ - jkR}"  →  "e^{-jkR}"
         *         "\dfrac{ - a}{b}" → "\dfrac{-a}{b}"  */
        INPLACE_SUB(result, "{ - ", "{-");

#undef INPLACE_SUB
        return result;
    }
}