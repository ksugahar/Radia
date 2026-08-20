/*
 * gvml_clip.cpp -- see gvml_clip.h for why this format exists at all.
 *
 * Two things: a minimal ZIP writer, and the three XML parts of the package.
 *
 * The XML was read off what PowerPoint itself puts on the clipboard, and only
 * what is load-bearing was kept.  In particular the size does NOT have to be
 * repeated on every <m:r> the way PowerPoint writes it -- a paragraph default
 * is enough, measured by pasting both into a real slide and reading the run
 * size back through PowerPoint's own object model.  That is what lets the OMML
 * writer stay untouched.
 */
#include "gvml_clip.h"

#include "mtef_omml.h"

#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

namespace mtef {

const double kPasteSizePt = 24.0;

namespace {

/* ---- a ZIP with stored entries ------------------------------------------ */

uint32_t crc32_of(const std::string& s) {
    static uint32_t table[256];
    static bool built = false;
    if (!built) {
        for (uint32_t i = 0; i < 256; ++i) {
            uint32_t c = i;
            for (int k = 0; k < 8; ++k)
                c = (c & 1) ? (0xEDB88320u ^ (c >> 1)) : (c >> 1);
            table[i] = c;
        }
        built = true;
    }
    uint32_t c = 0xFFFFFFFFu;
    for (unsigned char ch : s) c = table[(c ^ ch) & 0xFF] ^ (c >> 8);
    return c ^ 0xFFFFFFFFu;
}

void put16(std::string& out, uint32_t v) {
    out.push_back(char(v & 0xFF));
    out.push_back(char((v >> 8) & 0xFF));
}

void put32(std::string& out, uint32_t v) {
    put16(out, v & 0xFFFF);
    put16(out, (v >> 16) & 0xFFFF);
}

struct Entry {
    std::string name, data;
    uint32_t crc = 0, offset = 0;
};

std::string zip_of(std::vector<Entry> entries) {
    std::string out;
    for (Entry& e : entries) {
        e.crc = crc32_of(e.data);
        e.offset = uint32_t(out.size());
        put32(out, 0x04034B50);              /* local file header      */
        put16(out, 20);                      /* version needed         */
        put16(out, 0);                       /* flags                  */
        put16(out, 0);                       /* method: stored         */
        put16(out, 0);                       /* time                   */
        put16(out, 0x21);                    /* date: 1980-01-01       */
        put32(out, e.crc);
        put32(out, uint32_t(e.data.size()));
        put32(out, uint32_t(e.data.size()));
        put16(out, uint32_t(e.name.size()));
        put16(out, 0);                       /* extra field length     */
        out += e.name;
        out += e.data;
    }

    const uint32_t dir_at = uint32_t(out.size());
    for (const Entry& e : entries) {
        put32(out, 0x02014B50);              /* central directory      */
        put16(out, 20);                      /* version made by        */
        put16(out, 20);                      /* version needed         */
        put16(out, 0);
        put16(out, 0);
        put16(out, 0);
        put16(out, 0x21);
        put32(out, e.crc);
        put32(out, uint32_t(e.data.size()));
        put32(out, uint32_t(e.data.size()));
        put16(out, uint32_t(e.name.size()));
        put16(out, 0);                       /* extra                  */
        put16(out, 0);                       /* comment                */
        put16(out, 0);                       /* disk number            */
        put16(out, 0);                       /* internal attributes    */
        put32(out, 0);                       /* external attributes    */
        put32(out, e.offset);
        out += e.name;
    }
    const uint32_t dir_size = uint32_t(out.size()) - dir_at;

    put32(out, 0x06054B50);                  /* end of central dir     */
    put16(out, 0);
    put16(out, 0);
    put16(out, uint32_t(entries.size()));
    put16(out, uint32_t(entries.size()));
    put32(out, dir_size);
    put32(out, dir_at);
    put16(out, 0);                           /* comment length         */
    return out;
}

/* ---- the parts ---------------------------------------------------------- */

const char* kXmlHead =
    "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\r\n";

const char* kMathNs =
    "http://schemas.openxmlformats.org/officeDocument/2006/math";
const char* kDrawNs =
    "http://schemas.openxmlformats.org/drawingml/2006/main";
const char* kA14Ns = "http://schemas.microsoft.com/office/drawing/2010/main";
const char* kCanvasNs =
    "http://schemas.openxmlformats.org/drawingml/2006/lockedCanvas";

std::string content_types() {
    return std::string(kXmlHead) +
        "<Types xmlns=\"http://schemas.openxmlformats.org/package/2006/"
        "content-types\">"
        "<Default Extension=\"rels\" ContentType=\"application/vnd."
        "openxmlformats-package.relationships+xml\"/>"
        "<Default Extension=\"xml\" ContentType=\"application/xml\"/>"
        "<Override PartName=\"/clipboard/drawings/drawing1.xml\" ContentType="
        "\"application/vnd.openxmlformats-officedocument.drawing+xml\"/>"
        "</Types>";
}

std::string root_rels() {
    return std::string(kXmlHead) +
        "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/"
        "relationships\"><Relationship Id=\"rId1\" Type=\"http://schemas."
        "openxmlformats.org/officeDocument/2006/relationships/drawing\" "
        "Target=\"clipboard/drawings/drawing1.xml\"/></Relationships>";
}

std::string drawing(const std::string& omml, int half_points) {
    /* The canvas is a nominal box; PowerPoint autofits the shape to the text
     * on paste, so these numbers only have to be sane.  EMU: 914400 per inch. */
    const char* cx = "2286000";
    const char* cy = "762000";
    char sz[16];
    std::snprintf(sz, sizeof(sz), "%d", half_points);

    std::string s = kXmlHead;
    s += std::string("<a:graphic xmlns:a=\"") + kDrawNs + "\">";
    s += std::string("<a:graphicData uri=\"") + kCanvasNs + "\">";
    s += std::string("<lc:lockedCanvas xmlns:lc=\"") + kCanvasNs + "\">";
    s += "<a:nvGrpSpPr><a:cNvPr id=\"0\" name=\"\"/><a:cNvGrpSpPr/>"
         "</a:nvGrpSpPr>";
    s += std::string("<a:grpSpPr><a:xfrm><a:off x=\"0\" y=\"0\"/><a:ext cx=\"") +
         cx + "\" cy=\"" + cy + "\"/><a:chOff x=\"0\" y=\"0\"/><a:chExt cx=\"" +
         cx + "\" cy=\"" + cy + "\"/></a:xfrm></a:grpSpPr>";
    s += "<a:sp><a:nvSpPr><a:cNvPr id=\"1\" name=\"Equation\"/>"
         "<a:cNvSpPr txBox=\"1\"/></a:nvSpPr>";
    s += std::string("<a:spPr><a:xfrm><a:off x=\"0\" y=\"0\"/><a:ext cx=\"") +
         cx + "\" cy=\"" + cy + "\"/></a:xfrm><a:prstGeom prst=\"rect\">"
         "<a:avLst/></a:prstGeom><a:noFill/></a:spPr>";
    s += "<a:txSp><a:txBody>";
    s += "<a:bodyPr wrap=\"none\" rtlCol=\"0\"><a:spAutoFit/></a:bodyPr>";
    /* The size, in the two places that carry it: the level-1 default and the
     * paragraph.  Repeating it on every run, the way PowerPoint writes it, was
     * measured to make no difference. */
    s += std::string("<a:lstStyle><a:lvl1pPr><a:defRPr sz=\"") + sz +
         "\"/></a:lvl1pPr></a:lstStyle>";
    s += std::string("<a:p><a:pPr><a:defRPr sz=\"") + sz + "\"/></a:pPr>";
    s += std::string("<a14:m xmlns:a14=\"") + kA14Ns + "\">" + omml + "</a14:m>";
    s += std::string("<a:endParaRPr lang=\"en-US\" sz=\"") + sz + "\"/></a:p>";
    s += "</a:txBody><a:useSpRect/></a:txSp>";
    s += "</a:sp></lc:lockedCanvas></a:graphicData></a:graphic>";
    return s;
}

}  // namespace

std::string tex_to_gvml(const std::string& latex, double size_pt,
                        bool display) {
    OmmlOptions opt;
    opt.display = display;
    /* The a14:m wrapper does not declare the math namespace, so the fragment
     * has to carry its own. */
    opt.declare_namespace = true;
    const std::string omml = tex_to_omml(latex, opt);
    if (omml.empty()) return std::string();

    int half_points = int(size_pt * 100.0 + 0.5);   /* DrawingML: 1/100 pt */
    if (half_points < 100) half_points = 100;

    std::vector<Entry> parts;
    parts.push_back({"[Content_Types].xml", content_types(), 0, 0});
    parts.push_back({"_rels/.rels", root_rels(), 0, 0});
    parts.push_back({"clipboard/drawings/drawing1.xml",
                     drawing(omml, half_points), 0, 0});
    return zip_of(parts);
}

}  // namespace mtef
