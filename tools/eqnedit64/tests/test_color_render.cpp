#include "equation_render.h"
#include "tex_parser.h"
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <windows.h>

int main() {
    const char* isolated = std::getenv("EQNEDIT64_ISOLATED_TEST_SESSION");
    if (!isolated || std::string(isolated) != "1") {
        std::cerr << "Colour rendering requires a disposable isolated font session\n";
        return 2;
    }
    auto root = eqnedit::parse_latex(
        "\\textcolor{red}{\\frac{x}{y}}+\\textcolor{blue}{\\sqrt{z}}+w");
    const auto svg = eqnedit::render_svg(*root);
    if (svg.find("fill=\"#FF0000\"") == std::string::npos ||
        svg.find("fill=\"#0000FF\"") == std::string::npos)
        throw std::runtime_error("SVG glyph/rule colour missing");
    BITMAPINFO info{};
    info.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
    info.bmiHeader.biWidth = 800; info.bmiHeader.biHeight = -300;
    info.bmiHeader.biPlanes = 1; info.bmiHeader.biBitCount = 32;
    info.bmiHeader.biCompression = BI_RGB;
    void* pixels = nullptr;
    HDC dc = CreateCompatibleDC(nullptr);
    HBITMAP bitmap = CreateDIBSection(dc, &info, DIB_RGB_COLORS, &pixels, nullptr, 0);
    HGDIOBJ old = SelectObject(dc, bitmap);
    auto check = [&]() {
        GdiFlush();
        unsigned red = 0, blue = 0;
        auto bytes = static_cast<unsigned char*>(pixels);
        for (int i = 0; i < 800 * 300; ++i) {
            const auto b = bytes[4*i], g = bytes[4*i+1], r = bytes[4*i+2];
            if (r > 150 && g < 80 && b < 80) ++red;
            if (b > 150 && g < 80 && r < 80) ++blue;
        }
        if (red < 10 || blue < 10) throw std::runtime_error("DIB/EMF lost colour ink");
    };
    RECT rect{0,0,800,300};
    FillRect(dc, &rect, static_cast<HBRUSH>(GetStockObject(WHITE_BRUSH)));
    eqnedit::draw_equation_gdi(*root, dc, 5, 5, 3, {}, nullptr, -1,
                             nullptr, -1, -1, false, false);
    check();
    HDC emf = CreateEnhMetaFileW(dc, nullptr, nullptr, nullptr);
    eqnedit::draw_equation_gdi(*root, emf, 5, 5, 3, {}, nullptr, -1,
                             nullptr, -1, -1, false, false, true);
    HENHMETAFILE image = CloseEnhMetaFile(emf);
    FillRect(dc, &rect, static_cast<HBRUSH>(GetStockObject(WHITE_BRUSH)));
    RECT emfRect{0,0,800,300};
    PlayEnhMetaFile(dc, image, &emfRect);
    check();
    DeleteEnhMetaFile(image);
    SelectObject(dc, old); DeleteObject(bitmap); DeleteDC(dc);
    std::cout << "PASS: SVG, GDI DIB and outline EMF colour ink\n";
}
