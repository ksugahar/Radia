// Visual evidence harness. Include the exact application source
// so neither font selection nor owner-draw is reimplemented by the reviewer.
#include "../src/eqnedt64_app.cpp"

bool save_proof(HBITMAP bitmap, const wchar_t* path) {
    IWICImagingFactory* factory = nullptr;
    IWICStream* stream = nullptr;
    IWICBitmap* source = nullptr;
    IWICBitmapEncoder* encoder = nullptr;
    IWICBitmapFrameEncode* frame = nullptr;
    HRESULT hr = CoCreateInstance(CLSID_WICImagingFactory, nullptr,
        CLSCTX_INPROC_SERVER, IID_PPV_ARGS(&factory));
    if (SUCCEEDED(hr)) hr = factory->CreateStream(&stream);
    if (SUCCEEDED(hr)) hr = stream->InitializeFromFilename(path, GENERIC_WRITE);
    if (SUCCEEDED(hr)) hr = factory->CreateBitmapFromHBITMAP(bitmap, nullptr,
        WICBitmapIgnoreAlpha, &source);
    if (SUCCEEDED(hr)) hr = factory->CreateEncoder(GUID_ContainerFormatPng, nullptr, &encoder);
    if (SUCCEEDED(hr)) hr = encoder->Initialize(stream, WICBitmapEncoderNoCache);
    if (SUCCEEDED(hr)) hr = encoder->CreateNewFrame(&frame, nullptr);
    if (SUCCEEDED(hr)) hr = frame->Initialize(nullptr);
    if (SUCCEEDED(hr)) hr = frame->WriteSource(source, nullptr);
    if (SUCCEEDED(hr)) hr = frame->Commit();
    if (SUCCEEDED(hr)) hr = encoder->Commit();
    if (frame) frame->Release();
    if (encoder) encoder->Release();
    if (source) source->Release();
    if (stream) stream->Release();
    if (factory) factory->Release();
    return SUCCEEDED(hr);
}

int wmain(int argc, wchar_t** argv) {
    DWORD session = 99;
    const char* isolated = std::getenv("EQNEDIT64_ISOLATED_TEST_SESSION");
    if (!ProcessIdToSessionId(GetCurrentProcessId(), &session) ||
        (session != 0 && (!isolated || std::strcmp(isolated, "1") != 0)) || argc != 2) return 90;
    if (FAILED(CoInitializeEx(nullptr, COINIT_APARTMENTTHREADED))) return 91;
    for (int dpi : {96, 144, 192}) {
        HFONT font = pick_button_font(MulDiv(17, dpi, 96));
        if (!font_resolves_to_face(font, EQNEDIT64_MATH_FONT_FACE)) return 243;
        auto proofPalettes = eqnedit::palettes();
        eqnedit::Palette selectors{"Selectors and persistent styles", "", 4, {}};
        for (const auto& palette : eqnedit::palettes())
            selectors.items.push_back({"", palette.face, palette.title, "", ""});
        for (const char* face : {"R x", "I x", "B x"})
            selectors.items.push_back({"", face, face, "", ""});
        proofPalettes.push_back(selectors);
        int index = 0;
        for (const auto& palette : proofPalettes) {
          const bool selector = size_t(index) == eqnedit::palettes().size();
          for (bool hot : {false, true}) {
            const int cellW = MulDiv(selector ? 52 : 34, dpi, 96);
            const int cellH = MulDiv(selector ? 30 : 28, dpi, 96);
            const int rows = (int(palette.items.size()) + palette.columns - 1) / palette.columns;
            const int width = palette.columns * cellW;
            const int height = rows * cellH;
            HDC dc = CreateCompatibleDC(nullptr);
            BITMAPINFO bi{};
            bi.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
            bi.bmiHeader.biWidth = width;
            bi.bmiHeader.biHeight = -height;
            bi.bmiHeader.biPlanes = 1;
            bi.bmiHeader.biBitCount = 32;
            void* bits = nullptr;
            HBITMAP bitmap = CreateDIBSection(dc, &bi, DIB_RGB_COLORS, &bits, nullptr, 0);
            if (!bitmap || !bits) return 92;
            HGDIOBJ old = SelectObject(dc, bitmap);
            RECT bounds{0, 0, width, height};
            FillRect(dc, &bounds, GetSysColorBrush(hot ? COLOR_HIGHLIGHT : COLOR_MENU));
            int itemIndex = 0;
            for (const auto& item : palette.items) {
                const int x = (itemIndex % palette.columns) * cellW;
                const int y = (itemIndex / palette.columns) * cellH;
                RECT cell{x, y, x + cellW, y + cellH};
                HFONT faceFont = font;
                if (selector && size_t(itemIndex) >= eqnedit::palettes().size())
                    faceFont = make_style_button_font(font, item.face == "I x",
                        item.face == "B x" ? FW_BOLD : FW_NORMAL);
                if (!palette_cell_draws_readably(faceFont, wide_utf8(item.face), dpi, item.command, hot, selector)) return 94;
                if (selector) draw_selector_face(dc, cell, faceFont, wide_utf8(item.face), dpi, hot);
                else draw_palette_cell(dc, cell, faceFont, wide_utf8(item.face), hot, item.command);
                if (faceFont != font) DeleteObject(faceFont);
                FrameRect(dc, &cell, GetSysColorBrush(COLOR_3DSHADOW));
                printf("%d,%d,%d,%s\n", dpi, index, itemIndex, item.command.c_str());
                ++itemIndex;
            }
            GdiFlush();
            SelectObject(dc, old);
            wchar_t path[1024];
            swprintf_s(path, L"%s/palette-%02d-%d-%s.png", argv[1], index, dpi, hot ? L"selected" : L"normal");
            const bool saved = save_proof(bitmap, path);
            DeleteObject(bitmap);
            DeleteDC(dc);
            if (!saved) return 93;
          }
            ++index;
        }
        DeleteObject(font);
    }
    CoUninitialize();
    return 0;
}
