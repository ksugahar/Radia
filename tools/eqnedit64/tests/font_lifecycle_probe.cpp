// Diagnostic only: never register fonts on shared LAB/100 sessions.
// Deliberately does not link the editor or renderer, cache fonts, or retry.
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <cstdio>
#include <cwchar>

namespace {
void record(const char* stage, long result) {
    FILETIME time{};
    GetSystemTimePreciseAsFileTime(&time);
    ULARGE_INTEGER ticks{};
    ticks.LowPart = time.dwLowDateTime;
    ticks.HighPart = time.dwHighDateTime;
    std::printf("{\"stage\":\"%s\",\"pid\":%lu,\"filetime_100ns\":%llu,\"result\":%ld}\n",
                stage, GetCurrentProcessId(), ticks.QuadPart, result);
    std::fflush(stdout);
}
bool enabled(const wchar_t* name, const wchar_t* expected) {
    wchar_t value[16]{};
    const DWORD n = GetEnvironmentVariableW(name, value, 16);
    return n > 0 && n < 16 && std::wcscmp(value, expected) == 0;
}
}

int wmain(int argc, wchar_t** argv) {
    if (argc != 5) {
        std::fprintf(stderr, "Usage: font_lifecycle_probe FONT exit|remove|hold none|measure HOLD_SECONDS(1..60)\n");
        return 2;
    }
    const bool remove = std::wcscmp(argv[2], L"remove") == 0;
    const bool hold = std::wcscmp(argv[2], L"hold") == 0;
    const bool measure = std::wcscmp(argv[3], L"measure") == 0;
    wchar_t* end = nullptr;
    const long seconds = std::wcstol(argv[4], &end, 10);
    if ((!remove && !hold && std::wcscmp(argv[2], L"exit") != 0) ||
        (!measure && std::wcscmp(argv[3], L"none") != 0) ||
        !end || *end || seconds < 1 || seconds > 60) return 2;
    // This guard is not proof of isolation: the workflow must use fresh hosted VMs.
    if (!enabled(L"GITHUB_ACTIONS", L"true") ||
        !enabled(L"EQNEDIT64_ISOLATED_TEST_SESSION", L"1")) {
        std::fprintf(stderr, "Refusing font registration outside disposable diagnostic CI.\n");
        return 3;
    }
    const DWORD attributes = GetFileAttributesW(argv[1]);
    if (attributes == INVALID_FILE_ATTRIBUTES || (attributes & FILE_ATTRIBUTE_DIRECTORY)) return 4;
    constexpr DWORD flags = FR_PRIVATE | FR_NOT_ENUM;
    record("register.begin", 0);
    const int count = AddFontResourceExW(argv[1], flags, nullptr);
    record("register.end", count); // API has no extended error information.
    if (count <= 0) return 5;
    int status = 0;
    if (measure) {
        record("measure.begin", 0);
        HDC dc = CreateCompatibleDC(nullptr);
        HFONT font = CreateFontW(-24, 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE,
            DEFAULT_CHARSET, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS,
            DEFAULT_QUALITY, DEFAULT_PITCH, L"Latin Modern Math");
        bool ok = false;
        if (dc && font) {
            HGDIOBJ old = SelectObject(dc, font);
            if (old && old != HGDI_ERROR) {
                wchar_t face[LF_FACESIZE]{};
                SIZE size{};
                ok = GetTextFaceW(dc, LF_FACESIZE, face) > 0 &&
                    std::wcscmp(face, L"Latin Modern Math") == 0 &&
                    GetTextExtentPoint32W(dc, L"(", 1, &size) && size.cx > 0;
                SelectObject(dc, old);
            }
        }
        if (font) DeleteObject(font);
        if (dc) DeleteDC(dc);
        record("measure.end", ok ? 1 : 0);
        if (!ok) status = 6;
    }
    if (hold) {
        record("hold.begin", seconds);
        Sleep(static_cast<DWORD>(seconds) * 1000);
        record("hold.end", 0);
    }
    if (remove) {
        record("remove.begin", 0);
        const BOOL ok = RemoveFontResourceExW(argv[1], flags, nullptr);
        record("remove.end", ok ? 1 : 0);
        if (!ok) status = 7;
    }
    record("return", status); // External observer must still record process exit.
    return status;
}
