#pragma once

// Opt-in diagnostics only. Never logs source text or changes font policy.
#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <cstdio>

namespace eqnedit {
inline void font_trace(const char* stage, long long result = 0) noexcept {
    const DWORD savedError = GetLastError();
    wchar_t directory[32768];
    const DWORD length = GetEnvironmentVariableW(L"EQNEDIT64_FONT_TRACE_DIR",
                                                directory, 32768);
    if (!length || length >= 32768) {
        SetLastError(savedError);
        return;
    }
    wchar_t filename[32768];
    if (swprintf_s(filename, L"%ls\\font-%lu.jsonl", directory,
                   GetCurrentProcessId()) < 0) {
        SetLastError(savedError);
        return;
    }
    HANDLE file = CreateFileW(filename, FILE_APPEND_DATA,
        FILE_SHARE_READ | FILE_SHARE_WRITE, nullptr, OPEN_ALWAYS,
        FILE_ATTRIBUTE_NORMAL, nullptr);
    if (file != INVALID_HANDLE_VALUE) {
        SYSTEMTIME now{};
        GetSystemTime(&now);
        DWORD session = MAXDWORD;
        ProcessIdToSessionId(GetCurrentProcessId(), &session);
        char line[512];
        const int count = sprintf_s(line,
            "{\"utc\":\"%04u-%02u-%02uT%02u:%02u:%02u.%03uZ\","
            "\"pid\":%lu,\"tid\":%lu,\"session\":%lu,\"stage\":\"%s\",\"result\":%lld}\n",
            now.wYear, now.wMonth, now.wDay, now.wHour, now.wMinute,
            now.wSecond, now.wMilliseconds, GetCurrentProcessId(),
            GetCurrentThreadId(), session, stage, result);
        DWORD written = 0;
        if (count > 0) WriteFile(file, line, static_cast<DWORD>(count), &written, nullptr);
        FlushFileBuffers(file);
        CloseHandle(file);
    }
    SetLastError(savedError);
}
struct FontTraceMainScope {
    FontTraceMainScope() { font_trace("app.enter"); }
    ~FontTraceMainScope() { font_trace("app.return"); }
};
} // namespace eqnedit
#else
namespace eqnedit {
inline void font_trace(const char*, long long = 0) noexcept {}
}
#endif
