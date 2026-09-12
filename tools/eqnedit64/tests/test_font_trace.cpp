#include "font_trace.h"
#include <filesystem>
#include <fstream>
#include <string>

int main(int argc, char** argv) {
    if (argc != 2) return 1;
    const auto directory = std::filesystem::absolute(argv[1]);
    std::filesystem::create_directories(directory);
    const auto file = directory / (L"font-" + std::to_wstring(GetCurrentProcessId()) + L".jsonl");
    if (std::filesystem::exists(file)) return 2;
    SetEnvironmentVariableW(L"EQNEDIT64_FONT_TRACE_DIR", nullptr);
    SetLastError(1234);
    eqnedit::font_trace("disabled");
    if (GetLastError() != 1234 || std::filesystem::exists(file)) return 3;
    SetEnvironmentVariableW(L"EQNEDIT64_FONT_TRACE_DIR", directory.c_str());
    SetLastError(1234);
    eqnedit::font_trace("test.begin", 17);
    if (GetLastError() != 1234) return 4;
    eqnedit::font_trace("test.end", -1);
    std::ifstream input(file);
    std::string first, second, extra;
    if (!std::getline(input, first) || !std::getline(input, second) ||
        std::getline(input, extra)) return 5;
    if (first.find("\"stage\":\"test.begin\",\"result\":17}") == std::string::npos ||
        second.find("\"stage\":\"test.end\",\"result\":-1}") == std::string::npos ||
        first.find("\"session\":") == std::string::npos) return 6;
    input.close();
    std::filesystem::remove(file);
    SetEnvironmentVariableW(L"EQNEDIT64_FONT_TRACE_DIR", nullptr);
    return 0;
}
