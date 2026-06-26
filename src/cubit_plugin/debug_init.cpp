// debug_init.cpp — DllMain override to diagnose static init failures
#ifdef _WIN32
#include <windows.h>
#include <cstdio>

// Write debug messages to a file since we're in DllMain context
static void dbg(const char* msg)
{
    FILE* f = fopen("C:/temp/ccm_dllmain.log", "a");
    if (f) {
        fprintf(f, "%s\n", msg);
        fflush(f);
        fclose(f);
    }
}

// Custom DllMain to wrap CRT initialization
BOOL APIENTRY DllMain(HMODULE hModule, DWORD reason, LPVOID lpReserved)
{
    switch (reason) {
    case DLL_PROCESS_ATTACH:
        dbg("DLL_PROCESS_ATTACH: start");
        break;
    case DLL_PROCESS_DETACH:
        dbg("DLL_PROCESS_DETACH");
        break;
    }
    return TRUE;
}
#endif
