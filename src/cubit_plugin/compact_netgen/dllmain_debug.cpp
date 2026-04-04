// DllMain debug: log DLL load/unload to file
#include <windows.h>
#include <cstdio>

BOOL APIENTRY DllMain(HMODULE hModule, DWORD reason, LPVOID lpReserved)
{
    FILE* f = fopen("C:\\compact_netgen_debug.log", "a");
    if (f) {
        switch (reason) {
        case DLL_PROCESS_ATTACH:
            fprintf(f, "DLL_PROCESS_ATTACH\n");
            break;
        case DLL_PROCESS_DETACH:
            fprintf(f, "DLL_PROCESS_DETACH\n");
            break;
        }
        fflush(f);
        fclose(f);
    }
    return TRUE;
}
