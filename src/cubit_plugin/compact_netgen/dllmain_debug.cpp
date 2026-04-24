// DllMain debug: log DLL load/unload to file
//
// Per-user filename under C:\temp so non-admin users can write their
// own log (shared C:\compact_netgen_debug.log at the root was
// effectively Admin-only — all non-admin fopen attempts silently
// returned NULL because C:\ root is Users:ReadAndExecute).  See
// CLAUDE.md Temp Directory Policy.
//
// NOTE: DllMain runs under the loader lock — avoid heavy Windows API
// calls here.  GetEnvironmentVariableA on USERNAME is safe (kernel32
// is already loaded).
#include <windows.h>
#include <cstdio>

BOOL APIENTRY DllMain(HMODULE hModule, DWORD reason, LPVOID lpReserved)
{
    char user[64];
    DWORD n = GetEnvironmentVariableA("USERNAME", user, sizeof(user));
    if (n == 0 || n >= sizeof(user)) {
        user[0] = 'u'; user[1] = 'n'; user[2] = 'k';
        user[3] = 'n'; user[4] = 'o'; user[5] = 'w';
        user[6] = 'n'; user[7] = '\0';
    } else {
        // Sanitize: remove path separators / special chars
        for (DWORD i = 0; i < n; ++i) {
            char c = user[i];
            if (!((c >= 'A' && c <= 'Z') ||
                  (c >= 'a' && c <= 'z') ||
                  (c >= '0' && c <= '9') ||
                  c == '-' || c == '_' || c == '.')) {
                user[i] = '_';
            }
        }
    }
    char path[MAX_PATH];
    snprintf(path, sizeof(path),
             "C:\\temp\\compact_netgen_debug_%s.log", user);
    FILE* f = fopen(path, "a");
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
