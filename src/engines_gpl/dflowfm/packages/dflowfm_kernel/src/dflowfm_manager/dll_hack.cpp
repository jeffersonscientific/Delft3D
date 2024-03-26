#include <windows.h>
#include <libloaderapi.h>
#include <errhandlingapi.h>
#include <iostream>

extern "C" {
    int nc_def_var(int ncid, const char* name, int xtype, int ndims, const int* dimidsp, int* varidp);
}

static void dll_hack_imp();

extern "C" {
    void dll_hack() {
        dll_hack_imp();
    }
}

void dll_hack_imp() {
    char path[MAX_PATH];
    HMODULE hm = NULL;

    if (GetModuleHandleEx(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS | 
            GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
            (LPCSTR) &nc_def_var, &hm) == 0) {
        int ret = GetLastError();
        std::cerr << "GetModuleHandle failed, error = " << ret;
        return;
    }
    if (GetModuleFileName(hm, path, sizeof(path)) == 0) {
        int ret = GetLastError();
        std::cerr << "GetModuleHandle failed, error = " << ret;
        return;
    }

    std::cout << "HARMEN: path is '" << path << "'.";
}
