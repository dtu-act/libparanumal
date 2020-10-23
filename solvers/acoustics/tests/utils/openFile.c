#include <acousticsTests.h>

int openFile(std::string filepath, FILE **iFP) {    

    char cwd[PATH_MAX];
    if (getcwd(cwd, sizeof(cwd)) == NULL) {
        perror("getcwd() error");
        return -1;
    }

    *iFP = fopen((char*)filepath.c_str(),"r");

    if (*iFP == NULL) {
        printf("ERROR: file could not be opened %s/%s)\n", cwd, (char*)filepath.c_str());
        return -1;
    }

    return 0;
}