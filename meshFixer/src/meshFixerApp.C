#include "fvCFD.H"
#include "meshFixer.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Foam::meshFixer fixer(mesh);
    fixer.fix();

    return 0;
}