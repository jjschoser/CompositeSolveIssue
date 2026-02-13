
#include <AMReX.H>
#include "CompositeSolveIssue.H"

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    {        
        amrex::Print() << "Running test with unigrid..." << std::endl;
        CompositeSolveIssue unigrid("test_unigrid", 128, 0);
        unigrid.solve();
        unigrid.writePlotfile();

        amrex::Print() << "Running test with AMR..." << std::endl;
        CompositeSolveIssue amr("test_AMR", 32, 2);
        amr.solve();
        amr.writePlotfile();
    }

    amrex::Finalize();

    return 0;
}
