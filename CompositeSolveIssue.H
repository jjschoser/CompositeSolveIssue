#ifndef COMPOSITE_SOLVE_ISSUE_H_
#define COMPOSITE_SOLVE_ISSUE_H_

#include <AMReX_MLMG.H>

#ifdef AMREX_USE_HYPRE
#include <AMReX_Hypre.H>
#endif

class CompositeSolveIssue
{
public:

    CompositeSolveIssue (const std::string& name, const int nCell, const int maxLevel);

    void solve ();
    void writePlotfile () const;

private:
    const std::string m_name;

    const amrex::RealArray m_probLo = {AMREX_D_DECL(0.0, -10e-3, 0.0)};
    const amrex::RealArray m_probHi = {AMREX_D_DECL(48e-3, 38e-3, 48e-3)};

    const amrex::Real m_plateThickness = 7e-3;
    const amrex::Real m_arcRadius = 3e-3;

    const amrex::Real m_arcSigma = 37439.37498536;
    const amrex::Real m_airSigma = 1e-12;
    const amrex::RealArray m_plateSigma = {AMREX_D_DECL(3e4, 0.1, 1.0)};

    const int m_refRatio = 2;
    const int m_maxGridSize = 64;
    const int m_verbose = 3;
    const int m_maxIter = 20000;
    const int m_maxOrder = 2;
    const int m_maxCoarseningLevel = 0;
    const amrex::Real m_relTol = 1e-10;
    const amrex::Real m_absTol = 0.0;

#ifdef AMREX_USE_HYPRE
    amrex::Hypre::Interface m_hypre_interface = amrex::Hypre::Interface::ij;
#endif

    amrex::Vector<amrex::Geometry> m_geom;
    amrex::Vector<amrex::BoxArray> m_grids;
    amrex::Vector<amrex::DistributionMapping> m_dmap;

    amrex::Vector<amrex::MultiFab> m_phi;
    amrex::Vector<amrex::MultiFab> m_rhs;
    amrex::Vector<amrex::MultiFab> m_sigmaCC;
    amrex::Vector<amrex::Array<amrex::MultiFab, amrex::SpaceDim>> m_sigmaFC;

    amrex::Real getSigma (const amrex::RealArray& pos, const int idim) const;
};

#endif
