#include "CompositeSolveIssue.H"

#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

CompositeSolveIssue::CompositeSolveIssue (const std::string& name, const int nCell, const int maxLevel) : m_name(name)
{
    const int nlevels = maxLevel + 1;
    m_geom.resize(nlevels);
    m_grids.resize(nlevels);
    m_dmap.resize(nlevels);

    m_phi.resize(nlevels);
    m_rhs.resize(nlevels);
    m_sigmaCC.resize(nlevels);
    m_sigmaFC.resize(nlevels);

    RealBox rb(m_probLo, m_probHi);
    Array<int,SpaceDim> is_periodic{AMREX_D_DECL(0,0,0)};
    Geometry::Setup(&rb, 0, is_periodic.data());
    Box domain0(IntVect{AMREX_D_DECL(0,0,0)}, IntVect{AMREX_D_DECL(nCell-1,nCell-1,nCell-1)});
    Box domain = domain0;
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        m_geom[ilev].define(domain);
        domain.refine(m_refRatio);
    }

    // Refine the entire domain
    domain = domain0;
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        m_grids[ilev].define(domain);
        m_grids[ilev].maxSize(m_maxGridSize);
        domain.refine(m_refRatio);
    }

    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        m_dmap[ilev].define(m_grids[ilev]);
        m_phi[ilev].define(m_grids[ilev], m_dmap[ilev], 1, 1);
        m_rhs[ilev].define(m_grids[ilev], m_dmap[ilev], 1, 0);
        m_sigmaCC[ilev].define(m_grids[ilev], m_dmap[ilev], SpaceDim, 1);
        for (int idim = 0; idim < SpaceDim; ++idim)
        {
            const BoxArray& ba = convert(m_grids[ilev],
                                                IntVect::TheDimensionVector(idim));
            m_sigmaFC[ilev][idim].define(ba, m_dmap[ilev], 1, 0);
        }
    }
}

void
CompositeSolveIssue::solve ()
{
    LPInfo info;
    info.setMaxCoarseningLevel(m_maxCoarseningLevel);

    const int nlevels = m_geom.size();

    MLABecLaplacian mlabec(m_geom, m_grids, m_dmap, info);

    mlabec.setMaxOrder(m_maxOrder);

    mlabec.setDomainBC({AMREX_D_DECL(LinOpBCType::Neumann,
                                     LinOpBCType::Neumann,
                                     LinOpBCType::Neumann)},
                        {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                      LinOpBCType::Dirichlet,
                                      LinOpBCType::Dirichlet)});

    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        m_phi[ilev].setVal(0.0);
        m_rhs[ilev].setVal(0.0);

        Box const& domain = m_geom[ilev].Domain();

        for (MFIter mfi(m_phi[ilev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            auto phifab = m_phi[ilev].array(mfi);

            if (bx.bigEnd(1) == domain.bigEnd(1))
            {
                Box const& bxhi = adjCellHi(bx, 1);
                ParallelFor(bxhi, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    phifab(i,j,k) = 1.0;
                });
            }
        }

        mlabec.setLevelBC(ilev, &(m_phi[ilev]));
    }

    mlabec.setScalars(0.0, 1.0);

    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        const Real* probLo = m_geom[ilev].ProbLo();
        const Real* dx = m_geom[ilev].CellSize();

        for (MFIter mfi(m_sigmaCC[ilev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox();
            auto sigmaCCfab = m_sigmaCC[ilev].array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const RealArray pos = {AMREX_D_DECL(probLo[0] + (i + 0.5) * dx[0],
                                                    probLo[1] + (j + 0.5) * dx[1],
                                                    probLo[2] + (k + 0.5) * dx[2])};
                AMREX_D_TERM(sigmaCCfab(i,j,k,0) = getSigma(pos, 0); ,
                             sigmaCCfab(i,j,k,1) = getSigma(pos, 1); ,
                             sigmaCCfab(i,j,k,2) = getSigma(pos, 2); )
            });
        }

        for (int idim = 0; idim < SpaceDim; ++idim)
        {
            const int iOffset = (idim == 0) ? 1 : 0;
            const int jOffset = (idim == 1) ? 1 : 0;
            const int kOffset = (idim == 2) ? 1 : 0;

            for (MFIter mfi(m_sigmaFC[ilev][idim], TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                auto sigmaCCfab = m_sigmaCC[ilev].array(mfi);
                auto sigmaFCfab = m_sigmaFC[ilev][idim].array(mfi);

                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    sigmaFCfab(i,j,k) = 2.0 * sigmaCCfab(i,j,k,idim) * sigmaCCfab(i-iOffset,j-jOffset,k-kOffset,idim) / (sigmaCCfab(i,j,k,idim) + sigmaCCfab(i-iOffset,j-jOffset,k-kOffset,idim));
                });
            }
        }
    }

    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        mlabec.setACoeffs(ilev, 0.0);
        mlabec.setBCoeffs(ilev, GetArrOfConstPtrs(m_sigmaFC[ilev]));
    }

    MLMG mlmg(mlabec);
    mlmg.setMaxIter(m_maxIter);
    mlmg.setVerbose(m_verbose);
#ifdef AMREX_USE_HYPRE
    mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
    mlmg.setHypreInterface(m_hypre_interface);
#endif

    mlmg.solve(GetVecOfPtrs(m_phi), GetVecOfConstPtrs(m_rhs), m_relTol, m_absTol);
}

void
CompositeSolveIssue::writePlotfile () const
{
    const int nlevels = m_geom.size();
    Vector<MultiFab> plotmf(nlevels);

    const int ncomp = 1 + SpaceDim;
    Vector<std::string> varname = {"phi", AMREX_D_DECL("sigma_xx", "sigma_yy", "sigma_zz")};

    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        plotmf[ilev].define(m_grids[ilev], m_dmap[ilev], ncomp, 0);
        MultiFab::Copy(plotmf[ilev], m_phi[ilev], 0, 0, 1, 0);
        MultiFab::Copy(plotmf[ilev], m_sigmaCC[ilev], 0, 1, SpaceDim, 0);
    }

    WriteMultiLevelPlotfile(m_name + "_plot", nlevels, amrex::GetVecOfConstPtrs(plotmf),
                            varname, m_geom, 0.0, Vector<int>(nlevels, 0),
                            Vector<IntVect>(nlevels, IntVect{m_refRatio}));
}

Real
CompositeSolveIssue::getSigma (const RealArray& pos, const int idim) const
{
    if (pos[1] < -m_plateThickness)
    {
        return m_airSigma;
    }
    else if (pos[1] < 0.0)
    {
        return m_plateSigma[idim];
    }

    Real r2 = 0.0;
    for (int d = 0; d < SpaceDim; ++d)
    {
        if (d != 1)
        {
            r2 += pos[d]*pos[d];
        }
    }
    return (r2 < m_arcRadius*m_arcRadius) ? m_arcSigma : m_airSigma;
}
