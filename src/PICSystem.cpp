///////////////////////////////////////////////////////////////////////////////
//
// File: PICSystem.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Image warping solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include "PICSystem.h"
#include <MultiRegions/ContField.h>
#include <cstdlib>
#include <math.h> // #include <numbers> in C++20
#include <random>
#include <iostream>

namespace Nektar
{

std::string PICSystem::className =
    GetEquationSystemFactory().RegisterCreatorFunction(
        "PICSystem", PICSystem::create,
        "System for an Electrostatic PIC code.");

PICSystem::PICSystem(
    const LibUtilities::SessionReaderSharedPtr& session,
    const SpatialDomains::MeshGraphSharedPtr& graph)
    : UnsteadySystem(session, graph),
      AdvectionSystem(session, graph),
      m_driftVel(2),
      m_particlePositions(2),
      m_particleVelocities(2)
{
}

void PICSystem::v_InitObject()
{
    AdvectionSystem::v_InitObject();

    // Since we are starting from a setup where each field is defined to be a
    // discontinuous field (and thus support DG), the first thing we do is to
    // recreate the phi field so that it is continuous to support the Poisson
    // solve. Note that you can still perform a Poisson solve using a
    // discontinuous field, which is done via the hybridisable discontinuous
    // Galerkin (HDG) approach.
    m_fields[2] = MemoryManager<MultiRegions::ContField>
        ::AllocateSharedPtr(
            m_session, m_graph, m_session->GetVariable(2), true, true);

    // Tell UnsteadySystem to only integrate first two fields in time
    // (i.e. vorticity and density), ignoring electrostatic potential since .
    m_intVariables = { 0, 1 };

    // Assign storage for drift velocity.
    for (int i = 0; i < 2; ++i)
    {
        m_driftVel[i] = Array<OneD, NekDouble>(m_fields[i]->GetNpoints());
    }

    // Load constant alpha from the session file.
    ASSERTL0(m_session->DefinesParameter("alpha"),
             "Session file should define parameter alpha.");
    m_session->LoadParameter("alpha", m_alpha, 1.0);

    // Load constant kappa from the session file.
    ASSERTL0(m_session->DefinesParameter("kappa"),
             "Session file should define parameter kappa.");
    m_session->LoadParameter("kappa", m_kappa, 1.0);

    m_fields[3] = MemoryManager<MultiRegions::ContField>
        ::AllocateSharedPtr(
            m_session, m_graph, m_session->GetVariable(3), true, true);

    // Assign storage for particle positions.
    for (int i = 0; i < 2; ++i)
    {
        m_particlePositions[i] = Array<OneD, NekDouble>(m_numberMacroParticles);
        m_particleVelocities[i] = Array<OneD, NekDouble>(m_numberMacroParticles);
    }

    // Load constant charge from the session file.
    ASSERTL0(m_session->DefinesParameter("charge"),
             "Session file should define parameter charge.");
    m_session->LoadParameter("charge", m_charge, 1.0);

    // Load constant mass from the session file.
    ASSERTL0(m_session->DefinesParameter("mass"),
             "Session file should define parameter mass.");
    m_session->LoadParameter("mass", m_mass, 1.0);

    // Load constant numberDensity from the session file.
    ASSERTL0(m_session->DefinesParameter("numberDensity"),
             "Session file should define parameter numberDensity.");
    m_session->LoadParameter("numberDensity", m_numberDensity, 1.0);

    // Load constant numberDensity from the session file.
    ASSERTL0(m_session->DefinesParameter("numberMacroparticles"),
             "Session file should define parameter numberMacroparticles.");
    m_session->LoadParameter("numberMacroparticles", m_numberMacroParticles, 1.0);


    // Type of advection class to be used. By default, we only support the
    // discontinuous projection, since this is the only approach we're
    // considering for this solver.
    switch(m_projectionType)
    {
        // Discontinuous field
        case MultiRegions::eDiscontinuous:
        {
            // Do not forwards transform initial condition.
            m_homoInitialFwd = false;

            // Define the normal velocity fields.
            if (m_fields[0]->GetTrace())
            {
                m_traceVn = Array<OneD, NekDouble>(GetTraceNpoints());
            }

            // The remainder of this code is fairly generic boilerplate for the
            // DG setup.
            std::string advName, riemName;

            // Load what type of advection we want to use -- in theory we also
            // support flux reconstruction for quad-based meshes, or you can use
            // a standard convective term if you were fully continuous in
            // space. Default is DG.
            m_session->LoadSolverInfo("AdvectionType", advName, "WeakDG");

            // Create an advection object of the type above using the factory
            // pattern.
            m_advObject = SolverUtils::
                GetAdvectionFactory().CreateInstance(advName, advName);

            // The advection object needs to know the flux vector being
            // calculated: this is done with a callback.
            m_advObject->SetFluxVector(&PICSystem::GetFluxVector, this);

            // Repeat the above for the Riemann solver: in this case we use an
            // upwind by default. The solver also needs to know the trace
            // normal, which we again implement using a callback.
            m_session->LoadSolverInfo("UpwindType", riemName, "Upwind");
            m_riemannSolver =
                GetRiemannSolverFactory().CreateInstance(riemName, m_session);
            m_riemannSolver->SetScalar(
                "Vn", &PICSystem::GetNormalVelocity, this);

            // Tell the advection object about the Riemann solver to use, and
            // then get it set up.
            m_advObject->SetRiemannSolver(m_riemannSolver);
            m_advObject->InitObject(m_session, m_fields);

            // Now fields are initialised with values it is safe to
            // initialise particles from rho


            break;
        }

        default:
        {
            ASSERTL0(false,
                     "Unsupported projection type: only discontinuous"
                     " projection supported.");
            break;
        }
    }

    ASSERTL0(m_explicitAdvection,
             "This solver only supports explicit-in-time advection.");

    // The m_ode object defines the timestepping to be used, and lives in the
    // SolverUtils::UnsteadySystem class. For explicit solvers, you need to
    // supply a right-hand side function, and a projection function (e.g. for
    // continuous Galerkin this would be an assembly-type operation to ensure
    // C^0 connectivity). These are done again through callbacks.
    m_ode.DefineOdeRhs    (&PICSystem::ExplicitTimeInt, this);
    m_ode.DefineProjection(&PICSystem::DoOdeProjection, this);

    this->DoInitialise(); // this populates the fields with ICs

    this->InitialiseParticles();

}

void PICSystem::InitialiseParticles() {

  std::mt19937_64 rng;
  std::uniform_real_distribution<NekDouble> f01(0, 1);
  auto x = Array<OneD, NekDouble>(2);
  const auto twopi = 2 * M_PI;
  const auto rho = m_fields[3]; // charge density rho
  const auto rhoPhys = rho->GetPhys(); // TODO: is this right?
  for (int j = 0; j < m_numberMacroParticles; ++j) {

    while (true) { // load this particle
      const auto f = f01(rng);
      int cellIndex = -1;
      while (cellIndex < 0) {
        for (int d = 0; d < 2; ++d) {
          x[d] = f01(rng) * 5 - 1; // TODO: distribute particles evenly in grid
        }
        cellIndex = rho->GetExpIndex(x);
      }
      // Chris has some code in his tensor product code that looks like this?
      const auto rho_x = rho->PhysEvaluate(x, rhoPhys); // TODO: how to do this?!

      if (f <= rho_x) { // accept the particle position
        for (int d = 0; d < 2; ++d) {
          m_particlePositions[d][j] = x[d];
        }
        const auto r1 = f01(rng); // 1st random number for Box-Muller transform
        const auto r2 = f01(rng); // 2nd random number for Box-Muller transform
        // two independent normally distributed random numbers
        const auto commonPart = std::sqrt(-2 * std::log(r1));
        m_particleVelocities[0][j] = commonPart * std::cos(twopi * r2);
        m_particleVelocities[1][j] = commonPart * std::sin(twopi * r2);
        break;
      }
    }
  }

  // m_fields[3] is rho
  const auto volume = 1.0; // TODO get proper answer for volume of grid
  const auto numberPhysicalParticles = this->m_numberDensity * volume;
  this->m_particleWeight = numberPhysicalParticles / this->m_numberMacroParticles;
}

void PICSystem::Deposit() {

  auto x = Array<OneD, NekDouble>(2);
  auto rho = m_fields[3];
  int currentCellIndex = -1;
  Array<OneD, NekDouble> cellVector(rho->GetNcoeffs());
  std::map<int, std::stack<int>> cellIndexToParticleIndex();
  // assume particles are sorted already by increasing cellIndex
  for (int j = 0; j < m_numberMacroParticles; ++j) {
    for (int d = 0; d < 2; ++d) {
      x[d] = m_particlePositions[d][j];
    }
    const auto particleCellIndex = rho->GetExpIndex(x);
    cellIndexToParticleIndex[particleCellIndex].push(j);
  }
  for (auto& cellIndexParticleIndicesPair : cellIndexToParticleIndex) {
    const auto cellIndex = cellIndexParticleIndicesPair.first;
    const auto particleIndices = cellIndexParticleIndicesPair.second;
    while (!particleIndices.empty()) {
      const auto particleIndex = particleIndices.top();
      for (int d = 0; d < 2; ++d) {
        x[d] = m_particlePositions[d][particleIndex];
      }
// here be dragons
// library/StdRegions/StdQuadExp.cpp is useful
// need to transform x in to this cell's reference space
// evaluate bases at barycentric coordinates(?)
      Array<OneD, NekDouble> basesAtX = rho->StdEvaluate(x, rho); // TODO: fix this. evaluate all bases in cell at position x
      if (particleCellIndex != currentCellIndex) {
        // this is a new cell
        if (currentCellIndex > -1) { // we've completed evaluating basis functions in the cell
          // MultiplyByInvMassMatrix(cellVector, chargeDensityCellCoefficients);
        }
        currentCellIndex = particleCellIndex;
        cellVector = basesAtX;
      } else {
        cellVector = cellVector + basesAtX;
      }
      particleIndices.pop();
    }

  }

}





/**
 * @brief Evaluate the right-hand side of the ODE system used to integrate in
 * time.
 *
 * This routine performs the bulk of the work in this class, and essentially
 * computes the right hand side term of the generalised ODE system
 *
 * \f\[ \frac{\partial \mathbf{u}}{\partial t} = \mathbf{R}(\mathbf{u}) \f\]
 *
 * The order of operations is as follows:
 *
 * - First, compute the electrostatic potential \f$ \phi \f$, given the 
 * - Using this, compute the drift velocity \f$ (\partial_y\phi,
 *   -\partial_x\phi).
 * - Then evaluate the \f$ \nabla\cdot\mathbf{F} \f$ operator using the
 *   advection object #m_advObject.
 * - Finally put this on the right hand side and evaluate the source terms for
 *   each field.
 *
 * The assumption here is that fields are ordered inside `m_fields` so that
 * field 0 is vorticity \f$ \zeta \f$, field 1 is number density \f$ n \f$, and
 * field 2 is electrostatic potential. Only \f$ \zeta \f$ and \f$ n \f$ are time
 * integrated.
 *
 * @param inarray    Array containing each field's current state.
 * @param outarray   The result of the right-hand side operator for each field
 *                   being time integrated.
 * @param time       Current value of time.
 */
void PICSystem::ExplicitTimeInt(
    const Array<OneD, const Array<OneD,NekDouble> > &inarray,
          Array<OneD,       Array<OneD,NekDouble> > &outarray,
    const NekDouble time)
{
    // nPts below corresponds to the total number of solution/integration
    // points: i.e. number of elements * quadrature points per element.
    int i, nPts = GetNpoints();

    // Set up factors for electrostatic potential solve. We support a generic
    // Helmholtz solve of the form (\nabla^2 - \lambda) u = f, so this sets
    // \lambda to zero.
    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorLambda] = 0.0;

    // Solve for phi. Output of this routine is in coefficient (spectral) space,
    // so backwards transform to physical space since we'll need that for the
    // advection step & computing drift velocity.
    // inarray[0] is vorticity
    // inarray[2] is phi
    m_fields[2]->HelmSolve(inarray[0], m_fields[2]->UpdateCoeffs(),
                           factors);
    m_fields[2]->BwdTrans (m_fields[2]->GetCoeffs(),
                           m_fields[2]->UpdatePhys());

    // Calculate drift velocity v_E: PhysDeriv takes input and computes spatial
    // derivatives.
    // this one is vx = d phi /dy & vy = d phi /dx. Sign is corrected for below
    m_fields[2]->PhysDeriv(m_fields[2]->GetPhys(), m_driftVel[1], m_driftVel[0]);

    // We frequently use vector math (Vmath) routines for one-line operations
    // like negating entries in a vector.
    // This is because v = [d phi /dy, -d phi /dx], so negate the vy term
    Vmath::Neg(nPts, m_driftVel[1], 1);

    // Do advection for zeta, n. The hard-coded '2' here indicates that we
    // should only advect the first two components of inarray.
    m_advObject->Advect(2, m_fields, m_driftVel, inarray, outarray, time);

    // Put advection term on the right hand side.
    for (i = 0; i < 2; ++i)
    {
        Vmath::Neg(nPts, outarray[i], 1);
    }

    // Add source term alpha*(phi - n) to right hand side.
    Array<OneD, NekDouble> sourceTerm(nPts);
    Vmath::Vsub(nPts, m_fields[2]->GetPhys(), 1, inarray[1], 1, sourceTerm, 1);
    Vmath::Smul(nPts, m_alpha, sourceTerm, 1, sourceTerm, 1);
    Vmath::Vadd(nPts, sourceTerm, 1, outarray[0], 1, outarray[0], 1);
    Vmath::Vadd(nPts, sourceTerm, 1, outarray[1], 1, outarray[1], 1);

    // Add source term -kappa * d(phi)/dy to n equation.
    Vmath::Svtvp(nPts, -m_kappa, m_driftVel[0], 1,
                 outarray[1], 1, outarray[1], 1);
}

/**
 * @brief Perform projection into correct polynomial space.
 *
 * This routine projects the @p inarray input and ensures the @p outarray output
 * lives in the correct space. Since we are hard-coding DG, this corresponds to
 * a simple copy from in to out, since no elemental connectivity is required and
 * the output of the RHS function is polynomial.
 */
void PICSystem::DoOdeProjection(
    const Array<OneD, const Array<OneD, NekDouble> > &inarray,
          Array<OneD,       Array<OneD, NekDouble> > &outarray,
    const NekDouble time)
{
    int nvariables = inarray.size(), npoints = GetNpoints();
    SetBoundaryConditions(time);

    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Vcopy(npoints, inarray[i], 1, outarray[i], 1);
    }
}

/**
 * @brief Compute the flux vector for this system.
 */
void PICSystem::GetFluxVector(
    const Array<OneD, Array<OneD, NekDouble> >               &physfield,
          Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux)
{
    ASSERTL1(flux[0].size() == m_driftVel.size(),
             "Dimension of flux array and velocity array do not match");

    int i , j;
    int nq = physfield[0].size();

    for (i = 0; i < flux.size(); ++i)
    {
        for (j = 0; j < flux[0].size(); ++j)
        {
            Vmath::Vmul(nq, physfield[i], 1, m_driftVel[j], 1,
                        flux[i][j], 1);
        }
    }
}

/**
 * @brief Compute the normal advection velocity for this system on the
 * trace/skeleton/edges of the 2D mesh.
 */
Array<OneD, NekDouble> &PICSystem::GetNormalVelocity()
{
    // Number of trace (interface) points
    int i;
    int nTracePts = GetTraceNpoints();

    // Auxiliary variable to compute the normal velocity
    Array<OneD, NekDouble> tmp(nTracePts);

    // Reset the normal velocity
    Vmath::Zero(nTracePts, m_traceVn, 1);

    // Compute dot product of velocity along trace with trace normals. Store in
    // m_traceVn.
    for (i = 0; i < m_driftVel.size(); ++i)
    {
        m_fields[0]->ExtractTracePhys(m_driftVel[i], tmp);

        Vmath::Vvtvp(nTracePts,
                     m_traceNormals[i], 1,
                     tmp,               1,
                     m_traceVn,         1,
                     m_traceVn,         1);
    }

    return m_traceVn;
}


}
