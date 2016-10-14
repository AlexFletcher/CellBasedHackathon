
#ifndef TESTLONGRANGESIGNALLINGEXAMPLE_HPP_
#define TESTLONGRANGESIGNALLINGEXAMPLE_HPP_

#include <cxxtest/TestSuite.h>
#include "SmartPointers.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "PottsMeshGenerator.hpp"
#include "WildTypeCellMutationState.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "NoCellCycleModel.hpp"
#include "CellLabel.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "OnLatticeSimulation.hpp"
#include "AveragedSourceParabolicPde.hpp"
#include "ParabolicBoxDomainPdeModifier.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "SurfaceAreaConstraintPottsUpdateRule.hpp"
#include "LongRangeSignallingExamplePde.hpp"
// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestLongRangeSignallingExample : public AbstractCellBasedWithTimingsTestSuite
{
public:

    void TestParabolicBoxDomainPdeModifierWithPottsBasedMonolayer() throw (Exception)
    {
        double end_time = 100;

        // Create a PottsMesh
        unsigned cell_width = 4;
        unsigned domain_width = 200;
        unsigned num_cells_across = 5;
        PottsMeshGenerator<2> generator(domain_width, num_cells_across, cell_width, domain_width, num_cells_across, cell_width);
        PottsMesh<2>* p_mesh = generator.GetMesh();
        p_mesh->Translate(-0.5*(domain_width-(num_cells_across-1)*cell_width),-0.5*(domain_width-(num_cells_across-1)*cell_width));
        p_mesh->Scale(0.25,0.25);

        // Create a vector of cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_transit_type);
        MAKE_PTR(CellLabel, p_label);
//        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());
        for (unsigned i=0; i<num_cells_across*num_cells_across; i++)
        {
            NoCellCycleModel* p_cycle_model = new NoCellCycleModel();
            p_cycle_model->SetDimension(2);
            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            p_cell->GetCellData()->SetItem("Morphogen", 1.0);
		    p_cell->AddCellProperty(p_label);
            cells.push_back(p_cell);
        }

        // Create a PottsBasedCellPopulation
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetTemperature(0.1);

        // Create an OnLatticeSimulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestLongRangeSignallingExample");
        simulator.SetDt(1.0);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetEndTime(end_time);

        // Create some update rules and add them to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        simulator.AddUpdateRule(p_volume_constraint_update_rule);
//        MAKE_PTR(SurfaceAreaConstraintPottsUpdateRule<2>, p_surface_area_update_rule);
//        simulator.AddUpdateRule(p_surface_area_update_rule);
        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
        simulator.AddUpdateRule(p_adhesion_update_rule);

        // Create a PDE modifier and pass it to the simulation

        // Create PDE and boundary condition
        double diffusion_coefficient = 1.0;
        double uptake_rate = 1.0;
        double sink_rate = 1.0;
        MAKE_PTR_ARGS(LongRangeSignallingExamplePde<2>, p_pde, (cell_population, 1.0, diffusion_coefficient, uptake_rate, sink_rate));
        // making neumann boundary condition
        MAKE_PTR_ARGS( ConstBoundaryCondition<2>, p_bc, (0.0) );

        // Todo: make a modifier to Change Cell Labels etc depending on concentration.
        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(-10.0, -10.0);
        ChastePoint<2> upper(10.0, 10.0);
        ChasteCuboid<2> cuboid(lower, upper);

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, true, &cuboid));
        p_pde_modifier->SetDependentVariableName("Morphogen");

        simulator.AddSimulationModifier(p_pde_modifier);

        // Run the simulation
        simulator.Solve();
    }
};

#endif /*TESTLONGRANGESIGNALLINGEXAMPLE_HPP_*/

