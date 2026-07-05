#include "PhysiMeSS_environment.h"
#include "PhysiMeSS_agent.h"
#include "PhysiMeSS_fibre.h"
#include "PhysiMeSS_cell.h"

#include "../../BioFVM/BioFVM_vector.h"

using namespace PhysiCell;
using namespace BioFVM;

PhysiMeSS_Environment physimess_environment;

void PhysiMeSS_Environment::compute_velocities(double dt)
{
    for (auto* pCell : *all_cells)
    {
		PhysiMeSS_Agent* agent = dynamic_cast<PhysiMeSS_Agent*>(pCell->get_mechanics_implementation()); 
            
        double movement_threshold = pCell->custom_data["fibre_stuck_threshold"];
        if (!isFibre(pCell) && agent->motility_data.is_motile) {
        
            // Here I changed this, because here we don't have access to the old position, and I didn't want to track the old position
            // So I'm using the previous velocity, which is not exactly the same (because of Adams-Bashforth), but is a good proxy
            // if (dist(pCell->old_position, pCell->position) < movement_threshold) {
            if (norm(agent->previous_velocity)*mechanics_dt < movement_threshold) {
                static_cast<PhysiMeSS_Cell*>(pCell)->get_physimess_agent()->stuck_counter++;
            } else {
                static_cast<PhysiMeSS_Cell*>(pCell)->get_physimess_agent()->stuck_counter = 0;
            }
        }
        
        if( pCell->functions.add_cell_basement_membrane_interactions )
        {
            pCell->functions.add_cell_basement_membrane_interactions(pCell, pCell->phenotype,dt);
        }
        
        agent->simple_pressure = 0.0;
        agent->neighbors.clear(); // new 1.8.0
        
        if (!isFibre(pCell)) {
            //First check the neighbors in my current voxel
            for (auto* neighbor: pCell->nearby_cells())
            {
                auto* neighbor_agent = dynamic_cast<PhysiMeSS_Agent*>(neighbor->get_mechanics_implementation());
                if (!isFibre(neighbor)) {
                    agent->add_potentials(neighbor_agent);
                }
            }
        } else {
            // Count crosslinks
            static_cast<PhysiMeSS_Fibre*>(pCell)->get_physimess_agent()->X_crosslink_count = 0;
            if (static_cast<PhysiMeSS_Fibre*>(pCell)->get_physimess_agent()->fibres_crosslinkers.size() > 0){
                static_cast<PhysiMeSS_Fibre*>(pCell)->get_physimess_agent()->X_crosslink_count = static_cast<PhysiMeSS_Fibre*>(pCell)->get_physimess_agent()->fibres_crosslinkers.size();
            }

        }

        // std::cout << " AGENT " << pCell->type_name << " " << pCell->ID << " has " ;
        //add potentials between pCell and its neighbors
        for (auto* neighbor : agent->physimess_neighbors)
        {
            // std::cout << neighbor->type_name << " " << neighbor->ID << " " ;
            
            // if( this->ID == other_agent->ID )
            if( pCell != neighbor->get_cell() )
            { 
                if (!isFibre(pCell) && !isFibre(neighbor->get_cell())) {
                    //Already done above
                    continue;
                } else 
                if (!isFibre(pCell) && isFibre(neighbor->get_cell())) {
                    dynamic_cast<PhysiMeSS_CellAgent*>(agent)->add_potentials_from_fibre(dynamic_cast<PhysiMeSS_FibreAgent*>(neighbor));
                } else  if (isFibre(pCell) && !isFibre(neighbor->get_cell())) {
                    dynamic_cast<PhysiMeSS_FibreAgent*>(agent)->add_potentials_from_cell(dynamic_cast<PhysiMeSS_CellAgent*>(neighbor));
                } else if (isFibre(pCell) && isFibre(neighbor->get_cell())) {
                    dynamic_cast<PhysiMeSS_FibreAgent*>(agent)->add_potentials_from_fibre(dynamic_cast<PhysiMeSS_FibreAgent*>(neighbor));
                } else {
                    // std::cout << " WARNING: interaction between errant cell-types has been called : " << pCell->type_name << ", " << neighbor->type_name << std::endl;
                    return;
                }
            }
        }
        // std::cout << std::endl;

        if (!isFibre(pCell)) {
            
            PhysiMeSS_CellAgent* cellAgent = dynamic_cast<PhysiMeSS_CellAgent*>(agent);
            int stuck_threshold = 10;
            int unstuck_threshold = 1;

            if (cellAgent->stuck_counter == stuck_threshold){
                // std::cout << "!HELP! cell " << pCell->ID << " gets stuck at time "
                // << PhysiCell_globals.current_time << std::endl;
                cellAgent->stuck_counter = 0;
                cellAgent->unstuck_counter = 1;
            }

            if (1 <= cellAgent->unstuck_counter && cellAgent->unstuck_counter < unstuck_threshold+1) {
                /*std::cout << " getting unstuck at time "
                << PhysiCell_globals.current_time << std::endl;*/
                cellAgent->unstuck_counter++;
                cellAgent->force_update_motility_vector(dt);
                cellAgent->velocity += cellAgent->motility_data.motility_vector;
            }
            else {
                cellAgent->update_motility_vector(dt);
                cellAgent->velocity += cellAgent->motility_data.motility_vector;
            }

            if(cellAgent->unstuck_counter == unstuck_threshold+1){
                cellAgent->unstuck_counter = 0;
            }
        }


    }
}
