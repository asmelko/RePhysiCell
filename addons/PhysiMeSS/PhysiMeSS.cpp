#include "PhysiMeSS.h"
#include "PhysiMeSS_agent.h"

#include <algorithm>
#include <iterator> 

#include "../../BioFVM/BioFVM_vector.h"


static double last_update_time = -mechanics_dt;

void remove_physimess_out_of_bounds_fibres()
{
    for (auto* cell : *all_cells) {
        if (isFibre(cell) && static_cast<PhysiMeSS_Fibre*>(cell)->get_physimess_agent()->fail_count >= 10)
        {
            // std::cout << "I failed to place " << cell->type_name << " " 
            //           << cell->get_ID() << " in the domain - I am deleting agent " 
            //           << std::endl;
            delete_cell(cell);
        }
    }
}

void physimess_update_cell_velocity( Cell* pCell, Phenotype& phenotype, double dt)
{
    // This function is a placeholder for the standard update of cell velocity.
    // The implementation of this function is in PhysiMeSS, and is called by the main PhysiCell update_velocity() function.
}

void physimess_mechanics( double dt ) 
{    
    if ( PhysiCell_globals.current_time >= last_update_time + dt - 0.5 * diffusion_dt)
    {
        last_update_time = PhysiCell_globals.current_time;
        
        //#pragma omp parallel for
        // This is not parallel because we are modifying the agend grid
        for( int i=0; i < (*all_cells).size(); i++ )
        {
            Cell* pC = (*all_cells)[i];
            PhysiMeSS_Agent* agent = dynamic_cast<PhysiMeSS_Agent*>(pC->get_mechanics_implementation());
            agent->physimess_voxels.clear();
            if( !pC->get_is_out_of_domain() )
            {
                agent->register_fibre_voxels();
            }
        }

        #pragma omp parallel for
        for( int i=0; i < (*all_cells).size(); i++ )
        {
            Cell* pC = (*all_cells)[i];
            PhysiMeSS_Agent* agent = dynamic_cast<PhysiMeSS_Agent*>(pC->get_mechanics_implementation());
            agent->physimess_neighbors.clear();
            if (isFibre(pC)) {
                dynamic_cast<PhysiMeSS_FibreAgent*>(agent)->fibres_crosslinkers.clear();
            }
            if( !pC->get_is_out_of_domain() )
            {
                agent->find_agent_neighbors();
            }
        }

        // #pragma omp parallel for
        // This is not parallel because we are modifying the agend grid
        for( int i=0; i < (*all_cells).size(); i++ )
        {
            Cell* pC = (*all_cells)[i];
            PhysiMeSS_Agent* agent = dynamic_cast<PhysiMeSS_Agent*>(pC->get_mechanics_implementation());
            if( !pC->get_is_out_of_domain() )
            {
                agent->deregister_fibre_voxels();
            }
        }
        
        // determine and add crosslinks
        #pragma omp parallel for
        for( int i=0; i < (*all_cells).size(); i++ )
        {
            Cell* pC = (*all_cells)[i];
            if (isFibre(pC)) {
                dynamic_cast<PhysiMeSS_FibreAgent*>(pC->get_mechanics_implementation())->add_crosslinks();
            }
        }
    }
}


void fibre_agent_SVG(std::ofstream& os, PhysiCell::Cell* pC, double z_slice, std::vector<std::string> (*cell_coloring_function)(Cell*), double X_lower, double Y_lower) {

	// place a rod if it's a fibre (note fibre already renamed here)
	if (isFibre(pC) ){
    
        PhysiMeSS_Fibre* pFibre = static_cast<PhysiMeSS_Fibre*>(pC);
        PhysiMeSS_FibreAgent* fibre_agent = pFibre->get_physimess_agent();
        
		int crosslinks = fibre_agent->X_crosslink_count;
        if (crosslinks >= 3){
			// if fibre has cross-links different colour than if not
			Write_SVG_line(os, (pC->position)[0] - (fibre_agent->mLength) * (pC->state.orientation)[0] - X_lower,
							(pC->position)[1] - (fibre_agent->mLength) * (pC->state.orientation)[1] - Y_lower,
							(pC->position)[0] + (fibre_agent->mLength) * (pC->state.orientation)[0] - X_lower,
							(pC->position)[1] + (fibre_agent->mLength) * (pC->state.orientation)[1] - Y_lower,
							4.0, "darkblue");
		}
		else if (crosslinks == 2){
			// if fibre has cross-links different colour than if not
			Write_SVG_line(os, (pC->position)[0] - (fibre_agent->mLength) * (pC->state.orientation)[0] - X_lower,
							(pC->position)[1] - (fibre_agent->mLength) * (pC->state.orientation)[1] - Y_lower,
							(pC->position)[0] + (fibre_agent->mLength) * (pC->state.orientation)[0] - X_lower,
							(pC->position)[1] + (fibre_agent->mLength) * (pC->state.orientation)[1] - Y_lower,
							4.0, "blue");
		}
		else if (crosslinks == 1){
			// if fibre has cross-links different colour than if not
			Write_SVG_line(os, (pC->position)[0] - (fibre_agent->mLength) * (pC->state.orientation)[0] - X_lower,
							(pC->position)[1] - (fibre_agent->mLength) * (pC->state.orientation)[1] - Y_lower,
							(pC->position)[0] + (fibre_agent->mLength) * (pC->state.orientation)[0] - X_lower,
							(pC->position)[1] + (fibre_agent->mLength) * (pC->state.orientation)[1] - Y_lower,
							4.0, "steelblue");
		}
		else {
    		Write_SVG_line(os, (pC->position)[0] - (fibre_agent->mLength) * (pC->state.orientation)[0] - X_lower,
						(pC->position)[1] - (fibre_agent->mLength) * (pC->state.orientation)[1] - Y_lower,
						(pC->position)[0] + (fibre_agent->mLength) * (pC->state.orientation)[0] - X_lower,
						(pC->position)[1] + (fibre_agent->mLength) * (pC->state.orientation)[1] - Y_lower,
						4.0, "lightskyblue");
		}

	}
	else{
        standard_agent_SVG(os, pC, z_slice, cell_coloring_function, X_lower, Y_lower);
	}
}

void fibre_agent_legend(std::ofstream& os, Cell_Definition* cell_definition, double& cursor_x, double& cursor_y, std::vector<std::string> (*cell_coloring_function)(Cell*), double temp_cell_radius) {
	
	// switch to the cell type 
	Cell C; 
	C.convert_to_cell_definition( *cell_definition );
	
	// get the colors using the current coloring function 
	std::vector<std::string> colors = cell_coloring_function(&C); 	

	// place the label 
	// place a rod if it's a fibre (note fibre not yet renamed)
	if (isFibre(&C)) {
		//Write_SVG_fibre(os, cursor_x, cursor_y , 0.5*temp_cell_radius , 1.0 , colors[1] , colors[0] );
		Write_SVG_line(os, cursor_x, cursor_y-20.0 , cursor_x , cursor_y+20.0 , 4.0 , "lightskyblue" );
	}
	else {
		standard_agent_legend(os, cell_definition, cursor_x, cursor_y, cell_coloring_function, temp_cell_radius);
	}
}