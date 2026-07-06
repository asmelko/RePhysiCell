#include "PhysiCell_mechanics_environment.h"
#include "../BioFVM/BioFVM_vector.h"
#include "PhysiCell_mechanics_agent.h"
#include "PhysiCell_mechanics_agent_interface.h"
#include "../core/PhysiCell_constants.h"
#include "../core/PhysiCell_cell.h"

#include "PhysiCell_mechanics_standard_models.h"

#include "../BioFVM/BioFVM_implementation.h"

#include <cmath>

using namespace BioFVM;

namespace PhysiCell {

// ============================================================================
// Initialization
// ============================================================================

mechanics_environment mech_environment; // global instance of the mechanics environment

void mechanics_environment::initialize(
    double x_start, double x_end,
    double y_start, double y_end,
    double z_start, double z_end,
    double voxel_size)
{
	initialize(x_start, x_end, y_start, y_end, z_start, z_end,
	           voxel_size, voxel_size, voxel_size);
}

void mechanics_environment::initialize(
    double x_start, double x_end,
    double y_start, double y_end,
    double z_start, double z_end,
    double dx, double dy, double dz)
{
	all_agents = (std::vector<Cell*> *) BioFVM_implementation::get_instance()->get_all_basic_agents();

	mechanics_mesh.resize(x_start, x_end, y_start, y_end, z_start, z_end, dx, dy, dz);
	agent_grid.resize(mechanics_mesh.voxels.size());
	max_cell_interactive_distance_in_voxel.resize(mechanics_mesh.voxels.size(), 0.0);
	agents_in_outer_voxels.resize(6);
}

void mechanics_environment::compute_velocities(double dt)
{
	#pragma omp parallel for 
	for( int i=0; i < (*all_agents).size(); i++ )
	{
		Mechanics_Agent* pCell = dynamic_cast<Mechanics_Agent*>(static_cast<Mechanics_Agent_PIMPL*>((*all_agents)[i])->get_mechanics_implementation()); 

		if ((*all_agents)[i]->functions.update_velocity != standard_update_cell_velocity || pCell->is_out_of_domain || !pCell->is_movable)
			continue;

		pCell->functions.add_cell_basement_membrane_interactions(dt);

		pCell->simple_pressure = 0.0; 
		pCell->neighbors.clear(); // new 1.8.0
		
		//First check the neighbors in my current voxel
		std::vector<Mechanics_Agent*>::iterator neighbor;
		std::vector<Mechanics_Agent*>::iterator end = agent_grid[pCell->get_current_mechanics_voxel_index()].end();
		for(neighbor = agent_grid[pCell->get_current_mechanics_voxel_index()].begin(); neighbor != end; ++neighbor)
		{
			pCell->add_potentials(*neighbor);
		}
		std::vector<int>::iterator neighbor_voxel_index;
		std::vector<int>::iterator neighbor_voxel_index_end = 
			mechanics_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].end();

		for( neighbor_voxel_index = 
			mechanics_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].begin();
			neighbor_voxel_index != neighbor_voxel_index_end; 
			++neighbor_voxel_index )
		{
			if(!is_neighbor_voxel(pCell, mechanics_mesh.voxels[pCell->get_current_mechanics_voxel_index()].center, mechanics_mesh.voxels[*neighbor_voxel_index].center, *neighbor_voxel_index))
				continue;
			end = agent_grid[*neighbor_voxel_index].end();
			for(neighbor = agent_grid[*neighbor_voxel_index].begin();neighbor != end; ++neighbor)
			{
				pCell->add_potentials(*neighbor);
			}
		}
		
		pCell->update_motility_vector(dt); 
		pCell->velocity += pCell->motility_data.motility_vector; 
	}
	return; 
}

void mechanics_environment::compute_spring_attachments(double dt)
{
	if( disable_automated_spring_adhesions == false )
	{
		#pragma omp parallel for 
		for( int i=0; i < (*all_agents).size(); i++ )
		{
			Mechanics_Agent* pC = dynamic_cast<Mechanics_Agent*>(static_cast<Mechanics_Agent_PIMPL*>((*all_agents)[i])->get_mechanics_implementation()); 
			Mechanics_Standard_Models().dynamic_spring_attachments(pC, dt); 
		}		
		#pragma omp parallel for 
		for( int i=0; i < (*all_agents).size(); i++ )
		{
			Mechanics_Agent* pC = dynamic_cast<Mechanics_Agent*>(static_cast<Mechanics_Agent_PIMPL*>((*all_agents)[i])->get_mechanics_implementation()); 
			if( pC->get_is_movable() )
			{
				for( int j=0; j < pC->spring_attachments.size(); j++ )
				{
					Mechanics_Agent* pC1 = dynamic_cast<Mechanics_Agent*>(pC->spring_attachments[j].first->get_mechanics_implementation()); 
					// standard_elastic_contact_function_confluent_rest_length(pC,pC->phenotype,pC1,pC1->phenotype,time_since_last_mechanics);  
					Mechanics_Standard_Models().standard_elastic_contact_function(pC,pC1,dt);  
				}
			}
		}	
	}
}

void mechanics_environment::update_positions(double dt)
{
	#pragma omp parallel for 
	for( int i=0; i < (*all_agents).size(); i++ )
	{
		Mechanics_Agent* pCell = dynamic_cast<Mechanics_Agent*>(static_cast<Mechanics_Agent_PIMPL*>((*all_agents)[i])->get_mechanics_implementation()); 
		if( pCell->get_is_out_of_domain() == false && pCell->get_is_movable() )
		{ pCell->update_position(dt); }
	}
}

void mechanics_environment::update_container()
{
	for( int i=0; i < (*all_agents).size(); i++ )
	{
		Mechanics_Agent* pCell = dynamic_cast<Mechanics_Agent*>(static_cast<Mechanics_Agent_PIMPL*>((*all_agents)[i])->get_mechanics_implementation()); 
		if(!pCell->get_is_out_of_domain() && pCell->get_is_movable())
			pCell->update_voxel_in_container();
	}
}

// ============================================================================
// Grid lifecycle helpers
// ============================================================================

void mechanics_environment::add_agent_to_voxel(
    Mechanics_Agent_PIMPL* agent, int voxel_index)
{
	add_agent_to_voxel(dynamic_cast<Mechanics_Agent*>(agent->get_mechanics_implementation()), voxel_index);
}

void mechanics_environment::add_agent_to_voxel(
    Mechanics_Agent* agent, int voxel_index)
{
	agent_grid[voxel_index].push_back(agent);
}

void mechanics_environment::remove_agent_from_voxel(
    Mechanics_Agent_PIMPL* agent, int voxel_index)
{
	remove_agent_from_voxel(dynamic_cast<Mechanics_Agent*>(agent->get_mechanics_implementation()), voxel_index);
}

void mechanics_environment::remove_agent_from_voxel(
    Mechanics_Agent* agent, int voxel_index)
{
	if (voxel_index < 0)
	{
		return; 
	}
	int delete_index = 0; 
	while( agent_grid[voxel_index][ delete_index ] != agent )
		delete_index++;
	
	// move last item to index location
    agent_grid[voxel_index][delete_index] = agent_grid[voxel_index][agent_grid[voxel_index].size()-1 ];
    // shrink the vector
    agent_grid[voxel_index].pop_back();
    
	return; 
}

int mechanics_environment::find_escaping_face_index(Mechanics_Agent* agent)
{
	if(agent->get_position()[0] <= mechanics_mesh.bounding_box[PhysiCell_constants::mesh_min_x_index])
	{ return PhysiCell_constants::mesh_lx_face_index; }
	if(agent->get_position()[0] >= mechanics_mesh.bounding_box[PhysiCell_constants::mesh_max_x_index])
	{ return PhysiCell_constants::mesh_ux_face_index; }
	if(agent->get_position()[1] <= mechanics_mesh.bounding_box[PhysiCell_constants::mesh_min_y_index])
	{ return PhysiCell_constants::mesh_ly_face_index; }
	if(agent->get_position()[1] >= mechanics_mesh.bounding_box[PhysiCell_constants::mesh_max_y_index])
	{ return PhysiCell_constants::mesh_uy_face_index; }
	if(agent->get_position()[2] <= mechanics_mesh.bounding_box[PhysiCell_constants::mesh_min_z_index])
	{ return PhysiCell_constants::mesh_lz_face_index; }
	if(agent->get_position()[2] >= mechanics_mesh.bounding_box[PhysiCell_constants::mesh_max_z_index])
	{ return PhysiCell_constants::mesh_uz_face_index; }
	return -1; 
}

void mechanics_environment::add_agent_to_outer_voxel(Mechanics_Agent_PIMPL* agent)
{
	add_agent_to_outer_voxel(dynamic_cast<Mechanics_Agent*>(agent->get_mechanics_implementation()));
}

void mechanics_environment::add_agent_to_outer_voxel(Mechanics_Agent* agent)
{
	int escaping_face= find_escaping_face_index(agent);
	agents_in_outer_voxels[escaping_face].push_back(agent);
	agent->is_out_of_domain=true;
	return; 
}


void mechanics_environment::register_agent( Mechanics_Agent_PIMPL* agent )
{
	register_agent(dynamic_cast<Mechanics_Agent*>(agent->get_mechanics_implementation()));
}

void mechanics_environment::register_agent( Mechanics_Agent* agent )
{
	agent_grid[agent->get_current_mechanics_voxel_index()].push_back(agent);
	return; 
}

void mechanics_environment::remove_agent(Mechanics_Agent_PIMPL* agent )
{
	remove_agent(dynamic_cast<Mechanics_Agent*>(agent->get_mechanics_implementation()));
}

void mechanics_environment::remove_agent(Mechanics_Agent* agent )
{
	remove_agent_from_voxel(agent, agent->get_current_mechanics_voxel_index());
	return; 
}

bool mechanics_environment::contain_any_cell(int voxel_index)
{
	// Let's replace this with clearer statements. 
	return agent_grid[voxel_index].size() > 0;
}

bool mechanics_environment::is_neighbor_voxel(Mechanics_Agent_PIMPL* pCellI, std::vector<double> my_voxel_center, std::vector<double> other_voxel_center, int other_voxel_index)
{
	return is_neighbor_voxel(dynamic_cast<Mechanics_Agent*>(pCellI->get_mechanics_implementation()), my_voxel_center, other_voxel_center, other_voxel_index);
}

bool mechanics_environment::is_neighbor_voxel(Mechanics_Agent* pCell, std::vector<double> my_voxel_center, std::vector<double> other_voxel_center, int other_voxel_index)
{
	double max_interactive_distance = pCell->mechanics_data.relative_maximum_adhesion_distance * pCell->radius_data.radius 
		+ max_cell_interactive_distance_in_voxel[other_voxel_index];
	
	int comparing_dimension = -1, comparing_dimension2 = -1;
	if(my_voxel_center[0] == other_voxel_center[0] && my_voxel_center[1] == other_voxel_center[1])
	{
		comparing_dimension = 2;
	}
	else if(my_voxel_center[0] == other_voxel_center[0] && my_voxel_center[2] == other_voxel_center[2])
	{
		comparing_dimension = 1;
	}
	else if(my_voxel_center[1] == other_voxel_center[1] && my_voxel_center[2] == other_voxel_center[2])
	{
		comparing_dimension = 0;
	}
	
	if(comparing_dimension != -1) 
	{ //then it is an immediate neighbor (through side faces)
		double surface_coord= 0.5*(my_voxel_center[comparing_dimension] + other_voxel_center[comparing_dimension]);
		if(std::fabs(pCell->get_position()[comparing_dimension] - surface_coord) > max_interactive_distance)
		{ return false; }
		return true;
	}
	comparing_dimension=-1;
	
	if(my_voxel_center[0] == other_voxel_center[0])
	{
		comparing_dimension = 1; comparing_dimension2 = 2;
	}
	else if(my_voxel_center[1] == other_voxel_center[1])
	{
		comparing_dimension=0; comparing_dimension2 = 2;
	}
	else if(my_voxel_center[2] == other_voxel_center[2])
	{
		comparing_dimension = 0; comparing_dimension2=1;
	}
	if(comparing_dimension != -1)
	{
		double line_coord1= 0.5*(my_voxel_center[comparing_dimension] + other_voxel_center[comparing_dimension]);
		double line_coord2= 0.5*(my_voxel_center[comparing_dimension2] + other_voxel_center[comparing_dimension2]);
		double distance_squared= std::pow( pCell->get_position()[comparing_dimension] - line_coord1,2)+ std::pow( pCell->get_position()[comparing_dimension2] - line_coord2,2);
		if(distance_squared > max_interactive_distance * max_interactive_distance)
		{ return false; }
		return true;
	}
	std::vector<double> corner_point= 0.5*(my_voxel_center+other_voxel_center);
	double distance_squared= (corner_point[0]-pCell->get_position()[0])*(corner_point[0]-pCell->get_position()[0])
		+(corner_point[1]-pCell->get_position()[1])*(corner_point[1]-pCell->get_position()[1]) 
		+(corner_point[2]-pCell->get_position()[2]) * (corner_point[2]-pCell->get_position()[2]);
	if(distance_squared > max_interactive_distance * max_interactive_distance)
	{ return false; }
	return true;
}


void mechanics_environment::attach_cells_as_spring( Mechanics_Agent_PIMPL* pCell_1I, Mechanics_Agent_PIMPL* pCell_2I, bool attacking_spring )
{
	attach_cells_as_spring(
		dynamic_cast<Mechanics_Agent*>(pCell_1I->get_mechanics_implementation()),
		dynamic_cast<Mechanics_Agent*>(pCell_2I->get_mechanics_implementation()),
		attacking_spring
	);
}

void mechanics_environment::attach_cells_as_spring( Mechanics_Agent* pCell_1, Mechanics_Agent* pCell_2, bool attacking_spring )
{
	pCell_1->attach_cell_as_spring( pCell_2, attacking_spring );
	pCell_2->attach_cell_as_spring( pCell_1, attacking_spring );
	return; 
}

void mechanics_environment::detach_cells_as_spring( Mechanics_Agent_PIMPL* pCell_1I, Mechanics_Agent_PIMPL* pCell_2I )
{
	detach_cells_as_spring(
		dynamic_cast<Mechanics_Agent*>(pCell_1I->get_mechanics_implementation()),
		dynamic_cast<Mechanics_Agent*>(pCell_2I->get_mechanics_implementation())
	);
}

void mechanics_environment::detach_cells_as_spring( Mechanics_Agent* pCell_1, Mechanics_Agent* pCell_2 )
{
	pCell_1->detach_cell_as_spring( pCell_2 );
	pCell_2->detach_cell_as_spring( pCell_1 );
	return; 
}

void mechanics_environment::standard_add_basement_membrane_interactions( Mechanics_Agent_PIMPL* pCellI, double dt )
{
	models.standard_add_basement_membrane_interactions(pCellI, dt);
}

void mechanics_environment::standard_domain_edge_avoidance_interactions( Mechanics_Agent_PIMPL* pCellI, double dt )
{
	models.standard_domain_edge_avoidance_interactions(pCellI, dt);
}

double mechanics_environment::distance_to_domain_edge(Mechanics_Agent_PIMPL* pCellI)
{
	return models.distance_to_domain_edge(pCellI);
}

void mechanics_environment::dynamic_spring_attachments( Mechanics_Agent_PIMPL* pCellI , double dt )
{
	models.dynamic_spring_attachments(pCellI, dt);
}

void mechanics_environment::dynamic_spring_attachments( Mechanics_Agent* pCell , double dt )
{
	models.dynamic_spring_attachments(pCell, dt);
}

void mechanics_environment::standard_elastic_contact_function( Mechanics_Agent_PIMPL* pC1I, Mechanics_Agent_PIMPL* pC2I, double dt )
{
	models.standard_elastic_contact_function(pC1I, pC2I, dt);
}

void mechanics_environment::standard_elastic_contact_function( Mechanics_Agent* pC1, Mechanics_Agent* pC2, double dt )
{
	models.standard_elastic_contact_function(pC1, pC2, dt);
}

void mechanics_environment::standard_elastic_contact_function_confluent_rest_length( Mechanics_Agent_PIMPL* pC1I, Mechanics_Agent_PIMPL* pC2I, double dt )
{
	models.standard_elastic_contact_function_confluent_rest_length(pC1I, pC2I, dt);
}

void mechanics_environment::chemotaxis_function( Mechanics_Agent_PIMPL* pCellI, double dt )
{
	models.chemotaxis_function(pCellI, dt);
}

void mechanics_environment::advanced_chemotaxis_function_normalized( Mechanics_Agent_PIMPL* pCellI, double dt )
{
	models.advanced_chemotaxis_function_normalized(pCellI, dt);
}

void mechanics_environment::advanced_chemotaxis_function( Mechanics_Agent_PIMPL* pCellI, double dt )
{
	models.advanced_chemotaxis_function(pCellI, dt);
}

std::vector<Mechanics_Agent_PIMPL*> mechanics_environment::agents_in_my_container( Mechanics_Agent_PIMPL* pCellI )
{
	return agents_in_my_container(dynamic_cast<Mechanics_Agent*>(pCellI->get_mechanics_implementation()));
}

std::vector<Mechanics_Agent_PIMPL*> mechanics_environment::agents_in_my_container( Mechanics_Agent* pCell )
{
	std::vector<Mechanics_Agent_PIMPL*> result;
	for( Mechanics_Agent* a : agent_grid[pCell->get_current_mechanics_voxel_index()] )
	{ result.push_back(a->pOwner); }
	return result;
}

std::vector<Mechanics_Agent_PIMPL*> mechanics_environment::nearby_agents( Mechanics_Agent_PIMPL* pCellI )
{
	return nearby_agents(dynamic_cast<Mechanics_Agent*>(pCellI->get_mechanics_implementation()));
}

std::vector<Mechanics_Agent_PIMPL*> mechanics_environment::nearby_agents( Mechanics_Agent* pCell )
{
	std::vector<Mechanics_Agent_PIMPL*> result;

	std::vector<Mechanics_Agent*>::iterator neighbor;
	std::vector<Mechanics_Agent*>::iterator end =
		agent_grid[pCell->get_current_mechanics_voxel_index()].end();
	for( neighbor = agent_grid[pCell->get_current_mechanics_voxel_index()].begin(); neighbor != end; ++neighbor)
	{ result.push_back((*neighbor)->pOwner); }

	std::vector<int>::iterator neighbor_voxel_index;
	std::vector<int>::iterator neighbor_voxel_index_end =
		mechanics_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].end();

	for( neighbor_voxel_index =
		mechanics_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].begin();
		neighbor_voxel_index != neighbor_voxel_index_end;
		++neighbor_voxel_index )
	{
		if(!is_neighbor_voxel(pCell, mechanics_mesh.voxels[pCell->get_current_mechanics_voxel_index()].center, mechanics_mesh.voxels[*neighbor_voxel_index].center, *neighbor_voxel_index))
			continue;
		end = agent_grid[*neighbor_voxel_index].end();
		for(neighbor = agent_grid[*neighbor_voxel_index].begin(); neighbor != end; ++neighbor)
		{ result.push_back((*neighbor)->pOwner); }
	}

	return result;
}

std::vector<Mechanics_Agent_PIMPL*> mechanics_environment::nearby_interacting_agents( Mechanics_Agent_PIMPL* pCellI )
{
	return nearby_interacting_agents(dynamic_cast<Mechanics_Agent*>(pCellI->get_mechanics_implementation()));
}

std::vector<Mechanics_Agent_PIMPL*> mechanics_environment::nearby_interacting_agents( Mechanics_Agent* pCell )
{
	std::vector<Mechanics_Agent_PIMPL*> result;

	std::vector<Mechanics_Agent*>::iterator neighbor;
	std::vector<Mechanics_Agent*>::iterator end = agent_grid[pCell->get_current_mechanics_voxel_index()].end();
	for( neighbor = agent_grid[pCell->get_current_mechanics_voxel_index()].begin(); neighbor != end; ++neighbor)
	{
		std::vector<double> displacement = (*neighbor)->get_position() - pCell->get_position();
		double distance = norm( displacement );
		if( distance <= pCell->mechanics_data.relative_maximum_adhesion_distance * pCell->radius_data.radius
			+ (*neighbor)->mechanics_data.relative_maximum_adhesion_distance * (*neighbor)->radius_data.radius
			&& (*neighbor) != pCell )
		{ result.push_back((*neighbor)->pOwner); }
	}

	std::vector<int>::iterator neighbor_voxel_index;
	std::vector<int>::iterator neighbor_voxel_index_end =
		mechanics_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].end();

	for( neighbor_voxel_index =
		mechanics_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].begin();
		neighbor_voxel_index != neighbor_voxel_index_end;
		++neighbor_voxel_index )
	{
		if(!is_neighbor_voxel(pCell, mechanics_mesh.voxels[pCell->get_current_mechanics_voxel_index()].center, mechanics_mesh.voxels[*neighbor_voxel_index].center, *neighbor_voxel_index))
			continue;
		end = agent_grid[*neighbor_voxel_index].end();
		for(neighbor = agent_grid[*neighbor_voxel_index].begin(); neighbor != end; ++neighbor)
		{
			std::vector<double> displacement = (*neighbor)->get_position() - pCell->get_position();
			double distance = norm( displacement );
			if( distance <= pCell->mechanics_data.relative_maximum_adhesion_distance * pCell->radius_data.radius
				+ (*neighbor)->mechanics_data.relative_maximum_adhesion_distance * (*neighbor)->radius_data.radius
				&& (*neighbor) != pCell )
			{ result.push_back((*neighbor)->pOwner); }
		}
	}

	return result;
}




} // namespace PhysiCell
