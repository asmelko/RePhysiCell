#include "PhysiCell_mechanics_agent.h"

#include "../BioFVM/BioFVM_vector.h"
#include "../core/PhysiCell_utilities.h"
#include "PhysiCell_mechanics_agent_interface.h"
#include "PhysiCell_mechanics_environment.h"
#include "../BioFVM/BioFVM_microenvironment_interface.h"
#include "../core/PhysiCell_cell.h"

#include "../BioFVM/BioFVM_implementation.h"
#include "PhysiCell_mechanics_functions.h"

using namespace BioFVM;

namespace PhysiCell {

Mechanics_Agent::Mechanics_Agent(Cell* pCell) : Mechanics_Agent::Mechanics_Agent(BioFVM::BioFVM_implementation::get_instance()->create_basic_agent(), pCell) {}

Mechanics_Agent::Mechanics_Agent(BioFVM::Basic_Agent_Interface* pBasicAgent, Cell* pCell) : Basic_Agent_PIMPL(pBasicAgent, false), functions(pCell) 
{
	pOwner = static_cast<Mechanics_Agent_PIMPL*>(pCell);
	velocity.resize(3, 0.0);
	previous_velocity.resize(3, 0.0);
	springs.clear();
	neighbors.clear();
	neighbors2.clear();

	mechanics_data.sync_to_cell_definitions();
	motility_data.sync_to_current_microenvironment();
}

void Mechanics_Agent::update_motility_vector( double dt_ )
{
	if( motility_data.is_motile == false )
	{
		motility_data.motility_vector.assign( 3, 0.0 ); 
		return; 
	}
	
	if( UniformRandom() < dt_ / motility_data.persistence_time || motility_data.persistence_time < dt_ )
	{
		/*
		// choose a uniformly random unit vector 
		double temp_angle = 6.28318530717959*UniformRandom();
		double temp_phi = 3.1415926535897932384626433832795*UniformRandom();
		
		double sin_phi = sin(temp_phi);
		double cos_phi = cos(temp_phi);
		
		if( phenotype.motility.restrict_to_2D == true )
		{ 
			sin_phi = 1.0; 
			cos_phi = 0.0;
		}
		
		std::vector<double> randvec; 
		randvec.resize(3,sin_phi); 
		
		randvec[0] *= cos( temp_angle ); // cos(theta)*sin(phi)
		randvec[1] *= sin( temp_angle ); // sin(theta)*sin(phi)
		randvec[2] = cos_phi; //  cos(phi)
		*/
		std::vector<double> randvec(3,0.0);
		if( motility_data.restrict_to_2D == true )
		{ randvec = UniformOnUnitCircle(); }
		else
		{ randvec = UniformOnUnitSphere(); }

		// if the update_bias_vector function is set, use it  
		functions.update_migration_bias( dt_ ); 
		
		motility_data.motility_vector = motility_data.migration_bias_direction; // motiltiy = bias_vector
		motility_data.motility_vector *= motility_data.migration_bias; // motility = bias*bias_vector 
		
		double one_minus_bias = 1.0 - motility_data.migration_bias; 
		
		axpy( &motility_data.motility_vector, one_minus_bias, randvec ); // motility = (1-bias)*randvec + bias*bias_vector
		
		normalize( &motility_data.motility_vector ); 
		
		motility_data.motility_vector *= motility_data.migration_speed; 
	}	
	return; 
} 

bool Mechanics_Agent::assign_position(const std::vector<double>& new_position)
{
	return assign_position(new_position[0], new_position[1], new_position[2]);
}

void Mechanics_Agent::set_previous_velocity(double xV, double yV, double zV)
{
	get_previous_velocity()[0] = xV;
	get_previous_velocity()[1] = yV;
	get_previous_velocity()[2] = zV;

	return; 
}

bool Mechanics_Agent::assign_position(double x, double y, double z)
{
	get_position_internal()[0] = x;
	get_position_internal()[1] = y;
	if ( !get_microenvironment_i()->simulate_2D() )
	{ get_position_internal()[2] = z; }
	
	// update microenvironment current voxel index
	update_voxel_index();
	// update current_mechanics_voxel_index
	current_mechanics_voxel_index= mech_environment.mechanics_mesh.nearest_voxel_index( get_position() );

    // Since it is most likely our first position, we update the max_cell_interactive_distance_in_voxel
	// which was not initialized at cell creation
	if( mech_environment.max_cell_interactive_distance_in_voxel[get_current_mechanics_voxel_index()] < 
		radius * mechanics_data.relative_maximum_adhesion_distance )
	{
		// get_container()->max_cell_interactive_distance_in_voxel[get_current_mechanics_voxel_index()]= phenotype.geometry.radius*parameters.max_interaction_distance_factor;
		mech_environment.max_cell_interactive_distance_in_voxel[get_current_mechanics_voxel_index()] = radius
			* mechanics_data.relative_maximum_adhesion_distance;
	}

	mech_environment.register_agent(this);
	
	if( !mech_environment.mechanics_mesh.is_position_valid(x,y,z) )
	{	
		is_out_of_domain = true; 
		set_is_active(false); 
		is_movable = false; 
		
		return false;
	}
	
	return true;
}

void Mechanics_Agent::update_position( double dt )
{
	// BioFVM Basic_Agent::update_position(dt) returns without doing anything. 
	// So we remove this to avoid any future surprises. 
	// 
	// Basic_Agent::update_position(dt);
		
	// use Adams-Bashforth 
	static double d1; 
	static double d2; 
	static bool constants_defined = false; 
	if( constants_defined == false )
	{
		d1 = dt; 
		d1 *= 1.5; 
		d2 = dt; 
		d2 *= -0.5; 
		constants_defined = true; 
	}
	
	// new AUgust 2017
	if( get_microenvironment_i()->simulate_2D() == true )
	{ velocity[2] = 0.0; }
	
	int dims = get_microenvironment_i()->simulate_2D() ? 2 : 3;

	// std::vector<double> old_position = position;
	for ( int i = 0 ; i < dims ; i++ )
	{
		get_position_internal()[i] += 
			( d1 * get_velocity()[i] + d2 * get_previous_velocity()[i] );
	}
	// overwrite previous_velocity for future use 
	// if(sqrt(dist(old_position, position))>3* phenotype.geometry.radius)
		// std::cout<<sqrt(dist(old_position, position))<<"old_position: "<<old_position<<", new position: "<< position<<", velocity: "<<velocity<<", previous_velocity: "<< previous_velocity<<std::endl;
	
		
	previous_velocity = velocity; 
	
	velocity[0]=0; velocity[1]=0; velocity[2]=0;
	if(mech_environment.mechanics_mesh.is_position_valid(get_position()[0],get_position()[1],get_position()[2]))
	{
		updated_current_mechanics_voxel_index=mech_environment.mechanics_mesh.nearest_voxel_index( get_position() );
	}
	else
	{
		updated_current_mechanics_voxel_index=-1;
		
		is_out_of_domain = true; 
		set_is_active(false); 
		is_movable = false; 
	}
	return; 
}

int Mechanics_Agent::get_current_mechanics_voxel_index()
{
	return current_mechanics_voxel_index;
}

void Mechanics_Agent::update_voxel_in_container()
{
	// call the method from BioFVM_basic_agent to update microenvironment's voxel index
	update_voxel_index();
	// int temp_current_voxel_index;
	// Check to see if we need to remove agents that are pushed out of boundary
	// if(!get_container()->underlying_mesh.is_position_valid(position[0],position[1],position[2]))	
		
	if(updated_current_mechanics_voxel_index==-1)// updated_current_mechanics_voxel_index is updated in update_position
	{
		// check if this agent has a valid voxel index, if so, remove it from previous voxel
		if( get_current_mechanics_voxel_index() >= 0)
		{
			{mech_environment.remove_agent_from_voxel(this, get_current_mechanics_voxel_index());}
		}
		{mech_environment.add_agent_to_outer_voxel(this);}
		// std::cout<<"cell out of boundary..."<< __LINE__<<" "<<ID<<std::endl;
		current_mechanics_voxel_index=-1;
		is_out_of_domain=true;
		set_is_active(false);
		return;
	}
	
	// temp_current_voxel_index= get_current_mechanics_voxel_index();
	// updated_current_mechanics_voxel_index=get_container()->underlying_mesh.nearest_voxel_index( position );
	
	// update mesh indices (if needed)
	if(updated_current_mechanics_voxel_index!= get_current_mechanics_voxel_index())
	{
		{
			mech_environment.remove_agent_from_voxel(this, get_current_mechanics_voxel_index());
			mech_environment.add_agent_to_voxel(this, updated_current_mechanics_voxel_index);
		}
		current_mechanics_voxel_index=updated_current_mechanics_voxel_index;
	}
	
	return; 
}

void Mechanics_Agent::add_potentials(Mechanics_Agent* other_agent)
{
	// if( this->ID == other_agent->ID )
	if( this == other_agent )
	{ return; }
/*
	// new April 2022: don't interact with cells with 0 volume 
	// does not seem to really help 
	if( other_agent->phenotype.volume.total < 1e-15 )
	{ std::cout << "zero size cell in mechanics!" << std::endl; return; }
*/
	// 12 uniform neighbors at a close packing distance, after dividing out all constants
	static double simple_pressure_scale = 0.027288820670331; // 12 * (1 - sqrt(pi/(2*sqrt(3))))^2 
	// 9.820170012151277; // 12 * ( 1 - sqrt(2*pi/sqrt(3)))^2

	double distance = 0; 
	for( int i = 0 ; i < 3 ; i++ ) 
	{ 
		displacement[i] = get_position()[i] - (*other_agent).get_position()[i]; 
		distance += displacement[i] * displacement[i]; 
	}
	// Make sure that the distance is not zero
	
	distance = std::max(sqrt(distance), 0.00001); 
	
	//Repulsive
	double R = radius+ (*other_agent).radius; 
	
	// double RN = phenotype.geometry.nuclear_radius + (*other_agent).phenotype.geometry.nuclear_radius;	
	double temp_r, c;
	if( distance > R ) 
	{
		temp_r=0;
	}
	// else if( distance < RN ) 
	// {
		// double M = 1.0; 
		// c = 1.0 - RN/R; 
		// c *= c; 
		// c -= M; 
		// temp_r = ( c*distance/RN  + M  ); 
	// }
	else
	{
		// temp_r = 1 - distance/R;
		temp_r = -distance; // -d
		temp_r /= R; // -d/R
		temp_r += 1.0; // 1-d/R
		temp_r *= temp_r; // (1-d/R)^2 
		
		// add the relative pressure contribution 
		simple_pressure += ( temp_r / simple_pressure_scale ); // New July 2017 
	}
	
	// August 2017 - back to the original if both have same coefficient 

	double effective_repulsion = sqrt( mechanics_data.cell_cell_repulsion_strength * other_agent->mechanics_data.cell_cell_repulsion_strength ); 
	temp_r *= effective_repulsion; 
	
	// temp_r *= phenotype.mechanics.cell_cell_repulsion_strength; // original 
	//////////////////////////////////////////////////////////////////
	
	// Adhesive
	//double max_interactive_distance = parameters.max_interaction_distance_factor * phenotype.geometry.radius + 
	//	(*other_agent).parameters.max_interaction_distance_factor * (*other_agent).phenotype.geometry.radius;
		
	double max_interactive_distance = mechanics_data.relative_maximum_adhesion_distance * radius + 
		(*other_agent).mechanics_data.relative_maximum_adhesion_distance * (*other_agent).radius;
		
	if(distance < max_interactive_distance ) 
	{	
		// double temp_a = 1 - distance/max_interactive_distance; 
		double temp_a = -distance; // -d
		temp_a /= max_interactive_distance; // -d/S
		temp_a += 1.0; // 1 - d/S 
		temp_a *= temp_a; // (1-d/S)^2 
		// temp_a *= phenotype.mechanics.cell_cell_adhesion_strength; // original 
		
		// August 2017 - back to the original if both have same coefficient 
		// May 2022 - back to oriinal if both affinities are 1
		int ii = find_cell_definition_index( this->get_type() ); 
		int jj = find_cell_definition_index( other_agent->get_type() ); 

		double adhesion_ii = mechanics_data.cell_cell_adhesion_strength * mechanics_data.cell_adhesion_affinities[jj]; 
		double adhesion_jj = other_agent->mechanics_data.cell_cell_adhesion_strength * other_agent->mechanics_data.cell_adhesion_affinities[ii]; 

		// double effective_adhesion = sqrt( phenotype.mechanics.cell_cell_adhesion_strength() * other_agent->phenotype.mechanics.cell_cell_adhesion_strength() ); 
		double effective_adhesion = sqrt( adhesion_ii*adhesion_jj ); 
		temp_a *= effective_adhesion; 
		
		temp_r -= temp_a;

		neighbors.push_back(other_agent); // move here in 1.10.2 so non-adhesive cells also added. 
	}
	/////////////////////////////////////////////////////////////////
	if( fabs(temp_r) < 1e-16 )
	{ return; }
	temp_r /= distance;
	// for( int i = 0 ; i < 3 ; i++ ) 
	// {
	//	velocity[i] += displacement[i] * temp_r; 
	// }
	axpy( &(velocity) , temp_r , displacement ); 
	
	
	// state.neighbors.push_back(other_agent); // new 1.8.0
	
	return;
}


int Mechanics_Agent::number_of_attached_cells( void )
{ return attached_cells.size(); } 

void Mechanics_Agent::attach_cell( Mechanics_Agent* pAddMe )
{
	#pragma omp critical
	{
		bool already_attached = false; 
		for( int i=0 ; i < attached_cells.size() ; i++ )
		{
			if( attached_cells[i] == pAddMe )
			{ already_attached = true; }
		}
		if( already_attached == false )
		{ attached_cells.push_back( pAddMe ); }
	}
	// pAddMe->attach_cell( this ); 
	return; 
}

void Mechanics_Agent::attach_cell_as_spring( Mechanics_Agent* pAddMe, bool attacking_spring )
{
	#pragma omp critical
	{
		bool already_attached = false; 
		for( int i=0 ; i < spring_attachments.size() ; i++ )
		{
			if( spring_attachments[i].first == pAddMe )
			{ already_attached = true; }
		}
		if( already_attached == false )
		{ spring_attachments.emplace_back( pAddMe, attacking_spring ); }
	}
	// pAddMe->attach_cell( this ); 
	return; 
}

void Mechanics_Agent::detach_cell( Mechanics_Agent* pRemoveMe )
{
	#pragma omp critical
	{
		bool found = false; 
		int i = 0; 
		while( !found && i < attached_cells.size() )
		{
			// if pRemoveMe is in the cell's list, remove it
			if( attached_cells[i] == pRemoveMe )
			{
				int n = attached_cells.size(); 
				// copy last entry to current position 
				attached_cells[i] = attached_cells[n-1]; 
				// shrink by one 
				attached_cells.pop_back(); 
				found = true; 
			}
			i++; 
		}
	}
	return; 
}

void Mechanics_Agent::detach_cell_as_spring( Mechanics_Agent* pRemoveMe )
{
	#pragma omp critical
	{
		bool found = false; 
		int i = 0; 
		while( !found && i < spring_attachments.size() )
		{
			// if pRemoveMe is in the cell's list, remove it
			if( spring_attachments[i].first == pRemoveMe )
			{
				int n = spring_attachments.size(); 
				// copy last entry to current position 
				spring_attachments[i] = spring_attachments[n-1]; 
				// shrink by one 
				spring_attachments.pop_back(); 
				found = true; 
			}
			i++; 
		}
	}
	return; 
}

void Mechanics_Agent::remove_self_from_all_neighbors( void )
{
	Mechanics_Agent* pCell = this; 
	// go through all neighbors (pN) of this (pC)

	for( int j = 0 ; j < pCell->neighbors.size(); j++ )
	{
	 	Mechanics_Agent* pN = pCell->neighbors[j]; 

		// for each pN, remove pC from list of neighbors 
			// find pC in neighbors 


			auto SearchResult = std::find( 
				pN->neighbors.begin(),pN->neighbors.end(),pCell );  		

			// if pC is indeed found, remove it  
			// erase pC from neighbors 
			if( SearchResult != pN->neighbors.end() )
			{
				// if the target is found, set the appropriate rate 
				pN->neighbors.erase( SearchResult ); 
			}
			else
			{ /* future error message */  }
	}

	return; 
}

void Mechanics_Agent::remove_all_attached_cells( void )
{
	{
		// remove self from any attached cell's list. 
		for( int i = 0; i < attached_cells.size() ; i++ )
		{
			attached_cells[i]->detach_cell( this ); 
		}
		// clear my list 
		attached_cells.clear(); 
	}
	return; 
}

void Mechanics_Agent::remove_all_spring_attachments( void )
{
	{
		// remove self from any attached cell's list. 
		for( int i = 0; i < spring_attachments.size() ; i++ )
		{
			spring_attachments[i].first->detach_cell_as_spring( this ); 
		}
		// clear my list 
		spring_attachments.clear(); 
	}
	return; 
}

}
