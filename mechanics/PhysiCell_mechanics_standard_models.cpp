#include "PhysiCell_mechanics_standard_models.h"
#include "PhysiCell_mechanics_agent.h"
#include "../core/PhysiCell_utilities.h"
#include "../core/PhysiCell_cell.h"

#include "../BioFVM/BioFVM_vector.h"
#include "../BioFVM/BioFVM_microenvironment_interface.h"
#include "PhysiCell_mechanics_agent.h"
#include "PhysiCell_mechanics_agent_interface.h"
#include "PhysiCell_mechanics_environment.h"

using namespace BioFVM;

namespace PhysiCell {

void Mechanics_Standard_Models::standard_add_basement_membrane_interactions( Mechanics_Agent_PIMPL* pCellI, double dt )
{
	Mechanics_Agent* pCell = dynamic_cast<Mechanics_Agent*>(pCellI->get_mechanics_implementation());
	
	double max_interactive_distance = pCell->mechanics_data.relative_maximum_adhesion_distance * pCell->radius_data.radius;
	double distance = pCell->functions.calculate_distance_to_membrane(dt); 
	//Note that the distance_to_membrane function must set displacement values (as a normal vector)
		
	double temp_a=0;
	// Adhesion to basement membrane
	if(distance< max_interactive_distance)
	{
		temp_a= (1- distance/max_interactive_distance);
		temp_a*=temp_a;
		temp_a*=-pCell->mechanics_data.cell_BM_adhesion_strength;
	}
	// Repulsion from basement membrane
	double temp_r = 0;
	if(distance < pCell->radius_data.radius)
	{
		temp_r = (1- distance/pCell->radius_data.radius);
		temp_r *= temp_r;
		temp_r *= pCell->mechanics_data.cell_BM_repulsion_strength;
	}
	temp_r += temp_a;
	if( fabs( temp_r ) < 1e-16 )
	{ return; }
	
	axpy( &( pCell->velocity ) , temp_r , pCell->displacement ); 
	return;	
}

void Mechanics_Standard_Models::standard_domain_edge_avoidance_interactions( Mechanics_Agent_PIMPL* pCellI, double dt )
{
	Mechanics_Agent* pCell = dynamic_cast<Mechanics_Agent*>(pCellI->get_mechanics_implementation());
	pCell->mechanics_data.cell_BM_repulsion_strength = 100;  
		
	double distance = pCell->functions.calculate_distance_to_membrane(dt); 
	//Note that the distance_to_membrane function must set displacement values (as a normal vector)
		
	// Repulsion from basement membrane
	double temp_r = 0;
	if(distance < pCell->radius_data.radius)
	{
		temp_r = (1- distance/pCell->radius_data.radius);
		temp_r *= temp_r;
		temp_r *= pCell->mechanics_data.cell_BM_repulsion_strength;
	}
	if( fabs( temp_r ) < 1e-16 )
	{ return; }
	
	axpy( &( pCell->velocity ) , temp_r , pCell->displacement ); 
	return;
}

double Mechanics_Standard_Models::distance_to_domain_edge(Mechanics_Agent_PIMPL* pCellI)
{
	Mechanics_Agent* pCell = dynamic_cast<Mechanics_Agent*>(pCellI->get_mechanics_implementation());
	
	static double tolerance = 1e-7;
	static double one_over_sqrt_2 = 0.70710678118;
	static double one_over_sqrt_3 = 0.57735026919; 
	
		
	double min_distance = 9e99; 
	int nearest_boundary = -1; 
	
	// check against xL and xU
	double temp_distance = pCell->get_position()[0] - get_microenvironment_i()->get_mesh().bounding_box[0]; 
	if( temp_distance < min_distance )
	{
		min_distance = temp_distance; 
		nearest_boundary = 0; 
	}
	temp_distance = get_microenvironment_i()->get_mesh().bounding_box[3] - pCell->get_position()[0]; 
	if( temp_distance < min_distance )
	{
		min_distance = temp_distance; 
		nearest_boundary = 1; 
	}
	
	// check against yL and yU
	temp_distance = pCell->get_position()[1] - get_microenvironment_i()->get_mesh().bounding_box[1]; 
	if( temp_distance < min_distance )
	{
		min_distance = temp_distance; 
		nearest_boundary = 2; 
	}
	temp_distance = get_microenvironment_i()->get_mesh().bounding_box[4] - pCell->get_position()[1]; 
	if( temp_distance < min_distance )
	{
		min_distance = temp_distance; 
		nearest_boundary = 3; 
	}	
	
	if( get_microenvironment_i()->simulate_2D() == false )
	{
		// if in 3D, check against zL and zU
		temp_distance = pCell->get_position()[2] - get_microenvironment_i()->get_mesh().bounding_box[2]; 
		if( temp_distance < min_distance )
		{
			min_distance = temp_distance; 
			nearest_boundary = 4; 
		}
		temp_distance = get_microenvironment_i()->get_mesh().bounding_box[5] - pCell->get_position()[2]; 
		if( temp_distance < min_distance )
		{
			min_distance = temp_distance; 
			nearest_boundary = 5; 
		}			
		
		// check for 3D exceptions 
		
		// lines 
		if( fabs( (pCell->get_position()[0]) - (pCell->get_position()[1]) ) < tolerance && 
			fabs( (pCell->get_position()[1]) - (pCell->get_position()[2]) ) < tolerance && 
			fabs( (pCell->get_position()[0]) - (pCell->get_position()[2]) ) < tolerance )
		{
			if( pCell->get_position()[0] > 0 )
			{
				if( pCell->get_position()[0] > 0 && pCell->get_position()[1] > 0 )
				{ pCell->displacement = { -one_over_sqrt_3 , -one_over_sqrt_3 , -one_over_sqrt_3 }; }
				if( pCell->get_position()[0] < 0 && pCell->get_position()[1] > 0 )
				{ pCell->displacement = { one_over_sqrt_3 , -one_over_sqrt_3 , -one_over_sqrt_3 }; }
				
				if( pCell->get_position()[0] > 0 && pCell->get_position()[1] < 0 )
				{ pCell->displacement = { -one_over_sqrt_3 , one_over_sqrt_3 , -one_over_sqrt_3 }; }
				if( pCell->get_position()[0] < 0 && pCell->get_position()[1] < 0 )
				{ pCell->displacement = { one_over_sqrt_3 , one_over_sqrt_3 , -one_over_sqrt_3 }; }
			} 
			else
			{
				if( pCell->get_position()[0] > 0 && pCell->get_position()[1] > 0 )
				{ pCell->displacement = { -one_over_sqrt_3 , -one_over_sqrt_3 , one_over_sqrt_3 }; }
				if( pCell->get_position()[0] < 0 && pCell->get_position()[1] > 0 )
				{ pCell->displacement = { one_over_sqrt_3 , -one_over_sqrt_3 , one_over_sqrt_3 }; }
				
				if( pCell->get_position()[0] > 0 && pCell->get_position()[1] < 0 )
				{ pCell->displacement = { -one_over_sqrt_3 , one_over_sqrt_3 , one_over_sqrt_3 }; }
				if( pCell->get_position()[0] < 0 && pCell->get_position()[1] < 0 )
				{ pCell->displacement = { one_over_sqrt_3 , one_over_sqrt_3 , one_over_sqrt_3 }; }				
			}
			return min_distance; 
		}
		
		// planes - let's not worry for today 
		
	}
	else
	{
		// check for 2D  exceptions 
		
		if( fabs( (pCell->get_position()[0]) - (pCell->get_position()[1]) ) < tolerance )
		{
			if( pCell->get_position()[0] > 0 && pCell->get_position()[1] > 0 )
			{ pCell->displacement = { -one_over_sqrt_2 , -one_over_sqrt_2 , 0 }; }
			if( pCell->get_position()[0] < 0 && pCell->get_position()[1] > 0 )
			{ pCell->displacement = { one_over_sqrt_2 , -one_over_sqrt_2 , 0 }; }
			
			if( pCell->get_position()[0] > 0 && pCell->get_position()[1] < 0 )
			{ pCell->displacement = { -one_over_sqrt_2 , one_over_sqrt_2 , 0 }; }
			if( pCell->get_position()[0] < 0 && pCell->get_position()[1] < 0 )
			{ pCell->displacement = { one_over_sqrt_2 , one_over_sqrt_2 , 0 }; }
			return min_distance; 
		}
	}
	
	// no exceptions 
	switch(nearest_boundary)
	{
		case 0:
			pCell->displacement = {1,0,0}; 
			return min_distance; 
		case 1:
			pCell->displacement = {-1,0,0}; 
			return min_distance;
		case 2:
			pCell->displacement = {0,1,0}; 
			return min_distance; 
		case 3: 
			pCell->displacement = {0,-1,0}; 
			return min_distance; 
		case 4: 
			pCell->displacement = {0,0,1}; 
			return min_distance; 
		case 5: 
			pCell->displacement = {0,0,-1}; 
			return min_distance; 
		default:
			pCell->displacement = {0,0,0};
			return 9e99; 
	}
	
	pCell->displacement = {0,0,0};
	return 9e99; 
}	


void Mechanics_Standard_Models::chemotaxis_function( Mechanics_Agent_PIMPL* pCellI, double dt )
{
	Mechanics_Agent* pCell = dynamic_cast<Mechanics_Agent*>(pCellI->get_mechanics_implementation());

	// bias direction is gradient for the indicated substrate 
	pCell->motility_data.migration_bias_direction = get_microenvironment_i()->nearest_gradient_vector(pCell->pos_entity->position)[pCell->motility_data.chemotaxis_index];
	// move up or down gradient based on this direction 
	pCell->motility_data.migration_bias_direction *= pCell->motility_data.chemotaxis_direction; 

	// normalize 
	normalize( &( pCell->motility_data.migration_bias_direction ) );
	
	return;
}

void Mechanics_Standard_Models::advanced_chemotaxis_function_normalized( Mechanics_Agent_PIMPL* pCellI, double dt )
{
	Mechanics_Agent* pCell = dynamic_cast<Mechanics_Agent*>(pCellI->get_mechanics_implementation());

	// We'll work directly on the migration bias direction 
	std::vector<double>* pVec = &(pCell->motility_data.migration_bias_direction);  
	// reset to zero. use memset to be faster??
	pVec->assign( 3, 0.0 ); 
	
	// a place to put each gradient prior to normalizing it 
	std::vector<double> temp(3,0.0); 

	// weighted combination of the gradients 
	for( int i=0; i < pCell->motility_data.chemotactic_sensitivities.size(); i++ )
	{
		// get and normalize ith gradient 
		temp = get_microenvironment_i()->nearest_gradient_vector(pCell->pos_entity->position)[i]; 
		normalize( &temp ); 
		axpy( pVec , pCell->motility_data.chemotactic_sensitivities[i] , temp ); 
	}
	// normalize that 
	normalize( pVec ); 
	
	return;
}

void Mechanics_Standard_Models::advanced_chemotaxis_function( Mechanics_Agent_PIMPL* pCellI, double dt )
{
	Mechanics_Agent* pCell = dynamic_cast<Mechanics_Agent*>(pCellI->get_mechanics_implementation());

	// We'll work directly on the migration bias direction 
	std::vector<double>* pVec = &(pCell->motility_data.migration_bias_direction);  
	// reset to zero. use memset to be faster??
	pVec->assign( 3, 0.0 ); 

	// weighted combination of the gradients 
	for( int i=0; i < pCell->motility_data.chemotactic_sensitivities.size(); i++ )
	{
		// get and normalize ith gradient 
		axpy( pVec , pCell->motility_data.chemotactic_sensitivities[i] , get_microenvironment_i()->nearest_gradient_vector(pCell->pos_entity->position)[i] ); 
	}
	// normalize that 
	normalize( pVec );

/*
 #pragma omp critical
 {
	std::cout << "\t\ttype: " << pCell->type_name 
	<< " bias: " << pCell->motility_data.migration_bias 
	<< " speed: " << pCell->motility_data.migration_speed 
	<< " direction: " << pCell->motility_data.migration_bias_direction << std::endl; 
 */

	return;
}

void Mechanics_Standard_Models::standard_elastic_contact_function( Mechanics_Agent_PIMPL* pC1I, Mechanics_Agent_PIMPL* pC2I, double dt )
{
	standard_elastic_contact_function(
		dynamic_cast<Mechanics_Agent*>(pC1I->get_mechanics_implementation()),
		dynamic_cast<Mechanics_Agent*>(pC2I->get_mechanics_implementation()),
		dt
	);
}

void Mechanics_Standard_Models::standard_elastic_contact_function( Mechanics_Agent* pC1, Mechanics_Agent* pC2, double dt )
{
	if( pC1->get_position().size() != 3 || pC2->get_position().size() != 3 )
	{ return; }
	
	std::vector<double> displacement = pC2->get_position();
	displacement -= pC1->get_position(); 

	// update May 2022 - effective adhesion 
	int ii = find_cell_definition_index( pC1->get_type() ); 
	int jj = find_cell_definition_index( pC2->get_type() ); 

	double adhesion_ii = pC1->mechanics_data.attachment_elastic_constant * pC1->mechanics_data.cell_adhesion_affinities[jj]; 
	double adhesion_jj = pC2->mechanics_data.attachment_elastic_constant * pC2->mechanics_data.cell_adhesion_affinities[ii]; 

	double effective_attachment_elastic_constant = sqrt( adhesion_ii*adhesion_jj ); 

	// axpy( &(pC1->velocity) , p1.mechanics.attachment_elastic_constant , displacement ); 
	axpy( &(pC1->velocity) , effective_attachment_elastic_constant , displacement ); 
	return; 
}

void Mechanics_Standard_Models::standard_elastic_contact_function_confluent_rest_length( Mechanics_Agent_PIMPL* pC1I, Mechanics_Agent_PIMPL* pC2I, double dt )
{
	Mechanics_Agent* pC1 = dynamic_cast<Mechanics_Agent*>(pC1I->get_mechanics_implementation());
	Mechanics_Agent* pC2 = dynamic_cast<Mechanics_Agent*>(pC2I->get_mechanics_implementation());

	if( pC1->get_position().size() != 3 || pC2->get_position().size() != 3 )
	{ return; }
	
	std::vector<double> displacement = pC2->get_position();
	displacement -= pC1->get_position(); 

	// update May 2022 - effective adhesion 
	int ii = find_cell_definition_index( pC1->get_type() ); 
	int jj = find_cell_definition_index( pC2->get_type() ); 

	double adhesion_ii = pC1->mechanics_data.attachment_elastic_constant * pC1->mechanics_data.cell_adhesion_affinities[jj]; 
	double adhesion_jj = pC2->mechanics_data.attachment_elastic_constant * pC2->mechanics_data.cell_adhesion_affinities[ii]; 

	double effective_attachment_elastic_constant = sqrt( adhesion_ii*adhesion_jj ); 
	// axpy( &(pC1->velocity) , effective_attachment_elastic_constant , displacement ); 

	// have the adhesion strength taper away at this rest lenght
	// set the rest length = confluent cell-cell spacing 
	// 
	double rest_length = ( pC1->radius_data.radius + pC2->radius_data.radius ) * 0.9523809523809523;  

	double strength = ( norm(displacement) - rest_length )*effective_attachment_elastic_constant;
	normalize( &displacement );
	axpy( &(pC1->velocity) , strength , displacement ); 

	return; 
}

void Mechanics_Standard_Models::dynamic_spring_attachments( Mechanics_Agent_PIMPL* pCellI , double dt )
{
	dynamic_spring_attachments(dynamic_cast<Mechanics_Agent*>(pCellI->get_mechanics_implementation()), dt);
}

void Mechanics_Standard_Models::dynamic_spring_attachments( Mechanics_Agent* pCell , double dt )
{
    if( !pCell ) { return; }

    // check for detachments 
    double detachment_probability = pCell->mechanics_data.detachment_rate * dt; 
	// detach_cells_as_spring swaps the detached cell with the last cell in the vector, so we need to iterate backwards
    for( int j=pCell->spring_attachments.size()-1; j >= 0; j-- )
    {
        std::pair<Mechanics_Agent_PIMPL*, bool> spring_pair = pCell->spring_attachments[j];
		Mechanics_Agent* pTest = dynamic_cast<Mechanics_Agent*>(spring_pair.first->get_mechanics_implementation());
		bool atacking_cell = spring_pair.second;
		if (atacking_cell == true) // do not let attackers detach randomly
		{ continue; }
        if( UniformRandom() <= detachment_probability )
        { get_mechanics_environment().detach_cells_as_spring( pCell , pTest ); }
    }

    // check if I have max number of attachments 
				if( pCell->spring_attachments.size() >= pCell->mechanics_data.maximum_number_of_attachments )
    { return; }

    // check for new attachments; 
    double attachment_probability = pCell->mechanics_data.attachment_rate * dt; 
    bool done = false; 
    int j = 0; 
    while( done == false && j < pCell->neighbors.size() )
    {
        Mechanics_Agent* pTest = dynamic_cast<Mechanics_Agent*>(pCell->neighbors[j]->get_mechanics_implementation());
		if( pTest->spring_attachments.size() < pTest->mechanics_data.maximum_number_of_attachments )
		{
			// std::string search_string = "adhesive affinity to " + pTest->type_name; 
			// double affinity = get_single_behavior( pCell , search_string );
			double affinity = pCell->mechanics_data.cell_adhesion_affinities[cell_definition_indices_by_type[pTest->get_type()]];

            double prob = attachment_probability * affinity; 
            if( UniformRandom() <= prob )
            {
                // attempt the attachment. testing for prior connection is already automated 
                get_mechanics_environment().attach_cells_as_spring( pCell, pTest, false ); 
                if( pCell->spring_attachments.size() >= pCell->mechanics_data.maximum_number_of_attachments )
                { done = true; }
            }
        }
        j++; 
    }
    return; 
}

}
