#include "PhysiCell_mechanics_functions.h"

#include "../core/PhysiCell_cell.h"
#include "PhysiCell_mechanics_standard_models.h"

using namespace PhysiCell;

void Mechanics_Functions::add_cell_basement_membrane_interactions( double dt )
{
    if ( pCell->functions.add_cell_basement_membrane_interactions != nullptr )
    { pCell->functions.add_cell_basement_membrane_interactions( pCell, pCell->phenotype, dt ); }
    else
    { return; }
}

double Mechanics_Functions::calculate_distance_to_membrane( double dt )
{
    if ( pCell->functions.calculate_distance_to_membrane != nullptr )
    { return pCell->functions.calculate_distance_to_membrane( pCell, pCell->phenotype, dt ); }
    else
    { return Mechanics_Standard_Models().distance_to_domain_edge(pCell); }
}

void Mechanics_Functions::update_migration_bias( double dt )
{
    if ( pCell->functions.update_migration_bias != nullptr )
    { pCell->functions.update_migration_bias( pCell, pCell->phenotype, dt ); }
    else
    { return; }
}
