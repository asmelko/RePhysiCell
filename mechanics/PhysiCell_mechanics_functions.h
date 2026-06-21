#ifndef __PhysiCell_mechanics_functions_h__
#define __PhysiCell_mechanics_functions_h__

namespace PhysiCell {

class Cell;

class Mechanics_Functions
{
    Cell* pCell;
public:
    Mechanics_Functions( Cell* _pCell ) : pCell(_pCell) {}
    
	void add_cell_basement_membrane_interactions( double dt );
	double calculate_distance_to_membrane( double dt );
	void update_migration_bias( double dt );
};

}

#endif // __PhysiCell_mechanics_functions_h__
