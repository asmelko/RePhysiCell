#ifndef __PhysiCell_mechanics_standard_models_interface_h__
#define __PhysiCell_mechanics_standard_models_interface_h__

#include "../core_interface/PhysiCell_core_entity_PIMPL.h"

namespace PhysiCell {

class Mechanics_Standard_Models_Interface
{
public:
    virtual void standard_add_basement_membrane_interactions( CoreEntityPIMPL* pCell, double dt ) = 0;
    virtual void standard_domain_edge_avoidance_interactions( CoreEntityPIMPL* pCell, double dt ) = 0;
    virtual double distance_to_domain_edge(CoreEntityPIMPL* pCell) = 0;

    virtual void dynamic_spring_attachments( CoreEntityPIMPL* pCell , double dt ) = 0;
    virtual void standard_elastic_contact_function( CoreEntityPIMPL* pC1, CoreEntityPIMPL* pC2, double dt ) = 0;
    virtual void standard_elastic_contact_function_confluent_rest_length( CoreEntityPIMPL* pC1I, CoreEntityPIMPL* pC2I, double dt ) = 0;

    virtual void chemotaxis_function( CoreEntityPIMPL* pCellI, double dt ) = 0;
    virtual void advanced_chemotaxis_function_normalized( CoreEntityPIMPL* pCellI, double dt ) = 0;
    virtual void advanced_chemotaxis_function( CoreEntityPIMPL* pCellI, double dt ) = 0;
};

}

#endif // __PhysiCell_mechanics_standard_models_interface_h__