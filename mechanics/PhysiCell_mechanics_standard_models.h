#ifndef __PhysiCell_mechanics_standard_models_h__
#define __PhysiCell_mechanics_standard_models_h__

#include "PhysiCell_mechanics_standard_models_interface.h"
#include "PhysiCell_mechanics_agent.h"

namespace PhysiCell {

class Mechanics_Standard_Models : public Mechanics_Standard_Models_Interface
{
public:
    void standard_add_basement_membrane_interactions( CoreEntityPIMPL* pCell, double dt ) override;
    void standard_domain_edge_avoidance_interactions( CoreEntityPIMPL* pCell, double dt ) override;
    double distance_to_domain_edge(CoreEntityPIMPL* pCell) override;

    void dynamic_spring_attachments( CoreEntityPIMPL* pCell , double dt ) override;
    void dynamic_spring_attachments( Mechanics_Agent* pCell , double dt );
    void standard_elastic_contact_function( CoreEntityPIMPL* pC1, CoreEntityPIMPL* pC2, double dt ) override;
    void standard_elastic_contact_function( Mechanics_Agent* pC1, Mechanics_Agent* pC2, double dt );
    void standard_elastic_contact_function_confluent_rest_length( CoreEntityPIMPL* pC1I, CoreEntityPIMPL* pC2I, double dt ) override;

    void chemotaxis_function( CoreEntityPIMPL* pCellI, double dt ) override;
    void advanced_chemotaxis_function_normalized( CoreEntityPIMPL* pCellI, double dt ) override;
    void advanced_chemotaxis_function( CoreEntityPIMPL* pCellI, double dt ) override;
};

}

#endif // __PhysiCell_mechanics_standard_models_h__