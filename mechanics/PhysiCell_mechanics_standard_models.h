#ifndef __PhysiCell_mechanics_standard_models_h__
#define __PhysiCell_mechanics_standard_models_h__

#include "PhysiCell_mechanics_standard_models_interface.h"
#include "PhysiCell_mechanics_agent.h"

namespace PhysiCell {

class Mechanics_Standard_Models : public Mechanics_Standard_Models_Interface
{
public:
    void standard_add_basement_membrane_interactions( Mechanics_Agent_PIMPL* pCell, double dt ) override;
    void standard_domain_edge_avoidance_interactions( Mechanics_Agent_PIMPL* pCell, double dt ) override;
    double distance_to_domain_edge(Mechanics_Agent_PIMPL* pCell) override;

    void dynamic_spring_attachments( Mechanics_Agent_PIMPL* pCell , double dt ) override;
    void dynamic_spring_attachments( Mechanics_Agent* pCell , double dt );
    void standard_elastic_contact_function( Mechanics_Agent_PIMPL* pC1, Mechanics_Agent_PIMPL* pC2, double dt ) override;
    void standard_elastic_contact_function( Mechanics_Agent* pC1, Mechanics_Agent* pC2, double dt );
    void standard_elastic_contact_function_confluent_rest_length( Mechanics_Agent_PIMPL* pC1I, Mechanics_Agent_PIMPL* pC2I, double dt ) override;

    void chemotaxis_function( Mechanics_Agent_PIMPL* pCellI, double dt ) override;
    void advanced_chemotaxis_function_normalized( Mechanics_Agent_PIMPL* pCellI, double dt ) override;
    void advanced_chemotaxis_function( Mechanics_Agent_PIMPL* pCellI, double dt ) override;
};

}

#endif // __PhysiCell_mechanics_standard_models_h__