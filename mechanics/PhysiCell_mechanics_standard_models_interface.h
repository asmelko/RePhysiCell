#ifndef __PhysiCell_mechanics_standard_models_interface_h__
#define __PhysiCell_mechanics_standard_models_interface_h__

#include "PhysiCell_mechanics_agent_PIMPL.h"

namespace PhysiCell {

class Mechanics_Standard_Models_Interface
{
public:
    virtual void standard_add_basement_membrane_interactions( Mechanics_Agent_PIMPL* pCell, double dt ) = 0;
    virtual void standard_domain_edge_avoidance_interactions( Mechanics_Agent_PIMPL* pCell, double dt ) = 0;
    virtual double distance_to_domain_edge(Mechanics_Agent_PIMPL* pCell) = 0;

    virtual void dynamic_spring_attachments( Mechanics_Agent_PIMPL* pCell , double dt ) = 0;
    virtual void standard_elastic_contact_function( Mechanics_Agent_PIMPL* pC1, Mechanics_Agent_PIMPL* pC2, double dt ) = 0;
    virtual void standard_elastic_contact_function_confluent_rest_length( Mechanics_Agent_PIMPL* pC1I, Mechanics_Agent_PIMPL* pC2I, double dt ) = 0;

    virtual void chemotaxis_function( Mechanics_Agent_PIMPL* pCellI, double dt ) = 0;
    virtual void advanced_chemotaxis_function_normalized( Mechanics_Agent_PIMPL* pCellI, double dt ) = 0;
    virtual void advanced_chemotaxis_function( Mechanics_Agent_PIMPL* pCellI, double dt ) = 0;
};

}

#endif // __PhysiCell_mechanics_standard_models_interface_h__