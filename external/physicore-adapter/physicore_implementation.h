#ifndef __physicore_implementation_h__
#define __physicore_implementation_h__

#include "../../BioFVM/BioFVM_implementation.h"

namespace BioFVM {

class physicore_implementation : public BioFVM::BioFVM_implementation 
{
public:
    Microenvironment_Interface* get_microenvironment() override;
    Basic_Agent_Interface* create_basic_agent() override;
    std::vector<Basic_Agent_Interface*>* get_all_basic_agents() override;
};

} // namespace BioFVM

#endif
