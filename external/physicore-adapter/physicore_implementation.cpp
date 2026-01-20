#include "physicore_implementation.h"
#include "microenvironment_adapter.h"

#include <memory>

namespace BioFVM {

Microenvironment_Interface* physicore_implementation::get_microenvironment() {
    static std::unique_ptr<physicore_microenvironment_adapter> global_adapter = std::make_unique<physicore_microenvironment_adapter>();

	return global_adapter.get();
}

Basic_Agent_Interface* physicore_implementation::create_basic_agent() {
    return new physicore_basic_agent_adapter();
}


std::vector<Basic_Agent_Interface*>* physicore_implementation::get_all_basic_agents(){
    static std::vector<Basic_Agent_Interface*> all_basic_agents;
    return &all_basic_agents;
}

} // namespace BioFVM
