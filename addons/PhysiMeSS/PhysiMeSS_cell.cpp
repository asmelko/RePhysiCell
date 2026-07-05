#include "PhysiMeSS_cell.h"
#include "PhysiMeSS_fibre.h"
#include "PhysiMeSS_environment.h"

#include "../../BioFVM/BioFVM_vector.h"
#include "../../core/PhysiCell_cell.h"

using namespace PhysiCell;
using namespace BioFVM;

PhysiMeSS_CellAgent::PhysiMeSS_CellAgent(Cell* pCell)
    : PhysiMeSS_Agent(pCell)
{
    stuck_counter = 0;
    unstuck_counter = 0;
}


void PhysiMeSS_CellAgent::register_fibre_voxels()
{
    //a cell will be in one voxel
    int voxel = physimess_environment.mechanics_mesh.nearest_voxel_index(get_position());
    physimess_voxels.push_back(voxel);

}

void PhysiMeSS_CellAgent::deregister_fibre_voxels() {

    //only do this for fibres
    return;
}


void PhysiMeSS_CellAgent::add_potentials_from_fibre(PhysiMeSS_FibreAgent* pFibre)
{
    Cell* pCell = get_cell();
    double distance = 0.0;
    pFibre->nearest_point_on_fibre(get_position(), displacement);
    for (int index = 0; index < 3; index++) {
        distance += displacement[index] * displacement[index];
    }
    distance = std::max(sqrt(distance), 0.00001);
    /*if( this->phenotype.motility.is_motile) {
        std::cout << " determining distance from " << this->type_name << " " << this->ID << " to "
                    << (*other_agent).type_name << " " << (*other_agent).ID
                    << "   the distance is " << distance << std::endl;
    }*/

    // as per PhysiCell
    static double simple_pressure_scale = 0.027288820670331;

    // check distance relative repulsion and adhesion distances
    // cell should repel from a fibre if it comes within cell radius plus fibre radius (note fibre radius ~2 micron)
    double R = radius_data.radius + pFibre->mRadius;
    // cell should feel adhesion over
    double max_interactive_distance =
            mechanics_data.relative_maximum_adhesion_distance * radius_data.radius+
            pFibre->mechanics_data.relative_maximum_adhesion_distance *
            pFibre->mRadius;

    // First Repulsion as per PhysiCell
    double temp_r = 0;
    if (distance > R) {
        temp_r = 0;
    } else {
        // temp_r = 1 - distance/R;
        temp_r = -distance; // -d
        temp_r /= R; // -d/R
        temp_r += 1.0; // 1-d/R
        temp_r *= temp_r; // (1-d/R)^2

        // add the relative pressure contribution NOT SURE IF NEEDED
        simple_pressure += (temp_r / simple_pressure_scale);

        double effective_repulsion = sqrt(mechanics_data.cell_cell_repulsion_strength *
                                         pFibre->mechanics_data.cell_cell_repulsion_strength);
        temp_r *= effective_repulsion;
    }

    if (fabs(temp_r) < 1e-16) { return; }
    temp_r /= distance;

    axpy(&velocity, temp_r, displacement);

    //Then additional repulsion/adhesion as per Cicely's code
    double fibre_adhesion = 0;
    double fibre_repulsion = 0;
    if (distance < max_interactive_distance) {
        const std::vector<double> previous_velocity = this->previous_velocity;
        double cell_velocity_dot_fibre_direction = 0.;
        for (unsigned int j = 0; j < 3; j++) {
            cell_velocity_dot_fibre_direction += pFibre->get_cell()->state.orientation[j] * previous_velocity[j];
        }
        double cell_velocity = 0;
        for (unsigned int j = 0; j < previous_velocity.size(); j++) {
            cell_velocity += previous_velocity[j] * previous_velocity[j];
        }
        cell_velocity = std::max(sqrt(cell_velocity), 1e-8);

        double p_exponent = 1.0;
        double q_exponent = 1.0;
        double xi = fabs(cell_velocity_dot_fibre_direction) / (cell_velocity);
        double xip = pow(xi, p_exponent);
        double xiq = pow((1 - xi * xi), q_exponent);

        fibre_adhesion = pCell->custom_data["vel_adhesion"] * xip *
                         (1 - cell_velocity / pCell->custom_data["cell_velocity_max"]);

        fibre_repulsion = pCell->custom_data["vel_contact"] * xiq;

        axpy(&velocity, fibre_adhesion, pFibre->get_cell()->state.orientation);
        naxpy(&velocity, fibre_repulsion, previous_velocity);

        degrade_fibre(pFibre);
    }
}


void PhysiMeSS_CellAgent::degrade_fibre(PhysiMeSS_FibreAgent* pFibre)
{
    Cell* pCell = get_cell();
    double distance = 0.0;
    pFibre->nearest_point_on_fibre(get_position(), displacement);
    for (int index = 0; index < 3; index++) {
        distance += displacement[index] * displacement[index];
    }
    distance = std::max(sqrt(distance), 0.00001);

    // Fibre degradation by cell - switched on by flag fibre_degradation
    double stuck_threshold = pCell->custom_data["fibre_stuck_time"];
    if (pCell->custom_data["fibre_degradation"] > 0.5 && stuck_counter >= stuck_threshold) {
        // if (stuck_counter >= stuck_threshold){
        //     std::cout << "Cell " << ID << " is stuck at time " << PhysiCell::PhysiCell_globals.current_time
        //                 << " near fibre " << pFibre->ID  << std::endl;;
        // }
        displacement *= -1.0/distance;
        double dotproduct = dot_product(displacement, motility_data.motility_vector);
        if (dotproduct >= 0) {
            double rand_degradation = PhysiCell::UniformRandom();
            double prob_degradation = pCell->custom_data["fibre_degradation_rate"];
            if (rand_degradation <= prob_degradation) {
                //std::cout << " --------> fibre " << (*other_agent).ID << " is flagged for degradation " << std::endl;
                // (*other_agent).parameters.degradation_flag = true;
                pFibre->get_cell()->flag_for_removal();
                // std::cout << "Degrading fibre agent " << pFibre->ID << " using flag for removal !!" << std::endl;
                stuck_counter = 0;
            }
        }
    }
}

void PhysiMeSS_CellAgent::force_update_motility_vector(double dt_) {
    Cell* pCell = get_cell();
    if (!motility_data.is_motile) {
        motility_data.motility_vector.assign(3, 0.0);
        return;
    }

    // force cell to update its motility because it is stuck
    // choose a uniformly random unit vector
    double temp_angle = 6.28318530717959 * PhysiCell::UniformRandom();
    double temp_phi = 3.1415926535897932384626433832795 * PhysiCell::UniformRandom();

    double sin_phi = sin(temp_phi);
    double cos_phi = cos(temp_phi);

    if (motility_data.restrict_to_2D)
    {
        sin_phi = 1.0;
        cos_phi = 0.0;
    }

    std::vector<double> randvec;
    randvec.resize(3, sin_phi);

    randvec[0] *= cos(temp_angle); // cos(theta)*sin(phi)
    randvec[1] *= sin(temp_angle); // sin(theta)*sin(phi)
    randvec[2] = cos_phi; //  cos(phi)

    // if the update_bias_vector function is set, use it
    /*if (functions.update_migration_bias) {
    functions.update_migration_bias(this, phenotype, dt_);
    }*/

    //phenotype.motility.motility_vector *= -1.0;//phenotype.motility.migration_bias_direction; // motiltiy = bias_vector
    //phenotype.motility.motility_vector *= phenotype.motility.migration_bias; // motility = bias*bias_vector

    double one_minus_bias = 1.0;// - phenotype.motility.migration_bias;

    axpy(&motility_data.motility_vector, one_minus_bias, randvec); // motility = (1-bias)*randvec + bias*bias_vector

    normalize(&motility_data.motility_vector);

    motility_data.motility_vector *= motility_data.migration_speed;

    return;
}


PhysiMeSS_Cell::PhysiMeSS_Cell()
    : Cell()
{
    // Swap pImpl from default Mechanics_Agent to PhysiMeSS_CellAgent
    delete Mechanics_Agent_PIMPL::pImpl;
    Mechanics_Agent_PIMPL::pImpl = new PhysiMeSS_CellAgent(this);
    Mechanics_Agent_PIMPL::pImpl->bind_position_entity(this);
}

PhysiMeSS_CellAgent* PhysiMeSS_Cell::get_physimess_agent() const
{
    return const_cast<PhysiMeSS_CellAgent*>(
        static_cast<const PhysiMeSS_CellAgent*>(get_mechanics_implementation()));
}
