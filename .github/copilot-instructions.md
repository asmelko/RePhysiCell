# PhysiCell AI Agent Instructions (RePhysiCell Fork)

## Project Overview
**RePhysiCell** is a fork of **PhysiCell** (v1.14.2), an open-source physics-based cell simulator for 3-D multicellular systems. This is a C++ framework for agent-based modeling of cells with rich phenotypic behaviors, diffusible substrates (via BioFVM), and optional intracellular models.

**Key Citation**: Ghaffarizadeh et al., PLoS Comput. Biol. 14(2): e1005991, 2018

**Branch**: `biofvm-interface-x` introduces a critical architectural change: **abstraction layer between PhysiCell and BioFVM**.

## Architecture Overview

### Key Architectural Innovation: Interface/Adapter Pattern
**CRITICAL**: This fork introduces an abstraction layer separating PhysiCell from BioFVM:

- **`BioFVM/BioFVM_microenvironment_interface.h`**: Abstract interface defining microenvironment operations
  - Pure virtual interface with ~50 methods for microenvironment access
  - Allows swapping BioFVM for alternative solvers without modifying PhysiCell core
  - All PhysiCell code accesses microenvironment via `get_microenvironment_i()`

- **`BioFVM/microenvironment_adapter.{h,cpp}`**: Concrete adapter wrapping BioFVM
  - Implements interface by delegating to `BioFVM::Microenvironment`
  - Global singleton `global_adapter` initialized via `BioFVM::initialize_microenvironment_interface()`
  - Manages ownership of wrapped BioFVM instance

**Usage Pattern**:
```cpp
// In main.cpp - MUST call first
BioFVM::initialize_microenvironment_interface();

// Access microenvironment throughout codebase
get_microenvironment_i()->simulate_diffusion_decay(dt);
get_microenvironment_i()->find_density_index("oxygen");
((Cell_Container*)get_microenvironment_i()->get_agent_container())->update_all_cells(t);
```

### Core Components
- **BioFVM/**: Diffusion-reaction PDE solver for microenvironment (substrates like oxygen, nutrients)
  - `BioFVM_microenvironment.*` - spatial discretization, diffusion solver (wrapped by adapter)
  - `BioFVM_basic_agent.*` - base class for all agents
  - `microenvironment_adapter.*` - **NEW**: adapter implementation
- **bio_interface/**: **NEW**: Abstract interface layer
  - `Bio_microenvironment_interface.h` - interface definition
  - Header-only library defined in CMake
- **core/**: PhysiCell cell simulation engine
  - `PhysiCell_cell.*` - `Cell` class (inherits from `BioFVM::Basic_Agent`)
  - `PhysiCell_phenotype.*` - `Phenotype` class with cycle, death, volume, mechanics, motility, secretion, interactions
  - `PhysiCell_cell_container.*` - spatial binning for efficient neighbor searches
  - `PhysiCell_standard_models.*` - built-in cycle/death models (Ki67, live, apoptosis, necrosis)
  - `PhysiCell_rules.*` - rules-based modeling (Cell Behavior Hypothesis Grammar v3)
  - `PhysiCell_signal_behavior.*` - signal/behavior dictionaries for intracellular coupling
- **modules/**: I/O, visualization, settings
  - `PhysiCell_settings.*` - XML config file parsing
  - `PhysiCell_SVG.*` - SVG snapshot generation
  - `PhysiCell_MultiCellDS.*` - data output (MultiCellDS format)
- **addons/**: Optional modules (PhysiBoSS, libRoadrunner, dFBA, PhysiMeSS)
  - **Note**: `libRoadrunner` addon updated to use interface (`get_microenvironment_i()`)
- **sample_projects/**: Example models demonstrating features (all updated with interface call)
- **custom_modules/**: User-defined C++ code (custom phenotype functions, setup)

### Data Flow
1. **Initialization**: 
   - `main.cpp` → **`BioFVM::initialize_microenvironment_interface()`** (wraps global BioFVM microenvironment)
   - → `load_PhysiCell_config_file()` → `setup_microenvironment()` → `create_cell_types()` → `setup_tissue()`
2. **Time loop**: 
   - `get_microenvironment_i()->simulate_diffusion_decay()` (BioFVM diffusion via interface)
   - → cell phenotype updates → cell mechanics → secretion/uptake → SVG/data output
3. **Configuration**: XML files in `./config/` define microenvironment, cell definitions, rules, user parameters

## Critical Build & Test Workflows

### Build Systems
**Two parallel build systems** (currently transitioning to CMake):
1. **Makefile** (legacy, still primary):
   ```bash
   make                    # Compile current project
   make template           # Load template project
   make clean              # Remove .o files
   make data-cleanup       # Clear output/
   make reset              # Reset to clean state
   ```
2. **CMake** (modern, on `cmake` branch):
   ```bash
   mkdir build && cd build
   cmake .. -DPHYSICELL_BUILD_PHYSIBOSS=ON  # Enable addons
   cmake --build . --target template
   ```
   CMake options: `PHYSICELL_BUILD_PHYSIBOSS`, `PHYSICELL_BUILD_LIBROADDRUNNER`, `PHYSICELL_BUILD_DFBA`
   - **Note**: `bio_interface` compiled as header-only INTERFACE library
   - BioFVM includes `microenvironment_adapter.cpp` in its OBJECT library

### Testing
- **Unit tests**: `./unit_tests/` - basic functionality tests (compile with `make` in that dir)
- **System tests**: `./tests/system/` - integration tests via Python scripts
- **CI/CD**: `.github/workflows/tests_cmake.yml` - automated cross-platform testing (Windows/Ubuntu/macOS)
- **Sample validation**: `beta/test_run_sample.py` runs projects, `beta/test_diff_svg.py` compares outputs

### Platform-Specific Notes
- **macOS**: Must define `PHYSICELL_CPP` env variable pointing to OpenMP-enabled g++ (e.g., `g++-15`)
- **Windows**: Requires MinGW-w64, setup via `beta/setup_windows_dep.py`
- **OpenMP**: Required for parallelization (`-fopenmp` flag)

## Project-Specific Conventions

### Cell Type Definition Pattern
**Modern approach** (v1.7+): Define cell types in XML (`./config/PhysiCell_settings.xml`):
```xml
<cell_definitions>
  <cell_definition name="cancer cell" ID="0">
    <phenotype>
      <cycle code="5" name="live">...</cycle>
      <death>...</death>
      <volume>...</volume>
      <mechanics>...</mechanics>
      <motility>...</motility>
      <secretion>...</secretion>
      <cell_interactions>...</cell_interactions>
      <cell_transformations>...</cell_transformations>
    </phenotype>
  </cell_definition>
</cell_definitions>
```
**Legacy approach**: C++ in `custom_modules/custom.cpp` via `create_cell_types()` function.

### Rules-Based Modeling (CBHG v3)
PhysiCell supports CSV-based hypothesis specification (`./config/cell_rules.csv`):
```
cell_type,signal,direction,behavior,max_response_value,half_max,Hill_power,applies_to_dead
cancer cell,oxygen,increases,cycle entry,0.01,21.5,4,0
```
This auto-generates Hill response functions at runtime—**no C++ recompilation needed**.

### Contact Functions & Cell Interactions
- Attach cells: `attach_cells(pCell1, pCell2)` or `pCell->attach_cell(pOther)`
- Contact function signature: `void my_contact(Cell* pMe, Phenotype&, Cell* pOther, Phenotype&, double dt)`
- Built-in interactions: `standard_elastic_contact_function`, phagocytosis (`cell_interactions.live_phagocytosis_rates`), attack (`cell_interactions.attack_rates`), fusion
- Neighbor search: `pCell->nearby_interacting_cells()` returns vector of cells within interaction distance

### Custom Modules Structure
Sample projects organize code as:
- `main.cpp`: Entry point, loads config, runs time loop
- `custom_modules/custom.cpp`: User phenotype functions
  - `create_cell_types()` - define cell definitions in C++
  - `setup_tissue()` - place initial cells
  - Custom phenotype functions (e.g., `tumor_cell_phenotype_with_signaling()`)
- `config/PhysiCell_settings.xml`: Simulation parameters
- `config/cell_rules.csv`: Optional rules-based behaviors

### Signal & Behavior Dictionaries
PhysiCell auto-generates dictionaries at runtime for intracellular model coupling:
- **Signals**: substrate concentrations/gradients, contact with cell types, pressure, damage
- **Behaviors**: secretion rates, cycle entry, death rates, motility parameters, adhesion affinities
- Access: `get_single_signal(pCell, "oxygen")`, `set_single_behavior(pCell, "cycle entry", 0.001)`

### Sample Project Pattern
All sample projects follow this structure:
```
sample_projects/my_project/
├── main.cpp
├── Makefile
├── custom_modules/
│   ├── custom.cpp
│   └── custom.h
└── config/
    ├── PhysiCell_settings.xml
    ├── cell_rules.csv (optional)
    └── initial_cells.csv (optional)
```
Load with: `make my-project-sample`, then `make` to compile.

## Critical Code Patterns

### DO: Initialize Interface Before Any PhysiCell Operations
```cpp
// FIRST THING in main.cpp before loading config or creating cells
BioFVM::initialize_microenvironment_interface();
```
**Why**: ALL PhysiCell components access microenvironment via `get_microenvironment_i()`. Without initialization, this returns `nullptr` and crashes.

### DO: Access Microenvironment via Interface
```cpp
// Correct: Use interface accessor
get_microenvironment_i()->simulate_diffusion_decay(dt);
get_microenvironment_i()->find_density_index("oxygen");
Voxel& v = get_microenvironment_i()->nearest_voxel(position);

// WRONG: Direct access to BioFVM (breaks abstraction)
microenvironment.simulate_diffusion_decay(dt); // Not available on this branch!
```

### DO: Update Cell Container via Interface
```cpp
// Correct pattern in time loop
((Cell_Container*)get_microenvironment_i()->get_agent_container())->update_all_cells(t);

// Cast is necessary because interface returns Basic_Agent_Container*,
// but PhysiCell uses Cell_Container (derived class)
```

### DO: Access Phenotype Parameters Correctly
```cpp
// Correct: Use accessor functions
pCell->phenotype.death.rates[apoptosis_model_index] = 0.001;
pCell->phenotype.secretion.secretion_rate("oxygen") = 10.0;

// Correct: Use dictionaries for rules-based models
set_single_behavior(pCell, "cycle entry", 0.01);
```

### DON'T: Assume Linear Data Layouts
```cpp
// WRONG: BioFVM vectors may not be contiguous
std::vector<double>& data = microenvironment.gradient_vector(voxel_idx);

// Correct: Access via proper indices/functions via interface
get_microenvironment_i()->gradient_vector(voxel_idx)[substrate_idx]
```

### DO: Follow Time Step Conventions
- `dt_diffusion` (default 0.01 min) - BioFVM solver
- `dt_mechanics` (default 0.1 min) - cell position updates
- `dt_phenotype` (default 6 min) - phenotype function calls
Set in XML `<overall>` section.

### DON'T: Modify Core Directly for Custom Behaviors
Use phenotype functions instead:
```cpp
// Correct: In create_cell_types()
pCD->functions.update_phenotype = my_custom_phenotype;
pCD->functions.custom_cell_rule = my_cell_contacts;
pCD->functions.contact_function = my_contact_function;
```

## Key Integration Points

### Adding Custom Substrates
1. Define in `<microenvironment_setup>` in XML
2. Access via `get_substrate_index("my_substrate")`
3. Use in rules: signals automatically include new substrates

### Intracellular Model Integration
PhysiCell supports three types via addons:
- **PhysiBoSS**: Boolean networks (MaBoSS .bnd/.cfg files)
- **libRoadrunner**: SBML ODE models
- **dFBA**: Flux balance analysis (SBML)

Specify in `<cell_definition><phenotype><intracellular>` XML element.

### CSV Cell Initialization
Load pre-positioned cells:
```cpp
// In setup_tissue()
load_cells_csv("./config/initial_cells.csv");
```
Format: `x,y,z,cell type` (uses cell type name, not ID).

## Common Gotchas

1. **Interface initialization**: **MUST call `BioFVM::initialize_microenvironment_interface()` as first line in `main()`** before any PhysiCell operations. Forgetting causes `nullptr` dereference crashes.
2. **Microenvironment access**: On this branch, ALWAYS use `get_microenvironment_i()->method()` instead of direct `microenvironment.method()` calls. Old PhysiCell examples won't compile.
3. **Adapter pattern understanding**: Interface (`Bio_microenvironment_interface.h`) defines contract, Adapter (`microenvironment_adapter.*`) implements by wrapping BioFVM. To add alternative solvers, implement new adapter—don't modify interface.
4. **Random seed**: Set `<random_seed>` in XML `<options>` section for reproducibility. Wrote to `output/random_seed.txt`.
5. **Neighbor list updates**: `state.neighbors` auto-updated by mechanics, but `state.attached_cells` requires manual management.
6. **Cell volume**: Use `set_target_volume()` or `set_target_radius()` to change size while preserving nuclear:cytoplasmic ratio.
7. **Death models**: Use `pCell->start_death(death_model_index)`, NOT `phenotype.death.trigger_death(index)`.
8. **PhysiBoSS projects**: Always enable `-DPHYSICELL_BUILD_PHYSIBOSS=ON` in CMake or use appropriate Makefile rules.
9. **Dirichlet conditions**: Per-boundary, per-substrate control via `<Dirichlet_options>` in XML.

## Documentation & Resources
- Primary docs: `./documentation-deprecated/` (being updated)
- Quick reference: `README.md` (comprehensive changelog)
- Sample projects: Best source of working examples
- Community: PhysiCell Slack (see README for invite)
- Key files to understand: `PhysiCell_cell.h`, `PhysiCell_phenotype.h`, `PhysiCell_settings.cpp`

## Version-Specific Notes
- **v1.14.x**: Cell Behavior Hypothesis Grammar v3, `Cell_Integrity` for damage modeling, asymmetric division
- **v1.10.x**: Introduced cell interactions (phagocytosis/attack/fusion), transformations, signal/behavior dictionaries
- **v1.7.x**: XML-based cell definitions became standard (prefer over C++ definitions)
