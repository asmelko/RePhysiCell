/*
#############################################################################
# If you use BioFVM in your project, please cite BioFVM and the version     #
# number, such as below:                                                    #
#                                                                           #
# We solved the diffusion equations using BioFVM (Version 1.1.7) [1]        #
#                                                                           #
# [1] A. Ghaffarizadeh, S.H. Friedman, and P. Macklin, BioFVM: an efficient #
#    parallelized diffusive transport solver for 3-D biological simulations,#
#    Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730 #
#                                                                           #
#############################################################################
#                                                                           #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)   #
#                                                                           #
# Copyright (c) 2015-2017, Paul Macklin and the BioFVM Project              #
# All rights reserved.                                                      #
#                                                                           #
# Redistribution and use in source and binary forms, with or without        #
# modification, are permitted provided that the following conditions are    #
# met:                                                                      #
#                                                                           #
# 1. Redistributions of source code must retain the above copyright notice, #
# this list of conditions and the following disclaimer.                     #
#                                                                           #
# 2. Redistributions in binary form must reproduce the above copyright      #
# notice, this list of conditions and the following disclaimer in the       #
# documentation and/or other materials provided with the distribution.      #
#                                                                           #
# 3. Neither the name of the copyright holder nor the names of its          #
# contributors may be used to endorse or promote products derived from this #
# software without specific prior written permission.                       #
#                                                                           #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       #
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED #
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A           #
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER #
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  #
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       #
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        #
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    #
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      #
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        #
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              #
#                                                                           #
#############################################################################
*/

#include "BioFVM_solvers.h" 
#include "BioFVM_vector.h" 

#include <iostream>
#include <omp.h>

namespace BioFVM{

// do I even need this? 
void diffusion_decay_solver__constant_coefficients_explicit( Microenvironment& M, double dt )
{
	static bool precomputations_and_constants_done = false; 
	if( !precomputations_and_constants_done )
	{
		std::cout	<< std::endl << "Using solver: " << __FUNCTION__ << std::endl 
					<< "     (constant diffusion coefficient with explicit stepping, implicit decay) ... " << std::endl << std::endl;  

		if( M.mesh.uniform_mesh == true )
		{
			std::cout << "Uniform mesh detected! Consider switching to a more efficient method, such as " << std::endl  
			<< "     diffusion_decay_solver__constant_coefficients_explicit_uniform_mesh" << std::endl  
			<< std::endl; 
		}

		precomputations_and_constants_done = true; 
	}

	return; 
}

void diffusion_decay_solver__constant_coefficients_explicit_uniform_mesh( Microenvironment& M, double dt )
{
	static bool precomputations_and_constants_done = false; 
	if( !precomputations_and_constants_done )
	{
		std::cout	<< std::endl << "Using solver: " << __FUNCTION__ << std::endl 
					<< "     (constant diffusion coefficient with explicit stepping, implicit decay, uniform mesh) ... " << std::endl << std::endl;  

		if( M.mesh.regular_mesh == false )
		{ std::cout << "Error. This code is only supported for regular meshes." << std::endl; }

		precomputations_and_constants_done = true; 
	}

	return; 
}

void diffusion_decay_solver__constant_coefficients_LOD_3D( Microenvironment& M, double dt )
{
	if( M.mesh.regular_mesh == false || M.mesh.Cartesian_mesh == false )
	{
		std::cout << "Error: This algorithm is written for regular Cartesian meshes. Try: other solvers!" << std::endl << std::endl; 
	return; 
	}

	// define constants and pre-computed quantities 
	
	if( !M.diffusion_solver_setup_done )
	{
		std::cout << std::endl << "Using method " << __FUNCTION__ << " (implicit 3-D LOD with Thomas Algorithm) ... " 
		<< std::endl << std::endl;  

		M.setup_thomas_constants(dt, 3);
	}

	unsigned int d_dim = M.number_of_densities();

	// x-diffusion 
	
	M.apply_dirichlet_conditions();
	#pragma omp parallel for 
	for( unsigned int k=0; k < M.mesh.z_coordinates.size() ; k++ )
	{
		for( unsigned int j=0; j < M.mesh.y_coordinates.size() ; j++ )
		{
			// Thomas solver, x-direction

			// remaining part of forward elimination, using pre-computed quantities 
			int n = M.voxel_index(0,j,k);
			for (int d = 0; d < d_dim; d++)
				M.p_density_vectors[n * d_dim + d] /= M.thomas_denom[d]; 

			for( unsigned int i=1; i < M.mesh.x_coordinates.size() ; i++ )
			{
				n = M.voxel_index(i,j,k); 

				for (int d = 0; d < d_dim; d++)
					M.p_density_vectors[n * d_dim + d] = 
						(M.p_density_vectors[n * d_dim + d] + 
						 M.thomas_constant1[d] * M.p_density_vectors[(n-M.thomas_i_jump) * d_dim + d]) /
						 M.thomas_denom[i* d_dim + d];
			}

			for( int i = M.mesh.x_coordinates.size()-2 ; i >= 0 ; i-- )
			{
				n = M.voxel_index(i,j,k); 
				for (int d = 0; d < d_dim; d++)
					M.p_density_vectors[n * d_dim + d] += 
						M.thomas_c[i* d_dim + d] * M.p_density_vectors[(n+M.thomas_i_jump) * d_dim + d];
			}

		}
	}

	// y-diffusion 

	M.apply_dirichlet_conditions();
	#pragma omp parallel for 
	for( unsigned int k=0; k < M.mesh.z_coordinates.size() ; k++ )
	{
		for( unsigned int i=0; i < M.mesh.x_coordinates.size() ; i++ )
		{
			// Thomas solver, y-direction

			// remaining part of forward elimination, using pre-computed quantities 
			int n = M.voxel_index(i,0,k);
			for (int d = 0; d < d_dim; d++)
				M.p_density_vectors[n * d_dim + d] /= M.thomas_denom[d]; 

			for( unsigned int j=1; j < M.mesh.y_coordinates.size() ; j++ )
			{
				n = M.voxel_index(i,j,k); 

				for (int d = 0; d < d_dim; d++)
					M.p_density_vectors[n * d_dim + d] = 
						(M.p_density_vectors[n * d_dim + d] + 
						 M.thomas_constant1[d] * M.p_density_vectors[(n-M.thomas_j_jump) * d_dim + d]) /
						 M.thomas_denom[j* d_dim + d];
			}

			for( int j = M.mesh.y_coordinates.size()-2 ; j >= 0 ; j-- )
			{
				n = M.voxel_index(i,j,k); 
				for (int d = 0; d < d_dim; d++)
					M.p_density_vectors[n * d_dim + d] += 
						M.thomas_c[j* d_dim + d] * M.p_density_vectors[(n+M.thomas_j_jump) * d_dim + d];
			}

		}
	}

 // z-diffusion 

	M.apply_dirichlet_conditions();
	#pragma omp parallel for 
	for( unsigned int j=0; j < M.mesh.y_coordinates.size() ; j++ )
	{
		for( unsigned int i=0; i < M.mesh.x_coordinates.size() ; i++ )
		{
			// Thomas solver, z-direction

			// remaining part of forward elimination, using pre-computed quantities 
			int n = M.voxel_index(i,j,0);
			for (int d = 0; d < d_dim; d++)
				M.p_density_vectors[n * d_dim + d] /= M.thomas_denom[d]; 

			for( unsigned int k=1; k < M.mesh.z_coordinates.size() ; k++ )
			{
				n = M.voxel_index(i,j,k); 

				for (int d = 0; d < d_dim; d++)
					M.p_density_vectors[n * d_dim + d] = 
						(M.p_density_vectors[n * d_dim + d] + 
						 M.thomas_constant1[d] * M.p_density_vectors[(n-M.thomas_k_jump) * d_dim + d]) /
						 M.thomas_denom[k* d_dim + d];
			}

			for( int k = M.mesh.z_coordinates.size()-2 ; k >= 0 ; k-- )
			{
				n = M.voxel_index(i,j,k); 
				for (int d = 0; d < d_dim; d++)
					M.p_density_vectors[n * d_dim + d] += 
						M.thomas_c[k* d_dim + d] * M.p_density_vectors[(n+M.thomas_k_jump) * d_dim + d];
			}

		}
	}
 
	M.apply_dirichlet_conditions();
	
	// reset gradient vectors 
//	M.reset_all_gradient_vectors(); 

	return; 
}

void diffusion_decay_solver__constant_coefficients_LOD_2D( Microenvironment& M, double dt )
{
	if( M.mesh.regular_mesh == false )
	{
		std::cout << "Error: This algorithm is written for regular Cartesian meshes. Try: something else." << std::endl << std::endl; 
		return; 
	}
	
	// constants for the linear solver (Thomas algorithm) 
	
	if( !M.diffusion_solver_setup_done )
	{
		std::cout << std::endl << "Using method " << __FUNCTION__ << " (2D LOD with Thomas Algorithm) ... " << std::endl << std::endl;  
		
		M.setup_thomas_constants(dt, 2);
	}

	// set the pointer

	unsigned int d_dim = M.number_of_densities();
	
	M.apply_dirichlet_conditions();

	// x-diffusion 
	#pragma omp parallel for 
	for( unsigned int j=0; j < M.mesh.y_coordinates.size() ; j++ )
	{
		// Thomas solver, x-direction

		// remaining part of forward elimination, using pre-computed quantities 
		unsigned int n = M.voxel_index(0,j,0);
		for (int d = 0; d < d_dim; d++)
			M.p_density_vectors[n * d_dim + d] /= M.thomas_denom[d]; 

		n += M.thomas_i_jump; 
		for( unsigned int i=1; i < M.mesh.x_coordinates.size() ; i++ )
		{
			for (int d = 0; d < d_dim; d++)
				M.p_density_vectors[n * d_dim + d] = 
					(M.p_density_vectors[n * d_dim + d] + 
						M.thomas_constant1[d] * M.p_density_vectors[(n-M.thomas_i_jump) * d_dim + d]) /
						M.thomas_denom[i* d_dim + d];
			n += M.thomas_i_jump; 
		}

		// back substitution 
		n = M.voxel_index( M.mesh.x_coordinates.size()-2 ,j,0); 

		for( int i = M.mesh.x_coordinates.size()-2 ; i >= 0 ; i-- )
		{
			for (int d = 0; d < d_dim; d++)
				M.p_density_vectors[n * d_dim + d] += 
					M.thomas_c[i* d_dim + d] * M.p_density_vectors[(n+M.thomas_i_jump) * d_dim + d];
			n -= M.thomas_i_jump; 
		}
	}

	// y-diffusion 

	M.apply_dirichlet_conditions();
	#pragma omp parallel for 
	for( unsigned int i=0; i < M.mesh.x_coordinates.size() ; i++ )
	{
		// Thomas solver, y-direction

		// remaining part of forward elimination, using pre-computed quantities 

		int n = M.voxel_index(i,0,0);
		for (int d = 0; d < d_dim; d++)
			M.p_density_vectors[n * d_dim + d] /= M.thomas_denom[d]; 

		n += M.thomas_j_jump; 
		for( unsigned int j=1; j < M.mesh.y_coordinates.size() ; j++ )
		{
			for (int d = 0; d < d_dim; d++)
				M.p_density_vectors[n * d_dim + d] = 
					(M.p_density_vectors[n * d_dim + d] + 
						M.thomas_constant1[d] * M.p_density_vectors[(n-M.thomas_j_jump) * d_dim + d]) /
						M.thomas_denom[j* d_dim + d];
			n += M.thomas_j_jump; 
		}

		// back substitution 
		n = M.voxel_index( i,M.mesh.y_coordinates.size()-2, 0); 

		for( int j = M.mesh.y_coordinates.size()-2 ; j >= 0 ; j-- )
		{
			for (int d = 0; d < d_dim; d++)
				M.p_density_vectors[n * d_dim + d] += 
					M.thomas_c[j* d_dim + d] * M.p_density_vectors[(n+M.thomas_j_jump) * d_dim + d];
			n -= M.thomas_j_jump; 
		}
	}

	M.apply_dirichlet_conditions();
	
	// reset gradient vectors 
//	M.reset_all_gradient_vectors(); 
	
	return; 
}

void diffusion_decay_explicit_uniform_rates( Microenvironment& M, double dt )
{
	using std::vector; 
	using std::cout; 
	using std::endl; 

	// static int n_jump_i = 1; 
	// static int n_jump_j = M.mesh.x_coordinates.size(); 
	// static int n_jump_k = M.mesh.x_coordinates.size() * M.mesh.y_coordinates.size(); 

	if( !M.diffusion_solver_setup_done )
	{	
		M.thomas_i_jump = 1; 
		M.thomas_j_jump = M.mesh.x_coordinates.size(); 
		M.thomas_k_jump = M.thomas_j_jump * M.mesh.y_coordinates.size(); 
	
		M.diffusion_solver_setup_done = true; 
	}
	
	if( M.mesh.uniform_mesh == false )
	{
		cout << "Error: This algorithm is written for uniform Cartesian meshes. Try: something else" << endl << endl; 
		return; 
	}

	// double buffering to reduce memory copy / allocation overhead 

	double* pNew = M.temporary_density_vectors1.data();
	double* pOld = M.temporary_density_vectors2.data();

	// swap the buffers 

	std::swap(pNew, pOld);
	M.p_density_vectors = pNew; 

	// static bool reaction_diffusion_shortcuts_are_set = false; 

	static vector<double> constant1 = (1.0 / ( M.mesh.dx * M.mesh.dx )) * M.diffusion_coefficients; 
	static vector<double> constant2 = dt * constant1; 
	static vector<double> constant3 = M.one + dt * M.decay_rates;

	static vector<double> constant4 = M.one - dt * M.decay_rates;

	unsigned int d_dim = M.number_of_densities();

	#pragma omp parallel for
	for( unsigned int i=0; i < M.number_of_voxels() ; i++ )
	{
		unsigned int number_of_neighbors = M.mesh.connected_voxel_indices[i].size(); 

		double d1 = -1.0 * number_of_neighbors; 

		for (int d = 0; d < d_dim; d++)
			pNew[i * d_dim + d] = pOld[i * d_dim + d] * constant4[d];

		for( unsigned int j=0; j < number_of_neighbors ; j++ )
		{
			int neighbor_idx = M.mesh.connected_voxel_indices[i][j];

			for (int d = 0; d < d_dim; d++)
				pNew[i * d_dim + d] += 
					constant2[d] * pOld[neighbor_idx * d_dim + d];
		}
		
		for (int d = 0; d < d_dim; d++)
			pNew[i * d_dim + d] += 
				constant2[d] * d1 * pOld[i * d_dim + d];
	}
	
	// reset gradient vectors 
//	M.reset_all_gradient_vectors(); 

	return; 
}

void diffusion_decay_solver__constant_coefficients_LOD_1D( Microenvironment& M, double dt )
{
	if( M.mesh.regular_mesh == false )
	{
		std::cout << "Error: This algorithm is written for regular Cartesian meshes. Try: something else." << std::endl << std::endl; 
		return; 
	}
	
	// constants for the linear solver (Thomas algorithm) 
	
	if( !M.diffusion_solver_setup_done )
	{
		std::cout << std::endl << "Using method " << __FUNCTION__ << " (2D LOD with Thomas Algorithm) ... " << std::endl << std::endl;  
		
		M.setup_thomas_constants(dt, 1);
	}

	unsigned int d_dim = M.number_of_densities();

	// set the pointer
	
	M.apply_dirichlet_conditions();

	// x-diffusion 
	#pragma omp parallel for 
	for( unsigned int j=0; j < M.mesh.y_coordinates.size() ; j++ )
	{
		// Thomas solver, x-direction

		// remaining part of forward elimination, using pre-computed quantities 
		unsigned int n = M.voxel_index(0,j,0);
		for (int d = 0; d < d_dim; d++)
			M.p_density_vectors[n * d_dim + d] /= M.thomas_denom[d]; 

		n += M.thomas_i_jump; 
		for( unsigned int i=1; i < M.mesh.x_coordinates.size() ; i++ )
		{
			for (int d = 0; d < d_dim; d++)
				M.p_density_vectors[n * d_dim + d] = 
					(M.p_density_vectors[n * d_dim + d] + 
						M.thomas_constant1[d] * M.p_density_vectors[(n-M.thomas_i_jump) * d_dim + d]) /
						M.thomas_denom[i* d_dim + d];
			n += M.thomas_i_jump; 
		}

		// back substitution 
		n = M.voxel_index( M.mesh.x_coordinates.size()-2 ,j,0); 

		for( int i = M.mesh.x_coordinates.size()-2 ; i >= 0 ; i-- )
		{
			for (int d = 0; d < d_dim; d++)
				M.p_density_vectors[n * d_dim + d] += 
					M.thomas_c[i* d_dim + d] * M.p_density_vectors[(n+M.thomas_i_jump) * d_dim + d];
			n -= M.thomas_i_jump; 
		}
	}

	M.apply_dirichlet_conditions();
	
	// reset gradient vectors 
//	M.reset_all_gradient_vectors(); 
	
	return; 
}


};
