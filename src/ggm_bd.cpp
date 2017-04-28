// ----------------------------------------------------------------------------|
//     Copyright (C) 2012-2016 Mohammadi A. and Wit C. E.
//
//     This file is part of BDgraph package.
//
//     BDgraph is free software: you can redistribute it and/or modify it under 
//     the terms of the GNU General Public License as published by the Free 
//     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.
//
//     Maintainer:
//     Abdolreza Mohammadi: a.mohammadi@rug.nl or a.mohammadi@uvt.nl
// ----------------------------------------------------------------------------|
#include <sstream>
#include <string>        // std::string, std::to_string
#include <vector>        // for using vector
#include <math.h>        // isinf, sqrt
#include <R.h>
#include <Rmath.h>
#include <algorithm>     // for transform function
#include <functional>    // for transform function

#include "matrix.h"
#include "rgwish.h"

using namespace std;

extern "C" {
// ----------------------------------------------------------------------------|
// birth-death MCMC for Gaussian Graphical models  
// for case D = I_p 
// it is for Bayesian model averaging (MA)
// ----------------------------------------------------------------------------|
void ggm_bdmcmc_ma( int *iter, int *burnin, int G[], double Ts[], double K[], int *p, 
			 double K_hat[], double p_links[],
			 int *b, int *b_star, double Ds[], int *print )
{
	int print_c = *print, iteration = *iter, burn_in = *burnin;
	int index_selected_edge, selected_edge_i, selected_edge_j, selected_edge_ij;
	int i, j, ij, one = 1, dim = *p, pxp = dim * dim, p1 = dim - 1, p2x2 = ( dim - 2 ) * 2;

	double Dsij, sum_weights = 0.0, weight_C, sum_rates;
	
	vector<double> p_links_Cpp( pxp, 0.0 ); 
	vector<double> K_hat_Cpp( pxp, 0.0 ); 

	//-- allocation for birth/death rates 
	vector<double> K121( 4 ); 
	vector<double> Kj12( p1 );              // K[j, -j]
	vector<double> sigmaj12( p1 );          // sigma[-j, j]  
	vector<double> sigmaj22( p1 * p1 );     // sigma[-j, -j]
	vector<double> Kj12xK22_inv( p1 ); 

	vector<double> K21( p2x2 );                          // K[-e, e]
	vector<double> sigma21( p2x2 );                      // sigma[-e, e]
	vector<double> sigma22( ( dim - 2 ) * ( dim - 2 ) ); // sigma[-e, -e]
	vector<double> sigma11_inv( 4 );                     // inv( sigma[e, e] )
	vector<double> sigma21xsigma11_inv( p2x2 ); 
	vector<double> K12xK22_inv( p2x2 );   
	// -- allocation for rgwish_sigma 
	vector<double> sigma_start( pxp ); 
	vector<double> inv_C( pxp ); 
	vector<double> beta_star( dim ); 
	vector<double> sigma_i( dim ); 
	vector<double> sigma_start_N_i( dim );   // For dynamic memory used
	vector<double> sigma_N_i( pxp );         // For dynamic memory used
	vector<int> N_i( dim );                  // For dynamic memory used
	// ----------------------------
	
	vector<double> sigma( pxp ); 
	vector<double> copyK( pxp ); 
	memcpy( &copyK[0], K, sizeof( double ) * pxp );
	inverse( &copyK[0], &sigma[0], &dim );			

	// Counting size of notes
	int ip;
	vector<int> size_node( dim, 0 );
	for( i = 0; i < dim; i++ )
	{
		ip = i * dim;
		for( j = 0; j < dim; j++ ) size_node[i] += G[ip + j];
	}
	
	// For finding the index of rates 
	int qp = dim * ( dim - 1 ) / 2;
	vector<double> rates( qp );
	vector<int> index_rates_row( qp );
	vector<int> index_rates_col( qp );
	int counter = 0;
	vector<double> Dsijj( pxp ); 
	for( j = 1; j < dim; j++ )
		for( i = 0; i < j; i++ )
		{
			index_rates_row[counter] = i;
			index_rates_col[counter] = j;
			counter++;
			
			// for calculating the birth/death rates
			ij        = j * dim + i;
			Dsij      = Ds[ij];
			Dsijj[ij] = Dsij * Dsij / Ds[j * dim + j]; 
		}

//-- Main loop for birth-death MCMC -------------------------------------------| 
	GetRNGstate();
	for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc++ )
	{
		if( ( i_mcmc + 1 ) % print_c == 0 ) Rprintf( " Iteration  %d                 \n", i_mcmc + 1 ); 
		
//----- STEP 1: calculating birth and death rates -----------------------------|		

		rates_bdmcmc( &rates[0], G, Ds, &Dsijj[0],
                      &sigma[0], &sigma21[0], &sigma22[0], &sigmaj12[0], &sigmaj22[0],    
                      &K[0], &K21[0], &K121[0], &Kj12[0], 
                      &K12xK22_inv[0], &Kj12xK22_inv[0], &sigma11_inv[0], &sigma21xsigma11_inv[0],  
                      b, &dim );

		// Selecting an edge based on birth and death rates
		select_edge( &rates[0], &index_selected_edge, &sum_rates, &qp );
		selected_edge_i = index_rates_row[ index_selected_edge ];
		selected_edge_j = index_rates_col[ index_selected_edge ];

//----- saving result ---------------------------------------------------------|	
		if( i_mcmc >= burn_in )
		{
			weight_C = 1.0 / sum_rates;
			
			// K_hat_Cpp[i] += K[i] / sum_rates;
			F77_NAME(daxpy)( &pxp, &weight_C, &K[0], &one, &K_hat_Cpp[0], &one );
			
			for( i = 0; i < pxp ; i++ )
				if( G[i] ) p_links_Cpp[i] += weight_C;
			
			sum_weights += weight_C;
		} 
//----- End of saving result --------------------------------------------------|	
			
		// Updating G (graph) based on selected edge
		selected_edge_ij    = selected_edge_j * dim + selected_edge_i;
		G[selected_edge_ij] = 1 - G[selected_edge_ij];
		G[selected_edge_i * dim + selected_edge_j] = G[selected_edge_ij];

		if( G[selected_edge_ij] )
		{ 
			++size_node[selected_edge_i]; 
			++size_node[selected_edge_j]; 
		}
		else
		{ 
			--size_node[selected_edge_i]; 
			--size_node[selected_edge_j]; 
		}

//------ STEP 2: Sampling from G-Wishart for new graph ------------------------|
		rgwish_sigma( G, &size_node[0], Ts, K, &sigma[0], b_star, &dim, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		
	}  
	PutRNGstate();
//-- End of main loop for birth-death MCMC ------------------------------------| 

	for( i = 0; i < pxp; i++ )
	{	
		p_links[i] = p_links_Cpp[i] / sum_weights;
		K_hat[i]   = K_hat_Cpp[i] / sum_weights;
	}
}
       
// ----------------------------------------------------------------------------|
// birth-death MCMC for Gaussian Graphical models  
// for case D = I_p 
// it is for maximum a posterior probability estimation (MAP)
// ----------------------------------------------------------------------------|
void ggm_bdmcmc_map( int *iter, int *burnin, int G[], double Ts[], double K[], int *p, 
			 int all_graphs[], double all_weights[], double K_hat[], 
			 char *sample_graphs[], double graph_weights[], int *size_sample_g,
			 int *b, int *b_star, double Ds[], int *print )
{
	int print_c = *print, iteration = *iter, burn_in = *burnin, count_all_g = 0;
	int index_selected_edge, selected_edge_i, selected_edge_j, selected_edge_ij, size_sample_graph = *size_sample_g;
	int i, j, ij, one = 1, dim = *p, pxp = dim * dim, p1 = dim - 1, p2x2 = ( dim - 2 ) * 2;

	double Dsij, sum_weights = 0.0, weight_C, sum_rates;
	
	bool this_one;

	string string_g;
	vector<string> sample_graphs_C( iteration - burn_in );

	//-- allocation for birth/death rates 
	vector<double> K121( 4 ); 
	vector<double> Kj12( p1 );              // K[j, -j]
	vector<double> sigmaj12( p1 );          // sigma[-j, j]  
	vector<double> sigmaj22( p1 * p1 );     // sigma[-j, -j]
	vector<double> Kj12xK22_inv( p1 ); 

	vector<double> K21( p2x2 );                          // K[-e, e]
	vector<double> sigma21( p2x2 );                      // sigma[-e, e]
	vector<double> sigma22( ( dim - 2 ) * ( dim - 2 ) ); // sigma[-e, -e]
	vector<double> sigma11_inv( 4 );                     // inv( sigma[e, e] )
	vector<double> sigma21xsigma11_inv( p2x2 ); 
	vector<double> K12xK22_inv( p2x2 );   
	// -- allocation for rgwish_sigma 
	vector<double> sigma_start( pxp ); 
	vector<double> inv_C( pxp ); 
	vector<double> beta_star( dim ); 
	vector<double> sigma_i( dim ); 
	vector<double> sigma_start_N_i( dim );   // For dynamic memory used
	vector<double> sigma_N_i( pxp );         // For dynamic memory used
	vector<int> N_i( dim );                  // For dynamic memory used
	// ----------------------------

	vector<double> sigma( pxp ); 
	vector<double> copyK( pxp ); 
	memcpy( &copyK[0], K, sizeof( double ) * pxp );
	inverse( &copyK[0], &sigma[0], &dim );			

	int qp = dim * ( dim - 1 ) / 2;
	vector<char> char_g( qp );              // char string_g[pp];

	// Counting size of notes
	int ip;
	vector<int> size_node( dim, 0 );
	for( i = 0; i < dim; i++ )
	{
		ip = i * dim;
		for( j = 0; j < dim; j++ ) size_node[i] += G[ip + j];
	}
	
	// For finding the index of rates 
	vector<double> rates( qp );
	vector<int> index_rates_row( qp );
	vector<int> index_rates_col( qp );
	int counter = 0;
	vector<double> Dsijj( pxp ); 
	for( j = 1; j < dim; j++ )
		for( i = 0; i < j; i++ )
		{
			index_rates_row[counter] = i;
			index_rates_col[counter] = j;
			counter++;
			
			// for calculating the birth/death rates
			ij        = j * dim + i;
			Dsij      = Ds[ij];
			Dsijj[ij] = Dsij * Dsij / Ds[j * dim + j]; 
		}

//-- Main loop for birth-death MCMC -------------------------------------------| 
	GetRNGstate();
	for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc++ )
	{
		if( ( i_mcmc + 1 ) % print_c == 0 ) Rprintf( " Iteration  %d                 \n", i_mcmc + 1 ); 
		
//----- STEP 1: calculating birth and death rates -----------------------------|		

		rates_bdmcmc( &rates[0], G, Ds, &Dsijj[0],
                      &sigma[0], &sigma21[0], &sigma22[0], &sigmaj12[0], &sigmaj22[0],    
                      &K[0], &K21[0], &K121[0], &Kj12[0], 
                      &K12xK22_inv[0], &Kj12xK22_inv[0], &sigma11_inv[0], &sigma21xsigma11_inv[0],  
                      b, &dim );
		
		// Selecting an edge based on birth and death rates
		select_edge( &rates[0], &index_selected_edge, &sum_rates, &qp );
		selected_edge_i = index_rates_row[ index_selected_edge ];
		selected_edge_j = index_rates_col[ index_selected_edge ];

//----- saving result ---------------------------------------------------------|	
		counter = 0;	
		for( j = 1; j < dim; j++ )
			for( i = 0; i < j; i++ )
				char_g[counter++] = G[ j * dim + i ] + '0'; 

		if( i_mcmc >= burn_in )
		{
			weight_C = 1.0 / sum_rates;
			
			//for( i = 0; i < pxp; i++ ) K_hat[i] += K[i] / sum_rates;
			F77_NAME(daxpy)( &pxp, &weight_C, &K[0], &one, &K_hat[0], &one );			

			string_g = string( char_g.begin(), char_g.end() );	
			all_weights[count_all_g] = weight_C;
			
			this_one = false;
			for( i = 0; i < size_sample_graph; i++ )
				if( sample_graphs_C[i] == string_g )
				{
					graph_weights[i] += all_weights[count_all_g];
					all_graphs[count_all_g] = i;
					this_one = true;
					break;
				} 
			
			if( !this_one || size_sample_graph == 0 )
			{
				sample_graphs_C[size_sample_graph] = string_g;
				graph_weights[size_sample_graph]   = all_weights[count_all_g];
				all_graphs[count_all_g]            = size_sample_graph; 
				size_sample_graph++;				
			}
			
			count_all_g++; 
			sum_weights += weight_C;
		} 
//----- End of saving result --------------------------------------------------|	
    			
		// Updating G (graph) based on selected edge
		selected_edge_ij    = selected_edge_j * dim + selected_edge_i;
		G[selected_edge_ij] = 1 - G[selected_edge_ij];
		G[selected_edge_i * dim + selected_edge_j] = G[selected_edge_ij];

		if( G[selected_edge_ij] )
		{ 
			++size_node[selected_edge_i]; 
			++size_node[selected_edge_j]; 
		}
		else
		{ 
			--size_node[selected_edge_i]; 
			--size_node[selected_edge_j]; 
		}

//------ STEP 2: Sampling from G-Wishart for new graph ------------------------|
		rgwish_sigma( G, &size_node[0], Ts, K, &sigma[0], b_star, &dim, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		
	}  
	PutRNGstate();
//-- End of main loop for birth-death MCMC ------------------------------------| 

	for( i = 0; i < size_sample_graph; i++ ) 
	{
		sample_graphs_C[i].copy( sample_graphs[i], qp, 0 );
		sample_graphs[i][qp] = '\0';
	}
	
	*size_sample_g = size_sample_graph;
	
	for( i = 0; i < pxp; i++ )
		K_hat[i] /= sum_weights;
}
        
// ----------------------------------------------------------------------------|
// Multiple birth-death MCMC for Gaussian Graphical models  
// for D = I_p 
// it is for Bayesian model averaging
// ----------------------------------------------------------------------------|
void ggm_bdmcmc_ma_multi_update( int *iter, int *burnin, int G[], double Ts[], double K[], int *p, 
			 double K_hat[], double p_links[],
			 int *b, int *b_star, double Ds[], int *multi_update, int *print )
{
	int print_c = *print, iteration = *iter, burn_in = *burnin, multi_update_C = *multi_update;
	int selected_edge_i, selected_edge_j, selected_edge_ij;
	int i, j,ij, one = 1, dim = *p, pxp = dim * dim, p1 = dim - 1, p2x2 = ( dim - 2 ) * 2;

	double Dsij, sum_weights = 0.0, weight_C, sum_rates;

	vector<double> p_links_Cpp( pxp, 0.0 ); 
	vector<double> K_hat_Cpp( pxp, 0.0 ); 

	//-- allocation for birth/death rates 
	vector<double> K121( 4 ); 
	vector<double> Kj12( p1 );              // K[j, -j]
	vector<double> sigmaj12( p1 );          // sigma[-j, j]  
	vector<double> sigmaj22( p1 * p1 );     // sigma[-j, -j]
	vector<double> Kj12xK22_inv( p1 ); 

	vector<double> K21( p2x2 );                          // K[-e, e]
	vector<double> sigma21( p2x2 );                      // sigma[-e, e]
	vector<double> sigma22( ( dim - 2 ) * ( dim - 2 ) ); // sigma[-e, -e]
	vector<double> sigma11_inv( 4 );                     // inv( sigma[e, e] )
	vector<double> sigma21xsigma11_inv( p2x2 ); 
	vector<double> K12xK22_inv( p2x2 );   
	// -- allocation for rgwish_sigma 
	vector<double> sigma_start( pxp ); 
	vector<double> inv_C( pxp ); 
	vector<double> beta_star( dim ); 
	vector<double> sigma_i( dim ); 
	vector<double> sigma_start_N_i( dim );   // For dynamic memory used
	vector<double> sigma_N_i( pxp );         // For dynamic memory used
	vector<int> N_i( dim );                  // For dynamic memory used
	// ----------------------------

	vector<double> sigma( pxp ); 
	vector<double> copyK( pxp ); 
	memcpy( &copyK[0], K, sizeof( double ) * pxp );
	inverse( &copyK[0], &sigma[0], &dim );			

	// Count size of notes
	int ip;
	vector<int> size_node( dim, 0 );
	for( i = 0; i < dim; i++ )
	{
		ip = i * dim;
		for( j = 0; j < dim; j++ ) size_node[i] += G[ip + j];
	}

	// For finding the index of rates 
	int qp = dim * ( dim - 1 ) / 2;
	vector<double> rates( qp );
	vector<int> index_rates_row( qp );
	vector<int> index_rates_col( qp );
	int counter = 0;
	vector<double> Dsijj( pxp ); 
	for( j = 1; j < dim; j++ )
		for( i = 0; i < j; i++ )
		{
			index_rates_row[counter] = i;
			index_rates_col[counter] = j;
			counter++;
			
			// for calculating the birth/death rates
			ij        = j * dim + i;
			Dsij      = Ds[ij];
			Dsijj[ij] = Dsij * Dsij / Ds[j * dim + j]; 
		}

	int size_index = multi_update_C;
	vector<int> index_selected_edges( multi_update_C );

//-- Main loop for birth-death MCMC -------------------------------------------| 
	GetRNGstate();
	for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc += size_index )
	{
		if( ( i_mcmc + 1 ) % print_c < multi_update_C ) Rprintf( " Iteration  %d                 \n", i_mcmc + 1 ); 
		
//----- STEP 1: calculating birth and death rates -----------------------------|		

		rates_bdmcmc( &rates[0], G, Ds, &Dsijj[0],
                      &sigma[0], &sigma21[0], &sigma22[0], &sigmaj12[0], &sigmaj22[0],    
                      &K[0], &K21[0], &K121[0], &Kj12[0], 
                      &K12xK22_inv[0], &Kj12xK22_inv[0], &sigma11_inv[0], &sigma21xsigma11_inv[0],  
                      b, &dim );

		// Selecting multiple edges based on birth and death rates
		select_multi_edges( &rates[0], &index_selected_edges[0], &size_index, &sum_rates, &multi_update_C, &qp );

//----- Saving result ---------------------------------------------------------|	
		if( i_mcmc >= burn_in )
		{
			weight_C = 1.0 / sum_rates;
			
			// K_hat_Cpp[i] += K[i] / sum_rates;
			F77_NAME(daxpy)( &pxp, &weight_C, &K[0], &one, &K_hat_Cpp[0], &one );
			
			for( i = 0; i < pxp ; i++ )
				if( G[i] ) p_links_Cpp[i] += weight_C;
			
			sum_weights += weight_C;
		} 
//----- End of saving result --------------------------------------------------|	

		// Updating graph based on selected edges
		for ( i = 0; i < size_index; i++ )
		{
			selected_edge_i = index_rates_row[ index_selected_edges[i] ];
			selected_edge_j = index_rates_col[ index_selected_edges[i] ];
			
			selected_edge_ij    = selected_edge_j * dim + selected_edge_i;
			G[selected_edge_ij] = 1 - G[selected_edge_ij];
			G[selected_edge_i * dim + selected_edge_j] = G[selected_edge_ij];
		
			if( G[selected_edge_ij] )
			{ 
				++size_node[selected_edge_i]; 
				++size_node[selected_edge_j]; 
			}
			else
			{ 
				--size_node[selected_edge_i]; 
				--size_node[selected_edge_j]; 
			}		
		}

//------ STEP 2: Sampling from G-Wishart for new graph ------------------------|
		rgwish_sigma( G, &size_node[0], Ts, K, &sigma[0], b_star, &dim, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		
	}  
	PutRNGstate();
//-- End of main loop for birth-death MCMC ------------------------------------| 

	for( i = 0; i < pxp; i++ )
	{	
		p_links[i] = p_links_Cpp[i] / sum_weights;
		K_hat[i]   = K_hat_Cpp[i] / sum_weights;
	}
}
    
// ----------------------------------------------------------------------------|
// Multiple birth-death MCMC for Gaussian Graphical models  
// for D = I_p 
// it is for maximum a posterior probability estimation (MAP)
// ----------------------------------------------------------------------------|
void ggm_bdmcmc_map_multi_update( int *iter, int *burnin, int G[], double Ts[], double K[], int *p, 
			 int all_graphs[], double all_weights[], double K_hat[], 
			 char *sample_graphs[], double graph_weights[], int *size_sample_g, int *counter_all_g,
			 int *b, int *b_star, double Ds[], int *multi_update, int *print )
{
	int print_c = *print, multi_update_C = *multi_update, iteration = *iter, burn_in = *burnin;
	int count_all_g = *counter_all_g, selected_edge_i, selected_edge_j, selected_edge_ij, size_sample_graph = *size_sample_g;
	int i, j,ij, one = 1, dim = *p, pxp = dim * dim, p1 = dim - 1, p2x2 = ( dim - 2 ) * 2;
	
	double Dsij, sum_weights = 0.0, weight_C, sum_rates;
	bool this_one;

	string string_g;
	vector<string> sample_graphs_C( iteration - burn_in );

	//-- allocation for birth/death rates 
	vector<double> K121( 4 ); 
	vector<double> Kj12( p1 );              // K[j, -j]
	vector<double> sigmaj12( p1 );          // sigma[-j, j]  
	vector<double> sigmaj22( p1 * p1 );     // sigma[-j, -j]
	vector<double> Kj12xK22_inv( p1 ); 

	vector<double> K21( p2x2 );                          // K[-e, e]
	vector<double> sigma21( p2x2 );                      // sigma[-e, e]
	vector<double> sigma22( ( dim - 2 ) * ( dim - 2 ) ); // sigma[-e, -e]
	vector<double> sigma11_inv( 4 );                     // inv( sigma[e, e] )
	vector<double> sigma21xsigma11_inv( p2x2 ); 
	vector<double> K12xK22_inv( p2x2 );   
	// -- allocation for rgwish_sigma 
	vector<double> sigma_start( pxp ); 
	vector<double> inv_C( pxp ); 
	vector<double> beta_star( dim ); 
	vector<double> sigma_i( dim ); 
	vector<double> sigma_start_N_i( dim );   // For dynamic memory used
	vector<double> sigma_N_i( pxp );         // For dynamic memory used
	vector<int> N_i( dim );                  // For dynamic memory used
	// ----------------------------

	vector<double> sigma( pxp ); 
	vector<double> copyK( pxp ); 
	memcpy( &copyK[0], K, sizeof( double ) * pxp );
	inverse( &copyK[0], &sigma[0], &dim );			

	int qp = dim * ( dim - 1 ) / 2;
	vector<char> char_g( qp );              // char string_g[pp];

	// Counting size of notes
	int ip;
	vector<int> size_node( dim, 0 );
	for( i = 0; i < dim; i++ )
	{
		ip = i * dim;
		for( j = 0; j < dim; j++ ) size_node[i] += G[ip + j];
	}
	
	// For finding the index of rates 
	vector<double> rates( qp );
	vector<int> index_rates_row( qp );
	vector<int> index_rates_col( qp );
	int counter = 0;
	vector<double> Dsijj( pxp ); 
	for( j = 1; j < dim; j++ )
		for( i = 0; i < j; i++ )
		{
			index_rates_row[counter] = i;
			index_rates_col[counter] = j;
			counter++;
			
			// for calculating the birth/death rates
			ij        = j * dim + i;
			Dsij      = Ds[ij];
			Dsijj[ij] = Dsij * Dsij / Ds[j * dim + j]; 
		}

	int size_index = multi_update_C;
	vector<int> index_selected_edges( multi_update_C );

//-- Main loop for birth-death MCMC -------------------------------------------| 
	GetRNGstate();
	for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc += size_index )
	{
		if( ( i_mcmc + 1 ) % print_c < multi_update_C ) Rprintf( " Iteration  %d                 \n", i_mcmc + 1 ); 
		
//----- STEP 1: calculating birth and death rates -----------------------------|		

		rates_bdmcmc( &rates[0], G, Ds, &Dsijj[0],
                      &sigma[0], &sigma21[0], &sigma22[0], &sigmaj12[0], &sigmaj22[0],    
                      &K[0], &K21[0], &K121[0], &Kj12[0], 
                      &K12xK22_inv[0], &Kj12xK22_inv[0], &sigma11_inv[0], &sigma21xsigma11_inv[0],  
                      b, &dim );
				
		// Selecting multiple edges based on birth and death rates
		select_multi_edges( &rates[0], &index_selected_edges[0], &size_index, &sum_rates, &multi_update_C, &qp );

//----- Saving result ---------------------------------------------------------|	
		counter = 0;	
		for( j = 1; j < dim; j++ )
			for( i = 0; i < j; i++ )
				char_g[counter++] = G[ j * dim + i ] + '0'; 

		if( i_mcmc >= burn_in )
		{
			weight_C = 1.0 / sum_rates;
			
			//for( i = 0; i < pxp; i++ ) K_hat[i] += K[i] / sum_rates;
			F77_NAME(daxpy)( &pxp, &weight_C, &K[0], &one, &K_hat[0], &one );			

			string_g = string( char_g.begin(), char_g.end() );	
			all_weights[count_all_g] = weight_C;
			
			this_one = false;
			for( i = 0; i < size_sample_graph; i++ )
				if( sample_graphs_C[i] == string_g )
				{
					graph_weights[i] += all_weights[count_all_g];
					all_graphs[count_all_g] = i;
					this_one = true;
					break;
				} 
			
			if( !this_one || size_sample_graph == 0 )
			{
				sample_graphs_C[size_sample_graph] = string_g;
				graph_weights[size_sample_graph]   = all_weights[count_all_g];
				all_graphs[count_all_g]          = size_sample_graph; 
				size_sample_graph++;				
			}
			
			count_all_g++; 
			sum_weights += weight_C;
		} 
//----- End of saving result --------------------------------------------------|	
			
		// Updating graph based on selected edges
		for ( i = 0; i < size_index; i++ )
		{
			selected_edge_i = index_rates_row[ index_selected_edges[i] ];
			selected_edge_j = index_rates_col[ index_selected_edges[i] ];
			
			selected_edge_ij    = selected_edge_j * dim + selected_edge_i;
			G[selected_edge_ij] = 1 - G[selected_edge_ij];
			G[selected_edge_i * dim + selected_edge_j] = G[selected_edge_ij];
		
			if( G[selected_edge_ij] )
			{ 
				++size_node[selected_edge_i]; 
				++size_node[selected_edge_j]; 
			}
			else
			{ 
				--size_node[selected_edge_i]; 
				--size_node[selected_edge_j]; 
			}		
		}

//------ STEP 2: Sampling from G-Wishart for new graph ------------------------|
		rgwish_sigma( G, &size_node[0], Ts, K, &sigma[0], b_star, &dim, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		
	}  
	PutRNGstate();
//-- End of main loop for birth-death MCMC ------------------------------------| 

	for( i = 0; i < ( iteration - burn_in ); i++ ) 
	{
		sample_graphs_C[i].copy( sample_graphs[i], qp, 0 );
		sample_graphs[i][qp] = '\0';
	}
	
	*size_sample_g = size_sample_graph;
	*counter_all_g = count_all_g;
	
	for( i = 0; i < pxp; i++ ) K_hat[i] /= sum_weights;		
}
              
} // End of exturn "C"
