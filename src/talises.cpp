// This file is part of TALISES.
//
// TALISES is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// TALISES is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with TALISES.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright (C) 2020 Sascha Vowe
// Copyright (C) 2017 Želimir Marojević, Ertan Göklü, Claus Lämmerzahl - Original implementation in ATUS2

#include <iomanip>
#include <cmath>
#include <omp.h>
#include "cft_1d.h"
#include "cft_2d.h"
#include "cft_3d.h"
#include "muParser.h"
#include "ParameterHandler.h"
#include "CRT_Base_IF.h"

using namespace std;

namespace RT_Solver
{
  template<class T, int dim, int internal_dim>
  class Raman_single : public CRT_Base_IF<T,dim,internal_dim> //TODO number of internal states sollte erkannt werden.
  {
  public:
    Raman_single( ParameterHandler * );
    virtual ~Raman_single() {};

  protected:
    bool run_custom_sequence( const sequence_item & );
  };

  template<class T, int dim, int internal_dim>
  Raman_single<T,dim,internal_dim>::Raman_single( ParameterHandler *p ) : CRT_Base_IF<T,dim,internal_dim>( p )
  {
  }

  template<class T, int dim, int internal_dim>
  bool Raman_single<T,dim,internal_dim>::run_custom_sequence( const sequence_item &item )
  {
    // return true if a custom sequence is found or else
    return false;
  }
}

int main( int argc, char *argv[] ){
  if ( argc != 2 )
  {
    printf( "No parameter xml file specified.\n" );
    return EXIT_FAILURE;
  }

  ParameterHandler params(argv[1]);
  int dim=0;
  int internal_dim = 0;
  int no_of_threads = 1;

  try
  {
    std::string tmp = params.Get_simulation("DIM");
    dim = std::stod(tmp);
    tmp = params.Get_simulation("INTERNAL_DIM");
    internal_dim = std::stod(tmp);
  }
  catch (mu::Parser::exception_type &e)
  {
    cout << "Message:  " << e.GetMsg() << "\n";
    cout << "Formula:  " << e.GetExpr() << "\n";
    cout << "Token:    " << e.GetToken() << "\n";
    cout << "Position: " << e.GetPos() << "\n";
    cout << "Errc:     " << e.GetCode() << "\n";
  }

  no_of_threads = std::stod(params.Get_simulation("N_THREADS"));

  char *envstr = getenv( "MY_NO_OF_THREADS" );
  if ( envstr != nullptr ) no_of_threads = atoi( envstr );

  fftw_init_threads();
  fftw_plan_with_nthreads( no_of_threads );
  omp_set_num_threads( no_of_threads );

  std::cout << "FYI: Number of threads : " << no_of_threads << "\n";

  try //TODO hardcode more options for internal levels lol
  {
    if ( dim == 1 )
    {
		if (internal_dim == 1)
    	{
		  RT_Solver::Raman_single<Fourier::cft_1d,1,1> rtsol( &params );
		  rtsol.run_sequence();
    	}
    	if (internal_dim == 2)
    	{
		  RT_Solver::Raman_single<Fourier::cft_1d,1,2> rtsol( &params );
		  rtsol.run_sequence();
    	}
    	else if (internal_dim == 3)
    	{
		  RT_Solver::Raman_single<Fourier::cft_1d,1,3> rtsol( &params );
		  rtsol.run_sequence();
    	}
    	else if (internal_dim == 4)
    	{
		  RT_Solver::Raman_single<Fourier::cft_1d,1,4> rtsol( &params );
		  rtsol.run_sequence();
    	}
    	else if (internal_dim == 5)
    	{
		  RT_Solver::Raman_single<Fourier::cft_1d,1,5> rtsol( &params );
		  rtsol.run_sequence();
    	}
    	else if (internal_dim == 6)
    	{
		  RT_Solver::Raman_single<Fourier::cft_1d,1,6> rtsol( &params );
		  rtsol.run_sequence();
    	}
    	else if (internal_dim == 7)
    	{
		  RT_Solver::Raman_single<Fourier::cft_1d,1,7> rtsol( &params );
		  rtsol.run_sequence();
    	}
    	else if (internal_dim == 8)
    	{
		  RT_Solver::Raman_single<Fourier::cft_1d,1,8> rtsol( &params );
		  rtsol.run_sequence();
    	}
    }
    else if ( dim == 2 )
    {
		if (internal_dim == 1)
    	{
		  RT_Solver::Raman_single<Fourier::cft_2d,2,1> rtsol( &params );
		  rtsol.run_sequence();
    	}
    	if (internal_dim == 2)
    	{
		  RT_Solver::Raman_single<Fourier::cft_2d,2,2> rtsol( &params );
		  rtsol.run_sequence();
    	}
    	else if (internal_dim == 3)
    	{
		  RT_Solver::Raman_single<Fourier::cft_2d,2,3> rtsol( &params );
		  rtsol.run_sequence();
    	}
    	else if (internal_dim == 4)
    	{
		  RT_Solver::Raman_single<Fourier::cft_2d,2,4> rtsol( &params );
		  rtsol.run_sequence();
    	}
    	else if (internal_dim == 5)
    	{
		  RT_Solver::Raman_single<Fourier::cft_2d,2,5> rtsol( &params );
		  rtsol.run_sequence();
    	}
    	else if (internal_dim == 6)
    	{
		  RT_Solver::Raman_single<Fourier::cft_2d,2,6> rtsol( &params );
		  rtsol.run_sequence();
    	}
    	else if (internal_dim == 7)
    	{
		  RT_Solver::Raman_single<Fourier::cft_2d,2,7> rtsol( &params );
		  rtsol.run_sequence();
    	}
    	else if (internal_dim == 8)
    	{
		  RT_Solver::Raman_single<Fourier::cft_2d,2,8> rtsol( &params );
		  rtsol.run_sequence();
    	}
    }
    else if ( dim == 3 )
    {
		if (internal_dim == 1)
    	{
		  RT_Solver::Raman_single<Fourier::cft_3d,3,1> rtsol( &params );
		  rtsol.run_sequence();
    	}
    	if (internal_dim == 2)
    	{
		  RT_Solver::Raman_single<Fourier::cft_3d,3,2> rtsol( &params );
		  rtsol.run_sequence();
    	}
    	else if (internal_dim == 3)
    	{
		  RT_Solver::Raman_single<Fourier::cft_3d,3,3> rtsol( &params );
		  rtsol.run_sequence();
    	}
    	else if (internal_dim == 4)
    	{
		  RT_Solver::Raman_single<Fourier::cft_3d,3,4> rtsol( &params );
		  rtsol.run_sequence();
    	}
    	else if (internal_dim == 5)
    	{
		  RT_Solver::Raman_single<Fourier::cft_3d,3,5> rtsol( &params );
		  rtsol.run_sequence();
    	}
    	else if (internal_dim == 6)
    	{
		  RT_Solver::Raman_single<Fourier::cft_3d,3,6> rtsol( &params );
		  rtsol.run_sequence();
    	}
    	else if (internal_dim == 7)
    	{
		  RT_Solver::Raman_single<Fourier::cft_3d,3,7> rtsol( &params );
		  rtsol.run_sequence();
    	}
    	else if (internal_dim == 8)
    	{
		  RT_Solver::Raman_single<Fourier::cft_3d,3,8> rtsol( &params );
		  rtsol.run_sequence();
    	}
    }
    else
    {
      cout << "You have found a new dimension!" << endl;
    }
  }
  catch (mu::Parser::exception_type &e)
  {
    cout << "Message:  " << e.GetMsg() << "\n";
    cout << "Formula:  " << e.GetExpr() << "\n";
    cout << "Token:    " << e.GetToken() << "\n";
    cout << "Position: " << e.GetPos() << "\n";
    cout << "Errc:     " << e.GetCode() << "\n";
  }
  catch (std::string &str)
  {
    cout << str << endl;
  }

  fftw_cleanup_threads();
  return EXIT_SUCCESS;
}
