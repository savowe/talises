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

#include <ostream>
#include <fstream>
#include <string>
#include <cstring>
#include <array>

#include "CRT_Base.h"
#include "ParameterHandler.h"
#include "gsl/gsl_complex_math.h"
#include "gsl/gsl_eigen.h"
#include "gsl/gsl_blas.h"
#include "muParser.h"

using namespace std;

#ifndef __class_CRT_Base_IF__
#define __class_CRT_Base_IF__

/** Template class for interferometry in <B>dim</B> dimensions with <B>no_int_states</B> internal states
  *
  * In this template class functions for the interaction of a BEC with a light field are defined.
  */
template <class T, int dim, int no_int_states>
class CRT_Base_IF : public CRT_Base<T,dim,no_int_states>
{
public:
  CRT_Base_IF( ParameterHandler * );
  virtual ~CRT_Base_IF();

  void run_sequence();

protected:
  using CRT_Base<T,dim,no_int_states>::m_header;
  using CRT_Base<T,dim,no_int_states>::m_params;
  using CRT_Base<T,dim,no_int_states>::m_fields;
  using CRT_Base<T,dim,no_int_states>::m_custom_fct;
  using CRT_shared::m_no_of_pts;

  CPoint<dim> x;
  std::vector<double> psi_real;
  std::vector<double> psi_imag;

  double psi_real_array[50];
  double psi_imag_array[50];

  bool position_dependent;
  bool time_dependent;
  bool nonlinear;

  mu::Parser* V_parser;

  static void Do_NL_Step_Wrapper(void *,sequence_item &);
  static void Numerical_Diagonalization_Wrapper(void *,sequence_item &);

  void Do_NL_Step();
  void Numerical_Diagonalization();

  void UpdateParams();

  /// Define custom sequences
  virtual bool run_custom_sequence( const sequence_item & )=0;

};


/** Calls UpdateParams() and defines stepfunctions
  *
  * @param Pointer to ParameterHandler object to read from xml files
  */
template <class T, int dim, int no_int_states>
CRT_Base_IF<T,dim,no_int_states>::CRT_Base_IF( ParameterHandler *params ) : CRT_Base<T,dim,no_int_states>(params)
{
  // Map between "freeprop" and Do_NL_Step
  this->m_map_stepfcts["freeprop"] = &Do_NL_Step_Wrapper;
  this->m_map_stepfcts["interact"] = &Numerical_Diagonalization_Wrapper;

  UpdateParams();
}

/// Destructor
template <class T, int dim, int no_int_states>
CRT_Base_IF<T,dim,no_int_states>::~CRT_Base_IF()
{
}

/** Set values to interferometer variables from xml (m_params)
  *
  * Including:
  *   - the laser beams (amplitude, laser_k, ...)
  *   - gravitation
  *   - rabi_threshold
  */
template <class T, int dim, int no_int_states>
void CRT_Base_IF<T,dim,no_int_states>::UpdateParams()
{
  char s[100];
  //for ( int i=0; i<dim; i++)
  //  beta[i] = m_params->Get_VConstant("Beta",i);
}




/** Wrapper function for Do_NL_Step()
  * @param ptr Function pointer to be set to Do_NL_Step()
  * @param seq Additional information about the sequence (for example file names if a file has to be read)
  */
template <class T, int dim, int no_int_states>
void CRT_Base_IF<T,dim,no_int_states>::Do_NL_Step_Wrapper ( void *ptr, sequence_item &seq )
{
  CRT_Base_IF<T,dim,no_int_states> *self = static_cast<CRT_Base_IF<T,dim,no_int_states>*>(ptr);
  self->Do_NL_Step();
}


/** Wrapper function for Numerical_n()
  * @param ptr Function pointer to be set to Numerical_Diagonalization()
  * @param seq Additional information about the sequence (for example file names if a file has to be read)
  */
template <class T, int dim, int no_int_states>
void CRT_Base_IF<T,dim,no_int_states>::Numerical_Diagonalization_Wrapper ( void *ptr, sequence_item &seq )
{
  CRT_Base_IF<T,dim,no_int_states> *self = static_cast<CRT_Base_IF<T,dim,no_int_states>*>(ptr);
  self->Numerical_Diagonalization();
}

/** Solves the potential part without any external fields but
  * gravity.
  */
template <class T, int dim, int no_int_states>
void CRT_Base_IF<T,dim,no_int_states>::Do_NL_Step()
{
  const double dt = -m_header.dt;
  double re1, im1, tmp1, phi[no_int_states];
  CPoint<dim> x;
  int nNum = this->V_parser->GetNumResults();
  vector<fftw_complex *> Psi;
  //Vector for the components of the wavefunction
  for ( int i=0; i<no_int_states; i++ )
    Psi.push_back(m_fields[i]->Getp2In());


  for ( int l=0; l<this->m_no_of_pts; l++ )
  {
      this->x = this->m_fields[0]->Get_x(l);
      double *V_ptr = this->V_parser->Eval(nNum);

      for ( int i=0; i<no_int_states; i++ )
      {
		  double V_real = *(V_ptr+(2*i));
		  phi[i] = V_real*dt;
      }
    /*/Loop through column
    for ( int i=0; i<no_int_states; i++ )
    {
      double tmp_density = Psi[i][l][0]*Psi[i][l][0] + Psi[i][l][1]*Psi[i][l][1];
      if (tmp_density <= 0.0)
      {
    	 phi[i] = 0.0;
      } else
      {
    	  phi[i] = -this->m_b*log( tmp_density );
    	  phi[i] += this->m_gs[i+no_int_states*i]*tmp_density;
      }
      x = m_fields[0]->Get_x(l);
      //phi[i] += beta[0]*x[0]-DeltaL[i];
      phi[i] *= dt;
    }*/

    //Compute exponential: exp(V)*Psi
    for ( int i=0; i<no_int_states; i++ )
    {
      sincos( phi[i], &im1, &re1 );

      tmp1 = Psi[i][l][0];
      Psi[i][l][0] = Psi[i][l][0]*re1 - Psi[i][l][1]*im1;
      Psi[i][l][1] = Psi[i][l][1]*re1 + tmp1*im1;
    }
  }
}


/** Solves the potential part in the presence of light fields with a numerical method
  *
  * In this function \f$ \exp(V)\Psi \f$ is calculated. The matrix exponential is computed
  * with the help of a numerical diagonalisation which uses the gsl library
  */
template <class T, int dim, int no_int_states>
void CRT_Base_IF<T,dim,no_int_states>::Numerical_Diagonalization()
{
  int nNum = this->V_parser->GetNumResults();
  this->V_parser->Eval(nNum); // initializes nNum
  const long long int N_V_eval = this->m_no_of_pts*no_int_states*nNum ;
  double V_eval[N_V_eval];
  if (this->position_dependent == true) //Calculate V(r,t) at t for all r
  {
    for ( int l=0; l<this->m_no_of_pts; l++ ) //TODO parallelizing this would be good
    {
      for ( int i=0; i<no_int_states; i++ )
      {
        vector<fftw_complex *> Psi;
        Psi.push_back(m_fields[i]->Getp2In());
        double *psi_real_ptr = *(Psi[0]);
        double *psi_imag_ptr = *(Psi[0]+1);
        this->psi_real_array[i] = *psi_real_ptr;
        this->psi_imag_array[i] = *psi_imag_ptr;
      }
      this->x = this->m_fields[0]->Get_x(l);
      double *V_ptr = this->V_parser->Eval(nNum);
      for (int j=0; j<nNum; j++)
      {
        V_eval[l*nNum+j] = *(V_ptr+(j));
      }
    }
  }
  #pragma omp parallel
  {
	  double re1, im1;
    const double dt = -m_header.dt;

    vector<fftw_complex *> Psi;
    for ( int i=0; i<no_int_states; i++ )
      Psi.push_back(m_fields[i]->Getp2In());

    gsl_matrix_complex *A = gsl_matrix_complex_calloc(no_int_states,no_int_states);
    gsl_matrix_complex *B = gsl_matrix_complex_calloc(no_int_states,no_int_states);
    gsl_eigen_hermv_workspace *w = gsl_eigen_hermv_alloc(no_int_states);
    gsl_vector *eval = gsl_vector_alloc(no_int_states);
    gsl_vector_complex *Psi_1 = gsl_vector_complex_alloc(no_int_states);
    gsl_vector_complex *Psi_2 = gsl_vector_complex_alloc(no_int_states);
    gsl_matrix_complex *evec = gsl_matrix_complex_alloc(no_int_states,no_int_states);

    #pragma omp for
    for ( int l=0; l<this->m_no_of_pts; l++ )
    {
      gsl_matrix_complex_set_zero(A);
      gsl_matrix_complex_set_zero(B);

      int m = 0;
      for ( int i=0; i<no_int_states; i++ )
      {
          for ( int j=i; j<no_int_states; j++ )
          {
        	  double V_real = V_eval[l*nNum+2*m];
        	  double V_imag = V_eval[l*nNum+2*m+1];
            if (i != j) //nondiagonal elements
            {
				gsl_matrix_complex_set(A,i,j, {V_real,V_imag});
				gsl_matrix_complex_set(A,j,i, {V_real,-V_imag});
            }
            else
            { //diagonal elements
				gsl_matrix_complex_set(A,i,i, {V_real,0});
            }
            m += 1;
          }
      }

      /*
      //Diagonal elements + Nonlinear part: \Delta+g|\Phi|^2
      for ( int i=0; i<no_int_states; i++ )
      {
      	double tmp_density = Psi[i][l][0]*Psi[i][l][0] + Psi[i][l][1]*Psi[i][l][1];
      	if (tmp_density <= 0.0)
      	{
      		phi[i] = 0.0;
      	} else
      	{
      		phi[i] = -this->m_b*log( tmp_density );
      	}
        x = this->m_fields[0]->Get_x(l);
        phi[i] += DeltaL[i];
        gsl_matrix_complex_set(A,i,i, {phi[i],0});
      }
	*/


      //Compute Eigenvalues + Eigenvector
      gsl_eigen_hermv(A,eval,evec,w);

      // exp(Eigenvalues)
      for ( int i=0; i<no_int_states; i++ )
      {
        sincos( dt*gsl_vector_get(eval,i), &im1, &re1 );
        gsl_matrix_complex_set(B,i,i, {re1,im1});
      }

      // H_new = Eigenvector * exp(Eigenvalues) * conjugate(Eigenvector)
      gsl_blas_zgemm(CblasNoTrans,CblasConjTrans,GSL_COMPLEX_ONE,B,evec,GSL_COMPLEX_ZERO,A);
      gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,evec,A,GSL_COMPLEX_ZERO,B);

      for ( int i=0; i<no_int_states; i++)
      {
        gsl_vector_complex_set(Psi_1,i, {Psi[i][l][0],Psi[i][l][1]});
      }

      // H_new * Psi
      gsl_blas_zgemv(CblasNoTrans,GSL_COMPLEX_ONE,B,Psi_1,GSL_COMPLEX_ZERO,Psi_2);

      for ( int i=0; i<no_int_states; i++)
      {
        Psi[i][l][0] = gsl_vector_complex_get(Psi_2,i).dat[0];
        Psi[i][l][1] = gsl_vector_complex_get(Psi_2,i).dat[1];
      }
    }
    gsl_matrix_complex_free(A);
    gsl_matrix_complex_free(B);
    gsl_eigen_hermv_free(w);
    gsl_vector_free(eval);
    gsl_vector_complex_free(Psi_1);
    gsl_vector_complex_free(Psi_2);
    gsl_matrix_complex_free(evec);
  }
}

/** Run all the sequences defined in the xml file
  *
  * For furher information about the sequences see sequence_item
  */
template <class T, int dim, int no_int_states>
void CRT_Base_IF<T,dim,no_int_states>::run_sequence()
{
  if ( m_fields.size() != no_int_states )
  {
    std::cerr << "Critical Error: m_fields.size() != no_int_states\n";
    exit(EXIT_FAILURE);
  }

  StepFunction step_fct=nullptr;
  StepFunction half_step_fct=nullptr;
  StepFunction full_step_fct=nullptr;
  char filename[1024];

  std::cout << "FYI: Found " << m_params->m_sequence.size() << " sequences." << std::endl;


  try
  {
    half_step_fct = this->m_map_stepfcts.at("half_step");
    full_step_fct = this->m_map_stepfcts.at("full_step");
  }
  catch (const std::out_of_range &oor)
  {
    std::cerr << "Critical Error: Invalid fct ptr to ft_half_step or ft_full_step ()" << oor.what() << ")\n";
    exit(EXIT_FAILURE);
  }

  int seq_counter=1;

  for ( auto seq : m_params->m_sequence )
  {
    if ( run_custom_sequence(seq) )
    {
      seq_counter++;
      continue;
    }

  
    double max_duration = 0;
    for ( int i = 0; i < seq.duration.size(); i++)
      if (seq.duration[i] > max_duration )
        max_duration = seq.duration[i];

    int subN = int(max_duration / seq.dt);
    int Nk = seq.Nk;
    int Na = subN / seq.Nk;

    /* Definitions for the Hamiltonian parser */
    this->V_parser = new mu::Parser; // Parser in heap
    /** Read in Hamiltonian strings from XML */
    std::string V_expression = "";
    V_expression += seq.V_real[0];
    V_expression += ",";
    V_expression += seq.V_imag[0];
    for (int i = 1; i<seq.V_real.size(); i++)
    {
      V_expression += ",";
      V_expression += seq.V_real[i];
      V_expression += ",";
      V_expression += seq.V_imag[i];
    }


    /* Set final Hamiltonian*/
    this->V_parser->SetExpr(V_expression);
    // Get the map with the used variables
    const std::map<std::__cxx11::basic_string<char>, double*> variables =  this->V_parser->GetUsedVar();
    // Get the number of variables 
    std::map<std::__cxx11::basic_string<char>, double*>::const_iterator item = variables.begin();
    // Query the variables
    position_dependent = false;
    time_dependent = false;
    nonlinear = false;

    for (; item!=variables.end(); ++item)
    {
      if ((item->first == "x" ) or (item->first == "y" ) or (item->first == "z" ))
      {
        position_dependent = true;
        cout << "Name: " << item->first << " Address: [0x" << item->second << "]\n";
      }
      if (item->first == "t")
      {
        time_dependent = true;
        cout << "Name: " << item->first << " Address: [0x" << item->second << "]\n";
      }
      if (item->first.rfind("psi_", 0) == 0) 
      {
        nonlinear = true;
        cout << "Name: " << item->first << " Address: [0x" << item->second << "]\n";
      }
    }
    /* Define Variables and Constants*/
    // self-defined constants
    std::map<std::string, double>::iterator it = this->m_params->m_map_constants.begin();
    while(it != this->m_params->m_map_constants.end())
    {
      this->V_parser->DefineConst(it->first, (double)it->second);
      it++;
    }
    // constants
    this->V_parser->DefineConst("pi", (double)M_PI);
    this->V_parser->DefineConst("e", (double)M_E);
    // variables
    if (time_dependent == true) {this->V_parser->DefineVar("t", &this->Get_t());}
    if (position_dependent == true) 
    {
      this->V_parser->DefineVar("x", &this->x[0]);
      if (dim >=2) {this->V_parser->DefineVar("y", &this->x[1]);}
      if (dim == 3) {this->V_parser->DefineVar("z", &this->x[2]);}
    }

    if (nonlinear == true)
    {
      for (int i = 0; i < this->m_fields.size(); i++)
      {
        std::string tmp_str = "psi_real_";
        tmp_str += std::to_string(i+1);
        this->V_parser->DefineVar(tmp_str, &this->psi_real_array[i] );
        tmp_str = "psi_imag_";
        tmp_str += std::to_string(i+1);
        this->V_parser->DefineVar(tmp_str, &this->psi_imag_array[i] );
      }
    }




    std::cout << "FYI: started new sequence " << seq.name << "\n";
    std::cout << "FYI: sequence no : " << seq_counter << "\n";
    std::cout << "FYI: duration    : " << max_duration << "\n";
    std::cout << "FYI: dt          : " << seq.dt << "\n";
    //std::cout << "FYI: Na          : " << Na << "\n";
    //std::cout << "FYI: Nk          : " << Nk << "\n";
    //std::cout << "FYI: Na*Nk*dt    : " << double(Na*Nk)*seq.dt << "\n";


    if ( double(Na*Nk)*seq.dt != max_duration )
      std::cout << "FYI: double(Na*Nk)*seq.dt != max_duration\n";

    if ( this->Get_dt() != seq.dt )
      this->Set_dt(seq.dt);

    try
    {
      step_fct = this->m_map_stepfcts.at(seq.name);
    }
    catch (const std::out_of_range &oor)
    {
      std::cerr << "Critical Error: Invalid sequence name " << seq.name << "\n(" << oor.what() << ")\n";
      exit(EXIT_FAILURE);
    }

    if ( seq.name == "freeprop" )
    double backup_t = m_header.t;
    double backup_end_t = m_header.t;
    for ( int k=0; k<no_int_states; k++ ) // Delete old packed Sequence
    {
      sprintf( filename, "Seq_%d_%d.bin", seq_counter, k+1 );
      std::remove(filename);
    }

      for ( int i=1; i<=Na; i++ )
      {
        (*half_step_fct)(this,seq);
        for ( int j=2; j<=Nk; j++ )
        {
          (*step_fct)(this,seq);
          (*full_step_fct)(this,seq);
        }
        (*step_fct)(this,seq);
        (*half_step_fct)(this,seq);

        std::cout << "t = " << to_string(m_header.t) << std::endl;

        if ( seq.output_freq == freq::each )
        {
          for ( int k=0; k<no_int_states; k++ )
          {
            sprintf( filename, "%.3f_%d.bin", this->Get_t(), k+1 );
            this->Save_Phi( filename, k );
          }
        }

        if ( seq.output_freq == freq::packed )
        {
          for ( int k=0; k<no_int_states; k++ )
          {
            sprintf( filename, "Seq_%d_%d.bin", seq_counter, k+1 );
            this->Append_Phi( filename, k );
          }
        }

        if ( seq.compute_pn_freq == freq::each )
        {
          for ( int c=0; c<no_int_states; c++ )
            std::cout << "N[" << c << "] = " << this->Get_Particle_Number(c) << std::endl;
        }

        if ( seq.custom_freq == freq::each && m_custom_fct != nullptr )
        {
          (*m_custom_fct)(this,seq);
        }
      }

      if (seq.output_freq == freq::last )
      {
        for ( int k=0; k<no_int_states; k++ )
        {
          sprintf( filename, "%.3f_%d.bin", this->Get_t(), k+1 );
          this->Save_Phi( filename, k );
        }
      }

      if ( seq.compute_pn_freq == freq::last )
      {
        for ( int c=0; c<no_int_states; c++ )
          std::cout << "N[" << c << "] = " << this->Get_Particle_Number(c) << std::endl;
      }

      if ( seq.custom_freq == freq::last && m_custom_fct != nullptr )
      {
        (*m_custom_fct)(this,seq);
      }

    seq_counter++;
  } // end of sequence loop
}
#endif
