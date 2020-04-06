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
#include "pugixml.hpp"
#include "muParser.h"
#include "my_structs.h"
#include <string>
#include <map>
#include <vector>

#ifndef __class_ParameterHandler__
#define __class_ParameterHandler__

/** \file Parameterhandler.h */

enum freq { none=0, each=1, last=2, packed=3 };

/** Contains elements for controlling a sequence */
struct sequence_item
{
  std::string name; ///< sequence identifier
  std::string content;

  std::vector<std::string> V_real;
  std::vector<std::string> V_imag;
  //std::string H;
  std::vector<double> duration; ///< duration of the sequence
  double dt; ///< dt for this sequence

  double chirp_min; ///< lower bound of interval for a phasescan
  double chirp_max; ///< upper bound of interval for a phasescan
  int no_of_chirps; ///< number of chirps for a phasescan
  int chirp_mode; ///< chirp mode
  int rabi_output_freq; ///< output frequency for rabi oscillations
  int output_freq; ///< output frequency of data files
  int custom_freq; ///< set output frequency of custom data files (e.g. slices)
  int comp; ///< component index of the wave function
  int compute_pn_freq; ///< set frequency for computing particle numbers
  int analyze; ///< output frequency for analyzing tools
  int Nk; ///< number of intermediate steps
  double time;
};

struct analyze_item
{
  std::string name;
  bool content;
  bool separate;
};

/// Class for reading from xml files to control the simulation
class ParameterHandler
{
public:
  /// Constructor which needs the name of the xml file as an input parameter
  ParameterHandler( const std::string );
  ~ParameterHandler();

  void Get_Header( generic_header&, bool=true );

  std::string Get_simulation( const std::string );

  /** Returns value of a tag <string> in the Constant section of the xml file */
  double Get_Constant( const std::string );
  /** Returns value of element int of a tag <string> in the VConstant section of the xml file */
  double Get_VConstant( const std::string, const int );
  /** Returns values as a vector of a tag <string> in the VConstant section of the xml file */
  std::vector<double> Get_VConstant( const std::string );
  /** Returns all momentum states in a vector */
  bool Get_MomentumStates( std::vector<double>&, const int );

  double Get_t();
  double Get_t_scale();
  double Get_dt();
  double Get_epsilon();
  double Get_stepsize();
  double Get_xMin();
  double Get_xMax();
  double Get_yMin();
  double Get_yMax();
  double Get_zMin();
  double Get_zMax();
  double Get_L();
  double Get_T();
  double Get_M();

  int Get_NX();
  int Get_NY();
  int Get_NZ();
  int Get_NA();
  int Get_NK();
  int Get_MaxIter();

  void Setup_muParser( mu::Parser& );
  std::map<std::string,double> m_map_constants; ///< xml -> double (for constant scalar values)
  std::vector<sequence_item> m_sequence; ///< vector of sequence_items ( Elements of the sequence )
  std::vector<analyze_item> m_analyze; ///< vector of analyze_items (for ana_tools)
protected:
  void populate_constants(); ///< Read constant values from xml and populate m_map_constants
  void populate_vconstants(); ///< Read vconstants values from xml and populate m_map_vconstants
  void populate_algorithm(); ///< Read functions and save in m_map_algorithm
  void populate_simulation(); ///< populate m_map_simulation
  void populate_sequence(); ///< Read from sequence node and save in m_sequence
  void populate_analyze();

  pugi::xml_document m_xml_doc; ///< load document here
  std::map<std::string,int> m_map_ai_type;
  std::map<std::string,int> m_map_freq; ///< xml -> int (options none, each and last e.g. for the frequency of computing particle numbers)
  std::map<std::string,std::vector<double>> m_map_vconstants; ///< xml -> double (for constant vectors)
  std::map<std::string,std::string> m_map_algorithm; ///< xml -> string (function)
  std::map<std::string,std::string> m_map_simulation;
};

#endif
