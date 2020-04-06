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

#ifndef CFT_1D_H
#define CFT_1D_H

#include "cft_base.h"


/// Contains classes for a Fourier transform in one, two or three dimensions with complex and real valued data
namespace Fourier
{
  /// Class for Fourier transform in one dimension with complex valued data
  class cft_1d : public Fourier::cft_base<1>
  {
  public:
    cft_1d( const generic_header&, bool=true, bool=false );

    void ft( int isign ); // -1 (forward) oder +1 (backward)
    void D1();
    void D2();

    CPoint<1> Get_k(const int64_t) final;
    CPoint<1> Get_x(const int64_t) final;
  protected:

    void fix( fftw_complex* data, double d );
    void scale( fftw_complex* data, double sx );
  };
}
#endif
