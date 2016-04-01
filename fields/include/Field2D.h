/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert & Johannes Rauch

   This file is part of GENFIT.

   GENFIT is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   GENFIT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with GENFIT.  If not, see <http://www.gnu.org/licenses/>.
*/
/** @addtogroup genfit
 * @{
 */

#ifndef genfit_Field2D_h
#define genfit_Field2D_h

#include "AbsBField.h"

#include "boost/tuple/tuple.hpp"
#include "boost/tuple/tuple_comparison.hpp"
#include <map>


namespace genfit {

/** @brief 2D Field
 *
 *  @author Haiwang Yu (New Mexico State University)
 * 
 */
class Field2D : public AbsBField {
 public:
  //! define the constant field in this ctor
  Field2D()
  { field_map_.clear(); }

  Field2D(std::string inname)
  { initialize(inname); }

  bool initialize(std::string inname = "");

  void plot(std::string option = "");

  //! return value at position
  TVector3 get(const TVector3& pos) const;
  void get(const double& posX, const double& posY, const double& posZ, double& Bx, double& By, double& Bz) const;

 private:
  typedef std::map< boost::tuple<double,double>, boost::tuple<double,double> > BMAP2D;
  typedef std::pair< boost::tuple<double,double>, boost::tuple<double,double> > BPAIR2D;
  BMAP2D field_map_;
};

} /* End of namespace genfit */
/** @} */

#endif // genfit_Field2D_h
