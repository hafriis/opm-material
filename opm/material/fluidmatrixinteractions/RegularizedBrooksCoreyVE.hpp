// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Opm::RegularizedBrooksCoreyVE
 */
#ifndef REGULARIZED_BROOKS_COREY_VE_HPP
#define REGULARIZED_BROOKS_COREY_VE_HPP

#include "RegularizedBrooksCorey.hpp"
#include "RegularizedBrooksCoreyParamsVE.hpp"

#include <opm/material/common/Spline.hpp>

namespace Opm {
/*!
 * \ingroup FluidMatrixInteractions
 * \brief Implementation of the regularized Brooks-Corey capillary
 *        pressure / relative permeability <-> saturation relation.
 *
 * This class bundles the "raw" curves as static members and doesn't
 * concern itself converting absolute to effective saturations and
 * vice versa.
 *
 * In order to avoid very steep gradients the marginal values are
 * "regularized".  This means that in stead of following the curve of
 * the material law in these regions, some linear approximation is
 * used.  Doing this is not worse than following the material
 * law. E.g. for very low wetting phase values the material laws
 * predict infinite values for \f$p_c\f$ which is completely
 * unphysical. In case of very high wetting phase saturations the
 * difference between regularized and "pure" material law is not big.
 *
 * Regularizing has the additional benefit of being numerically
 * friendly: Newton's method does not like infinite gradients.
 *
 * The implementation is accomplished as follows:
 * - check whether we are in the range of regularization
 *   - yes: use the regularization
 *   - no: forward to the standard material law.
 *
 * \see BrooksCorey
 */
template <class TraitsT, class ParamsT = RegularizedBrooksCoreyParamsVE<TraitsT> >
class RegularizedBrooksCoreyVE : public TraitsT
{


public:
    typedef Opm::RegularizedBrooksCorey<TraitsT, ParamsT> RegularizedBrooksCorey;

    typedef TraitsT Traits;
    typedef ParamsT Params;
    typedef typename Traits::Scalar Scalar;

    //! The number of fluid phases
    static const int numPhases = Traits::numPhases;
    static_assert(numPhases == 2,
                  "The regularized Brooks-Corey capillary pressure law only "
                  "applies to the case of two fluid phases");

    //! Specify whether this material law implements the two-phase
    //! convenience API
    static const bool implementsTwoPhaseApi = true;

    //! Specify whether this material law implements the two-phase
    //! convenience API which only depends on the phase saturations
    static const bool implementsTwoPhaseSatApi = true;

    //! Specify whether the quantities defined by this material law
    //! are saturation dependent
    static const bool isSaturationDependent = true;

    //! Specify whether the quantities defined by this material law
    //! are dependent on the absolute pressure
    static const bool isPressureDependent = false;

    //! Specify whether the quantities defined by this material law
    //! are temperature dependent
    static const bool isTemperatureDependent = false;

    //! Specify whether the quantities defined by this material law
    //! are dependent on the phase composition
    static const bool isCompositionDependent = false;

    static_assert(Traits::numPhases == 2,
                  "The number of fluid phases must be two if you want to use "
                  "this material law!");

    /*!
     * \brief The capillary pressure-saturation curves depending on absolute saturations.
     *
     * \param values A random access container which stores the
     *               relative pressure of each fluid phase.
     * \param params The parameter object expressing the coefficients
     *               required by the material law.
     * \param fs The fluid state for which the capillary pressure
     *           ought to be calculated
     */
    template <class Container, class FluidState>
    static void capillaryPressures(Container& values, const Params& params, const FluidState& fs)
    {
#warning TODO --- NOTE: We first treat the case with zero fine-scale capillary pressure.
        RegularizedBrooksCorey::capillaryPressures(values, params, fs);
        Scalar density_n = fs.density(Traits::nonWettingPhaseIdx);//OK???
        Scalar density_w = fs.density(Traits::wettingPhaseIdx);//OK???
        Scalar g = 9.80665; //gravity acceleration
        const auto& S = fs.saturation(Traits::nonWettingPhaseIdx);
        Scalar Smax = fs.getSmax();
        Scalar h = params.compute_h_FromSandSmax(S, Smax);

        values[Traits::wettingPhaseIdx] = 0.0; // reference phase
        values[Traits::nonWettingPhaseIdx] = (density_w - density_n)*g*h;
    }

    /*!
     * \brief Calculate the saturations of the phases starting from
     *        their pressure differences.
     */
    template <class Container, class FluidState>
    static void saturations(Container& values, const Params& params, const FluidState& fs)
    {
#warning TODO: not really required! --- What to do here in this case??? --- Should actually anything be done here???
        RegularizedBrooksCorey::saturations(values, params, fs);
    }

    /*!
     * \brief The relative permeability-saturation curves depending on absolute saturations.
     *
     * \param values A random access container which stores the
     *               relative permeability of each fluid phase.
     * \param params The parameter object expressing the coefficients
     *               required by the material law.
     * \param fs The fluid state for which the relative permeabilities
     *           ought to be calculated
     */
    template <class Container, class FluidState>
    static void relativePermeabilities(Container& values, const Params& params, const FluidState& fs)
    {
        //RegularizedBrooksCorey::relativePermeabilities(values, params, fs);
        //values[Traits::wettingPhaseIdx] = krw<FluidState, Evaluation>(params, fs);
        //values[Traits::nonWettingPhaseIdx] = krn<FluidState, Evaluation>(params, fs);
        
#warning TODO --- Seems OK now...
         const auto& S = fs.saturation(Traits::nonWettingPhaseIdx);
         Scalar Smax = fs.getSmax();
         Scalar h = params.compute_h_FromSandSmax(S, Smax);
         Scalar hmax = params.compute_hmax_FromSandSmax(S, Smax);
         Scalar viscosity_w = fs.viscosity(Traits::wettingPhaseIdx);//OK???

         values[Traits::wettingPhaseIdx] = params.computeWettingPhaseRelPerm(h, hmax, viscosity_w);
         values[Traits::nonWettingPhaseIdx] = params.computeNonWettingPhaseRelPerm(h, hmax);
    }
    
};
} // namespace Opm

#endif
