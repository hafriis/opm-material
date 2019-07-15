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
 * \copydoc Opm::RegularizedBrooksCoreyParamsVE
 */
#ifndef OPM_REGULARIZED_BROOKS_COREY_PARAMS_VE_HPP
#define OPM_REGULARIZED_BROOKS_COREY_PARAMS_VE_HPP

#include "RegularizedBrooksCorey.hpp"
#include "RegularizedBrooksCoreyParams.hpp"

#include <cassert>

#include <opm/material/common/EnsureFinalized.hpp>

namespace Opm {
/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief   Parameters that are necessary for the \em regularization of
 *          the Brooks-Corey capillary pressure model.
 */
template <class TraitsT>
class RegularizedBrooksCoreyParamsVE : public Opm::RegularizedBrooksCoreyParams<TraitsT>
{
    typedef Opm::RegularizedBrooksCoreyParams<TraitsT> RegularizedBrooksCoreyParams;
    typedef typename TraitsT::Scalar Scalar;

public:
    typedef TraitsT Traits;

    RegularizedBrooksCoreyParamsVE()
        : RegularizedBrooksCoreyParams()
        , srw_(0.0)
        , srn_(0.0)
        , krnEndPoint_(0.01)
        , krwEndPoint_(0.01)
    {
    }

    RegularizedBrooksCoreyParamsVE(Scalar entryPressure, Scalar lambda)
        : RegularizedBrooksCoreyParams(entryPressure, lambda)
        , srw_(0.0)
        , srn_(0.0)
        , krnEndPoint_(0.01)
        , krwEndPoint_(0.01)
    { finalize(); }

    /*!
     * \brief Calculate all dependent quantities once the independent
     *        quantities of the parameter object have been set.
     */
    void finalizePlain()
    {
        RegularizedBrooksCoreyParams::finalize();
    }

    void finalize()
    {
    }

    /*!
     * \brief Return the threshold saturation below which the
     *        capillary pressure is regularized.
     */
    Scalar krnEndPoint() const
    { EnsureFinalized::check(); return krnEndPoint_; }

    /*!
     * \brief Return the capillary pressure at the low threshold
     *        saturation of the wetting phase.
     */
    Scalar krwEndPoint() const
    { EnsureFinalized::check(); return krwEndPoint_; }

    /*!
     * \brief Set the threshold saturation below which the capillary
     *        pressure is regularized.
     */
    void setKrnEndPoint(Scalar value)
    { krnEndPoint_ = value; }

    /*!
     * \brief Set the threshold saturation below which the capillary
     *        pressure is regularized.
     */
    void setKrwEndPoint(Scalar value)
    { krwEndPoint_ = value; }

    //*********************HAF********************************************
    void setResidualSaturation(const int phaseIdx, Scalar value)
    {
        if (phaseIdx == Traits::wettingPhaseIdx)
        {
            srw_ = value;
        }
        else if (phaseIdx == Traits::nonWettingPhaseIdx)
        {
            srn_ = value;
        }
    }

    Scalar getResidualSaturation(const int phaseIdx)
    {
        Scalar ret = 0.0;
        if (phaseIdx == Traits::wettingPhaseIdx)
        {
            ret = srw_;
        }
        else if (phaseIdx == Traits::nonWettingPhaseIdx)
        {
            ret = srn_;
        }

        return ret;
    }

    Scalar compute_h_FromSandSmax(const Scalar& S, const Scalar& Smax, const Scalar& H) const
    {
        return (H*(S*(1.0-srw_)-Smax*srn_))/((1.0-srw_)*(1.0-srw_-srn_));
    }

    Scalar compute_hmax_FromSandSmax(const Scalar& S, const Scalar& Smax, const Scalar& H) const
    {
        return (H*Smax)/(1.0-srw_);
    }

    Scalar computeNonWettingPhaseRelPerm(const Scalar& h, const Scalar& hmax, const Scalar& H) const
    {
        return krnEndPoint_*(h/H);
    }

    Scalar computeWettingPhaseRelPerm(const Scalar& h, const Scalar& hmax, const Scalar& H, const Scalar& viscosity_w) const
    {
        return (((H - hmax)/H) + (viscosity_w*krwEndPoint_)*((hmax-h)/H));
    }

private:

    Scalar srw_;
    Scalar srn_;
    Scalar krnEndPoint_;
    Scalar krwEndPoint_;
    //Legge in permeabilitet!!!!!!!!!!!!!!!! Ting bør ikke endre seg i tiden!!!
    //hva med tetthet???
    //hva med viskositet pga. computeWettingPhaseRelPerm???
};
} // namespace Opm

#endif
