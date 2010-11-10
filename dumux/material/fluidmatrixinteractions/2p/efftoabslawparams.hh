// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief A default implementation of the parameters for the adapter
 *        class to convert material laws from effective to absolute
 *        saturations.
 */
#ifndef DUMUX_EFF_TO_ABS_LAW_PARAMS_HH
#define DUMUX_EFF_TO_ABS_LAW_PARAMS_HH

namespace Dumux
{
/*!
 * \ingroup material
 *
 * \brief A default implementation of the parameters for the adapter
 *        class to convert material laws from effective to absolute
 *        saturations.
 */
template <class EffLawParamsT>
class EffToAbsLawParams : public EffLawParamsT
{
    typedef EffLawParamsT EffLawParams;
public:
    typedef typename EffLawParams::Scalar Scalar;

    EffToAbsLawParams()
        : EffLawParams()
    { Swr_ = Snr_ = 0; }

    /*!
     * \brief Return the residual wetting saturation.
     */
    Scalar Swr() const
    { return Swr_; }

    /*!
     * \brief Set the residual wetting saturation.
     */
    void setSwr(Scalar v)
    { Swr_ = v; }

    /*!
     * \brief Return the residual non-wetting saturation.
     */
    Scalar Snr() const
    { return Snr_; }

    /*!
     * \brief Set the residual non-wetting saturation.
     */
    void setSnr(Scalar v)
    { Snr_ = v; }

private:
    Scalar Swr_;
    Scalar Snr_;
};

}

#endif