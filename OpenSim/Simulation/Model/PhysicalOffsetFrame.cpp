/* -------------------------------------------------------------------------- *
 *                  OpenSim:  PhysicalOffsetFrame.cpp                         *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2017 Stanford University and the Authors                *
 * Author(s): Ajay Seth                                                       *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

//=============================================================================
// INCLUDES
//=============================================================================
#include "PhysicalOffsetFrame.h"

using namespace OpenSim;

void PhysicalOffsetFrame::extendAddToSystem(SimTK::MultibodySystem& system) const
{
    Super::extendAddToSystem(system);

    // note: this frame might be attached to other `PhysicalOffsetFrame`s, or dependent
    //       frames, so use `findBaseFrame` to ensure that we "dig" all the way to the
    //       underlying mobilized body (e.g. an `OpenSim::Body`)
    setMobilizedBodyIndex(dynamic_cast<const PhysicalFrame&>(findBaseFrame()).getMobilizedBodyIndex());
}

void PhysicalOffsetFrame::extendSetMobilizedBodyIndex(const SimTK::MobilizedBodyIndex& mbix) const
{
    // the `MobilizedBodyIndex` has been set on this frame (by `setMobilizedBodyIndex`), but
    // this extension point also ensures that the `MobilizedBodyIndex` is also transitively
    // assigned to the parent frame (which may, itself, be a `PhysicalOffsetFrame`, or a
    // `Body`, etc.)
    PhysicalFrame::setMobilizedBodyIndexOf(getParentFrame(), mbix);
}
