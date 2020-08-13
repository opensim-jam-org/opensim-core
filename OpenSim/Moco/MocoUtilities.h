#ifndef OPENSIM_MOCOUTILITIES_H
#define OPENSIM_MOCOUTILITIES_H
/* -------------------------------------------------------------------------- *
 * OpenSim: MocoUtilities.h                                                   *
 * -------------------------------------------------------------------------- *
 * Copyright (c) 2017 Stanford University and the Authors                     *
 *                                                                            *
 * Author(s): Christopher Dembia, Nicholas Bianco                             *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0          *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

#include "MocoTrajectory.h"
#include "osimMocoDLL.h"
#include <regex>

#include <OpenSim/Common/GCVSplineSet.h>
#include <OpenSim/Common/Logger.h>
#include <OpenSim/Common/PiecewiseLinearFunction.h>
#include <OpenSim/Common/Reporter.h>
#include <OpenSim/Common/Storage.h>
#include <OpenSim/Common/CommonUtilities.h>
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/StatesTrajectory.h>

namespace OpenSim {

class StatesTrajectory;
class Model;
class MocoTrajectory;
class MocoProblem;

/// Calculate the requested outputs using the model in the problem and the
/// states and controls in the MocoTrajectory.
/// The output paths can be regular expressions. For example,
/// ".*activation" gives the activation of all muscles.
///
/// The output paths must correspond to outputs that match the type provided in
/// the template argument, otherwise they are not included in the report.
///
/// @note The provided trajectory is not modified to satisfy kinematic
/// constraints, but SimTK::Motions in the Model (e.g., PositionMotion) are
/// applied. Therefore, this function expects that you've provided a trajectory
/// that already satisfies kinematic constraints. If your provided trajectory
/// does not satisfy kinematic constraints, many outputs will be incorrect.
/// For example, in a model with a patella whose location is determined by a
/// CoordinateCouplerConstraint, the length of a muscle that crosses the patella
/// will be incorrect.
///
/// @note Parameters and Lagrange multipliers in the MocoTrajectory are **not**
///       applied to the model.
/// @ingroup mocoutil
template <typename T>
TimeSeriesTable_<T> analyze(Model model, const MocoTrajectory& trajectory,
        std::vector<std::string> outputPaths) {

    // Initialize the system so we can access the outputs.
    model.initSystem();
    // Create the reporter object to which we'll add the output data to create
    // the report.
    auto* reporter = new TableReporter_<T>();
    // Loop through all the outputs for all components in the model, and if
    // the output path matches one provided in the argument and the output type
    // agrees with the template argument type, add it to the report.
    for (const auto& comp : model.getComponentList()) {
        for (const auto& outputName : comp.getOutputNames()) {
            const auto& output = comp.getOutput(outputName);
            auto thisOutputPath = output.getPathName();
            for (const auto& outputPathArg : outputPaths) {
                if (std::regex_match(
                            thisOutputPath, std::regex(outputPathArg))) {
                    // Make sure the output type agrees with the template.
                    if (dynamic_cast<const Output<T>*>(&output)) {
                        log_debug("Adding output {} of type {}.",
                                output.getPathName(), output.getTypeName());
                        reporter->addToReport(output);
                    } else {
                        log_warn("Ignoring output {} of type {}.",
                                output.getPathName(), output.getTypeName());
                    }
                }
            }
        }
    }
    model.addComponent(reporter);
    model.initSystem();

    // Get states trajectory.
    TimeSeriesTable states = trajectory.exportToStatesTable();
    auto statesTraj = StatesTrajectory::createFromStatesTable(model, states);

    // Loop through the states trajectory to create the report.
    for (int i = 0; i < (int)statesTraj.getSize(); ++i) {
        // Get the current state.
        auto state = statesTraj[i];

        // Enforce any SimTK::Motion's included in the model.
        model.getSystem().prescribe(state);

        // Create a SimTK::Vector of the control values for the current state.
        SimTK::RowVector controlsRow =
                trajectory.getControlsTrajectory().row(i);
        SimTK::Vector controls(controlsRow.size(),
                controlsRow.getContiguousScalarData(), true);

        // Set the controls on the state object.
        model.realizeVelocity(state);
        model.setControls(state, controls);

        // Generate report results for the current state.
        model.realizeReport(state);
    }

    return reporter->getTable();
}

/// Given a MocoTrajectory and the associated OpenSim model, return the model
/// with a prescribed controller appended that will compute the control values
/// from the MocoTrajectory. This can be useful when computing state-dependent
/// model quantities that require realization to the Dynamics stage or later.
/// The function used to fit the controls can either be GCVSpline or
/// PiecewiseLinearFunction.
/// @ingroup mocoutil
OSIMMOCO_API void prescribeControlsToModel(const MocoTrajectory& trajectory,
        Model& model, std::string functionType = "GCVSpline");

/// Use the controls and initial state in the provided trajectory to simulate
/// the model using an ODE time stepping integrator (OpenSim::Manager), and
/// return the resulting states and controls. We return a MocoTrajectory (rather
/// than a StatesTrajectory) to facilitate comparing optimal control solutions
/// with time stepping. Use integratorAccuracy to override the default setting.
///
/// @note This function expects all Actuator%s in the model to be in the Model's
/// ForceSet.
/// @ingroup mocoutil
OSIMMOCO_API MocoTrajectory simulateTrajectoryWithTimeStepping(
        const MocoTrajectory& trajectory, Model model,
        double integratorAccuracy = SimTK::NaN);

/// The map provides the index of each state variable in
/// SimTK::State::getY() from its each state variable path string.
/// Empty slots in Y (e.g., for quaternions) are ignored.
/// @ingroup mocoutil
OSIMMOCO_API
std::vector<std::string> createStateVariableNamesInSystemOrder(
        const Model& model);

#ifndef SWIG
/// Same as above, but you can obtain a map from the returned state variable
/// names to the index in SimTK::State::getY() that accounts for empty slots
/// in Y.
/// @ingroup mocoutil
OSIMMOCO_API
std::vector<std::string> createStateVariableNamesInSystemOrder(
        const Model& model, std::unordered_map<int, int>& yIndexMap);

/// The map provides the index of each state variable in
/// SimTK::State::getY() from its state variable path string.
/// @ingroup mocoutil
OSIMMOCO_API
std::unordered_map<std::string, int> createSystemYIndexMap(const Model& model);
#endif

/// Create a vector of control names based on the actuators in the model for
/// which appliesForce == True. For actuators with one control (e.g.
/// ScalarActuator) the control name is simply the actuator name. For actuators
/// with multiple controls, each control name is the actuator name appended by
/// the control index (e.g. "/actuator_0"); modelControlIndices has length equal
/// to the number of controls associated with actuators that apply a force
/// (appliesForce == True). Its elements are the indices of the controls in the
/// Model::updControls() that are associated with actuators that apply a force.
/// @ingroup mocoutil
OSIMMOCO_API
std::vector<std::string> createControlNamesFromModel(
        const Model& model, std::vector<int>& modelControlIndices);
/// Same as above, but when there is no mapping to the modelControlIndices.
/// @ingroup mocoutil
OSIMMOCO_API
std::vector<std::string> createControlNamesFromModel(const Model& model);
/// The map provides the index of each control variable in the SimTK::Vector
/// returned by Model::getControls(), using the control name as the
/// key.
/// @throws Exception if the order of actuators in the model does not match
///     the order of controls in Model::getControls(). This is an internal
///     error, but you may be able to avoid the error by ensuring all Actuator%s
///     are in the Model's ForceSet.
/// @ingroup mocoutil
OSIMMOCO_API
std::unordered_map<std::string, int> createSystemControlIndexMap(
        const Model& model);

/// Throws an exception if the order of the controls in the model is not the
/// same as the order of the actuators in the model.
/// @ingroup mocoutil
OSIMMOCO_API void checkOrderSystemControls(const Model& model);

/// Throws an exception if the same label appears twice in the list of labels.
/// The argument copies the provided labels since we need to sort them to check
/// for redundancies.
/// @ingroup mocoutil
OSIMMOCO_API void checkRedundantLabels(std::vector<std::string> labels);

/// Throws an exception if any label in the provided list does not match any
/// state variable names in the model.
OSIMMOCO_API void checkLabelsMatchModelStates(const Model& model,
        const std::vector<std::string>& labels);

/// Convert a trajectory covering half the period of a symmetric motion into a
/// trajectory over the full period. This is useful for simulations of half a
/// gait cycle.
/// This converts time, states, controls, and derivatives; all other quantities
/// from the input trajectory are ignored.
/// If a column in the trajectory does not match addPatterns, negatePatterns,
/// negateAndShiftPatterns, or symmetryPatterns, then the second half of the
/// period contains the same data as the first half.
///
/// @param halfPeriodTrajectory The input trajectory covering half a period.
/// @param addPatterns If a column label matches an addPattern, then the second
/// half of the period for that column is (first_half_trajectory +
/// half_period_value - initial_value).
/// @param negatePatterns If a column label matches a negatePattern, then the
/// second half of the period for that column is (-first_half_trajectory).
/// This is usually relevant for only 3D models.
/// @param negateAndShiftPatterns If a column label matches a
/// negateAndShiftPattern, then the second half of the period for that column is
/// (-first_half_trajectory + 2 * half_period_value). This is usually relevant
/// for only 3D models.
/// @param symmetryPatterns This argument is a list of pairs, where the first
/// element of the pair is a pattern to match, and the second is a substitution
/// to convert the column label into the opposite column label of the symmetric
/// pair. If a column label matches a symmetryPattern, then its first
/// half-period is copied into the second half of the period for the column
/// identified by the substitution.
///
/// The default values for the patterns are intended to handle the column labels
/// for typical 2D or 3D OpenSim gait models.
/// The default values for negatePatterns, negateAndShiftPatterns, and
/// symmetryPatterns warrant an explanation. The string pattern before the
/// regex "(?!/value)" is followed by
/// anything except "/value" since it is contained in the negative lookahead
/// "(?!...)".  R"()" is a string literal that permits us to not escape
/// backslash characters. The regex "_r(\/|_|$)" matches "_r" followed by either
/// a forward slash (which is escaped), an underscore, OR the end of the string
/// ($). Since the forward slash and end of the string are within parentheses,
/// whatever matches this is captured and is available in the substitution (the
/// second element of the pair) as $1. The default symmetry patterns cause the
/// following replacements:
/// - "/jointset/hip_r/hip_flexion_r/value" becomes "/jointset/hip_l/hip_flexion_l/value"
/// - "/forceset/soleus_r" becomes "/forceset/soleus_l"
/// @ingroup mocoutil
OSIMMOCO_API MocoTrajectory createPeriodicTrajectory(
        const MocoTrajectory& halfPeriodTrajectory,
        std::vector<std::string> addPatterns = {".*pelvis_tx/value"},
        std::vector<std::string> negatePatterns = {
                                            ".*pelvis_list(?!/value).*",
                                            ".*pelvis_rotation.*",
                                            ".*pelvis_tz(?!/value).*",
                                            ".*lumbar_bending(?!/value).*",
                                            ".*lumbar_rotation.*"},
        std::vector<std::string> negateAndShiftPatterns = {
                                                   ".*pelvis_list/value",
                                                   ".*pelvis_tz/value",
                                                   ".*lumbar_bending/value"},
        std::vector<std::pair<std::string, std::string>> symmetryPatterns =
                {{R"(_r(\/|_|$))", "_l$1"}, {R"(_l(\/|_|$))", "_r$1"}});

/// This obtains the value of the OPENSIM_MOCO_PARALLEL environment variable.
/// The value has the following meanings:
/// - 0: run in series (not parallel).
/// - 1: run in parallel using all cores.
/// - greater than 1: run in parallel with this number of parallel jobs.
/// If the environment variable is not set, this function returns -1.
///
/// This variable does not indicate which calculations are parallelized
/// or how the parallelization is achieved. Moco may even ignore or override
/// the setting from the environment variable. See documentation elsewhere
/// (e.g., from a specific MocoSolver) for more information.
/// @ingroup mocoutil
OSIMMOCO_API int getMocoParallelEnvironmentVariable();

/// Thrown by FileDeletionThrower::throwIfDeleted().
/// @ingroup mocoutil
class FileDeletionThrowerException : public Exception {
public:
    FileDeletionThrowerException(const std::string& file, size_t line,
            const std::string& func, const std::string& deletedFile)
            : Exception(file, line, func) {
        addMessage("File '" + deletedFile + "' deleted.");
    }
};
/// This class helps a user cause an exception within the code. The constructor
/// writes a file, and the destructor deletes the file. The programmer can call
/// throwIfDeleted() to throw the FileDeletionThrowerException exception if the
/// file is deleted (by a user) before the object is destructed. If the file
/// could not be written by the constructor, then throwIfDeleted() does not
/// throw an exception.
/// @ingroup mocoutil
class FileDeletionThrower {
public:
    FileDeletionThrower()
            : FileDeletionThrower(
                      "OpenSimMoco_delete_this_to_throw_exception_" +
                      getFormattedDateTime() + ".txt") {}
    FileDeletionThrower(std::string filepath)
            : m_filepath(std::move(filepath)) {
        std::ofstream f(m_filepath);
        m_wroteInitialFile = f.good();
        f.close();
    }
    ~FileDeletionThrower() {
        if (m_wroteInitialFile) {
            std::ifstream f(m_filepath);
            if (f.good()) {
                f.close();
                std::remove(m_filepath.c_str());
            }
        }
    }
    void throwIfDeleted() const {
        if (m_wroteInitialFile) {
            OPENSIM_THROW_IF(!std::ifstream(m_filepath).good(),
                    FileDeletionThrowerException, m_filepath);
        }
    }

private:
    bool m_wroteInitialFile = false;
    const std::string m_filepath;
};

/// Obtain the ground reaction forces, centers of pressure, and torques
/// resulting from Force elements (e.g., SmoothSphereHalfSpaceForce), using a
/// model and states trajectory. Forces and torques are expressed in the ground
/// frame with respect to the ground origin. Hence, the centers of pressure are
/// at the origin. Paths to Force elements should be provided separately for
/// elements of the right and left feet. The output is a table formatted for use
/// with OpenSim tools; the labels of the columns distinguish between right
/// ("<>_r") and left ("<>_l") forces, centers of pressure, and torques. The
/// forces and torques used are taken from the first six outputs of
/// getRecordValues(); this order is of use for, for example, the
/// SmoothSphereHalfSpaceForce contact model but might have a different meaning
/// for different contact models.
/// @ingroup mocoutil
OSIMMOCO_API
TimeSeriesTable createExternalLoadsTableForGait(Model model,
        const StatesTrajectory& trajectory,
        const std::vector<std::string>& forcePathsRightFoot,
        const std::vector<std::string>& forcePathsLeftFoot);

/// Same as above, but with a MocoTrajectory instead of a StatesTrajectory.
/// @ingroup mocoutil
OSIMMOCO_API
TimeSeriesTable createExternalLoadsTableForGait(Model model,
        const MocoTrajectory& trajectory,
        const std::vector<std::string>& forcePathsRightFoot,
        const std::vector<std::string>& forcePathsLeftFoot);

} // namespace OpenSim

#endif // OPENSIM_MOCOUTILITIES_H