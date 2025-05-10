/* -------------------------------------------------------------------------- *
 *                          OpenSim:  testForces.cpp                          *
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

//==============================================================================
//
//  Tests Include:
//      1. PointToPointSpring
//      2. BushingForce
//      3. ElasticFoundationForce
//      4. HuntCrossleyForce
//      5. SmoothSphereHalfSpaceForce
//      6. CoordinateLimitForce
//      7. RotationalCoordinateLimitForce
//      8. ExternalForce
//      9. PathSpring
//     10. ExpressionBasedPointToPointForce
//     11. Blankevoort1991Ligament
//     12. Smith2018ArticularContactForce
//
//     Add tests here as Forces are added to OpenSim
//
//==============================================================================
#include "SimTKcommon/internal/Xml.h"
#include <ctime> // clock(), clock_t, CLOCKS_PER_SEC

#include <OpenSim/Analyses/osimAnalyses.h>
#include <OpenSim/Auxiliary/auxiliaryTestFunctions.h>
#include <OpenSim/Simulation/osimSimulation.h>
#include <OpenSim/Actuators/osimActuators.h>
 #include<math.h>

using namespace OpenSim;
using namespace std;

//==============================================================================
// Common Parameters for the simulations are just global.
const static double integ_accuracy = 1.0e-4;
const static SimTK::Vec3 gravity_vec = SimTK::Vec3(0, -9.8065, 0);
//==============================================================================

void testPathSpring();
void testExternalForce();
void testSpringMass();
void testBushingForce();
void testTwoFrameLinkerUpdateFromXMLNode();
void testFunctionBasedBushingForce();
void testExpressionBasedBushingForceTranslational();
void testExpressionBasedBushingForceRotational();
void testElasticFoundation();
void testHuntCrossleyForce();
void testSmoothSphereHalfSpaceForce();
void testCoordinateLimitForce();
void testCoordinateLimitForceRotational();
void testExpressionBasedPointToPointForce();
void testExpressionBasedCoordinateForce();
void testSerializeDeserialize();
void testTranslationalDampingEffect(Model& osimModel, Coordinate& sliderCoord,
        double start_h, Component& componentWithDamping);
void testBlankevoort1991Ligament();
void testSmith2018ArticularContactForce();
void testSmith2018ContactDetection();

int main() {
    SimTK::Array_<std::string> failures;
    /*
    try { testPathSpring(); }
    catch (const std::exception& e){
        cout << e.what() <<endl; failures.push_back("testPathSpring");
    }

    try { testExternalForce(); }
    catch (const std::exception& e){
        cout << e.what() <<endl; failures.push_back("testExternalForce");
    }

    try { testSpringMass(); }
    catch (const std::exception& e){
        cout << e.what() <<endl; failures.push_back("testP2PSpringMass");
    }

    try { testBushingForce(); }
    catch (const std::exception& e){
        cout << e.what() <<endl; failures.push_back("testBushingForce");
    }

    try { testTwoFrameLinkerUpdateFromXMLNode(); }
    catch (const std::exception& e){
        cout << e.what() <<endl;
        failures.push_back("testTwoFrameLinkerUpdateFromXMLNode");
    }

    try { testFunctionBasedBushingForce(); }
    catch (const std::exception& e){
        cout << e.what() <<endl;
        failures.push_back("testFunctionBasedBushingForce");
    }

    try { testExpressionBasedBushingForceTranslational(); }
    catch (const std::exception& e){
        cout << e.what() <<endl;
        failures.push_back("testExpressionBasedBushingForceTranslational");
    }

    try { testExpressionBasedBushingForceRotational(); }
    catch (const std::exception& e) {
        cout << e.what() << endl;
        failures.push_back("testExpressionBasedBushingForceRotational");
    }

    try { testElasticFoundation(); }
    catch (const std::exception& e){
        cout << e.what() <<endl; failures.push_back("testElasticFoundation");
    }

    try { testHuntCrossleyForce(); }
    catch (const std::exception& e){
        cout << e.what() <<endl; failures.push_back("testHuntCrossleyForce");
    }

    try { testSmoothSphereHalfSpaceForce(); }
    catch (const std::exception& e){
        cout << e.what() <<endl;
        failures.push_back("testSmoothSphereHalfSpaceForce");
    }

    try { testCoordinateLimitForce(); }
    catch (const std::exception& e){
        cout << e.what() <<endl; failures.push_back("testCoordinateLimitForce");
    }

    try { testCoordinateLimitForceRotational(); }
    catch (const std::exception& e){
        cout << e.what() <<endl;
        failures.push_back("testCoordinateLimitForceRotational");
    }

    try { testExpressionBasedPointToPointForce(); }
    catch (const std::exception& e){
        cout << e.what() <<endl;
        failures.push_back("testExpressionBasedPointToPointForce");
    }

    try { testExpressionBasedCoordinateForce(); }
    catch (const std::exception& e){
        cout << e.what() <<endl;
        failures.push_back("testExpressionBasedCoordinateForce");
    }

    try { testSerializeDeserialize(); }
    catch (const std::exception& e){
        cout << e.what() <<endl;
        failures.push_back("testSerializeDeserialize");
    }
    
    try {
        testBlankevoort1991Ligament();
    } catch (const std::exception& e) {
        cout << e.what() << endl;
        failures.push_back("testBlankevoort1991Ligament");
    }
    */
    try {
        testSmith2018ArticularContactForce();
    }
    /*
    try {
        testSmith2018ContactDetection();
    }*/
    catch (const std::exception& e) {
        cout << e.what() << endl;
        failures.push_back("testSmith2018ArticularContactForce");
    }

    if (!failures.empty()) {
        cout << "Done, with failure(s): " << failures << endl;
        return 1;
    }

    cout << "Done. All cases passed." << endl;

    return 0;
}

//==============================================================================
// Test Cases
//==============================================================================

void testExpressionBasedCoordinateForce() {
    using namespace SimTK;

    double mass = 1;
    double stiffness = 10;
    double damp_coeff = 5;
    double start_h = 0.5;
    double start_v = 0;
    double ball_radius = 0.25;

    double omega = sqrt(stiffness / mass);
    // note: test case designed for 0 <= zeta < 1 (under damped system)
    double zeta = damp_coeff / (2 * sqrt(mass * stiffness));
    double damp_freq = omega * sqrt(1 - pow(zeta, 2));

    double dh = mass * gravity_vec(1) / stiffness;

    // Setup OpenSim model
    Model osimModel{};
    osimModel.setName("SpringMass");
    // OpenSim bodies
    const Ground& ground = osimModel.getGround();
    ;
    OpenSim::Body ball(
            "ball", mass, Vec3(0), mass * SimTK::Inertia::sphere(0.1));
    ball.attachGeometry(new Sphere(0.1));
    ball.scale(Vec3(ball_radius), false);

    // Add joints
    SliderJoint slider("slider", ground, Vec3(0), Vec3(0, 0, Pi / 2), ball,
            Vec3(0), Vec3(0, 0, Pi / 2));

    double positionRange[2] = {-10, 10};
    // Rename coordinate for the slider joint
    auto& sliderCoord = slider.updCoordinate();
    sliderCoord.setName("ball_h");
    sliderCoord.setRange(positionRange);

    osimModel.addBody(&ball);
    osimModel.addJoint(&slider);

    osimModel.setGravity(gravity_vec);

    // ode for basic mass-spring-dampener system
    ExpressionBasedCoordinateForce spring("ball_h", "-10*q-5*qdot");

    osimModel.addForce(&spring);

    // Create the force reporter
    ForceReporter* reporter = new ForceReporter(&osimModel);
    osimModel.addAnalysis(reporter);

    SimTK::State& osim_state = osimModel.initSystem();

    // move ball to initial conditions
    sliderCoord.setValue(osim_state, start_h);
    osimModel.getMultibodySystem().realize(osim_state, Stage::Position);

    //==========================================================================
    // Compute the force at the specified times.

    double final_t = 2.0;
    double nsteps = 10;
    double dt = final_t / nsteps;

    Manager manager(osimModel);
    manager.setIntegratorAccuracy(1e-7);
    osim_state.setTime(0.0);
    manager.initialize(osim_state);

    for (int i = 1; i <= nsteps; i++) {
        osim_state = manager.integrate(dt * i);
        osimModel.getMultibodySystem().realize(osim_state, Stage::Acceleration);
        Vec3 pos = ball.findStationLocationInGround(osim_state, Vec3(0));

        double height =
                exp(-1 * zeta * omega * osim_state.getTime()) *
                        ((start_h - dh) *
                                        cos(damp_freq * osim_state.getTime()) +
                                ((1 / damp_freq) *
                                        (zeta * omega * (start_h - dh) +
                                                start_v) *
                                        sin(damp_freq *
                                                osim_state.getTime()))) +
                dh;

        ASSERT_EQUAL(height, pos(1), 1e-6);
    }

    // Test copying
    ExpressionBasedCoordinateForce* copyOfSpring = spring.clone();

    ASSERT(*copyOfSpring == spring);

    osimModel.print("ExpressionBasedCoordinateForceModel.osim");

    osimModel.disownAllComponents();
}

void testExpressionBasedPointToPointForce() {
    using namespace SimTK;

    double mass = 100;
    double ball_radius = 0.25;

    Random::Uniform rand;
    Vec3 p1(rand.getValue(), rand.getValue(), rand.getValue());
    Vec3 p2(rand.getValue(), rand.getValue(), rand.getValue());

    // Setup OpenSim model
    Model model{};
    model.setName("ExpressionBasedPointToPointForce");
    // OpenSim bodies
    const Ground& ground = model.getGround();
    OpenSim::Body ball(
            "ball", mass, Vec3(0), mass * SimTK::Inertia::sphere(ball_radius));
    ball.attachGeometry(new Sphere(ball_radius));
    ball.scale(Vec3(ball_radius), false);

    // define body's joint
    FreeJoint free("free", ground, Vec3(0), Vec3(0, 0, Pi / 2), ball, Vec3(0),
            Vec3(0, 0, Pi / 2));

    model.addBody(&ball);
    model.addJoint(&free);

    string expression("2/(d^2)-3.0*(d-0.2)*(1+0.0123456789*ddot)");

    ExpressionBasedPointToPointForce* p2pForce =
            new ExpressionBasedPointToPointForce(
                    "ground", p1, "ball", p2, expression);
    p2pForce->setName("P2PTestForce");

    model.addForce(p2pForce);

    // Create the force reporter
    ForceReporter* reporter = new ForceReporter(&model);
    model.addAnalysis(reporter);

    // model.setUseVisualizer(true);
    SimTK::State& state = model.initSystem();

    model.print("ExpressionBasedPointToPointForceModel.osim");

    Vector& q = state.updQ();
    Vector& u = state.updU();

    for (int i = 0; i < state.getNU(); ++i) {
        q[i] = rand.getValue();
        u[i] = rand.getValue();
    }

    //==========================================================================
    // Compute the force and torque at the specified times.
    Manager manager(model);
    manager.setIntegratorAccuracy(1e-6);
    state.setTime(0.0);
    manager.initialize(state);

    double final_t = 1.0;
    state = manager.integrate(final_t);

    // manager.getStateStorage().print("testExpressionBasedPointToPointForce.sto");

    // force is only velocity dependent but is only compute in Dynamics
    model.getMultibodySystem().realize(state, Stage::Dynamics);

    // Now check that the force reported by spring
    double model_force = p2pForce->getForceMagnitude(state);

    // Save the forces
    // reporter->getForceStorage().print("path_spring_forces.mot");
    double d = (p1 - ball.findStationLocationInGround(state, p2)).norm();
    const MobilizedBody& b1 = ground.getMobilizedBody();
    const MobilizedBody& b2 = ball.getMobilizedBody();

    double ddot =
            b1.calcStationToStationDistanceTimeDerivative(state, p1, b2, p2);

    // string expression("2/(d^2)-3.0*(d-0.2)*(1+0.0123456789*ddot)");
    double analytical_force =
            2 / (d * d) - 3.0 * (d - 0.2) * (1 + 0.0123456789 * ddot);

    // something is wrong if the block does not reach equilibrium
    ASSERT_EQUAL(analytical_force, model_force, 1e-5);

    // Before exiting lets see if copying the P2P force works
    ExpressionBasedPointToPointForce* copyOfP2pForce = p2pForce->clone();
    ASSERT(*copyOfP2pForce == *p2pForce);

    model.disownAllComponents();
}

void testPathSpring() {
    using namespace SimTK;

    double mass = 1;
    double stiffness = 10;
    double restlength = 0.5;
    double dissipation = 0.1;
    double start_h = 0.5;

    // Setup OpenSim model
    Model osimModel{};
    osimModel.setName("PathSpring");
    // OpenSim bodies
    const Ground& ground = osimModel.getGround();
    ;
    OpenSim::Body pulleyBody("PulleyBody", mass, Vec3(0),
            mass * SimTK::Inertia::brick(0.1, 0.1, 0.1));
    OpenSim::Body block("block", mass, Vec3(0),
            mass * SimTK::Inertia::brick(0.2, 0.1, 0.1));
    block.attachGeometry(new Brick(Vec3(0.2, 0.1, 0.1)));
    block.scale(Vec3(0.2, 0.1, 0.1), false);

    WrapCylinder* pulley = new WrapCylinder();
    pulley->set_radius(0.1);
    pulley->set_length(0.05);

    // Add the wrap object to the body, which takes ownership of it
    pulleyBody.addWrapObject(pulley);

    // Add joints
    WeldJoint weld("pulley", ground, Vec3(0, 1.0, 0), Vec3(0), pulleyBody,
            Vec3(0), Vec3(0));
    SliderJoint slider("slider", ground, Vec3(0), Vec3(0, 0, Pi / 2), block,
            Vec3(0), Vec3(0, 0, Pi / 2));

    double positionRange[2] = {-10, 10};
    // Rename coordinate for the slider joint
    auto& sliderCoord = slider.updCoordinate();
    sliderCoord.setName("block_h");
    sliderCoord.setRange(positionRange);

    osimModel.addBody(&block);
    osimModel.addJoint(&weld);

    osimModel.addBody(&pulleyBody);
    osimModel.addJoint(&slider);

    osimModel.setGravity(gravity_vec);

    PathSpring spring("spring", restlength, stiffness, dissipation);
    spring.updGeometryPath().appendNewPathPoint(
            "origin", block, Vec3(-0.1, 0.0, 0.0));

    int N = 10;
    for (int i = 1; i < N; ++i) {
        double angle = i * Pi / N;
        double x = 0.1 * cos(angle);
        double y = 0.1 * sin(angle);
        spring.updGeometryPath().appendNewPathPoint(
                "", pulleyBody, Vec3(-x, y, 0.0));
    }

    spring.updGeometryPath().appendNewPathPoint(
            "insertion", block, Vec3(0.1, 0.0, 0.0));

    // BUG in defining wrapping API requires that the Force containing the
    // GeometryPath be connected to the model before the wrap can be added
    osimModel.addForce(&spring);

    // Create the force reporter
    ForceReporter* reporter = new ForceReporter(&osimModel);
    osimModel.addAnalysis(reporter);

    // osimModel.setUseVisualizer(true);
    SimTK::State& osim_state = osimModel.initSystem();

    sliderCoord.setValue(osim_state, start_h);
    osimModel.getMultibodySystem().realize(osim_state, Stage::Position);

    //==========================================================================
    // Compute the force and torque at the specified times.
    Manager manager(osimModel);
    manager.setIntegratorAccuracy(1e-6);
    osim_state.setTime(0.0);
    manager.initialize(osim_state);

    double final_t = 10.0;
    osim_state = manager.integrate(final_t);

    // tension should only be velocity dependent
    osimModel.getMultibodySystem().realize(osim_state, Stage::Velocity);

    // Now check that the force reported by spring
    double model_force = spring.getTension(osim_state);

    // get acceleration of the block
    osimModel.getMultibodySystem().realize(osim_state, Stage::Acceleration);
    double hddot =
            osimModel.getCoordinateSet().get("block_h").getAccelerationValue(
                    osim_state);

    // the tension should be half the weight of the block
    double analytical_force = -0.5 * (gravity_vec(1) - hddot) * mass;

    // Save the forces
    reporter->getForceStorage().print("path_spring_forces.mot");

    // something is wrong if the block does not reach equilibrium
    ASSERT_EQUAL(analytical_force, model_force, 1e-3);

    // Before exiting lets see if copying the spring works
    PathSpring* copyOfSpring = spring.clone();
    ASSERT(*copyOfSpring == spring);

    osimModel.disownAllComponents();
}

void testSpringMass() {
    using namespace SimTK;

    double mass = 1;
    double stiffness = 10;
    double restlength = 1.0;
    double start_h = 0.5;
    double ball_radius = 0.25;

    double omega = sqrt(stiffness / mass);

    double dh = mass * gravity_vec(1) / stiffness;

    // Setup OpenSim model
    Model osimModel{};
    osimModel.setName("SpringMass");
    // OpenSim bodies
    const Ground& ground = osimModel.getGround();
    ;
    OpenSim::Body ball(
            "ball", mass, Vec3(0), mass * SimTK::Inertia::sphere(0.1));
    ball.attachGeometry(new Sphere(0.1));
    ball.scale(Vec3(ball_radius), false);

    // Add joints
    SliderJoint slider("slider", ground, Vec3(0), Vec3(0, 0, Pi / 2), ball,
            Vec3(0), Vec3(0, 0, Pi / 2));

    double positionRange[2] = {-10, 10};
    // Rename coordinate for the slider joint
    auto& sliderCoord = slider.updCoordinate();
    sliderCoord.setName("ball_h");
    sliderCoord.setRange(positionRange);

    osimModel.addBody(&ball);
    osimModel.addJoint(&slider);

    osimModel.setGravity(gravity_vec);

    PointToPointSpring spring(osimModel.updGround(), Vec3(0., restlength, 0.),
            ball, Vec3(0.), stiffness, restlength);

    osimModel.addForce(&spring);

    // Create the force reporter
    ForceReporter* reporter = new ForceReporter(&osimModel);
    osimModel.addAnalysis(reporter);

    SimTK::State& osim_state = osimModel.initSystem();

    sliderCoord.setValue(osim_state, start_h);
    osimModel.getMultibodySystem().realize(osim_state, Stage::Position);

    //==========================================================================
    // Compute the force and torque at the specified times.

    Manager manager(osimModel);
    manager.setIntegratorAccuracy(1e-7);
    osim_state.setTime(0.0);
    manager.initialize(osim_state);

    double final_t = 2.0;
    double nsteps = 10;
    double dt = final_t / nsteps;

    for (int i = 1; i <= nsteps; i++) {
        osim_state = manager.integrate(dt * i);
        osimModel.getMultibodySystem().realize(osim_state, Stage::Acceleration);
        Vec3 pos = ball.findStationLocationInGround(osim_state, Vec3(0));

        double height = (start_h - dh) * cos(omega * osim_state.getTime()) + dh;
        ASSERT_EQUAL(height, pos(1), 1e-5);

        // Now check that the force reported by spring
        Array<double> model_force = spring.getRecordValues(osim_state);

        // get the forces applied to the ground and ball
        double analytical_force = -stiffness * height;
        // analytical force corresponds in direction to the force on the ball Y
        // index = 7
        ASSERT_EQUAL(analytical_force, model_force[7], 1e-4);
    }

    // Save the forces
    osimModel.disownAllComponents();

    // Before exiting lets see if copying the spring works
    PointToPointSpring* copyOfSpring = spring.clone();

    ASSERT(*copyOfSpring == spring);

    // Verify that the PointToPointSpring is correctly deserialized from
    // previous major version of OpenSim.
    Model bouncer("bouncing_block_30000.osim");
    SimTK::State& s = bouncer.initSystem();
    bouncer.realizeAcceleration(s);

    /*Vec3 comA =*/bouncer.calcMassCenterAcceleration(s);
}

void testBushingForce() {
    using namespace SimTK;

    double mass = 1;
    double stiffness = 10;
    double start_h = 0.5;
    double ball_radius = 0.25;

    double omega = sqrt(stiffness / mass);

    double dh = mass * gravity_vec(1) / stiffness;

    // Setup OpenSim model
    Model osimModel{};
    osimModel.setName("BushingTest");
    // OpenSim bodies
    const Ground& ground = osimModel.getGround();
    ;
    auto* ball = new OpenSim::Body(
            "ball", mass, Vec3(0), mass * SimTK::Inertia::sphere(0.1));
    ball->attachGeometry(new Sphere{0.1});
    ball->scale(Vec3(ball_radius), false);

    // Add joints
    auto* slider = new SliderJoint("slider", ground, Vec3(0),
            Vec3(0, 0, Pi / 2), *ball, Vec3(0), Vec3(0, 0, Pi / 2));

    double positionRange[2] = {-10, 10};
    // Rename coordinate for the slider joint
    auto& sliderCoord = slider->updCoordinate();
    sliderCoord.setName("ball_h");
    sliderCoord.setRange(positionRange);

    osimModel.addBody(ball);
    osimModel.addJoint(slider);

    Vec3 rotStiffness(0);
    Vec3 transStiffness(stiffness);
    Vec3 rotDamping(0);
    Vec3 transDamping(0);

    osimModel.setGravity(gravity_vec);

    auto* spring = new BushingForce("bushing", ground, *ball);
    spring->set_translational_stiffness(transStiffness);
    spring->set_rotational_stiffness(rotStiffness);
    spring->set_translational_damping(transDamping);
    spring->set_rotational_damping(rotDamping);

    osimModel.addForce(spring);
    const BushingForce& bushingForce =
            osimModel.getComponent<BushingForce>("forceset/bushing");

    // To print (serialize) the latest connections of the model, it is
    // necessary to finalizeConnections() first.
    osimModel.finalizeConnections();
    osimModel.print("BushingForceModel.osim");

    Model previousVersionModel("BushingForceModel_30000.osim");
    previousVersionModel.print("BushingForceModel_30000_in_Latest.osim");

    const BushingForce& bushingForceFromPrevious =
            previousVersionModel.getComponent<BushingForce>("forceset/bushing");

    ASSERT(bushingForce == bushingForceFromPrevious, __FILE__, __LINE__,
            "current bushing force FAILED to match bushing force from previous "
            "model.");

    // Create the force reporter
    ForceReporter* reporter = new ForceReporter(&osimModel);
    osimModel.addAnalysis(reporter);

    SimTK::State& osim_state = osimModel.initSystem();

    sliderCoord.setValue(osim_state, start_h);
    osimModel.getMultibodySystem().realize(osim_state, Stage::Position);

    //==========================================================================
    // Compute the force and torque at the specified times.
    Manager manager(osimModel);
    manager.setIntegratorAccuracy(1e-6);
    osim_state.setTime(0.0);
    manager.initialize(osim_state);

    double final_t = 2.0;
    double nsteps = 10;
    double dt = final_t / nsteps;

    for (int i = 1; i <= nsteps; i++) {
        osim_state = manager.integrate(dt * i);
        osimModel.getMultibodySystem().realize(osim_state, Stage::Acceleration);
        Vec3 pos = ball->findStationLocationInGround(osim_state, Vec3(0));

        double height = (start_h - dh) * cos(omega * osim_state.getTime()) + dh;
        ASSERT_EQUAL(height, pos(1), 1e-4);

        // Now check that the force reported by spring
        Array<double> model_force = spring->getRecordValues(osim_state);

        // get the forces applied to the ground and ball
        double analytical_force = -stiffness * height;
        // analytical force corresponds in direction to the force on the ball Y
        // index = 7
        ASSERT_EQUAL(analytical_force, model_force[7], 2e-4);
    }

    manager.getStateStorage().print("bushing_model_states.sto");

    // Save the forces
    reporter->getForceStorage().print("bushing_forces.mot");

    // Before exiting lets see if copying the spring works
    BushingForce* copyOfSpring = spring->clone();

    ASSERT(*copyOfSpring == *spring);
}

// testBushingForce() performs similar checks as does this test, but this test
// ensures intermediate offset frames are created correctly. This test still
// uses BushingForce to test the TwoFrameLinker.
void testTwoFrameLinkerUpdateFromXMLNode() {
    using namespace SimTK;

    double mass = 1;
    double start_h = 0.5;
    double ball_radius = 0.25;

    // Setup OpenSim model
    Model osimModel{};
    osimModel.setName("TwoFrameLinkerUpdateFromXMLNodeTest");
    // OpenSim bodies
    const Ground& ground = osimModel.getGround();
    ;
    auto* ball = new OpenSim::Body(
            "ball", mass, Vec3(0), mass * SimTK::Inertia::sphere(0.1));
    ball->attachGeometry(new Sphere{0.1});
    ball->scale(Vec3(ball_radius), false);

    // Add joints
    auto* slider = new SliderJoint("slider", ground, Vec3(0),
            Vec3(0, 0, Pi / 2), *ball, Vec3(0), Vec3(0, 0, Pi / 2));

    double positionRange[2] = {-10, 10};
    // Rename coordinate for the slider joint
    auto& sliderCoord = slider->updCoordinate();
    sliderCoord.setName("ball_h");
    sliderCoord.setRange(positionRange);

    osimModel.addBody(ball);
    osimModel.addJoint(slider);

    Vec3 rotStiffness(15, 21, 30);
    Vec3 transStiffness(10, 10, 10);
    Vec3 rotDamping(0.4, 0.5, 0.6);
    Vec3 transDamping(0.2, 0.3, 0.4);

    osimModel.setGravity(gravity_vec);

    auto* spring = new BushingForce("bushing", ground,
            Transform(Rotation(BodyRotationSequence, -0.5, XAxis, 0, YAxis, 0.5,
                              ZAxis),
                    Vec3(1, 2, 3)),
            *ball,
            Transform(Rotation(BodyRotationSequence, 0.1, XAxis, 0.2, YAxis,
                              0.3, ZAxis),
                    Vec3(4, 5, 6)),
            transStiffness, rotStiffness, transDamping, rotDamping);

    spring->print("bushingForceAPICreated.xml");

    osimModel.addForce(spring);
    const BushingForce& bushingForce =
            osimModel.getComponent<BushingForce>("./forceset/bushing");

    // It's necessary to correct the connectee paths in the BushingForce, which
    // we can do with finalizeConnections() (they are incorrect otherwise
    // because `spring` is initially orphaned).
    osimModel.finalizeConnections();
    osimModel.print("BushingForceOffsetModel.osim");

    Model previousVersionModel("BushingForceOffsetModel_30000.osim");
    previousVersionModel.finalizeConnections();
    previousVersionModel.print("BushingForceOffsetModel_30000_in_Latest.osim");

    const BushingForce& bushingForceFromPrevious =
            previousVersionModel.getComponent<BushingForce>(
                    "./forceset/bushing");

    ASSERT(bushingForce == bushingForceFromPrevious, __FILE__, __LINE__,
            "current bushing force FAILED to match bushing force from previous "
            "model.");
}

void testFunctionBasedBushingForce() {
    using namespace SimTK;

    double mass = 1;
    double stiffness = 10;
    double start_h = 0.5;
    double ball_radius = 0.25;

    double omega = sqrt(stiffness / mass);

    double dh = mass * gravity_vec(1) / stiffness;

    // Setup OpenSim model
    Model osimModel{};
    osimModel.setName("FunctionBasedBushingTest");
    // OpenSim bodies
    const Ground& ground = osimModel.getGround();
    ;
    OpenSim::Body ball(
            "ball", mass, Vec3(0), mass * SimTK::Inertia::sphere(0.1));
    ball.attachGeometry(new Sphere(0.1));
    ball.scale(Vec3(ball_radius), false);

    // Add joints
    SliderJoint slider("slider", ground, Vec3(0), Vec3(0, 0, Pi / 2), ball,
            Vec3(0), Vec3(0, 0, Pi / 2));

    double positionRange[2] = {-10, 10};
    // Rename coordinate for the slider joint
    auto& sliderCoord = slider.updCoordinate();
    sliderCoord.setName("ball_h");
    sliderCoord.setRange(positionRange);

    osimModel.addBody(&ball);
    osimModel.addJoint(&slider);

    Vec3 rotStiffness(0);
    Vec3 transStiffness(stiffness);
    Vec3 rotDamping(0);
    Vec3 transDamping(0);

    osimModel.setGravity(gravity_vec);

    FunctionBasedBushingForce spring("linear_bushing", ground, Vec3(0), Vec3(0),
            ball, Vec3(0), Vec3(0), transStiffness, rotStiffness, transDamping,
            rotDamping);

    osimModel.addForce(&spring);
    osimModel.finalizeConnections(); // fix warning on write that results in
                                     // invalid model file
    osimModel.print("FunctionBasedBushingForceModel.osim");

    // Create the force reporter
    ForceReporter* reporter = new ForceReporter(&osimModel);
    osimModel.addAnalysis(reporter);

    SimTK::State& osim_state = osimModel.initSystem();

    sliderCoord.setValue(osim_state, start_h);
    osimModel.getMultibodySystem().realize(osim_state, Stage::Position);

    //==========================================================================
    // Compute the force and torque at the specified times.
    Manager manager(osimModel);
    manager.setIntegratorAccuracy(1e-6);
    osim_state.setTime(0.0);
    manager.initialize(osim_state);

    double final_t = 2.0;
    double nsteps = 10;
    double dt = final_t / nsteps;

    for (int i = 1; i <= nsteps; i++) {
        osim_state = manager.integrate(dt * i);
        osimModel.getMultibodySystem().realize(osim_state, Stage::Acceleration);
        Vec3 pos = ball.findStationLocationInGround(osim_state, Vec3(0));

        double height = (start_h - dh) * cos(omega * osim_state.getTime()) + dh;
        ASSERT_EQUAL(height, pos(1), 1e-4);

        // Now check that the force reported by spring
        Array<double> model_force = spring.getRecordValues(osim_state);

        // get the forces applied to the ground and ball
        double analytical_force = -stiffness * height;
        // analytical force corresponds in direction to the force on the ball Y
        // index = 7
        ASSERT_EQUAL(analytical_force, model_force[7], 2e-4);
    }

    manager.getStateStorage().print("function_based_bushing_model_states.sto");

    // Save the forces
    reporter->getForceStorage().print("function_based_bushing_forces.mot");
    // Now add damping to the bushing force and make sure energy is
    // monotonically decreasing
    testTranslationalDampingEffect(osimModel, sliderCoord, start_h, spring);
    // The following line is BAD but necessary due to mixing stack and heap
    // allocation
    osimModel.disownAllComponents();

    // Before exiting lets see if copying the spring works
    FunctionBasedBushingForce* copyOfSpring = spring.clone();

    ASSERT(*copyOfSpring == spring);
}

void testExpressionBasedBushingForceTranslational() {
    using namespace SimTK;

    double mass = 1;
    double stiffness = 10;
    double start_h = 0.5;
    double ball_radius = 0.25;

    double omega = sqrt(stiffness / mass);

    double dh = mass * gravity_vec(1) / stiffness;

    // Setup OpenSim model
    Model osimModel{};
    osimModel.setName("ExpressionBasedBushingTranslationTest");
    osimModel.setGravity(gravity_vec);

    // Create ball body and attach it to ground
    // with a vertical slider

    const Ground& ground = osimModel.getGround();

    OpenSim::Body ball(
            "ball", mass, Vec3(0), mass * SimTK::Inertia::sphere(0.1));
    ball.attachGeometry(new Sphere(0.1));
    ball.scale(Vec3(ball_radius), false);

    SliderJoint sliderY("slider", ground, Vec3(0), Vec3(0, 0, Pi / 2), ball,
            Vec3(0), Vec3(0, 0, Pi / 2));

    double positionRange[2] = {-10, 10};
    // Rename coordinate for the slider joint
    auto& sliderCoord = sliderY.updCoordinate();
    sliderCoord.setName("ball_h");
    sliderCoord.setRange(positionRange);

    osimModel.addBody(&ball);
    osimModel.addJoint(&sliderY);

    // Create base body and attach it to ground with a weld

    OpenSim::Body base(
            "base_body", mass, Vec3(0), mass * SimTK::Inertia::sphere(0.1));
    base.attachGeometry(new Sphere(0.1));
    base.scale(Vec3(ball_radius), false);

    WeldJoint weld("weld", ground, Vec3(0), Vec3(0), base, Vec3(0), Vec3(0));
    osimModel.addBody(&base);
    osimModel.addJoint(&weld);

    // create an ExpressionBasedBushingForce that represents an
    // uncoupled, linear bushing between the ball body and welded base body
    Vec3 rotStiffness(0);
    Vec3 transStiffness(stiffness);
    Vec3 rotDamping(0);
    Vec3 transDamping(0);

    ExpressionBasedBushingForce spring("linear_bushing", base, Vec3(0), Vec3(0),
            ball, Vec3(0), Vec3(0), transStiffness, rotStiffness, transDamping,
            rotDamping);

    spring.setName("translational_linear_bushing");

    osimModel.addForce(&spring);
    osimModel.finalizeConnections();
    osimModel.print("ExpressionBasedBushingForceTranslationalModel.osim");

    // Create the force reporter
    ForceReporter* reporter = new ForceReporter(&osimModel);
    osimModel.addAnalysis(reporter);

    SimTK::State& osim_state = osimModel.initSystem();

    // set the initial height of the ball on slider
    sliderCoord.setValue(osim_state, start_h);
    osimModel.getMultibodySystem().realize(osim_state, Stage::Position);

    //==========================================================================
    // Compute the force and torque at the specified times.
    Manager manager(osimModel);
    manager.setIntegratorAccuracy(1e-6);
    osim_state.setTime(0.0);
    manager.initialize(osim_state);

    double final_t = 2.0;
    double nsteps = 10;
    double dt = final_t / nsteps;

    for (int i = 1; i <= nsteps; ++i) {
        osim_state = manager.integrate(dt * i);
        osimModel.getMultibodySystem().realize(osim_state, Stage::Acceleration);
        Vec3 pos = ball.findStationLocationInGround(osim_state, Vec3(0));

        // compute the height based on the analytic solution for 1-D spring-mass
        // system with zero-velocity at initial offset.
        double height = (start_h - dh) * cos(omega * osim_state.getTime()) + dh;

        // check that the simulated solution is equivalent to the analytic
        // solution
        ASSERT_EQUAL(height, pos(1), 1e-4);

        // get the forces applied to the base and ball
        Array<double> model_force = spring.getRecordValues(osim_state);

        // compute the expected force on the ball
        double analytical_force = -stiffness * height;

        // check analytical force corresponds to the force on the ball
        // in the Y direction, index = 7
        ASSERT_EQUAL(analytical_force, model_force[7], 2e-4);
    }

    manager.getStateStorage().print(
            "expression_based_bushing_translational_model_states.sto");

    // Save the forces
    reporter->getForceStorage().print(
            "expression_based_bushing_translational_model_forces.mot");

    // Now add damping to the bushing force and make sure energy is
    // monotonically decreasing
    testTranslationalDampingEffect(osimModel, sliderCoord, start_h, spring);

    // The following line is BAD but necessary due to mixing stack and heap
    // allocation
    osimModel.disownAllComponents();

    // Before exiting lets see if copying the spring works
    ExpressionBasedBushingForce* copyOfSpring = spring.clone();

    ASSERT(*copyOfSpring == spring);
}

void testExpressionBasedBushingForceRotational() {
    using namespace SimTK;

    double mass = 5;
    double stiffness = 2;
    double start_theta = Pi / 6;
    double ball_radius = 0.25;

    // Setup OpenSim model
    Model osimModel{};
    osimModel.setName("ExpressionBasedBushingRotationalTest");
    const Ground& ground = osimModel.getGround();

    // Create base body and attach it to ground with a weld

    OpenSim::Body base("base_body", mass, Vec3(0),
            mass * SimTK::Inertia::sphere(ball_radius));

    WeldJoint weld("weld", ground, Vec3(0), Vec3(0), base, Vec3(0), Vec3(0));
    osimModel.addBody(&base);
    osimModel.addJoint(&weld);

    // Create ball body and attach it to ground
    // with a pin joint

    OpenSim::Body ball(
            "ball", mass, Vec3(0), mass * SimTK::Inertia::sphere(ball_radius));

    PinJoint pin("pin", ground, Vec3(0), Vec3(Pi / 2, 0, 0), ball, Vec3(0),
            Vec3(Pi / 2, 0, 0));

    double thetaRange[2] = {-2 * Pi, 2 * Pi};
    // Rename coordinate for the pin joint
    auto& pinCoord = pin.updCoordinate();
    pinCoord.setName("ball_theta");
    pinCoord.setRange(thetaRange);

    osimModel.addBody(&ball);
    osimModel.addJoint(&pin);

    // create an ExpressionBasedBushingForce that represents an
    // uncoupled, linear bushing between the ball body and welded base body
    Vec3 rotStiffness(stiffness);
    Vec3 transStiffness(0);
    Vec3 rotDamping(0);
    Vec3 transDamping(0);

    ExpressionBasedBushingForce spring("rotatinal_spring", base, Vec3(0),
            Vec3(0), ball, Vec3(0), Vec3(0), transStiffness, rotStiffness,
            transDamping, rotDamping);

    spring.setName("rotational_linear_bushing");

    osimModel.addForce(&spring);
    osimModel.finalizeConnections();
    osimModel.print("ExpressionBasedBushingForceRotationalModel.osim");

    // Create the force reporter
    ForceReporter* reporter = new ForceReporter(&osimModel);
    osimModel.addAnalysis(reporter);

    SimTK::State& osim_state = osimModel.initSystem();

    // set the initial pin joint angle
    pinCoord.setValue(osim_state, start_theta);
    osimModel.getMultibodySystem().realize(osim_state, Stage::Position);

    //=========================================================================
    // Compute the force and torque at the specified times.
    Manager manager(osimModel);
    manager.setIntegratorAccuracy(1e-6);
    osim_state.setTime(0.0);
    manager.initialize(osim_state);

    double final_t = 2.0;
    double nsteps = 10;
    double dt = final_t / nsteps;

    double I_y = ball.getInertia().getMoments()[1];

    double omega = sqrt(stiffness / I_y);

    for (int i = 1; i <= nsteps; ++i) {
        osim_state = manager.integrate(dt * i);
        osimModel.getMultibodySystem().realize(osim_state, Stage::Acceleration);

        // compute the current rotation about the y axis
        double simulated_theta = pin.getCoordinate(PinJoint::Coord::RotationZ)
                                         .getValue(osim_state);

        // compute the rotation about y-axis from the analytic solution
        // for 1-D spring-mass system with zero-velocity at initial offset.
        double analytical_theta =
                start_theta * cos(omega * osim_state.getTime());

        // check that the simulated solution is
        //  equivalent to the analytic solution
        ASSERT_EQUAL(analytical_theta, simulated_theta, 1e-4);

        // get the forces applied to the base and ball
        Array<double> model_forces = spring.getRecordValues(osim_state);

        // compute the expected force on the ball
        double analytical_moment = -stiffness * analytical_theta;

        // check analytical moment corresponds to the moment on the ball
        // in the Y direction, index = 4
        ASSERT_EQUAL(analytical_moment, model_forces[4], 2e-4);
    }

    manager.getStateStorage().print(
            "expression_based_bushing_rotational_model_states.sto");

    // Save the forces
    reporter->getForceStorage().print(
            "expression_based_bushing_rotational_model_forces.mot");

    // Now add damping to the bushing force and make sure energy is
    // monotonically decreasing
    spring.set_rotational_damping(Vec3(100.0));

    SimTK::State& osim_state2 = osimModel.initSystem();
    // set the initial pin joint angle
    pinCoord.setValue(osim_state2, start_theta);
    osimModel.getMultibodySystem().realize(osim_state2, Stage::Position);

    //==========================================================================
    // Compute the force and torque at the specified times.
    Manager manager2(osimModel);
    manager2.setIntegratorAccuracy(1e-6);
    osim_state2.setTime(0.0);
    manager2.initialize(osim_state2);

    double lastEnergy = 1E20; // Large
    for (int i = 1; i <= nsteps; ++i) {
        osim_state2 = manager2.integrate(dt * i);
        osimModel.getMultibodySystem().realize(
                osim_state2, Stage::Acceleration);
        double newEnergy = osimModel.calcKineticEnergy(osim_state2) +
                           osimModel.calcPotentialEnergy(osim_state2);
        ASSERT(newEnergy < lastEnergy);
        lastEnergy = newEnergy;
    }
    osimModel.disownAllComponents();

    // Before exiting lets see if copying the spring works
    ExpressionBasedBushingForce* copyOfSpring = spring.clone();

    ASSERT(*copyOfSpring == spring);
}

// Test our wrapping of elastic foundation in OpenSim
// Simple simulation of bouncing ball with dissipation should generate contact
// forces that settle to ball weight.
void testElasticFoundation() {
    using namespace SimTK;

    double start_h = 0.5;

    // Setup OpenSim model
    Model osimModel{"BouncingBallModelEF.osim"};

    // Create the force reporter
    ForceReporter* reporter = new ForceReporter(&osimModel);
    osimModel.addAnalysis(reporter);

    SimTK::State& osim_state = osimModel.initSystem();

    osimModel.getCoordinateSet().get("ball_ty").setValue(osim_state, start_h);
    osimModel.getMultibodySystem().realize(osim_state, Stage::Position);

    const OpenSim::Body& ball = osimModel.getBodySet().get("ball");

    //==========================================================================
    // Compute the force and torque at the specified times.
    Manager manager(osimModel);
    manager.setIntegratorAccuracy(1e-6);
    osim_state.setTime(0.0);
    manager.initialize(osim_state);

    double final_t = 2.0;

    // start timing
    clock_t startTime = clock();

    osim_state = manager.integrate(final_t);

    // end timing
    cout << "Elastic Foundation simulation time = "
         << 1.e3 * (clock() - startTime) / CLOCKS_PER_SEC << "ms" << endl;
    ;

    // make sure we can access dynamic variables
    osimModel.getMultibodySystem().realize(osim_state, Stage::Acceleration);

    // Print out the motion for visualizing/debugging
    manager.getStateStorage().print("bouncing_ball_states.sto");

    // Save the forces
    reporter->getForceStorage().print("elastic_contact_forces.mot");

    // Bouncing ball should have settled to rest on ground due to dissipation
    // In that case the force generated by contact should be identically body
    // weight in vertical and zero else where.
    OpenSim::ElasticFoundationForce& contact =
            (OpenSim::ElasticFoundationForce&)osimModel.getForceSet().get(
                    "contact");

    Array<double> contact_force = contact.getRecordValues(osim_state);
    ASSERT_EQUAL(
            contact_force[0], 0.0, 1e-4); // no horizontal force on the ball
    ASSERT_EQUAL(contact_force[1], -ball.getMass() * gravity_vec[1],
            2e-3); // vertical is weight
    ASSERT_EQUAL(
            contact_force[2], 0.0, 1e-4); // no horizontal force on the ball
    ASSERT_EQUAL(contact_force[3], 0.0, 1e-4); // no torque on the ball
    ASSERT_EQUAL(contact_force[4], 0.0, 1e-4); // no torque on the ball
    ASSERT_EQUAL(contact_force[5], 0.0, 1e-4); // no torque on the ball

    // Before exiting lets see if copying the spring works
    OpenSim::ElasticFoundationForce* copyOfForce = contact.clone();

    bool isEqual = (*copyOfForce == contact);

    if (!isEqual) {
        contact.print("originalForce.xml");
        copyOfForce->print("copyOfForce.xml");
    }

    ASSERT(isEqual);
}

// Test our wrapping of Hunt-Crossley force in OpenSim
// Simple simulation of bouncing ball with dissipation should generate contact
// forces that settle to ball weight.
void testHuntCrossleyForce() {
    using namespace SimTK;

    double start_h = 0.5;

    // Setup OpenSim model
    Model osimModel{"BouncingBall_HuntCrossley.osim"};

    // Create the force reporter
    ForceReporter* reporter = new ForceReporter(&osimModel);
    osimModel.addAnalysis(reporter);

    SimTK::State& osim_state = osimModel.initSystem();

    osimModel.getCoordinateSet()[4].setValue(osim_state, start_h);
    osimModel.getMultibodySystem().realize(osim_state, Stage::Position);

    const OpenSim::Body& ball = osimModel.getBodySet().get("ball");

    //==========================================================================
    // Compute the force and torque at the specified times.
    Manager manager(osimModel);
    manager.setIntegratorAccuracy(1e-6);
    osim_state.setTime(0.0);
    manager.initialize(osim_state);

    double final_t = 2.0;

    // start timing
    clock_t startTime = clock();

    osim_state = manager.integrate(final_t);

    // end timing
    cout << "Hunt Crossley simulation time = "
         << 1.e3 * (clock() - startTime) / CLOCKS_PER_SEC << "ms" << endl;

    // make sure we can access dynamic variables
    osimModel.getMultibodySystem().realize(osim_state, Stage::Acceleration);

    // Print out the motion for visualizing/debugging
    manager.getStateStorage().print("bouncing_ball_HC_states.sto");

    // Save the forces
    reporter->getForceStorage().print("HuntCrossley_contact_forces.mot");

    // Bouncing ball should have settled to rest on ground due to dissipation
    // In that case the force generated by contact should be identically body
    // weight in vertical and zero else where.
    OpenSim::HuntCrossleyForce& contact =
            (OpenSim::HuntCrossleyForce&)osimModel.getForceSet().get("contact");

    Array<double> contact_force = contact.getRecordValues(osim_state);
    ASSERT_EQUAL(
            contact_force[0], 0.0, 1e-4); // no horizontal force on the ball
    ASSERT_EQUAL(contact_force[1], -ball.getMass() * gravity_vec[1],
            1e-3); // vertical is weight
    ASSERT_EQUAL(
            contact_force[2], 0.0, 1e-4); // no horizontal force on the ball
    ASSERT_EQUAL(contact_force[3], 0.0, 1e-4); // no torque on the ball
    ASSERT_EQUAL(contact_force[4], 0.0, 1e-4); // no torque on the ball
    ASSERT_EQUAL(contact_force[5], 0.0, 1e-4); // no torque on the ball

    // Before exiting lets see if copying the force works
    OpenSim::HuntCrossleyForce* copyOfForce = contact.clone();

    bool isEqual = (*copyOfForce == contact);

    if(!isEqual){
        contact.print("originalForce.xml");
        copyOfForce->print("copyOfForce.xml");
    }

    ASSERT(isEqual);
}

// Test our wrapping of SimTK::SmoothSphereHalfSpaceForce.
// Simple simulation of bouncing ball with dissipation should generate contact
// forces that settle to ball weight.
void testSmoothSphereHalfSpaceForce()
{
    using namespace SimTK;

    double start_h = 0.5;

    // Setup OpenSim model
    Model osimModel{"BouncingBall_SmoothSphereHalfSpace.osim"};

    // Create the force reporter
    ForceReporter* reporter = new ForceReporter(&osimModel);
    osimModel.addAnalysis(reporter);

    SimTK::State& osim_state = osimModel.initSystem();

    osimModel.getCoordinateSet()[4].setValue(osim_state, start_h);
    osimModel.getMultibodySystem().realize(osim_state, Stage::Position );

    const OpenSim::Body &ball = osimModel.getBodySet().get("ball");

    //==========================================================================
    // Compute the force and torque at the specified times.
    Manager manager(osimModel);
    manager.setIntegratorAccuracy(1e-6);
    osim_state.setTime(0.0);
    manager.initialize(osim_state);

    double final_t = 2.0;

    // start timing
    clock_t startTime = clock();

    osim_state = manager.integrate(final_t);

    // end timing
    cout << "SmoothSphereHalfSpace simulation time = "
            << 1.e3*(clock()-startTime)/CLOCKS_PER_SEC << "ms" << endl;

    //make sure we can access dynamic variables
    osimModel.getMultibodySystem().realize(osim_state, Stage::Acceleration);

    // Print out the motion for visualizing/debugging
    manager.getStateStorage().print(
            "bouncing_ball_SmoothSphereHalfSpace_states.sto");

    // Save the forces
    reporter->getForceStorage().print("SmoothSphereHalfSpace_contact_forces.mot");

    // Bouncing ball should have settled to rest on ground due to dissipation
    // In that case the force generated by contact should be identically body weight
    // in vertical and zero else where.
    auto& contact = osimModel.getComponent<OpenSim::SmoothSphereHalfSpaceForce>(
            "forceset/contact");

    Array<double> contact_force = contact.getRecordValues(osim_state);
    ASSERT_EQUAL(contact_force[0], 0.0, 1e-4); // no horizontal force on the ball
    ASSERT_EQUAL(contact_force[1], -ball.getMass()*gravity_vec[1], 1e-3); // vertical is weight
    ASSERT_EQUAL(contact_force[2], 0.0, 1e-4); // no horizontal force on the ball
    ASSERT_EQUAL(contact_force[3], 0.0, 1e-4); // no torque on the ball
    ASSERT_EQUAL(contact_force[4], 0.0, 1e-4); // no torque on the ball
    ASSERT_EQUAL(contact_force[5], 0.0, 1e-4); // no torque on the ball

    // Before exiting lets see if copying the force works
    OpenSim::SmoothSphereHalfSpaceForce* copyOfForce = contact.clone();

    bool isEqual = (*copyOfForce == contact);

    if (!isEqual) {
        contact.print("originalForce.xml");
        copyOfForce->print("copyOfForce.xml");
    }

    ASSERT(isEqual);
}

void testCoordinateLimitForce() {
    using namespace SimTK;

    double mass = 1;
    double ball_radius = 0.25;

    // Setup OpenSim model
    auto osimModel = std::unique_ptr<Model>{new Model};
    osimModel->setName("CoordinateLimitForceTest");
    // OpenSim bodies
    const Ground& ground = osimModel->getGround();
    ;
    OpenSim::Body ball(
            "ball", mass, Vec3(0), mass * SimTK::Inertia::sphere(0.1));
    ball.attachGeometry(new Sphere{0.1});
    ball.scale(Vec3(ball_radius), false);

    // Add joints
    SliderJoint slider("slider", ground, Vec3(0), Vec3(0, 0, Pi / 2), ball,
            Vec3(0), Vec3(0, 0, Pi / 2));

    double positionRange[2] = {0.1, 2};
    // Rename coordinate for the slider joint
    auto& sliderCoord = slider.updCoordinate();
    sliderCoord.setName("ball_h");
    sliderCoord.setRange(positionRange);

    osimModel->addBody(&ball);
    osimModel->addJoint(&slider);

    osimModel->setGravity(gravity_vec);

    // Define the parameters of the Coordinate Limit Force
    double K_upper = 200.0;
    double K_lower = 1000.0;
    double damping = 0.01;
    double trans = 0.05;
    CoordinateLimitForce limitForce("ball_h", positionRange[1], K_upper,
            positionRange[0], K_lower, damping, trans, true);

    osimModel->addForce(&limitForce);

    // BAD: have to set memoryOwner to false or program will crash when this
    // test is complete.
    osimModel->disownAllComponents();

    osimModel->finalizeConnections();
    osimModel->print("CoordinateLimitForceTest.osim");

    // Check serialization and deserialization
    Model loadedModel{"CoordinateLimitForceTest.osim"};

    ASSERT(loadedModel == *osimModel, "Deserialized CoordinateLimitForceTest "
                                      "failed to be equivalent to original.");

    // check copy
    auto copyModel = std::unique_ptr<Model>{osimModel->clone()};

    ASSERT(*copyModel == loadedModel, "Clone of CoordinateLimitForceTest "
                                      "failed to be equivalent to original.");

    copyModel->print("cloneCoordinateLimitForceTest.osim");

    osimModel = std::move(copyModel);

    // Create the force reporter
    ForceReporter* reporter = new ForceReporter(osimModel.get());
    osimModel->addAnalysis(reporter);

    SimTK::State& osim_state = osimModel->initSystem();

    double dh = 0.2;
    double start_h = positionRange[1];
    double start_v = 2.0;
    const Coordinate& q_h = osimModel->getCoordinateSet()[0];
    q_h.setValue(osim_state, start_h);
    q_h.setSpeedValue(osim_state, start_v);

    osimModel->getMultibodySystem().realize(osim_state, Stage::Acceleration);

    CoordinateLimitForce* clf =
            dynamic_cast<CoordinateLimitForce*>(&osimModel->getForceSet()[0]);

    // initial energy of the system;
    double clfPE = clf->computePotentialEnergy(osim_state);
    double constStiffnessPE = 0.5 * K_upper * dh * dh;
    ASSERT(clfPE < constStiffnessPE);
    double energy0 = clfPE + mass * (-gravity_vec[1]) * start_h +
                     0.5 * mass * start_v * start_v;
    // system KE + PE including strain energy in CLF
    double eSys0 = osimModel->getMultibodySystem().calcEnergy(osim_state);

    //==========================================================================
    // Compute the force and torque at the specified times.
    Manager manager(*osimModel);
    manager.setIntegratorAccuracy(1e-6);
    osim_state.setTime(0.0);
    manager.initialize(osim_state);

    double final_t = 1.0;
    double nsteps = 20;
    double dt = final_t / nsteps;

    for (int i = 1; i <= nsteps; i++) {
        osim_state = manager.integrate(dt * i);
        osimModel->getMultibodySystem().realize(
                osim_state, Stage::Acceleration);

        double h = q_h.getValue(osim_state);
        double v = q_h.getSpeedValue(osim_state);

        // Now check that the force reported by spring
        Array<double> model_force = clf->getRecordValues(osim_state);

        double ediss = clf->getDissipatedEnergy(osim_state);
        double clfE = clf->computePotentialEnergy(osim_state) + ediss;

        // EK + EM of mass alone
        double eMass = 0.5 * mass * v * v - mass * gravity_vec[1] * h;
        double e = eMass + clfE;
        // system KE + PE including strain energy in CLF
        double eSys =
                osimModel->getMultibodySystem().calcEnergy(osim_state) + ediss;

        ASSERT_EQUAL(1.0, e / energy0, integ_accuracy,
                "CoordinateLimitForce Failed to conserve energy");
        ASSERT_EQUAL(1.0, eSys / eSys0, integ_accuracy,
                "CoordinateLimitForce Failed to conserve system energy");

        // get the forces applied to the ball by the limit force
        if (h > (positionRange[1] + trans)) {
            ASSERT_EQUAL(-K_upper * (h - positionRange[1]) - damping * v,
                    model_force[0], 1e-4);
        } else if (h < (positionRange[0] - trans)) {
            ASSERT_EQUAL(K_lower * (positionRange[0] - h) - damping * v,
                    model_force[0], 1e-4);
        } else if ((h < positionRange[1]) && (h > positionRange[0])) {
            // Verify no force is being applied when within limits
            ASSERT_EQUAL(0.0, model_force[0], 1e-5);
        }
    }

    manager.getStateStorage().print("coordinte_limit_force_model_states.sto");

    // Save the forces
    reporter->getForceStorage().print("limit_forces.mot");
}

void testCoordinateLimitForceRotational() {
    using namespace SimTK;

    double mass = 1;
    double edge = 0.2;

    // Setup OpenSim model
    Model osimModel{};
    osimModel.setName("RotationalCoordinateLimitForceTest");
    // OpenSim bodies
    const Ground& ground = osimModel.getGround();
    ;
    OpenSim::Body block("block", mass, Vec3(0),
            mass * SimTK::Inertia::brick(edge, edge, edge));
    block.attachGeometry(new Brick(Vec3(edge, edge, edge)));
    block.scale(Vec3(edge), false);

    // Add joints
    PinJoint pin("pin", ground, Vec3(0), Vec3(0, 0, 0), block,
            Vec3(0, -edge, 0), Vec3(0, 0, 0));

    // NOTE: Angular limits are in degrees NOT radians
    double positionRange[2] = {-30, 90};
    // Rename coordinate for the pin joint
    auto& pinCoord = pin.updCoordinate();
    pinCoord.setName("theta");
    pinCoord.setRange(positionRange);

    osimModel.addBody(&block);
    osimModel.addJoint(&pin);

    osimModel.setGravity(Vec3(0));

    // Define the parameters of the Coordinate Limit Force
    // For rotational coordinates, these are in Nm/degree
    double K_upper = 10.0;
    double K_lower = 20.0;
    double damping = 0.1;
    double trans = 0.05;
    CoordinateLimitForce limitForce("theta", positionRange[1], K_upper,
            positionRange[0], K_lower, damping, trans, true);

    osimModel.addForce(&limitForce);

    // BAD: have to set memoryOwner to false or program will crash when this
    // test is complete.
    osimModel.disownAllComponents();

    // Create the force reporter
    ForceReporter* reporter = new ForceReporter(&osimModel);
    osimModel.addAnalysis(reporter);

    SimTK::State& osim_state = osimModel.initSystem();

    // Start 2 degrees beyond the upper limit
    double start_q = SimTK_DEGREE_TO_RADIAN * positionRange[1] + SimTK::Pi / 90;
    double start_v = 0.0;
    const Coordinate& coord = osimModel.getCoordinateSet()[0];
    coord.setValue(osim_state, start_q);
    coord.setSpeedValue(osim_state, start_v);

    osimModel.getMultibodySystem().realize(osim_state, Stage::Acceleration);

    CoordinateLimitForce* clf =
            dynamic_cast<CoordinateLimitForce*>(&osimModel.getForceSet()[0]);

    // Now check that the force reported by spring
    Array<double> model_force = clf->getRecordValues(osim_state);

    ASSERT_EQUAL(model_force[0] / (-2 * K_upper), 1.0, integ_accuracy);

    double clfPE = clf->computePotentialEnergy(osim_state);
    double constSpringPE = 0.5 * (K_upper * 2.0) * 2.0 * SimTK_DEGREE_TO_RADIAN;
    ASSERT_EQUAL(clfPE / constSpringPE, 1.0, 0.001,
            "Specified upper rotational stiffness not met.");
    ASSERT(clfPE < constSpringPE);

    // Now test lower bound
    start_q = SimTK_DEGREE_TO_RADIAN * positionRange[0] - SimTK::Pi / 90;
    coord.setValue(osim_state, start_q);
    osimModel.getMultibodySystem().realize(osim_state, Stage::Acceleration);
    model_force = clf->getRecordValues(osim_state);

    ASSERT_EQUAL(model_force[0] / (2 * K_lower), 1.0, integ_accuracy);

    clfPE = clf->computePotentialEnergy(osim_state);
    constSpringPE = 0.5 * (K_lower * 2.0) * 2.0 * SimTK_DEGREE_TO_RADIAN;
    ASSERT_EQUAL(clfPE / constSpringPE, 1.0, 0.001);
    ASSERT(clfPE < constSpringPE);

    // total system energy prior to simulation
    double eSys0 = osimModel.getMultibodySystem().calcEnergy(osim_state);

    //==========================================================================
    // Perform a simulation an monitor energy conservation
    Manager manager(osimModel);
    manager.setIntegratorAccuracy(1e-8);
    osim_state.setTime(0.0);
    manager.initialize(osim_state);

    double final_t = 1.0;
    double nsteps = 20;
    double dt = final_t / nsteps;

    for (int i = 1; i <= nsteps; i++) {
        osim_state = manager.integrate(dt * i);
        osimModel.getMultibodySystem().realize(osim_state, Stage::Acceleration);

        double ediss = clf->getDissipatedEnergy(osim_state);
        // system KE + PE including strain energy in CLF
        double eSys =
                osimModel.getMultibodySystem().calcEnergy(osim_state) + ediss;
        /*double EKsys = */ osimModel.getMultibodySystem().calcKineticEnergy(
                osim_state);

        ASSERT_EQUAL(eSys / eSys0, 1.0, integ_accuracy);
    }

    manager.getStateStorage().print(
            "rotational_coordinte_limit_force_states.sto");

    // Save the forces
    reporter->getForceStorage().print("rotational_limit_forces.mot");
}

void testExternalForce() {
    using namespace SimTK;

    // define a new model properties
    double mass = 1;
    double angRange[2] = {-Pi, Pi};
    double posRange[2] = {-1, 1};

    // construct a new OpenSim model
    Model model;
    model.setName("ExternalForceTest");
    // OpenSim bodies
    const Ground& ground = model.getGround();
    OpenSim::Body tower("tower", mass, Vec3(0),
            mass * SimTK::Inertia::brick(0.1, 1.0, 0.2));
    tower.attachGeometry(new Brick(Vec3(0.1, 1.0, 0.2)));
    tower.scale(Vec3(0.1, 1.0, 0.2));

    // Add joint connecting the tower to the ground and associate joint to tower
    // body
    FreeJoint freeJoint("groundTower", ground, Vec3(0), Vec3(0), tower,
            Vec3(0, -0.5, 0), Vec3(0));
    // Set range and default value for each Coordinate in freeJoint.
    for (int i = 0; i < freeJoint.numCoordinates(); ++i) {
        if (freeJoint.get_coordinates(i).getMotionType() ==
                Coordinate::Translational) {
            freeJoint.upd_coordinates(i).setRange(posRange);
        } else {
            freeJoint.upd_coordinates(i).setRange(angRange);
        }
        freeJoint.upd_coordinates(i).setDefaultValue(0);
    }

    // add the tower body to the model
    model.addBody(&tower);
    model.addJoint(&freeJoint);

    // Force is 10N in Y, Torque is 2Nm about Z and Point is 0.1m in X
    Vec3 force(0, 10, 0), torque(0, 0, 2), point(0.1, 0, 0);
    Storage forces("external_force_data.sto");
    forces.setName("ExternalForcesData");

    /***************************** CASE 1 ************************************/
    // Apply force with both force and point specified in ground and no torque
    ExternalForce xf(forces, "force", "point", "", "tower", "ground", "ground");

    model.addForce(&xf);
    model.finalizeConnections(); // Needed so sockets have correct absolute path
                                 // on print
    model.print("ExternalForceTest.osim");
    ForceReporter frp;
    model.addAnalysis(&frp);
    // Everything allocated on the stack, so no need for model to own to free
    model.disownAllComponents();

    SimTK::State& s = model.initSystem();

    // set the starting location of the tower to be right over the point
    freeJoint.getCoordinate(FreeJoint::Coord::TranslationX).setValue(s, 0.1);

    double accuracy = 1e-6;
    Manager manager(model);
    manager.setIntegratorAccuracy(accuracy);

    // Specify the initial and final times of the simulation.
    double tf = 2.0;
    s.setTime(0.0);
    manager.initialize(s);
    s = manager.integrate(tf);

    manager.getStateStorage().print("external_force_test_model_states.sto");
    ;

    Vec3 a_y = force / mass + model.getGravity();
    double d_y = 0.5 * (a_y.norm()) * (tf * tf);

    double y_sim = model.getCoordinateSet()[4].getValue(s);

    // Vertical displacement
    ASSERT_EQUAL(d_y, y_sim, 10 * accuracy);
    // all rotations should remain zero
    for (int i = 0; i < 3; i++) {
        double val = model.getCoordinateSet()[i].getValue(s);
        ASSERT_EQUAL(0.0, val, 10 * accuracy);
    }
    ASSERT_EQUAL(
            point[0], model.getCoordinateSet()[3].getValue(s), 10 * accuracy);

    model.updForceSet().setSize(0);

    /***************************** CASE 2 ************************************/
    // Apply force with both force and point specified in ground as well as
    // torque
    ExternalForce xf2(
            forces, "force", "point", "torque", "tower", "ground", "ground");

    model.addForce(&xf2);
    model.print("ExternalForceTest.osim");
    // Everything allocated on the stack, so no need for model to own to free
    model.disownAllComponents();

    // recreate a new underlying system with corresponding state
    SimTK::State& s2 = model.initSystem();

    // set the starting location of the tower to be offset as to counter-balance
    // the torque point is 0.1, by moving fwd to 0.3, force has -0.2m moment-arm
    // to generate -2Nm
    freeJoint.getCoordinate(FreeJoint::Coord::TranslationX).setValue(s2, 0.3);
    model.setPropertiesFromState(s2);

    Manager manager2(model);
    manager2.setIntegratorAccuracy(accuracy);
    s2.setTime(0.0);
    manager2.initialize(s2);
    s2 = manager2.integrate(tf);

    // all dofs should remain constant
    for (int i = 0; i < model.getCoordinateSet().getSize(); i++) {
        double val = model.getCoordinateSet()[i].getValue(s2);
        double def = model.getCoordinateSet()[i].getDefaultValue();
        if (i == 4) { // Y-direction
            ASSERT_EQUAL(d_y, val, 10 * accuracy);
        } else {
            ASSERT_EQUAL(def, val, 10 * accuracy);
        }
    }

    model.updForceSet().setSize(0);
    /***************************** CASE 3 ************************************/
    // Apply force with only force (and torque) in ground but point on body
    ExternalForce xf3(
            forces, "force", "point", "torque", "tower", "ground", "tower");
    // Also Apply force with both force (and torque) and point in ground
    ExternalForce xf4(
            forces, "force", "point", "", "tower", "ground", "ground");

    model.addForce(&xf3);
    model.addForce(&xf4);

    // Everything allocated on the stack, so no need for model to own to free
    model.disownAllComponents();

    // recreate a new underlying system with corresponding state
    SimTK::State& s3 = model.initSystem();

    // only xf4 is should be affected and set it to offset Tz+px*Fy = 2+0.1*10
    // = 3.
    freeJoint.getCoordinate(FreeJoint::Coord::TranslationX)
            .setValue(s3, 0.4); // yield -3Nm for force only
    model.setPropertiesFromState(s3);

    Manager manager3(model);
    manager3.setIntegratorAccuracy(accuracy);
    s3.setTime(0.0);
    manager3.initialize(s3);
    s3 = manager3.integrate(tf);

    // all dofs should remain constant except Y
    for (int i = 0; i < model.getCoordinateSet().getSize(); i++) {
        double val = model.getCoordinateSet()[i].getValue(s3);
        double def = model.getCoordinateSet()[i].getDefaultValue();
        if (i != 4) { // ignore Y-direction
            ASSERT_EQUAL(def, val, 10 * accuracy);
        }
    }

    model.updForceSet().setSize(0);
    /***************************** CASE 4 ************************************/
    // Add joint connecting a "sensor" reference to the ground in which to
    // describe the applied external force
    OpenSim::Body sensor(
            "sensor", 1, Vec3(0), SimTK::Inertia::brick(0.1, 0.1, 0.1));
    sensor.attachGeometry(new Brick(Vec3(0.1, 0.1, 0.1)));
    sensor.scale(Vec3(0.02, 0.1, 0.01));

    // locate joint at 0.3m above tower COM
    WeldJoint weldJoint("sensorWeld", ground, Vec3(0, 0.8, 0), Vec3(0), sensor,
            Vec3(0), Vec3(0, 0, Pi / 2));

    // add the sensor body to the model
    model.addBody(&sensor);
    model.addJoint(&weldJoint);

    // Apply force with both force and point in sensor body
    ExternalForce xf5(
            forces, "force", "point", "", "tower", "sensor", "sensor");
    // Counter-balance with torque only in tower body
    ExternalForce xf6(forces, "", "", "torque", "tower", "tower", "");

    model.addForce(&xf5);
    model.addForce(&xf6);
    // Everything allocated on the stack, so no need for model to own to free
    model.disownAllComponents();

    // recreate a new underlying system with corresponding state
    SimTK::State& s4 = model.initSystem();
    model.getGravityForce().disable(s4);

    // only xf4 is should be affected and set it to offset Tz+px*Fy = 2+0.1*10
    // = 3.
    freeJoint.getCoordinate(FreeJoint::Coord::TranslationX).setValue(s4, 0);
    model.setPropertiesFromState(s4);

    RungeKuttaMersonIntegrator integrator4(model.getMultibodySystem());
    integrator4.setAccuracy(accuracy);
    Manager manager4(model);
    manager4.setIntegratorAccuracy(accuracy);
    s4.setTime(0.0);
    manager4.initialize(s4);
    s4 = manager4.integrate(tf);

    // all dofs should remain constant except X-translation
    for (int i = 0; i < model.getCoordinateSet().getSize(); i++) {
        double val = model.getCoordinateSet()[i].getValue(s4);
        double def = model.getCoordinateSet()[i].getDefaultValue();
        if (i != 3) ASSERT_EQUAL(def, val, 10 * accuracy);
    }
}

void testSerializeDeserialize() {
    std::cout << "Test serialize & deserialize." << std::endl;

    std::string origModelFile{"PushUpToesOnGroundWithMuscles.osim"};
    std::string oldModelFile{"testForces_SerializeDeserialize_old.osim"};
    std::string newModelFile{"testForces_SerializeDeserialize_new.osim"};

    // Toggle of 'isDisabled' property for some muscles in model file.
    std::set<std::string> flippedMuscles{
            "glut_med1_r", "bifemlh_r", "add_mag2_r"};
    {
        auto xml = SimTK::Xml::Document{origModelFile};
        auto thelenMuscles = xml.getRootElement()
                                     .getRequiredElement("Model")
                                     .getRequiredElement("ForceSet")
                                     .getRequiredElement("objects")
                                     .getAllElements("Thelen2003Muscle");
        for (unsigned i = 0; i < thelenMuscles.size(); ++i) {
            const auto& muscleName =
                    thelenMuscles[i].getRequiredAttributeValue("name");
            if (flippedMuscles.find(muscleName) != flippedMuscles.end()) {
                auto elem = thelenMuscles[i].getRequiredElement("isDisabled");
                elem.setValue("true");
            }
        }
        xml.writeToFile(oldModelFile);
    }

    // Model with Force::isDisabled (version < 30508)
    Model oldModel{oldModelFile};
    oldModel.print(newModelFile);
    Model newModel{newModelFile};

    const auto& oldForceSet = oldModel.getForceSet();
    const auto& newForceSet = newModel.getForceSet();

    ASSERT(oldForceSet.getSize() == newForceSet.getSize());
    for (int i = 0; i < oldForceSet.getSize(); ++i) {
        ASSERT(oldForceSet.get(i).get_appliesForce() ==
                newForceSet.get(i).get_appliesForce());

        if (flippedMuscles.find(newForceSet.get(i).getName()) !=
                flippedMuscles.end())
            ASSERT(newForceSet.get(i).get_appliesForce() == false);
    }

    std::remove(oldModelFile.c_str());
    std::remove(newModelFile.c_str());
}

void testTranslationalDampingEffect(Model& osimModel, Coordinate& sliderCoord,
        double start_h, Component& componentWithDamping) {
    using namespace SimTK;
    ASSERT(componentWithDamping.hasProperty("translational_damping"));

    AbstractProperty& aProp =
            componentWithDamping.updPropertyByName("translational_damping");
    Property<SimTK::Vec3>& aPropVec3 =
            dynamic_cast<Property<SimTK::Vec3>&>(aProp);
    aPropVec3.setValue(Vec3(100.));
    SimTK::State& osim_state2 = osimModel.initSystem();

    // set the initial height of the ball on slider
    sliderCoord.setValue(osim_state2, start_h);
    osimModel.getMultibodySystem().realize(osim_state2, Stage::Position);

    //==========================================================================
    // Compute the Energy to make sure it goes down due to damping
    Manager manager2(osimModel);
    manager2.setIntegratorAccuracy(1e-6);
    osim_state2.setTime(0.0);
    manager2.initialize(osim_state2);

    double lastEnergy = 1E20; // Large
    for (int i = 1; i <= 10; ++i) {
        osim_state2 = manager2.integrate(0.2 * i);
        osimModel.getMultibodySystem().realize(
                osim_state2, Stage::Acceleration);
        double newEnergy = osimModel.calcKineticEnergy(osim_state2) +
                           osimModel.calcPotentialEnergy(osim_state2);
        ASSERT(newEnergy < lastEnergy);
        lastEnergy = newEnergy;
    }
}

/*
=============
Test Setup 1:
=============
This test is modified from testPathSpring. The ligament is fixed to ground,
wraps over a cylinder (pulley) and suspends a hanging mass against gravity.
            ____
           / __ \
          / /  \ \
         /  \__/  \
        /          \
        |          |
        |       -------
        |       |     |
        |       |     |
     ___|___    -------
     ///////
Test 1.1: The force in the ligament must be equal to the inertial and
gravitational forces acting on the block.

Test 1.2: With damping coefficient set to 0, energy must be conserved.

Test 1.3: With damping on, but block velocity = 0, the damping force should be
        zero

Test 1.4: With damping on and block velocity set to 1.0, the damping force
        should match an analytical value.

============
Test Setup 2
============
One end of the ligament is fixed to ground the other to a block. The motion of
the block is prescribed to check the spring force, damping force and potential
energy in the ligament in the slack, toe, and linear regions and at multiple
speeds.
/|               ______
/|              |      |
/|--------------|      |
/|              |______|
/|
*/
void testBlankevoort1991Ligament() {
    using namespace SimTK;

    //=========================================================================
    // Test Setup 1
    //=========================================================================
    double mass = 0.1;
    double stiffness = 10;
    double restlength = 1.2;
    double start_h = 0.1;

    // Setup OpenSim model
    Model osimModel{};
    osimModel.setName("LigamentTest");
    osimModel.setGravity(gravity_vec);

    // OpenSim bodies
    Ground& ground = osimModel.updGround();

    WrapCylinder* pulley1 = new WrapCylinder();
    pulley1->setName("pulley1");
    pulley1->set_radius(0.1);
    pulley1->set_length(0.05);
    pulley1->set_translation(Vec3(0.1, 0.4, 0));

    WrapCylinder* pulley2 = new WrapCylinder();
    pulley2->setName("pulley2");
    pulley2->set_radius(0.1);
    pulley2->set_length(0.05);
    pulley2->set_translation(Vec3(0.1, 0.4, 0));

    // Add the wrap object to the body, which takes ownership of it
    ground.addWrapObject(pulley1);
    ground.addWrapObject(pulley2);

    // Add Block to Model
    OpenSim::Body* block = new OpenSim::Body("block", mass, Vec3(0),
            mass * SimTK::Inertia::brick(0.05, 0.05, 0.05));
    block->attachGeometry(new Brick(Vec3(0.05, 0.05, 0.05)));

    SliderJoint* slider = new SliderJoint("slider", ground, Vec3(0.2, 0, 0),
            Vec3(0, 0, Pi / 2), *block, Vec3(0), Vec3(0, 0, Pi / 2));

    double positionRange[2] = {-10, 10};
    auto& sliderCoord = slider->updCoordinate();
    sliderCoord.setName("block_h");
    sliderCoord.setRange(positionRange);

    osimModel.addJoint(slider);
    osimModel.addBody(block);

    Blankevoort1991Ligament* ligament =
            new Blankevoort1991Ligament("ligament", stiffness, restlength);
    ligament->set_damping_coefficient(0.0);

    ligament->updGeometryPath().appendNewPathPoint(
            "origin", ground, Vec3(0.0, 0.0, 0.0));
    ligament->updGeometryPath().addPathWrap(*pulley1);
    ligament->updGeometryPath().appendNewPathPoint(
            "midpoint", ground, Vec3(0.1, 0.6, 0.0));
    ligament->updGeometryPath().addPathWrap(*pulley2);
    ligament->updGeometryPath().appendNewPathPoint(
            "insertion", *block, Vec3(0.0, 0.0, 0.0));

    osimModel.addForce(ligament);

    // Create the force reporter
    ForceReporter reporter = ForceReporter(&osimModel);
    osimModel.addAnalysis(&reporter);

    //osimModel.setUseVisualizer(true);
    SimTK::State& osim_state = osimModel.initSystem();

    sliderCoord.setValue(osim_state, start_h);
    osimModel.getMultibodySystem().realize(osim_state, Stage::Position);

    // Compute system energy at initial state
    osimModel.realizeReport(osim_state);
    double KE0 = osimModel.calcKineticEnergy(osim_state);
    double PE0 = osimModel.calcPotentialEnergy(osim_state);
    double E0 = KE0 + PE0;

    // Simulate
    Manager manager(osimModel);
    manager.setIntegratorAccuracy(1e-6);
    osim_state.setTime(0.0);
    manager.initialize(osim_state);

    double final_t = 3.0;
    osim_state = manager.integrate(final_t);

    // tension should only be velocity dependent
    osimModel.getMultibodySystem().realize(osim_state, Stage::Velocity);

    // Force in ligament
    double model_force = ligament->getTotalForce(osim_state);

    // get acceleration of the block
    osimModel.getMultibodySystem().realize(osim_state, Stage::Acceleration);
    double hddot =
            osimModel.getCoordinateSet().get("block_h").getAccelerationValue(
                    osim_state);

    // the ligament tension should be the weight of the block
    double analytical_force = -(gravity_vec(1) - hddot) * mass;

    // Save the forces
    reporter.getForceStorage().print("ligament_forces.mot");

    // something is wrong if the block does not reach equilibrium
    ASSERT_EQUAL(analytical_force, model_force, 1e-3, __FILE__, __LINE__,
        "Expected Blankevoort1991Ligament to force to be equal to the "
        "inertial and gravitational forces acting on block as is "
        "necessary for dynamic equilibrium.");

    // Check that Energy is conserved
    double KE1 = osimModel.calcKineticEnergy(osim_state);
    double PE1 = osimModel.calcPotentialEnergy(osim_state);
    double E1 = KE1 + PE1;

    ASSERT_EQUAL(E0, E1, 1e-3, __FILE__, __LINE__,
        "Expected Blankevoort1991Ligament with damping set to zero to "
        "conserve energy in a forward dynamic simulation.");

    // Test damping force
    double damping_coeff = 0.001;
    ligament->set_damping_coefficient(damping_coeff);

    osim_state = osimModel.initSystem();
    osimModel.realizeVelocity(osim_state);
    double damp_force0 = ligament->getDampingForce(osim_state);

    ASSERT_EQUAL(damp_force0, 0.0000, 1e-3, __FILE__, __LINE__,
        "Expected Blankevoort1991Ligament damping force to be zero when all "
        "generalized speeds were zero.");

    double block_velocity = -1.0;
    sliderCoord.setSpeedValue(osim_state, block_velocity);

    osimModel.realizeReport(osim_state);
    double damp_force1 = ligament->getDampingForce(osim_state);

    double analytical_damping_force = -damping_coeff * block_velocity;

    ASSERT_EQUAL(damp_force1, analytical_damping_force, 1e-3, __FILE__,
        __LINE__,
        "Expected Blankevoort1991Ligament damping force to be equal to "
        "analytical value.");

    //=========================================================================
    // Test Setup 2
    //=========================================================================
    Model model;

    OpenSim::Body* brick = new OpenSim::Body(
            "brick", 1.0, Vec3(0), SimTK::Inertia::brick(0.05, 0.05, 0.05));
    block->attachGeometry(new Brick(Vec3(0.05, 0.05, 0.05)));
    model.addBody(brick);

    SliderJoint* slot = new SliderJoint("slot", model.updGround(), Vec3(0),
            Vec3(0), *brick, Vec3(0), Vec3(0));

    double range[2] = {-10, 10};
    auto& slotCoord = slot->updCoordinate();
    slotCoord.setName("displacement");
    slotCoord.setRange(range);
    model.addJoint(slot);

    double lig_stiffness = 10;
    double lig_slack_length = 1.0;

    Blankevoort1991Ligament* lig = new Blankevoort1991Ligament("ligament",
            model.updGround(), SimTK::Vec3(0), *brick, SimTK::Vec3(0),
            lig_stiffness, lig_slack_length);
    model.addForce(lig);

    SimTK::State state = model.initSystem();

    // Strech the ligament
    int nSteps = 25;
    double disp = 0.95;
    double time = 0.0;
    double disp_step = 0.01;
    double time_step = 1.0;

    std::vector<double> ind_col;

    std::vector<std::string> outputs;
    outputs.push_back("strain");
    outputs.push_back("strain_rate");
    outputs.push_back("length");
    outputs.push_back("lengthening_speed");
    outputs.push_back("spring_force");
    outputs.push_back("damping_force");
    outputs.push_back("total_force");
    outputs.push_back("potential_energy");

    SimTK::Matrix output_data(nSteps, (int)outputs.size());

    for (int i = 0; i < nSteps; ++i) {
        slotCoord.setValue(state, disp);
        slotCoord.setSpeedValue(state, disp_step / time_step);
        state.setTime(time);
        model.realizeReport(state);

        ind_col.push_back(time);
        for (int j = 0; j < (int)outputs.size(); ++j) {
            output_data(i, j) = lig->getOutputValue<double>(state, outputs[j]);
        }

        disp += disp_step;
        time += time_step;
    }
    TimeSeriesTable results(ind_col,output_data,outputs);
    STOFileAdapter::write(results, "ligament_strain_test.sto");

    //Check that potential energy and spring and damping forces are zero
    //when the ligament is slack
    ASSERT(results.getDependentColumn("strain").getElt(0, 0) < 0.0,
            __FILE__, __LINE__,
        "Expected Blankevoort1991Ligament to be slack at first time step of "
        "test case.");

    ASSERT_EQUAL(results.getDependentColumn("potential_energy").getElt(0, 0),
            0.0, 1e-3,
        __FILE__, __LINE__,
        "Expected potential energy in Blankevoort1991Ligament to be "
        "equal to zero when the ligament is slack");

    ASSERT_EQUAL(results.getDependentColumn("spring_force").getElt(0, 0),
            0.0, 1e-3,
        __FILE__, __LINE__,
        "Expected spring_force in Blankevoort1991Ligament to be"
        "equal to zero when the ligament is slack");

    ASSERT_EQUAL(results.getDependentColumn("damping_force").getElt(0, 0),
            0.0, 1e-3,
        __FILE__, __LINE__,
        "Expected damping_force in Blankevoort1991Ligament to be"
        "equal to zero when the ligament is slack");

    double damping_frc = results.getDependentColumn("damping_force").getElt(12, 0);
    double spring_frc = results.getDependentColumn("spring_force").getElt(12, 0);

    ASSERT_EQUAL(results.getDependentColumn("total_force").getElt(12, 0),
            damping_frc + spring_frc, 1e-3,
        __FILE__, __LINE__,
        "Expected total_force in Blankevoort1991Ligament to be"
        "equal to the damping_force + spring_force");

    //Check that the spring_force and potential_energy are greater when the
    //ligment crosses the transition from the toe region to linear region

    int toe_index = 10;
    int linear_index = 12;

    double transition_strain = lig->get_transition_strain();

    ASSERT(results.getDependentColumn("strain").getElt(toe_index, 0) <
        transition_strain, __FILE__, __LINE__,
        "Expected strain at the toe_index to be less than the "
        "transition_strain property in Blankevoort1991Ligament test.");

    ASSERT(results.getDependentColumn("strain").getElt(linear_index, 0) >
        transition_strain, __FILE__, __LINE__,
        "Expected strain at the linear_index to be greater than the "
        "transition_strain property in Blankevoort1991Ligament test.");

    ASSERT(results.getDependentColumn("potential_energy")
                   .getElt(linear_index, 0) >
                   results.getDependentColumn("potential_energy")
                           .getElt(toe_index, 0),
        __FILE__, __LINE__,
        "Expexted potential_energy in the Blankevoort1991Ligament to be "
        "greater in the linear region compared to the toe region");

    ASSERT(results.getDependentColumn("spring_force").getElt(linear_index, 0) >
        results.getDependentColumn("spring_force").getElt(toe_index, 0),
        __FILE__, __LINE__,
        "Expected the spring_force in the Blankevoort1991Ligament to be "
        " greater in the linear region compared to the toe region");

    //Check that damping is nonzero if ligament is lengthening
    slotCoord.setSpeedValue(state, 1.0);
    model.realizeReport(state);
    double damping_lengthening =
        lig->getOutputValue<double>(state, "damping_force");

    ASSERT(damping_lengthening > 0.0, __FILE__, __LINE__,
        "Expected the damping force in Blankevoort1991Ligament to be greater "
        "than zero when the ligament is streched beyond the slack length "
        "and the lengthening_speed is positive.");

    //Check that damping is zero if ligament is shortening
    slotCoord.setSpeedValue(state, -1.0);
    model.realizeReport(state);
    double damping_shortening =
        lig->getOutputValue<double>(state, "damping_force");

    ASSERT_EQUAL(damping_shortening, 0.0, 1e-3, __FILE__, __LINE__,
        "Expected the damping force in Blankevoort1991Ligament to be "
        " zero when the ligament is streched beyond the slack "
        "length, but the lengthening_speed is negative.");

    //Check linear stiffness in force/length
    slotCoord.setSpeedValue(state, 0.0);
    slotCoord.setValue(state, 1.06);
    model.realizePosition(state);
    double F1 = lig->getOutputValue<double>(state, "spring_force");
    double L1 = lig->getOutputValue<double>(state, "length");

    slotCoord.setValue(state, 1.08);
    model.realizePosition(state);
    double F2 = lig->getOutputValue<double>(state, "spring_force");
    double L2 = lig->getOutputValue<double>(state, "length");

    double calc_stiff = (F2 - F1) / (L2 - L1);
    double lig_stiff = lig->getLinearStiffnessForcePerLength();

    ASSERT_EQUAL(calc_stiff, lig_stiff, 1e-3, __FILE__, __LINE__,
        "Expected the calculated linear_stiffness in force/length in the "
        "Blankevoort1991Ligament to be equal to the value returned by "
        "getLinearStiffnessForcePerLength().");

    //Check setting slack_length through reference force
    double ref_force = 2.5;
    double coord_value = 1.1;
    sliderCoord.setValue(state, coord_value);
    lig->setSlackLengthFromReferenceForce(ref_force, state);
    state = model.initSystem();
    sliderCoord.setValue(state, coord_value);
    model.realizePosition(state);
    double reported_force = lig->getSpringForce(state);

    ASSERT_EQUAL(ref_force, reported_force, 1e-3, __FILE__, __LINE__,
        "Expected the force in the Blankevoort1991Ligament at the input "
        "reference state be equal to the force value input "
        "to setSlackLengthFromReferenceForce().");

    //Check setting slack_length through reference strain
    double ref_strain = 0.05;
    sliderCoord.setValue(state, coord_value);
    lig->setSlackLengthFromReferenceStrain(ref_strain, state);
    state = model.initSystem();
    sliderCoord.setValue(state, coord_value);
    model.realizePosition(state);
    double reported_strain = lig->getStrain(state);

    ASSERT_EQUAL(ref_strain, reported_strain, 1e-3, __FILE__, __LINE__,
        "Expected the strain in the Blankevoort1991Ligament at the input "
        "reference state be equal to the strain value input "
        "to setSlackLengthFromReferenceStrain().");
}

/*
* testSmith2018ArticularContactForce
*

Test 1
Half sphere is collided with plane

 o                  o   
  o                o
    o            o 
       o  o  o
--------------------

Test 2 
Cylinder pin (flat circular end head) is collided with plane

    |        |
    |        |
    |        |
    |        |
    |________|
-------------------
T2.1: cylinder vertically penetrates plane
T2.2: cylinder slides horizontally
T2.3: cylinder twists (rotation around vertical axis)


*/
void testSmith2018ArticularContactForce() {
    Model model = Model();
    
    double gravity_acc = 9.8067;
    model.setGravity(SimTK::Vec3(0, -gravity_acc, 0));

    double mass = 1.0;
    double radius = 0.1;
    
    OpenSim::Body* plane = new OpenSim::Body("plane", 0.0, SimTK::Vec3(0),
        SimTK::Inertia(0));
        //mass * SimTK::Inertia::brick(0.05, 0.05, 0.05));

    OpenSim::Body* indenter = new OpenSim::Body("indenter", mass,
        SimTK::Vec3(0), mass * SimTK::Inertia::brick(0.05, 0.05, 0.05));

    //indenter->attachGeometry(new Mesh("half_sphere_10cm_radius.stl"));

    model.addBody(plane);
    model.addBody(indenter);

    SpatialTransform indenter_transform = SpatialTransform();
    indenter_transform.upd_rotation1().set_coordinates(0, "plane_indenter_rx");
    indenter_transform.upd_rotation1().set_axis(SimTK::Vec3(1, 0, 0));
    indenter_transform.upd_rotation1().set_function(LinearFunction());

    indenter_transform.upd_rotation2().set_coordinates(0, "plane_indenter_ry");
    indenter_transform.upd_rotation2().set_axis(SimTK::Vec3(0, 1, 0));
    indenter_transform.upd_rotation2().set_function(LinearFunction());

    indenter_transform.upd_rotation3().set_coordinates(0, "plane_indenter_rz");
    indenter_transform.upd_rotation3().set_axis(SimTK::Vec3(0, 0, 1));
    indenter_transform.upd_rotation3().set_function(LinearFunction());

    indenter_transform.upd_translation1().set_coordinates(0, "plane_indenter_tx");
    indenter_transform.upd_translation1().set_axis(SimTK::Vec3(1, 0, 0));
    indenter_transform.upd_translation1().set_function(LinearFunction());

    indenter_transform.upd_translation2().set_coordinates(0, "plane_indenter_ty");
    indenter_transform.upd_translation2().set_axis(SimTK::Vec3(0, 1, 0));
    indenter_transform.upd_translation2().set_function(LinearFunction());

    indenter_transform.upd_translation3().set_coordinates(0, "plane_indenter_tz");
    indenter_transform.upd_translation3().set_axis(SimTK::Vec3(0, 0, 1));
    indenter_transform.upd_translation3().set_function(LinearFunction());

    WeldJoint* ground_plane = new WeldJoint("ground_plane", model.getGround(), *plane);
    CustomJoint* plane_indenter = new CustomJoint("plane_indenter", *plane, *indenter, indenter_transform);

    model.addJoint(ground_plane);
    model.addJoint(plane_indenter);

    Smith2018ContactMesh* plane_mesh = new Smith2018ContactMesh("plane", "x_z_plane.stl", *plane);
    Smith2018ContactMesh* indenter_mesh = new Smith2018ContactMesh("indenter", "half_sphere_10cm_radius.stl", *indenter);
    
    //indenter_mesh->set_mesh_back_file("x_z_disk_10cm_radius.stl");
    indenter_mesh->set_mesh_back_file("half_sphere_10cm_radius.stl");
    indenter_mesh->set_use_variable_thickness(true);
    indenter_mesh->set_max_thickness(0.5);

    double E_indenter = 1e6;
    double v_indenter = 0.46;
    double E_plane = 1e6;
    double v_plane = 0.46;

    plane_mesh->set_elastic_modulus(E_plane);
    plane_mesh->set_poissons_ratio(v_plane);
    plane_mesh->set_thickness(0.2);

    indenter_mesh->set_elastic_modulus(E_indenter);
    indenter_mesh->set_poissons_ratio(v_indenter);
    indenter_mesh->set_thickness(0.1);

    model.addContactGeometry(plane_mesh);
    model.addContactGeometry(indenter_mesh);

    Smith2018ArticularContactForce* contact = new Smith2018ArticularContactForce("contact", *plane_mesh, *indenter_mesh);
    contact->set_max_proximity(0.2);
    contact->set_elastic_foundation_formulation("linear");
    contact->set_use_lumped_contact_model(true);
    model.addForce(contact);

    model.setUseVisualizer(true);
    SimTK::State& state = model.initSystem();
    contact->setModelingOption(state, "flip_meshes", 1);

    model.print("contact_test.osim");

    Coordinate& tx = model.updCoordinateSet().get("plane_indenter_tx");
    Coordinate& ty = model.updCoordinateSet().get("plane_indenter_ty");
    Coordinate& tz = model.updCoordinateSet().get("plane_indenter_tz");
    Coordinate& rx = model.updCoordinateSet().get("plane_indenter_rx");
    Coordinate& ry = model.updCoordinateSet().get("plane_indenter_ry");
    Coordinate& rz = model.updCoordinateSet().get("plane_indenter_rz");

    model.realizeReport(state);

    const ModelVisualizer& viz = model.getVisualizer();
    viz.show(state);

    // Verify Variable Thickness
    double sphere_radius = 0.1;
    double thickness = 0;
    double max_thickness = 0;
    for (int i = 0; i < indenter_mesh->getNumMeshBackVertices(); ++i) {
        thickness = indenter_mesh->getMeshBackVertexThickness(i);
        //std::cout << thickness << std::endl;
        if (thickness > max_thickness) {
            max_thickness = thickness;
        }
    }
    //std::cout << max_thickness << " " << sphere_radius << std::endl;
/*    ASSERT_EQUAL(sphere_radius, max_thickness, 1e-3, __FILE__,
        __LINE__,
        "Expected the maximum thickness to be equal to the radius.");
        */
   //indenter_mesh->set_use_variable_thickness(false);

    // Verify Contact Area
    double expected_casting_area = 2 * SimTK::Pi * SimTK::square(radius);
    double reported_casting_area = contact->getOutputValue<double>(
            state, "casting_total_contact_area");

     ASSERT_EQUAL(expected_casting_area, reported_casting_area, 1e-3,
     __FILE__, __LINE__,
            "Expected the contact area of the half sphere to be equal to the "
            "analytical formula.");

    double expected_target_area = SimTK::Pi * SimTK::square(radius);
    double reported_target_area =
            contact->getOutputValue<double>(state, "target_total_contact_area");

     ASSERT_EQUAL(expected_target_area, reported_target_area, 1e-3, __FILE__,
     __LINE__,
            "Expected the contact area of the half sphere to be equal to the "
            "analytical formula.");
       
    //Move half sphere to half depth
     double depth = 0.05;
    ty.setValue(state, depth);
    model.realizeReport(state);

    // Max Proximity
    double reported_casting_max_prox = contact->getOutputValue<double>(
            state, "casting_total_max_proximity");

    double reported_target_max_prox = contact->getOutputValue<double>(
            state, "target_total_max_proximity");

    double expected_max_prox = radius/2;

    ASSERT_EQUAL(expected_max_prox, reported_casting_max_prox, 1e-3, __FILE__,
            __LINE__,
            "Expected the maximum casting proximity "
            "to be equal to the half sphere radius.");

    ASSERT_EQUAL(expected_max_prox, reported_target_max_prox, 1e-3, __FILE__,
            __LINE__,
            "Expected the maximum target proximity "
            "to be equal to the half sphere radius.");
    
    //Test 2
    //===================================================================
    indenter_mesh->set_mesh_file("x_z_disk_10cm_radius.stl");    
    //indenter_mesh->set_mesh_file("half_sphere_10cm_radius.stl");
    model.finalizeFromProperties();

    state = model.initSystem();
    contact->setModelingOption(state, "flip_meshes", 1);
    model.realizeReport(state);
    
    depth = 0.01;
    ty.setValue(state, (0.0 - depth));
    model.realizeReport(state);


    //Proximity

    double reported_max_prox = contact->getOutputValue<double>(
            state, "casting_total_max_proximity");

    double reported_mean_prox = contact->getOutputValue<double>(
            state, "casting_total_mean_proximity");

    ASSERT_EQUAL(depth, reported_max_prox, 1e-3, __FILE__,
            __LINE__,
            "Expected the max proximity to be equal to indentation depth "
            "of flat cylinder bottom.");
    
    ASSERT_EQUAL(reported_mean_prox, reported_max_prox, 1e-3, __FILE__,
            __LINE__,
            "Expected the mean and max proximity to be equal because "
            "the cylinder bottom is flat.");

    // Pressure
    double reported_casting_total_max_pressure = contact->getOutputValue<double>(
            state, "casting_total_max_pressure");

    double reported_mean_pressure = contact->getOutputValue<double>(
            state, "casting_total_mean_pressure");

    double area = contact->getOutputValue<double>(
        state, "casting_total_contact_area");

    SimTK::Vec3 reported_force = contact->getOutputValue<SimTK::Vec3>(
            state, "casting_total_contact_force");

    double expected_pressure = reported_force.norm() / area;

    ASSERT_EQUAL(expected_pressure, reported_mean_pressure, 1e-3, __FILE__, __LINE__,
            "Expected the mean pressure to be equal to force/area "
            "for flat disk.");

    ASSERT_EQUAL(reported_mean_pressure, reported_casting_total_max_pressure, 1e-3, __FILE__,
            __LINE__,
            "Expected the mean and max pressure to be equal because "
            "the disk is flat.");

    // Force
    SimTK::Vec3 reported_moment = contact->getOutputValue<SimTK::Vec3>(
        state, "casting_total_contact_moment");

    SimTK::Vec3 reported_COprox = contact->getOutputValue<SimTK::Vec3>(
        state, "casting_total_center_of_proximity");

    SimTK::Vec3 reported_COpressure = contact->getOutputValue<SimTK::Vec3>(
        state, "casting_total_center_of_pressure");

    ASSERT_EQUAL(0., reported_force(0), 1e-3, __FILE__, __LINE__,
        "Expected the x (horizontal) contact force component to be "
        "zero for flat disk contacting plane.");

    ASSERT_EQUAL(0., reported_force(2), 1e-3, __FILE__, __LINE__,
        "Expected the z (horizontal) contact force component to be "
        "zero for flat disk contacting plane.");

    ASSERT_EQUAL(0., reported_moment.norm(), 1e-3, __FILE__, __LINE__,
        "Expected the contact moment to be "
        "zero for flat disk contacting plane.");

    ASSERT_EQUAL(reported_COprox, reported_COpressure, 1e-3, __FILE__, __LINE__,
        "Expected the center of proximity and pressure to be the same "
        "flat disk contacting plane.");

    // Test Regional Metrics

    SimTK::Vec3 reported_target_force = contact->getOutputValue<SimTK::Vec3>(
        state, "target_total_contact_force");

    SimTK::Vec3 reported_target_moment = contact->getOutputValue<SimTK::Vec3>(
        state, "target_total_contact_moment");

    SimTK::Vec3 reported_target_COprox = contact->getOutputValue<SimTK::Vec3>(
        state, "target_total_center_of_proximity");

    SimTK::Vec3 reported_target_COpressure = contact->getOutputValue<SimTK::Vec3>(
        state, "target_total_center_of_pressure");

    reported_target_area = contact->getOutputValue<double>(
        state, "target_total_contact_area");

    SimTK::Vector_<SimTK::Vec3> reported_target_regional_force = contact->getOutputValue<SimTK::Vector_<SimTK::Vec3>>(
        state, "target_regional_contact_force");

    SimTK::Vector_<SimTK::Vec3> reported_target_regional_moment = contact->getOutputValue<SimTK::Vector_<SimTK::Vec3>>(
        state, "target_regional_contact_moment");

    SimTK::Vector_<SimTK::Vec3> reported_target_regional_COprox = contact->getOutputValue<SimTK::Vector_<SimTK::Vec3>>(
        state, "target_regional_center_of_proximity");

    SimTK::Vector_<SimTK::Vec3> reported_target_regional_COpressure = contact->getOutputValue<SimTK::Vector_<SimTK::Vec3>>(
        state, "target_regional_center_of_pressure");

    SimTK::Vector reported_target_regional_area = contact->getOutputValue<SimTK::Vector>(
        state, "target_regional_contact_area");

    //regional area
    ASSERT_EQUAL(
        reported_target_regional_area(0) + reported_target_regional_area(1),
        reported_target_area, 1e-3, __FILE__, __LINE__,
        "Expected the regional contact area in the +x and -x regions (0,1) "
        "to sum to the total contact area.");

    ASSERT_EQUAL(
        reported_target_regional_area(4) + reported_target_regional_area(5),
        reported_target_area, 1e-3, __FILE__, __LINE__,
        "Expected the regional contact area in the +z and -z regions (4,5) "
        "to sum to the total contact area.");

    ASSERT_EQUAL(
        reported_target_regional_force(0) + reported_target_regional_force(1),
        reported_target_force, 1e-3, __FILE__, __LINE__,
        "Expected the regional contact force in the +x and -x regions (0,1) "
        "to sum to the total contact force.");

    ASSERT_EQUAL(
        reported_target_regional_force(4) + reported_target_regional_force(5),
        reported_target_force, 1e-3, __FILE__, __LINE__,
        "Expected the regional contact force in the +z and -z regions (4,5) "
        "to sum to the total contact force.");

    ASSERT_EQUAL(
        reported_target_regional_moment(0) + reported_target_regional_moment(1),
        reported_target_moment, 1e-3, __FILE__, __LINE__,
        "Expected the regional contact moment in the +x and -x regions (0,1) "
        "to sum to the total contact moment.");

    ASSERT_EQUAL(
        reported_target_regional_moment(4) + reported_target_regional_moment(5),
        reported_target_moment, 1e-3, __FILE__, __LINE__,
        "Expected the regional contact moment in the +z and -z regions (4,5) "
        "to sum to the total contact moment.");
    
    // Translate +x +z to test different regions
    SimTK::Vector tx_values = SimTK::Vector(2, -1);
    tx_values.set(0, -0.2);
    tx_values.set(1, 0.2);

    SimTK::Vector tz_values = SimTK::Vector(2, -1);
    tz_values.set(0, -0.2);
    tz_values.set(1, 0.2);


    for (int i = 0; i < tx_values.size(); ++i) {
        for (int j = 0; j < tz_values.size(); ++j) {
            tx.setValue(state, tx_values(i));
            tz.setValue(state, tz_values(j));
            model.realizeReport(state);

            reported_target_force = contact->getOutputValue<SimTK::Vec3>(
                state, "target_total_contact_force");

            reported_target_moment = contact->getOutputValue<SimTK::Vec3>(
                state, "target_total_contact_moment");

            reported_target_COprox = contact->getOutputValue<SimTK::Vec3>(
                state, "target_total_center_of_proximity");

            reported_target_COpressure = contact->getOutputValue<SimTK::Vec3>(
                state, "target_total_center_of_pressure");

            ASSERT_EQUAL(tx_values(i), reported_target_COpressure(0), 
                1e-3, __FILE__, __LINE__,
                "Expected the plane COP x component to be located at the "
                "center ofthe flat disk.");

            ASSERT_EQUAL(tz_values(j), reported_target_COpressure(2), 
                1e-3, __FILE__, __LINE__,
                "Expected the plane COP z component to be located at the "
                "center of the flat disk.");

            SimTK::Vec3 expected_moment = reported_target_COpressure % reported_target_force;

            ASSERT_EQUAL(expected_moment, reported_target_moment, 1e-3, __FILE__, __LINE__,
                "Expected the contact moment to be equal to COP x Force .");

            // Regional Metrics

            reported_target_regional_force = contact->getOutputValue<SimTK::Vector_<SimTK::Vec3>>(
                state, "target_regional_contact_force");

            reported_target_regional_moment = contact->getOutputValue<SimTK::Vector_<SimTK::Vec3>>(
                state, "target_regional_contact_moment");

            reported_target_regional_COprox = contact->getOutputValue<SimTK::Vector_<SimTK::Vec3>>(
                state, "target_regional_center_of_proximity");

            reported_target_regional_COpressure = contact->getOutputValue<SimTK::Vector_<SimTK::Vec3>>(
                state, "target_regional_center_of_pressure");

            ASSERT_EQUAL(reported_target_regional_COprox(0),
                reported_target_regional_COpressure(0), 1e-3, __FILE__, __LINE__,
                "Expected the regional center of pressure and proximity "
                "to be equal for flat disk.");

            if (tx_values(i) < 0.0) {
                ASSERT_EQUAL(reported_target_regional_COpressure(0)(0), tx_values(i), 1e-3, __FILE__, __LINE__,
                    "Expected the region 0 plane COP x component to be "
                    "located at the center of the flat disk when contact is "
                    "in -x half space");
                ASSERT_EQUAL(reported_target_regional_COpressure(1), SimTK::Vec3(-1), 1e-3, __FILE__, __LINE__,
                    "Expected the region 1 plane COP to be -1 when contact is "
                    "in -x half space.");
                ASSERT_EQUAL(reported_target_regional_force(0), reported_target_force, 1e-3, __FILE__, __LINE__,
                    "Expected the region 0 plane force to be equal to the "
                    "total force when contact is only in -x half space");
                ASSERT_EQUAL(reported_target_regional_force(1), SimTK::Vec3(0), 1e-3, __FILE__, __LINE__,
                    "Expected the region 1 plane force to be 0 when contact is "
                    "in -x half space.");
                ASSERT_EQUAL(reported_target_regional_moment(0), reported_target_moment, 1e-3, __FILE__, __LINE__,
                    "Expected the region 0 plane moment to be equal to the "
                    "total force when contact is only in -x half space");
                ASSERT_EQUAL(reported_target_regional_moment(1), SimTK::Vec3(0), 1e-3, __FILE__, __LINE__,
                    "Expected the region 1 plane moment to be 0 when contact is "
                    "in -x half space.");
            }
            else {
                ASSERT_EQUAL(reported_target_regional_COpressure(1)(0), tx_values(i), 1e-3, __FILE__, __LINE__,
                    "Expected the region 1 plane COP x component to be "
                    "located at the center of the flat disk when contact is "
                    "in +x half space");
                ASSERT_EQUAL(reported_target_regional_COpressure(0), SimTK::Vec3(-1), 1e-3, __FILE__, __LINE__,
                    "Expected the region 0 plane COP to be -1 when contact is "
                    "in +x half space.");
                ASSERT_EQUAL(reported_target_regional_force(1), reported_target_force, 1e-3, __FILE__, __LINE__,
                    "Expected the region 0 plane force to be equal to the "
                    "total force when contact is only in +x half space");
                ASSERT_EQUAL(reported_target_regional_force(0), SimTK::Vec3(0), 1e-3, __FILE__, __LINE__,
                    "Expected the region 1 plane force to be 0 when contact is "
                    "in +x half space.");
                ASSERT_EQUAL(reported_target_regional_moment(1), reported_target_moment, 1e-3, __FILE__, __LINE__,
                    "Expected the region 0 plane moment to be equal to the "
                    "total force when contact is only in +x half space");
                ASSERT_EQUAL(reported_target_regional_moment(0), SimTK::Vec3(0), 1e-3, __FILE__, __LINE__,
                    "Expected the region 1 plane moment to be 0 when contact is "
                    "in +x half space.");
            }

            if (tz_values(j) < 0.0) {
                ASSERT_EQUAL(reported_target_regional_COpressure(4)(2), 
                    tz_values(j), 1e-3, __FILE__, __LINE__,
                    "Expected the region 4 plane COP z component to be "
                    "located at the center of the flat disk when contact is "
                    "in -z half space");
                ASSERT_EQUAL(reported_target_regional_COpressure(5), 
                    SimTK::Vec3(-1), 1e-3, __FILE__, __LINE__,
                    "Expected the region 5 plane COP to be -1 when contact is "
                    "in -z half space.");
                ASSERT_EQUAL(reported_target_regional_force(4), 
                    reported_target_force, 1e-3, __FILE__, __LINE__,
                    "Expected the region 4 plane force to be equal to the "
                    "total force when contact is only in -z half space");
                ASSERT_EQUAL(reported_target_regional_force(5), 
                    SimTK::Vec3(0), 1e-3, __FILE__, __LINE__,
                    "Expected the region 5 plane force to be 0 when contact is"
                    " in -z half space.");
                ASSERT_EQUAL(reported_target_regional_moment(4), 
                    reported_target_moment, 1e-3, __FILE__, __LINE__,
                    "Expected the region 4 plane moment to be equal to the "
                    "total force when contact is only in -z half space");
                ASSERT_EQUAL(reported_target_regional_moment(5), 
                    SimTK::Vec3(0), 1e-3, __FILE__, __LINE__,
                    "Expected the region 5 plane moment to be 0 when contact "
                    "is in -z half space.");
            }
            else {
                ASSERT_EQUAL(reported_target_regional_COpressure(5)(2), tz_values(j), 1e-3, __FILE__, __LINE__,
                    "Expected the region 5 plane COP z component to be "
                    "located at the center of the flat disk when contact is "
                    "in +z half space");
                ASSERT_EQUAL(reported_target_regional_COpressure(4), SimTK::Vec3(-1), 1e-3, __FILE__, __LINE__,
                    "Expected the region 4 plane COP to be -1 when contact is "
                    "in +z half space.");
                ASSERT_EQUAL(reported_target_regional_force(5), reported_target_force, 1e-3, __FILE__, __LINE__,
                    "Expected the region 4 plane force to be equal to the "
                    "total force when contact is only in +z half space");
                ASSERT_EQUAL(reported_target_regional_force(4), SimTK::Vec3(0), 1e-3, __FILE__, __LINE__,
                    "Expected the region 5 plane force to be 0 when contact is "
                    "in +z half space.");
                ASSERT_EQUAL(reported_target_regional_moment(5), reported_target_moment, 1e-3, __FILE__, __LINE__,
                    "Expected the region 4 plane moment to be equal to the "
                    "total force when contact is only in +z half space");
                ASSERT_EQUAL(reported_target_regional_moment(4), SimTK::Vec3(0), 1e-3, __FILE__, __LINE__,
                    "Expected the region 5 plane moment to be 0 when contact is "
                    "in +z half space.");
            }                
        }
    }

    // Pressure - depth relationship
    double depth1 = 0;
    double depth2 = 0.01;
    double depth3 = 0.02;

    tx.setValue(state, 0.0);
    tz.setValue(state, 0.0);
    ty.setValue(state, -depth1);

    model.realizeReport(state);

    double reported_casting_total_max_pressure1 = contact->getOutputValue<double>(
        state, "casting_total_max_pressure");
    double reported_target_total_max_pressure1 = contact->getOutputValue<double>(
        state, "target_total_max_pressure");

    ty.setValue(state, -depth2);
    model.realizeReport(state);

    double reported_casting_total_max_pressure2 = contact->getOutputValue<double>(
        state, "casting_total_max_pressure");
    double reported_target_total_max_pressure2 = contact->getOutputValue<double>(
        state, "target_total_max_pressure");

    ty.setValue(state, -depth3);
    model.realizeReport(state);

    double reported_casting_total_max_pressure3 = contact->getOutputValue<double>(
        state, "casting_total_max_pressure");
    double reported_target_total_max_pressure3 = contact->getOutputValue<double>(
        state, "target_total_max_pressure");

    double casting_pressure_diff_21 = 
        reported_casting_total_max_pressure2 - reported_casting_total_max_pressure1;
    double target_pressure_diff_21 = 
        reported_target_total_max_pressure2 - reported_target_total_max_pressure1;

    double casting_pressure_diff_32 =
        reported_casting_total_max_pressure3 - reported_casting_total_max_pressure2;
    double target_pressure_diff_32 =
        reported_target_total_max_pressure3 - reported_target_total_max_pressure2;

    ASSERT_EQUAL(casting_pressure_diff_21, casting_pressure_diff_32, 
        1e-3, __FILE__, __LINE__,
        "Expected the pressure depth relationship to be linear when "
        "elastic_foundation_formulation is 'linear'.");

    // Run Forward Simulation
    
    // two disks - different area
    // settled contact force = gravity
    // settled depth of model1 > model2

    for (Coordinate& coord : model.updComponentList<Coordinate>()) {
        coord.setValue(state, 0);
        coord.set_locked(true);
    }
    ty.set_locked(false);
    
    indenter->setMass(10.0);
    
    state = model.initSystem();
    ty.setValue(state, -0.01);
   
    SimTK::CPodesIntegrator integrator(model.getSystem(), 
        SimTK::CPodes::BDF, SimTK::CPodes::Newton);
    integrator.setAccuracy(1e-8);
    integrator.setMaximumStepSize(1e-3);
    SimTK::TimeStepper timestepper(model.getSystem(), integrator);
    timestepper.initialize(state);
    
    state.setTime(0.0);
    model.realizeReport(state);    

    double KE0 = model.calcKineticEnergy(state);
    double PE0 = model.calcPotentialEnergy(state);
    double E0 = KE0 + PE0;

    // Integrate Forward in Time
    double dt = 0.01;
    double stop_time = 2.0;
    int nSteps = (int)lround(stop_time / dt);

    StatesTrajectory states_trj;
    states_trj.append(state);

    for (int i = 0; i <= nSteps; ++i) {
        double t = (i + 1) * dt;

        timestepper.stepTo(t);
        state = timestepper.getState();        

        // Record Result
        model.realizeReport(state);
        states_trj.append(state);
    }

    double KE1 = model.calcKineticEnergy(state);
    double PE1 = model.calcPotentialEnergy(state);
    double E1 = KE1 + PE1;

    ASSERT_EQUAL(E0, E1,
        1e-3, __FILE__, __LINE__,
        "Expected energy to be conserved with contact force and no damping.");

    model.print("test.osim");
    STOFileAdapter().write(states_trj.exportToTable(model),"test_states.sto");

    //Add damper to settle contact
    SpringGeneralizedForce* sgf = new SpringGeneralizedForce("plane_indenter_ty");
    sgf->setName("damper");
    sgf->setViscosity(200.0);
    model.addForce(sgf);

    state = model.initSystem();
    //contact->setModelingOption(state, "flip_meshes", 1);
    ty.setValue(state, 0);
    SimTK::CPodesIntegrator integrator1(model.getSystem(),
        SimTK::CPodes::BDF, SimTK::CPodes::Newton);
    integrator1.setAccuracy(1e-4);
    integrator1.setMaximumStepSize(1e-3);
    SimTK::TimeStepper timestepper1(model.getSystem(), integrator1);
    timestepper1.initialize(state);

    StatesTrajectory states_trj1;
    states_trj1.append(state);

    for (int i = 0; i <= nSteps; ++i) {
        double t = (i + 1) * dt;

        timestepper1.stepTo(t);
        state = timestepper1.getState();

        // Record Result
        model.realizeReport(state);
        states_trj1.append(state);
    }
    model.print("test.osim");
    STOFileAdapter().write(states_trj1.exportToTable(model), "test_states1.sto");

    double settle_depth1 = ty.getValue(state);

    SimTK::Vec3 reported_casting_force = contact->getOutputValue<SimTK::Vec3>(
        state, "casting_total_contact_force");

    double weight = indenter->getMass() * gravity_acc;
    double KE2 = model.calcKineticEnergy(state);


    ASSERT_EQUAL(reported_casting_force.norm(), weight,
        1e-3, __FILE__, __LINE__,
        "Expected the settled indenter contact force to be equal to the "
        "weight of the indenter.");

    // TEST 3 - Ray intersection tests
    // ========================================================================
    // half sphere radius 10 cm, center at (0,0,0), top half (+y) deleted

    Smith2018ContactMesh mesh = 
        Smith2018ContactMesh("half_sphere","half_sphere_10cm_radius.stl");

    mesh.finalizeFromProperties();

    SimTK::Vec3 origin(0, -0.2, 0);
    SimTK::UnitVec3 ray(0, 0.5, 0);
    double min_proximity = 0.0;
    double max_proximity = 0.1;

    int tri_index;
    SimTK::Vec3 intersect_pt;
    double distance;
    
    bool ray_hit = mesh.rayIntersectMesh(origin, ray, 
            min_proximity, max_proximity, tri_index, intersect_pt, distance);
    
    ASSERT(ray_hit, __FILE__, __LINE__,
        "Expected the ray to intersect the half sphere. ");

    double expected_distance = 0.1; // origin y - sphere radius
    ASSERT_EQUAL(distance, expected_distance,
        1e-3, __FILE__, __LINE__,
        "The distance along the ray to contact with sphere to be equal to the "
        "ray origin minus the sphere radius.");

    // Flip ray
    ray = SimTK::UnitVec3(0,-0.5,0);

    ray_hit = mesh.rayIntersectMesh(origin, ray,
        min_proximity, max_proximity, tri_index, intersect_pt, distance);
    
    ASSERT(!ray_hit, __FILE__, __LINE__,
        "Expected the ray pointed away from half sphere to fail the "
        "rayIntersectMesh test. ");

    // move origin
    origin = SimTK::Vec3(0, 0, 0);
    ray_hit = mesh.rayIntersectMesh(origin, ray,
        min_proximity, max_proximity, tri_index, intersect_pt, distance);
    
    ASSERT(!ray_hit, __FILE__, __LINE__,
        "Expected the ray pointed from inside from half sphere to hit "
        "a back face and fail the rayIntersectMesh test. ");
        
    // Test 4 - test constructor from vertices and tri

    SimTK::PolygonalMesh ply_mesh;
    ply_mesh.loadFile("half_sphere_10cm_radius.stl");

    SimTK::Vector_<SimTK::Vec3> vertices(ply_mesh.getNumVertices());
    for (int v = 0; v < ply_mesh.getNumVertices(); ++v) {
        vertices.set(v, ply_mesh.getVertexPosition(v));
    }
    
    //SimTK::Vector_<SimTK::Vec3> triangles;
    std::vector<std::vector<int>> triangles;
    for (int t = 0; t < ply_mesh.getNumFaces(); ++t) {
        std::vector<int> face_vertex(3);
        for (int j = 0; j < 3; j++) {
             face_vertex[j] = ply_mesh.getFaceVertex(t, j);
        }
        triangles.push_back(face_vertex);
    }

    Smith2018ContactMesh t4_mesh = Smith2018ContactMesh("half_sphere",vertices,triangles);
    
    t4_mesh.finalizeFromProperties();

    origin = SimTK::Vec3 (0, -0.2, 0);
    ray = SimTK::UnitVec3 (0, 0.5, 0);

    ray_hit = t4_mesh.rayIntersectMesh(origin, ray,
        min_proximity, max_proximity, tri_index, intersect_pt, distance);
    
    ASSERT(ray_hit, __FILE__, __LINE__,
        "Expected the ray to intersect the half sphere. ");
}

void testSmith2018ContactDetection() {
    Model model = Model("C:\\Users\\csmith\\research\\dod1\\simulation\\DoD1_003\\model\\scaled_model_with_contact.osim");

    double gravity_acc = 9.8067;
    model.setGravity(SimTK::Vec3(0, -gravity_acc, 0));

    double mass = 1.0;
    double radius = 0.1;

    OpenSim::Body* sphere1 = new OpenSim::Body("sphere1", 0.0, SimTK::Vec3(0),
        mass * SimTK::Inertia::brick(0.05, 0.05, 0.05));

    OpenSim::Body* sphere2 = new OpenSim::Body("sphere2", mass, SimTK::Vec3(0),
        mass * SimTK::Inertia::brick(0.05, 0.05, 0.05));

    //indenter->attachGeometry(new Mesh("half_sphere_10cm_radius.stl"));

    model.addBody(sphere1);
    model.addBody(sphere2);

    SpatialTransform indenter_transform = SpatialTransform();
    indenter_transform.upd_rotation1().set_coordinates(0, "plane_indenter_rx");
    indenter_transform.upd_rotation1().set_axis(SimTK::Vec3(1, 0, 0));
    indenter_transform.upd_rotation1().set_function(LinearFunction());

    indenter_transform.upd_rotation2().set_coordinates(0, "plane_indenter_ry");
    indenter_transform.upd_rotation2().set_axis(SimTK::Vec3(0, 1, 0));
    indenter_transform.upd_rotation2().set_function(LinearFunction());

    indenter_transform.upd_rotation3().set_coordinates(0, "plane_indenter_rz");
    indenter_transform.upd_rotation3().set_axis(SimTK::Vec3(0, 0, 1));
    indenter_transform.upd_rotation3().set_function(LinearFunction());

    indenter_transform.upd_translation1().set_coordinates(0, "plane_indenter_tx");
    indenter_transform.upd_translation1().set_axis(SimTK::Vec3(1, 0, 0));
    indenter_transform.upd_translation1().set_function(LinearFunction());

    indenter_transform.upd_translation2().set_coordinates(0, "plane_indenter_ty");
    indenter_transform.upd_translation2().set_axis(SimTK::Vec3(0, 1, 0));
    indenter_transform.upd_translation2().set_function(LinearFunction());

    indenter_transform.upd_translation3().set_coordinates(0, "plane_indenter_tz");
    indenter_transform.upd_translation3().set_axis(SimTK::Vec3(0, 0, 1));
    indenter_transform.upd_translation3().set_function(LinearFunction());

    WeldJoint* ground_sphere1 = new WeldJoint("ground_sphere1", model.getGround(), *sphere1);
    CustomJoint* sphere1_sphere2 = new CustomJoint("sphere1_sphere2", *sphere1, *sphere2, indenter_transform);

    model.addJoint(ground_sphere1);
    model.addJoint(sphere1_sphere2);

   // Smith2018ContactMesh* sphere1_mesh = new Smith2018ContactMesh("sphere1", "sphere_10cm_radius.stl", *sphere1);
   // Smith2018ContactMesh* sphere2_mesh = new Smith2018ContactMesh("sphere2", "sphere_10cm_radius.stl", *sphere2);

    std::string fem_file = "C:\\Users\\csmith\\research\\hip_microinstability\\opensim\\1\\Geometry\\femur_cart.stl";
    std::string pel_file = "C:\\Users\\csmith\\research\\hip_microinstability\\opensim\\1\\Geometry\\pelvis_cart.stl";

    Smith2018ContactMesh* sphere1_mesh = new Smith2018ContactMesh("sphere1", fem_file, *sphere1);
    Smith2018ContactMesh* sphere2_mesh = new Smith2018ContactMesh("sphere2", pel_file, *sphere2);

    double E_indenter = 1e6;
    double v_indenter = 0.46;
    double E_plane = 1e6;
    double v_plane = 0.46;

    sphere1_mesh->set_elastic_modulus(E_plane);
    sphere1_mesh->set_poissons_ratio(v_plane);
    sphere1_mesh->set_thickness(0.2);

    sphere2_mesh->set_elastic_modulus(E_indenter);
    sphere2_mesh->set_poissons_ratio(v_indenter);
    sphere2_mesh->set_thickness(0.1);

    model.addContactGeometry(sphere1_mesh);
    model.addContactGeometry(sphere2_mesh);

    Smith2018ArticularContactForce* contact = new Smith2018ArticularContactForce("contact", *sphere1_mesh, *sphere2_mesh);
    contact->set_max_proximity(0.2);
    contact->set_elastic_foundation_formulation("linear");
    contact->set_use_lumped_contact_model(true);
    contact->set_use_fast_contact_detection(false);
    model.addForce(contact);

    model.setUseVisualizer(true);
    SimTK::State& state = model.initSystem();
    contact->setModelingOption(state, "flip_meshes", 1);

    model.print("contact_dectection_test.osim");

    Coordinate& tx = model.updCoordinateSet().get("plane_indenter_tx");
    Coordinate& ty = model.updCoordinateSet().get("plane_indenter_ty");
    Coordinate& tz = model.updCoordinateSet().get("plane_indenter_tz");
    Coordinate& rx = model.updCoordinateSet().get("plane_indenter_rx");
    Coordinate& ry = model.updCoordinateSet().get("plane_indenter_ry");
    Coordinate& rz = model.updCoordinateSet().get("plane_indenter_rz");

    model.realizeReport(state);

    const ModelVisualizer& viz = model.getVisualizer();
    viz.show(state);
    // Verify Contact Area
   /* double expected_casting_area = 2 * SimTK::Pi * SimTK::square(radius);
    double reported_casting_area = contact->getOutputValue<double>(
        state, "casting_total_contact_area");

    ASSERT_EQUAL(expected_casting_area, reported_casting_area, 1e-3,
        __FILE__, __LINE__,
        "Expected the contact area of the half sphere to be equal to the "
        "analytical formula.");

    double expected_target_area = SimTK::Pi * SimTK::square(radius);
    double reported_target_area =
        contact->getOutputValue<double>(state, "target_total_contact_area");

    ASSERT_EQUAL(expected_target_area, reported_target_area, 1e-3, __FILE__,
        __LINE__,
        "Expected the contact area of the half sphere to be equal to the "
        "analytical formula.");
        */
    //Move half sphere to half depth
    /*double depth = 0.05;
    ty.setValue(state, depth);
    model.realizeReport(state);
    viz.show(state);*/

    for (double d = 0; d < 0.01; d += 0.001) {
        std::cout << d << std::endl;
        tz.setValue(state, d);
        model.realizeReport(state);
        viz.show(state);
    }

    // Max Proximity
   /* double reported_casting_max_prox = contact->getOutputValue<double>(
        state, "casting_total_max_proximity");

    double reported_target_max_prox = contact->getOutputValue<double>(
        state, "target_total_max_proximity");

    double expected_max_prox = radius / 2;

    ASSERT_EQUAL(expected_max_prox, reported_casting_max_prox, 1e-3, __FILE__,
        __LINE__,
        "Expected the maximum casting proximity "
        "to be equal to the half sphere radius.");

    ASSERT_EQUAL(expected_max_prox, reported_target_max_prox, 1e-3, __FILE__,
        __LINE__,
        "Expected the maximum target proximity "
        "to be equal to the half sphere radius.");*/
}
