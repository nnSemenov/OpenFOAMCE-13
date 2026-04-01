lnInclude(${CMAKE_CURRENT_SOURCE_DIR})

set(target_name rigidBodyDynamics)

add_library(${target_name} SHARED
    bodies/rigidBody/rigidBody.C
    bodies/masslessBody/masslessBody.C
    bodies/jointBody/jointBody.C
    bodies/compositeBody/compositeBody.C
    bodies/subBody/subBody.C
    bodies/sphere/sphere.C
    bodies/cuboid/cuboid.C
    bodies/pointMasses/pointMasses.C

    joints/joint/joint.C
    joints/null/nullJoint.C
    joints/composite/compositeJoint.C
    joints/floating/floatingJoint.C

    joints/Rx/Rx.C
    joints/Ry/Ry.C
    joints/Rz/Rz.C
    joints/Ra/Ra.C
    joints/Rs/Rs.C
    joints/Rzyx/Rzyx.C
    joints/Rxyz/Rxyz.C
    joints/Ryxz/Ryxz.C

    joints/Px/Px.C
    joints/Py/Py.C
    joints/Pz/Pz.C
    joints/Pa/Pa.C
    joints/Pxyz/Pxyz.C

    joints/rigid/rigid.C
    joints/rotating/rotating.C
    joints/function/function.C
    joints/functionDot/functionDot.C

    restraints/restraint/rigidBodyRestraint.C
    restraints/restraint/rigidBodyRestraintNew.C
    restraints/linearSpring/linearSpring.C
    restraints/linearDamper/linearDamper.C
    restraints/linearAxialAngularSpring/linearAxialAngularSpring.C
    restraints/sphericalAngularDamper/sphericalAngularDamper.C
    restraints/externalForce/externalForce.C

    rigidBodyModel/rigidBodyModel.C
    rigidBodyModel/forwardDynamics.C

    rigidBodyModelState/rigidBodyModelState.C
    rigidBodyModelState/rigidBodyModelStateIO.C

    rigidBodyMotion/rigidBodyMotion.C
    rigidBodyMotion/rigidBodyMotionIO.C

    rigidBodySolvers/rigidBodySolver/rigidBodySolver.C
    rigidBodySolvers/rigidBodySolver/rigidBodySolverNew.C
    rigidBodySolvers/symplectic/symplectic.C
    rigidBodySolvers/Newmark/Newmark.C
    rigidBodySolvers/CrankNicolson/CrankNicolson.C
)

target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude"
)

target_include_directories(${target_name} PRIVATE

)

target_link_libraries(${target_name} PRIVATE
    OpenFOAM_Defines
    OpenFOAM
    OSspecific
)

target_compile_definitions(${target_name} PRIVATE
    LIB_NAME=lib${target_name}.so
)

set(FOAM_regular_libraries "${FOAM_regular_libraries};${target_name}" PARENT_SCOPE)


