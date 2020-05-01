#include "Common.h"
#include "QuadControl.h"

#include "Utility/SimpleConfig.h"

#include "Utility/StringUtils.h"
#include "Trajectory.h"
#include "BaseController.h"
#include "Math/Mat3x3F.h"

#ifdef __PX4_NUTTX
#include <systemlib/param/param.h>
#endif

enum RotorIndices
{
    ROTOR_FL = 0,
    ROTOR_FR,
    ROTOR_RL,
    ROTOR_RR,
    ROTOR_NUM
};

void QuadControl::Init()
{
    BaseController::Init();

    // variables needed for integral control
    integratedAltitudeError = 0;

#ifndef __PX4_NUTTX
    // Load params from simulator parameter system
    ParamsHandle config = SimpleConfig::GetInstance();

    // Load parameters (default to 0)
    kpPosXY = config->Get(_config+".kpPosXY", 0);
    kpPosZ = config->Get(_config + ".kpPosZ", 0);
    KiPosZ = config->Get(_config + ".KiPosZ", 0);

    kpVelXY = config->Get(_config + ".kpVelXY", 0);
    kpVelZ = config->Get(_config + ".kpVelZ", 0);

    kpBank = config->Get(_config + ".kpBank", 0);
    kpYaw = config->Get(_config + ".kpYaw", 0);

    kpPQR = config->Get(_config + ".kpPQR", V3F());

    maxDescentRate = config->Get(_config + ".maxDescentRate", 100);
    maxAscentRate = config->Get(_config + ".maxAscentRate", 100);
    maxSpeedXY = config->Get(_config + ".maxSpeedXY", 100);
    maxAccelXY = config->Get(_config + ".maxHorizAccel", 100);

    maxTiltAngle = config->Get(_config + ".maxTiltAngle", 100);

    minMotorThrust = config->Get(_config + ".minMotorThrust", 0);
    maxMotorThrust = config->Get(_config + ".maxMotorThrust", 100);
#else
    // load params from PX4 parameter system
    //TODO
    param_get(param_find("MC_PITCH_P"), &Kp_bank);
    param_get(param_find("MC_YAW_P"), &Kp_yaw);
#endif
}

VehicleCommand QuadControl::GenerateMotorCommands(float collThrustCmd, V3F momentCmd)
{
    // Convert a desired 3-axis moment and collective thrust command to
    //   individual motor thrust commands
    // INPUTS:
    //   collThrustCmd: desired collective thrust [N]
    //   momentCmd: desired rotation moment about each axis [N m]
    // OUTPUT:
    //   set class member variable cmd (class variable for graphing) where
    //   cmd.desiredThrustsN[0..3]: motor commands, in [N]

    // HINTS:
    // - you can access parts of momentCmd via e.g. momentCmd.x
    // You'll need the arm length parameter L, and the drag/thrust ratio kappa

    ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

    //cmd.desiredThrustsN[0] = mass * 9.81f / 4.f; // front left
    //cmd.desiredThrustsN[1] = mass * 9.81f / 4.f; // front right
    //cmd.desiredThrustsN[2] = mass * 9.81f / 4.f; // rear left
    //cmd.desiredThrustsN[3] = mass * 9.81f / 4.f; // rear right

    const float Ln = this->L / std::sqrt(2.0F);

    const float Fc_N{ collThrustCmd };
    const float Fp_N{ momentCmd.x / Ln };
    const float Fq_N{ momentCmd.y / Ln };
    const float Fr_N{ -momentCmd.z / this->kappa };

    // rotors 3 and 4 are swapped compared to the lessons.
    float F_N[ROTOR_NUM];
    F_N[ROTOR_FL] = 0.25F * (Fc_N + Fp_N + Fq_N + Fr_N);
    F_N[ROTOR_FR] = 0.25F * (Fc_N - Fp_N + Fq_N - Fr_N);
    F_N[ROTOR_RL] = 0.25F * (Fc_N + Fp_N - Fq_N - Fr_N);
    F_N[ROTOR_RR] = 0.25F * (Fc_N - Fp_N - Fq_N + Fr_N);

    cmd.desiredThrustsN[ROTOR_FL] = CONSTRAIN(F_N[ROTOR_FL], this->minMotorThrust, this->maxMotorThrust);
    cmd.desiredThrustsN[ROTOR_FR] = CONSTRAIN(F_N[ROTOR_FR], this->minMotorThrust, this->maxMotorThrust);
    cmd.desiredThrustsN[ROTOR_RL] = CONSTRAIN(F_N[ROTOR_RL], this->minMotorThrust, this->maxMotorThrust);
    cmd.desiredThrustsN[ROTOR_RR] = CONSTRAIN(F_N[ROTOR_RR], this->minMotorThrust, this->maxMotorThrust);

    /////////////////////////////// END STUDENT CODE ////////////////////////////

    return cmd;
}

V3F QuadControl::BodyRateControl(V3F pqrCmd, V3F pqr)
{
    // Calculate a desired 3-axis moment given a desired and current body rate
    // INPUTS:
    //   pqrCmd: desired body rates [rad/s]
    //   pqr: current or estimated body rates [rad/s]
    // OUTPUT:
    //   return a V3F containing the desired moments for each of the 3 axes

    // HINTS:
    //  - you can use V3Fs just like scalars: V3F a(1,1,1), b(2,3,4), c; c=a-b;
    //  - you'll need parameters for moments of inertia Ixx, Iyy, Izz
    //  - you'll also need the gain parameter kpPQR (it's a V3F)

    V3F momentCmd;

    ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

    // p_dot_dot_cmd = up_bar = Kp_p * (p_target - p_measured)
    // tau_x_cmd = p_dot_dot_cmd * I_xx

    // q_dot_dot_cmd = uq_bar = Kp_q * (q_target - q_measured)
    // tau_y_cmd = q_dot_dot_cmd * I_yy

    // r_dot_dot_cmd = ur_bar = Kp_r * (r_target - r_measured)
    // tau_z_cmd = r_dot_dot_cmd * I_zz

    V3F I(this->Ixx, this->Iyy, this->Izz);
    V3F PQRError = pqrCmd - pqr;
    V3F PQRDotDot = this->kpPQR * PQRError;

    momentCmd = PQRDotDot * I;

    /////////////////////////////// END STUDENT CODE ////////////////////////////

    return momentCmd;
}

// returns a desired roll and pitch rate
V3F QuadControl::RollPitchControl(V3F accelCmd, Quaternion<float> attitude, float collThrustCmd)
{
    // Calculate a desired pitch and roll angle rates based on a desired global
    //   lateral acceleration, the current attitude of the quad, and desired
    //   collective thrust command
    // INPUTS:
    //   accelCmd: desired acceleration in global XY coordinates [m/s2]
    //   attitude: current or estimated attitude of the vehicle
    //   collThrustCmd: desired collective thrust of the quad [N]
    // OUTPUT:
    //   return a V3F containing the desired pitch and roll rates. The Z
    //     element of the V3F should be left at its default value (0)

    // HINTS:
    //  - we already provide rotation matrix R: to get element R[1,2] (python) use R(1,2) (C++)
    //  - you'll need the roll/pitch gain kpBank
    //  - collThrustCmd is a force in Newtons! You'll likely want to convert it to acceleration first

    V3F pqrCmd;
    Mat3x3F R = attitude.RotationMatrix_IwrtB();

    ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

    float R11{ R(0, 0) };
    float R12{ R(0, 1) };
    float R13{ R(0, 2) };
    float R21{ R(1, 0) };
    float R22{ R(1, 1) };
    float R23{ R(1, 2) };
    float R33{ R(2, 2) };

    float bx_m = R13;
    float by_m = R23;

    float c_mps2 = -collThrustCmd / this->mass;
    float bx_c = CONSTRAIN(accelCmd.x / c_mps2, -this->maxTiltAngle, this->maxTiltAngle);
    float by_c = CONSTRAIN(accelCmd.y / c_mps2, -this->maxTiltAngle, this->maxTiltAngle);

    if (collThrustCmd < 0)
    {
        bx_c = 0.0F;
        by_c = 0.0F;
    }

    float bxError = bx_c - bx_m;
    float byError = by_c - by_m;

    float bxDot_c = this->kpBank * bxError;
    float byDot_c = this->kpBank * byError;

    // convert the commanded orientation rates from inertial frame to body frame

    pqrCmd.x = (R21 * bxDot_c) - (R11 * byDot_c);
    pqrCmd.x /= R33;

    pqrCmd.y = (R22 * bxDot_c) - (R12 * byDot_c);
    pqrCmd.y /= R33;

    pqrCmd.z = 0.0F;

    /////////////////////////////// END STUDENT CODE ////////////////////////////

    return pqrCmd;
}

float QuadControl::AltitudeControl(float posZCmd, float velZCmd, float posZ, float velZ, Quaternion<float> attitude, float accelZCmd, float dt)
{
    // Calculate desired quad thrust based on altitude setpoint, actual altitude,
    //   vertical velocity setpoint, actual vertical velocity, and a vertical
    //   acceleration feed-forward command
    // INPUTS:
    //   posZCmd, velZCmd: desired vertical position and velocity in NED [m]
    //   posZ, velZ: current vertical position and velocity in NED [m]
    //   accelZCmd: feed-forward vertical acceleration in NED [m/s2]
    //   dt: the time step of the measurements [seconds]
    // OUTPUT:
    //   return a collective thrust command in [N]

    // HINTS:
    //  - we already provide rotation matrix R: to get element R[1,2] (python) use R(1,2) (C++)
    //  - you'll need the gain parameters kpPosZ and kpVelZ
    //  - maxAscentRate and maxDescentRate are maximum vertical speeds. Note they're both >=0!
    //  - make sure to return a force, not an acceleration
    //  - remember that for an upright quad in NED, thrust should be HIGHER if the desired Z acceleration is LOWER

    Mat3x3F R = attitude.RotationMatrix_IwrtB();
    float thrust = 0;

    ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

    float bz{ R(2, 2) };

    float zPosError = posZCmd - posZ;
    float zVelError = velZCmd - velZ;
    this->integratedAltitudeError += zPosError * dt;

    float termP = this->kpPosZ * zPosError;
    float termI = this->KiPosZ * this->integratedAltitudeError;
    float termD = this->kpVelZ * zVelError;

    float u1Bar = termP + termI + termD + accelZCmd;
    u1Bar = (u1Bar - CONST_GRAVITY) / bz;
    u1Bar = CONSTRAIN(u1Bar, -this->maxAscentRate / dt, this->maxDescentRate / dt);

    thrust = -(this->mass * u1Bar);

    /////////////////////////////// END STUDENT CODE ////////////////////////////

    return thrust;
}

// returns a desired acceleration in global frame
V3F QuadControl::LateralPositionControl(V3F posCmd, V3F velCmd, V3F pos, V3F vel, V3F accelCmdFF)
{
    // Calculate a desired horizontal acceleration based on
    //  desired lateral position/velocity/acceleration and current pose
    // INPUTS:
    //   posCmd: desired position, in NED [m]
    //   velCmd: desired velocity, in NED [m/s]
    //   pos: current position, NED [m]
    //   vel: current velocity, NED [m/s]
    //   accelCmdFF: feed-forward acceleration, NED [m/s2]
    // OUTPUT:
    //   return a V3F with desired horizontal accelerations.
    //     the Z component should be 0
    // HINTS:
    //  - use the gain parameters kpPosXY and kpVelXY
    //  - make sure you limit the maximum horizontal velocity and acceleration
    //    to maxSpeedXY and maxAccelXY

    // make sure we don't have any incoming z-component
    accelCmdFF.z = 0;
    velCmd.z = 0;
    posCmd.z = pos.z;

    // we initialize the returned desired acceleration to the feed-forward value.
    // Make sure to _add_, not simply replace, the result of your controller
    // to this variable
    V3F accelCmd = accelCmdFF;

    ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

    velCmd.x = CONSTRAIN(velCmd.x, -this->maxSpeedXY, this->maxSpeedXY);
    velCmd.y = CONSTRAIN(velCmd.y, -this->maxSpeedXY, this->maxSpeedXY);

    float xAccelCmd_mps2 = this->kpPosXY * (posCmd.x - pos.x) + this->kpVelXY * (velCmd.x - vel.x);
    float yAccelCmd_mps2 = this->kpPosXY * (posCmd.y - pos.y) + this->kpVelXY * (velCmd.y - vel.y);

    // accumulate to consider the feedforward term from initialization

    accelCmd.x += xAccelCmd_mps2;
    accelCmd.y += yAccelCmd_mps2;

    accelCmd.x = CONSTRAIN(accelCmd.x, -this->maxAccelXY, this->maxAccelXY);
    accelCmd.y = CONSTRAIN(accelCmd.y, -this->maxAccelXY, this->maxAccelXY);

    accelCmd.z = 0.0F;

    /////////////////////////////// END STUDENT CODE ////////////////////////////

    return accelCmd;
}

double constrainAngle(double x)
{
    x = fmodf(x + M_PI, 2 * M_PI);
    if (x < 0)
        x += 2 * M_PI;
    return x - M_PI;
}

// returns desired yaw rate
float QuadControl::YawControl(float yawCmd, float yaw)
{
    // Calculate a desired yaw rate to control yaw to yawCmd
    // INPUTS:
    //   yawCmd: commanded yaw [rad]
    //   yaw: current yaw [rad]
    // OUTPUT:
    //   return a desired yaw rate [rad/s]
    // HINTS:
    //  - use fmodf(foo,b) to unwrap a radian angle measure float foo to range [0,b].
    //  - use the yaw control gain parameter kpYaw

    float yawRateCmd=0;
    ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

    yawCmd = constrainAngle(yawCmd);

    float yawError_rad = yawCmd - yaw;

    yawError_rad = constrainAngle(yawError_rad);

    yawRateCmd = this->kpYaw * yawError_rad;

    /////////////////////////////// END STUDENT CODE ////////////////////////////

    return yawRateCmd;
}

VehicleCommand QuadControl::RunControl(float dt, float simTime)
{
    curTrajPoint = GetNextTrajectoryPoint(simTime);

    float collThrustCmd = AltitudeControl(curTrajPoint.position.z, curTrajPoint.velocity.z, estPos.z, estVel.z, estAtt, curTrajPoint.accel.z, dt);

    // reserve some thrust margin for angle control
    float thrustMargin = .1f*(maxMotorThrust - minMotorThrust);
    collThrustCmd = CONSTRAIN(collThrustCmd, (minMotorThrust+ thrustMargin)*4.f, (maxMotorThrust-thrustMargin)*4.f);

    V3F desAcc = LateralPositionControl(curTrajPoint.position, curTrajPoint.velocity, estPos, estVel, curTrajPoint.accel);

    V3F desOmega = RollPitchControl(desAcc, estAtt, collThrustCmd);
    desOmega.z = YawControl(curTrajPoint.attitude.Yaw(), estAtt.Yaw());

    V3F desMoment = BodyRateControl(desOmega, estOmega);

    return GenerateMotorCommands(collThrustCmd, desMoment);
}
