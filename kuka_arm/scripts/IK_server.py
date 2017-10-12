#!/usr/bin/env python

# Copyright (C) 2017 Electric Movement Inc.
#
# This file is part of Robotic Arm: Pick and Place project for Udacity
# Robotics nano-degree program
#
# All Rights Reserved.

# Author: Harsh Pandya

# import modules
import rospy
import tf
from kuka_arm.srv import *
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint
from geometry_msgs.msg import Pose
from mpmath import *
from sympy import *


def handle_calculate_IK(req):
    rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
    if len(req.poses) < 1:
        print "No valid poses received"
        return -1
    else:
		
        ### Your FK code here
        # Create symbols
	q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8') # theta_i
	d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8')
	a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7')
	alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7')
	#
	#   
	# Create Modified DH parameters
	s = {alpha0:     0, a0:         0, d1:  0.75,
	     alpha1: -pi/2, a1:      0.35, d2:     0, q2: q2-pi/2,
	     alpha2:     0, a2:      1.25, d3:     0,
	     alpha3: -pi/2, a3:    -0.054, d4:   1.5,
	     alpha4:  pi/2, a4:         0, d5:     0,
	     alpha5: -pi/2, a5:         0, d6:     0,}

	#
	#            
	# Define Modified DH Transformation matrix
        def TF_Matrix(alpha, a, d, q):
            TF = Matrix([[ cos(q),             -sin(q),            0,            a],
                         [ sin(q)*cos(alpha), cos(q)*cos(alpha), -sin(alpha), -sin(alpha)*d],
                         [ sin(q)*sin(alpha), cos(q)*sin(alpha), cos(alpha), cos(alpha)*d],
                         [ 0,                   0,                   0,           1]])
       	    return TF
	#
	#
	# Create individual transformation matrices
        T0_1 = TF_Matrix(alpha0, a0, d1, q1).subs(s)
        T1_2 = TF_Matrix(alpha1, a1, d2, q2).subs(s)
        T2_3 = TF_Matrix(alpha2, a2, d3, q3).subs(s)
        T3_4 = TF_Matrix(alpha3, a3, d4, q4).subs(s)
        T4_5 = TF_Matrix(alpha4, a4, d5, q5).subs(s)
        T5_6 = TF_Matrix(alpha5, a5, d6, q6).subs(s)
        T6_EE = TF_Matrix(alpha6, a6, d7, q7).subs(s)
    
        T0_EE = T0_1 * T1_2 * T2_3 * T3_4 * T4_5 * T5_6 * T6_EE
	#
	#
	# Extract rotation matrices from the transformation matrices
	#
	#
        ###

        # Initialize service response
        joint_trajectory_list = []
        for x in xrange(0, len(req.poses)):
            # IK code starts here
            joint_trajectory_point = JointTrajectoryPoint()

	    # Extract end-effector position and orientation from request
	    # px,py,pz = end-effector position
	    # roll, pitch, yaw = end-effector orientation
            px = req.poses[x].position.x
            py = req.poses[x].position.y
            pz = req.poses[x].position.z

            (roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
                [req.poses[x].orientation.x, req.poses[x].orientation.y,
                    req.poses[x].orientation.z, req.poses[x].orientation.w])
     
            ### Your IK code here 
	    # Compensate for rotation discrepancy between DH parameters and Gazebo
            r, p, y = symbols('r p y')
            
            # Find end-effector rotation matrix
            # Roll
            ROT_x = Matrix([[1,      0,       0],
                            [0, cos(r), -sin(r)],
              	            [0, sin(r),  cos(r)]]) 
            # Pitch
            ROT_y = Matrix([[cos(p),  0, sin(p)],
                            [0,       1,      0],
            	            [-sin(p), 0, cos(p)]]) 
            
            # Yaw
            ROT_z = Matrix([[cos(y),  -sin(y), 0],
                            [sin(y),   cos(y), 0],
            	            [     0,        0, 1]]) 
            # Correct for rotational differences between DH parameters and URDF definitions
            Rotation_ERR = ROT_z.subs(y, radians(180)) * ROT_y.subs(p, radians(-90))
	    # Get the desired orientation for the end effector
            ROT_EE = ROT_z * ROT_y * ROT_x
            ROT_EE = ROT_EE * Rotation_ERR
            ROT_EE = ROT_EE.subs({'r': roll, 'p': pitch, 'y': yaw})
            
            EE = Matrix([[px], [py], [pz]])
            # Get the vector for the wrist center wrt the base frame
	    # Take the end effector vector and subtract the distance from end effect to wrist
	    # center, rotated into place with the rotation matrix from frame 0 to frame 6
	    # r_wc = r_ee - d_ee * R
            WC = EE - (0.303) * ROT_EE[:,2]
            #
            #
            # Calculate joint angles using Geometric IK method

            # Get joint angles with atan2
            # theta1 controls movement in the yaw direction, so compute 
            # theta1 = tan^-1(y, x)
            theta1 = atan2(WC[1], WC[0])
            
            # Get theta2, theta3 with SSS triangle
            # side_a comes from a triangle, with the hypotenuse as the straight
            # line distance from joint 3 to joint 5
            side_a = 1.501
            # side_b is hardest to compute - it is obtained from
            # its components.
            # wc[2] is the wrist center z-component from the base frame, so
            # subtract 0.75 to get the z-distance starting from joint 2
            # wc[0] and w[1] are the x and y components of the
            # wrist center (as seen from top-down view). By using
            # the pythagorean theorem, we can calculate the
            # hypotenuse. However, joint 2 juts out from joint 1,
            # so subtract this constant offset (0.35).
            # Then, with the x and y components of side B, we
            # compute it with the pythagorean theorem.
            side_b = sqrt(pow((sqrt(WC[0] * WC[0] + WC[1] * WC[1]) - 0.35), 2) + pow((WC[2] - 0.75), 2))
            # side_c is the link length between joint 2 and joint 3
            side_c = 1.25
            
            # By the Law of Cosines:
            # cos(A) = (b^2 + c^2 - a^2) / 2*b*c
            # cos(B) = (c^2 + a^2 - b^2) / 2*c*a
            # cos(C) = (a^2 + b^2 - c^2) / 2*a*b
            angle_a = acos((side_b * side_b + side_c * side_c - side_a * side_a) / (2 * side_b * side_c))
            angle_b = acos((side_a * side_a + side_c * side_c - side_b * side_b) / (2 * side_a * side_c))
            angle_c = acos((side_a * side_a + side_b * side_b - side_c * side_c) / (2 * side_a * side_b))
            
            # We know theta2, angle_a, the angle for side B add
            # up to 90 degrees. Therefore, compute the side B
            # angle with atan2 and then solve for theta2.
            theta2 = pi / 2 - angle_a - atan2(WC[2] - 0.75, sqrt(WC[0] *
            WC[0] + WC[1] * WC[1]) - 0.35)
            # Above the sag angle caused by the dip between joint 3 and 4, you
            # can obtain an upper right triangle with the horizontal. The angle
            # between this horizontal and the link between joints 2 and 3 is 90
            # degrees when theta2 is zero. Therefore, you can obtain theta3 by
            # taking 90 degrees and subtracting that angle and angle_b.
            theta3 = pi / 2 - (angle_b + 0.036) 
            
            # Extract the rotation matrix from the individual transforms and
            # multiply them up to get the rotation matrix from frame 3 to the
            # base frame
            R0_3 = T0_1[0:3, 0:3] * T1_2[0:3, 0:3] * T2_3[0:3, 0:3]
            R0_3 = R0_3.evalf(subs={q1: theta1, q2: theta2, q3: theta3})
            # ROT_EE is the rotation matrix from base frame to gripper obtained
            # from desired roll, pitch, yaw specified in request.
            #
            # For this calculation, just use the tranpose, as the inverse is
            # equal to the transpose for rotation matrices.
            #R3_6 = R0_3.inv("LU") * ROT_EE
            R3_6 = R0_3.T * ROT_EE
            
            # Print the symbolic matrix and solve the linear equations to
            # obtain solutions for thetas 4-6
            theta4 = atan2(R3_6[2,2], -R3_6[0,2])
            theta5 = atan2(sqrt(R3_6[0,2]*R3_6[0,2] + R3_6[2,2] *
            R3_6[2,2]), R3_6[1,2])
            theta6 = atan2(-R3_6[1,1], R3_6[1,0])
            #
            #
            ###
        		
            # Populate response for the IK request
            # In the next line replace theta1,theta2...,theta6 by your joint angle variables
	    joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5, theta6]
	    joint_trajectory_list.append(joint_trajectory_point)

        rospy.loginfo("length of Joint Trajectory List: %s" % len(joint_trajectory_list))
        return CalculateIKResponse(joint_trajectory_list)


def IK_server():
    # initialize node and declare calculate_ik service
    rospy.init_node('IK_server')
    s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
    print "Ready to receive an IK request"
    rospy.spin()

if __name__ == "__main__":
    IK_server()
