# -*- coding: utf-8 -*-

# Nodo ROS para el brazo articulado.
# --------------------------------------------------------
# autor: Gabriel Montoiro Varela
# ver:   0.1

import roslib
import rospy
from sensor_msgs.msg import JointState
import numpy as np
from my_robot_arm_controller import MyRoboticArm

class MyRoboticArmNode():
    """Crea un nodo ROS para publicar el estado del robot."""

    def __init__(self, arm):
        """Crea el nodo ROS para la publicación del estado del brazo articulado."""

        self.arm = arm

	# iniciamos el nodo ROS
        rospy.init_node('my_robotic_arm')

	# creamos el publicador del estado de las articulaciones
	self.state_pub = rospy.Publisher('/joint_states', JointState, queue_size=10)

    def publish_state(self):
        """Publica un nuevo mensaje ROS con el estado actual del robot."""

        msg = JointState()

        msg.header.stamp = rospy.Time.now()
        msg.header.frame_id = 'base_link'
	msg.name = []
	"""
        msg.name = ['base_link__link_01', \
                    'link_01__link_02', \
                    'link_02__link_03', \
                    'link_03__link_04', \
                    'link_04__link_05', \
                    'link_05__gripper']
	"""
        msg.position = list(self.arm.get_theta())
        msg.velocity = []
        msg.effort = []
        
        self.state_pub.publish(msg)
        


