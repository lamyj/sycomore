import numpy
from .units import *

from . import gamma

def pulse(angle, phase):
    M = numpy.identity(4)
    
    Rz = axis_angle_to_matrix([0,0,1], -phase.convert_to(rad))
    Rx = axis_angle_to_matrix([1,0,0], angle.convert_to(rad))
    M[:3,:3] = numpy.matmul(numpy.linalg.inv(Rz), numpy.matmul(Rx, Rz))
    
    return M

def time_interval(
        species, duration, delta_omega=0*Hz, 
        gradient_amplitude=0*T/m, position=0*m):
    
    E = relaxation(species, duration)
    
    delta_omega = (
        # Field-related dephasing
        delta_omega 
        # Species-related dephasing, e.g. chemical shift or susceptibility
        + species.delta_omega 
        # Gradient-related dephasing
        + gamma*numpy.dot(gradient_amplitude, position)
    )
    F = phase_accumulation(duration * 2*numpy.pi*rad * delta_omega)
    
    return numpy.matmul(F, E)

def relaxation(species, duration):
    E_1 = numpy.exp((-duration*species.R1).magnitude)
    E_2 = numpy.exp((-duration*species.R2).magnitude)
    
    M = numpy.diag([E_2, E_2, E_1, 1])
    M[2,3] = 1-E_1
    
    return M

def phase_accumulation(angle):
    M = numpy.identity(4)
    M[:3,:3] = axis_angle_to_matrix([0,0,1], angle.convert_to(rad))
    return M

def axis_angle_to_matrix(axis, angle) :
    r""" Convert an (axis, angle) to a rotation matrix.
    
         This formula comes from Rodrigues' rotation formula,
         :math:`R = I + \hat{\omega} \sin \theta + \hat{\omega}^2 (1-\cos \theta)`
         where :math:`\hat{}` gives the antisymmetric matrix equivalent of the cross product
        
         .. math ::
            
             \hat{\omega} = \begin{matrix}
                                        0 & -\omega_z &  \omega_y \\
                                 \omega_z &         0 & -\omega_x \\
                                -\omega_y &  \omega_x &         0 \\
                            \end{matrix} 
        
         Diagonal terms can be rewritten :
         
         .. math ::
             
             \begin{matrix}
                 1+(1-\cos \theta)*(\omega_x^2-1) & = & 1+(1-\cos \theta)*\omega_x^2-(1-\cos \theta) \\
                                                  & = & \cos \theta+\omega_x^2*(1-\cos \theta)
             \end{matrix}
    """
    
    result = numpy.ndarray((3,3))
    
    cos = numpy.cos(angle)
    sin = numpy.sin(angle)
    one_minus_cos = 1.-cos
    
    result[0][0] = cos+axis[0]**2*(one_minus_cos)
    result[1][1] = cos+axis[1]**2*(one_minus_cos)
    result[2][2] = cos+axis[2]**2*(one_minus_cos)
    
    result[0][1] = -axis[2]*sin+axis[0]*axis[1]*one_minus_cos
    result[1][0] = +axis[2]*sin+axis[0]*axis[1]*one_minus_cos
    
    result[0][2] = +axis[1]*sin+axis[0]*axis[2]*one_minus_cos
    result[2][0] = -axis[1]*sin+axis[0]*axis[2]*one_minus_cos
    
    result[1][2] = -axis[0]*sin+axis[1]*axis[2]*one_minus_cos
    result[2][1] = +axis[0]*sin+axis[1]*axis[2]*one_minus_cos
    
    return result
