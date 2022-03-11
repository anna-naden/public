""" classes and functions for entanglement """
from sympy import Matrix
import sympy as sp
import numpy as np
from scipy import linalg
def get_states(rho):
    """ Extract the state of a particle from its density matrix
    by getting its eigenvectors. Only the states that have
    nonzero probability are retained. """
    matrix = Matrix(rho)
    states = []
    for triple in matrix.eigenvects():
        prob = triple[0]
        if prob != 0:
            evects = triple[2]
            for evec in evects:
                evec1 = []
                for component in evec:
                    evec1.append(component)
                states.append(np.array(evec1))
    return np.array(states)
def initial_state():
    """ Create initial state for teleportation """
    asym, bsym = sp.symbols('a b')
    state = [asym,bsym,0,0,0,0,0,0]
    return Qstate(np.array(state))
def initial_state2():
    """ Create initial state for teleportation """
    asym, bsym = sp.symbols('a b')
    state = np.zeros(32, dtype=type(asym))
    state[0]=asym
    state[1]=bsym
    return Qstate(state)
def partial_trace_first_subsystem(rho,first,second):
    """
    Partial trace over second subsystem
    first=#qubits in first subsyste
    second=# qubits in second subsystem
    first+second=dim
    """
    first=2**first
    second=2**second
    return np.trace(rho.reshape(second,first,second,first), axis1=1, axis2=3)
def partial_trace_second_subsystem(rho,first,second):
    """
    Partial trace over first subsystem
    first=#qubits in first subsyste
    second=# qubits in second subsystem
    first+second=dim

    """
    first = 2**first
    second = 2**second
    return np.trace(rho.reshape(second,first,second,first), axis1=0, axis2=2)
class Tensor():
    """ Like qutip's Qobj """
    def __init__(self, val):
        self.val = val
    def __mul__(self,other):
        if isinstance(other, Tensor):
            return Tensor(linalg.kron(self.val, other.val))
        return Qstate(np.matmul(other.val,np.conj(self.val)))
    def __add__(self,other):
        return Tensor(self.val+other.val)
class Qstate():
    """ Like qutip's qobj """
    def __init__(self,val):
        self.val = val
    def __rmul__(self, other):
        return Qstate(np.matmul(self,other.val))
    def braket(self,other):
        """ Like qutip's ket2dm"""
        return Tensor(np.outer(other.val,np.conj(self.val)))
