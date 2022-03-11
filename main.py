""" demonstration of entanglement """
import numpy as np
import sympy as sp
from classes import Tensor, Qstate, initial_state, initial_state2, \
    get_states, partial_trace_first_subsystem

 # Quantum gates
identity = Tensor(np.array([[1,0],[0,1]]))
x_gate = Tensor(np.array([[0,1],[1,0]]))
y_gate = Tensor(np.array([[0,-complex(0,1)],[complex(0,1),0]], dtype=complex))
z_gate = Tensor(np.array([[1,0],[0,-1]]))

# Bell states in the Bell basis (diagnonal)
phi_plus = Qstate(np.array([1,0,0,0]))
phi_minus = Qstate(np.array([0,1,0,0]))
psi_plus = Qstate(np.array([0,0,1,0]))
psi_minus = Qstate(np.array([0,0,0,1]))

# Bell states in the computational basis
phi_plus_c = Qstate(np.array([1,0,0,1]/np.sqrt(2)))
phi_minus_c = Qstate(np.array([1,0,0,-1]/np.sqrt(2)))
psi_plus_c = Qstate(np.array([0,1,1,0]/np.sqrt(2)))
psi_minus_c = Qstate(np.array([0,-1,1,0]/np.sqrt(2)))

# The quantum operation used by Alice to perform her measurements in the Bell basis
gate00 = phi_plus.braket(phi_plus_c)
gate10 = phi_minus.braket(phi_minus_c)
gate01 = psi_plus.braket(psi_plus_c)
gate11 = psi_minus.braket(psi_minus_c)
bell_gate = gate00+gate01+gate10+gate11

# The projection operator Alice uses to place particles two and three in the Psi minus Bell state

entangle_gate = psi_minus_c.braket(phi_plus)*identity
state = entangle_gate * initial_state()

# Apply the Bell operator to the particles one and two,
# leaving particle three unaffected
# In forming the Tensor product, list the operators
# in reverse order (last subspace/last qubits first)
state = identity * bell_gate * state

# Alice's measurement with four possible outcomes
for i,psi in enumerate([phi_plus, phi_minus, psi_plus, psi_minus]): #The four outcomes
    proj = psi.braket(psi)
    proj = identity * proj
    state1 = proj * state

    # Bob examines the state of third particle.
    # Simulate by forming density matrix and taking partial trace
    rho1= state1.braket(state1)
    rho1 = partial_trace_first_subsystem(rho1.val,2,1)
    state1 = get_states(rho1)[0]
    asym,bsym = sp.symbols('a b')

    # Depending on the measurement outcome (i), Bob selects and applies a quantum operation
    if i == 0: #phi_plus
        gate = y_gate
    elif i == 1: #phi_minus
        gate = x_gate
    elif i == 2: #psi_plus
        gate = z_gate
    else: # psi_minus
        gate = identity
    state2 = Qstate(get_states(rho1)[0])
    state2 = gate * state2

    # Verify that the final state is parallel to the initial state (cross prduct is zero)
    assert state2.val[0]*bsym-state2.val[1]*asym==0
   # Quantum repeater measures in the Bell basis.
   # Simlate by constructing the transformation from the Bell
# (diagonal) basis to the computational basis
gate00 = phi_plus_c.braket(phi_plus)
gate10 = phi_minus_c.braket(phi_minus)
gate01 = psi_plus_c.braket(psi_plus)
gate11 = psi_minus_c.braket(psi_minus)
bell_gate = gate00+gate10+gate01+gate11

# entangle qubits 1 and 2, leave 0, 3 and 4 unchanged
state = identity*identity*bell_gate*identity*initial_state2()

# entangle qubits 3 and 4, leave 0,1 and 2 unchanged
state = bell_gate*identity*identity*identity*state

# All possible outcomes for the measurement of qubits 1 and 2
for i1,psi in enumerate([phi_plus_c, phi_minus_c, psi_plus_c, psi_minus_c]):

    proj = identity*psi.braket(psi)*identity*identity
    state1 = proj*state

    # All possible outcomes for the measurement of qubits 0 and 1
    for i2, psi2 in enumerate([phi_plus_c, phi_minus_c, psi_plus_c, psi_minus_c]):

        proj = identity*identity*identity*psi2.braket(psi2)
        state2 = proj * state1

        # Form density matrix and extract state
        rho1= state2.braket(state2).val
        rho1 = partial_trace_first_subsystem(rho1,4,1)
        state2 = Qstate(get_states(rho1)[0])

        # Depending on the measurement outcomes, the receiver applies a gate to the link qubit.
        if i1==0:
            if i2 == 0:
                GATE = None
            elif i2 == 1:
                GATE = z_gate
            elif i2 == 2:
                GATE = x_gate
            else:
                GATE = y_gate
        elif i1==1:
            if i2 == 0:
                GATE = z_gate
            elif i2 == 1:
                GATE = None
            elif i2 == 2:
                GATE = y_gate
            else:
                GATE = x_gate

        elif i1 == 2:
            if i2 == 0:
                GATE = x_gate
            elif i2 == 1:
                GATE = y_gate
            elif i2 == 2:
                GATE = None
            else:
                GATE = z_gate
        elif i1 == 3:
            if i2 == 0:
                GATE = y_gate
            elif i2 == 1:
                GATE = x_gate
            elif i2 == 2:
                GATE = z_gate
            else:
                GATE = None
        if GATE is not None:
            state2 = GATE * state2

        # Verify that the intial qubit and the reconstructed,
        # repeated, qubit are parallel - zero cross product
        asym,bsym = sp.symbols('a b')
        assert state2.val[0]*bsym-state2.val[1]*asym==0
