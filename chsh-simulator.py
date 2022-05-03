import numpy as np
from qutip import Bloch
def bloch_to_Hilbert (b, show=False):
    if b=='z':
        theta=0
        Phi=0
    elif b=='x':
        theta=90
        Phi=0
    elif b== 'z-x':
        theta=45
        Phi=180
    else:
        assert b=='z+x'
        theta=45
        Phi=0
    theta=np.radians(theta)
    phi=np.radians(Phi)
    psi1=[np.cos(theta/2),np.exp(1j*phi)*np.sin(theta/2)]
    psi2=[psi1[1],-psi1[0]]

    if show:
        b=Bloch(figsize=[8,8])
        vect=np.array([np.cos(theta)*np.cos(phi),np.cos(theta)*np.sin(phi),np.sin(theta)])
        b.add_vectors([vect,-vect])
        b.show()
        b.fig.suptitle("Basis")
        print('hello')
    return np.array(psi1),np.array(psi2)

wins1=[[0,0],[1,1]],[[0,0],[1,1]],[[0,0],[1,1]],[[0,1],[1,0]]
  
alice_bases=['z','x']
bob_bases=['z+x','z-x']
for ref_alice in 0,1:
  alice_basis = bloch_to_Hilbert(alice_bases [ref_alice])
  for ref_bob in 0,1:
    if ref_bob==1:
        bob_basis = bloch_to_Hilbert (bob_bases [ref_bob], show=True)
    else:
        bob_basis = bloch_to_Hilbert (bob_bases [ref_bob], show=False)
    wins = wins1[2*ref_alice+ref_bob]
    p=0
    for win in wins:
        alice_resp = win[0]
        bob_resp = win[1]
        alice_state = alice_basis[alice_resp]
        bob_state = bob_basis[bob_resp]
        p += (1/2)*np.matmul(np.conj(alice_state),bob_state)**2
    p = np.round(p,4)
    assert np.imag(p)==0
    print(f'win probability {np.real(p)}')