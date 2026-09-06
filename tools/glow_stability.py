#!/usr/bin/env python3
"""Small-signal stability of the Grubert glow against a lumped ballast.

MODEL. The discharge carries I = G(n) V with G proportional to n_e, and the
density obeys  dn/dt = (nu_iz(V) - nu_loss) n.  Linearising about a state with
net growth rate gamma = nu_iz - nu_loss:

    dI = G dV + (I/n) dn ,      s dn = gamma dn + n nu' dV
    =>  Y_d(s) = G + I*nu'/(s - gamma)          nu' = d(nu_iz)/dV

The gap adds C_gap in parallel; the circuit adds Z(s) in series. The loop is

    1 + Z(s) * [ Y_d(s) + s C_gap ] = 0

Multiplying by (s - gamma) clears the pole and gives a polynomial whose roots
must all lie in the left half plane. Routh-Hurwitz on that polynomial is the
criterion.

gamma is the LINEARISATION POINT and it matters: gamma = 0 at the operating
point (steady state, by definition), gamma = 1/tau_g during the runaway.
"""
import math

# ---- measured operating-point and runaway quantities, 2026-09-06 -----------
Cg    = 1.77083756352e-16     # gap capacitance [F], from dischargeCurrent
Vop   = 200.0                 # gap voltage at the runaway [V]
Irun  = 1.25e-4               # conduction current at the runaway [A]
G     = Irun/Vop              # chord conductance [S]
tau_g = 0.68e-9               # measured ionisation e-folding [s]
gamma_run = 1.0/tau_g

# nu' = d(nu_iz)/dV, from the case's own tables:
#   1056 Td (meanE 16.03) -> k=1.694e-14 -> nu=4.09e8
#   2610 Td (meanE 37.16) -> k=6.075e-14 -> nu=1.467e9
# uniform-field V equivalents 190 V and 470 V
nu_p  = (1.467e9 - 4.09e8)/(470.0-190.0)

Iop   = 1.022e-6

def coeffs(R, L, Cext, gamma):
    """Polynomial coefficients, highest power first, of
       (s-gamma) + (R+sL)*[G(s-gamma) + I*nu' + sC(s-gamma)] = 0 ."""
    C = Cg + Cext
    I = Irun
    # (R + sL) * [ C s^2 + (G - gamma C) s + (I nu' - gamma G) ]  +  (s - gamma)
    A2, A1, A0 = C, G - gamma*C, I*nu_p - gamma*G
    # R * quadratic
    c3 = L*A2
    c2 = R*A2 + L*A1
    c1 = R*A1 + L*A0 + 1.0
    c0 = R*A0 - gamma
    return [c3, c2, c1, c0]

def stable(c):
    """Routh-Hurwitz. c = [c3,c2,c1,c0]; drop leading zeros for lower order."""
    while c and abs(c[0]) < 1e-300:
        c = c[1:]
    if len(c) == 3:                       # quadratic
        return all(x > 0 for x in c)
    if len(c) == 4:                       # cubic: all>0 AND c2*c1 > c3*c0
        return all(x > 0 for x in c) and c[1]*c[2] > c[0]*c[3]
    if len(c) == 2:
        return all(x > 0 for x in c)
    return False

def why(c):
    while c and abs(c[0]) < 1e-300: c = c[1:]
    bad=[f"c{len(c)-1-i}={v:.4g}" for i,v in enumerate(c) if v <= 0]
    if len(c)==4 and not (c[1]*c[2] > c[0]*c[3]):
        bad.append(f"c2*c1={c[1]*c[2]:.4g} !> c3*c0={c[0]*c[3]:.4g}")
    return ", ".join(bad) if bad else "-"

print(f"G = I/V = {G:.4g} S,  nu' = {nu_p:.4g} 1/(V s),  gamma_run = {gamma_run:.4g} 1/s")
print(f"C_gap = {Cg:.4g} F\n")

print("=== STEP 2: does the model REPRODUCE the two observed instabilities? ===")
print("    (linearised at the RUNAWAY state, gamma = 1/tau_g)")
for lbl,R,Cext in (("fix_n11  R=1e8, C=5e-14", 1e8, 5e-14),
                   ("fix_fast R=1e6, C=5e-14", 1e6, 5e-14)):
    c = coeffs(R, 0.0, Cext, gamma_run)
    s = stable(c)
    print(f"  {lbl:<26} -> {'STABLE  *** MODEL WRONG ***' if s else 'UNSTABLE  (observed: UNSTABLE)'}   [{why(c)}]")

print("\n=== STEP 3: can ANY series R / RL / RLC stabilise the RUNAWAY state? ===")
found=[]
for R in (1e6,1e7,1e8,1e9,1e10):
    for L in (0.0,1e-6,1e-3,1e-1,1e1,1e3):
        c = coeffs(R, L, 0.0, gamma_run)
        if stable(c): found.append((R,L))
print(f"  searched 5 R x 6 L with no shunt C: {len(found)} stable combinations")
if not found:
    c = coeffs(1e8, 1e-3, 0.0, gamma_run)
    print(f"  the blocking coefficient is the CONSTANT term:")
    print(f"    c0 = R*(I*nu' - gamma*G) - gamma")
    print(f"       I*nu'   = {Irun*nu_p:.4g}")
    print(f"       gamma*G = {gamma_run*G:.4g}")
    print(f"    I*nu' < gamma*G by {gamma_run*G/(Irun*nu_p):.2f}x, so c0 < 0 for ANY R and L.")
    print(f"    nu' > gamma/V is required: {nu_p:.4g} vs {gamma_run/Vop:.4g}  -> fails by {gamma_run/Vop/nu_p:.2f}x")
    print(f"    *** that condition contains NO circuit parameter. ***")

print("\n=== STEP 4: the same criteria AT THE OPERATING POINT (gamma = 0) ===")
for lbl,R,L,Cext in (("R=1e8, C=5e-14 (as run)", 1e8, 0.0, 5e-14),
                     ("R=1e8, C=0",              1e8, 0.0, 0.0),
                     ("R=1e8, L=1mH, C=0",       1e8, 1e-3, 0.0),
                     ("R=1e6, C=0",              1e6, 0.0, 0.0)):
    c = coeffs(R, L, Cext, 0.0)
    print(f"  {lbl:<26} -> {'STABLE' if stable(c) else 'UNSTABLE'}   [{why(c)}]")
print("\n  -> the TARGET state is stable. It is the PATH through gamma >> 0 that is not.")
