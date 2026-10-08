# Petri net reachability analysis of the PLC controlled pneumatic milling cell.
# Standard library only. A marking is a tuple of token counts over PLACES.
from collections import deque
CTRL = ['Idle','Clamp','SpinUp','Advance','Retract','SpinDown','Release','Fault']
PLANT = ['C_open','C_moving','C_closed','F_home','F_moving','F_fwd',
         'S_stopped','S_running','G_closed','G_open','P_ok','P_low']
PLACES = CTRL + PLANT
PERM = ('G_closed','P_ok')      # permissives tested by every energizing step
T = {}                          # name -> (input places, output places, test places)
def add(name, pre, post, test=()):
    T[name] = (tuple(pre), tuple(post), tuple(test))
# controller: normal progression
add('t_start',    ['Idle'],     ['Clamp'],    PERM + ('C_open','F_home','S_stopped'))
add('t_clamp_ok', ['Clamp'],    ['SpinUp'],   PERM + ('C_closed',))
add('t_spin_ok',  ['SpinUp'],   ['Advance'],  PERM + ('S_running',))
add('t_fwd_ok',   ['Advance'],  ['Retract'],  PERM + ('F_fwd',))
add('t_home_ok',  ['Retract'],  ['SpinDown'], PERM + ('F_home',))
add('t_stop_ok',  ['SpinDown'], ['Release'],  PERM + ('S_stopped',))
add('t_complete', ['Release'],  ['Idle'],     PERM + ('C_open',))
# controller: fault family from every active state
for s in CTRL[1:7]:
    add('f_%s_guard' % s, [s], ['Fault'], ('G_open',))
    add('f_%s_press' % s, [s], ['Fault'], ('P_low',))
    add('f_%s_abort' % s, [s], ['Fault'])          # movement deadline or operator Stop
add('f_Advance_spin', ['Advance'], ['Fault'], ('S_stopped',))
add('t_reset', ['Fault'], ['Idle'], PERM + ('C_open','F_home','S_stopped'))
# plant: clamp cylinder (command tested, completion uncontrolled)
add('t_ce1', ['C_open'],   ['C_moving'], PERM + ('Clamp',))
add('t_ce2', ['C_moving'], ['C_closed'])
add('t_cr2', ['C_moving'], ['C_open'])
add('t_cr1', ['C_closed'], ['C_moving'], PERM + ('Release',))
add('t_cr1_rec', ['C_closed'], ['C_moving'], PERM + ('Fault','F_home','S_stopped'))
# plant: feed carriage
add('t_fa1', ['F_home'],   ['F_moving'], PERM + ('Advance','C_closed','S_running'))
add('t_fa2', ['F_moving'], ['F_fwd'])
add('t_fb2', ['F_moving'], ['F_home'])
add('t_fr1', ['F_fwd'],    ['F_moving'], PERM + ('Retract',))
add('t_fr1_rec', ['F_fwd'], ['F_moving'], PERM + ('Fault',))
# plant: spindle
add('t_ss',     ['S_stopped'], ['S_running'], PERM + ('SpinUp',))
add('t_sc',     ['S_running'], ['S_stopped'], ('SpinDown',))
add('t_sc_rec', ['S_running'], ['S_stopped'], ('Fault',))
add('t_strip',  ['S_running'], ['S_stopped'])   # drive trips or coasts without command
# environment
add('t_g_open', ['G_closed'], ['G_open']);  add('t_g_close', ['G_open'], ['G_closed'])
add('t_p_drop', ['P_ok'], ['P_low']);       add('t_p_rise', ['P_low'], ['P_ok'])

idx = {p: i for i, p in enumerate(PLACES)}
def fire(m, t):
    pre, post, test = T[t]
    if any(m[idx[p]] < 1 for p in pre + test): return None
    n = list(m)
    for p in pre: n[idx[p]] -= 1
    for p in post: n[idx[p]] += 1
    return tuple(n)
M0 = tuple(1 if p in ('Idle','C_open','F_home','S_stopped','G_closed','P_ok') else 0 for p in PLACES)
seen = {M0: 0}; queue = deque([M0]); edges = []
while queue:
    m = queue.popleft()
    for t in T:
        n = fire(m, t)
        if n is None: continue
        if n not in seen: seen[n] = len(seen); queue.append(n)
        edges.append((m, t, n))
R = list(seen)
def has(m, *ps): return all(m[idx[p]] == 1 for p in ps)
def ctrl(m): return [p for p in CTRL if m[idx[p]]][0]
print('places', len(PLACES), 'transitions', len(T))
print('reachable markings', len(R), 'arcs', len(edges), 'bound', max(max(m) for m in R))
# P-invariants: a weight vector x with x.C = 0 keeps x.M constant in every reachable marking
C = {t: {p: 0 for p in PLACES} for t in T}
for t, (pre, post, test) in T.items():
    for p in pre:  C[t][p] -= 1
    for p in post: C[t][p] += 1
groups = {'controller': CTRL, 'clamp': PLANT[0:3], 'feed': PLANT[3:6],
          'spindle': PLANT[6:8], 'guard': PLANT[8:10], 'pressure': PLANT[10:12]}
for g, ps in groups.items():
    structural = all(sum(C[t][p] for p in ps) == 0 for t in T)
    print('P-invariant', g, 'x.C = 0:', structural, 'token sum', {sum(m[idx[p]] for p in ps) for m in R})
# T-invariant: the firing count vector of one nominal cycle satisfies C.y = 0
cycle = ['t_start','t_ce1','t_ce2','t_clamp_ok','t_ss','t_spin_ok','t_fa1','t_fa2','t_fwd_ok',
         't_fr1','t_fb2','t_home_ok','t_sc','t_stop_ok','t_cr1','t_cr2','t_complete']
print('T-invariant nominal cycle C.y = 0:', all(sum(C[t][p] for t in cycle) == 0 for p in PLACES),
      len(cycle), 'firings')
m = M0
for t in cycle: m = fire(m, t)
print('nominal cycle fires from M0 and returns to M0:', m == M0)
# safety properties checked over every reachable marking or arc
S1 = all(has(m,'C_closed') for m in R if not has(m,'F_home'))
S2 = all(has(m,'C_closed') for m in R if has(m,'S_running'))
S3 = all(has(m,'S_stopped') for m in R if has(m,'Release'))
S4 = all(has(m,'C_closed','S_running','G_closed','P_ok') for (m,t,n) in edges if t == 't_fa1')
S5 = all(has(m,'C_closed') for m in R if ctrl(m) in ('SpinUp','Advance','Retract','SpinDown'))
S6 = all(ctrl(m) in ('Advance','Retract','Fault') for m in R if not has(m,'F_home'))
print('S1 feed away from home => clamp closed', S1)
print('S2 spindle running => clamp closed', S2)
print('S3 Release => spindle stopped', S3)
print('S4 feed advance fires only with guard, pressure, clamp and spindle', S4)
print('S5 SpinUp, Advance, Retract, SpinDown => clamp closed', S5)
print('S6 feed away from home => Advance, Retract or Fault', S6)
# one step response: a permissive lost in an active state enables a fault transition at once,
# and no energizing step (controller progression or commanded movement) is enabled
lost = [m for m in R if ctrl(m) not in ('Idle','Fault') and (has(m,'G_open') or has(m,'P_low'))]
fault_ready = all(any(t.startswith('f_') for (mm,t,n) in edges if mm == m) for m in lost)
energizing = [t for t,(pre,post,test) in T.items() if 'G_closed' in test]
leak = sorted({t for m in lost for (mm,t,n) in edges if mm == m and t in energizing})
print('markings with a permissive lost in an active state', len(lost),
      '| fault transition enabled in all:', fault_ready, '| energizing steps enabled:', leak)
# deadlock freedom, dead transitions, reversibility (M0 is a home state, which implies liveness)
dead = [m for m in R if not any(fire(m, t) for t in T)]
fired = {t for (m,t,n) in edges}
rev = {}
for (m,t,n) in edges: rev.setdefault(n, []).append(m)
back = {M0}; q = deque([M0])
while q:
    x = q.popleft()
    for y in rev.get(x, []):
        if y not in back: back.add(y); q.append(y)
print('dead markings', len(dead), '| transitions never fired', sorted(set(T) - fired))
print('M0 reachable from every reachable marking (reversible):', len(back) == len(R))
print('markings per controller place', {s: sum(1 for m in R if ctrl(m) == s) for s in CTRL})
