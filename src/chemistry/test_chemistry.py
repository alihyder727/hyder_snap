#! /usr/bin/env python3
#import matplotlib
#matplotlib.use('Agg')
from pylab import *

L = 1.

def svp(T):
  return exp(T);

def reaction_rate(q, params):
  k1, k2, k3, k4 = params[0], params[1], params[2], params[3]
  if q[1] >= svp(q[0]):
    r0 = (q[1] - svp(q[0]))*L  # temperature
    r1 = -k1*(q[1] - svp(q[0]))
    r2 = k1*(q[1] - svp(q[0])) - k3*q[2]*q[3] - k2*q[2];
    r3 = k3*q[2]*q[3] + k2*q[2]
  else:
    r0 = (q[1] - svp(q[0]))*L # temperature
    r1 =  k1*q[2]*(svp(q[0]) - q[1]) + k4*q[3]*(svp(q[0]) - q[1])
    r2 = -k1*q[2]*(svp(q[0]) - q[1]) - k3*q[2]*q[3] - k2*q[2]
    r3 = -k4*q[3]*(svp(q[0]) - q[1]) + k3*q[2]*q[3] + k2*q[2]
  return array([r0,r1,r2,r3])

def jacobian_matrix(q, params):
  jac = zeros((len(q),len(q)))
  k1, k2, k3, k4 = params[0], params[1], params[2], params[3]
  if q[1] >= svp(q[0]):
    jac = array([[-L*svp(q[0]),     L,          0.,       0.],
                 [ k1*svp(q[0]),  -k1,          0.,       0.],
                 [-k1*svp(q[0]),   k1, -k2-k3*q[3], -k3*q[2]],
                 [ 0.,             0.,  k2+k3*q[3],  k3*q[2]]])
  else:
    jac = array([[-L*svp(q[0]),     L,          0.,       0.],
                 [(k1*q[2]+k4*q[3])*svp(q[0]), -k4*q[3]-k1*q[2],  k1*(svp(q[0]) - q[1]), k4*(svp(q[0]) - q[1])],
                 [-k1*q[2]*svp(q[0]),   k1*q[2],  -k1*(svp(q[0]) - q[1])-k3*q[3]-k2,  -k3*q[2]],
                 [-k4*q[3]*svp(q[0]),   k4*q[3],   k3*q[3]+k2,  -k4*(svp(q[0]) - q[1])+k3*q[2]]])
  return jac

def solve_trbdf2(q, dt, params):
  gamma = 2. - sqrt(2.)
  jac = jacobian_matrix(q, params)
  A = eye(len(q)) - gamma/2.*dt*jac
  r0 = reaction_rate(q, params)
  b = dt*(1. - gamma/2.)*r0 + dt*gamma/2.*dot(A, r0)
  dq = linalg.solve(dot(A,A), b)
  return q+dq

def solve_trbdf2_new(q0, q1, dt, params):
  gamma = 2. - sqrt(2.)
  gamma3 = 1./(gamma*(2. - gamma))
  jac = jacobian_matrix(q1, params)
  A = eye(len(q0)) - gamma/2.*dt*jac
  r0 = reaction_rate(q0, params)
  r1 = reaction_rate(q1, params)
  b = dt*(1. - gamma/2.)*r0 + dt*gamma/2.*dot(A, r1) \
      + gamma3*(q0 - q1) + (1. - gamma3)*dot(A, q0 - q1)
  dq = linalg.solve(dot(A,A), b)
  return q1+dq

def solve_bdf1(q, dt, params):
  jac = jacobian_matrix(q, params)
  rate = reaction_rate(q, params)
  A = eye(len(q))/dt - jac
  dq = linalg.solve(A, rate)
  return q+dq

def truth_q0(T, q, t, k):
  return svp(T) + (q - svp(T))*exp(-k*t)

def truth_q1(t, k1, k2, k3):
  return k1/(k2 - k3)*(exp(-k3*t) - exp(-k2*t))

def truth_q2(t, k1, k2, k3):
  return k1/k2*1./(k2 - k3)*(k2*(1. - exp(-k3*t)) - k3*(1. - exp(-k2*t)))

if __name__ == '__main__':
  params = [1./0.1, 1./10., 0., 1./1.]
  #params = [1./0.1, 0., 0., 0.]
  time = logspace(-2, 3, 20)
  sol1, sol2, sol3 = [], [], []

  # condensation
  q = array([1., 3., 0., 0.])
  for t in time:
    #q1 = q.copy()
    #q2 = q.copy()
    q1 = solve_bdf1(q, t, params)
    #q2 = solve_trbdf2(q, t, params)
    print('1', q1)
    q2 = solve_trbdf2_new(q, q, t, params)
    print('2', q2)

    '''
    if q[1] > svp(q[0]):
      qsat = svp(q2[0])
      if q2[1] < qsat:
        dq = (q2[1] - qsat)/(1. + svp(q2[0])*L)
        q2[2] += dq
        q2[1] -= dq
        q2[0] += L*dq
        print(dq, q2[0], q2[1], svp(q2[0]))
    '''

    alpha = 1.;
    for i in range(len(q1)):
      if q2[i] < 0.:
        alpha = min(alpha, q1[i]/(q1[i] - q2[i]))
    #print(alpha, q1, q2)
    q3 = (1. - alpha)*q1 + alpha*q2

    sol1.append(q1)
    sol2.append(q2)
    sol3.append(q3)
  sol1 = array(sol1)
  sol2 = array(sol2)
  sol3 = array(sol3)

  fig, axs = subplots(2, 1, figsize = (12, 10), sharex = True)
  ax = axs[0]
  ax.plot([min(time), max(time)], [svp(q[0]), svp(q[0])], '-', color = '0.7', linewidth = 2)
  ax.plot(time, truth_q0(q[0], q[1], time, params[0]), 'k-', linewidth = 2)
  ax.plot(time, sol1[:,1], 'C1:', linewidth = 2, label = 'vapor, BDF')
  ax.plot(time, sol2[:,1], 'C1--', linewidth = 2, label = 'vapor, TR-BDF2')
  ax.plot(time, sol3[:,1], 'C1-', linewidth = 2, label = 'vapor, Blend')
  ax.legend()
  ax2 = ax.twinx()
  ax2.plot(time, sol1[:,0], 'C0:', linewidth = 2, label = 'T, BDF')
  ax2.plot(time, sol2[:,0], 'C0--', linewidth = 2, label = 'T, TR-BDF2')
  ax2.plot(time, sol3[:,0], 'C0-', linewidth = 2, label = 'T, Blend')


  ax = axs[1]
  ax.plot([min(time), max(time)], [0., 0.], '-', color = '0.7', linewidth = 2)
  ax.plot(time, truth_q1(time, params[0]*(q[1] - svp(q[0])), params[0], params[1]), 'k-', linewidth = 2)
  ax.plot(time, truth_q2(time, params[0]*(q[1] - svp(q[0])), params[0], params[1]), 'k-', linewidth = 2)
  ax.plot(time, sol1[:,2], 'C2:', linewidth = 2, label = 'cloud, BDF')
  ax.plot(time, sol2[:,2], 'C2--', linewidth = 2, label = 'cloud, TR-BDF2')
  ax.plot(time, sol3[:,2], 'C2-', linewidth = 2, label = 'cloud, Blend')
  ax.plot(time, sol1[:,3], 'C3:', linewidth = 2, label = 'precip, BDF')
  ax.plot(time, sol2[:,3], 'C3--', linewidth = 2, label = 'precip, TR-BDF2')
  ax.plot(time, sol3[:,3], 'C3-', linewidth = 2, label = 'precip, Blend')
  ax.legend()
  ax.set_xscale('log')
  ax.set_xlim([min(time), max(time)])

  show()
  #savefig('chemistry_integrator_latent.png', bbox_inches = 'tight')
