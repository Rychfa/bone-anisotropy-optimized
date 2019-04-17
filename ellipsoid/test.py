import numpy as np
import vector as vec
import matplotlib.pyplot as plt
import plotHelper as ph
import cvxpy as cvx

# construct set of points
N = 10

# randomly generate lengths
ls = np.random.rand(3)

# randomly generate rotation matrix
q = np.random.rand(3)*2*np.pi
R = vec.Rotation.from_rotation_vector(q)

# generate ellipsoid points
thetas = np.random.rand(N)*2*np.pi
phis   = np.random.rand(N)*2*np.pi

Bpoints = np.array([ls[0]*np.cos(thetas)*np.cos(phis),
                   ls[1]*np.cos(thetas)*np.sin(phis),
                   ls[2]*np.sin(thetas)]) # 3 x N

Ipoints = np.bmat([ R.body_to_inertial( vec.Vec3(Bpoints[:,i]) ).to_array() for i in range(N) ])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(Ipoints[0,:], Ipoints[1,:], Ipoints[2,:])


# algorithm
Q = np.matrix(np.eye(3))
alpha = 0.25
beta  = 0.5

def _cost(Q):
    res = 0
    for i in range(N):
        res = res + (Ipoints[:,i].T*Q*Ipoints[:,i] - 1)**2
    return res

while True:
    # determine step direction
    grad = np.matrix(np.zeros((3,3)))

    for i in range(N):
        residual_i = Ipoints[:,i].T*Q*Ipoints[:,i] - 1
        jacobian_i = Ipoints[:,i]*Ipoints[:,i].T

        grad = grad + 2*residual_i[0,0]*jacobian_i

    step = - grad

    # backtracking line search
    t = 1
    while _cost(Q+t*step) > _cost(Q) + alpha*t*np.trace(grad*step):
        t = t*beta

    step = step * t

    if np.linalg.norm(step) < 1e-3:
        break

    Q = Q + step

print("[Gradient Descent] Minimizer at ", Q)

evals, evecs = np.linalg.eig(Q)
for i in range(3):
    ax.plot([0, 1 / np.sqrt(evals[i]) * evecs[0, i]], [0, 1 / np.sqrt(evals[i]) * evecs[1, i]],
            [0, 1 / np.sqrt(evals[i]) * evecs[2, i]], c='r', linestyle='-')


Q = cvx.Variable((3,3))
error =   0
for i in range(N):
    residual = Ipoints[:,i].T*Q*Ipoints[:,i] - 1
    error = error + residual**2
objective = cvx.Minimize(error)
prob = cvx.Problem(objective)
prob.solve()

print("[CVXPY] Minimizer at ", Q.value)


evals, evecs = np.linalg.eig(Q.value)
for i in range(3):
    ax.plot([0, 1/np.sqrt(evals[i])*evecs[0,i]], [0,1/np.sqrt(evals[i])*evecs[1,i]], [0, 1/np.sqrt(evals[i])*evecs[2,i]], c='g', linestyle='--')

ph.axis_equal_3d(ax)