import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt

#[ Given coordinates for a divertor plate, parametrize it using
#[ cublic split interpolation to give Gkeyll a function f(t)
#[ where t\in[0,1] that describes the divertor plate.

def parametrize_curve(points):
  #[ Parametrizes a 2D curve using cubic spline interpolation.
  #[
  #[ Args:
  #[   points: A list or numpy array of 2D points, e.g., [(x1, y1), (x2, y2), ...].
  #[
  #[ Returns:
  #[   A tuple containing:
  #[     - t: An array of parameter values (from 0 to 1).
  #[     - x_t: An array of x-coordinates corresponding to t.
  #[     - y_t: An array of y-coordinates corresponding to t.

  points = np.array(points)
  distances = np.cumsum(np.sqrt(np.sum(np.diff(points, axis=0)**2, axis=1)))
  distances = np.insert(distances, 0, 0)
  
  #[ Normalize distances to the range [0, 1].
  t = distances / distances[-1]

  #[ Create cubic spline interpolations for x and y coordinates.
  cs_x = CubicSpline(t, points[:, 0])
  cs_y = CubicSpline(t, points[:, 1])

  #[ Evaluate the spline at more finely spaced t values.
  t_fine = np.linspace(0, 1, 100)
  x_t = cs_x(t_fine)
  y_t = cs_y(t_fine)
  print(t_fine)
  print(x_t)
  print(y_t)

  return t_fine, x_t, y_t
#[ End of parametrize_curve .........................]#

#[ Load the divertor plate coordinates.
points = np.loadtxt(open('./LTX_VV_103955_468ms.csv'),delimiter=',')
##[ Sometimes the divertor plate coordinates need sorting.
#points = np.sort(points, axis=1)
#points[:,0] = points[-1::-1,0]
#points[:,1] = points[-1::-1,1]
t, x_t, y_t = parametrize_curve(points)

#[ Plot the data and the parametrized curve.
plt.plot(np.array(points)[:, 0], np.array(points)[:, 1], 'o', label='Data')
plt.plot(x_t, y_t, '-', label='Parametrized Curve')
plt.xlabel(r'$R$ (m)')
plt.ylabel(r'$Z$ (m)')
plt.legend()
plt.grid(True)
plt.show()
