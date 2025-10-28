import math
from matplotlib import pyplot as plt
from main import get_atoms 

atoms = get_atoms()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for a in atoms:
    # The s parameter is the area of the point, so using the a.rad we must calculate for the area
    point_area = (a.rad*2)*math.pi
    # point_area = a.rad
    ax.scatter(xs=a.coords.x, ys=a.coords.y, zs=a.coords.z, s=point_area, color=a.color)

ax.set_xlabel('X Axis')
ax.set_ylabel('Y Axis')
ax.set_zlabel('Z Axis')

plt.legend()
plt.show()