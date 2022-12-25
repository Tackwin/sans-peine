import numpy as np
import matplotlib.pyplot as plt
import magpylib as magpy
import math

fig, [ax1, ax2, ax3] = plt.subplots(1, 3, figsize=(10,5))

# create an observer grid in the xz-symmetry plane
ts = np.linspace(-2, 2, 60)
grid = np.array([[(x,0,z) for x in ts] for z in ts])

rings = 50
radia = 50
stars = 2
radius = 1
height = 1
r_in = 0.5
r_out = 1.0

# compute B- and H-fields of a cuboid magnet on the grid
cyl = magpy.magnet.CylinderSegment(magnetization=(0,0,1), dimension=(0.5, 1, 1, 0, 360))
B = cyl.getB(grid)
H = cyl.getH(grid)
cyl.rotate()

m_B = np.linalg.norm(B, axis=2)
# m_B = np.log(m_B)
# display field with Pyplot
# ax1.streamplot(grid[:,:,0], grid[:,:,2], B[:,:,0], B[:,:,2], density=2,
#     color=np.log(np.linalg.norm(B, axis=2)), linewidth=1, cmap='autumn')
ax1.imshow(m_B, cmap="viridis")

# ax2.streamplot(grid[:,:,0], grid[:,:,2], H[:,:,0], H[:,:,2], density=2,
#     color=np.log(np.linalg.norm(B, axis=2)), linewidth=1, cmap='winter')

print(grid.shape)
dipole_B = np.zeros(grid.shape)
dipole_H = np.zeros(grid.shape)

for r in range(rings):
	for s in range(stars):
		n = 1
		for n in range(radia):
			R = r_in + ((n + 0.5) / radia) * (r_out - r_in)
			p = (
				R * math.cos(2 * math.pi * s / stars),
				R * math.sin(2 * math.pi * s / stars),
				-height / 2 + height * (r + 0.5) / rings
			)
			dip = magpy.misc.Dipole(moment=(0, 0, 0.01 / (rings * stars * radia)), position=p)

			dipole_B += dip.getB(grid)
	# subcy = magpy.magnet.CylinderSegment(
	# 	magnetization=(0, 0, 1),
	# 	position=(0, 0, (r - (rings - 1) / 2) / rings),
	# 	dimension=(0.5, 1, 1 / rings, 0, 360)
	# )

	# dipole_B += subcy.getB(grid)
	# dipole_H += subcy.getH(grid)


mdipole_B = np.linalg.norm(dipole_B, axis=2)
# mdipole_B = np.log(mdipole_B)
# display field with Pyplot
# ax3.streamplot(grid[:,:,0], grid[:,:,2], dipole_B[:,:,0], dipole_B[:,:,2], density=2,
#     color=np.log(np.linalg.norm(B, axis=2)), linewidth=1, cmap='autumn')
ax2.imshow(mdipole_B, cmap="viridis")

# ax4.streamplot(grid[:,:,0], grid[:,:,2], dipole_H[:,:,0], dipole_H[:,:,2], density=2,
#     color=np.log(np.linalg.norm(B, axis=2)), linewidth=1, cmap='winter')

diff = np.linalg.norm(dipole_B - B, axis=2)
diff = np.abs(diff)

ax3.imshow(np.log(diff), cmap="viridis")


plt.tight_layout()
plt.show()