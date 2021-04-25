import pylab
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from matplotlib import cm
f = open('DataMain.txt', 'r')
try:
	# переправить
	nums = f.read().splitlines()

	#  Считали размерность
	number_x = int(nums[0])
	number_y = int(nums[1])

	x = []
	#  Считали границы по Х
	xBorder1 = float(nums[2])
	xBorder2 = float(nums[3])
	StepX = float(nums[4])
	x = np.arange(xBorder1, xBorder2 + StepX/100, StepX)
	# print(len(x))
	# print(x)

	y = [] 
	#  Считали границы по Y
	yBorder1 = int(nums[5])
	yBorder2 = int(nums[6])
	StepY = float(nums[7])
	y = np.arange(yBorder1, yBorder2 + StepY/100, StepY)
	# print(len(y))
	# print(y)

	grid = []
	for i in range(8,8 + int(number_x)*int(number_y)):
		grid.append(nums[i])

finally:
   f.close()

fig = pylab.figure()
axes = Axes3D(fig)

#  Создали из массивов двумерные сетки
X,Y = np.meshgrid(x, y)

#  Привели тип массива значений сетки к флоату
for i, elem in enumerate(grid):
    grid[i] = float(elem)

#  Преобразовали массив значений сетки из строк 
z = np.array(grid)

#  Задали ему размерность
Z = z.reshape(X.shape)

axes.plot_surface(X, Y, Z, rstride= 3, cstride=3, cmap = 'inferno' ) #  'coolwarm' cm.jet
axes.set_xlabel('X')
axes.set_ylabel('Y')
axes.set_zlabel('Z')
fig.canvas.set_window_title('Задача Дирихле для уравнения Пуассона - МВР')

pylab.show()
