import csv
from matplotlib import pyplot as plt

x = []
y = []
reader = csv.DictReader(open('D:/Projects/git/K_NNG/K_NNG/points.csv'))

for row in reader:
    x.append(float(row['x']))
    y.append(float(row['y']))

plt.plot(x, y)

plt.xlabel('point.x')
plt.ylabel('point.y')
plt.title('line segments')

plt.show()