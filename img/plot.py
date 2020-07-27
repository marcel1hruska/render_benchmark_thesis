import matplotlib.pyplot as plt
import csv

x = []
y = []
z = []
w = []

with open('color_matching_functions.txt','r') as csvfile:
    plots = csv.reader(csvfile, delimiter='\t')
    for row in plots:
        x.append(float(row[0]))
        y.append(float(row[1]))
        z.append(float(row[2]))
        w.append(float(row[3]))

plt.plot(x,y, label=r'$\bar{x}$')
plt.plot(x,z, label=r'$\bar{y}$')
plt.plot(x,w, label=r'$\bar{z}$')


plt.xlabel('Wavelength (nm)')
plt.title('CIE Color Matching Functions')
plt.legend()
plt.show()