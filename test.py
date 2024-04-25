import numpy as np

q1v4 = [-1.381815, -1.017865, 0.315646, 1.670802]
q2v4 = [-1.781527, -1.017865, 0.389931, 1.597989]
i = [-11.351995, -1.017865, 2.168544, 1.597989]
div = [-7.103927, -0.636967]

q1 = q1v4[:3]
q2 = q2v4[:3]
q1w = q1v4[3]
q2w = q2v4[3]

# print(q1)

p = [-q1w, 0, q1[2]]
n = np.cross([0, 1, 0], p)

d1 = np.dot(n, np.subtract(q1, p))
d2 = np.dot(n, np.subtract(q2, p))

t = d1 / (d1 - d2)

dir = np.subtract(q2, q1)

i = np.add(q1, dir*t)

print(i, end="")
print(" " + str(q2w))
print(i[0]/q2w, end=" ")
print(i[1]/q2w)
