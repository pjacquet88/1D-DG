import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

print('Animation is beginning')

if (len(sys.argv)==4):
    file_name=str(sys.argv[1])
    print_step=int(sys.argv[2])
    max_step=int(sys.argv[3])
    
path="../Files/"+file_name
frames=int((max_step/print_step)+1)

if file_name=='VP':
    xmin = 0
    xmax = 100
    ymin = 0.5
    ymax = 2.0
    interval=100

else:
    xmin = 0
    xmax = 1
    ymin = -0.5
    ymax = 0.5
    interval=20

fig = plt.figure() # initialise la figure
line, = plt.plot([],[]) 
plt.xlim(xmin, xmax)
plt.ylim(ymin,ymax)

# fonction à définir quand blit=True
# crée l'arrière de l'animation qui sera présent sur chaque image
def init():
    line.set_data([],[])
    return line,

def animate(i):
    current_file_name=path+str(i*print_step)+".dat"
    if file_name=='VP':
        plt.legend(['Model update iter n°'+str(i)])
        plt.xlabel('x')
        plt.ylabel('Velocity Model')

    else:
        plt.legend(['Pressure iter n°'+str(i)])
        plt.xlabel('x')
        plt.ylabel('Pressure')

    file=np.loadtxt(current_file_name)
    x = file[:,0]
    y = file[:,1]
    line.set_data(x, y)

    # if file_name=='VP':
    #     plt.legend("Velocity model "+str(i))
    # else:
    #     plt.legend("Pressure modeling "+str(i))
    # return line,


print('Animation created')
ani = animation.FuncAnimation(fig, animate, init_func=init, frames=frames, blit=False, interval=interval, repeat=True)

# print('Gif creation in progress')
# ani.save('test.gif',dpi=80,writer='pillow')
# print('Gif created')

plt.show()
