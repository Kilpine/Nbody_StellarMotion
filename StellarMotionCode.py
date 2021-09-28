# Original author Alexi Musick 9_20_21

# Check the readme file to have an overview of the programs function.

# Importing packages and libraries.

print ('Importing packages...')
import sys, os, gc, traceback
import pylab
import numpy as np
import matplotlib.pyplot as plt
import imageio # be sure that imageio is installed. Use 'pip install imageio'.

# We write out some functions to help with the processing of the movement of the stars.

def update(state,dt,masses):
    '''
    The function works to keep track of the past, current, and future states
    of the positions and velocities of the stars. The function utilizes basic
    kinematic equations to calculate the new postions and velocities from the
    old ones.
    '''
    r_old = state[0] # Old position pulled from the state tuple.
    v_old = state[1] # Old velocity pulled from the state tuple.
    r_new = (v_old * dt) + r_old
    acc = acceleration(r_old,masses) # Need to call acceleration in order to get the new velocity.
    v_new = 0.5 * acc * dt + v_old
    new_state = (r_new,v_new) # Make a new tuple that contains the new positions and velocities.
    return new_state

# We will also need the acceleration of the stars.

def acceleration(positions,masses):
    '''
    The function calculates the acceleration of the stars from the use
    of newton's gravity equation.
    '''
    accelerations = [] # Open list to append our accelerations to.
    for i,star_pos in enumerate(positions):
        total_acc = np.zeros(3) # A numpy array is made for the accelerations of the system.
        for j,star_pos2 in enumerate(positions):
            if i == j:
                continue
            r = star_pos - star_pos2
            rlen = np.linalg.norm(r)
            if rlen == 0.0:
                continue
            a = -G * masses[j] * r * (rlen ** -3.0)
            total_acc += a
        accelerations.append(total_acc)
    return np.array(accelerations)

# Lets write a function to calculate the Kinetic energy of the stars.

def kineticenergy(masses,velocities):
    '''
    The function calculates the kinetic energy of the
    system. The function is based off the classic 
    kinematic equation to calculate KE.
    '''
    mass_master = np.array([masses,masses,masses])
    # m = np.hstack(mass_master)
    K = 0.5 * np.sum(np.sum(mass_master.T*velocities**2))
    return K
            
# We also write a function to calculate the gravitational potential energy of the system.

def potentialenergy(positions,masses):
    '''
    The function calculates the potential energy of the indivual stars
    before adding them into a total. Each star feels a gravitational
    potential due to the effects of the other orbitting stars.
    '''
    P = 0 # Starting with an empty value which gets added to through the nested for loop.
    for i,star_pos in enumerate(positions):
        for j,star_pos2 in enumerate(positions):
            if i == j:
                continue
            r = star_pos - star_pos2
            rlen = np.linalg.norm(r)
            pe_ij = G * masses[i] * masses[j] / rlen # Gravitational potential.
            P += pe_ij  
    return P 

# A function to turn the .pngs of the plots into a singular .gif file.

# Now we define some of our global intital conditions and variables. These will be soon changed in simulation() based on the users input.

ti = 0.0 # Start time of the simulation.
tf = 0.0 # Final time of the simulation.
dt = 0.01 # Timestep through the simulation.
t = tf - ti # Total time of simulation.
Nt = int(np.ceil(t/dt)) # Number of time steps.
all_time = np.arange(Nt+2)*dt # All of the time steps. Plus 2 is needed to match shape of arrays. 
G = 39.42 # Gravitational Constant in terms of unit AU^3/M(solar)*yr^2

# Making a list of string names of the .png files to be saved.
    
png_to_convert = []

# We finally get to building the main 'simulation' and using the functions and variables we have written out above.

def simulation():
    '''
    The 'simulation' makes an array of the stellar history of the stars (including: position, velocity).
    With this array we then are easily able to plot out the values onto an xy graph.
    We also reference the total energy of the system to make sure there is no oddities.
    '''
    
    # Now we create a numpy array for the stars postions in terms of AU. These can be easily added to later. 
    
    positions = np.array([
        [-1,9,-1],
        [9,-1,1],
        [-11,-11,4],
        [4,-1,-6],
        [-1,4,4]
    ]) 
    
    # We writeout the x,y,z intial velocitites in terms of AU/yr of the stars into a numpy array.
    
    velocities = np.array([
        [-0.7,0.1,0],
        [0.3,1.1,0],
        [0.8,-0.4,0],
        [-0.7,0.1,0],
        [0.3,-0.9,0]
    ])
    
    # Let us let the user choose the timestep now.
    
    choosen_time = float(input('Decide time step of the system: '))
    tf = choosen_time
    t = tf - ti
    Nt = int(np.ceil(t/dt)) # Number of time steps.
    all_time = np.arange(Nt+2)*dt 
    
    # Since according to the homework mass is arbitrary we allow it to be choosen at will from an input function. These can easily be added to later.
    
    Nstar = len(positions) # The total number of stars.
    mass = float(input("Enter total solar mass of the system: ")) # Get user input for total mass of the system.
    masses = mass*np.ones(Nstar)/Nstar # Divides total mass into indivual masses for each star.
    
    # Here we write out some more useful information regarding the KE and PE energies of the system.
    
    KE = kineticenergy(masses,velocities)
    PE = potentialenergy(positions, masses)
    KE_history = [KE]
    PE_history = [PE]
    
    # Making a list that contains two tuples, postions and velocity. This will be updated with our update() function.
    
    star_history = [(positions,velocities)]
    
    # In the while loop below we want it to run throughout the time of the system from the beginning.
    
    current_time = ti
    
    # We have a while loop that updates the positions and velocities of the system and also the KE and PE.
    
    while current_time <= tf:
        
        current_state = star_history[-1]
        next_state = update(current_state,dt,masses)
        star_history.append(next_state)
        r_new = next_state[0]
        v_new = next_state[1]
        
        KE_history.append(kineticenergy(masses,v_new))
        PE_history.append(potentialenergy(r_new,masses))
        current_time += dt 
    
    # Now converting those KE and PE list into numpy arrays.
    
    KE_history = np.array(KE_history)
    PE_history = np.array(PE_history)
    TotalE_history = PE_history + KE_history
    
    # Now lets make the figure to be plotted onto.
    
    f = plt.figure(figsize = (20,20))
    ax1 = plt.subplot(3,1,1)
    ax2 = plt.subplot(3,1,2)
    ax3 = plt.subplot(3,1,3, projection = '3d')

    # Data for a three-dimensional line.
    
    zline = np.linspace(0, 15, 1000)

    # Setting aspect and limits on our graphs.
    
    ax1.set_aspect('equal', 'box')
    ax2.set_aspect('equal', 'box')
    ax3.set_aspect('auto','box')
    ax1.set_title('Stellar Motion', fontsize = 20)
    limits = (-55, 55)
    ax1.set(xlim = (limits),ylim = (limits),xlabel = ('X position'), ylabel = ('Y position'))
    ax2.set(xlim = (limits),ylim = (limits),xlabel = ('X position'), ylabel = ('Z position'))
    ax3.set_xlim3d(-100,100)
    ax3.set_ylim3d(-60,60)
    ax3.set_zlim3d(-50,50)
    
    # These numpy arrays exist to address an issue I was encounting in the shape of my normal arrays.
    
    histories = [{'x':[],'y':[],'z':[]} for i in range(Nstar)]
    colors = np.array([])
    size = np.array([])
    markers = np.array([])
    
    # The for loop plots the points from the star_history array and from KE_history and PE_history.
    
    for i,state in enumerate(star_history):
        
        # Getting the postions out of the earlier declared state tuple.
        
        positions = state[0]
        x = positions[:,0]
        y = positions[:,1]
        z = positions[:,2]
        
        for j in range(len(x)):
            [histories[j][n].append(m[j]) for n,m in zip(['x','y','z'],[x,y,z])]
        
    
    # Making the scatter plot of the stellar motion.
    
    colors = ['black','red','blue','purple','black']
    markers = ['o','s','^','^','s']
    facecolors = ['white','white','blue','white','black']
    
    for k in range(len(histories)):
        
        histories_k = histories[k]
        
        ax1.scatter(histories_k['x'][:-1],histories_k['y'][:-1],color = colors[k], s = 1)
        ax2.scatter(histories_k['x'][:-1],histories_k['z'][:-1],color = colors[k], s = 1)
        
        ax1.plot(histories_k['x'][-1],histories_k['y'][-1],color = colors[k], marker = markers[k],
                 markerfacecolor = facecolors[k], markeredgecolor = colors[k], markersize = 15)
        ax2.plot(histories_k['x'][-1],histories_k['z'][-1],color = colors[k], marker = markers[k],
                 markerfacecolor = facecolors[k], markeredgecolor = colors[k], markersize = 15)
        
        # We have to if statements to remove an issue with Nt not scaling correctly on certain timestep arrays.
        
        if t < 20.00 or t > 40.00:
            ax3.scatter3D(histories_k['x'], histories_k['y'], histories_k['z'], c = range(Nt+2), cmap='magma_r')
        if t >= 20.00 and t <= 40.00:
            ax3.scatter3D(histories_k['x'], histories_k['y'], histories_k['z'], c = range(Nt+1), cmap='magma_r')

    timestep_png = int(t * 100) # We declare a variable in order to quickly label our .png files.
    to_save = 'fivestar{}.png'.format(timestep_png)
    plt.savefig(to_save, dpi = 240)
    png_to_convert.append(to_save)
    print('Now creating plots at ' + str(timestep_png) + ' time steps:')
    plt.show()
    
    # Now making the scatter plots of the energies. Uncomment to run. 
    
    #f, ax = plt.subplots(1, 1, figsize=(6,4))
    #print(PE_history.shape, KE_history.shape, colors.shape, size.shape)
    #plt.scatter(all_time, KE_history, color = 'blue')
    #plt.scatter(all_time, PE_history, color = 'red')
    #plt.scatter(all_time,TotalE_history, color = 'green')
    #plt.show()

# Finally we get to running the simulation with a couple of different options for the user to choose from.

answer = str(input('Welcome to the program!\nBegin calculating?(Y/N): '))
if answer == 'Y' or answer == 'y':
    print('There are currently 5 stars in the simulation. Change base code to add or remove stars.')
    simulation()
else:
    print('Goodbye!')

# Now we want to give the user an oppurtunity to plot at different time steps.

multiplotting_answer =  str(input('Plot at different time steps? (Y/N): '))
if multiplotting_answer == 'Y' or multiplotting_answer == 'y':
        number_plots = int(input('How many more plots would you like to create?: '))
        
        for i in list(range(number_plots)):
            simulation()
else:
    print('Goodbye!')

# Now lets give the option 

png_to_gif_answer = str(input('Would you like to turn your .png files into a GIF? (Y/N): '))
if answer == 'Y' or answer == 'y':
    duration = 0.5
    name = 'fivestar.gif'
    create_gif(png_to_convert, duration, name)
    print('GIF Created in directory!\nHave a great day!')

def create_gif(filenames, durations, name):
    images = []
    for filename in filenames:
        images.append(imageio.imread(filename))
    output_file = name
    imageio.mimsave(output_file, images, format = 'GIF', duration = durations)