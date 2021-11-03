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
    Original Author Alexi Musick
    
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
    Original Author Alexi Musick
    
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
    Original Author Alexi Musick
    
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
    Original Author Alexi Musick
    
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
            pe_ij = - G * masses[i] * masses[j] / rlen # Gravitational potential.
            P += pe_ij  
    return P 

# A function to turn the .pngs of the plots into a singular .gif file.

'''
Stating intial values for our simulation including the intial time, final time, timesteps, number of steps, and
gravitational constant. For Homework 3 we need to put these time steps into intervals. There are two methods that come to mind
when thinking of how best to do this, slice indices or creating a mask. I think both can be used to make this problem very easy. 
I will also be making tf = 200.0 yr (tdyn(100 timesteps per tdyn)) to include all potential timesteps (20,000) the problem requires. This mask can be found in the
main() method.

'''
ti = 0.0 #t time of the simulation.
tf = 200.0 # Final time of the simulation.
dt = 0.01 # Timestep through the simulation (tdyn).
t = tf - ti # Total time of simulation.
Nt = int(np.ceil(t/dt)) # Number of time steps.
all_time = np.arange(Nt+2)*dt # All of the time steps. Plus 2 is needed to match shape of arrays. 
G = 39.42 # Gravitational Constant in terms of unit AU^3/M(solar)*yr^2

# Making a list of string names of the .png files to be saved.
png_to_convert = []

# We finally get to building the main 'main' and using the functions and variables we have written out above.
def main():
    '''
    Original Author Alexi Musick
    
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
    
    # Let us let the user choose the minimum and maximum for the mask.
    minimum = int(input('Decide a minimum for Tdyn interval: '))
    maximum = int(input('Decide a maximum for Tdyn interval: '))
    min_mask = minimum * 1000 
    max_mask = (maximum * 1000) - 1 # The subtraction here is to allow the last point of the graph to be the marker points.
    mask_length = max_mask - min_mask
    all_time_mask = range(mask_length + 1)
    
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
    limits = (-20,20)
    ax1.set(xlim = (limits),ylim = (limits),xlabel = ('X position [AU]'), ylabel = ('Y position [AU]'))
    ax2.set(xlim = (limits),ylim = (limits),xlabel = ('X position [AU]'), ylabel = ('Z position [AU]'))
    ax3.set_xlim3d(-20,20)
    ax3.set_ylim3d(-20,20)
    ax3.set_zlim3d(-20,20)
    ax3.set_xlabel('X position [AU]')
    ax3.set_ylabel('Y position [AU]')
    ax3.set_zlabel('Z position [AU]')
    
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
        
        ax1.scatter(histories_k['x'][min_mask:max_mask],histories_k['y'][min_mask:max_mask],color = colors[k], s = 1)
        ax2.scatter(histories_k['x'][min_mask:max_mask],histories_k['z'][min_mask:max_mask],color = colors[k], s = 1)
        
        ax1.plot(histories_k['x'][max_mask + 1],histories_k['y'][max_mask + 1],color = colors[k], marker = markers[k],
        markerfacecolor = facecolors[k], markeredgecolor = colors[k], markersize = 15)
        ax2.plot(histories_k['x'][max_mask + 1],histories_k['z'][max_mask + 1],color = colors[k], marker = markers[k],
        markerfacecolor = facecolors[k], markeredgecolor = colors[k], markersize = 15)
        
        # We have two if statements to remove an issue with Nt not scaling correctly on certain timestep arrays.
        ax3.scatter3D(histories_k['x'][min_mask:max_mask], histories_k['y'][min_mask:max_mask], 
        histories_k['z'][min_mask:max_mask], c = range(mask_length), cmap='hot_r')


    # A few lines to help save the figures we created.
    to_save = 'Stellar_Intervals_' + str(minimum) + '_to_' + str(maximum)
    plt.savefig(to_save, dpi = 300, bbox_inches='tight')
    png_to_convert.append(to_save)
    print('Now creating plots between timesteps: ' + str(minimum) + ' and ' + str(maximum))
    plt.show()
    
    # Energy plotting. Potential energy is perhaps too high but by all means looks to me to be written correctly using kinematic equations.
    f, ax4 = plt.subplots(1, 1, figsize=(6,4))
    plt.scatter(all_time_mask, KE_history[min_mask:max_mask + 1], color = 'plum', label = 'KE')
    plt.scatter(all_time_mask, PE_history[min_mask:max_mask + 1], color = 'mediumorchid', label = 'PE')
    plt.scatter(all_time_mask, TotalE_history[min_mask:max_mask + 1], color = 'sienna', label = 'Etotal')
    plt.xlabel('Time [yr]')
    plt.ylabel('Energy ' + r'$[AU^{3}/M_\odot*yr^{2}]$')
    plt.legend(loc="upper right")
    
    # A few lines to help save the energy graphs.
    to_save = 'Stellar_Energy_Intervals_' + str(minimum) + '_to_' + str(maximum)
    plt.savefig(to_save, dpi = 300, bbox_inches='tight')
    print('Creating energy graphs betweem timestps: ' + str(minimum) + ' and ' + str(maximum))
    plt.show()

# Finally we get to running the simulation with a couple of different options for the user to choose from.
answer = str(input('Welcome to the stellar motion program! V.2\nStart plotting?(Y/N): '))
if answer == 'Y' or answer == 'y':
    print('There are currently 5 stars in the simulation. Change base code to add or remove stars.')
    main()
else:
    print('Goodbye!')

# Now we want to give the user an oppurtunity to plot at different time steps.
multiplotting_answer =  str(input('Would you like to enter the stars at\ndifferent time intervals? (Y/N): '))
if multiplotting_answer == 'Y' or multiplotting_answer == 'y':
        number_plots = int(input('How many more plots would you like to create?: '))
        
        for i in list(range(number_plots)):
            main()
        png_to_gif_answer = str(input('Would you like to turn your .png files into a GIF? (Y/N): '))
        if png_to_gif_answer == 'Y' or png_to_gif_answer == 'y':
            duration = 0.5
            name = 'stellar_time_intervals.gif'
            create_gif(png_to_convert, duration, name)
            print('GIF Created in directory!\nHave a great day!')
else:
    print('Goodbye!')
